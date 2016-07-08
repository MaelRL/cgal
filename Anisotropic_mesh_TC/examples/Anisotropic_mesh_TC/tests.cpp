#define ANISO_COMPUTE_DUAL_FROM_TANGENT_PLANE
#define ANISO_TC_DEBUG

#include <CGAL/Epick_d.h>

#include <CGAL/Metric_field.h>
#include <CGAL/Full_cell_refine_queue.h>
#include <CGAL/Star.h>
#include <CGAL/Starset.h>

#include <CGAL/IO/Starset_output.h>

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <fstream>

const int d = 2;
const int D = 5;

using namespace CGAL::Anisotropic_mesh_TC;

typedef CGAL::Epick_d< CGAL::Dimension_tag<d> >         Kd;
typedef CGAL::Epick_d< CGAL::Dimension_tag<D> >         KD;

typedef Tangent_star<Kd, KD>                            Star;
typedef Star*                                           Star_handle;
typedef typename Star::FT                               FT;
typedef typename Star::Point_d                          Point_d;
typedef typename Star::Point_D                          Point_D;

typedef Full_cell_refine_queue<Kd, KD>                  Queue;
typedef typename Queue::Rfcell_set_iterator             Queue_iterator;

void test_transformations(const Star_handle star)
{
  std::vector<typename Star::FT> c(d, 1);
  typename Star::Point_d p = Kd().construct_point_d_object()(d, c.begin(), c.begin() + d);
  typename Star::WPoint_D ps = star->to_S(p);
  typename Star::WPoint_d pt = star->to_T(ps);

  std::cout << "p: " << p[0] << " " << p[1] << std::endl;
  std::cout << "ps: " << ps.point()[0] << " " << ps.point()[1] << " " << ps.point()[2] << " " << ps.point()[3] << " " << ps.point()[4] << std::endl;
  std::cout << "pt: " << pt.point()[0] << " " << pt.point()[1] << std::endl;
}

void read_points(std::vector<Point_d>& points)
{
  std::ifstream in("aniso_regular.mesh");
  std::string word;
  int useless, nv, dd;
  FT x;

  in >> word >> useless; //MeshVersionFormatted i
  in >> word >> dd; //Dimension d
  in >> word >> nv;

  assert(dd == Kd::Dimension::value);

  for(int i=0; i<nv; ++i)
  {
    std::vector<FT> ids;
    for(int j=0; j<dd; ++j)
    {
      in >> x;
      ids.push_back(x);
    }
    in >> useless;
    points.push_back(Kd().construct_point_d_object()(ids.begin(), ids.begin() + d));

    if(points.size() == 7)
      break;
  }

  std::cout << points.size() << " points" << std::endl;
}

bool is_simplex_degenerated_in_Rd(const typename Star::Full_cell_handle fch)
{
  const int d = Star::dDim::value;
  Eigen::Matrix<double, d, d> m;

  for(int i=0; i<d; ++i)
    for(int j=0; j<d; ++j)
      m(i,j) = fch->vertex(i+1)->point().point()[j]-fch->vertex(0)->point().point()[j];

  return (std::abs(m.determinant()) < 1e-5); //fixme not scalable
}

template<typename Starset>
void test_fch(const Starset& starset,
              const Star_handle star,
              const typename Starset::Full_cell_handle fch,
              Full_cell_refine_queue<Kd, KD>& refine_queue)
{
  std::cout << "test fch: ";
  std::cout << fch->vertex(0)->data() << " ";
  std::cout << fch->vertex(1)->data() << " ";
  std::cout << fch->vertex(2)->data() << std::endl;

  //if you want to use other criteria, fix the criteria computations first ! (todo)

  if(!starset.is_consistent(fch))
  {
    std::cout << "push inconsistency" << std::endl;
    FT vol = starset.compute_volume(fch);
    refine_queue.push(star, fch, vol, 5);
  }
  std::cout << "end test fch" << std::endl;
}

void fill_refinement_queue(const Starset<Kd, KD>& starset,
                           Full_cell_refine_queue<Kd, KD>& queue,
                           const int id = -1)
{
  typedef Starset<Kd, KD> Starset;
  typename Starset::const_iterator si = starset.begin();
  typename Starset::const_iterator siend = starset.end();
  for (; si != siend; si++)
  {
    if(id != -1 && (*si)->index() != id)
      continue;

    Star_handle star = *si;
    std::cout << "fill @ " << star->index() << std::endl;

    typename Starset::Full_cell_handle_iterator fchi = star->finite_incident_full_cells_begin();
    typename Starset::Full_cell_handle_iterator fend = star->finite_incident_full_cells_end();
    for(; fchi!=fend; ++fchi)
    {
      typename Starset::Full_cell_handle fch = *fchi;
      if((*fchi)->maximal_dimension() != d)
        continue;

      if(!is_simplex_degenerated_in_Rd(fch)
         //&& star->is_inside(fch, starset.stars())
         )
        test_fch(starset, star, fch, queue);
    }
  }

  std::cout << "nstars: " << starset.size() << std::endl;
  queue.print_queues();
}

bool next_refine_cell(Queue_iterator& it,
                      typename Star::Full_cell_handle& fch,
                      Queue& queue,
                      const Starset<Kd,KD>& starset)
{
  std::cout << "next refine face" << std::endl;
  while(true)
  {
    if(!queue.top(it))
    {
      std::cout << "empty queue at face pop time..." << std::endl;
      return false;
    }

    if(it->star->has_cell(fch, it->full_cell.vertices()))
    {
      //recompute the refinement point for fch, overkill, but to be safe
      if(it->star->compute_dual(fch, starset.stars()))
        return true;
      else
        queue.pop();
    }
    else
      queue.pop();
  }
}

bool refine(Starset<Kd, KD>& starset)
{
  Queue queue;
  fill_refinement_queue(starset, queue);

  while(!queue.empty())
  {
    Queue_iterator rfsit;
    typename Star::Full_cell_handle fch;
    next_refine_cell(rfsit, fch, queue, starset);
    assert(fch->data().second); // make sure the cell's circumcenter has been computed
    const Point_d& ref = fch->data().first;
    const Star_handle star = rfsit->star;

// VERBOSE ---------------------------------------------------------------------
    std::cout << "INSERT IN STARS: " << ref[0] << " " << ref[1] << std::endl;
    for(int i=0; i<=d; ++i)
    {
      typename Star::E_Vector_d v;
      Point_d pi = starset[fch->vertex(i)->data()]->m_center;
      for(int j=0; j<d; ++j)
        v(j) = ref[j]-pi[j];

      v = starset[fch->vertex(i)->data()]->metric().get_transformation() * v;
      FT dist  = std::sqrt(v.transpose() * v);
      std::cout << "metric dist: " << fch->vertex(i)->data() << " " << dist << std::endl;
    }
// VERBOSE ---------------------------------------------------------------------

    starset.insert_in_stars(ref);
    if(rfsit->star->has_cell(fch, rfsit->full_cell.vertices()))
    {
      std::cout << star->index() << " star has cell unbroken by refinement point ";
      std::cout << fch->vertex(0)->data() << " ";
      std::cout << fch->vertex(1)->data() << " ";
      std::cout << fch->vertex(2)->data() << std::endl;
      //assert(0);
    }
    else
    {
      std::cout << "broke the face" << std::endl;
      queue.pop();
    }

    std::ofstream outm("aniso_TC.mesh");
    output_medit(starset, outm);

    fill_refinement_queue(starset, queue, star->index());
    std::cout << "queue size: " << queue.count() << std::endl;
  }

  return true;
}

int main(int, char **)
{
  std::freopen("log_aniso_TC.txt", "w", stdout); // redirect std::cout

  //Custom_metric_field<Kd>* mf = new Custom_metric_field<Kd>();
  Euclidean_metric_field<Kd>* mf = new Euclidean_metric_field<Kd>();

  std::vector<Point_d> points;
  read_points(points);

  Starset<Kd, KD> starset(mf);

  for(std::size_t i=0; i<points.size(); ++i)
  {
    Star_handle si = new Star(points[i], i, mf->compute_metric(points[i]));
    starset.push_back(si);
  }

  starset.rebuild();
  refine(starset);

  std::ofstream out("aniso_TC.off");
  output_off(starset, out);
  std::ofstream outm("aniso_TC.mesh");
  output_medit(starset, outm);

  std::cout << std::endl;
}

#include <CGAL/Epick_d.h>

#include <CGAL/Metric_field.h>
#include <CGAL/Star.h>
#include <CGAL/Starset.h>

#include <CGAL/IO/Starset_output.h>

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <fstream>

#define dim_3
const int d = 3;
const int D = 9;

using namespace CGAL::Anisotropic_mesh_TC;

typedef CGAL::Epick_d< CGAL::Dimension_tag<d> >         Kd;
typedef CGAL::Epick_d< CGAL::Dimension_tag<D> >         KD;

typedef Tangent_star<Kd, KD>                            Star;
typedef Star*                                           Star_handle;
typedef typename Star::FT                               FT;
typedef typename Star::Point_d                          Point_d;
typedef typename Star::Point_D                          Point_D;

void test_transformations(const Star_handle star)
{
  std::vector<typename Star::FT> c(d, 1);
  typename Star::Point_d p = Kd().construct_point_d_object()(d, c.begin(), c.begin() + d);
  typename Star::WPoint_D ps = star->to_S(p);
  typename Star::WPoint_d pt = star->to_T(p);

  std::cout << "p: " << p[0] << " " << p[1] << std::endl;
  std::cout << "ps: " << ps.point()[0] << " " << ps.point()[1] << " " << ps.point()[2] << " " << ps.point()[3] << " " << ps.point()[4] << std::endl;
  std::cout << "pt: " << pt.point()[0] << " " << pt.point()[1] << std::endl;
}

void read_points(std::vector<Point_d>& points)
{
  std::ifstream in("aniso_regular.mesh");
  std::string word;
  int useless, nv, d;
  FT x;

  in >> word >> useless; //MeshVersionFormatted i
  in >> word >> d; //Dimension d
  in >> word >> nv;

  assert(d == Kd::Dimension::value);

  for(int i=0; i<nv; ++i)
  {
    std::vector<FT> ids;
    for(int j=0; j<d; ++j)
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

int main(int, char **)
{
  //std::freopen("log.txt", "w", stdout); // redirect std::cout to log.txt

  Custom_metric_field<Kd>* mgarbf = new Custom_metric_field<Kd>();
  Euclidean_metric_field<Kd>* mf = new Euclidean_metric_field<Kd>();

  std::vector<Point_d> points;
  read_points(points);

  Starset<Kd, KD> starset;

  for(std::size_t i=0; i<points.size(); ++i)
  {
    Star_handle si = new Star(points[i], i, mf->compute_metric(points[i]));
    starset.push_back(si);
  }

  starset.rebuild();

  std::ofstream out("aniso_TC.off");
  output_off(starset, out);
  std::ofstream outm("aniso_TC.mesh");
  output_medit(starset, outm);

  std::cout << std::endl;
}

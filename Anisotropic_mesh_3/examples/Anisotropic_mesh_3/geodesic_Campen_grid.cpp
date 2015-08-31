// this is based on the paper :
// Campen et al, Practical Anisotropic Geodesy EG 2013

// The canvas used is an orthogonal grid

#include <CGAL/Starset.h>
#include <CGAL/Stretched_Delaunay_3.h>

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <CGAL/Voronoi_painter.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <boost/array.hpp>
#include <Eigen/Dense>
#include <omp.h>

#include <ctime>
#include <iostream>
#include <limits>
#include <set>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K>
class Campen_canvas_point :
    public Canvas_point<K>
{
public:
  typedef boost::array<Campen_canvas_point*, 6>                      Neighbors;
  typedef Canvas_point<K>                                            Base;

  typedef typename Base::FT                                          FT;
  typedef typename Base::Point_3                                     Point_3;
  typedef typename Base::Metric                                      Metric;
  typedef typename Base::Vector3d                                    Vector3d;

  Neighbors neighbors;

  bool compute_closest_seed(const Base* anc)
  {
    // returns true if we improved the distance
    CGAL_assertion(anc->state == KNOWN);

    const int k = 8; // depth of the ancestor edge
    FT d = FT_inf;

    // the path is 'this' to ancestor1, ancestor1 to ancestor2, etc.
    // stored as 'this', ancestor1, ancestor2, etc.
    boost::array<const Base*, k+1> ancestor_path;
    for(int i=1; i<k+1; ++i)
      ancestor_path[i] = NULL;
    ancestor_path[0] = this;

    const Base* curr_ancestor = anc;
    for(int i=1; i<=k; ++i)
    {
      // add the new segment to the ancestor path
      ancestor_path[i] = curr_ancestor;

      Vector3d ancestor_edge;
      ancestor_edge(0) = this->point.x() - curr_ancestor->point.x();
      ancestor_edge(1) = this->point.y() - curr_ancestor->point.y();
      ancestor_edge(2) = this->point.z() - curr_ancestor->point.z();
      FT ancestor_edge_length = ancestor_edge.norm();
      Vector3d normalized_anc_edge = ancestor_edge / ancestor_edge_length;

      // compute the distance for the current depth (i)
      FT dist_to_ancestor = 0.;
      for(int j=0; j<i; ++j) // we add a part for each edge in the path
      {
        // get the metric for the current edge
        const Base* e0 = ancestor_path[j];
        const Base* e1 = ancestor_path[j+1];

        CGAL_assertion(e0 && e1);

        const Metric& m0 = e0->metric;
        const Metric& m1 = e1->metric;

        Vector3d curr_edge;
        curr_edge(0) = e0->point.x() - e1->point.x();
        curr_edge(1) = e0->point.y() - e1->point.y();
        curr_edge(2) = e0->point.z() - e1->point.z();

        // interpolate between both metric and transform the normalized edge
        // then we have (transformed_edge).norm() = || e ||_M = sqrt(e^t M e)
        Eigen::Matrix3d f = 0.5*(m0.get_transformation() + m1.get_transformation());
        Vector3d transformed_curr_edge = f*normalized_anc_edge;

        FT sp = curr_edge.dot(normalized_anc_edge);
        FT l = transformed_curr_edge.norm(); // length of the normalized anc edge in the metric

        dist_to_ancestor += sp * l;
      }
      dist_to_ancestor = (std::max)(dist_to_ancestor, 0.);

      // add ancestor edge length to the distance at that ancestor
      FT dist_at_anc = curr_ancestor->distance_to_closest_seed;
      FT new_d = dist_at_anc + dist_to_ancestor;

      if(new_d < d)
        d = new_d;

      if(!curr_ancestor->ancestor) // can't go any farther up in the ancestor tree
        break;

      curr_ancestor = curr_ancestor->ancestor;
    }

    if(d < this->distance_to_closest_seed)
    {
      this->ancestor = anc;
      this->distance_to_closest_seed = d;
      this->closest_seed_id = anc->closest_seed_id;
      return true;
    }

    return false;
  }

  PQ_state update_neighbors_distances(std::vector<Campen_canvas_point*>& trial_pq) const
  {
    // consider all the neighbors of a KNOWN point and compute their distance to 'this'
#if (verbose > 10)
    std::cout << "update neighbors of " << this->index << std::endl;
#endif
    CGAL_assertion(this->state == KNOWN);

    PQ_state pqs_ret = NOTHING_TO_DO;
    typename Neighbors::const_iterator it = neighbors.begin(),
                                       iend = neighbors.end();
    for(; it!=iend; ++it)
    {
      Campen_canvas_point* cp = *it;
      if(!cp)
        continue;
      else if(cp->state == KNOWN)
        continue; // dual_shenanigans(cp);
      else if(cp->state == TRIAL)
      {
        // note that we don't insert in trial_pq since it's already in
        if(cp->compute_closest_seed(this))
          pqs_ret = REBUILD_TRIAL;
      }
      else // cp->state == FAR
      {
        CGAL_assertion(cp->state == FAR);

        // note that cp->distance_to_closest_seed is not necessarily FT_inf here :
        // if we're refining, we've assigned FAR to all points after inserting a new
        // seed, therefore we must verify that compute_closest_seed is an update
        // before inserting it in the trial_queue
        if(cp->compute_closest_seed(this))
        {
          cp->state = TRIAL;
          trial_pq.push_back(cp);
          std::push_heap(trial_pq.begin(), trial_pq.end(),
                         Canvas_point_comparer<Campen_canvas_point>());
        }
      }
    }
    return pqs_ret;
  }

  template<typename MF>
  Campen_canvas_point(const Point_3& _point,
                      const std::size_t _index,
                      const MF _mf)
    :
      Base(_point, _index, _mf),
      neighbors()
  {
    for(std::size_t i=0; i<neighbors.size(); ++i)
      neighbors[i] = NULL;
  }
};

template<typename K, typename KExact, typename Canvas_point, typename Metric_field>
class Campen_canvas :
    public Canvas<K, KExact, Canvas_point, Metric_field>
{
public:
  typedef Canvas<K, KExact, Canvas_point, Metric_field>             Base;

  typedef typename Base::FT                                         FT;
  typedef typename Base::Point_3                                    Point_3;

  typedef typename Base::Simplex                                    Simplex;
  typedef typename Base::Edge                                       Edge;
  typedef typename Base::Tri                                        Tri;
  typedef typename Base::Tet                                        Tet;

  typedef typename Base::Vector3d                                   Vector3d;

  const Point_3 center;
  const FT side;
  const std::size_t n; // number of points per side

  FT offset_x; // offset is the bottom left point
  FT offset_y;
  FT offset_z;
  const FT step;
  const FT sq_n;

  void locate_and_initialize(const Point_3& p,
                             const std::size_t seed_id)
  {
    int index_x = std::floor((p.x()-offset_x)/step);
    int index_y = std::floor((p.y()-offset_y)/step);
    int index_z = std::floor((p.z()-offset_z)/step);

    Canvas_point* cp = &(this->points[index_z*sq_n + index_y*n + index_x]);

#if (verbose > 10)
    std::cout << "looking for p: " << p.x() << " " << p.y() << " " << p.z() << std::endl;
    std::cout << "found cp: " << cp->index << " [" << cp->point.x()
              << ", " << cp->point.y() << " " << cp->point.z() << "]" << std::endl;
    std::cout << "index: " << index_x << " " << index_y << " " << index_z << std::endl;
    std::cout << "offset: " << offset_x << " " << offset_y << " " << offset_z << std::endl;
#endif

    Vector3d v;
    v(0) = p.x() - cp->point.x();
    v(1) = p.y() - cp->point.y();
    v(2) = p.z() - cp->point.z();
    const Eigen::Matrix3d& m = cp->metric.get_mat();
    FT d = std::sqrt(v.transpose() * m * v);

    Base::initialize_canvas_point(cp, d, seed_id);
  }

  void initialize()
  {
#if (verbose > 5)
    std::cout << "canvas initialization" << std::endl;
#endif

    // create the canvas points
    for(unsigned int k=0; k<n; ++k)
    {
      for(unsigned int j=0; j<n; ++j)
      {
        for(unsigned int i=0; i<n; ++i)
        {
          Point_3 p(offset_x+i*step, offset_y+j*step, offset_z+k*step); // fill from bot left to top right
          Canvas_point cp(p, i + j*n + k*sq_n, this->mf);
          this->points.push_back(cp);
        }
      }
    }

    // assign the neighbors
    for(unsigned int k=0; k<n; ++k)
    {
      for(unsigned int j=0; j<n; ++j)
      {
        for(unsigned int i=0; i<n; ++i)
        {
          std::size_t curr_id = i + j*n + k*sq_n;
          if(k != n-1) // there is a neighbor above
            this->points[curr_id].neighbors[0] = &(this->points[i + j*n + (k+1)*sq_n]);
          if(i != 0) // there is a neighbor left
            this->points[curr_id].neighbors[1] = &(this->points[i-1 + j*n + k*sq_n]);
          if(j != n-1) // there is a neighbor back
            this->points[curr_id].neighbors[2] = &(this->points[i + (j+1)*n + k*sq_n]);
          if(i != n-1) // there is a neighbor right
            this->points[curr_id].neighbors[3] = &(this->points[i+1 + j*n + k*sq_n]);
          if(j != 0) // there is a neighbor front
            this->points[curr_id].neighbors[4] = &(this->points[i + (j-1)*n + k*sq_n]);
          if(k != 0) // there is a neighbor below
            this->points[curr_id].neighbors[5] = &(this->points[i + j*n + (k-1)*sq_n]);
        }
      }
    }
#if (verbose > 5)
    std::cout << "neighbors assigned" << std::endl;
#endif

    Base::compute_canvas_bbox();
    this->seeds.initialize_seeds();
    Base::locate_seeds_on_canvas();

#if (verbose > 5)
    std::cout << "canvas initialized" << std::endl;
#endif
  }

  void cosphericity_sweeper() const
  {
    // we define a tiny (2k+1)*(2k+1) grid ('sweeper'), and we look at the colors
    // if we have more than 5 colors at a time, we have a (quasi-) cosphericity

    std::set<boost::array<std::size_t, 5> > quasi_cosphericities;

    const int h = (std::max)(static_cast<std::size_t>(1), n/100);
    const std::size_t sweeper_size = 2*h+1;
    if(sweeper_size > n)
    {
      std::cerr << "WARNING: sweeper too big for that canvas !" << std::endl;
      return;
    }

    // we move the center of the sweeper in the canvas
    for(std::size_t i=h; i<n-h; ++i)
    {
      for(std::size_t j=h; j<n-h; ++j)
      {
        for(std::size_t k=h; k<n-h; ++k)
        {
          //const Canvas_point& sweeper_center = points[j*n + i];

          // gather the pts
          std::set<std::size_t> colors;
          for(int lx=-h; lx<=h; ++lx)
          {
            for(int ly=-h; ly<=h; ++ly)
            {
              for(int lz=-h; lz<=h; ++lz)
              {
                const Canvas_point& cp = this->points[(k+lz)*sq_n + (j+ly)*n + (i+lx)];
                colors.insert(cp.closest_seed_id);
              }
            }
          }

          if(colors.size() > 4)
          {
            // we have at least 5 pts in a quasi cospherical configuration...
            // get them all with combinations!

            std::vector<boost::array<std::size_t, 5> > combis =
                combinations<5>(colors, colors.begin(), 5);

            quasi_cosphericities.insert(combis.begin(), combis.end());
          }
        }
      }
    }

    std::cout << "Sweeper no sweeping ! (with size: " << h << ")" << std::endl;
    std::cout << "Found " << quasi_cosphericities.size() << " (quasi-)cosphericities" << std::endl;

    std::set<boost::array<std::size_t, 5> >::const_iterator sit =
                                                   quasi_cosphericities.begin();
    std::set<boost::array<std::size_t, 5> >::const_iterator send =
                                                     quasi_cosphericities.end();
    for(; sit!=send; ++sit)
    {
      const boost::array<std::size_t, 5>& quasi_cosphericity = *sit;
      for(std::size_t i=0; i<sit->size(); ++i)
        std::cout << quasi_cosphericity[i] << " ";
      std::cout << std::endl;
    }
    std::cout << "End of quasi cosphericities" << std::endl;
  }

  void add_tet_to_simplices(std::set<Tet>& simplices,
                            const Canvas_point& cp,
                            const Canvas_point* n1,
                            const Canvas_point* n2,
                            const Canvas_point* n3) const
  {
    Tet t;
    if(n1 && n2 && n3)
    {
      t[0] = cp.index;
      t[1] = n1->index;
      t[2] = n2->index;
      t[3] = n3->index;
      std::sort(t.begin(), t.end()); // is that really needed ?
      simplices.insert(t);
    }
  }

  void add_tet_to_simplices(std::set<Tet>& simplices, const Canvas_point& cp,
                            std::size_t i1, std::size_t i2, std::size_t i3) const
  {
    return add_tet_to_simplices(simplices, cp, cp.neighbors[i1], cp.neighbors[i2], cp.neighbors[i3]);
  }

  void output_canvas(const std::string str_base) const
  {
#if (verbose > 0)
    std::cout << "Output canvas" << std::endl;
#endif

    std::ofstream out((str_base + ".mesh").c_str());
    std::ofstream out_bb((str_base + ".bb").c_str());
    out << "MeshVersionFormatted 1" << std::endl;
    out << "Dimension 3" << std::endl;
    out << "Vertices" << std::endl;
    out << this->points.size() << std::endl;
    out_bb << "3 1 " << this->points.size() << " 2" << std::endl;

    std::set<Tet> simplices;
    for(std::size_t i=0; i<this->points.size(); ++i)
    {
      const Canvas_point& cp = this->points[i];
      out << cp.point.x() << " " << cp.point.y() << " " << cp.point.z() << " " << i+1 << std::endl;

//      out_bb << cp.closest_seed_id << std::endl;
      out_bb << cp.distance_to_closest_seed << std::endl;

      // compute the tets...
      const boost::array<Canvas_point*, 6>& ns = cp.neighbors; // just to get a shorter name...
      for(std::size_t i=0; i<ns.size(); ++i)
      {
#if 0
        // lazily : the eight corner tetrahedra (creates self intersections,
        // but it doesn't really matter except visually)
        add_tet_to_simplices(simplices, cp, 0, 1, 2);
        add_tet_to_simplices(simplices, cp, 0, 2, 3);
        add_tet_to_simplices(simplices, cp, 0, 3, 4);
        add_tet_to_simplices(simplices, cp, 0, 4, 1);
        add_tet_to_simplices(simplices, cp, 5, 1, 2);
        add_tet_to_simplices(simplices, cp, 5, 2, 3);
        add_tet_to_simplices(simplices, cp, 5, 3, 4);
        add_tet_to_simplices(simplices, cp, 5, 4, 1);
#else
        // a bit smarter : a cube is decomposed in 5 tets and we pick them
        // correctly so we have no overlap
        add_tet_to_simplices(simplices, cp, 0, 2, 3);
        add_tet_to_simplices(simplices, cp, 0, 1, 4);
        add_tet_to_simplices(simplices, cp, 3, 4, 5);
        add_tet_to_simplices(simplices, cp, 1, 2, 5);

        // the last one is a bit annoying since we have to grab neighbors of neighbors
        const Canvas_point* n4 = cp.neighbors[4];
        const Canvas_point* n3 = cp.neighbors[3];
        if(n4 && n3)
          add_tet_to_simplices(simplices, cp, n4->neighbors[3],
                               n4->neighbors[0], n3->neighbors[0]);
#endif
      }
    }

    out << "Tetrahedra" << std::endl;
    out << simplices.size() << std::endl;
    for(typename std::set<Tet>::iterator it = simplices.begin();
                                             it != simplices.end(); ++it)
    {
      const Tet& tet = *it;
      std::set<std::size_t> materials;

      for(std::size_t i=0; i<tet.size(); ++i)
      {
        out << tet[i] + 1 << " ";
        materials.insert(this->points[tet[i]].closest_seed_id);
      }
      std::size_t mat = (materials.size()==1)?(*(materials.begin())):(this->seeds.size());
      out << mat << std::endl;
    }

    std::set<Tri> triangles;
    for(typename std::set<Tet>::iterator it = simplices.begin();
                                         it != simplices.end(); ++it)
    {
      const Tet& tet = *it;
      Tri tr; // could use a permutation to make it look nicer, I suppose...
      tr[0] = tet[0]; tr[1] = tet[1]; tr[2] = tet[2]; triangles.insert(tr);
      tr[0] = tet[0]; tr[1] = tet[1]; tr[2] = tet[3]; triangles.insert(tr);
      tr[0] = tet[0]; tr[1] = tet[2]; tr[2] = tet[3]; triangles.insert(tr);
      tr[0] = tet[1]; tr[1] = tet[2]; tr[2] = tet[3]; triangles.insert(tr);
    }

    out << "Triangles" << std::endl;
    out << triangles.size() << std::endl;
    for(typename std::set<Tri>::iterator it = triangles.begin();
                                         it != triangles.end(); ++it)
    {
      const Tri& tr = *it;
      std::set<std::size_t> materials;
      for(std::size_t i=0; i<tr.size(); ++i)
      {
        out << tr[i] + 1 << " ";
        materials.insert(this->points[tr[i]].closest_seed_id);
      }
      std::size_t mat = (materials.size()==1)?(*(materials.begin())):(this->seeds.size());
      out << mat << std::endl;
    }
  }

  void compute_dual() const
  {
#if (verbose > 5)
    std::cout << "dual computations" << std::endl;
#endif

    for(std::size_t i=0; i<this->points.size(); ++i)
    {
      const Canvas_point& cp = this->points[i];
      Simplex dual_simplex;
      dual_simplex.insert(cp.closest_seed_id);

      typename Canvas_point::Neighbors::const_iterator it = cp.neighbors.begin(),
                                                       iend = cp.neighbors.end();
      for(; it!=iend; ++it)
      {
        const Canvas_point* cq = *it;
        if(!cq)
          continue;

        dual_simplex.insert(cq->closest_seed_id);
      }

      std::cout << "dual simp: " << dual_simplex.size() << std::endl;
      Base::add_simplex_to_triangulation(&cp, dual_simplex);
    }
  }

  Campen_canvas(const std::string& _canvas_str,
                const std::string& _seeds_str,
                const Point_3& _center,
                const FT _side,
                const std::size_t points_per_side,
                const std::size_t _max_seeds_n,
                const std::size_t _n_refine,
                const Metric_field _mf)
    :
      Base(_canvas_str, _seeds_str, _max_seeds_n, _n_refine, _mf),
      center(_center),
      side(_side),
      n(points_per_side),
      offset_x(center.x() - side/2.),
      offset_y(center.y() - side/2.),
      offset_z(center.z() - side/2.),
      step(side / (n-1)),
      sq_n(n*n)
  { }
};

} // namespace Anisotropic_mesh_3
} // namespace CGAL

using namespace CGAL::Anisotropic_mesh_3;

int main(int, char**)
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
  typedef CGAL::Exact_predicates_exact_constructions_kernel    KExact;

  typedef typename K::FT                                       FT;
  typedef typename K::Point_3                                  Point_3;

  typedef typename CGAL::Anisotropic_mesh_3::Euclidean_metric_field<K>* MF;
  //typedef typename CGAL::Anisotropic_mesh_3::Custom_metric_field<K>* MF;

  typedef Campen_canvas_point<K>                               Campen_canvas_point;
  typedef Campen_canvas<K, KExact, Campen_canvas_point, MF>    Canvas;

  std::cout.precision(17);
//  std::freopen("log.txt", "w", stdout);

  double duration;
  std::clock_t start = std::clock();

  MF mf = new Euclidean_metric_field<K>(1., 1., 3.);
//  MF mf = new Custom_metric_field<K>();

  const std::string canvas_str = "grid";
  const std::string seeds_str = "base_mesh.mesh";
  std::size_t max_seeds_n = 1;
  std::size_t n_refine = 0;

  // canvas geometry
  Point_3 center(1., 1., 1.);
  const FT canvas_side = 2.;
  FT points_per_side = 50.; // number of points per side of the canvas

  Canvas canvas(canvas_str, seeds_str,
                center, canvas_side, points_per_side,
                max_seeds_n, n_refine,
                mf);

  canvas.initialize();
  canvas.paint();
  canvas.refine();
  canvas.output_canvas_data_and_dual(canvas_str + "_tr");
//  canvas.cosphericity_sweeper();

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cerr << "Total: " << duration << std::endl;
}

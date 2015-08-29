// this is based on the paper :
// Campen et al, Practical Anisotropic Geodesy EG 2013

// the canvas used is a triangulation

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <CGAL/Voronoi_painter.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <Eigen/Dense>
#include <omp.h>
#include <boost/array.hpp>
#include <boost/tuple/tuple.hpp>

#include <iostream>
#include <map>
#include <vector>
#include <set>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

// The following classes implement some functions specific to the Campen algorithm

template<typename K>
class Campen_canvas_point :
    public Canvas_point<K>
{
public:
  typedef boost::unordered_set<Campen_canvas_point<K> *>             Neighbors;
  typedef Canvas_point<K>                                            Base;

  typedef typename Base::FT                                          FT;
  typedef typename Base::Point_3                                     Point_3;
  typedef typename Base::Metric                                      Metric;
  typedef typename Base::Vector3d                                    Vector3d;

  Neighbors neighbors;

  // this function is the heart of the painter
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
        cp->compute_closest_seed(this); // always returns true here so no need to test it
        cp->state = TRIAL;
        trial_pq.push_back(cp);
        std::push_heap(trial_pq.begin(), trial_pq.end(),
                       Canvas_point_comparer<Campen_canvas_point>());
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
  { }
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

  void initialize()
  {
#if (verbose > 5)
    std::cout << "canvas initialization" << std::endl;
#endif

    // read and create the canvas points
    std::ifstream in((this->canvas_str + ".mesh").c_str());
    std::string word;
    std::size_t useless, nv, ne, nt, ntet, dim;
    FT r_x, r_y, r_z;
    int i1, i2, i3, i4;

    in >> word >> useless; // MeshVersionFormatted i
    in >> word >> dim; // Dimension d
    CGAL_assertion(dim == 3);

    while(in >> word)
    {
      if(word == "End")
        break;
      else if(word == "Vertices")
      {
        in >> nv;
        std::cout << "canvas made of " << nv << " vertices" << std::endl;

        for(std::size_t i=0; i<nv; ++i)
        {
          in >> r_x >> r_y >> r_z >> useless;
          Point_3 p(r_x, r_y, r_z);
          Canvas_point cp(p, i, this->mf);
          this->points.push_back(cp);
        }
      }
      else if(word == "Edges")
      {
        //ignoring edges for now... could find a use for them... ? todo
        in >> ne;
        for(std::size_t i=0; i<ne; ++i)
          in >> useless >> useless >> useless;
      }
      else if(word == "Triangles")
      {
        //ignoring triangles for now... could find a use for them... ? todo
        in >> nt;
        std::cout << "canvas made of " << nt << " triangles" << std::endl;

        for(std::size_t i=0; i<nt; ++i)
          in >> useless >> useless >> useless >> useless;
      }
      else if(word == "Tetrahedra")
      {
        // read the tetrahedra and assign the neighbors
        in >> ntet;
        std::cout << "canvas made of " << ntet << " tets" << std::endl;

        for(std::size_t i=0; i<ntet; ++i)
        {
          in >> i1 >> i2 >> i3 >> i4 >> useless;
          --i1; --i2; --i3; --i4; // '--' because we're reading medit data...

          this->points[i1].neighbors.insert(&this->points[i2]);
          this->points[i1].neighbors.insert(&this->points[i3]);
          this->points[i1].neighbors.insert(&this->points[i4]);

          this->points[i2].neighbors.insert(&this->points[i1]);
          this->points[i2].neighbors.insert(&this->points[i3]);
          this->points[i2].neighbors.insert(&this->points[i4]);

          this->points[i3].neighbors.insert(&this->points[i1]);
          this->points[i3].neighbors.insert(&this->points[i2]);
          this->points[i3].neighbors.insert(&this->points[i4]);

          this->points[i4].neighbors.insert(&this->points[i1]);
          this->points[i4].neighbors.insert(&this->points[i2]);
          this->points[i4].neighbors.insert(&this->points[i3]);

          Tet tet;
          tet[0] = i1; tet[1] = i2; tet[2] = i3; tet[3] = i4;
          this->tetrahedra.push_back(tet);
        }
      }
    }

    Base::compute_canvas_bbox();
    this->seeds.initialize_seeds();
    Base::locate_seeds_on_canvas();

#if (verbose > 5)
    std::cout << "canvas initialized" << std::endl;
#endif
  }

  void check_edelsbrunner()
  {
    std::cout << "it's Edel's time" << std::endl;

    // for each simplex, associate all the tetrahedra (vector of size_t, the size_t
    // is the index of the tet in 'tetrahedra') of the base canvas that correspond
    // to that simplex...
    // Then check that this tet set satisfies the closed ball property
    typedef std::map<Simplex, std::vector<std::size_t> > Dual_map;

    Dual_map dual_map;
    int failed_blob_counter = 0;

    for(std::size_t i=0; i<this->tetrahedra.size(); ++i)
    {
      Simplex dual_simplex;
      const Tet& tet = this->tetrahedra[i];

      for(int j=0; j<4; ++j)
        dual_simplex.insert(this->points[tet[j]].closest_seed_id);

      std::pair<typename Dual_map::iterator, bool> is_insert_successful;
      std::vector<std::size_t> vec;
      vec.push_back(i);
      is_insert_successful = dual_map.insert(std::make_pair(dual_simplex, vec));

      if(!is_insert_successful.second)
        (is_insert_successful.first)->second.push_back(i);
    }

    // now, for each simplex, we want the tetrahedra corresponding to that simplex
    // to be a single blob

    typename Dual_map::const_iterator dmit = dual_map.begin();
    typename Dual_map::const_iterator dmend = dual_map.end();
    for(; dmit!=dmend; ++dmit)
    {
      const Simplex& simplex = dmit->first;
      std::cout << "looking at the simplex : ";
      typename Simplex::iterator sit = simplex.begin();
      typename Simplex::iterator send = simplex.end();
      for(; sit!=send; ++sit)
        std::cout << *sit << " ";
      std::cout << std::endl;

      const std::vector<std::size_t>& canvas_tetrahedra = dmit->second;

      // collect all the triangles that appear once
      std::set<Tri> triangles;

      for(std::size_t i=0; i<canvas_tetrahedra.size(); ++i)
      {
        const Tet& tr = this->tetrahedra[canvas_tetrahedra[i]];
        for(int j=0; j<4; ++j)
        {
          Tri triangle;
          triangle[0] = tr[(j+1)%4];
          triangle[1] = tr[(j+2)%4];
          triangle[2] = tr[(j+3)%4];
          // gotta sort since we'll use a map with edge keys
          std::sort(triangle.begin(), triangle.end());

          std::pair<typename std::set<Tri>::iterator, bool> is_insert_successful;
          is_insert_successful = triangles.insert(triangle);

          // already in 'triangles', thus triangle is internal & ignore it
          if(!is_insert_successful.second)
            triangles.erase(is_insert_successful.first);
        }
      }

      std::cout << triangles.size() << " triangles" << std::endl;
      CGAL_assertion(!triangles.empty());

      // check that all the border triangles of the blob actually form a blob
      // note: 'bow-ties' configurations could appear... what to do ?
      std::map<const Tri*, bool> visited_status;
      std::map<Edge, std::deque<const Tri*> > incident_tris; // edge to tri

      typename std::set<Tri>::iterator it = triangles.begin();
      typename std::set<Tri>::iterator iend = triangles.end();
      for(; it!=iend; ++it)
      {
        const Tri& tr = *it;

        std::pair<typename  std::map<Edge, std::deque<const Tri*> >::iterator, bool>
                                                           is_insert_successful;

        std::deque<const Tri*> vec;
        vec.push_back(&tr);

        Edge e;
        e[0] = tr[0]; e[1] = tr[1];
        is_insert_successful = incident_tris.insert(std::make_pair(e, vec));
        if(!is_insert_successful.second)
          is_insert_successful.first->second.push_back(&tr);

        e[0] = tr[0]; e[1] = tr[2];
        is_insert_successful = incident_tris.insert(std::make_pair(e, vec));
        if(!is_insert_successful.second)
          is_insert_successful.first->second.push_back(&tr);

        e[0] = tr[1]; e[1] = tr[2];
        is_insert_successful = incident_tris.insert(std::make_pair(e, vec));
        if(!is_insert_successful.second)
          is_insert_successful.first->second.push_back(&tr);

        visited_status[&tr] = false;
      }

      CGAL_assertion(!incident_tris.empty());

      std::deque<const Tri*> triangles_to_visit;
      triangles_to_visit.push_back(visited_status.begin()->first);
      while(!triangles_to_visit.empty())
      {
//        std::cout << "triangles to visit: " << triangles_to_visit.size() << std::endl;
        const Tri& current_tri = *(triangles_to_visit.front());
        triangles_to_visit.pop_front();

        visited_status[&current_tri] = true;

        // add all the triangles incident to the edges of current tri
        for(int i=0; i<3; ++i)
        {
          Edge e;
          if(i==2) // since current_tri is sorted and we need to find e in the map
          {
            e[0] = current_tri[0];
            e[1] = current_tri[2];
          }
          else
          {
            e[0] = current_tri[i];
            e[1] = current_tri[i+1];
          }

          typename std::deque<const Tri*>::iterator dit = incident_tris[e].begin();
          for(; dit!=incident_tris[e].begin(); ++dit)
          {
            if(!visited_status[*dit])
              triangles_to_visit.push_back(*dit);
          }
        }
      }

      typename std::map<const Tri*, bool>::const_iterator mit = visited_status.begin();
      for(; mit != visited_status.end(); ++mit)
      {
        if(!mit->second)
        {
          std::cout << "Blob not achieved" << std::endl;
          failed_blob_counter++;
          break;
        }
      }
    }
    std::cout << "end edelsbrunner" << std::endl;
    std::cout << failed_blob_counter << " out of " << dual_map.size()
              << " failed" << std::endl;
  }

  void compute_dual() const
  {
#if (verbose > 5)
    std::cout << "dual computations" << std::endl;
#endif

    for(std::size_t i=0; i<this->tetrahedra.size(); ++i)
    {
      const Tet& tet = this->tetrahedra[i];
      Simplex dual_simplex;

      for(std::size_t j=0; j<Tet::size(); ++j)
        dual_simplex.insert(this->points[tet[j]].closest_seed_id);

#if (verbose>20)
      std::cout << "tet: " << std::endl;
      for(std::size_t j=0; j<tet.size(); ++j)
        std::cout << tet[j] << " (" << this->points[tet[j]].point << ")" << std::endl;
      std::cout << " witnesses : ";

      typename Simplex::const_iterator sit = dual_simplex.begin();
      typename Simplex::const_iterator send = dual_simplex.end();
      for(; sit!=send; ++sit)
        std::cout << *sit << " ";
      std::cout << std::endl;
#endif

      Base::add_simplex_to_triangulation(tet, dual_simplex);
    }
  }

  void output_canvas(const std::string str_base) const
  {
    std::ofstream out((str_base + ".mesh").c_str());
    std::ofstream out_bb((str_base + ".bb").c_str());

    out << "MeshVersionFormatted 1" << std::endl;
    out << "Dimension 3" << std::endl;
    out << "Vertices" << std::endl;
    out << this->points.size() << std::endl;
    out_bb << "3 1 " << this->points.size() << " 2" << std::endl;

    for(std::size_t i=0; i<this->points.size(); ++i)
    {
      const Canvas_point& cp = this->points[i];
      out << cp.point.x() << " " << cp.point.y() << " " << cp.point.z() << " " << i+1 << std::endl;

//      out_bb << cp.closest_seed_id << std::endl;
      out_bb << cp.distance_to_closest_seed << std::endl;
    }

    out << "Tetrahedra" << std::endl;
    out << this->tetrahedra.size() << std::endl;
    for(std::size_t i=0; i<this->tetrahedra.size(); ++i)
    {
      std::set<std::size_t> materials;
      for(std::size_t j=0; j<Tet::size(); ++j)
      {
        materials.insert(this->points[this->tetrahedra[i][j]].closest_seed_id);
        out << this->tetrahedra[i][j]+1 << " ";
      }
      std::size_t mat = (materials.size()==1)?(*(materials.begin())):(this->seeds.size());
      out << mat << std::endl;
    }

    out_bb << "End" << std::endl;
    out << "End" << std::endl;
  }

  Campen_canvas(const std::string& _canvas_str,
                const std::string& _seeds_str,
                const std::size_t _max_seeds_n,
                const std::size_t _n_refine,
                const Metric_field _mf)
    :
      Base(_canvas_str, _seeds_str, _max_seeds_n, _n_refine, _mf)
  { }
};

} // namespace Anisotropic_mesh_3
} // namespace CGAL

using namespace CGAL::Anisotropic_mesh_3;

int main(int, char**)
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
  typedef CGAL::Exact_predicates_exact_constructions_kernel    KExact;

  typedef Euclidean_metric_field<K>*                           MF;
//  typedef Custom_metric_field<K>*                              MF;

  typedef Campen_canvas_point<K>                               Campen_canvas_point;
  typedef Campen_canvas<K, KExact, Campen_canvas_point, MF>    Canvas;

  std::cout.precision(17);
//  std::freopen("log.txt", "w", stdout);

  double duration;
  std::clock_t start = std::clock();
  std::srand(0);

  MF mf = new Euclidean_metric_field<K>(1., 1., 3.);
//  MF mf = new Custom_metric_field<K>();

  const std::string canvas_str = "mild_base_mesh";
  const std::string seeds_str = "base_mesh.mesh";
  std::size_t max_seeds_n = 1;
  std::size_t n_refine = 2;

//  CGAL::generate_canvas<K>();
//  exit(0);

  Canvas canvas(canvas_str, seeds_str, max_seeds_n, n_refine, mf);

  canvas.initialize();
  canvas.paint();
  canvas.refine();
  canvas.output_canvas_data_and_dual(canvas_str + "_tr");

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cerr << "duration: " << duration << std::endl;
}

#ifndef CGAL_ANISOTROPIC_MESH_3_MVPQ_HELPER_H
#define CGAL_ANISOTROPIC_MESH_3_MVPQ_HELPER_H

//mutable vertex priority queue :
//vertices with a null metric are sorted according to the number of vertices with a non null metric in the first ring.

#include <Eigen/Dense>

#include <CGAL/assertions.h>
#include <CGAL/Modifiable_priority_queue.h>
#include <CGAL/squared_distance_3.h>

#include <CGAL/helpers/metric_helper.h>

#include <iostream>
#include <functional>
#include <set>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template <typename K>
typename K::FT angle_in_radian(const typename K::Vector_3& u,
                               const typename K::Vector_3& v,
                               K k = K())
{
  typedef typename K::FT       FT;
  typedef typename K::Vector_3 Vector_3;

  typename K::Compute_scalar_product_3 scalar_product =
    k.compute_scalar_product_3_object();
  typename K::Construct_cross_product_vector_3 cross_product =
    k.construct_cross_product_vector_3_object();
  typename K::Compute_squared_length_3 sq_length =
    k.compute_squared_length_3_object();

  // -------------------------------------
  // Angle between two vectors (in rad)
  // uv = |u||v| cos(u,v)
  // u^v  = w
  // |w| = |u||v| |sin(u,v)|
  // -------------------------------------
  FT product = CGAL::sqrt(sq_length(u) * sq_length(v));

  // Check
  if ( product == FT(0) )
    return FT(0);

  // Sine
  Vector_3 w = cross_product(u,v);
  FT abs_sin = CGAL::sqrt(sq_length(w)) / product;

  if ( abs_sin < FT(-1) ) { abs_sin = FT(-1); }
  if ( abs_sin > FT(1) ) { abs_sin = FT(1); }

  // We just need cosine sign
  FT cosine_sign = scalar_product(u,v);

  if ( cosine_sign >= FT(0) )
    return FT(std::asin(abs_sin));
  else
    return FT(CGAL_PI) - FT(std::asin(abs_sin));
}


template<typename Colored_polyhedron>
struct More_v
{
  typedef std::pair<typename Colored_polyhedron::Vertex_handle, int> mvpq_type;

  bool operator()(const mvpq_type* t1,const mvpq_type* t2) const
  {
    return t1->second > t2->second;
  }
};

//property map
template<typename Colored_polyhedron>
struct First_of_pair_v
{
  typedef std::pair<typename Colored_polyhedron::Vertex_handle, int> mvpq_type;
  typedef mvpq_type*                                                 key_type;
  typedef int                                                        value_type;
  typedef boost::readable_property_map_tag                           category;
};

//get() function for property map.
template<typename Colored_polyhedron>
typename First_of_pair_v<Colored_polyhedron>::value_type
get(const First_of_pair_v<Colored_polyhedron>&,
    const typename First_of_pair_v<Colored_polyhedron>::key_type& k)
{
  return k->first->tag();
}

template<typename Colored_polyhedron>
class Colored_modifiable_vertex_priority_queue :
  public Modifiable_priority_queue< std::pair<typename Colored_polyhedron::Vertex_handle, int>*,
                                    More_v<Colored_polyhedron>,
                                    First_of_pair_v<Colored_polyhedron> >
  {
  public:
    typedef Modifiable_priority_queue<
             std::pair<typename Colored_polyhedron::Vertex_handle,
                       int>*,
             More_v<Colored_polyhedron>,
             First_of_pair_v<Colored_polyhedron> >                         Base;
    typedef typename Base::IndexedType                                     IndexedType;
    typedef typename Base::Compare                                         Compare;
    typedef typename Base::ID                                              ID;
    typedef typename Base::size_type                                       size_type;

    typedef typename std::pair<typename Colored_polyhedron::Vertex_handle, int> Cmvpq_type;

    typedef typename Colored_polyhedron::Traits::Kernel                    K;
    typedef typename K::FT                                                 FT;
    typedef typename K::Vector_3                                           Vector_3;
    typedef typename Colored_polyhedron::Point_3                           Point_3;
    typedef typename Colored_polyhedron::Vertex_handle                     Vertex_handle;
    typedef typename Colored_polyhedron::Facet_handle                      Facet_handle;
    typedef typename Colored_polyhedron::Vertex_iterator                   Vertex_iterator;
    typedef typename Colored_polyhedron::Halfedge_around_vertex_circulator HV_circulator;

    struct Vertex_Compare
    {
    public:
      bool operator()(const Vertex_handle& v1, const Vertex_handle& v2) const
      {
        return v1 < v2;
      }

      Vertex_Compare(){}
    };

    FT compute_coeff_at_vertex(const Vertex_handle& v0,
                               HV_circulator h)
    {
      //this can probably be easily optimized noting that pi+1's prev is pi
      //and other similar relations...
      //naive for now

      Point_3 p0 = v0->point();
      Point_3 pi = h->opposite()->vertex()->point();
      FT p0pi = std::sqrt(CGAL::squared_distance(p0, pi));
      Point_3 prev, next;

      HV_circulator init = h;
      //std::cout << "first h " << h->opposite()->vertex()->point() << std::endl;
      //std::cout << "first init " << init->opposite()->vertex()->point() << std::endl;

      //find previous
      while(1)
      {
        h--;
        //std::cout << "now testing -- : " << h->opposite()->vertex()->point() << std::endl;
        if(h->opposite()->vertex()->is_colored())
          break;
      }

      prev = h->opposite()->vertex()->point();
      if(p0 == prev) //only one point is colored
      {
        std::cout << "1 colored"<< std::endl;
        return 1;
      }

      //find next
      h = init;
      //std::cout << "second h " << h->opposite()->vertex()->point() << std::endl;
      //std::cout << "second init " << init->opposite()->vertex()->point() << std::endl;

      while(1)
      {
        h++;
        //std::cout << "now testing ++ : " << h->opposite()->vertex()->point() << std::endl;
        if(h->opposite()->vertex()->is_colored())
          break;
      }

      next = h->opposite()->vertex()->point();
      if(prev == next) //only two points are colored
      {
        std::cout << "2 colored : " << 1./p0pi << std::endl;
        return 1./p0pi;
      }

      //general case
      FT alpha_i = 0.5 * angle_in_radian<K>(Vector_3(p0,pi),Vector_3(p0,prev));
      FT alpha_j = 0.5 * angle_in_radian<K>(Vector_3(p0,next),Vector_3(p0,pi));

      //coeff
      FT coeff = (std::tan(alpha_i) + std::tan(alpha_j)) / p0pi;

      /*
      std::cout << "compute coeff at point : " << p0 << std::endl;
      std::cout << "opposite point is ---- : " << pi << std::endl;
      std::cout << "prev and next : " << prev << " " << next << std::endl;
      std::cout << "angles : "<< alpha_i << " " << alpha_j << std::endl;
      std::cout << "coeff in fine : "<< coeff << std::endl;
      std::cout << "check h at the end : " << pi << " should be equal to " << h->opposite()->vertex()->point() << std::endl;
      */

      //re-init h because it's needed outside of this function
      h = init;

      return coeff;
    }

    void get_first_ring(const Vertex_handle& v,
                        std::set<Vertex_handle, Vertex_Compare>& vs,
                        std::map<Vertex_handle, FT>& coeffs,
                        bool compute_interpolation_coeffs)
    {
      FT sum = 0;

      HV_circulator h = v->vertex_begin();
      HV_circulator hend = h;
      do
      {
        Vertex_handle vi = h->opposite()->vertex();
        vs.insert(vi);

        if(compute_interpolation_coeffs)
          if(vi->is_colored())
          {
            FT cvi = compute_coeff_at_vertex(v, h);
            coeffs[vi] = cvi;
            sum += cvi;
          }
        h++;
      }
      while(h != hend);

      if(compute_interpolation_coeffs)
      {
        typename std::map<Vertex_handle, FT>::iterator it = coeffs.begin();
        for(; it != coeffs.end(); ++it)
          it->second /= sum;

        FT sum_check = 0;
        it = coeffs.begin();
        for(; it != coeffs.end(); ++it)
           sum_check += it->second;

        //std::cout << "sum check should be 1 : " << sum_check << std::endl;
      }
    }

    void get_first_ring(const Vertex_handle& v,
                        std::set<Vertex_handle, Vertex_Compare>& vs)
    {
      std::map<Vertex_handle, FT> useless;
      get_first_ring(v, vs, useless, false);
    }

    int count_colored_vertices_in_first_ring(const Vertex_handle& v)
    {
      int ret = 0;
      Vertex_Compare vc;
      std::set<Vertex_handle, Vertex_Compare> neigh_vertices(vc);

      get_first_ring(v, neigh_vertices);

      typename std::set<Vertex_handle, Vertex_Compare>::iterator it = neigh_vertices.begin();
      typename std::set<Vertex_handle, Vertex_Compare>::iterator itend = neigh_vertices.end();
      for(; it!=itend; ++it)
        if((*it)->is_colored())
          ret++;

      return ret;
    }

    void initialize_cmvpq(Colored_polyhedron& P)
    {
      std::cout << "initalizing cmvpq..." << std::endl;
      data.reserve(P.size_of_vertices());

      for(Vertex_iterator v = P.vertices_begin(); v != P.vertices_end(); ++v)
      {
        if(v->is_colored()) //queue is for uncolored vertices only
          continue;

        int n = count_colored_vertices_in_first_ring(v);
        data[v->tag()] = Cmvpq_type(v, n);
        this->push(&data[0]+v->tag());
      }
    }

    void update_vertex(const Vertex_handle& v, int value)
    {
      int index = v->tag();
      data[index].second = value;
      this->update(&data[0]+index, true);
    }

    void increase_vertex_value(const Vertex_handle& v)
    {
      int index = v->tag();
      data[index].second++;
      this->update(&data[0]+index, true);
    }

    void increase_neigh_vertices(std::set<Vertex_handle, Vertex_Compare>& neigh_vertices)
    {
      typename std::set<Vertex_handle, Vertex_Compare>::iterator it = neigh_vertices.begin();
      typename std::set<Vertex_handle, Vertex_Compare>::iterator itend = neigh_vertices.end();
      for(; it!=itend; ++it)
        if(!(*it)->is_colored()) //only increase value of non colored vertices
          increase_vertex_value(*it);
    }

    Eigen::Matrix3d compute_top_vertex_metric_blend(Vertex_handle& v)
    {
      std::set<Vertex_handle, Vertex_Compare> neigh_vertices;
      get_first_ring(v, neigh_vertices);

      Point_3 p = v->point();

      typename std::set<Vertex_handle, Vertex_Compare>::iterator it = neigh_vertices.begin();
      typename std::set<Vertex_handle, Vertex_Compare>::iterator itend = neigh_vertices.end();

      //std::cout << "*--*-*-**-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << std::endl;
      //std::cout << "coloring : " << top_vertex->first->tag() << " " << top_vertex->first->point() << std::endl;

      FT sum_dist = 0.;
      int counter = 0;
      std::vector<FT> sq_distances;
      for(; it!=itend; ++it)
      {
        if((*it)->is_colored())
        {
          //std::cout << "in : " << (*it)->tag() << " ranked : " << (*it)->colored_rank() << std::endl;
          sq_distances.push_back(CGAL::squared_distance(p, (*it)->point()));
          sum_dist += std::sqrt(sq_distances[counter]);
          counter++;
        }
      }

      FT rb = sum_dist / (counter);
      FT wpp_inv = 1./(rb * rb);

      // blend
      Eigen::Matrix3d m = Eigen::Matrix3d::Zero();
      FT wsum = 0.;
      it = neigh_vertices.begin();
      counter = 0;
      for(; it!=itend; ++it)
      {
        if((*it)->is_colored())
        {
          FT w = std::exp(-wpp_inv * sq_distances[counter]);
          wsum += w;
          for(int j = 0; j < 3; j++)
            for(int k = 0; k < 3; k++)
              m(j,k) = m(j,k) + w * ((*it)->metric())(j,k);
          counter++;
        }
      }

      FT wsum_inv = 1. / wsum;
      for(int j = 0; j < 3; j++)
        for(int k = 0; k < 3; k++)
          m(j,k) = m(j,k) * wsum_inv;

      increase_neigh_vertices(neigh_vertices);
      return m;
    }

    Eigen::Matrix3d compute_top_vertex_metric_interpolation(Vertex_handle& v)
    {
      std::map<Vertex_handle, FT> coeffs;
      std::set<Vertex_handle, Vertex_Compare> neigh_vertices;
      get_first_ring(v, neigh_vertices, coeffs, true /*compute coeffs*/);

      if(coeffs.empty())
      {
        std::cout << "surely you're joking mr Feynman" << std::endl;
        return Eigen::Matrix3d::Zero();
      }

      //std::cout << "check map size " << count_colored_vertices_in_first_ring(v) << " " << coeffs.size() << std::endl;

      std::vector<std::pair<Eigen::Matrix3d, typename K::FT> > w_metrics;

      typename std::set<Vertex_handle, Vertex_Compare>::iterator it = neigh_vertices.begin();
      typename std::set<Vertex_handle, Vertex_Compare>::iterator itend = neigh_vertices.end();
      for(; it!=itend; ++it)
        if((*it)->is_colored())
          w_metrics.push_back(std::make_pair((*it)->metric(), coeffs[*it]));

      Eigen::Matrix3d m = logexp_interpolate<K>(w_metrics);
      increase_neigh_vertices(neigh_vertices);
      return m;
    }

    Eigen::Matrix3d compute_top_vertex_metric_intersection(Vertex_handle& v)
    {
      std::set<Vertex_handle, Vertex_Compare> neigh_vertices;
      get_first_ring(v, neigh_vertices);

      Point_3 p = v->point();
      Eigen::Matrix3d m = Eigen::Matrix3d::Zero();

      int counter = 0;
      typename std::set<Vertex_handle, Vertex_Compare>::iterator it = neigh_vertices.begin();
      typename std::set<Vertex_handle, Vertex_Compare>::iterator itend = neigh_vertices.end();

      std::cout << "*--*-*-**-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << std::endl;
      std::cout << "coloring : " << v->tag() << " " << v->point() << std::endl;

      for(; it!=itend; ++it)
      {
        if((*it)->is_colored())
        {
          std::cout << "in : " << (*it)->tag() << " ranked : " << (*it)->colored_rank() << std::endl;
          counter++;
          Eigen::Matrix3d scaled_m = (*it)->metric(); //scale_matrix_to_point<K>((*it)->metric(),(*it)->point(), p);
          std::cout << "scaled is : " << std::endl << scaled_m << std::endl;
          if(counter == 1)
            m = scaled_m;
          else
            m = matrix_intersection<K>(m, scaled_m);
          std::cout << "count : " << counter << std::endl << m << std::endl;
        }
      }
      increase_neigh_vertices(neigh_vertices);
      return m;
    }

    void color_top_vertex(int rank)
    {
      Cmvpq_type* top_vertex = this->extract_top().get();

      if(top_vertex->second == 0)
      {
        std::cout << "top vertex has zero neighbors colored" << std::endl;
        return;
      }

      //Eigen::Matrix3d new_vertex_metric = compute_top_vertex_metric_intersection(top_vertex->first);
      //Eigen::Matrix3d new_vertex_metric = compute_top_vertex_metric_blend(top_vertex->first);
      Eigen::Matrix3d new_vertex_metric = compute_top_vertex_metric_interpolation(top_vertex->first);

      top_vertex->first->metric() = new_vertex_metric;
      top_vertex->first->metric_origin() = 3;
      top_vertex->first->colored_rank() = rank;
    }

    void color_all_vertices()
    {
      std::cout << "spreading colors..." << std::endl;
      int rank = 0;
      while(!this->empty())
        color_top_vertex(rank++);
    }

    Colored_modifiable_vertex_priority_queue(size_type largest_ID,
                                      Compare const& c,
                                      ID const& id )
      : Base(largest_ID, c, id)
    {}

  private:
    std::vector<Cmvpq_type> data;
  };

} //namespace Aniso
} //namespace CGAL

#endif


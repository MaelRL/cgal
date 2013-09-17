#ifndef mvpq_HELPER_H
#define mvpq_HELPER_H

//mutable vertex priority queue :
//vertices with a null metric are sorted according to the number of vertices with a non null metric in the first ring.

#include <Eigen/Dense>

#include <CGAL/assertions.h>
#include <CGAL/Modifiable_priority_queue.h>

#include <CGAL/helpers/metric_helper.h>

#include <iostream>
#include <functional>
#include <set>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

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

    void get_first_ring(const Vertex_handle& v,
                        std::set<Vertex_handle, Vertex_Compare>& vs)
    {
      HV_circulator h = v->vertex_begin();
      HV_circulator hend = h;
      do
      {
        vs.insert(h->opposite()->vertex());
        h++;
      }
      while(h != hend);
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
      update(&data[0]+index, true);
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

    void color_top_vertex()
    {
      Cmvpq_type* top_vertex = this->extract_top().get();

      if(top_vertex->second == 0)
      {
        std::cout << "top vertex has zero neighbors colored" << std::endl;
        return;
      }

      std::set<Vertex_handle, Vertex_Compare> neigh_vertices;
      get_first_ring(top_vertex->first, neigh_vertices);

      Point_3 p = top_vertex->first->point();
      Eigen::Matrix3d new_vertex_metric = Eigen::Matrix3d::Zero();

      int counter = 0;
      typename std::set<Vertex_handle, Vertex_Compare>::iterator it = neigh_vertices.begin();
      typename std::set<Vertex_handle, Vertex_Compare>::iterator itend = neigh_vertices.end();

      for(; it!=itend; ++it)
      {
        if((*it)->is_colored())
        {
          counter++;
          Eigen::Matrix3d scaled_m = scale_matrix_to_point<K>((*it)->metric(), p, (*it)->point());
          if(counter == 1)
            new_vertex_metric = scaled_m;
          else
            new_vertex_metric = matrix_intersection<K>(new_vertex_metric, scaled_m);
        }
      }

//      std::cout << counter << " " << top_vertex->second << std::endl;
      top_vertex->first->metric() = new_vertex_metric;

      increase_neigh_vertices(neigh_vertices);
    }

    void color_all_vertices()
    {
      std::cout << "spreading colors..." << std::endl;
      while(!this->empty())
        color_top_vertex();
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


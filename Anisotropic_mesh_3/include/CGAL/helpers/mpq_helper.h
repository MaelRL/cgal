#ifndef CGAL_ANISOTROPIC_MESH_3_MPQ_HELPER_H
#define CGAL_ANISOTROPIC_MESH_3_MPQ_HELPER_H

#include <CGAL/assertions.h>
#include <CGAL/Modifiable_priority_queue.h>

#include <iostream>
#include <functional>
#include <set>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename Colored_polyhedron>
struct More
{
  typedef std::pair<typename Colored_polyhedron::Facet_handle, int> Mpq_type;

  bool operator()(const Mpq_type* t1,const Mpq_type* t2) const
  {
    return t1->second > t2->second;
  }
};

//property map
template<typename Colored_polyhedron>
struct First_of_pair
{
  typedef std::pair<typename Colored_polyhedron::Facet_handle, int> Mpq_type;
  typedef Mpq_type*                                                 key_type;
  typedef int                                                       value_type;
  typedef boost::readable_property_map_tag                          category;
};

//get() function for property map.
template<typename Colored_polyhedron>
typename First_of_pair<Colored_polyhedron>::value_type
get(const First_of_pair<Colored_polyhedron>&,
    const typename First_of_pair<Colored_polyhedron>::key_type& k)
{
  return k->first->tag();
}

template<typename Colored_polyhedron>
class Colored_modifiable_priority_queue :
  public Modifiable_priority_queue< std::pair<typename Colored_polyhedron::Facet_handle, int>*,
                                    More<Colored_polyhedron>,
                                    First_of_pair<Colored_polyhedron> >
  {
  public:
    typedef Modifiable_priority_queue<
             std::pair<typename Colored_polyhedron::Facet_handle,
                       int>*,
             More<Colored_polyhedron>,
             First_of_pair<Colored_polyhedron> >                           Base;
    typedef typename Base::IndexedType                                     IndexedType;
    typedef typename Base::Compare                                         Compare;
    typedef typename Base::ID                                              ID;
    typedef typename Base::size_type                                       size_type;

    typedef typename std::pair<typename Colored_polyhedron::Facet_handle, int> Cmpq_type;

    typedef typename Colored_polyhedron::Vertex_handle                     Vertex_handle;
    typedef typename Colored_polyhedron::Facet_handle                      Facet_handle;
    typedef typename Colored_polyhedron::Facet_iterator                    Facet_iterator;
    typedef typename Colored_polyhedron::Halfedge_around_vertex_circulator HV_circulator;

    struct Facet_Compare
    {
    public:
      bool operator()(const Facet_handle& fh1, const Facet_handle& fh2) const
      {
        return fh1 < fh2;
      }

      Facet_Compare(){}
    };

    void get_first_ring(const Facet_handle& f,
                        std::set<Facet_handle, Facet_Compare>& fhs)
    {
      std::vector<Vertex_handle> v(3);
      v[0] = f->halfedge()->vertex();
      v[1] = f->halfedge()->next()->vertex();
      v[2] = f->halfedge()->next()->next()->vertex();

      for(int i=0; i<3; ++i)
      {
        HV_circulator h = v[i]->vertex_begin();
        HV_circulator hend = h;
        do
        {
          fhs.insert(h->facet());
          h++;
        }
        while(h != hend);
      }
    }

    int count_colored_facets_in_first_ring(const Facet_handle& f)
    {
      int ret = 0;
      Facet_Compare fc;
      std::set<Facet_handle, Facet_Compare> neigh_facets(fc);

      get_first_ring(f, neigh_facets);

      typename std::set<Facet_handle, Facet_Compare>::iterator it = neigh_facets.begin();
      typename std::set<Facet_handle, Facet_Compare>::iterator itend = neigh_facets.end();
      for(; it!=itend; ++it)
        if((*it)->color() > 0.)
          ret++;

      return ret;
    }

    void initialize_cmpq(Colored_polyhedron& P)
    {
      std::cout << "initalizing cmpq..." << std::endl;
      data.reserve(P.size_of_facets());

      for(Facet_iterator f = P.facets_begin(); f != P.facets_end(); ++f)
      {
        if(f->color() > 0.) //queue is for uncolored facets only
          continue;

        int n = count_colored_facets_in_first_ring(f);
        data[f->tag()] = Cmpq_type(f, n);
        this->push(&data[0]+f->tag());
      }
    }

    void update_facet(const Facet_handle& f, int value)
    {
      int index = f->tag();
      data[index].second = value;
      update(&data[0]+index, true);
    }

    void increase_facet_value(const Facet_handle& f)
    {
      int index = f->tag();
      data[index].second++;
      this->update(&data[0]+index, true);
    }

    void increase_neigh_facets(std::set<Facet_handle, Facet_Compare>& neigh_facets)
    {
      typename std::set<Facet_handle, Facet_Compare>::iterator it = neigh_facets.begin();
      typename std::set<Facet_handle, Facet_Compare>::iterator itend = neigh_facets.end();
      for(; it!=itend; ++it)
        if((*it)->color() == 0.) //only increase value of non colored facets
          increase_facet_value(*it);
    }

    void color_top_facet()
    {
      Cmpq_type* top_facet = this->extract_top().get();

      if(top_facet->second == 0)
      {
        std::cout << "top facet has zero neighbors colored" << std::endl;
        return;
      }
      std::set<Facet_handle, Facet_Compare> neigh_facets;
      get_first_ring(top_facet->first, neigh_facets);

      double new_facet_color = 0.;
      int counter = 0;
      typename std::set<Facet_handle, Facet_Compare>::iterator it = neigh_facets.begin();
      typename std::set<Facet_handle, Facet_Compare>::iterator itend = neigh_facets.end();
      for(; it!=itend; ++it)
      {
        if((*it)->color() > 0.)
        {
          new_facet_color += (*it)->color();
          counter++;
        }
      }
      new_facet_color /= top_facet->second;
      top_facet->first->color() = new_facet_color;
      increase_neigh_facets(neigh_facets);
    }

    void color_all_facets()
    {
      std::cout << "spreading colors..." << std::endl;
      while(!this->empty())
        color_top_facet();
    }

    Colored_modifiable_priority_queue(size_type largest_ID,
                                      Compare const& c,
                                      ID const& id )
      : Base(largest_ID, c, id)
    {}

  private:
    std::vector<Cmpq_type> data;
  };

} //namespace Aniso
} //namespace CGAL

#endif


#ifndef CGAL_ANISOTROPIC_MESH_3_FACET_REFINE_QUEUE_H
#define CGAL_ANISOTROPIC_MESH_3_FACET_REFINE_QUEUE_H

#include <CGAL/Stretched_Delaunay_3.h>
#include <CGAL/helpers/combinatorics_helper.h>

#ifndef ANISO_USE_BOOST_UNORDERED_SET
#define ANISO_USE_BOOST_UNORDERED_SET
#include <boost/functional/hash.hpp>
#include <boost/unordered_set.hpp>
#endif

#include <set>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K, typename KExact = K>
class Refine_facet_comparer;

template<typename K, typename KExact = K>
class Refine_facet_iterator_comparer;

template<typename K, typename KExact = K>
class Refine_facet_hash;

template<typename K, typename KExact = K>
class Refine_facet
{
  typedef Refine_facet<K, KExact>                           Self;

public:
  typedef typename K::FT                                    FT;
  typedef Stretched_Delaunay_3<K, KExact>                   Star;
  typedef Star*                                             Star_handle;
  typedef typename Star::Facet                              Facet;

  typedef Refine_facet_comparer<K, KExact>                  Rfacet_comparer;
  typedef Refine_facet_iterator_comparer<K, KExact>         Rfacet_it_comparer;
  typedef Refine_facet_hash<K, KExact>                      Rfacet_hash;

#ifdef ANISO_USE_BOOST_UNORDERED_SET
  typedef typename boost::unordered_set<Self,
                                        Rfacet_hash,
                                        Rfacet_comparer>    Rfacet_set;
#else
  typedef typename std::set<Self, Rfacet_comparer>          Rfacet_set;
#endif

  typedef typename Rfacet_set::iterator                     Rfacet_set_iterator;
  typedef typename std::multiset<Rfacet_set_iterator,
                                 Rfacet_it_comparer>        Queue;
  typedef typename Queue::iterator                          Queue_iterator;

public:
  Star_handle star;
  Facet_ijk facet;

  FT value;
  int queue_type; // type of queue
  mutable Queue_iterator queue_it; // position in queues[queue_type]
  bool prev_rejection;

public:
  Refine_facet() : star(NULL), value(0), queue_type(-1), prev_rejection(false) { }
  Refine_facet(Star_handle star_, const Facet& facet_, FT value_, int queue_type_)
    :
    star(star_),
    facet(facet_),
    value(value_),
    queue_type(queue_type_),
    prev_rejection(false)
  { }
};

template<typename K, typename KExact>
inline std::ostream& operator<<(std::ostream& os, const Refine_facet<K, KExact>& src)
{
  os << src.star->index_in_star_set() << " || ";
  os << src.facet[0] << " " << src.facet[1] << " " << src.facet[2];
  os << " || qt: " << src.queue_type;
  os << " || val: " << src.value;
  os << " || rej: " << src.prev_rejection << std::endl;
  return os;
}

template<typename K, typename KExact/* = K */>
class Refine_facet_comparer
{
public:
  Refine_facet_comparer() { }
  bool operator() (const Refine_facet<K, KExact>& left,
                   const Refine_facet<K, KExact>& right) const
  {
#ifdef ANISO_USE_BOOST_UNORDERED_SET
    return std::equal(left.facet.vertices().begin(),
                      left.facet.vertices().end(),
                      right.facet.vertices().begin());
#else
    return std::lexicographical_compare(left.facet.vertices().begin(),
                                        left.facet.vertices().end(),
                                        right.facet.vertices().begin(),
                                        right.facet.vertices().end());
#endif
  }
};

template<typename K, typename KExact>
class Refine_facet_hash
{
public:
  typedef Refine_facet<K, KExact>                           Rfacet;

  Refine_facet_hash(){}

  std::size_t operator()(const Rfacet& rf) const
  {
    return boost::hash_range(rf.facet.vertices().begin(),
                             rf.facet.vertices().end());
  }
};

template<typename K, typename KExact/* = K*/>
class Refine_facet_iterator_comparer
{
  typedef Refine_facet<K, KExact>                           Rfacet;
  typedef Refine_facet_comparer<K, KExact>                  Rfacet_comparer;

#ifdef ANISO_USE_BOOST_UNORDERED_SET
  typedef typename boost::unordered_set<Rfacet,
                                        Refine_facet_hash<K, KExact>,
                                        Rfacet_comparer>    Rfacet_set;
#else
  typedef typename std::set<Rfacet, Rfacet_comparer>        Rfacet_set;
#endif

  typedef typename Rfacet_set::iterator                     Rfacet_set_iterator;
  typedef typename Rfacet_set::const_iterator               Rfacet_set_c_iterator;

public:
  Refine_facet_iterator_comparer() { }

  bool operator() (Rfacet_set_iterator left,
                   Rfacet_set_iterator right) const
  {
    return left->value > right->value;
  }
};

template<typename K, typename KExact = K>
class Facet_refine_queue
{
public:
  typedef typename K::FT                                    FT;
  typedef typename K::Point_3                               Point_3;

  typedef Refine_facet<K, KExact>                           Rfacet;
  typedef Refine_facet_comparer<K, KExact>                  Rfacet_comparer;
  typedef Refine_facet_iterator_comparer<K, KExact>         Rfacet_it_comparer;

#ifdef ANISO_USE_BOOST_UNORDERED_SET
  typedef typename boost::unordered_set<Rfacet,
                                        Refine_facet_hash<K, KExact>,
                                        Rfacet_comparer>    Rfacet_set;
#else
  typedef typename std::set<Rfacet, Rfacet_comparer>        Rfacet_set;
#endif

  typedef typename Rfacet_set::iterator                     Rfacet_set_iterator;

  typedef typename std::multiset<Rfacet_set_iterator,
                                 Rfacet_it_comparer>        Queue;
  typedef typename Queue::iterator                          Queue_iterator;

  typedef typename Rfacet::Star_handle                      Star_handle;
  typedef typename Rfacet::Facet                            Facet;

public:
  static const int nb_queues = 6;
  static const int encroachment_queue = 0;
  static const int over_distortion_queue = 1;
  static const int over_circumradius_queue = 2;
  static const int over_approximation_queue = 3;
  static const int bad_shape_queue = 4;
  static const int start_pick_valid = 5;
  static const int inconsistent_queue = 5;

  //rfacets: set of the facets needing refinement. Avoids having the same facet (i,j,k)
  //         in different queues (either due to diff values, star or queue_type).

  //queues: multisets of iterators to the rfacets set, sorted by value. Multiset chosen
  //        over prio_queue to access/remove any element.

private:
  Rfacet_set rfacets;
  Queue *queues[nb_queues];

  Queue encroachments;
  Queue over_distortions;
  Queue over_circumradii;
  Queue over_approximation;
  Queue bad_shapes;
  Queue inconsistents;

private:
  void update_rfacet(Rfacet_set_iterator& rfacet_it,
                     Star_handle star,
                     FT value,
                     int queue_type,
                     bool prev_rejection)
  {
/*
    std::cout << "update rfacet " << *rfacet_it << std::endl;
    std::cout << "with new values" << std::endl;
    std::cout << star->index_in_star_set() << " " << value << " " << queue_type << std::endl;
*/

    Rfacet_set_iterator rfacet_hint = rfacets.end();
    Queue_iterator queue_hint = queues[queue_type]->end();

#ifndef ANISO_USE_BOOST_UNORDERED_SET
    if(rfacet_it == rfacets.begin())
      rfacet_hint = rfacets.end(); //can't use hint in that case, end() is safe
    else
      rfacet_hint = (--rfacet_it)++;
#endif

    //In the general case, we don't know where the insertion will be. In the rejection case
    //the insertion is done at the end.
    if(!prev_rejection || queues[queue_type]->size() <= 1)
      queue_hint = queues[queue_type]->end();
    else
      queue_hint = --(queues[queue_type]->end());

    Rfacet new_rfacet = *rfacet_it;
    new_rfacet.star = star;
    new_rfacet.value = value;
    new_rfacet.queue_type = queue_type;
    new_rfacet.prev_rejection = prev_rejection;

    //clean both queues of that element
    queues[rfacet_it->queue_type]->erase(rfacet_it->queue_it);

#ifdef ANISO_USE_BOOST_UNORDERED_SET
    rfacet_hint =
#endif
    rfacets.erase(rfacet_it);

    //re-insert with the new values
    rfacet_it = rfacets.insert(rfacet_hint, new_rfacet);
    Queue_iterator qit = queues[queue_type]->insert(queue_hint, rfacet_it);

    rfacet_it->queue_it = qit;
  }

  void update_rfacet_value(Rfacet_set_iterator& rfacet_it,
                           FT value)
  {
    return update_rfacet(rfacet_it, rfacet_it->star, value,
                         rfacet_it->queue_type, rfacet_it->prev_rejection);
  }

public:
  //reject the top of the queue n° queue_type to the end of the same queue
  void reject_rfacet(int queue_type)
  {
    Rfacet_set_iterator rfacet_it = *(queues[queue_type]->begin());
    FT min_value = queue_min_value(queue_type);
    FT epsilon = 1e-10;
    return update_rfacet(rfacet_it, rfacet_it->star, min_value - epsilon,
                         rfacet_it->queue_type, true /*rejection*/);
  }

  void push(Star_handle star,
            const Facet& facet,
            FT value,
            int queue_type,
            bool force_push = false) //update even if it means a lower priority
  {
    Rfacet new_rfacet(star, facet, value, queue_type);

    //std::cout << "rfacet created ready to push : " << new_rfacet << std::endl;

    std::pair<Rfacet_set_iterator, bool> is_insert_successful;
    is_insert_successful = rfacets.insert(new_rfacet);

    if(!is_insert_successful.second) // the facet already exists in one of the queues
    {
      Rfacet_set_iterator old_rfacet_it = is_insert_successful.first;
      if(old_rfacet_it->queue_type > queue_type || // the new facet has a higher prio queue_type
         old_rfacet_it->value < value || // new value has higher prio within the same queue
         force_push)
      {
        /*
        std::cout << "refine_facet needs update: " << std::endl;
        std::cout << old_rfacet_it->queue_type << " " << queue_type << std::endl;
        std::cout << old_rfacet_it->value << " " << value << std::endl;
        */
        update_rfacet(old_rfacet_it, star, value, queue_type, old_rfacet_it->prev_rejection);
      }
    }
    else
    {
      //std::cout << "insert succesful" << std::endl;
      is_insert_successful.first->queue_it = queues[queue_type]->insert(is_insert_successful.first);
    }
  }

  bool need_picking_valid(int queue_type) const
  {
    return (queue_type >= start_pick_valid);
  }

  bool empty(int start_id = 0, int end_id = nb_queues - 1) const
  {
    for (int i=start_id; i<=end_id; i++)
      if (!queues[i]->empty())
        return false;
    return true;
  }

  bool top(Rfacet_set_iterator& facet,
           int start_id = 0, int end_id = nb_queues - 1) const
  {
    for (int i=start_id; i<=end_id; i++)
      if(!queues[i]->empty())
      {
        facet = *(queues[i]->begin());
/*
        std::cout << "found top at: " << i << std::endl;
        std::cout << "star: " << facet->star->index_in_star_set() << std::endl;
        std::cout << facet->facet[0] << " " << facet->facet[1] << " " << facet->facet[2] << std::endl;
*/
        return true;
      }
    return false;
  }

  bool pop(int start_id = 0, int end_id = nb_queues - 1)
  {
    for (int i=start_id; i<=end_id; i++)
      if (!queues[i]->empty())
      {
        rfacets.erase(*(queues[i]->begin()));
        queues[i]->erase(queues[i]->begin());
        return true;
      }
    return false;
  }

  void clear(int start_id = 0, int end_id = nb_queues - 1)
  {
    for(int i=start_id; i<=end_id; i++)
      queues[i]->clear();
  }

  FT queue_min_value(int queue_type)
  {
    return (*(--queues[queue_type]->end()))->value;
  }

  void clean()
  {
    CGAL::Timer t;
    t.start();
    int size = rfacets.size();

    Rfacet_set_iterator rfit = rfacets.begin();
    Rfacet_set_iterator rfend = rfacets.end();
    for(; rfit != rfend; ++rfit)
    {
      Facet f;
      if(!rfit->star->has_facet(rfit->facet, f))
      {
        queues[rfit->queue_type]->erase(rfit->queue_it);
        rfacets.erase(rfit);
      }
    }

    t.stop();
    std::cout << "clean: " << size << " " << rfacets.size() << " in " << t.time() << std::endl;
  }

public:
  void print_rfacets()
  {
    std::cout << "rfacets: " << std::endl;
    Rfacet_set_iterator it = rfacets.begin();
    Rfacet_set_iterator itend = rfacets.end();
    for(;it!=itend; ++it)
      std::cout << *it;
  }

  void print_queue(int queue_type)
  {
    std::cout << "fqueue : " << queue_type << std::endl;
    Queue* si = queues[queue_type];
    Queue_iterator vit = si->begin();
    Queue_iterator viend = si->end();
    for(; vit!=viend; ++vit)
      std::cout << **vit;
  }

  void print_queues(int start_id = 0, int end_id = nb_queues - 1)
  {
    for(int i=start_id; i<= end_id; ++i)
      print_queue(i);
  }

  void print() const
  {
    std::cout << "FQueue : ( ";
    for(int i = 0; i < nb_queues; i++)
      std::cout << queues[i]->size() <<" ";
    std::cout << ") ";
    std::cout << rfacets.size() << std::endl;
  }

  std::size_t count(int start_id = 0, int end_id = nb_queues - 1) const
  {
    std::size_t count = 0;
    for(int i=start_id; i<= end_id; ++i)
      count += queues[i]->size();
    return count;
  }

  bool is_facet_in(Star_handle star, Facet facet, FT value, int queue_type)
  {
    bool is_facet_in, is_rest_identical;
    Rfacet_set_iterator it = rfacets.find(Rfacet(star, facet, value, queue_type));

    if((is_facet_in = (it != rfacets.end()))) // = and not == on purpose
    {
      is_rest_identical = (star->index_in_star_set() == it->star->index_in_star_set() &&
                           queue_type == it->queue_type &&
                           std::abs(value - it->value) < 1e-11);
    }
/*
    if(is_facet_in && !is_rest_identical)
    {
      std::cout << "facet is in, but from a different star: " << std::endl;
      std::cout << "facet check: ";
      std::cout << facet.first->vertex((facet.second+1)%4)->info() << " ";
      std::cout << facet.first->vertex((facet.second+2)%4)->info() << " ";
      std::cout << facet.first->vertex((facet.second+3)%4)->info() << std::endl;
      std::cout << it->facet[0] << " " << it->facet[1] << " " << it->facet[2] << std::endl;
      std::cout << "stars: " << star->index_in_star_set() << " " << it->star->index_in_star_set() << std::endl;
      std::cout << "values: " << value << " " << it->value << std::endl;
      std::cout << "queue type: " << queue_type << " " << it->queue_type << std::endl;
    }
*/
    return is_facet_in;
  }

public:
  Facet_refine_queue()
    :
      rfacets(),
      encroachments(Rfacet_it_comparer()),
      over_distortions(Rfacet_it_comparer()),
      over_circumradii(Rfacet_it_comparer()),
      over_approximation(Rfacet_it_comparer()),
      bad_shapes(Rfacet_it_comparer()),
      inconsistents(Rfacet_it_comparer())
  {
    queues[encroachment_queue] = &encroachments;
    queues[over_distortion_queue] = &over_distortions;
    queues[over_circumradius_queue] = &over_circumradii;
    queues[over_approximation_queue] = &over_approximation;
    queues[bad_shape_queue] = &bad_shapes;
    queues[inconsistent_queue] = &inconsistents;
  }
}; // Facet_refine_queue

} // Anisotropic_mesh_3
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_REFINE_QUEUE_SURFACE_H

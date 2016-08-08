#ifndef CGAL_ANISOTROPIC_MESH_2_FACE_REFINE_QUEUE_H
#define CGAL_ANISOTROPIC_MESH_2_FACE_REFINE_QUEUE_H

#include <CGAL/Stretched_Delaunay_2.h>
#include <CGAL/helpers/combinatorics_helper.h>

#ifndef ANISO_USE_BOOST_UNORDERED_SET
#define ANISO_USE_BOOST_UNORDERED_SET
#include <boost/functional/hash.hpp>
#include <boost/unordered_set.hpp>
#endif

#include <set>

namespace CGAL
{
namespace Anisotropic_mesh_2
{

template<typename K, typename KExact = K>
class Refine_face_comparer;

template<typename K, typename KExact = K>
class Refine_face_iterator_comparer;

template<typename K, typename KExact = K>
class Refine_face_hash;

template<typename K, typename KExact = K>
class Refine_face
{
  typedef Refine_face<K, KExact>                           Self;

public:
  typedef typename K::FT                                    FT;
  typedef Stretched_Delaunay_2<K, KExact>                   Star;
  typedef Star*                                             Star_handle;
  typedef typename Star::Face                               Face;
  typedef typename Star::Face_handle                        Face_handle;

  typedef Refine_face_comparer<K, KExact>                   Rface_comparer;
  typedef Refine_face_iterator_comparer<K, KExact>          Rface_it_comparer;
  typedef Refine_face_hash<K, KExact>                       Rface_hash;

#ifdef ANISO_USE_BOOST_UNORDERED_SET
  typedef boost::unordered_set<Self,
                               Rface_hash,
                               Rface_comparer>              Rface_set;
#else
  typedef std::set<Self, Rface_comparer>                    Rface_set;
#endif

  typedef typename Rface_set::iterator                      Rface_set_iterator;
  typedef std::multiset<Rface_set_iterator,
                        Rface_it_comparer>                  Queue;
  typedef typename Queue::iterator                          Queue_iterator;

public:
  Star_handle star;
  Facet_ijk face;

  FT value;
  int queue_type; // type of queue
  mutable Queue_iterator queue_it; // position in queues[queue_type]
  bool prev_rejection;

public:
  Refine_face() : star(NULL), value(0), queue_type(-1), prev_rejection(false) { }
  Refine_face(Star_handle star_, const Face_handle& face_, FT value_, int queue_type_)
    :
    star(star_),
    face(face_),
    value(value_),
    queue_type(queue_type_),
    prev_rejection(false)
  { }
};

template<typename K, typename KExact>
inline std::ostream& operator<<(std::ostream& os, const Refine_face<K, KExact>& src)
{
  os << src.star->index_in_star_set() << " || ";
  os << src.face[0] << " " << src.face[1] << " " << src.face[2];
  os << " || qt: " << src.queue_type;
  os << " || val: " << src.value;
  os << " || rej: " << src.prev_rejection << std::endl;
  return os;
}

template<typename K, typename KExact/* = K */>
class Refine_face_comparer
{
public:
  Refine_face_comparer() { }
  bool operator() (const Refine_face<K, KExact>& left,
                   const Refine_face<K, KExact>& right) const
  {
    CGAL_precondition(left.face.vertices().size() == right.face.vertices().size());
#ifdef ANISO_USE_BOOST_UNORDERED_SET
    return std::equal(left.face.vertices().begin(),
                      left.face.vertices().end(),
                      right.face.vertices().begin());
#else
    return std::lexicographical_compare(left.face.vertices().begin(),
                                        left.face.vertices().end(),
                                        right.face.vertices().begin(),
                                        right.face.vertices().end());
#endif
  }
};

template<typename K, typename KExact>
class Refine_face_hash
{
public:
  typedef Refine_face<K, KExact>                           Rface;

  Refine_face_hash(){}

  std::size_t operator()(const Rface& rf) const
  {
    return boost::hash_range(rf.face.vertices().begin(),
                             rf.face.vertices().end());
  }
};

template<typename K, typename KExact/* = K*/>
class Refine_face_iterator_comparer
{
  typedef Refine_face<K, KExact>                           Rface;
  typedef Refine_face_comparer<K, KExact>                  Rface_comparer;

#ifdef ANISO_USE_BOOST_UNORDERED_SET
  typedef boost::unordered_set<Rface,
                               Refine_face_hash<K, KExact>,
                               Rface_comparer>             Rface_set;
#else
  typedef std::set<Rface, Rface_comparer>                  Rface_set;
#endif

  typedef typename Rface_set::iterator                     Rface_set_iterator;
  typedef typename Rface_set::const_iterator               Rface_set_c_iterator;

public:
  Refine_face_iterator_comparer() { }

  bool operator() (Rface_set_iterator left,
                   Rface_set_iterator right) const
  {
    return left->value > right->value;
  }
};

template<typename K, typename KExact = K>
class Face_refine_queue
{
public:
  typedef typename K::FT                                    FT;
  typedef typename K::Point_2                               Point_2;

  typedef Refine_face<K, KExact>                            Rface;
  typedef Refine_face_comparer<K, KExact>                   Rface_comparer;
  typedef Refine_face_iterator_comparer<K, KExact>          Rface_it_comparer;

#ifdef ANISO_USE_BOOST_UNORDERED_SET
  typedef boost::unordered_set<Rface,
                               Refine_face_hash<K, KExact>,
                               Rface_comparer>               Rface_set;
#else
  typedef std::set<Rface, Rface_comparer>                    Rface_set;
#endif

  typedef typename Rface_set::iterator                       Rface_set_iterator;

  typedef std::multiset<Rface_set_iterator,
                        Rface_it_comparer>                   Queue;
  typedef typename Queue::iterator                           Queue_iterator;

  typedef typename Rface::Star_handle                        Star_handle;
  typedef typename Rface::Face                               Face;
  typedef typename Rface::Face_handle                        Face_handle;

public:
  static const int nb_queues = 6;
  static const int encroachment_queue = 0;
  static const int over_distortion_queue = 1;
  static const int over_circumradius_queue = 2;
  static const int bad_shape_queue = 3;
  static const int start_pick_valid = 4;
  static const int custom_inconsistent_queue = 4;
  static const int inconsistent_queue = 5;

  //rfaces: set of the faces needing refinement. Avoids having the same face (i,j,k)
  //         in different queues (either due to diff values, star or queue_type).

  //queues: multisets of iterators to the rfaces set, sorted by value. Multiset chosen
  //        over prio_queue to access/remove any element.

private:
  Rface_set rfaces;
  Queue *queues[nb_queues];

  Queue encroachments;
  Queue over_distortions;
  Queue over_circumradii;
  Queue over_approximation;
  Queue bad_shapes;
  Queue custom_inconsistent;
  Queue inconsistent;

private:
  void update_rface(Rface_set_iterator& rface_it,
                     Star_handle star,
                     FT value,
                     int queue_type,
                     bool prev_rejection)
  {
/*
    std::cout << "update rface " << *rface_it << std::endl;
    std::cout << "with new values" << std::endl;
    std::cout << star->index_in_star_set() << " " << value << " " << queue_type << std::endl;
*/

    Rface_set_iterator rface_hint = rfaces.end();
    Queue_iterator queue_hint = queues[queue_type]->end();

#ifndef ANISO_USE_BOOST_UNORDERED_SET
    if(rface_it == rfaces.begin())
      rface_hint = rfaces.end(); //can't use hint in that case, end() is safe
    else
      rface_hint = (--rface_it)++;
#endif

    //In the general case, we don't know where the insertion will be.
    // In the rejection case, the insertion is done at the end.
    if(!prev_rejection || queues[queue_type]->size() <= 1)
      queue_hint = queues[queue_type]->end();
    else
      queue_hint = --(queues[queue_type]->end());

    Rface new_rface = *rface_it;
    new_rface.star = star;
    new_rface.value = value;
    new_rface.queue_type = queue_type;
    new_rface.prev_rejection = prev_rejection;

    //clean both queues of that element
    queues[rface_it->queue_type]->erase(rface_it->queue_it);

#ifdef ANISO_USE_BOOST_UNORDERED_SET
    rface_hint =
#endif
    rfaces.erase(rface_it);

    //re-insert with the new values
    rface_it = rfaces.insert(rface_hint, new_rface);
    Queue_iterator qit = queues[queue_type]->insert(queue_hint, rface_it);

    rface_it->queue_it = qit;
  }

  void update_rface_value(Rface_set_iterator& rface_it,
                           FT value)
  {
    return update_rface(rface_it, rface_it->star, value,
                         rface_it->queue_type, rface_it->prev_rejection);
  }

public:
  //reject the top of the queue n° queue_type to the end of the same queue
  void reject_rface(int queue_type)
  {
    Rface_set_iterator rface_it = *(queues[queue_type]->begin());
    FT min_value = queue_min_value(queue_type);
    FT epsilon = 1e-10;
    return update_rface(rface_it, rface_it->star, min_value - epsilon,
                         rface_it->queue_type, true /*rejection*/);
  }

  void push(Star_handle star,
            const Face_handle& fh,
            FT value,
            int queue_type,
            bool force_push = false) //update even if it means a lower priority
  {
    Rface new_rface(star, fh, value, queue_type);

    //std::cout << "rface created ready to push : " << new_rface << std::endl;

    std::pair<Rface_set_iterator, bool> is_insert_successful;
    is_insert_successful = rfaces.insert(new_rface);

    if(!is_insert_successful.second) // the face already exists in one of the queues
    {
      Rface_set_iterator old_rface_it = is_insert_successful.first;
      if(old_rface_it->queue_type > queue_type || // the new face has a higher prio queue_type
         old_rface_it->value < value || // new value has higher prio within the same queue
         force_push)
      {
        /*
        std::cout << "refine_face needs update: " << std::endl;
        std::cout << old_rface_it->queue_type << " " << queue_type << std::endl;
        std::cout << old_rface_it->value << " " << value << std::endl;
        */
        update_rface(old_rface_it, star, value, queue_type, old_rface_it->prev_rejection);
      }
    }
    else
    {
      //std::cout << "insert successful" << std::endl;
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

  bool top(Rface_set_iterator& face,
           int start_id = 0, int end_id = nb_queues - 1) const
  {
    for (int i=start_id; i<=end_id; i++)
      if(!queues[i]->empty())
      {
        face = *(queues[i]->begin());
/*
        std::cout << "found top at: " << i << std::endl;
        std::cout << "star: " << face->star->index_in_star_set() << std::endl;
        std::cout << face->face[0] << " " << face->face[1] << " " << face->face[2] << std::endl;
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
        rfaces.erase(*(queues[i]->begin()));
        queues[i]->erase(queues[i]->begin());
        return true;
      }
    return false;
  }

  void clear(int start_id = 0, int end_id = nb_queues - 1)
  {
    rfaces.clear();
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
    int size = rfaces.size();

    Rface_set_iterator rfit = rfaces.begin();
    Rface_set_iterator rfend = rfaces.end();
    for(; rfit != rfend; ++rfit)
    {
      if(!rfit->star->has_face(rfit->face))
      {
        queues[rfit->queue_type]->erase(rfit->queue_it);
        rfaces.erase(rfit);
      }
    }

    t.stop();
    std::cout << "clean: " << size << " " << rfaces.size() << " in " << t.time() << std::endl;
  }

public:
  void print_rfaces()
  {
    std::cout << "rfaces: " << std::endl;
    Rface_set_iterator it = rfaces.begin();
    Rface_set_iterator itend = rfaces.end();
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
    std::cout << rfaces.size() << std::endl;
  }

  std::size_t count(int start_id = 0, int end_id = nb_queues - 1) const
  {
    std::size_t count = 0;
    for(int i=start_id; i<= end_id; ++i)
      count += queues[i]->size();
    return count;
  }

  bool is_face_in(Star_handle star, Face_handle face, FT value, int queue_type)
  {
    bool is_face_in, is_rest_identical;
    Rface_set_iterator it = rfaces.find(Rface(star, face, value, queue_type));

    if((is_face_in = (it != rfaces.end()))) // = and not == on purpose
    {
      is_rest_identical = (star->index_in_star_set() == it->star->index_in_star_set() &&
                           queue_type == it->queue_type &&
                           std::abs(value - it->value) < 1e-11);
    }

    if(is_face_in && !is_rest_identical)
    {
      std::cout << "face is in, but from a different star: " << std::endl;
      std::cout << "face check: ";
      std::cout << face->vertex(0)->info() << " ";
      std::cout << face->vertex(1)->info() << " ";
      std::cout << face->vertex(2)->info() << std::endl;
      std::cout << it->face[0] << " " << it->face[1] << " " << it->face[2] << std::endl;
      std::cout << "stars: " << star->index_in_star_set() << " " << it->star->index_in_star_set() << std::endl;
      std::cout << "values: " << value << " " << it->value << std::endl;
      std::cout << "queue type: " << queue_type << " " << it->queue_type << std::endl;
    }

    return is_face_in;
  }

public:
  Face_refine_queue()
    :
      rfaces(),
      encroachments(Rface_it_comparer()),
      over_distortions(Rface_it_comparer()),
      over_circumradii(Rface_it_comparer()),
      bad_shapes(Rface_it_comparer()),
      custom_inconsistent(Rface_it_comparer()),
      inconsistent(Rface_it_comparer())
  {
    queues[encroachment_queue] = &encroachments;
    queues[over_distortion_queue] = &over_distortions;
    queues[over_circumradius_queue] = &over_circumradii;
    queues[bad_shape_queue] = &bad_shapes;
    queues[custom_inconsistent_queue] = &custom_inconsistent;
    queues[inconsistent_queue] = &inconsistent;
  }
}; // face_refine_queue

} // Anisotropic_mesh_2
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_2_REFINE_QUEUE_SURFACE_H

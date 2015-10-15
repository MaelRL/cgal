#ifndef CGAL_ANISOTROPIC_MESH_3_CANVAS_REFINE_H
#define CGAL_ANISOTROPIC_MESH_3_CANVAS_REFINE_H

#include <CGAL/Canvas/primal_simplex.h>

#include <CGAL/Bbox_3.h>

#include <boost/unordered_set.hpp>
#include <boost/functional/hash.hpp>

#include <algorithm>
#include <cstddef>
#include <ostream>
#include <set>
#include <utility>

// the structure is the same as local Delaunay's refinement queues :
// we have a set of the simplices that are to be refined and the queues only contain
// iterators to those elements. This is made to avoid having the same elements multiple
// times in different queues or with different values.

// In local Delaunay, this structure is quite important to have because
// each simplex usually exists in multiple stars . Is it really that useful
// to have such a complex structure here ? Probably not.

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename Canvas, std::size_t size>
struct Queue_entry_hash;

template<typename Canvas, std::size_t size>
struct Queue_entry_comparer;

template<typename Iterator>
struct Queue_value_comparer;

template<typename Canvas, std::size_t size>
struct Queue_entry
{
private:
  typedef Queue_entry<Canvas, size>                                   Self;
public:
  typedef typename Canvas::FT                                         FT;
  typedef typename Canvas::Canvas_point                               Canvas_point;
  typedef typename Canvas::Point_3                                    Point_3;

  typedef Primal_simplex<Canvas, size>                                PSimplex;
  typedef typename PSimplex::BSimplex                                 Simplex;

  typedef boost::unordered_set<Self,
                               Queue_entry_hash<Canvas, size>,
                               Queue_entry_comparer<Canvas, size> >   Queues_entries;
  typedef typename Queues_entries::iterator                           Queue_entry_iterator;
  typedef std::multiset<Queue_entry_iterator,
                        Queue_value_comparer<Queue_entry_iterator> >  Queue;
  typedef typename Queue::iterator                                    Queue_iterator;

  PSimplex primal_simplex;
  FT value;
  std::size_t queue_type;
  mutable Queue_iterator queue_it;

  const Simplex& simplex() const { return primal_simplex.simplex(); }
  const Canvas_point* dual_point() const { return primal_simplex.dual_point(); }

  Queue_entry(const PSimplex& psimplex_, const FT value_, const std::size_t queue_type_)
    :
      primal_simplex(psimplex_),
      value(value_),
      queue_type(queue_type_)
  { }
};

template<typename Canvas, std::size_t size>
inline std::ostream& operator<<(std::ostream& os, const Queue_entry<Canvas, size>& qe)
{
  os << qe.primal_simplex;
  os << " || val: " << qe.value;
  os << " || qt: " << qe.queue_type << std::endl;
  return os;
}

template<typename Canvas, std::size_t size>
struct Queue_entry_hash
{
  typedef Queue_entry<Canvas, size>       Q_entry;

  std::size_t operator()(const Q_entry& e) const
  {
    return boost::hash_range(e.simplex().begin(), e.simplex().end());
  }

  Queue_entry_hash() { }
};

template<typename Canvas, std::size_t size>
struct Queue_entry_comparer
{
  typedef Queue_entry<Canvas, size>       Q_entry;

  bool operator()(const Q_entry& left, const Q_entry& right) const
  {
    const typename Q_entry::Simplex& s1 = left.simplex();
    const typename Q_entry::Simplex& s2 = right.simplex();

    CGAL_assertion(s1.size() == s2.size());

    return std::equal(s1.begin(), s1.end(), s2.begin());
  }

  Queue_entry_comparer() { }
};

template<typename Iterator>
struct Queue_value_comparer
{
  bool operator()(const Iterator& left, const Iterator& right) const
  {
    // order by decreasing value
    return left->value > right->value;
  }

  Queue_value_comparer() { }
};

template<typename Canvas, std::size_t size>
class Canvas_refinement_queue
{
public:
  typedef typename Canvas::Kernel                                      Kernel;
  typedef typename Kernel::FT                                          FT;
  typedef typename Canvas::Canvas_point                                Canvas_point;

  typedef Queue_entry<Canvas, size>                                    Q_entry;
  typedef Primal_simplex<Canvas, size>                                 PSimplex;
  typedef boost::unordered_set<Q_entry,
                               Queue_entry_hash<Canvas, size>,
                               Queue_entry_comparer<Canvas, size> >    Queues_entries;
  typedef typename Queues_entries::iterator                            Queue_entry_iterator;
  typedef std::multiset<Queue_entry_iterator,
                        Queue_value_comparer<Queue_entry_iterator> >   Queue;
  typedef typename Queue::iterator                                     Queue_value_iterator;

  static const int nb_queues = 4;
  static const int over_distortion_queue = 0;
  static const int over_circumradius_queue = 1;
  static const int bad_shape_queue = 2;
  static const int intersections_queue = 3;

  Queues_entries queues_entries;

  Queue *queues[nb_queues];

  Queue distortion_queue;
  Queue size_queue;
  Queue quality_queue;
  Queue intersection_queue;

  void update_queue_entry(Queue_entry_iterator& qe_it,
                          const PSimplex& psimplex,
                          const FT val, const std::size_t q_type)
  {
    Q_entry new_entry = *qe_it;
    new_entry.primal_simplex = psimplex;
    new_entry.value = val;
    new_entry.queue_type = q_type;

    // clean both queues of that element
    queues[qe_it->queue_type]->erase(qe_it->queue_it);
    Queue_entry_iterator entry_hint = queues_entries.erase(qe_it);

    // re-insert with the new values
    qe_it = queues_entries.insert(entry_hint, new_entry);
    Queue_value_iterator qit = queues[q_type]->insert(qe_it);

    qe_it->queue_it = qit;
  }

  void push(const PSimplex& psimplex,
            const FT val,
            const std::size_t q_type)
  {
    Q_entry new_entry(psimplex, val, q_type);

    std::pair<Queue_entry_iterator, bool> is_insert_successful;
    is_insert_successful = queues_entries.insert(new_entry);

    if(!is_insert_successful.second) // already exists
    {
      Queue_entry_iterator old_qe_it = is_insert_successful.first;
      if(old_qe_it->queue_type > q_type || // the new face has a higher prio queue_type
         old_qe_it->value < val) // new value has higher prio within the same queue
      {
        update_queue_entry(old_qe_it, psimplex, val, q_type);
      }
    }
    else
    {
      is_insert_successful.first->queue_it =
                             queues[q_type]->insert(is_insert_successful.first);
    }
  }

  bool empty(int start_id = 0, int end_id = nb_queues - 1) const
  {
    for (int i=start_id; i<=end_id; i++)
    {
      if (!queues[i]->empty())
        return false;
    }
    return true;
  }

  bool top(Queue_entry_iterator& e, int start_id = 0, int end_id = nb_queues - 1) const
  {
    for (int i=start_id; i<=end_id; i++)
    {
      if(!queues[i]->empty())
      {
        e = *(queues[i]->begin());
        return true;
      }
    }
    return false;
  }

  bool pop(int start_id = 0, int end_id = nb_queues - 1)
  {
    for (int i=start_id; i<=end_id; i++)
    {
      if (!queues[i]->empty())
      {
        queues_entries.erase(*(queues[i]->begin())); // got to be careful with the order
        queues[i]->erase(queues[i]->begin());
        return true;
      }
    }
    return false;
  }

  void clean()
  {
    // When we insert a point, we change the connectivity of the dual
    // and some simplices that were marked for refinement possibly disappears.
    // We have to clean them off the queues

    // Ideally it should be rather simplex: when we're spreading for a new point
    // over an already colored canvas, we detect which simplex we're overwriting
    // and delete them from the queue... Tedious to code so todo...
  }

  void clear(int start_id = 0, int end_id = nb_queues - 1)
  {
    queues_entries.clear();
    for(int i=start_id; i<=end_id; i++)
      queues[i]->clear();
  }

public:
  void print_entries()
  {
    std::cout << "rfaces: " << std::endl;
    Queue_entry_iterator qeit = queues_entries.begin();
    Queue_entry_iterator qeend = queues_entries.end();
    for(; qeit!=qeend; ++qeit)
      std::cout << *qeit;
  }

  void print_queue(int queue_type)
  {
    std::cout << "queue : " << queue_type
              << " (size: " << queues[queue_type]->size() << ")" << std::endl;
    Queue* si = queues[queue_type];
    Queue_value_iterator vit = si->begin();
    Queue_value_iterator viend = si->end();
    for(; vit!=viend; ++vit)
      std::cout << **vit;
  }

  void print_queues(int start_id = 0, int end_id = nb_queues - 1)
  {
    for(int i=start_id; i<= end_id; ++i)
      print_queue(i);
  }

  Canvas_refinement_queue()
    :
      distortion_queue(),
      size_queue(),
      quality_queue(),
      intersection_queue()
  {
    queues[over_distortion_queue] = &distortion_queue;
    queues[over_circumradius_queue] = &size_queue;
    queues[bad_shape_queue] = &quality_queue;
    queues[intersections_queue] = &intersection_queue;
  }
};

}  // namespace Anisotropic_mesh_3
}  // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_CANVAS_REFINE_H

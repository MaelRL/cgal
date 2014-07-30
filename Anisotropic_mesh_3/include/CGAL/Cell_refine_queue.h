#ifndef CGAL_ANISOTROPIC_MESH_3_CELL_REFINE_QUEUE_H
#define CGAL_ANISOTROPIC_MESH_3_CELL_REFINE_QUEUE_H

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
class Refine_cell_comparer;

template<typename K, typename KExact = K>
class Refine_cell_iterator_comparer;

template<typename K, typename KExact = K>
class Refine_cell_hash;

template<typename K, typename KExact = K>
class Refine_cell
{
  typedef Refine_cell<K, KExact>                            Self;

public:
  typedef typename K::FT                                    FT;
  typedef Stretched_Delaunay_3<K>                           Star;
  typedef Star*                                             Star_handle;
  typedef typename Star::Cell_handle                        Cell_handle;

  typedef Refine_cell_comparer<K, KExact>                   Rcell_comparer;
  typedef Refine_cell_iterator_comparer<K, KExact>          Rcell_it_comparer;
  typedef Refine_cell_hash<K, KExact>                       Rcell_hash;

#ifdef ANISO_USE_BOOST_UNORDERED_SET
  typedef typename boost::unordered_set<Self,
                                        Rcell_hash,
                                        Rcell_comparer>     Rcell_set;
#else
  typedef typename std::set<Self, Rcell_comparer>           Rcell_set;
#endif

  typedef typename Rcell_set::iterator                      Rcell_set_iterator;
  typedef typename std::multiset<Rcell_set_iterator,
                                 Rcell_it_comparer>         Queue;
  typedef typename Queue::iterator                          Queue_iterator;

public:
  Star_handle star;
  Cell_ijkl cell;

  FT  value;
  int queue_type; // type of queue
  mutable Queue_iterator queue_it; //position in queues[queue_type]
  bool prev_rejection;

public:
  Refine_cell() : star(NULL), value(0), queue_type(-1), prev_rejection(false) { }
  Refine_cell(Star_handle star_, Cell_handle cell_, FT value_, int queue_type_)
    :
    star(star_),
    cell(cell_),
    value(value_),
    queue_type(queue_type_),
    prev_rejection(false)
  { }
};

template<typename K, typename KExact>
inline std::ostream& operator<<(std::ostream& os, const Refine_cell<K, KExact>& src)
{
  os << src.star->index_in_star_set() << " || ";
  os << src.cell[0] << " " << src.cell[1] << " " << src.cell[2] << " " << src.cell[3];
  os << " || val: " << src.queue_type;
  os << " || val: " << src.value;
  os << " || rej: " << src.prev_rejection << std::endl;
  return os;
}

template<typename K, typename KExact/* = K*/>
class Refine_cell_comparer
{
public:
  Refine_cell_comparer() { }
  bool operator() (const Refine_cell<K, KExact>& left,
                   const Refine_cell<K, KExact>& right) const
  {
#ifdef ANISO_USE_BOOST_UNORDERED_SET
    return std::equal(left.cell.vertices().begin(), left.cell.vertices().end(), right.cell.vertices().begin());
#else
    return std::lexicographical_compare(left.cell.vertices().begin(),
                                        left.cell.vertices().end(),
                                        right.cell.vertices().begin(),
                                        right.cell.vertices().end());
#endif
  }
};

template<typename K, typename KExact>
class Refine_cell_hash
{
public:
  typedef Refine_cell<K, KExact>                           Rcell;

  Refine_cell_hash(){}

  std::size_t operator()(const Rcell& rc) const
  {
    return boost::hash_range(rc.cell.vertices().begin(),
                             rc.cell.vertices().end());
  }
};


template<typename K, typename KExact/* = K*/>
class Refine_cell_iterator_comparer
{
  typedef Refine_cell<K, KExact>                            Rcell;
  typedef Refine_cell_comparer<K, KExact>                   Rcell_comparer;
  typedef Refine_cell_hash<K, KExact>                       Rcell_hash;

#ifdef ANISO_USE_BOOST_UNORDERED_SET
  typedef typename boost::unordered_set<Rcell,
                                        Rcell_hash,
                                        Rcell_comparer>     Rcell_set;
#else
  typedef typename std::set<Rcell, Rcell_comparer>          Rcell_set;
#endif

  typedef typename Rcell_set::iterator                      Rcell_set_iterator;

public:
  Refine_cell_iterator_comparer() { }
  bool operator() (const Rcell_set_iterator left,
                   const Rcell_set_iterator right) const
  {
    return left->value > right->value;
  }
};

template<typename K, typename KExact = K>
class Cell_refine_queue
{
public:
  typedef typename K::FT                                        FT;
  typedef typename K::Point_3                                   Point_3;

  typedef Refine_cell<K, KExact>                                Rcell;
  typedef Refine_cell_comparer<K, KExact>                       Rcell_comparer;
  typedef Refine_cell_iterator_comparer<K, KExact>              Rcell_it_comparer;
  typedef Refine_cell_hash<K, KExact>                           Rcell_hash;

#ifdef ANISO_USE_BOOST_UNORDERED_SET
  typedef typename boost::unordered_set<Rcell,
                                        Rcell_hash,
                                        Rcell_comparer>         Rcell_set;
#else
  typedef typename std::set<Rcell, Rcell_comparer>              Rcell_set;
#endif

  typedef typename Rcell_set::iterator                          Rcell_set_iterator;
  typedef typename std::multiset<Rcell_set_iterator,
                                 Rcell_it_comparer>             Queue;
  typedef typename Queue::iterator                              Queue_iterator;

  typedef typename Rcell::Star_handle                           Star_handle;
  typedef typename Rcell::Cell_handle                           Cell_handle;

public:
  static const int nb_queues = 6;
  static const int encroachment_queue = 0;
  static const int over_distortion_queue = 1;
  static const int over_circumradius_queue = 2;
  static const int bad_shape_queue = 3;
  static const int start_pick_valid = 4;
  static const int sliver_queue = 4;
  static const int inconsistent_queue = 5;

  //rcells: set of the cells needing refinement. Avoids having the same cell (i,j,k)
  //         in different queues (either due to diff values, star or queue_type).

  //queues: multisets of iterators to the rcells set, sorted by value. Multiset chosen
  //        over prio_queue to access/remove any element.

private:
  Rcell_set rcells;
  Queue *queues[nb_queues];

  Queue encroachments;
  Queue over_distortions;
  Queue over_circumradii;
  Queue bad_shapes;
  Queue slivers;
  Queue inconsistents;

private:
  void update_rcell(Rcell_set_iterator& rcell_it,
                    Star_handle star,
                    FT value,
                    int queue_type,
                    bool prev_rejection)
  {
/*
    std::cout << "update rcell " << *rcell_it << std::endl;
    std::cout << "with new values" << std::endl;
    std::cout << star->index_in_star_set() << " " << value << " " << queue_type << std::endl;
*/

    Rcell_set_iterator rcell_hint = rcells.end();
    Queue_iterator queue_hint = queues[queue_type]->end();

#ifndef ANISO_USE_BOOST_UNORDERED_SET
    if(rcell_it == rcells.begin())
      rcell_hint = rcells.end(); //can't use hint in that case, end() is safe
    else
      rcell_hint = (--rcell_it)++;
#endif

    //In the general case, we don't know where the insertion will be. In the rejection case
    //the insertion is done at the end.
    if(!prev_rejection || queues[queue_type]->size() <= 1)
      queue_hint = queues[queue_type]->end();
    else
      queue_hint = --(queues[queue_type]->end());

    Rcell new_rcell = *rcell_it;
    new_rcell.star = star;
    new_rcell.value = value;
    new_rcell.queue_type = queue_type;
    new_rcell.prev_rejection = prev_rejection;

    //clean both queues of that element
    queues[rcell_it->queue_type]->erase(rcell_it->queue_it);

#ifdef ANISO_USE_BOOST_UNORDERED_SET
    rcell_hint =
#endif
    rcells.erase(rcell_it);

    //re-insert with the new values
    rcell_it = rcells.insert(rcell_hint, new_rcell);
    Queue_iterator qit = queues[queue_type]->insert(queue_hint, rcell_it);

    rcell_it->queue_it = qit;
  }

  void update_rcell_value(Rcell_set_iterator& rcell_it,
                          FT value)
  {
    return update_rcell(rcell_it, rcell_it->star, value,
                         rcell_it->queue_type, rcell_it->prev_rejection);
  }

public:
  //reject the top of the queue n° queue_type to the end of the same queue
  void reject_rcell(int queue_type)
  {
    Rcell_set_iterator rcell_it = *(queues[queue_type]->begin());
    FT min_value = queue_min_value(queue_type);
    FT epsilon = 1e-10;
    return update_rcell(rcell_it, rcell_it->star, min_value - epsilon,
                        rcell_it->queue_type, true /*rejection*/);
  }

  void push(Star_handle star,
            const Cell_handle cell,
            FT value,
            int queue_type,
            bool force_push = false) //update even if it means a lower priority
  {
    Rcell new_rcell(star, cell, value, queue_type);

    //std::cout << "rcell created ready to push: "  << new_rcell << std::endl;

    std::pair<Rcell_set_iterator, bool> is_insert_successful;
    is_insert_successful = rcells.insert(new_rcell);

    if(!is_insert_successful.second) // the cell already exists in one of the queues
    {
      //std::cout << "cell already in" << std::endl;
      Rcell_set_iterator old_rcell_it = is_insert_successful.first;
      if(old_rcell_it->queue_type > queue_type || // the new cell has a higher prio queue_type
         old_rcell_it->value < value || // new value has higher prio within the same queue
         force_push)
      {
        /*
        std::cout << "refine_cell needs update: " << std::endl;
        std::cout << old_rcell_it->queue_type << " " << queue_type << std::endl;
        std::cout << old_rcell_it->value << " " << value << std::endl;
        */
        update_rcell(old_rcell_it, star, value, queue_type, old_rcell_it->prev_rejection);
      }
    }
    else
    {
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

  bool top(Rcell_set_iterator& cell,
           int start_id = 0, int end_id = nb_queues - 1) const
  {
    for (int i=start_id; i<=end_id; i++)
      if (!queues[i]->empty())
      {
        cell = *(queues[i]->begin());
/*
        std::cout << "found cell top at: " << i << std::endl;
        std::cout << "star: " << cell->star->index_in_star_set() << std::endl;
        std::cout << cell->cell[0] << " " << cell->cell[1] << " " << cell->cell[2] << " " << cell->cell[3] << std::endl;
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
        rcells.erase(*(queues[i]->begin()));
        queues[i]->erase(queues[i]->begin());
        return true;
      }
    return false;
  }

  void clear(int start_id = 0, int end_id = nb_queues - 1)
  {
    rcells.clear();
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
    int size = rcells.size();

    Rcell_set_iterator rcit = rcells.begin();
    Rcell_set_iterator rcend = rcells.end();
    for(; rcit != rcend; ++rcit)
    {
      Cell_handle c;
      if(!rcit->star->has_cell(rcit->cell, c))
      {
        queues[rcit->queue_type]->erase(rcit->queue_it);
        rcells.erase(rcit);
      }
    }

    t.stop();
    std::cout << "clean: " << size << " " << rcells.size() << " in " << t.time() << std::endl;
  }

public:
  void print_rcells()
  {
    std::cout << "rcells: " << std::endl;
    Rcell_set_iterator it = rcells.begin();
    Rcell_set_iterator itend = rcells.end();
    for(;it!=itend; ++it)
      std::cout << *it;
  }

  void print_queue(int queue_type)
  {
    std::cout << "cqueue : " << queue_type << std::endl;
    Queue* si = queues[queue_type];
    Queue_iterator vit = si->begin();
    Queue_iterator viend = si->end();
    for(; vit!=viend; ++vit)
      std::cout << **vit;
  }

  void print_queues(int start_id = 0, int end_id = nb_queues - 1)
  {
    for(int i=start_id; i<=end_id; ++i)
      print_queue(i);
  }

  void print() const
  {
    std::cout << "CQueue : ( ";
    for(int i = 0; i < nb_queues; i++)
      std::cout << queues[i]->size() <<" ";
    std::cout << ") ";
    std::cout << rcells.size() << std::endl;
  }

  std::size_t count() const
  {
    std::size_t count = 0;
    for(int i=0; i< nb_queues; ++i)
      count += queues[i]->size();
    return count;
  }

  bool is_cell_in(Star_handle star, Cell_handle cell, FT value, int queue_type)
  {
    bool is_cell_in, is_rest_identical;

    Rcell_set_iterator it = rcells.find(Rcell(star, cell, value, queue_type));
    if((is_cell_in = (it != rcells.end()))) // = and not == on purpose
    {
      is_rest_identical = (star->index_in_star_set() == it->star->index_in_star_set() &&
                           queue_type == it->queue_type &&
                           std::abs(value  - it->value) < 1e-11);
    }

    if(is_cell_in && !is_rest_identical)
    {
      std::cout.precision(25);
      std::cout << "is_cell_in : unicity says in, but they differ: " << std::endl;
      std::cout << "cell check: ";
      std::cout << cell->vertex(0)->info() << " " << cell->vertex(1)->info() << " ";
      std::cout << cell->vertex(2)->info() << " " << cell->vertex(3)->info() << std::endl;
      std::cout << it->cell[0] << " " << it->cell[1] << " " << it->cell[2] << " " << it->cell[3] << std::endl;
      std::cout << "stars: " << star->index_in_star_set() << " " << it->star->index_in_star_set() << std::endl;
      std::cout << "values: " << value << " " << it->value << std::endl;
      std::cout << "queue type: " << queue_type << " " << it->queue_type << std::endl;
      std::cout.precision(8);
    }

    return is_cell_in;
  }

public:
  Cell_refine_queue()
    :
      rcells(),
      encroachments(Rcell_it_comparer()),
      over_distortions(Rcell_it_comparer()),
      over_circumradii(Rcell_it_comparer()),
      bad_shapes(Rcell_it_comparer()),
      slivers(Rcell_it_comparer()),
      inconsistents(Rcell_it_comparer())
  {
    queues[encroachment_queue] = &encroachments;
    queues[over_distortion_queue] = &over_distortions;
    queues[over_circumradius_queue] = &over_circumradii;
    queues[bad_shape_queue] = &bad_shapes;
    queues[sliver_queue] = &slivers;
    queues[inconsistent_queue] = &inconsistents;
  }
}; // Cell_refine_queue

} // Anisotropic_mesh_3
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_CELL_REFINE_QUEUE_H

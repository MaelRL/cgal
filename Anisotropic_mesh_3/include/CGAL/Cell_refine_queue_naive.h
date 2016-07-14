#ifndef CGAL_ANISOTROPIC_MESH_3_CELL_REFINE_QUEUE_H
#define CGAL_ANISOTROPIC_MESH_3_CELL_REFINE_QUEUE_H

#include <CGAL/Stretched_Delaunay_3.h>
#include <CGAL/helpers/combinatorics_helper.h>

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
class Refine_cell
{
  typedef Refine_cell<K, KExact>                            Self;

public:
  typedef typename K::FT                                    FT;
  typedef Stretched_Delaunay_3<K>                           Star;
  typedef Star*                                             Star_handle;
  typedef typename Star::Cell_handle                        Cell_handle;

  typedef Refine_cell_comparer<K, KExact>                   Rcell_comparer;

  typedef std::priority_queue<Self, 
                              typename std::vector<Self>,
                              Rcell_comparer>               Queue;

public:
  Star_handle star;
  Cell_ijkl cell;

  FT  value;
  int queue_type; // type of queue
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
    return left.value < right.value;
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

  typedef std::priority_queue<Rcell, 
                              typename std::vector<Rcell>,
                              Rcell_comparer>                   Queue;

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
  Queue *queues[nb_queues];

  Queue encroachments;
  Queue over_distortions;
  Queue over_circumradii;
  Queue bad_shapes;
  Queue slivers;
  Queue inconsistents;

public:
  //reject the top of the queue n° queue_type to the end of the same queue
  // void reject_rcell(int queue_type)
  // {
  //   Rcell_set_iterator rcell_it = *(queues[queue_type]->begin());
  //   FT min_value = queue_min_value(queue_type);
  //   FT epsilon = 1e-10;
  //   return update_rcell(rcell_it, rcell_it->star, min_value - epsilon,
  //                       rcell_it->queue_type, true /*rejection*/);
  // }

  void push(Star_handle star,
            const Cell_handle cell,
            FT value,
            int queue_type,
            bool force_push = false) //update even if it means a lower priority
  {
    Rcell new_rcell(star, cell, value, queue_type);
    queues[queue_type]->push(new_rcell);

//    std::cout << "pushed: star " << star->index_in_star_set() << " v: " << value << " in q: " << queue_type << std::endl;
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

  bool top(Rcell const *& cell,
           int start_id = 0, int end_id = nb_queues - 1) const
  {
    for (int i=start_id; i<=end_id; i++)
      if (!queues[i]->empty())
      {
        cell = &(queues[i]->top());

        // std::cout << "found cell top at: " << i << std::endl;
        // std::cout << "star: " << cell->star->index_in_star_set() << std::endl;
        // std::cout << cell->cell[0] << " " << cell->cell[1] << " " << cell->cell[2] << " " << cell->cell[3] << std::endl;

        return true;
      }
    return false;
  }

  bool pop(int start_id = 0, int end_id = nb_queues - 1)
  {
    for (int i=start_id; i<=end_id; i++)
      if (!queues[i]->empty())
      {
        queues[i]->pop();

//        std::cout << "popped queue: " << i << std::endl;

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
    return (--queues[queue_type]->end()).value;
  }

public:
  // void print_queue(int queue_type)
  // {
  //   std::cout << "cqueue : " << queue_type << std::endl;
  //   Queue* si = queues[queue_type];
  //   Queue_iterator vit = si->begin();
  //   Queue_iterator viend = si->end();
  //   for(; vit!=viend; ++vit)
  //     std::cout << *vit;
  // }

  // void print_queues(int start_id = 0, int end_id = nb_queues - 1)
  // {
  //   for(int i=start_id; i<=end_id; ++i)
  //     print_queue(i);
  // }

  void print() const
  {
    std::cout << "CQueue : ( ";
    for(int i = 0; i < nb_queues; i++)
      std::cout << queues[i]->size() <<" ";
    std::cout << ") ";
  }

  std::size_t count() const
  {
    std::size_t count = 0;
    for(int i=0; i< nb_queues; ++i)
      count += queues[i]->size();
    return count;
  }


//fixme there is no find() for priority queues so this needs a fix to work
//once it is fixed, the NAIVE macro can be removed in aniso_ref_cells_3
  // bool is_cell_in(Star_handle star, Cell_handle cell, FT value, int queue_type)
  // {
  //   bool is_cell_in, is_rest_identical;

  //   Rcell rc(star, cell, value, queue_type);

  //   Queue* q = queues[queue_type];
  //   Rcell* it = q.find(rc)

  //   if(it != q.end()) // = already exists
  //   {
  //     is_rest_identical = (star->index_in_star_set() == it->star->index_in_star_set() &&
  //                          queue_type == it->queue_type &&
  //                          std::abs(value  - it->value) < 1e-11);
  //   }

  //   if(is_cell_in && !is_rest_identical)
  //   {
  //     std::cout.precision(25);
  //     std::cout << "is_cell_in : unicity says in, but they differ: " << std::endl;
  //     std::cout << "cell check: ";
  //     std::cout << cell->vertex(0)->info() << " " << cell->vertex(1)->info() << " ";
  //     std::cout << cell->vertex(2)->info() << " " << cell->vertex(3)->info() << std::endl;
  //     std::cout << it->cell[0] << " " << it->cell[1] << " " << it->cell[2] << " " << it->cell[3] << std::endl;
  //     std::cout << "stars: " << star->index_in_star_set() << " " << it->star->index_in_star_set() << std::endl;
  //     std::cout << "values: " << value << " " << it->value << std::endl;
  //     std::cout << "queue type: " << queue_type << " " << it->queue_type << std::endl;
  //     std::cout.precision(8);
  //   }

  //   return is_cell_in;
  // }

public:
  Cell_refine_queue()
    :
      encroachments(Rcell_comparer()),
      over_distortions(Rcell_comparer()),
      over_circumradii(Rcell_comparer()),
      bad_shapes(Rcell_comparer()),
      slivers(Rcell_comparer()),
      inconsistents(Rcell_comparer())
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

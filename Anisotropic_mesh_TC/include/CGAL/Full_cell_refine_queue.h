#ifndef CGAL_ANISOTROPIC_MESH_TC_FULL_CELL_REFINE_QUEUE_H
#define CGAL_ANISOTROPIC_MESH_TC_FULL_CELL_REFINE_QUEUE_H

#include <CGAL/Star.h>
#include <CGAL/helpers/combinatorics_helper.h>

#ifndef ANISO_USE_BOOST_UNORDERED_SET
#define ANISO_USE_BOOST_UNORDERED_SET
#include <boost/functional/hash.hpp>
#include <boost/unordered_set.hpp>
#endif

#include <set>

namespace CGAL
{
namespace Anisotropic_mesh_TC
{

template<typename Kd, typename KD>
class Refine_full_cell_comparer;

template<typename Kd, typename KD>
class Refine_full_cell_iterator_comparer;

template<typename Kd, typename KD>
class Refine_full_cell_hash;

template<typename Kd, typename KD>
class Refine_full_cell
{
  typedef Refine_full_cell<Kd, KD>                          Self;

public:
  typedef typename Kd::FT                                   FT;
  typedef Tangent_star<Kd, KD>                              Star;
  typedef Star*                                             Star_handle;
  typedef typename Star::Full_cell                          Full_cell;
  typedef typename Star::Full_cell_handle                   Full_cell_handle;

  typedef Refine_full_cell_comparer<Kd, KD>                 Rfcell_comparer;
  typedef Refine_full_cell_iterator_comparer<Kd, KD>        Rfcell_it_comparer;
  typedef Refine_full_cell_hash<Kd, KD>                     Rfcell_hash;

  typedef typename Kd::Dimension                            dDim;
  typedef Ordered_simplex_base<dDim::value+1>               Simplex;

#ifdef ANISO_USE_BOOST_UNORDERED_SET
  typedef typename boost::unordered_set<Self,
                                        Rfcell_hash,
                                        Rfcell_comparer>    Rfcell_set;
#else
  typedef typename std::set<Self, Rfcell_comparer>          Rfcell_set;
#endif

  typedef typename Rfcell_set::iterator                     Rfcell_set_iterator;
  typedef typename std::multiset<Rfcell_set_iterator,
                                 Rfcell_it_comparer>        Queue;
  typedef typename Queue::iterator                          Queue_iterator;

public:
  Star_handle star;
  Simplex full_cell;

  FT value;
  int queue_type; // type of queue
  mutable Queue_iterator queue_it; // position in queues[queue_type]
  bool prev_rejection;

public:
  Refine_full_cell() : star(NULL), value(0), queue_type(-1), prev_rejection(false) { }
  Refine_full_cell(Star_handle star_,
                   const Full_cell_handle& fch_,
                   FT value_,
                   int queue_type_)
    :
    star(star_),
    value(value_),
    queue_type(queue_type_),
    prev_rejection(false)
  {
    //below is ugly, it should be done with an outside function to avoid overhead and copies
    std::vector<int> ids;
    for(int i=0; i<=dDim::value; ++i)
      ids.push_back(fch_->vertex(i)->data());
    full_cell = Simplex(ids.begin(), ids.end());
  }
};

template<typename Kd, typename KD>
inline std::ostream& operator<<(std::ostream& os,
                                const Refine_full_cell<Kd, KD>& src)
{
  os << "Star: " << src.star->index() << std::endl;
  os << src.full_cell;
  os << " || qt: " << src.queue_type;
  os << " || val: " << src.value;
  os << " || rej: " << src.prev_rejection << std::endl;
  return os;
}

template<typename Kd, typename KD>
class Refine_full_cell_comparer
{
public:
  Refine_full_cell_comparer() { }
  bool operator() (const Refine_full_cell<Kd, KD>& left,
                   const Refine_full_cell<Kd, KD>& right) const
  {
#ifdef ANISO_USE_BOOST_UNORDERED_SET
    return std::equal(left.full_cell.vertices().begin(),
                      left.full_cell.vertices().end(),
                      right.full_cell.vertices().begin());
#else
    return std::lexicographical_compare(left.full_cell.vertices().begin(),
                                        left.full_cell.vertices().end(),
                                        right.full_cell.vertices().begin(),
                                        right.full_cell.vertices().end());
#endif
  }
};

template<typename Kd, typename KD>
class Refine_full_cell_hash
{
public:
  typedef Refine_full_cell<Kd, KD>                           Rfcell;

  Refine_full_cell_hash(){}

  std::size_t operator()(const Rfcell& rf) const
  {
    return boost::hash_range(rf.full_cell.vertices().begin(),
                             rf.full_cell.vertices().end());
  }
};

template<typename Kd, typename KD>
class Refine_full_cell_iterator_comparer
{
  typedef Refine_full_cell<Kd, KD>                         Rfcell;
  typedef Refine_full_cell_comparer<Kd, KD>                Rfcell_comparer;

#ifdef ANISO_USE_BOOST_UNORDERED_SET
  typedef typename boost::unordered_set<Rfcell,
                                        Refine_full_cell_hash<Kd, KD>,
                                        Rfcell_comparer>   Rfcell_set;
#else
  typedef typename std::set<Rfcell, Rfcell_comparer>       Rfcell_set;
#endif

  typedef typename Rfcell_set::iterator                    Rfcell_set_iterator;
  typedef typename Rfcell_set::const_iterator              Rfcell_set_c_iterator;

public:
  Refine_full_cell_iterator_comparer() { }

  bool operator() (Rfcell_set_iterator left,
                   Rfcell_set_iterator right) const
  {
    return left->value > right->value;
  }
};

template<typename Kd, typename KD>
class Full_cell_refine_queue
{
public:
  typedef typename Kd::FT                                 FT;

  typedef Refine_full_cell<Kd, KD>                        Rfcell;
  typedef Refine_full_cell_comparer<Kd, KD>               Rfcell_comparer;
  typedef Refine_full_cell_iterator_comparer<Kd, KD>      Rfcell_it_comparer;

#ifdef ANISO_USE_BOOST_UNORDERED_SET
  typedef typename boost::unordered_set<Rfcell,
                                        Refine_full_cell_hash<Kd, KD>,
                                        Rfcell_comparer>  Rfcell_set;
#else
  typedef typename std::set<Rfcell, Rfcell_comparer>      Rfcell_set;
#endif

  typedef typename Rfcell_set::iterator                   Rfcell_set_iterator;

  typedef typename std::multiset<Rfcell_set_iterator,
                                 Rfcell_it_comparer>      Queue;
  typedef typename Queue::iterator                        Queue_iterator;

  typedef typename Rfcell::Star_handle                    Star_handle;
  typedef typename Rfcell::Full_cell                      Full_cell;
  typedef typename Rfcell::Full_cell_handle               Full_cell_handle;

public:
  static const int nb_queues = 6;
  static const int encroachment_queue = 0;
  static const int over_distortion_queue = 1;
  static const int over_circumradius_queue = 2;
  static const int bad_shape_queue = 3;
  static const int start_pick_valid = 4;
  static const int custom_inconsistent_queue = 4;
  static const int inconsistent_queue = 5;

  //rfull_cells: set of the faces needing refinement. Avoids having the same face (i,j,k)
  //         in different queues (either due to diff values, star or queue_type).

  //queues: multisets of iterators to the rfull_cells set, sorted by value. Multiset chosen
  //        over prio_queue to access/remove any element.

private:
  Rfcell_set rfull_cells;
  Queue *queues[nb_queues];

  Queue encroachments;
  Queue over_distortions;
  Queue over_circumradii;
  Queue over_approximation;
  Queue bad_shapes;
  Queue custom_inconsistent;
  Queue inconsistent;

private:
  void update_Rfcell(Rfcell_set_iterator& Rfcell_it,
                     Star_handle star,
                     FT value,
                     int queue_type,
                     bool prev_rejection)
  {
/*
    std::cout << "update Rfcell " << *Rfcell_it << std::endl;
    std::cout << "with new values" << std::endl;
    std::cout << star->index_in_star_set() << " " << value << " " << queue_type << std::endl;
*/

    Rfcell_set_iterator Rfcell_hint = rfull_cells.end();
    Queue_iterator queue_hint = queues[queue_type]->end();

#ifndef ANISO_USE_BOOST_UNORDERED_SET
    if(Rfcell_it == rfull_cells.begin())
      Rfcell_hint = rfull_cells.end(); //can't use hint in that case, end() is safe
    else
      Rfcell_hint = (--Rfcell_it)++;
#endif

    //In the general case, we don't know where the insertion will be. In the rejection case
    //the insertion is done at the end.
    if(!prev_rejection || queues[queue_type]->size() <= 1)
      queue_hint = queues[queue_type]->end();
    else
      queue_hint = --(queues[queue_type]->end());

    Rfcell new_Rfcell = *Rfcell_it;
    new_Rfcell.star = star;
    new_Rfcell.value = value;
    new_Rfcell.queue_type = queue_type;
    new_Rfcell.prev_rejection = prev_rejection;

    //clean both queues of that element
    queues[Rfcell_it->queue_type]->erase(Rfcell_it->queue_it);

#ifdef ANISO_USE_BOOST_UNORDERED_SET
    Rfcell_hint =
#endif
    rfull_cells.erase(Rfcell_it);

    //re-insert with the new values
    Rfcell_it = rfull_cells.insert(Rfcell_hint, new_Rfcell);
    Queue_iterator qit = queues[queue_type]->insert(queue_hint, Rfcell_it);

    Rfcell_it->queue_it = qit;
  }

  void update_Rfcell_value(Rfcell_set_iterator& Rfcell_it,
                           FT value)
  {
    return update_Rfcell(Rfcell_it, Rfcell_it->star, value,
                         Rfcell_it->queue_type, Rfcell_it->prev_rejection);
  }

public:
  //reject the top of the queue n° queue_type to the end of the same queue
  void reject_rfcell(int queue_type)
  {
    Rfcell_set_iterator Rfcell_it = *(queues[queue_type]->begin());
    FT min_value = queue_min_value(queue_type);
    FT epsilon = 1e-10;
    return update_Rfcell(Rfcell_it, Rfcell_it->star, min_value - epsilon,
                         Rfcell_it->queue_type, true /*rejection*/);
  }

  void push(Star_handle star,
            const Full_cell_handle& fhc,
            FT value,
            int queue_type,
            bool force_push = false) //update even if it means a lower priority
  {
    Rfcell new_Rfcell(star, fhc, value, queue_type);
    std::cout << "Rfcell created ready to push : " << new_Rfcell << std::endl;

    std::pair<Rfcell_set_iterator, bool> is_insert_successful;
    is_insert_successful = rfull_cells.insert(new_Rfcell);

    if(!is_insert_successful.second) // the face already exists in one of the queues
    {
      Rfcell_set_iterator old_Rfcell_it = is_insert_successful.first;
      if(old_Rfcell_it->queue_type > queue_type || // the new face has a higher prio queue_type
         old_Rfcell_it->value < value || // new value has higher prio within the same queue
         force_push)
      {
        /*
        std::cout << "Refine_full_cell needs update: " << std::endl;
        std::cout << old_Rfcell_it->queue_type << " " << queue_type << std::endl;
        std::cout << old_Rfcell_it->value << " " << value << std::endl;
        */
        update_Rfcell(old_Rfcell_it, star, value, queue_type, old_Rfcell_it->prev_rejection);
      }
    }
    else
    {
      std::cout << "push successful" << std::endl;
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

  bool top(Rfcell_set_iterator& face,
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
        rfull_cells.erase(*(queues[i]->begin()));
        queues[i]->erase(queues[i]->begin());
        return true;
      }
    return false;
  }

  void clear(int start_id = 0, int end_id = nb_queues - 1)
  {
    rfull_cells.clear();
    for(int i=start_id; i<=end_id; i++)
      queues[i]->clear();
  }

  FT queue_min_value(int queue_type)
  {
    return (*(--queues[queue_type]->end()))->value;
  }

  void clean()
  {
    Rfcell_set_iterator rfit = rfull_cells.begin();
    Rfcell_set_iterator rfend = rfull_cells.end();
    for(; rfit != rfend; ++rfit)
    {
      if(!rfit->star->has_face(rfit->face))
      {
        queues[rfit->queue_type]->erase(rfit->queue_it);
        rfull_cells.erase(rfit);
      }
    }
  }

public:
  void print_rfull_cells()
  {
    std::cout << "rfull_cells: " << std::endl;
    Rfcell_set_iterator it = rfull_cells.begin();
    Rfcell_set_iterator itend = rfull_cells.end();
    for(;it!=itend; ++it)
      std::cout << *it;
  }

  void print_queue(int queue_type)
  {
    std::cout << "fqueue : " << queue_type;
    std::cout << " (size: " << queues[queue_type]->size() << ")" << std::endl;
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
    std::cout << rfull_cells.size() << std::endl;
  }

  std::size_t count(int start_id = 0, int end_id = nb_queues - 1) const
  {
    std::size_t count = 0;
    for(int i=start_id; i<= end_id; ++i)
      count += queues[i]->size();
    return count;
  }

  bool is_full_cell_in(Star_handle star, Full_cell_handle fch,
                       FT value, int queue_type)
  {
    bool is_full_cell_in, is_rest_identical;
    Rfcell_set_iterator it = rfull_cells.find(Rfcell(star, fch, value, queue_type));

    if((is_full_cell_in = (it != rfull_cells.end()))) // = and not == on purpose
    {
      is_rest_identical = (star->index_in_star_set() == it->star->index_in_star_set() &&
                           queue_type == it->queue_type &&
                           std::abs(value - it->value) < 1e-11);
    }

    if(is_full_cell_in && !is_rest_identical)
    {
      std::cout << "face is in, but from a different star: " << std::endl;
      std::cout << "face check: ";
      std::cout << fch->vertex(0)->info() << " ";
      std::cout << fch->vertex(1)->info() << " ";
      std::cout << fch->vertex(2)->info() << std::endl;
      std::cout << it->face[0] << " " << it->face[1] << " " << it->face[2] << std::endl;
      std::cout << "stars: " << star->index_in_star_set() << " " << it->star->index_in_star_set() << std::endl;
      std::cout << "values: " << value << " " << it->value << std::endl;
      std::cout << "queue type: " << queue_type << " " << it->queue_type << std::endl;
    }

    return is_full_cell_in;
  }

public:
  Full_cell_refine_queue()
    :
      rfull_cells(),
      encroachments(Rfcell_it_comparer()),
      over_distortions(Rfcell_it_comparer()),
      over_circumradii(Rfcell_it_comparer()),
      bad_shapes(Rfcell_it_comparer()),
      custom_inconsistent(Rfcell_it_comparer()),
      inconsistent(Rfcell_it_comparer())
  {
    queues[encroachment_queue] = &encroachments;
    queues[over_distortion_queue] = &over_distortions;
    queues[over_circumradius_queue] = &over_circumradii;
    queues[bad_shape_queue] = &bad_shapes;
    queues[custom_inconsistent_queue] = &custom_inconsistent;
    queues[inconsistent_queue] = &inconsistent;
  }
}; // Full_cell_refine_queue

} // Anisotropic_mesh_2
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_2_REFINE_QUEUE_SURfcell_H

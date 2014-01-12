// Copyright (c) 2011  INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Author(s) : Kan-Le Shi

#ifndef CGAL_ANISOTROPIC_MESH_3_CELL_REFINE_QUEUE_H
#define CGAL_ANISOTROPIC_MESH_3_CELL_REFINE_QUEUE_H

#include <list>
#include <queue>
#include <vector>
#include <iostream>
#include <algorithm>
#include <CGAL/Stretched_Delaunay_3.h>

namespace CGAL {
namespace Anisotropic_mesh_3 {

template<typename K>
class Refine_cell
{
public:
  typedef typename K::FT                 FT;
  typedef Stretched_Delaunay_3<K>        Star;
  typedef Star                          *Star_handle;
  typedef typename Star::Cell_handle     Cell_handle;

public:
  Star_handle star;
  int vertices[4]; // the global ids, the first one is the tag
  FT  value;
public:
  Refine_cell() : star(NULL), value(0) { }
  Refine_cell(Star_handle star_, Cell_handle cell_, FT value_, int tag_) :
    star(star_), value(value_)
  {
    for (int i = 0; i < 4; i++)
      vertices[i] = cell_->vertex((tag_ + i) % 4)->info();
  }
  ~Refine_cell() { }
  bool operator==(Refine_cell<K> right)
  {
    return star == right.star &&
           vertices[0] == right.vertices[0] &&
           vertices[1] == right.vertices[1] &&
           vertices[2] == right.vertices[2] &&
           vertices[3] == right.vertices[3] &&
           value == right.value;
  }
};

template<typename K>
class Refine_cell_comparer
{
public:
  Refine_cell_comparer() { }
  bool operator() (Refine_cell<K> &left, Refine_cell<K> &right) const
  {
    return left.value < right.value;
  }
};

template<typename K>
class Cell_refine_queue
{
public:
  typedef CGAL::Anisotropic_mesh_3::Refine_cell<K>                  Refine_cell;
  typedef CGAL::Anisotropic_mesh_3::Refine_cell_comparer<K>         Refine_cell_comparer;
  typedef typename std::priority_queue<Refine_cell,
    std::vector<Refine_cell>, Refine_cell_comparer>                 Refine_cell_queue;
  typedef typename Refine_cell::Star_handle                         Star_handle;
  typedef typename Refine_cell::Cell_handle                         Cell_handle;
  typedef typename K::FT                                            FT;
  typedef typename K::Point_3                                       Point_3;

public:
  static const int nb_queues = 6;
  static const int encroachment_queue = 0;
  static const int over_distortion_queue = 1;
  static const int over_circumradius_queue = 2;
  static const int start_pick_valid = 3;
  static const int bad_shape_queue = 3;
  static const int sliver_queue = 4;
  static const int inconsistent_queue = 5;

private:
  Refine_cell_queue *queues[6];
public:
  std::list<std::pair<Point_3, Point_3> > edge_encroachments;
  Refine_cell_queue encroachments;
  Refine_cell_queue over_distortions;
  Refine_cell_queue over_circumradii;
  Refine_cell_queue bad_shapes;
  Refine_cell_queue slivers;
  Refine_cell_queue inconsistents;

#define DEFINE_PUSH(func_name, id) \
  void func_name(Star_handle star, Cell_handle cell, FT value, int tag) { \
    queues[id]->push(Refine_cell(star, cell, value, tag)); }

  DEFINE_PUSH(push_encroachment, 0)
  DEFINE_PUSH(push_over_distortion, 1)
  DEFINE_PUSH(push_over_circumradius, 2)
  DEFINE_PUSH(push_bad_shape, 3)
  DEFINE_PUSH(push_sliver, 4)
  DEFINE_PUSH(push_inconsistent, 5)

#undef DEFINE_PUSH

public:
  bool need_picking_valid(int queue_type) const
  {
    return (queue_type >= start_pick_valid);
  }

  bool empty() const
  {
    for (int i = 0; i < nb_queues; i++)
      if (!queues[i]->empty())
        return false;
    return true;
  }

  unsigned int size() const
  {
    unsigned int nb = 0;
    for (int i = 0; i < nb_queues; i++)
      nb += queues[i]->size();
    return nb;
  }

  Refine_cell top(int &queue_type) const
  {
    for (int i = 0; i < nb_queues; i++)
      if (!queues[i]->empty())
      {
        queue_type = i;
        return queues[i]->top();
      }
    return Refine_cell();
  }

  bool top(Refine_cell &cell, int &queue_type) const
  {
    for (int i = 0; i < nb_queues; i++)
      if (!queues[i]->empty())
      {
        queue_type = i;
        cell = queues[i]->top();
        return true;
      }
    return false;
  }

  bool pop()
  {
    for (int i = 0; i < nb_queues; i++)
      if (!queues[i]->empty())
      {
        queues[i]->pop();
        return true;
      }
    return false;
  }

  void clear()
  {
    for(int i = 0; i < nb_queues; i++)
      while(!queues[i]->empty())
        queues[i]->pop();
  }

public:
  void print() const
  {
    std::cout << "CQueue : ( ";
    for(int i = 0; i < nb_queues; i++)
      std::cout << queues[i]->size() <<" ";
    std::cout << ") " << std::endl;
  }

  std::size_t count() const
  {
    std::size_t count = 0;
    for(int i=0; i< nb_queues; ++i)
      count += queues[i]->size();
    return count;
  }

  //temp hack to access the priority queues underlying container
  template <class T, class S, class C>
  S& Container(std::priority_queue<T, S, C>& q)
  {
    struct HackedQueue : private std::priority_queue<T, S, C>
    {
      static S& Container(std::priority_queue<T, S, C>& q)
      {
        return q.*&HackedQueue::c;
      }
    };
    return HackedQueue::Container(q);
  }

  bool is_cell_in(Star_handle star, Cell_handle cell, FT value, int tag,
                  int queue_id)
  {
    Refine_cell rc(star, cell, value, tag);
    Refine_cell_queue* qi = queues[queue_id];
    std::vector<Refine_cell>& vi = Container(*qi);
    typename std::vector<Refine_cell>::iterator vit = vi.begin();
    typename std::vector<Refine_cell>::iterator viend = vi.end();
    for(; vit!=viend; ++vit)
      if(*vit == rc)
        return true;
    return false;
  }

  void print_queue(int queue_id)
  {
    Refine_cell_queue* qi = queues[queue_id];
    std::vector<Refine_cell>& vi = Container(*qi);
    typename std::vector<Refine_cell>::iterator vit = vi.begin();
    typename std::vector<Refine_cell>::iterator viend = vi.end();
    for(; vit!=viend; ++vit)
    {
      Refine_cell rci = *vit;
      std::cout << rci.star->index_in_star_set() << " || ";
      std::cout << rci.vertices[0] << " " << rci.vertices[1] << " " << rci.vertices[2] << " " << rci.vertices[3];
      std::cout << " || " << rci.value << std::endl;
    }
  }

public:
  Cell_refine_queue() :
    edge_encroachments(),
    encroachments(Refine_cell_comparer()),
    over_distortions(Refine_cell_comparer()),
    over_circumradii(Refine_cell_comparer()),
    bad_shapes(Refine_cell_comparer()),
    slivers(Refine_cell_comparer()),
    inconsistents(Refine_cell_comparer())
  {
    queues[encroachment_queue] = &encroachments;
    queues[over_distortion_queue] = &over_distortions;
    queues[over_circumradius_queue] = &over_circumradii;
    queues[bad_shape_queue] = &bad_shapes;
    queues[sliver_queue] = &slivers;
    queues[inconsistent_queue] = &inconsistents;
  }

  ~Cell_refine_queue() { }
};

}
}

#endif // CGAL_ANISOTROPIC_MESH_3_CELL_REFINE_QUEUE_H

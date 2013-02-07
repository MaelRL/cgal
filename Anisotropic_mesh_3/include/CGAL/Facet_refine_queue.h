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

#ifndef CGAL_ANISOTROPIC_MESH_3_REFINE_QUEUE_SURFACE_H
#define CGAL_ANISOTROPIC_MESH_3_REFINE_QUEUE_SURFACE_H

#include <iostream>
#include <fstream>
#include <utility>

#include <list>
#include <queue>
#include <vector>
#include <iostream>
#include <algorithm>
#include <CGAL/Stretched_Delaunay_3.h>

namespace CGAL {
  namespace Anisotropic_mesh_3 {

    template<typename K, typename KExact = K>
    class Refine_facet 
    {
    public:
      typedef Stretched_Delaunay_3<K, KExact> Star;
      typedef Star* Star_handle;
      typedef typename Star::Facet Facet;
      typedef typename K::FT FT;
    public:
      Star_handle star;
      int vertices[3]; // the global ids
      FT value;
    public:
      Refine_facet() : star(NULL), value(0) { }
      Refine_facet(Star_handle star_, const Facet &facet_, FT value_) :
        star(star_), 
        value(value_) 
      { 
        vertices[0] = facet_.first->vertex((facet_.second + 1) % 4)->info();
        vertices[1] = facet_.first->vertex((facet_.second + 2) % 4)->info();
        vertices[2] = facet_.first->vertex((facet_.second + 3) % 4)->info();
      }
      ~Refine_facet() { }
    };

    template<typename K, typename KExact = K>
    class Refine_facet_comparer 
    {
    public:
      Refine_facet_comparer() { };
      bool operator() (Refine_facet<K, KExact> &left, Refine_facet<K, KExact> &right) const 
      {
        return left.value < right.value;
      }
    };

    template<typename K, typename KExact = K>
    class Facet_refine_queue 
    {
    public:
      typedef CGAL::Anisotropic_mesh_3::Refine_facet<K, KExact>	Refine_facet;
      typedef CGAL::Anisotropic_mesh_3::Refine_facet_comparer<K, KExact> Refine_facet_comparer;
      typedef std::priority_queue<Refine_facet, 
        std::vector<Refine_facet>, Refine_facet_comparer> Refine_facet_queue;
      typedef typename Refine_facet::Star_handle Star_handle;
      typedef typename Refine_facet::Facet Facet;
      typedef typename K::FT FT;
      typedef typename K::Point_3 Point_3;

    public:
      static const int nb_queues = 6;
      static const int encroachment_queue = 0;
      static const int over_distortion_queue = 1;
      static const int over_circumradius_queue = 2;
      static const int bad_shape_queue = 3;
      static const int over_approximation_queue = 4;
      static const int start_pick_valid = 5;
      static const int inconsistent_queue = 5;
    private:
      Refine_facet_queue *queues[6];
    public:
      std::list<std::pair<Point_3, Point_3> > edge_encroachments;
      Refine_facet_queue encroachments;
      Refine_facet_queue over_distortions;
      Refine_facet_queue over_circumradii;
      Refine_facet_queue over_approximation;
      Refine_facet_queue bad_shapes;
      Refine_facet_queue slivers;
      Refine_facet_queue inconsistents;

#define DEFINE_PUSH(func_name, id) \
  void func_name(Star_handle star, const Facet &facet, FT value) { \
  queues[id]->push(Refine_facet(star, facet, value)); }

      DEFINE_PUSH(push_encroachment, encroachment_queue)
      DEFINE_PUSH(push_over_distortion, over_distortion_queue)
      DEFINE_PUSH(push_over_circumradius, over_circumradius_queue)
      DEFINE_PUSH(push_bad_approximation, over_approximation_queue)
      DEFINE_PUSH(push_bad_shape, bad_shape_queue)
      DEFINE_PUSH(push_inconsistent, inconsistent_queue)

#undef DEFINE_PUSH

    public:
      bool need_picking_valid(int queue_type) 
      {
        return (queue_type >= start_pick_valid);
      }
      bool empty() 
      {
        for (int i = 0; i < nb_queues; i++)
          if (!queues[i]->empty())
            return false;
        return true;
      }
      unsigned int size()
      {
        unsigned int nb = 0;
        for (int i = 0; i < nb_queues; i++)
          nb += queues[i]->size();
        return nb;
      }
      Refine_facet top(int &queue_type) 
      {
        for (int i = 0; i < nb_queues; i++)
          if (!queues[i]->empty()) 
          {
            queue_type = i;
            return queues[i]->top();
          }
          return Refine_facet();
      }
      bool top(Refine_facet &facet, int &queue_type) 
      {
        for (int i = 0; i < nb_queues; i++)
          if (!queues[i]->empty()) 
          {
            queue_type = i;
            facet = queues[i]->top();
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
      void print()
      {
        std::cout << "Queue : ( ";
        for(int i = 1; i < nb_queues; i++)
          std::cout << queues[i]->size() <<" ";
        std::cout << ") " << std::endl;
      }

    public:
      Facet_refine_queue() :
        edge_encroachments(),
        encroachments(Refine_facet_comparer()),
        over_distortions(Refine_facet_comparer()),
        over_circumradii(Refine_facet_comparer()),
        over_approximation(Refine_facet_comparer()),
        bad_shapes(Refine_facet_comparer()),
        inconsistents(Refine_facet_comparer())
      {
        queues[encroachment_queue] = &encroachments;
        queues[over_distortion_queue] = &over_distortions;
        queues[over_circumradius_queue] = &over_circumradii;
        queues[over_approximation_queue] = &over_approximation;
        queues[bad_shape_queue] = &bad_shapes;
        queues[inconsistent_queue] = &inconsistents;
      }
      ~Facet_refine_queue() {
      }
    };

  }
}

#endif // CGAL_ANISOTROPIC_MESH_3_REFINE_QUEUE_SURFACE_H

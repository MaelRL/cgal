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

#ifndef CGAL_ANISOTROPIC_MESH_3_STAR_SET_3_H
#define CGAL_ANISOTROPIC_MESH_3_STAR_SET_3_H

#include <CGAL/Stretched_Delaunay_3.h>
#include <CGAL/Constrain_surface_3.h>
#include <CGAL/Metric_field.h>
#include <CGAL/Criteria.h>
#include <CGAL/Cell_refine_queue.h>
#include <CGAL/Output_cells.h>
#include <iostream>
#include <fstream>
#include <utility>
#include "stdio.h"

#include <CGAL/helpers/statistics_helper.h>

#ifdef  CGAL_ANISOTROPIC_MESH_3_DEBUG_INFO
#define DEBUG_OUTPUT(x)	std::cout << x
#else
#define DEBUG_OUTPUT(x)
#endif

namespace CGAL{
  namespace Anisotropic_mesh_3{

    template<typename K>
    class Cell_star_set_3 {
    public:
      typedef typename Cell_star_set_3<K>         Self;
      typedef typename K::FT			  FT;
      typedef Stretched_delauney_3<K>		  Star;
      typedef Star*				  Star_handle;
      typedef std::vector<Star_handle>	          Star_vector;
      typedef typename Star_vector::iterator	  Star_iterator;
      typedef Constrain_surface_3<K>		  Constrain_surface;
      typedef Metric_field<K>			  Metric_field;
      typedef typename Metric_field::Metric	  Metric;
      typedef typename K::Point_3		  Point_3;
      typedef Criteria_base<K>		  Criteria;
      typedef typename Star::Vertex_handle	  Vertex_handle;
      typedef typename Star::Cell_handle	  Cell_handle;
      typedef typename Star::Facet_handle	  Facet_handle;
      typedef typename Star::Cell_handle_handle Cell_handle_handle;
      typedef typename Star::Facet		  Facet;
      typedef typename Cell_refine_queue<K>	  Refine_queue;
      typedef typename Refine_cell<K>		  Refine_cell;

      typedef typename std::pair<unsigned, int>	Global_triangulation_vertex_info;
      typedef typename CGAL::Triangulation_data_structure_3<
        CGAL::Triangulation_vertex_base_with_info_3<Global_triangulation_vertex_info, K> >
        Global_triangulation_ds;
      typedef typename CGAL::Delaunay_triangulation_3<K, Global_triangulation_ds>
        Global_triangulation;
      typedef typename Global_triangulation_ds::Vertex_handle	Global_triangulation_vertex_handle;

    public:
      Constrain_surface* m_pConstrain;
      const Metric_field &metric_field;
      const Criteria &criteria;
      Star_vector m_stars;
      Global_triangulation global_triangulation;
      std::vector<Point_3> m_points;
      std::vector<int> point_rounds;
      Refine_queue refine_queue;

      // the following global variables are only for debugging or benchmark
      int current_round;
      int pick_count;
      int vertex_with_picking_count;
      int vertex_without_picking_count;
      int intersection_count;
      time_t start_time;

    public:
      bool empty() const
      {
        return m_stars.empty();
      }
    public:
      Constrain_surface* constrain_surface() { return m_pConstrain; }
      void set_constrain_surface(Constrain_surface* pcs){ m_pConstrain = pcs; }

    public:
      bool is_valid_point(const Point_3 &p, const FT radius_bound, Star_handle &star)
      {
        int new_vertex_id = (int)m_points.size();
        Global_triangulation_vertex_handle gh = global_triangulation.insert(p);
        if (gh == Global_triangulation_vertex_handle())
          return false;
        else
          gh->info() = Global_triangulation_vertex_info((unsigned)new_vertex_id, 0);
        star = create_star(p, new_vertex_id, gh);
        global_triangulation.remove(gh);


        // only check adjacent stars for sliver
        Star::Vertex_handle_handle nsi = star->begin_neighboring_vertices();
        Star::Vertex_handle_handle nsiend = star->end_neighboring_vertices();
        for (; nsi != nsiend; nsi++) {
          if (star->is_infinite(*nsi))
            continue;
          if (m_stars[(*nsi)->info()]->can_form_sliver(p, radius_bound)) {
            DEBUG_OUTPUT("S");
            delete star;
            star = NULL;
            return false;
          }
        }

        // check conflicts
        bool is_cospherical = false;
        Cell_handle_handle ci = star->begin_finite_star_cells();
        Cell_handle_handle ciend = star->end_finite_star_cells();
        for (; ci != ciend; ci++) {
          if (!star->is_inside(*ci))
            continue;			

          bool is_in_bound = (star->criteria().compute_squared_circumradius(
            (*ci)->vertex(0)->point(), (*ci)->vertex(1)->point(), 
            (*ci)->vertex(2)->point(), (*ci)->vertex(3)->point()) < radius_bound);
          for (int i = 0; i < 4; i++) {
            Star::Vertex_handle v = (*ci)->vertex(i);
            if (v == star->center())
              continue;
            Star_handle foreign_star = m_stars[v->info()];
            if (foreign_star->metric().compute_distortion(star->metric()) > sqrt(criteria.distortion))
              continue;
            // if there is a point in foreign star falls into the circum-sphere of the 
            // four points under the metric of the foreign star, it is inconsistent
            Star::Vertex_handle_handle fvi = foreign_star->begin_neighboring_vertices();
            Star::Vertex_handle_handle fviend = foreign_star->end_neighboring_vertices();
            Point_3 foreign_cc = foreign_star->traits()->construct_circumcenter_3_object()(
              (*ci)->vertex(0)->point(), (*ci)->vertex(1)->point(), 
              (*ci)->vertex(2)->point(), (*ci)->vertex(3)->point());
            Star::Traits::Compute_squared_distance_3 o = foreign_star->traits()->compute_squared_distance_3_object();
            FT squared_foreign_cr = o(foreign_cc, (*ci)->vertex(0)->point());
            for (; fvi != fviend; fvi++) {
              if (foreign_star->is_infinite(*fvi)) continue;
              if ((*fvi)->info() == (*ci)->vertex(0)->info()) continue;
              if ((*fvi)->info() == (*ci)->vertex(1)->info()) continue;
              if ((*fvi)->info() == (*ci)->vertex(2)->info()) continue;
              if ((*fvi)->info() == (*ci)->vertex(3)->info()) continue;
              if (o(foreign_cc, (*fvi)->point()) < squared_foreign_cr) {
                if (!is_in_bound) {
                  for (int k = 0; k < 4; k++) {
                    if (star->criteria().compute_volume(
                      (*ci)->vertex((k + 0) % 4)->point(),
                      (*ci)->vertex((k + 1) % 4)->point(),
                      (*ci)->vertex((k + 2) % 4)->point(),
                      (*fvi)->point()) < 1e-6) {
                        continue;
                    }
                    if (star->criteria().compute_squared_circumradius(
                      (*ci)->vertex((k + 0) % 4)->point(),
                      (*ci)->vertex((k + 1) % 4)->point(),
                      (*ci)->vertex((k + 2) % 4)->point(),
                      (*fvi)->point()) < radius_bound) {
                        is_in_bound = true;
                        break;
                    }
                  }
                }
                if (!is_in_bound)
                  continue;
                is_cospherical = true;
                break;
              }
            }
            if (is_cospherical)
              break;
          }
          if (is_cospherical)
            break;
        }
        if (is_cospherical) {
          delete star;
          DEBUG_OUTPUT("C");
          star = NULL;
          return false;
        }

        return true;
      }

      Point_3 pick_valid(const Star_handle star, const Cell_handle &cell, Star_handle &retstar) 
      {
        Point_3 circumcenter = cell->circumcenter(*(star->traits()));
        FT squared_circum_radius = star->traits()->compute_squared_distance_3_object()
          (circumcenter, cell->vertex(0)->point());
        FT radius_bound = squared_circum_radius *  criteria.beta * criteria.beta;
        FT circum_radius = sqrt(squared_circum_radius) * criteria.delta;
        Star::Traits::Compute_random_point_3 random = 
          star->traits()->compute_random_point_3_object();

        int tried_times = 0;
        while (true) {
          pick_count++;
          Point_3 p = random(circumcenter, 
            ((FT)tried_times / criteria.max_times_to_try_in_picking_region) * circum_radius);

          if (is_valid_point(p, radius_bound, retstar))
            return p;
          if ((tried_times++) > criteria.max_times_to_try_in_picking_region) {
            DEBUG_OUTPUT("X\n");
            retstar = NULL;
            return circumcenter;
          }
        }
      }

      Point_3 pick_valid(const Star_handle star, const Facet &facet, Star_handle &retstar) 
      {
        Point_3 circumcenter = star->compute_circumcenter(facet);
        Point_3 cell_circumcenter = facet.first->circumcenter(*(star->traits()));
        FT squared_circum_radius = star->traits()->compute_squared_distance_3_object()
          (circumcenter, facet.first->vertex((facet.second + 1) % 4)->point());
        FT radius_bound = squared_circum_radius * criteria.beta * criteria.beta;
        FT circum_radius = sqrt(squared_circum_radius) * criteria.delta;
        Star::Traits::Compute_random_point_3 random = 
          star->traits()->compute_random_point_3_object();

        int tried_times = 0;
        while (true) {
          pick_count++;
          intersection_count++;
          Point_3 p = star->compute_exact_dual_intersection(random(circumcenter, 
            cell_circumcenter, ((FT)tried_times / criteria.max_times_to_try_in_picking_region) * circum_radius), facet);

          if (is_valid_point(p, radius_bound, retstar))
            return p;
          if ((tried_times++) > criteria.max_times_to_try_in_picking_region) {
            DEBUG_OUTPUT("X\n");
            retstar = NULL;
            return compute_insert_or_snap_point(star, facet);
          }
        }
      }

      Point_3 compute_insert_or_snap_point(const Star_handle star, const typename Constrain_surface::EdgeIterator &edge)
      {
        Point_3 p = edge->first;
        Point_3 q = edge->second;
        Point_3 c((p.x() + q.x()) / 2.0, (p.y() + q.y()) / 2.0, (p.z() + q.z()) / 2.0);
        m_pConstrain->edge_split(edge, c);
        return c;
      }

      Point_3 compute_insert_or_snap_point(const Star_handle star, const Facet &facet) 
      {
        intersection_count++;
        Point_3 p = star->compute_exact_dual_intersection(facet);

        // test encroachment
#ifdef CHECK_EDGE_ENCROACHMENT
        Star_iterator si = m_stars.begin();
        Star_iterator siend = m_stars.end();
        Constrain_surface::EdgeIterator encroached_edge;
        for (; si != siend; si++) {
          if ((*si)->is_encroached(p, encroached_edge))
            return compute_insert_or_snap_point(*si, encroached_edge);
        }
#endif
        return p;
      }

      Point_3 compute_insert_or_snap_point(const Star_handle star, const Cell_handle cell)
      {
        Point_3 p = cell->circumcenter(*(star->traits()));
        // test encroachment

#ifdef CHECK_EDGE_ENCROACHMENT
        {
          Star_iterator si = m_stars.begin();
          Star_iterator siend = m_stars.end();
          Constrain_surface::EdgeIterator encroached_edge;
          for (; si != siend; si++) {
            if ((*si)->is_encroached(p, encroached_edge))
              return compute_insert_or_snap_point(*si, encroached_edge);
          }
        }
#endif
#ifdef CHECK_FACE_ENCROACHMENT
        {
#ifdef ONLY_CONSIDER_NEIGHBORING_ENCROACHMENT
          Facet encroached_facet;
          if (m_stars.size() < 1000) {
            Star_iterator si = m_stars.begin();
            Star_iterator siend = m_stars.end();
            Facet encroached_facet;
            for (; si != siend; si++) {
              if ((*si)->is_encroached(p, encroached_facet))
                return compute_insert_or_snap_point(*si, encroached_facet);
            }
          } else {
            if (star->is_encroached(p, encroached_facet))
              return compute_insert_or_snap_point(star, encroached_facet);

            Star::Vertex_handle_handle nvi = star->begin_neighboring_vertices();
            Star::Vertex_handle_handle nviend = star->end_neighboring_vertices();
            for (; nvi != nviend; nvi++) {
              if (star->is_infinite(*nvi))
                continue;
              Star_handle css = m_stars[(*nvi)->info()];
              if (css->is_encroached(p, encroached_facet))
                return compute_insert_or_snap_point(css, encroached_facet);
            }
          }
#else
          Star_iterator si = m_stars.begin();
          Star_iterator siend = m_stars.end();
          Facet encroached_facet;
          for (; si != siend; si++) {
            if ((*si)->is_encroached(p, encroached_facet))
              return compute_insert_or_snap_point(*si, encroached_facet);
          }
#endif
        }
#endif
        return p;
      }

      Point_3 compute_insert_or_snap_valid_point(const Star_handle star, const Facet &facet, Star_handle &retstar)
      {
        Point_3 p = pick_valid(star, facet, retstar);
        // test encroachment

#ifdef CHECK_EDGE_ENCROACHMENT
        Star_iterator si = m_stars.begin();
        Star_iterator siend = m_stars.end();
        Constrain_surface::EdgeIterator encroached_edge;
        for (; si != siend; si++) {
          if ((*si)->is_encroached(p, encroached_edge)) {
            delete retstar;
            retstar = NULL;
            return compute_insert_or_snap_point(*si, encroached_edge);
          }
        }
#endif
        return p;
      }

      Point_3 compute_insert_or_snap_valid_point(const Star_handle star, const Cell_handle cell, Star_handle &retstar) 
      {
        Point_3 p = pick_valid(star, cell, retstar);
        // test encroachment

#ifdef CHECK_EDGE_ENCROACHMENT
        {
          Star_iterator si = m_stars.begin();
          Star_iterator siend = m_stars.end();
          Constrain_surface::EdgeIterator encroached_edge;
          for (; si != siend; si++) {
            if ((*si)->is_encroached(p, encroached_edge)) {
              delete retstar;
              retstar = NULL;
              return compute_insert_or_snap_point(*si, encroached_edge);
            }
          }
        }
#endif
#ifdef CHECK_FACE_ENCROACHMENT
        {

#ifdef ONLY_CONSIDER_NEIGHBORING_ENCROACHMENT
          Facet encroached_facet;
          {
            if (star->is_encroached(p, encroached_facet)) {
              delete retstar;
              retstar = NULL;
              return compute_insert_or_snap_valid_point(star, encroached_facet, retstar);
            }

            Star::Vertex_handle_handle nvi = star->begin_neighboring_vertices();
            Star::Vertex_handle_handle nviend = star->end_neighboring_vertices();
            for (; nvi != nviend; nvi++) {
              if (star->is_infinite(*nvi))
                continue;
              Star_handle css = m_stars[(*nvi)->info()];
              if (css->is_encroached(p, encroached_facet)) {
                delete retstar;
                retstar = NULL;
                return compute_insert_or_snap_valid_point(css, encroached_facet, retstar);
              }
            }
          }
#else
          Star_iterator si = m_stars.begin();
          Star_iterator siend = m_stars.end();
          Facet encroached_facet;
          for (; si != siend; si++) {
            if ((*si)->is_encroached(p, encroached_facet)) {
              delete retstar;
              retstar = NULL;
              return compute_insert_or_snap_valid_point(*si, encroached_facet, retstar);
            }
          }
#endif
        }
#endif
        return p;
      }

      int insert(const Point_3 &point, 
                 Star_vector &modified_stars, 
                 Global_triangulation_vertex_handle &global_handle) 
      {
        //global triangulation
        int this_id = (int)m_points.size();
        Global_triangulation_vertex_handle gh = global_triangulation.insert(point);
        if (gh == Global_triangulation_vertex_handle())
          return -1;
        else
          gh->info() = Global_triangulation_vertex_info((unsigned)this_id, 0);

        modified_stars.clear();

        // insert p into all stars
#ifdef	INSERT_POINT_ENHANCEMENT
        if (m_stars.size() > 1000) 
        {
          current_round++;
          std::list<Global_triangulation_vertex_handle> queue;
          std::back_insert_iterator<std::list<Global_triangulation_vertex_handle> > 
            queue_insertor(queue);
          gh->info().second = current_round;
          global_triangulation.adjacent_vertices(gh, queue_insertor);
          while (!queue.empty()) 
          {
            Global_triangulation_vertex_handle h = queue.front();
            if (h->info().second == current_round)// visited
            { 
              queue.pop_front();
              continue;
            }
            unsigned vid = h->info().first;
            h->info().second = current_round;
            if (m_stars[vid]->insert_to_star(point, (int)m_points.size()) != Star::Vertex_handle()) 
            {
              global_triangulation.adjacent_vertices(h, queue_insertor);
              modified_stars.push_back(m_stars[vid]);
            }
            queue.pop_front();
          }
        } 
        else 
          naive_insertion(point, modified_stars);

#else
        naive_insertion(point, modified_stars);

#endif
        m_points.push_back(point);
        point_rounds.push_back(0);
        global_handle = gh;
        return this_id;
      }

      void naive_insertion(const Point_3 &point, 
                           Star_vector &modified_stars)
      {
        Star_iterator si = m_stars.begin();
        Star_iterator siend = m_stars.end();
        for (; si != siend; si++) 
        {
          if ((*si)->insert_to_star(point, (int)m_points.size()) != Star::Vertex_handle())
            modified_stars.push_back(*si);
        }
      }

      void check_conflicts(Star_vector &modified_stars, int relative_point = -1) 
      {
        Star_iterator si = modified_stars.begin();
        Star_iterator siend = modified_stars.end();
        for (; si != siend; si++) {
          Star_handle star = *si;

          // check encroachment
          Facet_handle bi = star->begin_boundary_facets();
          Facet_handle biend = star->end_boundary_facets();
          for (; bi != biend; bi++) 
          {
            if ((relative_point >= 0) && (bi->first->vertex(bi->second)->info() != relative_point))
              continue;
#ifdef CHECK_FACE_ENCROACHMENT
            if (star->is_facet_encroached(bi->first->vertex(bi->second)->point(), *bi))
              refine_queue.push_encroachment(star, bi->first, star->compute_volume(bi->first), bi->second);
#endif
          }

          Cell_handle_handle ci = star->begin_finite_star_cells();
          Cell_handle_handle ciend = star->end_finite_star_cells();
          for (; ci != ciend; ci++) 
          {
            Cell_handle c = *ci;
            if (!star->is_inside(c))
              continue;
            assert(!star->is_infinite(c));

            if (relative_point >= 0) 
            {
              bool relative = false;
              for (int i = 0; i < 4; i++)
                if (relative_point == c->vertex(i)->info()) 
                {
                  relative = true;
                  break;
                }
                if (!relative)
                  continue;
            }

            // over distortion
            for (int i = 0; i < 4; i++)
              for (int j = i + 1; j < 4; j++) {
                FT over_distoration = 					
                  m_stars[c->vertex(i)->info()]->metric().compute_distortion(
                  m_stars[c->vertex(j)->info()]->metric()) - criteria.distortion;
                if (over_distoration > 0) {
                  refine_queue.push_over_distortion(star, c, over_distoration, 0);
                  goto next_cell;
                }
              }

              // too big
              FT over_circumradius =
                star->compute_circumradius_overflow(c);
              if (over_circumradius > 0) {
                refine_queue.push_over_circumradius(star, c, over_circumradius, 0);
                goto next_cell;
              }

              // bad shape
              FT over_radius_edge_ratio =
                star->compute_radius_edge_ratio_overflow(c);
              if (over_radius_edge_ratio > 0) {
                refine_queue.push_bad_shape(star, c, over_radius_edge_ratio, 0);
                goto next_cell;
              }

              // sliverity
              FT over_sliverity =
                star->compute_sliverity_overflow(c);
              if (over_sliverity > 0) {
                refine_queue.push_sliver(star, c, over_sliverity, 0);
                goto next_cell;
              }

              // consistency
              for (int i = 0; i < 4; i++) {
                Vertex_handle v = c->vertex(i);
                if (star->center() == v)
                  continue;
                if (!m_stars[v->info()]->has_cell(c)) {
                  refine_queue.push_inconsistent(star, c, star->compute_volume(c), 0);
                  goto next_cell;
                }
              }
next_cell:	;
          } // for cells
        } // for stars
      }

      void initialize_stars() 
      {
        DEBUG_OUTPUT("Initializing stars...");
        Constrain_surface::Pointset initial_points = m_pConstrain->initial_points();
        Constrain_surface::Pointset::iterator pi = initial_points.begin();
        Constrain_surface::Pointset::iterator pend = initial_points.end();
        for (int i = 0; pi != pend; pi++, i++) 
        {
          DEBUG_OUTPUT(".");
          Point_3 p = *pi;
          Global_triangulation_vertex_handle h = global_triangulation.insert(p);
          if (h == Global_triangulation_vertex_handle())
            continue;
          else
            h->info() = std::pair<unsigned, int>(i, 0);
        
          Star_handle star = new Star(criteria, p, i, metric_field.compute_metric(p), m_pConstrain);
          Constrain_surface::Pointset::iterator qj = initial_points.begin();
          Constrain_surface::Pointset::iterator qend = initial_points.end();
          for (int j = 0; qj != qend; qj++, j++)
            if (i != j)
              star->insert_to_star(*qj, j);
          m_stars.push_back(star);
          m_points.push_back(p);
          point_rounds.push_back(0);

        }
        DEBUG_OUTPUT("\nInitializing conflicts...\n");
        check_conflicts(m_stars);
        DEBUG_OUTPUT(" [ok]" << std::endl);
      }

      bool next_refine_cell(Refine_cell &refine_cell, Cell_handle &cell, bool &need_picking_valid, int &queue_type)
      {
        while (true) 
        {
          if (!refine_queue.top(refine_cell, queue_type))
            return false;
          refine_queue.pop();
          if (refine_cell.star->has_cell_ref(refine_cell.vertices, cell)) 
          {
            if (!refine_cell.star->is_inside(cell))
              continue;
            need_picking_valid = refine_queue.need_picking_valid(queue_type);
            return true;
          }
        }
      }

      Star_handle create_star(const Point_3 &p, int pid, Global_triangulation_vertex_handle &global_handle) 
      {
        Star_handle star = new Star(criteria, p, pid, metric_field.compute_metric(p), m_pConstrain);

#ifdef	CREATE_STAR_ENHANCEMENT
        if (m_stars.size() > 1000) {
          current_round++;
          std::list<Global_triangulation_vertex_handle> queue;
          std::back_insert_iterator<std::list<Global_triangulation_vertex_handle> > 
            queue_insertor(queue);
          global_handle->info().second = current_round;
          global_triangulation.adjacent_vertices(global_handle, queue_insertor);
          while (!queue.empty()) {
            Global_triangulation_vertex_handle h = queue.front();
            if (h->info().second == current_round)	{ // visited
              queue.pop_front();
              continue;
            }
            unsigned vid = h->info().first;
            h->info().second = current_round;
            if (star->insert_to_star(m_points[vid], vid) != Star::Vertex_handle())
              global_triangulation.adjacent_vertices(h, queue_insertor);
            queue.pop_front();
          }
        } else {
          for (int i = 0; i < pid; i++)
            star->insert_to_star(m_points[i], i);
        }
#else
        for (int i = 0; i < pid; i++)
          star->insert_to_star(m_points[i], i);
#endif
        return star;
      }

      bool refine()
      {
        int this_id = (int)m_points.size();

        Point_3 p;
        int queue_type = 0;
        Refine_cell cell;
        bool need_picking_valid;
        Cell_handle c;
        Star_handle star = NULL;		

repick_cell:
        {
          if (!next_refine_cell(cell, c, need_picking_valid, queue_type))
            return false;
#ifdef ONLY_CONSIDER_SURFACE
          if (!cell.star->is_boundary_star())
            return true;
#endif
          if (queue_type == 0) { // encroachment
            int tag = 0;
            for (int i = 0; i < 4; i++)
              if (c->vertex(i)->info() == cell.vertices[0]) {
                tag = i;
                break;
              }
              if (!cell.star->has_facet(Facet(c, tag))) {
                DEBUG_OUTPUT("[repick]");
                goto repick_cell;
              }
              p = compute_insert_or_snap_point(cell.star, Facet(c, tag));
              vertex_without_picking_count++;
          } else {
            if (need_picking_valid) {
              p = compute_insert_or_snap_valid_point(cell.star, c, star);
              vertex_with_picking_count++;
            } else {
              p = compute_insert_or_snap_point(cell.star, c);
              vertex_without_picking_count++;
            }
          }


#ifdef CHECK_EXISTING_VERTICES
          // check already exists
          for (int i = (int)m_points.size() - 1; i >= 0; i--) {
            Point_3 q = m_points[i];
            if (fabs(p.x() - q.x()) +
              fabs(p.y() - q.y()) +
              fabs(p.z() - q.z()) < 1e-4)
            {
              std::cout << "\nsame with " << i << std::endl;
              std::cout << p << " vs. " << q << std::endl;
              return true;
            }
          }
#endif

          DEBUG_OUTPUT(this_id << " refine: " << p << " [" << 
            c->vertex(0)->info() << "-" <<
            c->vertex(1)->info() << "-" <<
            c->vertex(2)->info() << "-" <<
            c->vertex(3)->info() << "] " <<
            "t:" << queue_type << " i:" << cell.vertices[0] << " s:" << cell.star->center()->info() <<
            std::endl);
        }

        Star_vector affected_stars;
        Global_triangulation_vertex_handle global_handle;
        int pid = insert(p, affected_stars, global_handle);
        if (pid < 0) {
          DEBUG_OUTPUT("[Already exists!]");
          return true;
        }

        // create star
        if (star == NULL)
          star = create_star(p, pid, global_handle);
        m_stars.push_back(star);
        check_conflicts(affected_stars, pid);
        return true;
      }

      void dump(const bool nb = false) 
      {
        std::ostringstream nbs;
        nbs << m_points.size();
        
        std::string title("dump");
        if(nb) title.append(nbs.str());
        title.append(".off");
        
        std::cout << "Write " << title << " (" << nbs.str() << " vertices)...";
        typename std::ofstream fx(title);

        fx << "OFF" << std::endl;
        fx << m_points.size() << " 0 0 " << std::endl;
        for (int i = 0; i < (int)m_points.size(); i++)
          fx << m_points[i] << std::endl;

        /*fx << 0 << std::endl;
        for (int i = 0; i < (int)m_stars.size(); i++) 
        {
          Star_handle star = m_stars[i];
          Star::Vertex_iterator vi = star->vertices_begin();
          Star::Vertex_iterator viend = star->vertices_end();
          for (; vi != viend; vi++) 
          {
            if (star->is_infinite(vi))
              continue;
            fx << (*vi).info() << " ";
          }
          fx << -1 << std::endl;
       }*/
       std::cout << "done.\n";
       fx.close();
      }

/*      void load_dump() 
      {
        m_stars.clear();
        m_points.clear();

        typename std::ifstream fx("dump.txt");
        DEBUG_OUTPUT("Loading...");
        int dimension;
        fx >> dimension;
        int point_count;
        fx >> point_count;
        for (int i = 0; i < point_count; i++) 
        {
          Point_3 p;
          fx >> p;
          m_points.push_back(p);
          point_rounds.push_back(0);
        }
        int cell_count_zero;
        fx >> cell_count_zero;
        for (int i = 0; i < point_count; i++) 
        {
          Star_handle star = new Star(
            criteria, m_points[i], i, metric_field.compute_metric(m_points[i]), m_pConstrain);
          while (true) 
          {
            int pid;
            fx >> pid;
            if (pid < 0)
              break;
            if (pid != i)
              star->insert_to_star(m_points[pid], pid);
          }
          m_stars.push_back(star);
          std::cout << ".";
        }
        DEBUG_OUTPUT("\nInitializing conflicts...");
        check_conflicts(m_stars);
        DEBUG_OUTPUT("\n\n");
      }
      */
      void report() 
      {
        typename std::ofstream fx("report.txt");

        fx << "[Parameters]" << std::endl << std::endl;
        criteria.report(fx);

        fx << std::endl << "[Metric field]" << std::endl << std::endl;
        metric_field.report(fx);

        fx << std::endl << "[Statistics]" << std::endl << std::endl;
        fx << "elapsed time:       " << time(NULL) - start_time << " sec." << std::endl;
        fx << "number of vertices: " << m_stars.size() << std::endl;
        fx << "picking times:      " << pick_count << std::endl;
        fx << "vertex via picking: " << vertex_with_picking_count << std::endl;
        fx << "vertex non-picking: " << vertex_without_picking_count << std::endl;
        fx << "picking rate:       " << (double)vertex_with_picking_count / (double)(vertex_without_picking_count + vertex_with_picking_count) << std::endl;
        fx << "picking time rate:  " << (double)pick_count / (double)(pick_count + vertex_without_picking_count) << std::endl;
        fx << "mean picking times: " << (double)pick_count / (double)vertex_with_picking_count << std::endl;
        fx << "intersection count: " << intersection_count << std::endl;
        fx.close();
      }

      void output() 
      {
        typename std::ofstream fx("mesh.volume.cgal");
        fx << 3 << std::endl;
        fx << m_points.size() << std::endl;
        for (int i = 0; i < (int)m_points.size(); i++)
          fx << m_points[i] << std::endl;

        Output_cells output_cells;
        Star_iterator si = m_stars.begin();
        Star_iterator siend = m_stars.end();
        for (; si != siend; si++) 
        {
          Star_handle star = *si;
          Star::Cell_handle_handle ci = star->begin_finite_star_cells();
          Star::Cell_handle_handle ciend = star->end_finite_star_cells();
          for (; ci != ciend; ci++) {
            if (!star->is_inside(*ci))
              continue;
            bool consistent = true;
            for (int i = 0; i < 4; i++)
              if (!m_stars[(*ci)->vertex(i)->info()]->has_cell(*ci)) 
              {
                consistent = false;
                break;
              }
              if (consistent)
                output_cells.insert(
                (*ci)->vertex(0)->info(), (*ci)->vertex(1)->info(), 
                (*ci)->vertex(2)->info(), (*ci)->vertex(3)->info());
          }
        }
        Output_cells::Cell_handle oci = output_cells.begin();
        Output_cells::Cell_handle ociend = output_cells.end();
        fx << output_cells.size() << std::endl;
        for (; oci != ociend; oci++) 
        {
          fx << oci->vertices[0] + 1 << " " << oci->vertices[1] + 1 << " "
            << oci->vertices[2] + 1 << " " << oci->vertices[3] + 1 << std::endl;
        }
        fx.close();		
      }

      void output_medit() 
      {
        unsigned int nb_inconsistent_stars = 0;
        typename std::ofstream fx("meshvolume.cgal.mesh");
        fx << "MeshVersionFormatted 1\n";
        fx << "Dimension 3\n";

        fx << "Vertices\n";
        fx << m_points.size() << std::endl;
        for (int i = 0; i < (int)m_points.size(); i++)
          fx << m_points[i] << " " << (i+1) << std::endl; // warning : indices start at 1 in Medit

        Output_cells output_cells;
        Star_iterator si = m_stars.begin();
        Star_iterator siend = m_stars.end();
        for (; si != siend; si++) 
        {
          Star_handle star = *si;
          Star::Cell_handle_handle ci = star->begin_finite_star_cells();
          Star::Cell_handle_handle ciend = star->end_finite_star_cells();
          for (; ci != ciend; ci++) 
          {
            if (!star->is_inside(*ci))
              continue;
            bool consistent = true;
            for (int i = 0; i < 4; i++)
              if (!m_stars[(*ci)->vertex(i)->info()]->has_cell(*ci)) 
              {
                consistent = false;
                nb_inconsistent_stars++;
                break;
              }
              if (consistent)
                output_cells.insert(
                (*ci)->vertex(0)->info(), (*ci)->vertex(1)->info(), 
                (*ci)->vertex(2)->info(), (*ci)->vertex(3)->info());
          }
        }
        fx << "Tetrahedra\n";
        fx << output_cells.size() << std::endl;
        Output_cells::Cell_handle oci = output_cells.begin();
        Output_cells::Cell_handle ociend = output_cells.end();
        for (; oci != ociend; oci++) 
        {
          fx << (oci->vertices[0] + 1) << " " << (oci->vertices[1] + 1) << " "
             << (oci->vertices[2] + 1) << " " << (oci->vertices[3] + 1) << " "
             << "1" /*color*/<< std::endl;
        }
        fx.close();	
        if(nb_inconsistent_stars > 0)
          std::cout << "\nWarning : there are " << nb_inconsistent_stars << " inconsistent stars in the ouput mesh.\n";
      }

      void refine_all(int max_count = INT_MAX) 
      {
        typename std::ofstream fe("report_times.txt");
        fe.close();

        pick_count = 0;
        vertex_with_picking_count = 0;
        vertex_without_picking_count = (int)m_points.size();
        intersection_count = 0;
        start_time = time(NULL);

        while ((int)m_points.size() < max_count) 
        {
//          int before_refinement = (int)m_points.size();
          if (!refine())
            break;
          /*if (m_points.size() != before_refinement) 
          {
            if (m_points.size() % 1000 == 0)
              dump();
          }*/
        }
        report(); // prints time
        output_medit();
        dump(true);
        histogram_vertices_per_star<Self>(*this);

      }

    public:

      Cell_star_set_3(const Criteria &criteria_, 
        const Metric_field &metric_field_, 
        Constrain_surface* pconstrain_, 
        bool init_new) :
      metric_field(metric_field_), 
        m_pConstrain(pconstrain_), 
        criteria(criteria_), 
        m_stars(), 
        m_points(), 
        refine_queue(),
        global_triangulation(), 
        current_round(0), 
        point_rounds() 
      {
        if (init_new)
          initialize_stars();
        //else
        //  load_dump();
      }
      Cell_star_set_3(const Cell_star_set_3& css) :
      metric_field(css.metric_field), 
        m_pConstrain(css.m_pConstrain), 
        criteria(css.criteria), 
        m_stars(), 
        m_points(), 
        refine_queue(),
        global_triangulation(), 
        current_round(0), 
        point_rounds() 
      {}
      Cell_star_set_3() :
      metric_field(Metric_field()),
        m_pConstrain(NULL), 
        criteria(Criteria()),
        m_stars(),
        m_points(), 
        refine_queue(),
        global_triangulation(), 
        current_round(0), 
        point_rounds() 
      {}


      ~Cell_star_set_3() 
      {
        Star_iterator si = m_stars.begin();
        Star_iterator siend = m_stars.end();
        for (; si != siend; si++)
          delete (*si);
      }
    };

  }
}

#undef DEBUG_OUTPUT

#endif // CGAL_ANISOTROPIC_MESH_3_STAR_SET_3_H

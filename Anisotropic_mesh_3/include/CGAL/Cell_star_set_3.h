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

#include <CGAL/Cell_refine_queue.h>
#include <CGAL/Constrain_surface_3.h>
#include <CGAL/Criteria.h>
#include <CGAL/Metric_field.h>
#include <CGAL/Stretched_Delaunay_3.h>
#include <CGAL/Output_facets.h>
#include <CGAL/Output_cells.h>

#include <CGAL/helpers/combinatorics_helper.h>
#include <CGAL/helpers/statistics_helper.h>
#include <CGAL/helpers/metric_helper.h>

#include <CGAL/aabb_tree/aabb_tree_bbox.h>
#include <CGAL/aabb_tree/aabb_tree_bbox_primitive.h>
#include <CGAL/kd_tree/Kd_tree_for_star_set.h>

#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <vector>
#include "stdio.h"


namespace CGAL{
namespace Anisotropic_mesh_3{

template<typename K>
class Cell_star_set_3 {
public:
  typedef Cell_star_set_3<K>                        Self;

  typedef typename K::FT                            FT;
  typedef typename K::Point_3                       Point_3;

  typedef Stretched_Delaunay_3<K>                   Star;
  typedef Star*                                     Star_handle;
  typedef std::set<Star_handle>                     Star_set;
  typedef std::vector<Star_handle>                  Star_vector;
  typedef typename Star::Base                       DT;
  typedef typename Star_vector::iterator            Star_iterator;
  typedef typename Star::Vertex_handle              Vertex_handle;
  typedef typename Star::Cell_handle                Cell_handle;
  typedef typename Star::Facet                      Facet;
  typedef typename Star::Facet_handle               Facet_handle;
  typedef typename Star::Cell_handle_handle         Cell_handle_handle;
  typedef typename Star::TPoint_3                   TPoint_3;
  typedef typename Star::Vector_3                   Vector_3;
  typedef typename Star::Index                      Index;
  typedef std::set<Index>                           Index_set;

  typedef CGAL::AABB_tree_bbox<K, Star>             AABB_tree;
  typedef CGAL::AABB_bbox_primitive<Star>           AABB_primitive;

  typedef CGAL::Kd_tree_for_star_set<K, Star_handle>          Kd_tree;
  typedef typename Kd_tree::Traits                            Kd_traits;
  typedef typename Kd_tree::Box_query                         Kd_Box_query;
  typedef typename Kd_tree::key_type                          Kd_point_info;

  typedef Constrain_surface_3<K>                    Constrain_surface;
  typedef CGAL::Anisotropic_mesh_3::Metric_field<K> Metric_field;
  typedef typename Metric_field::Metric             Metric;
  typedef Criteria_base<K>                          Criteria;

  typedef typename CGAL::Anisotropic_mesh_3::Cell_refine_queue<K>             Refine_queue;
  typedef typename CGAL::Anisotropic_mesh_3::Refine_cell<K>                   Refine_cell;

public:
  const Constrain_surface* const m_pConstrain;
  const Metric_field* m_metric_field;
  const Criteria* m_criteria;
  Star_vector m_stars;
  Refine_queue m_refine_queue;

  DT m_ch_triangulation;
  mutable AABB_tree m_aabb_tree; //bboxes of stars
  Kd_tree m_kd_tree;

  // the following global variables are only for debugging or benchmark
  int pick_count;
  int vertex_with_picking_count;
  int vertex_without_picking_count;
  int intersection_count;
  time_t start_time;

public:
  const Constrain_surface* constrain_surface() { return m_pConstrain; }
  const Constrain_surface* const constrain_surface() const { return m_pConstrain; }
  const Metric_field* metric_field() const { return m_metric_field; }
  const Criteria* criteria() const { return m_criteria; }
  void set_criteria(const Criteria* criteria_) { m_criteria = criteria_; }
  void set_metric_field(const Metric_field* metric_field_) { m_metric_field = metric_field_; }

public:
  std::size_t number_of_stars() const { return m_stars.size(); }
  bool empty() const { return m_stars.empty(); }
  std::size_t size() const { return m_stars.size(); }
  bool is_infinite_vertex(int i) const { return (i == Star::infinite_vertex_index()); }

private:
  Star_handle get_star(Star_handle s) const { return s; }
  Star_handle get_star(std::size_t i) const { return m_stars[i]; }
  Star_handle get_star(int i) const         { return m_stars[i]; }
  Star_handle get_star(typename Index_set::const_iterator it) const   { return m_stars[*it]; }
  Star_handle get_star(typename Star_set::const_iterator it) const    { return *it; }
  Star_handle get_star(typename Star_vector::const_iterator it) const { return *it; }

public:
  TPoint_3 transform_to_star_point(const Point_3& p, Star_handle star) const
  {
    return star->metric().transform(p);
  }
  Point_3 transform_from_star_point(const TPoint_3& p, Star_handle star) const
  {
    return star->metric().inverse_transform(p);
  }

public:
  double duration(const time_t& start) const
  {
    return ((clock() - start + 0.) / ((double)CLOCKS_PER_SEC));
  }

  void update_aabb_tree(Star_handle star) const
  {
    m_aabb_tree.update_primitive(AABB_primitive(star));
  }

  void build_aabb_tree()
  {
    m_aabb_tree.rebuild(m_stars.begin(), m_stars.end());
  }

public:
  Point_3 compute_circumcenter(const Point_3& p0, const Point_3& p1,
                               const Point_3& p2, const Point_3& p3,
                               Star_handle here) const
  {
#ifdef ANISO_USE_CC_EXACT
    return back_from_exact(
       transform_from_star_point(here->compute_circumcenter(
              transform_to_star_point(to_exact(p0), here),
              transform_to_star_point(to_exact(p1), here),
              transform_to_star_point(to_exact(p2), here),
              transform_to_star_point(to_exact(p3), here)), here));
#else
    return transform_from_star_point(here->compute_circumcenter(
              transform_to_star_point(p0, here),
              transform_to_star_point(p1, here),
              transform_to_star_point(p2, here),
              transform_to_star_point(p3, here)), here);
#endif
  }

  bool is_inside_domain(const Point_3& p) const
  {
    return (m_pConstrain->side_of_constraint(p) == ON_POSITIVE_SIDE);
  }

public:
  void update_bboxes() const
  {
#ifdef DEBUG_UPDATE_AABB_TREE
    std::cout << "updating bboxes. tree : " << m_aabb_tree.size() << " " << m_aabb_tree.m_insertion_buffer_size() << std::endl;
#endif
    std::size_t i;
    std::size_t N = m_stars.size();
    for(i = 0; i < N; i++)
    {
      m_stars[i]->update_bbox();
      if(m_aabb_tree.m_insertion_buffer_size() == 1) // tree has just been rebuilt
        m_stars[i]->set_bbox_needs_aabb_update(false);
      if(m_stars[i]->bbox_needs_aabb_update() && i<m_aabb_tree.size())
      {
        m_stars[i]->set_bbox_needs_aabb_update(false);
        update_aabb_tree(m_stars[i]);
      }
    }

    for(i = 0; i < N; i++)
      if(m_stars[i]->bbox_needs_aabb_update() && i<m_aabb_tree.size())
        std::cout << "forgot some stars in the update" << std::endl;
  }

  template<typename OutputIterator>
  void finite_stars_in_conflict(const Point_3& p,
                                OutputIterator oit) const
  {
    update_bboxes();
#ifndef NO_USE_AABB_TREE_OF_BBOXES
    //bigger set of stars
    m_aabb_tree.all_intersected_primitives(p, oit);
#else
    //exact set
    for(unsigned int i = 0; i < m_stars.size(); i++)
    {
      Star_handle s = m_stars[i];
      Cell_handle ch;
      if(s->is_conflicted(transform_to_star_point(p, s), ch))
        *oit++ = s;
    }
#endif
  }

  template<typename OutputIterator>
  void all_stars(OutputIterator oit) const
  {
    for(std::size_t i = 0; i < m_stars.size(); i++)
      *oit++ = m_stars[i];
  }

  template<typename OutputIterator>
  bool infinite_stars_in_conflict(const Point_3& p,
                                  OutputIterator oit,
                                  const bool collect = true) const
  {
    if(m_ch_triangulation.dimension() < 3)
    {
      all_stars(oit);
      return true;
    }

    int li, lj;
    typename DT::Locate_type lt;
    typename DT::Cell_handle c = m_ch_triangulation.locate(p, lt, li, lj);

    if(lt != DT::OUTSIDE_CONVEX_HULL && lt != DT::OUTSIDE_AFFINE_HULL)
      return false;

    std::vector<typename DT::Facet> bfacets;
    m_ch_triangulation.find_conflicts(p, c, std::back_inserter(bfacets), Emptyset_iterator());

    if(!collect)
       return !bfacets.empty();

    typename std::vector<typename DT::Facet>::iterator fit = bfacets.begin();
    typename std::vector<typename DT::Facet>::iterator fitend = bfacets.end();
    for(; fit != fitend; fit++)
    {
      typename DT::Facet f = *fit;
      if(m_ch_triangulation.is_infinite(f))
        continue;

      *oit++ = m_stars[f.first->vertex((f.second + 1) % 4)->info()];
      *oit++ = m_stars[f.first->vertex((f.second + 2) % 4)->info()];
      *oit++ = m_stars[f.first->vertex((f.second + 3) % 4)->info()];
    }
    return true;
  }

public:
  template<typename StarConstIterator>
  void remove_from_stars(const Index& id,
                         StarConstIterator begin,
                         StarConstIterator end)
  {
    StarConstIterator it;
    for(it = begin; it != end; it++)
    {
      Star_handle s = get_star(it);
      if(id != s->index_in_star_set())
        s->remove(id);
    }
  }

  void create_star(const Point_3 &p,
                   int pid,
                   const Index_set& modified_stars,//by p's insertion
                   Star_handle& star,
                   const bool surface_star = true) const
  {
    //todo if(surf)
    star->reset(p, pid, m_metric_field->compute_metric(p), surface_star);

    typename Index_set::const_iterator cit = modified_stars.begin();
    typename Index_set::const_iterator citend = modified_stars.end();
    for (; cit != citend; cit++)
    {
      Star_handle si = get_star(cit);
      star->insert_to_star(si->center_point(), si->index_in_star_set(), false);
    }

    typename Star::Bbox bbox = star->bbox(); // volume bbox + surface bbox when star is not a topo_disk
    Point_3 pmin(bbox.xmin(), bbox.ymin(), bbox.zmin());
    Point_3 pmax(bbox.xmax(), bbox.ymax(), bbox.zmax());
    Kd_Box_query query(pmin, pmax, /*3=dim,*/ 0./*epsilon*/, typename Kd_tree::Star_pmap(m_stars));
    std::set<Kd_point_info> indices;
    m_kd_tree.search(std::inserter(indices, indices.end()), query);

    if(indices.size() != modified_stars.size())// because 'indices' contains 'modified_stars'
    {
      Index_set diff;
      std::set_difference(indices.begin(), indices.end(),
        modified_stars.begin(), modified_stars.end(), std::inserter(diff, diff.end()));
      typename Index_set::iterator it = diff.begin();
      while(it != diff.end())
      {
        Star_handle si = get_star(it++);
        star->insert_to_star(si->center_point(), si->index_in_star_set(), true/*conditional*/);
      }
    }
  }

  Star_handle create_star(const Point_3 &p,
                          int pid,
                          const Index_set& modified_stars)
  {
    Star_handle star = new Star(m_criteria, m_pConstrain, true /*surface star*/);
    create_star(p, pid, modified_stars, star, true /*surface star*/);
    return star;
  }

  Star_handle create_inside_star(const Point_3 &p,
                                 int pid,
                                 const Index_set& modified_stars) const //by p's insertion
  {
    Star_handle star = new Star(m_criteria, m_pConstrain, false/*surface*/);
    create_star(p, pid, modified_stars, star, false/*surface*/);
    return star;
  }

  Index perform_insertions(const Point_3& p,
                           const Index& this_id,
                           const Star_set& target_stars,
                           Index_set& modified_stars,
                           const bool conditional,
                           const bool infinite_stars = false)
  {
    Vertex_handle v_in_ch = Vertex_handle();
    typename Star_set::const_iterator it = target_stars.begin();
    typename Star_set::const_iterator itend = target_stars.end();
    for(; it != itend; it++)
    {
      Star_handle si = *it;
      Vertex_handle vi = si->insert_to_star(p, this_id, conditional);

      if(vi == Vertex_handle())  //no conflict
        continue;
      else if(vi->info() < this_id) // already in star set (should not happen)
      {
        std::cout << "Warning! Insertion of p"<< this_id
          << " (" << p << ") in S" << si->index_in_star_set() << " failed. vi->info() :"<< vi->info() << std::endl;
        if(v_in_ch != Vertex_handle())
          m_ch_triangulation.remove(v_in_ch);
        remove_from_stars(this_id, target_stars.begin(), ++it);

        si->print_vertices(true);
        si->print_faces();
        std::cout << "Metric : \n" << si->metric().get_transformation() << std::endl;
        return vi->info();
      }
      else // inserted, standard configuration
      {
        modified_stars.insert(si->index_in_star_set());
        // update triangulation of convex hull
        if(infinite_stars)//(si->is_infinite() && v_in_ch != Vertex_handle())
        {
          v_in_ch = m_ch_triangulation.insert(p);
          v_in_ch->info() = this_id;
        }
      }
    }
    return this_id;
  }

  Index insert_to_stars(const Point_3& p,
                       Index_set& modified_stars,
                       const bool conditional)
  {
    Index this_id = static_cast<Index>(m_stars.size());

    // find conflicted stars & insert p to these stars
    Index id;
    Star_set target_stars;
    if(conditional)
      finite_stars_in_conflict(p, std::inserter(target_stars, target_stars.end())); // aabb tree/exhaustive
    else
      all_stars(std::inserter(target_stars, target_stars.end()));

    id = perform_insertions(p, this_id, target_stars, modified_stars, conditional, false);

    target_stars.clear();
    infinite_stars_in_conflict(p, std::inserter(target_stars, target_stars.end()));//convex hull
    id = perform_insertions(p, this_id, target_stars, modified_stars, conditional, true);

    return this_id;
  }

  Index simulate_insert_to_stars(const Point_3& p,
                                 Index_set& modified_stars) const
  {
    Index this_id = static_cast<Index>(m_stars.size());

    // find conflicted stars
    typename std::vector<Star_handle>::const_iterator it = m_stars.begin();
    typename std::vector<Star_handle>::const_iterator itend = m_stars.end();
    for(; it != itend; it++)
    {
      Star_handle si = get_star(it);
      int id = si->simulate_insert_to_star(p, this_id);

      if(id == -1)
        continue;
      else if(id < (int)this_id)
        return id;
      else
        modified_stars.insert(si->index_in_star_set());
    }
    return this_id;
  }

  Index insert(const Point_3& p, const bool conditional)
  {
    Index_set modified_stars;
    return insert(p, modified_stars, conditional);
  }

  Index insert(const Point_3 &p,
               Index_set& modified_stars,
               const bool conditional)
  {
    Index id = insert_to_stars(p, modified_stars, conditional);
    if(id < 0 || id < (int)m_stars.size())
      return id;

    Star_handle star = create_star(p, id, modified_stars);
    m_stars.push_back(star);
    modified_stars.insert(star->index_in_star_set());
    m_kd_tree.insert(star->index_in_star_set());
#ifndef NO_USE_AABB_TREE_OF_BBOXES
    m_aabb_tree.insert(AABB_primitive(star));
#endif

    return id;
  }

  Index insert_in_domain(const Point_3& p, const bool conditional)
  {
    Index_set modified_stars;
    return insert_in_domain(p, modified_stars, conditional);
  }

  Index insert_in_domain(const Point_3& p,
                         Index_set& modified_stars,
                         const bool conditional)
  {
    Index id = insert_to_stars(p, modified_stars, conditional);

    if(id < 0 || id < (int)m_stars.size())
      return id;

    Star_handle star = create_inside_star(p, id, modified_stars);
    m_stars.push_back(star);
    modified_stars.insert(star->index_in_star_set());
    m_kd_tree.insert(star->index_in_star_set());
#ifndef NO_USE_AABB_TREE_OF_BBOXES
    m_aabb_tree.insert(AABB_primitive(star));
#endif

    return id;
  }

  template<typename Element, typename OneMap>
  void add_to_map(const Element& e, OneMap& m) const
  {
    std::pair<typename OneMap::iterator, bool> is_insert_successful;
    is_insert_successful = m.insert(std::make_pair(e,1));
    if(!is_insert_successful.second)
      (is_insert_successful.first)->second += 1; // m[e] += 1
  }

  void cells_created(Star_handle star,
                      std::map<Cell_ijkl, int>& cells) const
  {
    Cell_handle_handle ci = star->begin_finite_star_cells();
    Cell_handle_handle ciend = star->end_finite_star_cells();
    for(; ci != ciend; ci++)
      add_to_map(Cell_ijkl(*ci), cells);
  }

  // finite cells(i,j,k,l) which would be created by the insertion
  // of 'p' in 'star' are added to 'cells'
  void cells_created(const Point_3& p, // new point, not yet in the star set
                      const int p_id,   // index of p in star set
                      Star_handle star,
                      std::map<Cell_ijkl, int>& cells) const
  {
    std::vector<Facet> bfacets;
    star->find_conflicts(p, std::back_inserter(bfacets));

    typename std::vector<Facet>::const_iterator fit = bfacets.begin();
    for( ; fit != bfacets.end(); fit++)
    {
      int id1 = fit->first->vertex((fit->second+1)%4)->info();
      int id2 = fit->first->vertex((fit->second+2)%4)->info();
      int id3 = fit->first->vertex((fit->second+3)%4)->info();
      int sid = star->index_in_star_set();

      if(id1 == star->infinite_vertex_index() ||
         id2 == star->infinite_vertex_index() ||
         id3 == star->infinite_vertex_index())
        continue;

      if(id1 != sid && id2 != sid && id3 != sid)
        continue;

      Point_3 c = compute_circumcenter(m_stars[id1]->center_point(),
                                       m_stars[id2]->center_point(),
                                       m_stars[id3]->center_point(),
                                       p, star);
      if(is_inside_domain(c))
        add_to_map(Cell_ijkl(id1, id2, id3, p_id), cells);
    }
  }

  bool check_consistency_and_sliverity(Star_handle to_be_refined,
                                       const Star_handle& new_star,
                                       const Index_set& modified_stars,
                                       FT sq_radius_bound) const
  {
    Point_3 p = new_star->center_point();
    int p_index = new_star->index_in_star_set();

    std::map<Cell_ijkl, int> cells;
    cells_created(new_star, cells);

    typename Index_set::iterator it;
    for(it = modified_stars.begin(); it != modified_stars.end(); ++it)
      cells_created(p, p_index, m_stars[(*it)], cells);

    typename std::map<Cell_ijkl, int>::iterator itc;
    for(itc = cells.begin(); itc != cells.end(); ++itc)
    {
      std::size_t nmax = number_of_stars();
      int c0 = (*itc).first.vertex(0);
      int c1 = (*itc).first.vertex(1);
      int c2 = (*itc).first.vertex(2);
      int c3 = (*itc).first.vertex(3);
      if((*itc).second != 4 )
      {
        if((*itc).first.is_infinite())
          return false;

        TPoint_3 tp0, tp1, tp2, tp3;
        TPoint_3 tp = transform_to_star_point(p, to_be_refined);
        tp0 = (c0 == nmax) ? tp : transform_to_star_point(m_stars[c0]->center_point(),
                                                          to_be_refined);
        tp1 = (c1 == nmax) ? tp : transform_to_star_point(m_stars[c1]->center_point(),
                                                          to_be_refined);
        tp2 = (c2 == nmax) ? tp : transform_to_star_point(m_stars[c2]->center_point(),
                                                          to_be_refined);
        tp3 = (c3 == nmax) ? tp : transform_to_star_point(m_stars[c3]->center_point(),
                                                          to_be_refined);

        double sqr = to_be_refined->compute_squared_circumradius(tp0, tp1, tp2, tp3);
        if(sqr < sq_radius_bound)
          return false; // small inconsistency (radius) is forbidden.
      }

      //check sliverity
      for(int i=0; i<4; ++i)
      {
        Star_handle star;
        if((*itc).first.vertex(i) == nmax)
          star = new_star;
        else
          star = m_stars[(*itc).first.vertex(i)];

        Point_3 p0, p1, p2, p3;
        p0 = (c0 == nmax) ? p : m_stars[c0]->center_point();
        p1 = (c1 == nmax) ? p : m_stars[c1]->center_point();
        p2 = (c2 == nmax) ? p : m_stars[c2]->center_point();
        p3 = (c3 == nmax) ? p : m_stars[c3]->center_point();

        FT sliver_overflow = star->compute_sliverity_overflow(p0, p1, p2, p3);
        if (sliver_overflow > 0.0)
          return false;

        //add squared_radius_bound & volume? todo
      }
    }
    return true;
  }

  bool is_valid_point(const Point_3 &p,
                      const FT sq_radiusbound,
                      const Star_handle& star,
                      Star_handle& new_star)
  {
    //simulate insert to stars
    Index_set modified_stars;
    int id = simulate_insert_to_stars(p, modified_stars);
    if(id < 0)
      return false;
    else if(id < (int)m_stars.size())
    {
      std::cout << "already in star set is_valid_point" << std::endl;
      new_star = m_stars[id];
      return false;
    }
    else
      create_star(p, id, modified_stars, new_star, new_star->is_surface_star());

    bool is = check_consistency_and_sliverity(star, new_star,
                                              modified_stars, sq_radiusbound);
    if(!is)
    {
      new_star->invalidate();
      std::cout << ("C");
    }

    return is;
  }

  Point_3 pick_valid(const Star_handle star, const Cell_handle &cell)
  {
    //std::cout << "pick valid : " << star->index_in_star_set() << std::endl;

    TPoint_3 circumcenter = cell->circumcenter(*(star->traits()));
    TPoint_3 tp2 = cell->vertex(0)->point();

    FT sq_circumradius = star->traits()->compute_squared_distance_3_object() (circumcenter, tp2);
    FT sq_radiusbound = sq_circumradius *  m_criteria->beta * m_criteria->beta;
    FT circumradius = sqrt(sq_circumradius) * m_criteria->delta;
    typename Star::Traits::Compute_random_point_3 random =
             star->traits()->compute_random_point_3_object();

    Star_handle new_star = new Star(m_criteria, m_pConstrain, false/*surface star*/);
    std::size_t tried_times = 0;

    while(true)
    {
      pick_count++;

      TPoint_3 tp = random(circumcenter,
        ((FT)tried_times / m_criteria->max_times_to_try_in_picking_region) * circumradius);
      Point_3 p = transform_from_star_point(tp, star);

      if(is_valid_point(p, sq_radiusbound, star, new_star))
        return p;
      if((tried_times++) > m_criteria->max_times_to_try_in_picking_region)
      {
        std::cout << ("X");
        return transform_from_star_point(circumcenter, star);
      }
    }
    delete new_star;
  }

  Point_3 compute_random_steiner_point(const Star_handle star,
                                       const Facet& facet,
                                       const TPoint_3& tccf,
                                       const FT& circumradius) const
  {
    typename Star::Traits::Compute_random_point_3 random =
             star->traits()->compute_random_point_3_object();

    bool found_acceptable_steiner_point = false;
    int failures_count = 0;
    Point_3 steiner_point;
    TPoint_3 random_point_within_sphere_1;
    TPoint_3 random_point_within_sphere_2;

    while(!found_acceptable_steiner_point)
    {
      if(++failures_count > 100)
        std::cout << "failures_count is getting big : " << failures_count << std::endl;

      random_point_within_sphere_1 = random(tccf, circumradius);
      random_point_within_sphere_2 = random(tccf, circumradius);

      found_acceptable_steiner_point =
            star->compute_steiner_dual_intersection(steiner_point,
                                                    random_point_within_sphere_1,
                                                    random_point_within_sphere_2,
                                                    tccf,
                                                    facet,
                                                    circumradius,
                                                    failures_count);
    }
    return steiner_point;
  }

  Point_3 pick_valid(const Star_handle star, const Facet &facet)
  {
    Point_3 center;
    star->compute_exact_dual_intersection(facet, center);
    TPoint_3 tccf = transform_to_star_point(center, star);
    TPoint_3 tp2 = facet.first->vertex((facet.second + 1) % 4)->point();

    FT sq_circumradius = star->traits()->compute_squared_distance_3_object()(tccf, tp2);
    FT sq_radiusbound = sq_circumradius * m_criteria->beta * m_criteria->beta;
    FT circumradius = sqrt(sq_circumradius) * m_criteria->delta;

    Star_handle new_star = new Star(m_criteria, m_pConstrain, true/*surface star*/);
    std::size_t tried_times = 0;

    while(true)
    {
      pick_count++;
      intersection_count++;
      Point_3 p = compute_random_steiner_point(star, facet, tccf, circumradius);

      if(is_valid_point(p, sq_radiusbound, star, new_star))
        return p;
      if((tried_times++) > m_criteria->max_times_to_try_in_picking_region)
      {
        std::cout << ("X\n");
        return compute_insert_or_snap_point(star, facet);
      }
    }
    delete new_star;
  }

  Point_3 compute_insert_or_snap_point(const typename Constrain_surface::EdgeIterator &edge)
  {
    //todo is this correct?
    Point_3 p = edge->first;
    Point_3 q = edge->second;
    Point_3 c((p.x() + q.x()) / 2.0, (p.y() + q.y()) / 2.0, (p.z() + q.z()) / 2.0);
    m_pConstrain->edge_split(edge, c);
    return c;
  }

  Point_3 compute_insert_or_snap_point(const Star_handle star, const Facet &facet)
  {
    intersection_count++;
    Point_3 p;
    star->compute_exact_dual_intersection(facet, p);
    //CHECK_EDGE_ENCROACHMENT
    return p;
  }

  Point_3 compute_insert_or_snap_point(const Star_handle star, const Cell_handle& cell)
  {
    Point_3 p = transform_from_star_point(star->compute_circumcenter(cell), star);
    /*
    std::cout << "in star " << star->index_in_star_set() << " cell : " << std::endl;
    std::cout << cell->vertex(0)->info() << " " << cell->vertex(0)->point() << std::endl;
    std::cout << cell->vertex(1)->info() << " " << cell->vertex(1)->point() << std::endl;
    std::cout << cell->vertex(2)->info() << " " << cell->vertex(2)->point() << std::endl;
    std::cout << cell->vertex(3)->info() << " " << cell->vertex(3)->point() << std::endl;
    std::cout << "insert or snap : " << p << std::endl;

    std::cout << m_stars[0]->bbox() << std::endl;
    Cell_handle useless;
    std::cout << "check conflicts : " << std::endl;
    std::cout << m_stars[cell->vertex(0)->info()]->is_conflicted(p, useless) << " "; //no transf coz iso met
    std::cout << m_stars[cell->vertex(1)->info()]->is_conflicted(p, useless) << " "; //no transf coz iso met
    std::cout << m_stars[cell->vertex(2)->info()]->is_conflicted(p, useless) << " "; //no transf coz iso met
    std::cout << m_stars[cell->vertex(3)->info()]->is_conflicted(p, useless) << std::endl; //no transf coz iso met
    */

    //CHECK_EDGE_ENCROACHMENT
    //CHECK_FACE_ENCROACHMENT
    //ONLY_CONSIDER_NEIGHBORING_ENCROACHMENT
    return p;
  }

  Point_3 compute_insert_or_snap_valid_point(const Star_handle star,
                                             const Facet &facet)
  {
    Point_3 p = pick_valid(star, facet);
    //CHECK_EDGE_ENCROACHMENT
    return p;
  }

  Point_3 compute_insert_or_snap_valid_point(const Star_handle& star,
                                             const Cell_handle& cell)
  {
    Point_3 p = pick_valid(star, cell);
    // test encroachment
    //CHECK_EDGE_ENCROACHMENT
    //CHECK_FACE_ENCROACHMENT
    //ONLY_CONSIDER_NEIGHBORING_ENCROACHMENT
    return p;
  }


  FT compute_distortion(const Cell_handle& c) const
  {
    FT distortion = 1.0;
    for(int i = 0; i < 4; i++)
    {
      for(int j = i + 1; j < 4; j++)
      {
        distortion = (std::max)(distortion,
          m_stars[c->vertex(i)->info()]->metric().compute_distortion(
          m_stars[c->vertex(j)->info()]->metric()));
      }
    }
    return distortion;
  }

  bool is_consistent(const Cell_handle& c,
                     std::vector<bool>& inconsistent_points,
                     const bool verbose = false) const
  {
    bool retval = true;
    bool cell_told = false;
    for(int i = 0; i < 4; i++)
    {
      int index = c->vertex(i)->info();
      if(is_infinite_vertex(index))
        continue;
      if(!m_stars[index]->has_cell(c))
      {
        if(verbose)
        {
          if(!cell_told)
          {
            m_stars[index]->cell_indices(c);
            std::cout << " inconsistent : ";
            cell_told = true;
          }
          std::cout << c->vertex(0)->info() << " " << c->vertex(1)->info() << " ";
          std::cout << c->vertex(2)->info() << " " << c->vertex(3)->info();
          std::cout << " not in S_" << index << ", ";
        }
        retval = false;
        inconsistent_points[i] = true;
      }
    }
    if(verbose && !retval)
      std::cout << "." << std::endl;
    return retval;
  }

  bool is_consistent(const Cell_handle& c,
                     const bool verbose = false) const
  {
    std::vector<bool> not_used(4);
    return is_consistent(c, not_used, verbose);
  }

  bool is_consistent(const bool verbose = true) const
  {
    std::size_t N = m_stars.size();
    for(std::size_t i = 0; i < N; i++)
    {
      Star_handle star = get_star(i);
      Cell_handle_handle cit = star->begin_finite_star_cells();
      Cell_handle_handle citend = star->end_finite_star_cells();
      for(; cit != citend; cit++)
      {
        if(!star->is_inside(*cit))
          continue;
        if(!is_consistent(*cit, verbose))
          return false;
      }
    }
    return true;
  }

  template<typename Stars>
  void fill_refinement_queue(const Stars &modified_stars, int relative_point = -1)
  {
    typename Stars::const_iterator si = modified_stars.begin();
    typename Stars::const_iterator siend = modified_stars.end();
    for(; si != siend; si++)
    {
      Star_handle star = get_star(si);

      /*
      typename Star::Facet_set_iterator fit = star->begin_restricted_facets();
      typename Star::Facet_set_iterator fend = star->end_restricted_facets();
      for(; fit != fend; fit++)
      {
        if((relative_point >= 0) && (fit->first->vertex(fit->second)->info() != relative_point))
          continue;
       //CHECK_FACE_ENCROACHMENT
      }
      */

      Cell_handle_handle ci = star->begin_finite_star_cells();
      Cell_handle_handle ciend = star->end_finite_star_cells();
      for(; ci != ciend; ci++)
      {
        Cell_handle c = *ci;
        if(!star->is_inside(c))
          continue;
        assert(!star->is_infinite(c));

        if(relative_point >= 0)
        {
          bool relative = false;
          for(int i = 0; i < 4; i++)
            if(relative_point == c->vertex(i)->info())
            {
              relative = true;
              break;
            }
            if(!relative)
              continue;
        }

        /*
        std::cout << "fill ";
        std::cout << c->vertex(0)->info() << " ";
        std::cout << c->vertex(1)->info() << " ";
        std::cout << c->vertex(2)->info() << " ";
        std::cout << c->vertex(3)->info() << std::endl;
        */

        // over distortion 1
        if(0 && m_criteria->distortion > 0.)
        {
          FT over_distortion = compute_distortion(c) - m_criteria->distortion;
          if(over_distortion > 0.)
          {
            m_refine_queue.push_over_distortion(star, c, over_distortion, 0);
            continue;
          }
        }

        // too big 2
        FT over_circumradius = star->compute_circumradius_overflow(c);
        if(over_circumradius > 0)
        {
          //std::cout << "over circum : " <<  over_circumradius << std::endl;
          m_refine_queue.push_over_circumradius(star, c, over_circumradius, 0);
          continue;
        }

        // bad shape 3
        FT over_radius_edge_ratio = star->compute_radius_edge_ratio_overflow(c);
        /*
        std::cout << "bad shape : " << star->criteria()->compute_squared_circumradius(c->vertex(0)->point(),
                                                                          c->vertex(1)->point(),
                                                                          c->vertex(2)->point(),
                                                                          c->vertex(3)->point());
        std::cout << " " << star->criteria()->compute_squared_shortest_edge(c->vertex(0)->point(),
                                                                            c->vertex(1)->point(),
                                                                            c->vertex(2)->point(),
                                                                            c->vertex(3)->point())
                 << std::endl;
        */
        if(over_radius_edge_ratio > 0)
        {
          m_refine_queue.push_bad_shape(star, c, over_radius_edge_ratio, 0);
          continue;
        }

        // sliverity 4
        FT over_sliverity = star->compute_sliverity_overflow(c);
        if(over_sliverity > 0)
        {
          m_refine_queue.push_sliver(star, c, over_sliverity, 0);
          continue;
        }

        // consistency 5
        if(!is_consistent(c))
        {
          m_refine_queue.push_inconsistent(star, c, star->compute_volume(c), 0);
          continue;
        }
      } // for cells
    } // for stars
  }

  bool next_refine_cell(Refine_cell &refine_cell,
                        Cell_handle &cell,
                        bool &need_picking_valid,
                        int &queue_type)
  {
    while(true)
    {
      //m_refine_queue.print();
      if(!m_refine_queue.top(refine_cell, queue_type))
        return false;
      m_refine_queue.pop();
      if(refine_cell.star->has_cell(refine_cell.vertices, cell))
      {
        if(!refine_cell.star->is_inside(cell))
        {
          std::cout << "out" << std::endl;
          continue;
        }
        need_picking_valid = m_refine_queue.need_picking_valid(queue_type);
        return true;
      }
      //else
      //  std::cout << "doesn't exist anymore" << std::endl;
    }
  }

  void clean_stars()
  {
    for(std::size_t i = 0; i < m_stars.size(); i++)
      m_stars[i]->clean();
  }

  bool refine()
  {
    Point_3 p;
    int queue_type = 0;
    Refine_cell bad_cell;
    bool need_picking_valid;
    Cell_handle c;

    if(!next_refine_cell(bad_cell, c, need_picking_valid, queue_type))
    {
      std::cout << "get next fail" << std::endl;
      return false;
    }
    if(need_picking_valid && compute_distortion(c) > m_criteria->distortion)
      need_picking_valid = false;

    Vertex_handle v1 = c->vertex(0);
    Vertex_handle v2 = c->vertex(1);
    Vertex_handle v3 = c->vertex(2);
    Vertex_handle v4 = c->vertex(3);

/*
#ifdef ONLY_CONSIDER_SURFACE
    if(!bad_cell.star->is_boundary_star())
      return true;
#endif
*/

    if(queue_type == 0)
    { // encroachment
      std::cerr << "Error : encroachment is not implemented.\n";
      return true;
      /*
      int tag = 0;
      for(int i = 0; i < 4; i++)
        if(c->vertex(i)->info() == bad_cell.vertices[0])
        {
          tag = i;
          break;
        }
      if(!bad_cell.star->has_facet(Facet(c, tag)))
      {
        std::cout << ("[repick]");
        return true; //todo check that this is correct
      }
      p = compute_insert_or_snap_point(bad_cell.star, Facet(c, tag));
      vertex_without_picking_count++;
      */
    }
    else
    {
      if(need_picking_valid)
      {
        p = compute_insert_or_snap_valid_point(bad_cell.star, c);
        vertex_with_picking_count++;
      }
      else
      {
        p = compute_insert_or_snap_point(bad_cell.star, c);
        vertex_without_picking_count++;
      }
    }

    std::cout << std::endl << number_of_stars() << " refine: " << p << " [" << c->vertex(0)->info() << "-";
    std::cout << c->vertex(1)->info() << "-" << c->vertex(2)->info() << "-";
    std::cout << c->vertex(3)->info() << "] " << "t:" << queue_type << " i:";
    std::cout << bad_cell.vertices[0] << " s:" << bad_cell.star->center()->info() << std::endl;

    Index_set modified_stars;
    Index pid = insert_in_domain(p, modified_stars, true);

    Cell_handle ctest;
    int i,j,k,l;
    if(bad_cell.star->is_cell(v1,v2,v3,v4,ctest,i,j,k,l))
    {
      std::cout << "bad cell still in..." << std::endl;
      return false;
    }

    if(!modified_stars.empty())
      fill_refinement_queue(modified_stars, pid);
    return true;
  }

  void print_stars()
  {
    Star_iterator it = m_stars.begin();
    Star_iterator itend = m_stars.end();

    for(; it!=itend; ++it)
    {
      Star_handle si = get_star(it);
      std::cout << "star " << si->index_in_star_set() << " || " << si->center_point() << " || ";
      std::cout << si->center()->info() << " || " << si->center()->point() << std::endl;
      std::cout << "neighbours : ";
      typename Star::Vertex_handle_handle nsi = si->begin_neighboring_vertices();
      typename Star::Vertex_handle_handle nsiend = si->end_neighboring_vertices();
      for(; nsi != nsiend; nsi++)
        std::cout << (*nsi)->info() << " ";
      std::cout << std::endl;

      /*
      std::cout << "double tap " << std::endl;
      Star_iterator it2 = m_stars.begin();
      Star_iterator it2end = m_stars.end();
      for(; it2!=it2end; ++it2)
      {
        Star_handle si2 = get_star(it2);
        if(si2->index_in_star_set() != si->index_in_star_set())
          si->insert_to_star(si2->center_point(), si2->index_in_star_set(), false);
      }

      nsi = si->begin_neighboring_vertices();
      for(; nsi != nsiend; nsi++)
        std::cout << (*nsi)->info() << " ";
      std::cout << std::endl;
      */

      Cell_handle_handle ci = si->begin_star_cells();
      Cell_handle_handle ciend = si->end_star_cells();
      for (; ci != ciend; ++ci)
      {
        Cell_handle c = *ci;
        if(!si->is_inside(c))
          continue;
        if(!is_consistent(c))
        {
          std::cout << c->vertex(0)->info() << " ";
          std::cout << c->vertex(1)->info() << " ";
          std::cout << c->vertex(2)->info() << " ";
          std::cout << c->vertex(3)->info() << " consist : " << is_consistent(c) << std::endl;
        }
      }
    }
  }

  void refine_all(std::size_t max_count = (std::size_t) -1)
  {
    pick_count = 0;
    intersection_count = 0;
    start_time = clock();

    std::cout << "\nInitializing conflicts...\n";
    fill_refinement_queue(m_stars, -1);
    m_refine_queue.print();
    std::cout << " [ok]" << std::endl;

    vertex_with_picking_count = 0;
    vertex_without_picking_count = (int)m_stars.size();

    std::size_t nbv = m_stars.size();
    while(nbv < max_count)
    {
      if(nbv % 100 == 0)
      {
        clean_stars();
        std::cout << " " << nbv << " vertices, ";
        std::cout << duration(start_time) << " sec.,\t";
        m_refine_queue.print();
        dump();
      }

      if(nbv % 1000 == 0)
      {
        std::ofstream fx("cell_mesher_partial.mesh");
        output_medit(fx);
      }

      if(!refine())
      {
        std::cout << "refine fail" << std::endl;
        clean_stars();
        break;
      }

      if(0 && number_of_stars() > 10000)
        break;

      nbv = m_stars.size();
    }

    //checking if really done!
    fill_refinement_queue(m_stars, -1);
    if(!m_refine_queue.empty())
      std::cout << "Stopped too early! Not done!" << std::endl;
    else
      std::cout << "Finished properly" << std::endl;
    m_refine_queue.print();

    //print_stars();
    std::cout << "triangulation is consistent : " << is_consistent() << std::endl;
    std::cout << duration(start_time) << " sec." << std::endl;
    std::cout << number_of_stars() << " stars" << std::endl;
    report();
    histogram_vertices_per_star<Self>(*this);
  }

  void initialize_stars(const int nb = 10)
  {
    std::cout << "Initializing stars...";
    int nbdone = 0;

    typename Constrain_surface::Pointset initial_points = m_pConstrain->get_surface_points(nb);
    typename Constrain_surface::Pointset::iterator pi = initial_points.begin();
    typename Constrain_surface::Pointset::iterator pend = initial_points.end();
    for(int i = 0; pi != pend; pi++, i++)
    {
      std::size_t this_id = m_stars.size();
      int id = -1;
      id = insert_in_domain(*pi, false/*under no condition*/); //points on the surface are marked as inside stars here...
      if(this_id == id)                                        //The problem is that if a point is a surface star, the eval
        nbdone++;                                              //of the the conflict zone is different in stretched delaunay.
    }                                                          //Need something better todo.
    std::cout << "done : " << nbdone << std::endl;
    clean_stars();
  }

  void dump()
  {
    std::ofstream fx("dump.txt");

    std::size_t ns = number_of_stars();
    fx << number_of_stars() << std::endl;
    for(std::size_t i=0; i<ns; ++i)
      fx << get_star(i)->center_point() << std::endl;

    for(std::size_t i=0; i<ns; ++i)
    {
      Star_handle star_i = get_star(i);
      typename Star::Vertex_handle_handle nsi = star_i->begin_neighboring_vertices();
      typename Star::Vertex_handle_handle nsiend = star_i->end_neighboring_vertices();
      fx << nsiend - nsi;
      for(; nsi != nsiend; nsi++)
        fx << " " << (*nsi)->info();
      fx << std::endl;
    }
  }

  void save_off(const bool nb = false)
  {
    std::ostringstream nbs;
    nbs << "aniso_3D_";
    if(nb)
      nbs << number_of_stars();
    nbs << ".off" << std::ends;
    std::ofstream fx(nbs.str().c_str());

    std::cout << "Saving as .off...";
    std::map<Index, int> match_indices;//because we won't use all of them
    int off_index = 0; // the corresponding index in .off
    std::vector<Point_3> points;
    Output_facets output_facets;
    Output_cells output_cells;

    Star_iterator it = m_stars.begin();
    Star_iterator itend = m_stars.end();
    for (; it != itend; ++it)
    {
      Star_handle star = *it;

      points.push_back(star->center_point());
      match_indices[star->index_in_star_set()] = off_index++;

      Cell_handle_handle ci = star->begin_star_cells();
      Cell_handle_handle ciend = star->end_star_cells();
      for (; ci != ciend; ++ci)
      {
        Cell_handle c = *ci;
        if(!star->is_infinite(c) && is_consistent(c))
        {
          output_cells.insert(c->vertex(0)->info(), c->vertex(1)->info(),
                              c->vertex(2)->info(), c->vertex(3)->info());
          output_facets.insert(c->vertex(0)->info(), c->vertex(1)->info(), c->vertex(2)->info());
          output_facets.insert(c->vertex(0)->info(), c->vertex(2)->info(), c->vertex(3)->info());
          output_facets.insert(c->vertex(0)->info(), c->vertex(1)->info(), c->vertex(3)->info());
          output_facets.insert(c->vertex(3)->info(), c->vertex(1)->info(), c->vertex(2)->info());
        }
      }
    }

    fx << "OFF" << std::endl;
    fx << points.size() << " " << output_facets.size() << " " << output_cells.size() << std::endl;
    for(unsigned int i = 0; i < points.size(); i++)
      fx << points[i] << std::endl;

    Output_facets::Facet_handle ofi = output_facets.begin();
    Output_facets::Facet_handle ofiend = output_facets.end();
    for (; ofi != ofiend; ofi++)
    {
      fx << "3  " << match_indices[ofi->vertices[0] ]
         << " "   << match_indices[ofi->vertices[1] ]
         << " "   << match_indices[ofi->vertices[2] ] << std::endl;
    }

    Output_cells::Cell_handle oci = output_cells.begin();
    Output_cells::Cell_handle ociend = output_cells.end();
    for (; oci != ociend; oci++)
    {
      fx << "4  " << match_indices[oci->vertices[0] ]
         << " "   << match_indices[oci->vertices[1] ]
         << " "   << match_indices[oci->vertices[2] ]
         << " "   << match_indices[oci->vertices[3] ] << std::endl;
    }

    std::cout << "done.\n";
  }

  void report()
  {
    typename std::ofstream fx("report.txt");

    fx << "[Parameters]" << std::endl << std::endl;
    m_criteria->report(fx);

    fx << std::endl << "[Metric field]" << std::endl << std::endl;
    m_metric_field->report(fx);

    fx << std::endl << "[Statistics]" << std::endl << std::endl;
    fx << "elapsed time:       " << duration(start_time) << " sec." << std::endl;
    fx << "number of vertices: " << m_stars.size() << std::endl;
    fx << "picking times:      " << pick_count << std::endl;
    fx << "vertex via picking: " << vertex_with_picking_count << std::endl;
    fx << "vertex non-picking: " << vertex_without_picking_count << std::endl;
    fx << "picking rate:       " << (double)vertex_with_picking_count / (double)(vertex_without_picking_count + vertex_with_picking_count) << std::endl;
    fx << "picking time rate:  " << (double)pick_count / (double)(pick_count + vertex_without_picking_count) << std::endl;
    fx << "mean picking times: " << (double)pick_count / (double)vertex_with_picking_count << std::endl;
    fx << "intersection count: " << intersection_count << std::endl;
  }

  void load_dump()
  {
    m_stars.clear();
    std::ifstream fx("dump.txt");
    assert(fx.is_open());
    std::cout << ("Loading dump file...");
    int ns, nn, id;
    fx >> ns;

    for(int i = 0; i < ns; i++)
    {
      Point_3 p;
      fx >> p;
      Star_handle star = new Star(m_criteria, m_pConstrain, false);
      star->reset(p, i, m_metric_field->compute_metric(p), false);
      m_stars.push_back(star);
    }

    for(int i = 0; i < ns; i++)
    {
      fx >> nn;
      Star_handle star_i = get_star(i);
      for(int i = 0; i < nn; i++)
      {
        fx >> id;
        if(id >= 0)
        {
          Star_handle star_id = get_star(id);
          star_i->insert_to_star(star_id->center_point(), star_id->index_in_star_set(), false);
        }
      }
    }

    std::cout << number_of_stars() << " stars" << std::endl;
  }

/*
  void output()
  {
    typename std::ofstream fx("mesh.volume.cgal");
    fx << 3 << std::endl;
    fx << m_points.size() << std::endl;
    for(int i = 0; i < (int)m_points.size(); i++)
      fx << m_points[i] << std::endl;

    Output_cells output_cells;
    Star_iterator it = m_stars.begin();
    Star_iterator itend = m_stars.end();
    for(; it != itend; it++)
    {
      Star_handle star = *it;
      Cell_handle_handle ci = star->begin_finite_star_cells();
      Cell_handle_handle ciend = star->end_finite_star_cells();
      for(; ci != ciend; ci++)
      {
        if(!star->is_inside(*ci))
          continue;
        bool consistent = true;
        for(int i = 0; i < 4; i++)
          if(!m_stars[(*ci)->vertex(i)->info()]->has_cell(*ci))
          {
            consistent = false;
            break;
          }
          if(consistent)
            output_cells.insert(
            (*ci)->vertex(0)->info(), (*ci)->vertex(1)->info(),
            (*ci)->vertex(2)->info(), (*ci)->vertex(3)->info());
      }
    }
    Output_cells::Cell_handle oci = output_cells.begin();
    Output_cells::Cell_handle ociend = output_cells.end();
    fx << output_cells.size() << std::endl;
    for(; oci != ociend; oci++)
    {
      fx << oci->vertices[0] + 1 << " " << oci->vertices[1] + 1 << " "
        << oci->vertices[2] + 1 << " " << oci->vertices[3] + 1 << std::endl;
    }
  }
*/

  void output_medit(std::ofstream& fx)
  {
    std::cout << "Saving medit" << std::endl;
    unsigned int nb_inconsistent_stars = 0;
    fx << "MeshVersionFormatted 1\n";
    fx << "Dimension 3\n";

    fx << "Vertices\n";
    fx << number_of_stars() << std::endl;
    for(std::size_t i = 0; i < number_of_stars(); i++)
      fx << get_star(i)->center_point() << " " << (i+1) << std::endl; // warning : indices start at 1 in Medit

    Output_cells output_cells;
    Star_iterator it = m_stars.begin();
    Star_iterator itend = m_stars.end();
    for(; it != itend; it++)
    {
      Star_handle star = *it;
      Cell_handle_handle ci = star->begin_finite_star_cells();
      Cell_handle_handle ciend = star->end_finite_star_cells();
      for(; ci != ciend; ci++)
      {
        if(!star->is_inside(*ci))
          continue;
        if(is_consistent(*ci))
          output_cells.insert((*ci)->vertex(0)->info(), (*ci)->vertex(1)->info(),
                              (*ci)->vertex(2)->info(), (*ci)->vertex(3)->info());
        else
          nb_inconsistent_stars++;
      }
    }
    fx << "Tetrahedra\n";
    fx << output_cells.size() << std::endl;
    Output_cells::Cell_handle oci = output_cells.begin();
    Output_cells::Cell_handle ociend = output_cells.end();
    for(; oci != ociend; oci++)
    {
      fx << (oci->vertices[0] + 1) << " " << (oci->vertices[1] + 1) << " "
         << (oci->vertices[2] + 1) << " " << (oci->vertices[3] + 1) << " "
         << "1" << std::endl;
    }
    if(nb_inconsistent_stars > 0)
      std::cout << "Warning : there are " << nb_inconsistent_stars << " inconsistent stars in the ouput mesh.\n";
    std::cout << "done" << std::endl;
  }

public:
  Cell_star_set_3(const Criteria* criteria_,
                  const Metric_field* metric_field_,
                  const Constrain_surface* const pconstrain_,
                  int nb = 10) :
    m_pConstrain(pconstrain_),
    m_metric_field(metric_field_),
    m_criteria(criteria_),
    m_stars(),
    m_refine_queue(),
    m_ch_triangulation(),
    m_aabb_tree(100/*insertion buffer size*/),
    m_kd_tree(m_stars)
  {
    if(nb)
      initialize_stars(nb);
    else
      load_dump();

    m_ch_triangulation.infinite_vertex()->info() = -10;

#ifndef NO_USE_AABB_TREE_OF_BBOXES
    build_aabb_tree();  // update bboxes and rebuild
#endif

  }

  ~Cell_star_set_3()
  {
    delete m_metric_field;
    delete m_criteria;
    Star_iterator it = m_stars.begin();
    Star_iterator itend = m_stars.end();
    for(; it != itend; it++)
      delete (*it);
  }
};

} //namespace Anisotropic_mesh_3
} //namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_STAR_SET_3_H

﻿// Copyright (c) 2011  INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
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

#ifndef CGAL_ANISOTROPIC_MESH_3_SURFACE_STAR_SET_3_H
#define CGAL_ANISOTROPIC_MESH_3_SURFACE_STAR_SET_3_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <stdio.h>
#include <map>
#include <algorithm>
#include <limits>
#include <set>
//#include <boost/tuple/tuple.hpp>

#include <CGAL/Timer.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/iterator.h>
#include <CGAL/Stretched_Delaunay_3.h>
#include <CGAL/Constrain_surface_3.h>
#include <CGAL/Metric_field.h>
#include <CGAL/Criteria.h>
#include <CGAL/Facet_refine_queue.h>
#include <CGAL/Output_facets.h>

#include <CGAL/helpers/statistics_helper.h>
#include <CGAL/helpers/timer_helper.h>
#include <CGAL/helpers/combinatorics_helper.h>

#include <CGAL/aabb_tree/aabb_tree_bbox.h>
#include <CGAL/aabb_tree/aabb_tree_bbox_primitive.h>

#include <CGAL/kd_tree/Kd_tree_for_star_set.h>

#include <CGAL/Profile_counter.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>

template<typename PointType>
struct No_condition
{
  bool operator()(const PointType& p) const
  {
    return true;
  }
};

namespace CGAL
{
  namespace Anisotropic_mesh_3
  {
    template<typename K, typename RefinementCondition = No_condition<typename K::Point_3> >
    class Surface_star_set_3 
    {
    public:
      typedef Surface_star_set_3<K, RefinementCondition> Self;

      //typedef typename CGAL::Exact_predicates_exact_constructions_kernel KExact;
      typedef  K KExact;

      typedef typename K::FT                      FT;
      typedef Stretched_Delaunay_3<K, KExact>     Star;
      typedef Star*                               Star_handle;
      typedef typename Star::Base                 DT; // DT_3 with vertex_base_with_info
      typedef std::vector<Star_handle>            Star_vector;
      typedef std::set<Star_handle>               Star_set;
      typedef typename Star_vector::iterator      Star_iterator;
      typedef typename Star::Index                Index;
      typedef std::set<Index>                     Index_set;
      
      typedef Constrain_surface_3<K>              Constrain_surface;
      typedef CGAL::Anisotropic_mesh_3::Metric_field<K> Metric_field;
      typedef typename Metric_field::Metric       Metric;
      typedef typename K::Point_3                 Point_3;
      typedef Criteria_base<K>                    Criteria;

      typedef typename Star::Vertex_handle        Vertex_handle;
      typedef typename Star::Cell_handle          Cell_handle;
      typedef typename Star::Facet_handle         Facet_handle;
      typedef typename Star::Facet_set_iterator   Facet_set_iterator;
      typedef typename Star::Cell_handle_handle   Cell_handle_handle;
      typedef typename Star::Facet                Facet;
      typedef typename Star::Edge                 Edge;
      typedef typename Star::TPoint_3             TPoint_3;
      typedef typename Star::Vector_3             Vector_3;

      typedef Facet_refine_queue<K, KExact>                     Refine_queue;
      typedef CGAL::Anisotropic_mesh_3::Refine_facet<K, KExact> Refine_facet;

      typedef CGAL::AABB_tree_bbox<K, Star>   AABB_tree;
      typedef CGAL::AABB_bbox_primitive<Star> AABB_primitive;

      typedef CGAL::Kd_tree_for_star_set<K, Star_handle> Kd_tree;
      typedef typename Kd_tree::Traits                            Kd_traits;         
      typedef typename Kd_tree::Box_query                         Kd_Box_query;
      typedef typename Kd_tree::key_type                          Kd_point_info;

    public :
      struct Facet_ijk;
      struct Edge_ij;

    public:
      std::set<Point_3> poles;
      typename Constrain_surface::Pointset initial_points;
      const Constrain_surface* const m_pConstrain;
      const Metric_field* m_metric_field;
      const Criteria* m_criteria;
      Star_vector m_stars;
      Refine_queue m_refine_queue;
      DT m_ch_triangulation;
      RefinementCondition m_refinement_condition;

      AABB_tree m_aabb_tree; //bboxes of stars
      Kd_tree m_kd_tree;     //stars* centers for box queries

      // the following global variables are only 
      // for debugging or benchmark
      mutable int vertex_with_picking_count;
      mutable int vertex_without_picking_count;
      mutable time_t start_time;
      mutable int vertex_with_smoothing_counter;
      mutable int vertex_without_smoothing_counter;
     
    private:
      typedef typename KExact::Point_3                      Exact_Point_3;
      typedef typename KExact::Point_3                      Exact_TPoint_3;
      typedef CGAL::Cartesian_converter<K, KExact> To_exact;
      typedef CGAL::Cartesian_converter<KExact, K> Back_from_exact;
      To_exact to_exact;
      Back_from_exact back_from_exact;

    private:
      mutable std::vector<Point_3> m_pick_valid_facet; // facet f trying to be refined (for visu)
      mutable std::vector<Point_3> m_pick_valid_cache; // the random points picked to refine f
      //the correspondance between a random point and the n problematic facets in pick_valid
      //the vector is made of n*2 points, since the third point is the key of the map
      mutable std::map<Point_3,std::vector<int> > pickvalid_problematic_facets;

      mutable std::vector<Point_3> red_points;
      mutable std::vector<Point_3> orange_points;
      mutable std::vector<Point_3> yellow_points;
      mutable std::vector<Point_3> green_points;

      mutable Metric metricA, metricB, metricP;

#ifdef USE_ANISO_TIMERS
private:
      mutable double m_pick_valid_timer;
      mutable double m_insert_timer;
      mutable double m_fill_queue_timer;
      mutable double m_compute_steiner_timer;
      mutable double m_create_star_timer;
#endif

    public:
      std::size_t pickvalid_problematic_facets_size() const {return pickvalid_problematic_facets.size();}
      const Constrain_surface* const constrain_surface() const { return m_pConstrain; }
      const Criteria* criteria() const { return m_criteria; }
      const Metric_field* metric_field() const { return m_metric_field; }
      void set_criteria(const Criteria* criteria_) { m_criteria = criteria_; }
      void set_metric_field(const Metric_field* metric_field_) { m_metric_field = metric_field_; }

    public:
      AABB_tree& aabb_tree() { return m_aabb_tree; }
      const AABB_tree& aabb_tree() const { return m_aabb_tree; }

    public:
      bool empty() const { return m_stars.empty(); }
      std::size_t size() const { return m_stars.size(); }

    public:
      TPoint_3 transform_to_star_point(const Point_3& p, Star_handle star) const {
        return star->metric().transform(p);
      }
      Point_3 transform_from_star_point(const TPoint_3& p, Star_handle star) const {
        return star->metric().inverse_transform(p);
      }
      //Exact_TPoint_3 transform_to_star_point(const Exact_Point_3& p, Star_handle star) const {
      //  return star->metric().transform(p);
      //}
      //Exact_Point_3 transform_from_star_point(const Exact_TPoint_3& p, Star_handle star) const {
      //  return star->metric().inverse_transform(p);
      //}
      
    private:
      Star_handle get_star(Star_handle s) const { return s; }
      Star_handle get_star(int i) const         { return m_stars[i]; }
      Star_handle get_star(typename Index_set::const_iterator it) const   { return m_stars[*it]; }
      Star_handle get_star(typename Star_set::const_iterator it) const    { return *it; }
      Star_handle get_star(typename Star_vector::const_iterator it) const { return *it; }

      template <typename StarIterator>
      Star_set get_stars(StarIterator begin, StarIterator end) const
      {
        Star_set stars;
        for(StarIterator it = begin; it != end; ++it)
          stars.insert(get_star(it));
        return stars;
      }

    public:
      FT compute_approximation_error(const bool verbose = false) const
      {
//        std::vector<int> histogram(100, 0);
        FT sq_approx = 0.;
        typename Star_vector::const_iterator it;
        for(it = m_stars.begin(); it != m_stars.end(); ++it)
        {
          Star_handle s = get_star(it);
          if(!s->is_surface_star())
            continue;                    
          typename Star::Facet_set_iterator fit = s->begin_restricted_facets();
          typename Star::Facet_set_iterator fend = s->end_restricted_facets();
          for(; fit != fend; ++fit)
          {
            Facet f = *fit;
            Point_3 cc;
            s->compute_dual_intersection(f, cc);
            if(!m_refinement_condition(cc))
              continue;

            FT sqd = sq_distance_to_surface(f, s);
            //FT d = std::sqrt(sqd);
            //histogram[std::floor(100*d)]++;
            sq_approx = (std::max)(sq_approx, sqd);   
          }
        }
        //std::cout << "Histogram : " << std::endl;
        //for(unsigned int i = 0; i < 100; i++)
        //  std::cout << "\t" << i << " : " << histogram[i] << std::endl;
        return std::sqrt(sq_approx);
      }

      bool is_surface_facet(const Facet& f) const
      {
        return m_stars[f.first->vertex((f.second + 1) % 4)->info()]->is_surface_star()
            && m_stars[f.first->vertex((f.second + 2) % 4)->info()]->is_surface_star()
            && m_stars[f.first->vertex((f.second + 3) % 4)->info()]->is_surface_star();
      }
      
      FT sq_distance_to_surface(const Facet& f, const Star_handle s) const
      {
        //Point_3 steiner;
        //s->compute_dual_intersection(f, steiner);
        //Point_3 cc = transform_from_star_point(s->compute_circumcenter(f), s);
        //return CGAL::squared_distance(steiner, cc);
        
        return m_pConstrain->compute_sq_approximation(barycenter(f));
      }

      Point_3 barycenter(const Facet& f) const
      {
        Point_3 p1 = m_stars[f.first->vertex((f.second + 1) % 4)->info()]->center_point();
        Point_3 p2 = m_stars[f.first->vertex((f.second + 2) % 4)->info()]->center_point();
        Point_3 p3 = m_stars[f.first->vertex((f.second + 3) % 4)->info()]->center_point();
        FT third = 1./3.;
        return Point_3(third * (p1.x() + p2.x() + p3.x()),
                       third * (p1.y() + p2.y() + p3.y()),
                       third * (p1.z() + p2.z() + p3.z()));
      }

    public:
      bool is_valid_point(const Point_3 &p,
                          const FT& sq_radius_bound, // in M_{star to_be_refined}
                          Star_handle to_be_refined,
                          Star_handle& new_star,
                          std::vector<int>& problematic_facets) const
      {
        Index_set modified_stars;
        int id = simulate_insert_to_stars(p, modified_stars);
        if(id < 0)
          return false;
        else if(id < (int)m_stars.size())
        {
          new_star = m_stars[id]; //already in star set
          return false;
        }
        else
          create_star(p, id, modified_stars, new_star);

        bool is = check_consistency(to_be_refined, new_star, modified_stars, sq_radius_bound, problematic_facets);

        if(!is)
          new_star->invalidate();

        return is;
      }

      static void sort(typename boost::tuple<int, int, int>& f, const bool verbose = false)
      {
        if(verbose)
          std::cout << "\n{ "<< f.get<0>() << " " << f.get<1>() << " " << f.get<2>() << "}"; 

        if(f.get<0>() < f.get<1>())    // a < b here
          if(f.get<0>() < f.get<2>())  // a < c     : a the smallest
            if(f.get<1>() < f.get<2>())// b < c     : a < b < c
              return; // they are in correct order
            else                   // c <= b        : a < c <= b
              f = boost::make_tuple(f.get<0>(), f.get<2>(), f.get<1>());
          else                     // a >= c        : c <= a < b
              f = boost::make_tuple(f.get<2>(), f.get<0>(), f.get<1>());
        else                       // b <= a here
          if(f.get<1>() < f.get<2>())  // b < c     : b the smallest
            if(f.get<0>() < f.get<2>())// a < c     : b <= a < c
              f = boost::make_tuple(f.get<1>(), f.get<0>(), f.get<2>());
            else                   // a >= c        : b < c <= a
              f = boost::make_tuple(f.get<1>(), f.get<2>(), f.get<0>());
        else                       // c <= b        : c <= b <= a
         f = boost::make_tuple(f.get<2>(), f.get<1>(), f.get<0>());
        
        if(verbose)
          std::cout << "\t --> { "<< f.get<0>() << " " << f.get<1>() << " " << f.get<2>() << "}";
      }

      struct Facet_ijk
      {
      private:
        typename boost::tuple<int,int,int> m_facet;
      public :
        Facet_ijk(const Facet& f)
        {
          m_facet = boost::make_tuple(
            f.first->vertex((f.second + 1) % 4)->info(),
            f.first->vertex((f.second + 2) % 4)->info(),
            f.first->vertex((f.second + 3) % 4)->info());
          sort(m_facet);
        }
        Facet_ijk(const int i = -1, const int j = -1, const int k = -1)
        {
          m_facet = boost::make_tuple(i,j,k);
          sort(m_facet);
        }
        int vertex(const int index) const
        {
          if(index == 0)      return m_facet.get<0>();
          else if(index == 1) return m_facet.get<1>();
          else if(index == 2) return m_facet.get<2>();
          else std::cerr << "Facet_ijk does not have vertex number " << index << "!\n";
          return -1;
        }
        bool has(const int index) const
        {
          return (vertex(0) == index) || (vertex(1) == index) || (vertex(2) == index);
        }
        bool has(const int index, int& other1, int& other2) const
        {
          for(int i = 0; i < 3; i++)
          {
            if(vertex(i) == index)
            {
              other1 = vertex((i+1)%3);
              other2 = vertex((i+2)%3);
              return true;
            }
          }
          return false;
        }
        bool operator==(const Facet_ijk& f) const
        {
          return vertex(0) == f.vertex(0)
            && vertex(1) == f.vertex(1)
            && vertex(2) == f.vertex(2);
        }
        bool operator<(const Facet_ijk& f) const
        {
          if (this->vertex(0) != f.vertex(0))   return (this->vertex(0) < f.vertex(0));
          if (this->vertex(1) != f.vertex(1))   return (this->vertex(1) < f.vertex(1));
          if (this->vertex(2) != f.vertex(2))   return (this->vertex(2) < f.vertex(2));
          return false;
        }
        Facet_ijk& operator=(const Facet_ijk& f)
        {
          m_facet = boost::make_tuple(f.vertex(0), f.vertex(1), f.vertex(2));
          return *this;
        }
        bool is_infinite() const
        {
          return has(-10);
        }
      }; // end struct Facet_ijk
      
      struct Edge_ij
      {
      private:
        std::pair<int, int> m_edge;
      public:
        Edge_ij(const Edge& e)
        {
          int i = e.first->vertex(e.second)->info();
          int j = e.first->vertex(e.third)->info();
          m_edge = std::pair<int, int>((std::min)(i,j), (std::max)(i,j));
        }
        Edge_ij(const int i = -1, const int j = -1)
          : m_edge(std::pair<int, int>((std::min)(i,j), (std::max)(i,j))) {}

        Edge_ij(const Edge_ij& e)
          : m_edge(std::pair<int, int>(e.vertex(0), e.vertex(1))) {}

        int vertex(const int index) const
        {
          if(index == 0)      return m_edge.first;
          else if(index == 1) return m_edge.second;
          std::cerr << "Edge_ij does not have vertex number " << index << "!\n";
          return -1;
        }
        bool has(const int i)
        {
          return (vertex(0) == i) || (vertex(1) == i);
        }
        bool has(const int& i, int& other)
        {
          if(vertex(0) == i)
          {
            other = vertex(1);
            return true;
          }
          if(vertex(1) == i)
          {
            other = vertex(0);
            return true;
          }
          return false;
        }
        bool operator==(const Edge_ij& e) const
        {
          return vertex(0) == e.vertex(0) && vertex(1) == e.vertex(1);
        }
        bool operator<(const Edge_ij& e) const
        {
          if (this->vertex(0) != e.vertex(0))   return (this->vertex(0) < e.vertex(0));
          if (this->vertex(1) != e.vertex(1))   return (this->vertex(1) < e.vertex(1));
          return false;
        }          
        Edge_ij& operator=(const Edge_ij& e)
        {
          m_edge = std::pair<int, int>(e.vertex(0), e.vertex(1));
          return *this;
        }
        bool is_infinite()
        {
          return has(-10);
        }
        bool is_valid()
        {
          return !has(-1);
        }
      };

      std::size_t count_restricted_facets() const
      {
        typename std::set<Facet_ijk> facets;
        for(unsigned int i = 0; i < m_stars.size(); i++)
        {
          typename Star::Facet_set_iterator fit = m_stars[i]->begin_restricted_facets();
          typename Star::Facet_set_iterator fitend = m_stars[i]->end_restricted_facets();
          for(; fit != fitend; fit++)
            facets.insert(Facet_ijk(*fit));          
        }
        return facets.size();
      }

      bool check_consistency_and_sliverity(Star_handle to_be_refined, //facet to be refined belongs to this star
                                           Star_handle new_star,     //the newly created star
                                           const Index_set& modified_stars,
                                           const double& sq_radius_bound,
                                           std::vector<int>& problematic_facets) const
      {
        return check_consistency(to_be_refined, new_star, modified_stars, sq_radius_bound, problematic_facets, true);
      }

      bool check_consistency(Star_handle to_be_refined, //facet to be refined belongs to this star
                             Star_handle new_star,     //the newly created star
                             const Index_set& modified_stars,
                             const double& sq_radius_bound,
                             std::vector<int>& problematic_facets,
                             const bool do_check_sliverity = false) const
      {
        //list all facets that would be created by p's insertion
        Point_3 p = new_star->center_point();
        int p_index = new_star->index_in_star_set();

        std::map<Facet_ijk, int> facets;
        std::vector<Index> bfacets_quads; //used for sliverity check, three point index  + star index
        facets_created(new_star, facets);

        typename Index_set::const_iterator it;
        for(it = modified_stars.begin(); it != modified_stars.end(); it++)
          facets_created(p, p_index, m_stars[(*it)], facets, bfacets_quads);

        typename std::map<Facet_ijk, int>::iterator itf;
        for(itf = facets.begin(); itf != facets.end(); itf++)
        {
          std::size_t nmax = m_stars.size();
          if( (*itf).second != 3) // the face is not there 3 times
          {
            if((*itf).first.is_infinite()) // should not happen
              return false;

            TPoint_3 tp0, tp1, tp2;
            TPoint_3 tp = transform_to_star_point(p, to_be_refined);

            tp0 = ((*itf).first.vertex(0) == nmax) ? tp 
              : transform_to_star_point(m_stars[(*itf).first.vertex(0)]->center_point(),to_be_refined);
            tp1 = ((*itf).first.vertex(1) == nmax) ? tp
              : transform_to_star_point(m_stars[(*itf).first.vertex(1)]->center_point(),to_be_refined);
            tp2 = ((*itf).first.vertex(2) == nmax) ? tp  
              : transform_to_star_point(m_stars[(*itf).first.vertex(2)]->center_point(),to_be_refined);

            double sqr = to_be_refined->compute_squared_circumradius(tp0, tp1, tp2);
            if(sqr < sq_radius_bound)
            {
              CGAL_PROFILER("[is_valid failure : small inconsistency]");
              // small inconsistency (radius) is forbidden. A big one is fine.
              for(int i=0; i<3; ++i)
                if( (*itf).first.vertex(i) != nmax )
                {
                  assert((*itf).first.vertex(i) < nmax);
                  problematic_facets.push_back( (*itf).first.vertex(i) );
                }
            }
          }
        }

        //consistency ok, check sliverity now
        if(do_check_sliverity && m_criteria->sliverity > 0.)
        {
          typename std::vector<Index>::iterator itp;
          for(itp = bfacets_quads.begin(); itp != bfacets_quads.end();)
          {
            CGAL::cpp11::array<TPoint_3, 4> ps;

            Star_handle quad_star = m_stars[*itp++];
            ps[0] = transform_to_star_point(p, quad_star);

            int index = 1;
            for(int i = 0; i < 3; ++i)
            {
              //std::cout << "index : " << index << " itp : " << *itp << std::endl;
              ps[index++] = transform_to_star_point(m_stars[*itp++]->center_point(), quad_star);
            }

            //todo : we should check sliverity in each star involved, not only in one star (?)
            if(quad_star->is_sliver(ps[0], ps[1], ps[2], ps[3]))
            {
              CGAL_PROFILER("[is_valid failure : sliverity]");
              return false;
            }
          }
        }

#ifdef ANISO_SLIVERITY_USE_BRUTE_FORCE
        if(do_check_sliverity && m_criteria->sliverity > 0.)
        {
          CGAL::cpp11::array<Point_3, 4> ps;

            //fill set with all points involved
          std::set<int> point_ids;
          for(it = modified_stars.begin(); it != modified_stars.end(); it++){
              point_ids.insert(*it);

              typename std::vector<Vertex_handle>::iterator nvi = m_stars[*it]->begin_neighboring_vertices();
              typename std::vector<Vertex_handle>::iterator nviend = m_stars[*it]->end_neighboring_vertices();
              for (; nvi != nviend; nvi++)
              {
                if (m_stars[*it]->is_infinite_vertex(*nvi))
                  continue;
                point_ids.insert((*nvi)->info());
              }
          }

            //pick three points from that set
          for(std::set<int>::iterator it1 = point_ids.begin(); it1 != point_ids.end(); ++it1)
          {
            for(std::set<int>::iterator it2 = it1; ++it2 != point_ids.end();)
            {
              for(std::set<int>::iterator it3 = it2; ++it3 != point_ids.end();)
              {
                  //build the 4 points and check sliverity for three points + p in each star metric
                //std::cout << "4 points : " << p_index << " " << *it1 << " " << *it2 << " " << *it3 << std::endl;
                assert(p_index != *it1 && p_index != *it2 && p_index != *it3 &&
                       *it1 != *it2 && *it1 != *it3 && it2 != *it3);

                  //using modified stars' metrics
                for(it = modified_stars.begin(); it != modified_stars.end(); it++)
                {
                  ps[0] = transform_to_star_point(p, m_stars[*it]);
                  ps[1] = transform_to_star_point(m_stars[*it1]->center_point(), m_stars[*it]);
                  ps[2] = transform_to_star_point(m_stars[*it2]->center_point(), m_stars[*it]);
                  ps[3] = transform_to_star_point(m_stars[*it3]->center_point(), m_stars[*it]);

                  if(m_stars[(*it)]->is_sliver(ps[0], ps[1], ps[2], ps[3]))
                  {
                    //std::cout << "failed for a modified star : " << *it << std::endl;
                    CGAL_PROFILER("[is_valid failure : sliverity1]");
                    return false;
                  }
                }

                  //using to_be_refined's metric
                ps[0] = transform_to_star_point(p, to_be_refined);
                ps[1] = transform_to_star_point(m_stars[*it1]->center_point(), to_be_refined);
                ps[2] = transform_to_star_point(m_stars[*it2]->center_point(), to_be_refined);
                ps[3] = transform_to_star_point(m_stars[*it3]->center_point(), to_be_refined);

                if(to_be_refined->is_sliver(ps[0], ps[1], ps[2], ps[3]))
                {
                  //std::cout << "failed for to_be_refined : " << to_be_refined->index_in_star_set() << std::endl;
                  CGAL_PROFILER("[is_valid failure : sliverity2]");
                  return false;
                }

                  //using new_star's metric
                ps[0] = transform_to_star_point(p, new_star);
                ps[1] = transform_to_star_point(m_stars[*it1]->center_point(), new_star);
                ps[2] = transform_to_star_point(m_stars[*it2]->center_point(), new_star);
                ps[3] = transform_to_star_point(m_stars[*it3]->center_point(), new_star);

                if(new_star->is_sliver(ps[0], ps[1], ps[2], ps[3]))
                {
                  //std::cout << "failed for new_star : " << new_star->index_in_star_set() << std::endl;
                  CGAL_PROFILER("[is_valid failure : sliverity3]");
                  return false;
                }

              }
            }
          }
        }
#endif
        if(problematic_facets.empty())
        {
          CGAL_PROFILER("[is_valid success]");
          return true;
        }
        else
          return false;
      }

      template<typename Map>
      void print_map(Map& m)
      {
        typename Map::iterator it;
        for(it = m.begin(); it != m.end(); it++)
        {
          std::cout << "\n\t[";
          std::cout << (*it).first.vertex(0) << " " << (*it).first.vertex(1) 
                << " " << (*it).first.vertex(2) << " \t- ";
          std::cout << (*it).second;
          std::cout << "]";
        }
      }

      void facets_created(Star_handle star,
                          std::map<Facet_ijk, int>& facets) const
      {
        typename Star::Facet_set_iterator it = star->begin_restricted_facets(); 
        typename Star::Facet_set_iterator iend = star->end_restricted_facets();
        for(; it != iend; it++)
          add_to_map(Facet_ijk(*it), facets);               
      }

      bool is_infinite_vertex(int i) const
      {
        return (i == Star::infinite_vertex_index());
      }

      // restricted facets(i,j,k) which would be created by the insertion 
      // of 'p' in 'star' are added to 'facets'
      void facets_created(const Point_3& p, // new point, not yet in the star set
                          const int p_id,   // index of p in star set
                          Star_handle star, 
                          std::map<Facet_ijk, int>& facets,
                          std::vector<Index>& bfacets_quads) const
      {
        int center_id = star->index_in_star_set();
        if(center_id == p_id)
          return facets_created(star, facets);
        
        // find boundary facets of the conflict zone
        std::vector<Facet> local_bfacets;
        int dim = star->dimension();
        // dimension 1
        if(dim == 1)
        {
          Edge_ij e(*(star->finite_edges_begin()));
          add_to_map(Facet_ijk(p_id, e.vertex(0), e.vertex(1)), facets);
          return;
        }

        std::set<Edge_ij> edges_around_c;
        // dimension 2
        if(dim == 2)
        {
          typename Star::Finite_edges_iterator eit = star->finite_edges_begin();
          typename Star::Finite_edges_iterator endit = star->finite_edges_begin();
          for(; eit != endit; eit++)
          {
            Edge_ij e(*eit);
            if(e.has(center_id))
              edges_around_c.insert(e);
          }
        }
        // dimension 3
        else if(dim == 3 && star->find_conflicts(p, std::back_inserter(local_bfacets)))
        { // get the circular set of edges/vertices around the center in 'boundary facets'
          get_cycle(local_bfacets, edges_around_c, center_id);
        }
        else std::cerr << "Warning : in 'facets_created', no conflict zone!\n";
        
        //add new facets to the list of sliverity test facets
        typename std::vector<Facet>::iterator it;
        for(it = local_bfacets.begin(); it != local_bfacets.end(); it++)
        {
          Index p1 = ((*it).first->vertex( ((*it).second + 1)%4 ))->info();
          Index p2 = ((*it).first->vertex( ((*it).second + 2)%4 ))->info();
          Index p3 = ((*it).first->vertex( ((*it).second + 3)%4 ))->info();

          if( p1 != star->infinite_vertex_index() &&
              p2 != star->infinite_vertex_index() &&
              p3 != star->infinite_vertex_index() )
          {
              //and center of the star to know the metric later on
            bfacets_quads.push_back( star->index_in_star_set() );
            assert(star->index_in_star_set() != -10);

              //three points for the facet
            bfacets_quads.push_back(p1);
            bfacets_quads.push_back(p2);
            bfacets_quads.push_back(p3);
          }
        }

        // get the set (circular if dim is 3) of vertices around the center
        std::set<int> vertices_around_c;
        get_vertices(edges_around_c, std::inserter(vertices_around_c, vertices_around_c.end()));

        // facets to be added are : (p, center_i, vertices of 'vertices_around_c')
        // add them only if they are restricted to the constrain surface
        typename std::set<int>::iterator iit;
        for(iit = vertices_around_c.begin(); iit != vertices_around_c.end(); iit++)
        {
          int v_id = *iit;
          if(is_infinite_vertex(v_id))
            continue; //it is an infinite vertex
                               
          Edge_ij e1, e2; // edges incident to vertex with 'v_id' in edges_around_c
          bool one_found = false;
          typename std::set<Edge_ij>::iterator eit;
          for(eit = edges_around_c.begin(); eit != edges_around_c.end(); eit++)
          {
            Edge_ij tmp = *eit;
            if(tmp.has(v_id)) //v_id is shared by e1 and e2
            {   
              if(!one_found)
              {
                e1 = tmp;
                one_found = true;
              }
              else
              {
                e2 = tmp;
                break;
              }
            }
          }

          //if(!e1.is_valid())
          //  std::cout << "Warning : e1 was not found!" << std::endl;
          if(!e2.is_valid())
          {
            std::cout << "Warning : e2 was not found!" << std::endl;
            continue;//todo : check this is correct
          }
          bool is_finite_e1 = !e1.is_infinite();
          bool is_finite_e2 = !e2.is_infinite();
          bool is_inside_c1 = false;
          bool is_inside_c2 = false;
          if(is_finite_e1)
          {
            Point_3 cc = compute_circumcenter(p, star->center_point(),
                                    m_stars[e1.vertex(0)]->center_point(),
                                    m_stars[e1.vertex(1)]->center_point(), star);
            is_inside_c1 = is_inside_domain(cc);
          }
          if(is_finite_e2)
          {
            Point_3 cc = compute_circumcenter(p, star->center_point(),
                                     m_stars[e2.vertex(0)]->center_point(),
                                     m_stars[e2.vertex(1)]->center_point(), star);
            is_inside_c2 = is_inside_domain(cc);
          }

          if((is_finite_e1 && is_finite_e2 && is_inside_c1 != is_inside_c2) //finite case
            || (!is_finite_e1 && is_finite_e2 && is_inside_c2)  // finite+infinite case
            || (!is_finite_e2 && is_finite_e1 && is_inside_c1)) // finite+infinite case
            add_to_map(Facet_ijk(center_id, p_id, v_id), facets);
        }
      }

      bool is_inside_domain(const Point_3& p) const
      {
        return (m_pConstrain->side_of_constraint(p) == ON_POSITIVE_SIDE);
      }

      template<typename FacetVector, typename EdgeSet>
      void get_cycle(const FacetVector& bfacets, 
                     EdgeSet& edges_around_c, 
                     const int& center_id) const
      {
        typename FacetVector::const_iterator fit = bfacets.begin();
        for( ; fit != bfacets.end(); fit++)
        {
          Facet_ijk f(*fit);
          int i1, i2;
          if(f.has(center_id, i1, i2))
            edges_around_c.insert(Edge_ij(i1, i2));
        }
        //debug_is_cycle(edges_around_c);
      }

      template<typename Edge_ijSet, typename IndexOutputIterator>
      void get_vertices(const Edge_ijSet& edges, IndexOutputIterator oit) const
      {
        typename Edge_ijSet::const_iterator eit;
        for(eit = edges.begin(); eit != edges.end(); eit++)
        {
          *oit++ = (*eit).vertex(0);
          *oit++ = (*eit).vertex(1);
        }
      }

      void debug_is_cycle(const std::set<Edge_ij>& edges)
      {
        std::map<int, int> m;
        typename std::set<Edge_ij>::const_iterator it;
        for(it = edges.begin(); it != edges.end(); it++)
        {
          int i = (*it).vertex(0);
          if(m.count(i) > 0) m.erase(i);
          else m[i] = 1;
          
          i = (*it).vertex(1);
          if(m.count(i) > 0) m.erase(i);
          else m[i] = 1;
        }
        if(!m.empty())
          std::cerr << "Error : edges_around_c is not a cycle!\n";
      }

      template<typename Element, typename OneMap>
      void add_to_map(const Element& e, OneMap& m) const
      {
        typename OneMap::iterator it = m.find(e);
        if(it != m.end())  
          m[e] = m[e] + 1;
        else 
          m.insert(std::pair<Element,int>(e,1));
      }

      Point_3 compute_circumcenter(const Facet& f, Star_handle here) const
      {
        return transform_from_star_point(here->compute_circumcenter(f), here);
      }
      Point_3 compute_circumcenter(Cell_handle cell, Star_handle here) const
      {
        return transform_from_star_point(here->compute_circumcenter(cell), here);
      }

      Point_3 compute_circumcenter(const Point_3& p0, const Point_3& p1, 
                                   const Point_3& p2, Star_handle here) const
      {
#ifdef ANISO_USE_CC_EXACT
        return back_from_exact(
           transform_from_star_point(here->compute_circumcenter(
                  transform_to_star_point(to_exact(p0), here),
                  transform_to_star_point(to_exact(p1), here),
                  transform_to_star_point(to_exact(p2), here)), here));
#else
        return transform_from_star_point(here->compute_circumcenter(
                  transform_to_star_point(p0, here),
                  transform_to_star_point(p1, here),
                  transform_to_star_point(p2, here)), here);
#endif
      }
      
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
      
      bool pick_valid(const Star_handle star, //to be refined
                      const Facet &facet,
                      Point_3& p,
                      const bool pick_valid_use_cube_probing = false) const //belongs to star and should be refined
      {
#ifdef USE_ANISO_TIMERS
        std::clock_t start_time = clock();
#endif
        // compute radius bound, in M_star
        Point_3 center;
#ifdef ANISO_USE_EXACT
        star->compute_exact_dual_intersection(facet, center);
#else
        star->compute_dual_intersection(facet, center);
#endif
        TPoint_3 tccf = transform_to_star_point(center, star);
        TPoint_3 tp2 = facet.first->vertex((facet.second + 1) % 4)->point();
        FT sq_circumradius = star->traits()->compute_squared_distance_3_object()(tccf, tp2);
        FT sq_radiusbound = m_criteria->beta * m_criteria->beta * sq_circumradius;

        // compute tools to pick random point in picking region
        FT circumradius = m_criteria->delta * std::sqrt(sq_circumradius);
        Cell_handle cell = (! star->is_infinite(facet.first)) ? facet.first
                                      : facet.first->neighbor(facet.second);
        // in M_star
#ifdef ANISO_USE_EXACT
        TPoint_3 tcccell = back_from_exact(star->compute_exact_circumcenter(cell));
#else
        TPoint_3 tcccell = star->compute_circumcenter(cell);
#endif
        typename Star::Traits::Compute_random_point_3 random = 
          star->traits()->compute_random_point_3_object();

        std::size_t tried_times = 1;
        bool success = false;

        pickvalid_problematic_facets.clear();
        m_pick_valid_cache.clear();
        m_pick_valid_facet.clear();
        m_pick_valid_facet.push_back(transform_from_star_point(facet.first->vertex((facet.second+1)%4)->point(), star));
        m_pick_valid_facet.push_back(transform_from_star_point(facet.first->vertex((facet.second+2)%4)->point(), star));
        m_pick_valid_facet.push_back(transform_from_star_point(facet.first->vertex((facet.second+3)%4)->point(), star));

        Star_handle newstar = new Star(m_criteria, m_pConstrain, true/*surface*/);

        //cube mode
        bool cube_probing = false;
        std::pair<int,int> A_n, B_n, C_n; // the three points forming ABC the triangle we want to study in cube mode

        // define the cube dimensions
        int cube_points_on_edge_n = 10;
        CGAL::Bbox_3 m_bbox = m_pConstrain->get_bbox();
        double cube_half_edge_size = 0.020*((std::max)((std::max)(m_bbox.xmax()-m_bbox.xmin(), m_bbox.ymax()-m_bbox.ymin()),m_bbox.zmax()-m_bbox.zmin()));

        // define the vector of points to be tested
        std::vector<TPoint_3> probing_points;
        probing_points.push_back(Point_3(1,1,1)); //dummy point

        while(true && !probing_points.empty())
        {
          if(cube_probing)
          {
            tried_times = 0;
            p = probing_points.back();
            probing_points.pop_back();
          }
          else
          {
            TPoint_3 next_point;
            next_point = random(tccf, tcccell, circumradius);
            p = star->compute_steiner_dual_intersection(next_point, facet);
          }

          std::vector<int> problematic_facets; //2*nb of facets : two points + p = one facet
          m_pick_valid_cache.push_back(p);

          if(is_valid_point(p, sq_radiusbound, star, newstar, problematic_facets))
          {
            success = true;
            CGAL_HISTOGRAM_PROFILER("[iterations for pick_valid]", tried_times);

            if(cube_probing)
              green_points.push_back(p);
            else
            {
              pickvalid_problematic_facets.clear();
              m_pick_valid_cache.clear();
              m_pick_valid_facet.clear();
              break;
            }
          }
          else
          {
            assert(!problematic_facets.empty());

            //check for 3 inconsistences within four points
            if(cube_probing)
            {
              //check which facets are involved in these inconsistences
              bool A_n_involved = false, B_n_involved = false, C_n_involved = false;
              bool others_involved = false;

              for(std::size_t i=0; i<problematic_facets.size();)
              {
                int i1 = problematic_facets[i];
                int i2 = problematic_facets[i+1];

                if(i1 == A_n.first && i2 == A_n.second)
                  A_n_involved = true;
                else if(i1 == B_n.first && i2 == B_n.second)
                  B_n_involved = true;
                else if(i1 == C_n.first && i2 == C_n.second)
                  C_n_involved = true;
                else
                  others_involved = true;

                i += 2;
              }

              if(A_n_involved && B_n_involved && C_n_involved && !others_involved)
                red_points.push_back(p);
              else if((A_n_involved || B_n_involved || C_n_involved) && others_involved)
                orange_points.push_back(p);
              else
                yellow_points.push_back(p);
            }
            else if(pick_valid_use_cube_probing &&
                    problematic_facets.size() == 6 && //three facets
                    pickvalid_problematic_facets.empty()) //first problematic facet (not really needed)
            {
              std::set<int> point_set(problematic_facets.begin(), problematic_facets.end());
              if(point_set.size() == 3)
              {
                //three facets & three points only are involved
                //enter cube_probing!
                cube_probing = true;
                p = compute_insert_or_snap_point(star, facet);

                //clear the containers up and deactivate tried_times
                 pickvalid_problematic_facets.clear();
                 m_pick_valid_cache.clear();
                 red_points.clear();
                 orange_points.clear();
                 yellow_points.clear();
                 green_points.clear();
                 tried_times = 0;

                //compute lower left back of the cube coordinates
                double x0 = p.x() - cube_half_edge_size;
                double y0 = p.y() - cube_half_edge_size;
                double z0 = p.z() - cube_half_edge_size;

                //fill probing points
                probing_points.clear(); //clear the dummy
                std::cout << "New cube :" << std::endl;
                double x, y, z;
                double step = 2*cube_half_edge_size/((double) cube_points_on_edge_n);
                for(int i=0; i<cube_points_on_edge_n; ++i) //x
                {
                  x = x0 + i*step;
                  for(int j=0; j<cube_points_on_edge_n; ++j) //y
                  {
                    y = y0 + j*step;
                    for(int k=0; k<cube_points_on_edge_n; ++k) //z
                    {
                      z = z0 + k*step;
                      probing_points.push_back(Point_3(x,y,z));
                    }
                  }
                }

                //fill A_n, B_n, C_n
                A_n = std::pair<int,int>(problematic_facets[0],problematic_facets[1]);
                B_n = std::pair<int,int>(problematic_facets[2],problematic_facets[3]);
                C_n = std::pair<int,int>(problematic_facets[4],problematic_facets[5]);

                std::cout << "entered cube probing in : " << star->index_in_star_set() << std::endl;
                std::cout << "ABC : " << A_n.first << " " << A_n.second << " " << B_n.first << " " << B_n.second << " " << C_n.first << " " << C_n.second << std::endl;
                std::cout << "The cube has " << probing_points.size() << " points to test" << std::endl;
                std::cout << "it is centered on : " << p << std::endl;
                std::cout << "points on edge : " << cube_points_on_edge_n << " half_edge : " << cube_half_edge_size;
                std::cout << " step : " << step << std::endl;
              }
            }
            pickvalid_problematic_facets[p] = problematic_facets;
          }
          if(++tried_times > m_criteria->max_times_to_try_in_picking_region || probing_points.empty())
          {
            p = compute_insert_or_snap_point(star, facet);
            CGAL_HISTOGRAM_PROFILER("[iterations for pick_valid]", tried_times);
            m_pick_valid_cache.push_back(p);

            if(cube_probing)
            {
              std::cout << red_points.size() << " red points" << std::endl;
              std::cout << orange_points.size() << " orange points" << std::endl;
              std::cout << yellow_points.size() << " yellow points" << std::endl;
              std::cout << green_points.size() << " green points" << std::endl;
              std::cout << "-------------------------------------" << std::endl;
            }

            break;
          }
        }
        delete newstar;
#ifdef USE_ANISO_TIMERS
        m_pick_valid_timer += duration(start_time);
#endif
        if(pick_valid_use_cube_probing)
          return !cube_probing;
        else
          return success;
      }

      bool exists(const Facet_ijk& f, unsigned int& star_id)
      {
        for(unsigned int i = 0; i < m_stars.size(); i++)
        {
          Star_handle s = m_stars[i];
          typename Star::Facet_set_iterator fit = s->begin_restricted_facets();
          typename Star::Facet_set_iterator fitend = s->end_restricted_facets();
          for(; fit != fitend; fit++)
          {
            Facet_ijk ftest(*fit);
            if(f == ftest)
            {
              star_id = i;
              return true;
            }
          }
        }
        return false;
      }

      Point_3 compute_insert_or_snap_point(const typename Constrain_surface::EdgeIterator &edge) 
      {
        Point_3 p = edge->first;
        Point_3 q = edge->second;
        Point_3 c((p.x() + q.x()) / 2.0, (p.y() + q.y()) / 2.0, (p.z() + q.z()) / 2.0);
        m_pConstrain->edge_split(edge, c);
        return c;
      }

      Point_3 compute_insert_or_snap_point(const Star_handle star, const Facet &facet) const
      {
        Point_3 p;
#ifdef ANISO_USE_EXACT
        star->compute_exact_dual_intersection(facet, p);
#else
        star->compute_dual_intersection(facet, p);
#endif
        // test encroachment
#ifdef CHECK_EDGE_ENCROACHMENT
        Star_iterator si = m_stars.begin();
        Star_iterator siend = m_stars.end();
        typename Constrain_surface::EdgeIterator encroached_edge;
        for (; si != siend; si++) {
          if ((*si)->is_encroached(p, encroached_edge))
            return compute_insert_or_snap_point(encroached_edge);
        }
#endif
        return p;
      }

public:
      // insert p to the star set
      int insert(const Point_3& p, const bool conditional)
      {
        Index_set modified_stars;
        return insert(p, modified_stars, conditional);
      }

private:
      Index insert(const Point_3& p,
                   Index_set& modified_stars, // should be empty, except in special cases
                   const bool conditional,
                   const bool smoothing = false)
      {         
#ifdef USE_ANISO_TIMERS
        std::clock_t start_time = clock();
#endif
        Index id = insert_to_stars(p, modified_stars, conditional); //performs insertion in CH
        
        if(id < 0 || id < (int)m_stars.size())
          return id;
        
        Star_handle star = create_star(p, id, modified_stars, smoothing);

        if(star->index_in_star_set() != m_stars.size())
          std::cout << "WARNING in insert..." << std::endl;

        m_stars.push_back(star);
        modified_stars.insert(star->index_in_star_set());
        m_kd_tree.insert(star->index_in_star_set());
#ifndef NO_USE_AABB_TREE_OF_BBOXES
        m_aabb_tree.insert(AABB_primitive(star));
#endif
#ifdef USE_ANISO_TIMERS
        m_insert_timer += duration(start_time);
#endif
        return id;
      }

//      Index insert(const Point_3& p,  
//                   Index_set& modified_stars, // should not be empty
//                   Star_handle new_star,
//                   const bool conditional)
//      {
//#ifdef USE_ANISO_TIMERS
//        std::clock_t start_time = clock();
//#endif
//        Star_set target_stars = get_stars(modified_stars.begin(), modified_stars.end());
//        Index id = perform_insertions(p, new_star->index_in_star_set(),
//                                      target_stars, modified_stars, conditional);
//        if(id != new_star->index_in_star_set())
//          return id;
//        
//        if(new_star->index_in_star_set() != m_stars.size())
//          std::cout << "WARNING in insert..." << std::endl;
//
//        m_stars.push_back(new_star);             
//        modified_stars.insert(new_star->index_in_star_set());
//        m_kd_tree.insert(new_star->index_in_star_set());
//#ifndef NO_USE_AABB_TREE_OF_BBOXES
//        m_aabb_tree.insert(AABB_primitive(new_star));
//#endif
//#ifdef USE_ANISO_TIMERS
//        m_insert_timer += duration(start_time);
//#endif
//        return id;
//      }

      Index insert_in_domain(const Point_3& p, const bool conditional)
      {
#ifdef USE_ANISO_TIMERS
        std::clock_t start_time = clock();
#endif
        Index_set modified_stars;
        Index id = insert_to_stars(p, modified_stars, conditional); //performs insertion in CH
        
        if(id < 0 || id < (int)m_stars.size())
          return id;
        
        Star_handle star = create_inside_star(p, id, modified_stars);

        if(star->index_in_star_set() != m_stars.size())
          std::cout << "WARNING in insert_in_domain..." << std::endl;

        m_stars.push_back(star);    
        modified_stars.insert(star->index_in_star_set());
        m_kd_tree.insert(star->index_in_star_set());
#ifndef NO_USE_AABB_TREE_OF_BBOXES
        m_aabb_tree.insert(AABB_primitive(star));
#endif
#ifdef USE_ANISO_TIMERS
        m_insert_timer += duration(start_time);
#endif
        return id;
      }

public:
      void print_stars() const
      {
        std::cout << "Star set :("<< m_stars.size() <<") \n";
        for(unsigned int i = 0; i < m_stars.size(); i++)
          m_stars[i]->print_vertices();
      }

      //void print_ch_triangulation()
      //{
      //  unsigned int N = m_ch_triangulation.number_of_vertices();
      //  std::cout << "CH Triangulation : " << N << " v.)\n";
      //  typename DT::Finite_vertices_iterator v = m_ch_triangulation.finite_vertices_begin();
      //  typename DT::Finite_vertices_iterator vend = m_ch_triangulation.finite_vertices_end();
      //  for(; v != vend; v++)
      //  {
      //    std::vector<Vertex_handle> inc;
      //    m_ch_triangulation.incident_vertices(v, std::back_inserter(inc));
      //    std::cout << "\t Star of p" << v->info() << " : [";
      //    for(unsigned int i = 0; i < inc.size(); i++)
      //    {
      //      std::cout << inc[i]->info() << " ";
      //    }
      //    std::cout << "]\n";
      //  }
      //}

      //void debug_check_topological_disks()
      //{
      //  std::cout << "Check topological disks...";
      //  bool all_ok = true;   
      //  for(unsigned int i = 0; i < m_stars.size(); i++)
      //  {
      //    if(!m_stars[i]->is_topological_disk())
      //    {
      //      std::cout << i << " ";
      //      all_ok = false;
      //    }
      //  }
      //  if(all_ok) std::cout << "Ok." << std::endl;
      //  else std::cout << " are not topological disks." << std::endl;
      //}

      int simulate_insert_to_stars(const Point_3& p,
                                   Index_set& modified_stars) const
      {
        Index this_id = static_cast<Index>(m_stars.size());
        
        // find conflicted stars
        Star_set stars;
        finite_stars_in_conflict(p, std::inserter(stars, stars.end())); //aabb tree
        infinite_stars_in_conflict(p, std::inserter(stars, stars.end()));//convex hull

        typename Star_set::iterator it = stars.begin();
        typename Star_set::iterator itend = stars.end();
        for(; it != itend; it++) 
        {
          Star_handle si = *it;
          int id = si->simulate_insert_to_star(p, this_id);
          
          if(id == -1)          // no conflict
            continue;
          else if(id < (int)this_id) // already in star set
            return id;
          else                  // to be inserted, std configuration
            modified_stars.insert(si->index_in_star_set());
        }
        return this_id;
      }
     
      // insert p inside existing stars
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
                
        return id;
      }

      // inserts p in target_stars ('conditional' to be in conflict or not)
      // modified_stars get filled with the stars in which p is actually inserted
      // triangulation of CH is updated accordingly
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

      //vertex with index 'id' should be removed from stars in [begin, end)
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
       
      void update_bboxes() const
      {
        std::size_t i;
        std::size_t N = m_stars.size();
        for(i = 0; i < N; i++)
          m_stars[i]->update_bbox();
      }

      template<typename OutputIterator>
      void finite_stars_in_conflict(const Point_3& p,
                                    OutputIterator oit) const
      {
//        update_bboxes();

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

      //bool infinite_stars_in_conflict(const Point_3& p) const
      //{
      //  std::vector<Star_handle> vec;
      //  return infinite_stars_in_conflict(p, std::back_inserter(vec), false);
      //}      

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

      //void debug_intersected_stars(const Point_3& p)
      //{
      //  std::set<Star_handle> with;
      //  intersected_stars(p, std::inserter(with, with.begin()), true);
      //  std::cout << "\n\tAABB : " << with.size() << " stars.";

      //  // check that p is inside all these bboxes
      //  typename std::set<Star_handle>::iterator it;
      //  for(it = with.begin(); it != with.end(); it++)
      //  {
      //    if(! (*it)->is_inside_bbox(p))
      //      std::cout << "-ERROR 1-";
      //  }
      //  std::cout << std::endl;

      //  // check that all stars of exhaustive search are in aabb tree search
      //  std::set<Star_handle> without;
      //  intersected_stars(p, std::inserter(without, without.begin()), false);
      //  std::cout << "\tnoAABB : " << without.size() << " stars.";
      //  for(it = without.begin(); it != without.end(); it++)
      //  {
      //    if(with.find(*it) == with.end())
      //      std::cout << "-ERROR 2-";
      //  }
      //  std::cout << std::endl;
      //}

      template<typename Stars>
      void fill_refinement_queue(const Stars& modified_stars, 
                                 int relative_point = -1) 
      {
#ifdef USE_ANISO_TIMERS
        std::clock_t start_time = clock();
#endif
        typename Stars::const_iterator si = modified_stars.begin();
        typename Stars::const_iterator siend = modified_stars.end();
        for (; si != siend; si++) 
        {
          Star_handle star = get_star(si);
          if(!star->is_surface_star())
            continue;

          Facet_set_iterator fi = star->begin_restricted_facets();
          Facet_set_iterator fiend = star->end_restricted_facets();
          for (; fi != fiend; fi++) 
          {
            Point_3 cc;
            star->compute_dual_intersection(*fi, cc);
            if(!m_refinement_condition(cc))
              continue;
            
            Cell_handle cell = fi->first;
            int offset = fi->second;
            if (relative_point >= 0) // we don't consider not-relative facets
            {
              bool relative = false;
              for (int i = 1; i <= 3; i++)
                if (relative_point == cell->vertex((offset + i) % 4)->info()) 
                {
                  relative = true;
                  break;
                }
              if (!relative)
                continue;
            }

            // over distortion : 1
            bool b_continue = false;
            if(m_criteria->distortion > 0.)
            {
              for (int i = 0; i < 3; i++) 
              {
                int index_1 = (offset + i + 1) % 4;
                int index_2 = (offset + (i + 1) % 3 + 1) % 4;
                FT over_distortion =
                  m_stars[cell->vertex(index_1)->info()]->metric().compute_distortion(
                  m_stars[cell->vertex(index_2)->info()]->metric()) - m_criteria->distortion;
                if (over_distortion > 0) 
                { // here, protect the edge
  #ifdef ANISO_DEBUG_REFINEMENT
                  Index im1 = cell->vertex(index_1)->info();
                  Index im2 = cell->vertex(index_2)->info();
                  typename Star::Metric m1 = m_stars[im1]->metric();
                  typename Star::Metric m2 = m_stars[im2]->metric();
  #endif

                  m_refine_queue.push_over_distortion(star, *fi, over_distortion);
                  b_continue = true;
                  break;
                }
              }
              if(b_continue) continue;
            }
            // too big : 2
            if(m_criteria->circumradius > 0.)
            {
              FT over_circumradius = star->compute_circumradius_overflow(*fi);
              if (over_circumradius > 0) 
              {
                m_refine_queue.push_over_circumradius(star, *fi, over_circumradius);
                continue;
              }
            }
            // bad shape : 3
            if(m_criteria->radius_edge_ratio > 0.)
            {
              FT over_radius_edge_ratio = star->compute_radius_edge_ratio_overflow(*fi);
              if (over_radius_edge_ratio > 0) 
              {
                m_refine_queue.push_bad_shape(star, *fi, over_radius_edge_ratio);
                continue;
              }
            }
            // bad approx : 4
            if(m_criteria->approximation > 0.)
            {
              FT over_approx = std::sqrt(sq_distance_to_surface(*fi, star)) - m_criteria->approximation;
              if(over_approx > 0.)
              {
                m_refine_queue.push_bad_approximation(star, *fi, over_approx);
                continue;
              }
            }

            // inconsistency : 5
            b_continue = false;
            for (int i = 1; i <= 3; i++) 
            {
              int vi = cell->vertex((offset + i) % 4)->info();
              if (star->index_in_star_set() == vi)
                continue;
              if (!m_stars[vi]->has_facet(*fi)) 
              {
                m_refine_queue.push_inconsistent(star, *fi, star->compute_volume(*fi));
                b_continue = true;
                break;
              }
            }
            if(b_continue) continue; //useless but safer if we add another criterion       
          } // facet
        } // star
#ifdef USE_ANISO_TIMERS
        m_fill_queue_timer += duration(start_time);
#endif
      }

      bool is_consistent(const Facet& f,
                         const bool verbose = false) const
      {
        bool retval = true;
        bool facet_told = false;
        for(int i = 1; i < 4; i++)
        {
          int index = f.first->vertex((f.second + i) % 4)->info();
          if(!m_stars[index]->has_facet(f))
          {
            if(verbose)
            {
              if(!facet_told)
              {
                std::cout << "f(" << f.first->vertex((f.second+1)%4)->info()
                          << " "  << f.first->vertex((f.second+2)%4)->info()
                          << " "  << f.first->vertex((f.second+3)%4)->info() << ") inconsistent : ";
                facet_told = true;
              }
              std::cout << "f not in S_" << index << ", ";
            }
            retval = false;
          }
        }
        if(verbose && !retval) std::cout << "." << std::endl;
        return retval;
      }

      bool is_consistent(const bool verbose = false) const
      {
        std::size_t N = m_stars.size();
        for(std::size_t i = 0; i < N; i++)
        {
          Star_handle star = m_stars[i];
          typename Star::Facet_set_iterator fit = star->begin_restricted_facets();
          typename Star::Facet_set_iterator fitend = star->end_restricted_facets();
          for(; fit != fitend; fit++)
            if(!is_consistent(*fit, verbose))
              return false;
        }
        return true;
      }

      void fill_refinement_queue()
      {
        m_refine_queue.clear();
        fill_refinement_queue(m_stars, -1);
#ifdef ANISO_VERBOSE
        m_refine_queue.print();
#endif
      }

      void initialize_stars(const int nb = 8) 
      {
#ifdef ANISO_VERBOSE
        std::cout << "Initialize "<< nb << " stars..." << std::endl;
#endif
        //typename Constrain_surface::Pointset initial_points = m_pConstrain->initial_points(nb);
        initial_points.clear();
        initial_points = m_pConstrain->get_surface_points(2*nb);

        initialize_medial_axis(); // poles

        std::cout << "now inserting initial points" << std::endl;

        typename Constrain_surface::Pointset::iterator pi = initial_points.begin();
        typename Constrain_surface::Pointset::iterator pend = initial_points.end();
#ifdef ANISO_VERBOSE
        std::cout << "(" << initial_points.size() << " initial points) ";
#endif
        int nbdone = 0;
        for (; pi != pend && nbdone < nb; pi++)
        {
          if(nbdone % 100 == 0)
            this->clean_stars();

          std::size_t this_id = m_stars.size();
          int id = -1;
          //if(m_refinement_condition(*pi))
            id = insert(*pi, false/*under no condition*/);
          if(this_id == id)
            nbdone++;
        }
#ifdef ANISO_VERBOSE
        std::cout << "done (" << m_stars.size() << " stars)." << std::endl;
#endif
      }

      void initialize_medial_axis()
      {
#ifdef ANISO_VERBOSE
        std::cout << "Initialize medial axis...";
        // compute poles (1/2) (1 per vertex : the furthest)
        // compute poles (2/2) (1 per vertex : the furthest on the other side)
        std::cout << "(compute poles...";
#endif
        poles.clear();
        poles = m_pConstrain->compute_poles();
#ifdef ANISO_VERBOSE
        std::cout << poles.size() << ")" << std::endl;

        //insert them all in all stars
        std::cout << "(insert poles...";
#endif
        unsigned int i = 1;
        unsigned int done = 0;
        typename std::set<Point_3>::const_iterator it;
        for(it = poles.begin(); it != poles.end(); ++it, ++i)
        {
          if(i % 500 == 0)
            this->clean_stars();

          bool conditional = (i % 10 != 0); //let's put 1/10 with no condition
                                            //(poles are sorted so geometry should be covered)
          //if(m_refinement_condition(*it))
          //{
            insert_in_domain(*it, conditional);
            ++done;
          //}
        }
        this->clean_stars();
#ifdef ANISO_VERBOSE
        std::cout << ") done ("<< done << " points)." << std::endl;
#endif
      }

      std::size_t number_of_stars() const
      {
        return m_stars.size();
      }

      unsigned int number_of_surface_stars() const
      {
        unsigned int nb = 0;
        for(unsigned int i = 0; i < m_stars.size(); ++i)
          if(m_stars[i]->is_surface_star())
            nb++;
        return nb;
      }

      std::size_t total_number_of_vertices() const
      {
        std::size_t nbv = 0;
        for(unsigned int i = 0; i < m_stars.size(); i++)
          nbv += m_stars[i]->number_of_vertices();
        return nbv;
      }
            
      void build_aabb_tree()
      {
#ifdef ANISO_VERBOSE
        std::cout << "Build the aabb tree...";
#endif
        m_aabb_tree.rebuild(m_stars.begin(), m_stars.end());
#ifdef ANISO_VERBOSE
        std::cout << " done." << std::endl;
#endif
      }

      bool next_refine_cell(Refine_facet &refine_facet, 
                            Facet &facet, 
                            bool &need_picking_valid, 
                            int &queue_type)
      {
        while(true) 
        {
          if (!m_refine_queue.top(refine_facet, queue_type))
            return false;
          m_refine_queue.pop();
          if (refine_facet.star->has_facet_ref(refine_facet.vertices, facet))
          {
            need_picking_valid = m_refine_queue.need_picking_valid(queue_type);
            return true;
          }
        }
      }

      Metric build_smoothed_metric(const Point_3& p, const Index_set& modified_stars) const
      {
        typename Index_set::const_iterator si = modified_stars.begin();
        typename Index_set::const_iterator siend = modified_stars.end();

        Star_handle star_i = get_star(si++);

        std::cout << "first star is : " << star_i->index_in_star_set() << std::endl;

        Metric metric_i = star_i->metric();
        Metric metric_p = m_metric_field->scale_metric_to_point(metric_i, star_i->center_point(), p);

        for (; si != siend; si++)
        {
          star_i = get_star(si);
          std::cout << "star : " << star_i->index_in_star_set() << "---------------------------------------------" << std::endl;

          metricA = metric_p;
          metric_i = star_i->metric();
          Metric scaled_metric_i = m_metric_field->scale_metric_to_point(metric_i, star_i->center_point(), p);
          metricB = star_i->metric();
          metric_p = m_metric_field->intersection(metric_p, scaled_metric_i);
          metricP = metric_p;
          std::cout << "intersected in buildsmooth : " << std::endl <<  metric_p.get_transformation() << std::endl;
        }

        Metric theoritical_at_p = m_metric_field->compute_metric(p);
        std::cout << "final metric at p : " << std::endl << metric_p.get_transformation() << std::endl << metric_p.get_inverse_transformation() << std::endl;
        std::cout << "theoritical at p was : " << std::endl << theoritical_at_p.get_transformation() << std::endl << theoritical_at_p.get_inverse_transformation() << std::endl;

        return metric_p;
      }

      Star_handle create_star(const Point_3 &p, 
                              int pid,
                              const Index_set& modified_stars,
                              const bool smoothing = false) const //by p's insertion
      {
        Star_handle star = new Star(m_criteria, m_pConstrain, true/*surface*/);
        create_star(p, pid, modified_stars, star, true/*surface*/, smoothing);
        return star;
      }

      Star_handle create_inside_star(const Point_3 &p, 
                                     int pid,
                                     const Index_set& modified_stars) const//by p's insertion
      {
        Star_handle star = new Star(m_criteria, m_pConstrain, false/*surface*/);
        create_star(p, pid, modified_stars, star, false/*surface*/);
        return star;
      }

      void create_star(const Point_3 &p, 
                       int pid,
                       const Index_set& modified_stars,//by p's insertion
                       Star_handle& star,
                       const bool surface_star = true, //o.w. center is inside volume
                       const bool smoothing = false) const
      {
#ifdef USE_ANISO_TIMERS
        std::clock_t start_time = clock();
#endif
        if(smoothing)
        {
          //check for poles in modified_stars
          Index_set modified_stars_without_poles;
          typename Index_set::const_iterator si = modified_stars.begin();
          typename Index_set::const_iterator siend = modified_stars.end();
          for (; si != siend; si++)
            if((get_star(si))->is_surface_star())
              modified_stars_without_poles.insert(*si);

          if(modified_stars_without_poles.empty())
          {
            std::cout << "can't smooth with an empty modified_stars_without_poles. ";
            std::cout << "Using default metric instead." << std::endl;
            std::cout << "modified_stars was empty : " << modified_stars.empty() << std::endl;
            star->reset(p, pid, m_metric_field->compute_metric(p), surface_star);
          }
          else
          {
            std::cout << "modified stars without poles : " << modified_stars_without_poles.size() << std::endl;
            std::cout << "modified stars : " << modified_stars.size() << std::endl;
            std::cout << "point P is : " << p << std::endl;
            Metric M = build_smoothed_metric(p, modified_stars_without_poles);
            star->reset(p, pid, M, surface_star);
            star->red_ellipsoid = true;
          }
        }
        else if(surface_star)
          star->reset(p, pid, m_metric_field->compute_metric(p), surface_star);
        else
          star->reset(p, pid, m_metric_field->uniform_metric(p), surface_star);

        typename Index_set::const_iterator si = modified_stars.begin();
        typename Index_set::const_iterator siend = modified_stars.end();
        for (; si != siend; si++) 
        {
          Star_handle star_i = get_star(si);
          star->insert_to_star(star_i->center_point(), star_i->index_in_star_set(), false);
            //no condition because they should be there for consistency
        }

        typename Star::Bbox bbox = star->bbox(); // volume bbox when star is not a topo_disk
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


#ifdef ANISO_DEBUG
        typename Star::Facet_set_iterator it = star->begin_restricted_facets();
        typename Star::Facet_set_iterator itend = star->end_restricted_facets();
        for(; it != itend; it++)
        {
          Facet f = *it;
          Vertex_handle v1 = f.first->vertex((f.second+1)%4);
          Vertex_handle v2 = f.first->vertex((f.second+2)%4);
          Vertex_handle v3 = f.first->vertex((f.second+3)%4);

          double deg_value = 1e-4;
          bool degenerated = ( (std::abs(v1->point().x()-v2->point().x()) < deg_value &&
                                std::abs(v1->point().y()-v2->point().y()) < deg_value &&
                                std::abs(v1->point().z()-v2->point().z()) < deg_value ) ||
                               (std::abs(v2->point().x()-v3->point().x()) < deg_value &&
                                std::abs(v2->point().y()-v3->point().y()) < deg_value &&
                                std::abs(v2->point().z()-v3->point().z()) < deg_value ) ||
                               (std::abs(v1->point().x()-v3->point().x()) < deg_value &&
                                std::abs(v1->point().y()-v3->point().y()) < deg_value &&
                                std::abs(v1->point().z()-v3->point().z()) < deg_value ) );

          if(degenerated){
              std::cout.precision(15);
              std::cout << "Building bad facet : " << pid << std::endl;
              std::cout << "\tp1 : " << v1->point() << std::endl;
              std::cout << "\tp2 : " << v2->point() << std::endl;
              std::cout << "\tp3 : " << v3->point() << std::endl;
              std::cout << "\tp was : " << p << std::endl;
          }
        }
#endif

#ifdef USE_ANISO_TIMERS
        m_create_star_timer += duration(start_time);
#endif
      }

      std::set<Facet_ijk> restricted_facets(Star_handle star)
      {
        std::set<Facet_ijk> facets;
        typename Star::Facet_set_iterator it = star->begin_restricted_facets();
        typename Star::Facet_set_iterator itend = star->end_restricted_facets();
        for(; it != itend; it++)
        {
          Facet f = *it;
          facets.insert(Facet_ijk(f));
        }
        return facets;
      }

      bool pick_valid_output(const bool need_picking_valid,
                             int & pick_valid_succeeded,
                             int & pick_valid_failed,
                             const bool pick_valid_causes_stop,
                             const int pick_valid_max_failures,
                             const bool success)
      {
        if(!success)
          pick_valid_failed++;
        else if(need_picking_valid)
          pick_valid_succeeded++;

        if((!success && pick_valid_failed % 100 == 0 && pick_valid_failed > 0) ||
           (success && need_picking_valid && pick_valid_succeeded % 100 == 0 && pick_valid_succeeded > 0))
        {
          std::cout << "pick_valid : ";
          std::cout << pick_valid_succeeded << " success and ";
          std::cout << pick_valid_failed << " failures" << std::endl;
        }
        if(pick_valid_causes_stop && pick_valid_failed >= pick_valid_max_failures)
        {
          //print the problematic facet map
          assert(m_pick_valid_cache.size() == pickvalid_problematic_facets.size());
          for(std::size_t i = 0; i<m_pick_valid_cache.size(); ++i)
          {
            Point_3 pickvalid_point = m_pick_valid_cache[i];
            std::cout << "id : " << i << " for the point " << pickvalid_point;
            std::cout << " , " << pickvalid_problematic_facets[pickvalid_point].size()/2 << " facets in play" << std::endl;

            std::vector<int> id_vector = pickvalid_problematic_facets[pickvalid_point];

            assert(id_vector.size()%2 == 0);
            for(std::size_t i=0; i<id_vector.size();)
            {
              std::cout << "(" << id_vector[i] << " " << id_vector[i+1] << "), ";
              i+=2;
            }
            std::cout << std::endl;
          }
          return false;
        }
        return true;
      }

      bool refine(int & pick_valid_succeeded,
                  int & pick_valid_failed,
                  const bool pick_valid_causes_stop = false,
                  const int pick_valid_max_failures = 100,
                  const bool pick_valid_use_cube_probing = false,
                  const bool metric_smoothing = false)
      {
        int queue_type = 0;
        Refine_facet bad_facet;
        bool need_picking_valid;
        Facet f; // to be refined
        if (!next_refine_cell(bad_facet, f, need_picking_valid, queue_type))
          return false;

        Vertex_handle v1 = f.first->vertex((f.second+1)%4);
        Vertex_handle v2 = f.first->vertex((f.second+2)%4);
        Vertex_handle v3 = f.first->vertex((f.second+3)%4);

#ifdef ANISO_DEBUG
        double deg_value = 1e-4;
        bool degenerated = ( (std::abs(v1->point().x()-v2->point().x()) < deg_value &&
                              std::abs(v1->point().y()-v2->point().y()) < deg_value &&
                              std::abs(v1->point().z()-v2->point().z()) < deg_value ) ||
                             (std::abs(v2->point().x()-v3->point().x()) < deg_value &&
                              std::abs(v2->point().y()-v3->point().y()) < deg_value &&
                              std::abs(v2->point().z()-v3->point().z()) < deg_value ) ||
                             (std::abs(v1->point().x()-v3->point().x()) < deg_value &&
                              std::abs(v1->point().y()-v3->point().y()) < deg_value &&
                              std::abs(v1->point().z()-v3->point().z()) < deg_value) );

        if(degenerated){
            std::cout << "trying to refine a degenerate facet" << std::endl;
        }
#endif

        Index_set modified_stars;// the ones that would be modified by p's insertion
        if (queue_type == 0) // encroachment
        { 
          std::cerr << "Error : encroachment is not implemented.\n";
          vertex_without_picking_count++;
          return true; //note false would stop refinement
        } 
         
        Point_3 steiner_point;
        bool success = compute_steiner_point(bad_facet.star, f,
                                             need_picking_valid, steiner_point,
                                             pick_valid_use_cube_probing);

        if(!pick_valid_output(need_picking_valid, pick_valid_succeeded, pick_valid_failed,
                              pick_valid_causes_stop, pick_valid_max_failures, success))
          return false;

#ifdef ANISO_DEBUG_REFINEMENT
        if(bad_facet.star->debug_steiner_point(steiner_point, f))
          std::cerr << "(Not exact)" << std::endl;
#endif
        if(!m_refinement_condition(steiner_point))
          return true; //note false would stop refinement

        //check if the facet trying to be refined is too small + success = false => enter metric smoothing
        bool smoothing = false;
        CGAL::Bbox_3 m_bbox = m_pConstrain->get_bbox();
        double max_sq_circumradius = 0.001*((std::max)((std::max)(m_bbox.xmax()-m_bbox.xmin(), m_bbox.ymax()-m_bbox.ymin()),m_bbox.zmax()-m_bbox.zmin()));
        max_sq_circumradius = max_sq_circumradius*max_sq_circumradius;

        if(metric_smoothing && !success && (bad_facet.star)->compute_squared_circumradius(f) < max_sq_circumradius)
        {
          std::cout << "smooth mode : " << (bad_facet.star)->compute_squared_circumradius(f) << " allowed : " << max_sq_circumradius << std::endl;
          smoothing = true;
          vertex_with_smoothing_counter++;
        }
        else if(metric_smoothing && !success)
        {
          std::cout << "compute_squared_circumradius : " << (bad_facet.star)->compute_squared_circumradius(f) << " allowed : " << max_sq_circumradius << std::endl;
          vertex_without_smoothing_counter++;
        }

        int pid = insert(steiner_point, modified_stars, true/*conditional*/, smoothing);

        //check if f has been destroyed
        // begin debug
        Cell_handle c;
        int i,j,k;
        if(bad_facet.star->is_facet(v1, v2, v3, c, i, j, k))
        {
          int index = 6 - i - j - k;
          Facet ff = bad_facet.star->make_canonical(Facet(c,index));
          if(bad_facet.star->is_restricted(ff))
          {
            modified_stars.clear();
            pop_back_star();
            if(bad_facet.star->is_facet(v1, v2, v3, c, i, j, k))
            {
              index = 6 - i - j - k;
              ff = bad_facet.star->make_canonical(Facet(c,index));

              steiner_point = compute_exact_steiner_point(bad_facet.star, ff, need_picking_valid);
              pid = insert(steiner_point, modified_stars, true/*conditional*/);
            }

#ifdef ANISO_DEBUG_REFINEMENT
            if(bad_facet.star->is_facet(v1, v2, v3, c, i, j, k))
            {
              index = 6 - i - j - k;
              ff = bad_facet.star->make_canonical(Facet(c,index));
              if(bad_facet.star->is_restricted(ff))
              {
                std::cout.precision(15);
                std::cout << "Bad facet still in : " << bad_facet.star << std::endl;
                std::cout << "\tp"<< v1->info() <<" : " << v1->point() << std::endl;
                std::cout << "\tp"<< v2->info() <<" : " << v2->point() << std::endl;
                std::cout << "\tp"<< v3->info() <<" : " << v3->point() << std::endl;
             
                std::cout << ", dim = " << bad_facet.star->dimension();
                std::cout<< ", nbv = " << bad_facet.star->number_of_vertices();
                std::cout << ", pid = " << pid;
                std::cout << ", p = " << steiner_point;
                std::cout << ", f(p) = " << m_pConstrain->side_of_constraint(steiner_point);
                std::cout << std::endl;
                bad_facet.star->debug_steiner_point(steiner_point, ff);
                std::cout << "(Exact)" << std::endl;
              }
              else
                std::cout << "bad facet still in but not restricted anymore" << std::endl;
            }
            else
              std::cout << "Switching to exact succeeded" << std::endl;
#endif
          }
        }
        else
        {
          //std::cout << "Facet " << v1->info() << " " << v2->info() << " " << v3->info() << " was dealt with" << std::endl;
        }
        //end debug

        if(!modified_stars.empty())
          fill_refinement_queue(modified_stars, pid);
        return true;
      }

      void pop_back_star()
      {
        Star_handle last = m_stars.back();
        delete last;

        m_stars.pop_back();
        Index id = static_cast<Index>(m_stars.size());
        remove_from_stars(id, m_stars.begin(), m_stars.end());
        m_kd_tree.remove_last();
        m_aabb_tree.remove_last();
      }

      bool compute_steiner_point(Star_handle to_be_refined,
                                 const Facet& f, //facet to be refined
                                 const bool need_picking_valid,
                                 Point_3& steiner,
                                 const bool pick_valid_use_cube_probing = false) const
      {
#ifdef USE_ANISO_TIMERS
        std::clock_t start_time = clock();
#endif
        bool success = true;
        if (need_picking_valid) 
        {            
          vertex_with_picking_count++;
          success = pick_valid(to_be_refined, f, steiner, pick_valid_use_cube_probing);
        } 
        else 
        {
          vertex_without_picking_count++;
          steiner = compute_insert_or_snap_point(to_be_refined, f);
        }
#ifdef USE_ANISO_TIMERS
        m_compute_steiner_timer += duration(start_time);
#endif
        return success;
      }

      Point_3 compute_exact_steiner_point(Star_handle to_be_refined,
                                    const Facet& f, //facet to be refined
                                    const bool need_picking_valid) const
      {
#ifdef USE_ANISO_TIMERS
        std::clock_t start_time = clock();
#endif
        Point_3 p;
        if (need_picking_valid) 
        {            
          vertex_with_picking_count++;
          pick_valid(to_be_refined, f, p);
        } 
        else 
        {
          vertex_without_picking_count++;
          to_be_refined->compute_exact_dual_intersection(f, p);
        }
#ifdef USE_ANISO_TIMERS
        m_compute_steiner_timer += duration(start_time);
#endif
        return p;
      }

      void dump() 
      {
        std::ostringstream nbs;
        nbs << (m_stars.size());
        
        std::string title("dump.off");        
        std::cout << "Write " << title << " (" << nbs.str() << " vertices)...";

        std::ofstream fx(title.c_str());
        fx << "OFF" << std::endl;
        fx << (m_stars.size()) << " 0 0" << std::endl;
        for (int i = 0; i < (int)m_stars.size(); i++)
          fx << m_stars[i]->center_point() << std::endl;
  
        std::cout << "done.\n";
        fx.close();
        //fx << 0 << std::endl;

        //for (int i = 0; i < (int)m_stars.size(); i++) 
        //{
        //  Star_handle star = m_stars[i];
        //  Star::Vertex_iterator vi = star->vertices_begin();
        //  Star::Vertex_iterator viend = star->vertices_end();
        //  for (; vi != viend; vi++) 
        //  {
        //    if (star->is_infinite(vi))
        //      continue;
        //    fx << (*vi).info() << " ";
        //  }
        //  fx << -1 << std::endl;
        //}
      }

      void report()
      {
        typename std::ofstream fx("report_surface.txt");
        fx << "[Parameters]" << std::endl << std::endl;
        m_criteria->report(fx);
        fx << std::endl << "[Metric field]" << std::endl << std::endl;
        m_metric_field->report(fx);
        fx << std::endl << "[Statistics]" << std::endl << std::endl;
        fx << "elapsed time:       " << time(NULL) - start_time << " sec." << std::endl;
        fx << "number of vertices: " << m_stars.size() << std::endl;
        fx << "vertices via picking: " << vertex_with_picking_count << std::endl;
        fx << "vertices non-picking: " << vertex_without_picking_count << std::endl;
        fx << "picking rate:       " << (double)vertex_with_picking_count / (double)(vertex_without_picking_count + vertex_with_picking_count) << std::endl;
        fx << "vertices with smoothing:  " << vertex_with_smoothing_counter << std::endl;
        fx << "vertices w/out smoothing: " << vertex_without_smoothing_counter << std::endl;
        fx.close();

        typename std::ofstream fe("report_times_surface.txt", std::ios_base::app);
        fe << m_stars.size() << "\t" << time(NULL) - start_time << std::endl;
        fe.close();
      }

      template<typename DualType, typename OutputIterator>
      void get_dual_segments(OutputIterator oit)
      {
        for(std::size_t i = 0; i < m_stars.size(); i++)
        {
          Star_handle si = m_stars[i];
          typename Star::Base::Finite_facets_iterator fit 
            = si->finite_facets_begin();
          typename Star::Base::Finite_facets_iterator fend 
            = si->finite_facets_end();
          for(; fit != fend; fit++)
          {
            typename K::Object_3 o = si->dual(*fit);
            typename K::Segment_3 s; 
            typename K::Ray_3 r;
            typename K::Line_3 l;
            if(CGAL::assign(s, o))
              *oit++ = DualType(s.source(), s.target());
            else if(CGAL::assign(r, o))
              *oit++ = DualType(r.source(), r.source() + r.to_vector() 
                * (m_pConstrain->get_bounding_radius() * 2.0 / std::sqrt(r.to_vector() * r.to_vector())));
            else if(CGAL::assign(l, o))
              *oit++ = DualType(Point_3(0.,0.,0.),Point_3(0.,0.,0.));            
          }
        }
      }

      void output_dual()
      {
        typename std::ofstream fx("dual.polylines.txt");
        
        typedef typename std::pair<Point_3, Point_3> PointsEdge;
        std::vector<PointsEdge> edges;
        get_dual_segments<PointsEdge>(std::back_inserter(edges));
        fx << (edges.size()) << std::endl;

        for(std::size_t i = 0; i < edges.size(); i++)
          fx << " " << edges[i].first << " " << edges[i].second << "\n";

        fx.close();
      }


      void output(const char* filename = "points.surface.off",
                  const bool consistent_only = true) 
      {
        std::cout << "Saving " << filename << "...";
        std::string file(filename);
        if(file.substr(file.length()-4).compare(".off") != 0)
          std::cout << "Warning : file name does not end with .off" << std::endl;

        std::ofstream fx(filename);
        output(fx,consistent_only);
      }

      void output(std::ofstream& fx, const bool consistent_only = true)
      {
        std::cout << "Saving " << "...";
        std::map<Index, int> match_indices;//because we won't use all of them
        int off_index = 0; // the corresponding index in .off 
        std::vector<Point_3> points;
        Output_facets output_facets;       

        Star_iterator si = m_stars.begin();
        Star_iterator siend = m_stars.end();
        for (; si != siend; ++si) 
        {
          Star_handle star = *si;
          if(!star->is_surface_star())
            continue;

          points.push_back(star->center_point());
          match_indices[star->index_in_star_set()] = off_index++;

          typename Star::Facet_set_iterator fi = star->begin_restricted_facets();
          typename Star::Facet_set_iterator fiend = star->end_restricted_facets();
          for (; fi != fiend; ++fi) 
          {
            Point_3 cc;
            star->compute_dual_intersection(*fi, cc);

            bool consistent = true;
            if(consistent_only)
            {
              for (int i = 1; i <= 3; i++) 
              {
                Vertex_handle v = fi->first->vertex((fi->second + i) % 4);
                if (v == star->center())
                  continue;
                if (!m_stars[v->info()]->is_surface_star())
                  continue;
                if (!m_stars[v->info()]->has_facet(*fi)) 
                {
                  consistent = false;
                  break;
                }
              }
            }
            if (consistent)
              output_facets.insert(fi->first->vertex((fi->second + 1) % 4)->info(),
                                   fi->first->vertex((fi->second + 2) % 4)->info(),
                                   fi->first->vertex((fi->second + 3) % 4)->info());
          }
        }

        fx << "OFF" << std::endl;
        fx << points.size() << " " << output_facets.size() << " " << 0 << std::endl;
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
        fx.close();
        std::cout << "done.\n";
      }

      // warning : this function empties the priority queue
      void output_refine_queue(int queue_type)
      {
        typename std::set<Facet_ijk> facets;
        Index_set vertices;
        while(true) 
        {
          Refine_facet refine_facet;
          Facet facet;
          if (!m_refine_queue.top(refine_facet, queue_type))
            break;
          m_refine_queue.pop();
          
          if (refine_facet.star->has_facet_ref(refine_facet.vertices, facet)) 
          {
            Facet_ijk f(facet);
            vertices.insert(f.vertex(0));
            vertices.insert(f.vertex(1));
            vertices.insert(f.vertex(2));
            facets.insert(f);
          }
        }

        typename std::ofstream fx("refine_queue.off");
        fx << "OFF" << std::endl;
        fx << vertices.size() << " " << facets.size() << " " << 0 << std::endl;
        
        typename Index_set::iterator iit;
        int i = 0;
        std::map<int, int> correspondances;
        for(iit = vertices.begin(); iit != vertices.end(); iit++, i++)
        {
          int index = *iit;
          fx << m_stars[index]->center_point() << std::endl;
          correspondances[index] = i;
        }
        typename std::set<Facet_ijk>::iterator fit;
        for(fit = facets.begin(); fit != facets.end(); fit++)
          fx << "3 " << correspondances[(*fit).vertex(0)] 
             << " "  << correspondances[(*fit).vertex(1)]
             << " "  << correspondances[(*fit).vertex(2)] << std::endl;
        
        fx.close();
      }

      //void output_metric_field_eigenvalues()
      //{
      //  typename std::ofstream fx("eigenvalues.txt");
      //  m_metric_field->report(fx);
      //  fx << std::endl;

      //  std::size_t i;
      //  for(i = 0; i < m_stars.size(); i++)
      //  {
      //    fx << "[" << m_stars[i]->center_point() << "]\n";
      //    typename CGAL::Anisotropic_mesh_3::Metric_base<K>::Aff_transformation_3
      //      t(m_stars[i]->metric().get_transformation());
      //    for(int j = 0; j < 3; j++)
      //    {
      //      for(int k = 0; k < 3; k++)
      //        fx << t.cartesian(j,k) << "\t";
      //      fx << std::endl;
      //    }
      //    fx << std::endl;
      //  }
      //  fx.close();
      //}

      void output_star_set()
      {
        typename std::ofstream fx("starset.mesh");
        fx << "MeshVersionFormatted 1\n\n";
        fx << "Dimension\n3" << std::endl;

        fx << "Vertices" << std::endl;
        fx << (m_stars.size()) << std::endl;
        for (int i = 0; i < (int)m_stars.size(); i++)
          fx << m_stars[i]->center_point() << " " << i << std::endl;
  
        fx << "Tetrahedra" << std::endl;
        fx << number_of_tets_in_star_set() << std::endl;
        for (int i = 0; i < (int)m_stars.size(); i++)
        {
          Star_handle si = m_stars[i];
          typename Star::Cell_handle_handle cit = si->begin_finite_star_cells();
          typename Star::Cell_handle_handle cend = si->end_finite_star_cells();
          for(; cit != cend; cit++)
          {
            Cell_handle c = *cit;
            //if(si->is_inside(c))
              fx << c->vertex(0)->info() << " " << c->vertex(1)->info() << " "
                 << c->vertex(2)->info() << " " << c->vertex(3)->info() << " " 
                 << i /*color*/ << std::endl;
          }
        }
      }
      
      int number_of_tets_in_star_set() const
      {
        int count = 0;
         for (int i = 0; i < (int)m_stars.size(); i++)
        {
          Star_handle si = m_stars[i];
          typename Star::Cell_handle_handle cit = si->begin_finite_star_cells();
          typename Star::Cell_handle_handle cend = si->end_finite_star_cells();
          for(; cit != cend; cit++)
           // if(si->is_inside(*cit))
              count++;
        }
        return count;
      }

public:
      void update_stars_criteria(){
        for(std::size_t i = 0; i < m_stars.size(); i++){
          m_stars[i]->set_criteria(m_criteria);
        }
      }

      void clean_stars() //remove useless vertices
      {
        for(std::size_t i = 0; i < m_stars.size(); i++)
          m_stars[i]->clean();
      }

public:
      void ouput_sq_radius_edge_ratios(const std::string& filename) const
      {
        std::ofstream fout(filename.c_str());
        for(std::size_t i = 0; i < m_stars.size(); i++)
        {
          Star_handle s = m_stars[i];
          if(!s->is_surface_star() || !s->is_topological_disk())
            continue;
          typename Star::Facet_set_iterator fi = s->begin_restricted_facets();
          typename Star::Facet_set_iterator fend = s->end_restricted_facets();
          for(; fi != fend; fi++)
          {
            FT rer = s->compute_squared_radius_edge_ratio(*fi);
            fout << rer << std::endl;
          }
        }
        fout.close();
      }
      
public:
      void refine_all(std::ofstream& fx, 
                      const double& starttime,
                      const int max_count = INT_MAX,
                      const bool pick_valid_causes_stop = false,
                      const bool pick_valid_use_cube_probing = false)
      {
        //if you modify this, do not forget to also modify the demo
        CGAL::Timer t;
        const int pick_valid_max_failures = 100;
        int pick_valid_failed_n = 0;
        int pick_valid_succeeded_n = 0;

        t.start();
        fill_refinement_queue();         
        int nbv = m_stars.size();
        while(nbv < max_count) 
        {
          if(nbv == 4) //dimension 3 reached
            update_bboxes();         
          if(nbv % 100 == 0)
          { 
            t.stop();
            fx << nbv << "\t" << (starttime + t.time()) << std::endl;
            t.start();
            clean_stars();//remove useless vertices
          }
          if(!refine(pick_valid_succeeded_n, pick_valid_failed_n,
                     pick_valid_causes_stop, pick_valid_max_failures,
                     pick_valid_use_cube_probing))
            break;
          nbv = m_stars.size();
        }
        t.stop();
        fx << nbv << "\t" << (starttime + t.time())  << std::endl;
      }

public:
      void refine_all(const int max_count = INT_MAX,
                      const bool pick_valid_causes_stop = false,
                      const bool pick_valid_use_cube_probing = false)
      {
        //if you modify this, do not forget to also modify the demo
#ifdef ANISO_VERBOSE
        std::cout << "\nRefine all...";
        std::clock_t start_time = clock();
#endif        
        const int pick_valid_max_failures = 100;
        int pick_valid_failed_n = 0;
        int pick_valid_succeeded_n = 0;

        fill_refinement_queue();

        vertex_with_picking_count = 0;
        vertex_without_picking_count = (int)m_stars.size();
#ifdef ANISO_VERBOSE
        std::cout << "There are " << count_restricted_facets() << " restricted facets.\n";
#endif
        std::size_t nbv = m_stars.size();
        while(nbv < max_count) 
        {
          if(nbv == 4) //dimension 3 reached
            update_bboxes();
          
          if(nbv % 100 == 0)
          { 
            clean_stars();//remove useless vertices
#ifdef ANISO_VERBOSE
            std::cerr << " " << nbv << " vertices, ";
            std::cerr << duration(start_time) << " sec.,\t";
            m_refine_queue.print();
#endif
#ifdef ANISO_DUMP_INTERMEDIATE_MESHES
            std::ostringstream oss;
            oss << "out_" << nbv << ".off";
            output(oss.str().c_str(), false/*consistent_only*/);
#endif
          }
          
          if(!refine(pick_valid_succeeded_n, pick_valid_failed_n,
                     pick_valid_causes_stop, pick_valid_max_failures,
                     pick_valid_use_cube_probing))
          {
            clean_stars();
            //debug_show_distortions();
            break;
          }
          nbv = m_stars.size();
        }

#ifdef ANISO_VERBOSE
        double time = duration(start_time);
        std::cout << "\nRefinement done (" << nbv << " vertices in " << time << " seconds)\n";
        if(pick_valid_failed())
          std::cout << "Pick valid failed and stopped mesher!" << std::endl;
        
        if(is_consistent(true/*verbose*/))
          std::cout << "Triangulation is consistent.\n";
        else
          std::cout << "Triangulation is not consistent.\n";

        std::cout << "Vertices via picking: " << vertex_with_picking_count << std::endl;
        std::cout << "Vertex non-picking:   " << vertex_without_picking_count << std::endl;
        std::cout << "picking rate:         " << (double)vertex_with_picking_count / 
          (double)(vertex_without_picking_count + vertex_with_picking_count) << std::endl;
        std::cout << "Approximation error : " << compute_approximation_error() << std::endl;
        std::cout << "Vertices with smoothing:  " << vertex_with_smoothing_counter << std::endl;
        std::cout << "Vertices w/out smoothing: " << vertex_without_smoothing_counter << std::endl;
        
        report();
        histogram_vertices_per_star<Self>(*this);
        output();
#endif
#ifdef USE_ANISO_TIMERS
        report_timers();
#endif
      }

      bool pick_valid_failed() const
      {
        return !m_pick_valid_cache.empty();
      }

      void debug_show_distortions() const
      {
        for(std::size_t i = 0; i < m_stars.size(); ++i)
        {
          Star_handle s = m_stars[i];
          if(!s->is_surface_star())
            continue;
          std::cout << "  " << i << " : ";
          typename Star::Facet_set_iterator fit = s->begin_restricted_facets();
          typename Star::Facet_set_iterator fend = s->end_restricted_facets();
          for(; fit != fend; ++fit)
          {
            std::cout << "(";
            Facet f = *fit;
            if(is_surface_facet(f))
            {
              for (int i = 0; i < 3; i++) 
              {
                int index_1 = (f.second + i + 1) % 4;
                int index_2 = (f.second + (i + 1) % 3 + 1) % 4;
                FT distortion =
                  m_stars[f.first->vertex(index_1)->info()]->metric().compute_distortion(
                  m_stars[f.first->vertex(index_2)->info()]->metric());

                FT costheta = std::abs(
                  m_stars[f.first->vertex(index_1)->info()]->metric().get_vmin()
                  * m_stars[f.first->vertex(index_2)->info()]->metric().get_vmin());

                std::cout << (distortion - 1./costheta) << "\t";
              }
            }
            std::cout << ")";
          }
          std::cout << std::endl;
        }
      }


public:
    double duration(const time_t& start) const
    {
      return ((clock() - start + 0.) / ((double)CLOCKS_PER_SEC));
    }

private:
     void cout_info()
     {
       //total nbv
       int count = 0;
       for(std::size_t i = 0; i < m_stars.size(); i++)
         count += m_stars[i]->number_of_vertices();
       std::cout << "\t\tTotal nb of vertices : " << count << std::endl;

       //stars sizes
       int step = 10;//count / (20 * m_stars.size());
       std::map<int, int> sizes;
       for(std::size_t i = 0; i < m_stars.size(); i++)
       {
         int box = (m_stars[i]->number_of_vertices()/step)*step;
         if(sizes.find(box) == sizes.end()) sizes[box] = 1;
         else                               sizes[box]++;
       }
       std::cout << "\t\tHistogram of nbv/star : " << std::endl;
       for(std::map<int, int>::iterator it = sizes.begin(); it != sizes.end(); it++)
         std::cout << "\t\t\t" << it->first << "\t" << it->second << std::endl;
     }

public:
      void gl_draw(const typename K::Plane_3& plane,
                   const bool draw_edges,
                   const int star_id = -1/*only this one*/) const 
      {
        if(star_id < 0) // draw them all
          for(std::size_t i = 0; i < m_stars.size(); i++)
            m_stars[i]->gl_draw(plane);//, draw_edges);
        else
          m_stars[star_id]->gl_draw(plane);//, colors, draw_edges);
      }

      void gl_draw_cell(const typename K::Plane_3& plane,
                   const int star_id = -1/*only this one*/) const
      {
        Point_3 debug_a(5, 5, 5);
        Vector_3 e1(1, 0, 0);
        Vector_3 e2(0, 1, 0);
        Vector_3 e3(0, 0, 1);
        Vector_3 e4(0.5*std::sqrt(2.), 0.5*std::sqrt(2.), 0);
        Vector_3 e5(-0.5*std::sqrt(2.), 0.5*std::sqrt(2.), 0);
        Vector_3 v1, v2, vn;

        Metric debug_ma(e3, e1, e2, 1.0, 0.4, 1.2, 0);
        Metric debug_mb(e3, e4, e5, 0.4, 0.1, 4, 0);
        Metric debug_mp = m_metric_field->intersection(debug_ma, debug_mb);

        debug_ma = metricA;
        debug_mb = metricB;
        debug_mp = metricP;

        FT ma_a = 1./debug_ma.get_max_eigenvalue();
        FT ma_b = 1./debug_ma.get_min_eigenvalue();
        FT ma_c = 1./debug_ma.get_third_eigenvalue();

        FT mb_a = 1./debug_mb.get_max_eigenvalue();
        FT mb_b = 1./debug_mb.get_min_eigenvalue();
        FT mb_c = 1./debug_mb.get_third_eigenvalue();

        FT mp_a = 1./debug_mp.get_max_eigenvalue();
        FT mp_b = 1./debug_mp.get_min_eigenvalue();
        FT mp_c = 1./debug_mp.get_third_eigenvalue();

        //A
        debug_ma.get_max_eigenvector(v1);
        debug_ma.get_min_eigenvector(v2);
        debug_ma.get_third_eigenvector(vn);

        ::GLdouble rot_mat[16];
        rot_mat[0] = v1.x(); rot_mat[4] = v2.x(); rot_mat[8] = vn.x();  rot_mat[12] = 5;
        rot_mat[1] = v1.y(); rot_mat[5] = v2.y(); rot_mat[9] = vn.y();  rot_mat[13] = 5;
        rot_mat[2] = v1.z(); rot_mat[6] = v2.z(); rot_mat[10] = vn.z(); rot_mat[14] = 5;
        rot_mat[3] = 0.; rot_mat[7] = 0.; rot_mat[11] = 0.; rot_mat[15] = 1.;

        ::glMatrixMode (GL_MODELVIEW);
        ::glPushMatrix();
        ::glMultMatrixd(rot_mat);
        gl_draw_ellipsoid<K>(CGAL::ORIGIN, 20, 20, ma_a, ma_b, ma_c, 240, 20, 20);
        ::glPopMatrix();

        //B
        debug_mb.get_max_eigenvector(v1);
        debug_mb.get_min_eigenvector(v2);
        debug_mb.get_third_eigenvector(vn);

        rot_mat[0] = v1.x(); rot_mat[4] = v2.x(); rot_mat[8] = vn.x();  rot_mat[12] = 5;
        rot_mat[1] = v1.y(); rot_mat[5] = v2.y(); rot_mat[9] = vn.y();  rot_mat[13] = 5;
        rot_mat[2] = v1.z(); rot_mat[6] = v2.z(); rot_mat[10] = vn.z(); rot_mat[14] = 5;
        rot_mat[3] = 0.; rot_mat[7] = 0.; rot_mat[11] = 0.; rot_mat[15] = 1.;

        ::glMatrixMode (GL_MODELVIEW);
        ::glPushMatrix();
        ::glMultMatrixd(rot_mat);
        gl_draw_ellipsoid<K>(CGAL::ORIGIN, 20, 20, mb_a, mb_b, mb_c, 20, 240, 20);
        ::glPopMatrix();

        //P
        debug_mp.get_max_eigenvector(v1);
        debug_mp.get_min_eigenvector(v2);
        debug_mp.get_third_eigenvector(vn);

        rot_mat[0] = v1.x(); rot_mat[4] = v2.x(); rot_mat[8] = vn.x();  rot_mat[12] = 5;
        rot_mat[1] = v1.y(); rot_mat[5] = v2.y(); rot_mat[9] = vn.y();  rot_mat[13] = 5;
        rot_mat[2] = v1.z(); rot_mat[6] = v2.z(); rot_mat[10] = vn.z(); rot_mat[14] = 5;
        rot_mat[3] = 0.; rot_mat[7] = 0.; rot_mat[11] = 0.; rot_mat[15] = 1.;

        ::glMatrixMode (GL_MODELVIEW);
        ::glPushMatrix();
        ::glMultMatrixd(rot_mat);
        gl_draw_ellipsoid<K>(CGAL::ORIGIN, 20, 20, mp_a, mp_b, mp_c, 20, 20, 240);
        ::glPopMatrix();

        ::glPointSize(5.);
        ::glBegin(GL_POINTS);
        ::glColor3f(0.87f, 0.14f, 0.14f);
        ::glVertex3d(debug_a.x(), debug_a.y(), debug_a.z());
        ::glEnd();


        if(star_id < 0) // draw them all
          for(std::size_t i = 0; i < m_stars.size(); i++)
            m_stars[i]->gl_draw_cell(plane);
        else
          m_stars[star_id]->gl_draw_cell(plane);
      }

      void gl_draw_picked_points(const typename K::Plane_3& plane,
                                 const int point_id = -1/*only this one*/) const
      {
        //facet trying to be refined
        if(!m_pick_valid_facet.empty())
        {
          gl_draw_triangle<K>(m_pick_valid_facet[0],
                              m_pick_valid_facet[1],
                              m_pick_valid_facet[2],
                              EDGES_AND_FACES, 84, 204, 204);
        }
        GLboolean light = (::glIsEnabled(GL_LIGHTING));
        if(light)
          ::glDisable(GL_LIGHTING);


        if(!red_points.empty() || !orange_points.empty() || !yellow_points.empty() || !green_points.empty())
        {//cube probing
          ::glPointSize(5.);
          ::glBegin(GL_POINTS);

          ::glColor3f(0.87f, 0.14f, 0.14f);
          for(std::size_t i=0; i<red_points.size();++i)
          {
            Point_3 red_point = red_points[i];
            ::glVertex3d(red_point.x(), red_point.y(), red_point.z());
          }

          ::glColor3f(0.87f, 0.55f, 0.10f);
          for(std::size_t i=0; i<orange_points.size();++i)
          {
            Point_3 orange_point = orange_points[i];
            ::glVertex3d(orange_point.x(), orange_point.y(), orange_point.z());
          }

          ::glColor3f(0.95f, 0.95f, 0.10f);
          for(std::size_t i=0; i<yellow_points.size();++i)
          {
            Point_3 yellow_point = yellow_points[i];
            ::glVertex3d(yellow_point.x(), yellow_point.y(), yellow_point.z());
          }

          ::glColor3f(0.14f, 0.87f, 0.14f);
          for(std::size_t i=0; i<green_points.size();++i)
          {
            Point_3 green_point = green_points[i];
            ::glVertex3d(green_point.x(), green_point.y(), green_point.z());
          }
          ::glEnd();

          if(light)
            ::glEnable(GL_LIGHTING);
        }

        float point_size = 5.;
        if(!red_points.empty() || !orange_points.empty() || !yellow_points.empty() || !green_points.empty())
         point_size = 0.01f;

        //all pickvalid points that were tried
        if(point_id < 0)
        {
          for(std::size_t i = 0; i<m_pick_valid_cache.size(); ++i)
          {
            Point_3 pickvalid_point = m_pick_valid_cache[i];

            //keeping all these points independant of the cut plane
            // if(!is_above_plane(plane, pickvalid_point))
            //   continue;

            //if(i != m_pick_valid_cache.size()-1) //point size depends on the number of problematic facets
            //  point_size += ((pickvalid_problematic_facets[pickvalid_point]).size())/2;
            ::glPointSize(point_size);

            ::glBegin(GL_POINTS);
            if( i == m_pick_valid_cache.size()-1 )
              ::glColor3f(0.87f, 0.14f, 0.14f); //circumcenter in red
            else
              ::glColor3f(0.14f, 0.87f, 0.14f); //others in green

            ::glVertex3d(pickvalid_point.x(), pickvalid_point.y(), pickvalid_point.z());
            ::glEnd();
          }

          if(light)
            ::glEnable(GL_LIGHTING);
        }
        else //only one
        {
          assert(point_id >= 0 && point_id < m_pick_valid_cache.size());
          Point_3 picked_point = m_pick_valid_cache[point_id];
          std::vector<int> picked_point_couples = pickvalid_problematic_facets[picked_point];

          //the pickvalid point
          ::glPointSize(5.);
          ::glBegin(GL_POINTS);
          ::glColor3f(0.14f, 0.87f, 0.14f);
          ::glVertex3d(picked_point.x(), picked_point.y(), picked_point.z());
          ::glEnd();

          if(light)
            ::glEnable(GL_LIGHTING);

          //all its problematic facets
          for(std::size_t i=0; i<picked_point_couples.size();)
          {
            gl_draw_triangle<K>(picked_point,
                                m_stars[picked_point_couples[i]]->center_point(),
                                m_stars[picked_point_couples[i+1]]->center_point(),
                                EDGES_AND_FACES, 166, 247, 170);
            i+=2;
          }
        }
      }

      void gl_draw_initial_points(const typename K::Plane_3& plane) const
      {
        typename Constrain_surface::Pointset::const_iterator pi = initial_points.begin();
        typename Constrain_surface::Pointset::const_iterator pend = initial_points.end();
        GLboolean light = (::glIsEnabled(GL_LIGHTING));

        ::glPointSize(10.);
        ::glColor3f(0.706f, 0.345f, 0.878f);
        if(light)
          ::glDisable(GL_LIGHTING);
        ::glBegin(GL_POINTS);
        for (; pi != pend; pi++)
        {
          // if(!is_above_plane(plane, *it))
          //   continue;

          ::glVertex3d((*pi).x(), (*pi).y(), (*pi).z());
        }
        ::glEnd();
        if(light)
          ::glEnable(GL_LIGHTING);
      }

      void gl_draw_poles(const typename K::Plane_3& plane) const
      {
        GLboolean light = (::glIsEnabled(GL_LIGHTING));
        ::glPointSize(5.);
        ::glColor3f(0.98f, 0.757f, 0.137f);
        if(light)
          ::glDisable(GL_LIGHTING);

        ::glBegin(GL_POINTS);
        typename std::set<Point_3>::const_iterator it;
        for(it = poles.begin(); it != poles.end(); ++it)
        {
          // if(!is_above_plane(plane, *it))
          //   continue;
          ::glVertex3d((*it).x(), (*it).y(), (*it).z());
        }
        ::glEnd();
        if(light)
          ::glEnable(GL_LIGHTING);
      }

      void gl_draw_metric(const typename K::Plane_3& plane,
                          double bbox_min, double eps,
                          const int star_id = -1/*only this one*/) const 
      {
        double coeff = bbox_min/20.;
        double glob_min = (std::max)(eps, m_pConstrain->global_min_curvature());
        coeff *= glob_min;

        if(star_id < 0) // draw them all
          for(std::size_t i = 0; i < m_stars.size(); i++)
            m_stars[i]->gl_draw_metric(plane, coeff);
        else
          m_stars[star_id]->gl_draw_metric(plane, coeff);
      }

      void gl_draw_dual(const typename K::Plane_3& plane,
                        const int star_id = -1) const
      {
        if(star_id < 0) // draw them all
          for(std::size_t i = 0; i < m_stars.size(); i++)
            m_stars[i]->gl_draw_dual(plane);
        else
          m_stars[star_id]->gl_draw_dual(plane);
      }

      void gl_draw_surface_delaunay_balls(const typename K::Plane_3& plane,
                                          const int star_id = -1) const
      {
        if(star_id < 0)
          for(std::size_t i = 0; i < m_stars.size(); i++)
            m_stars[i]->gl_draw_surface_delaunay_balls(plane);
        else
          m_stars[star_id]->gl_draw_surface_delaunay_balls(plane);
      }

      bool is_above_plane(const typename K::Plane_3& plane,
                          const typename K::Point_3& pa,
                          const typename K::Point_3& pb,
                          const typename K::Point_3& pc) const
      {
        typedef typename K::Oriented_side Side;
        using CGAL::ON_ORIENTED_BOUNDARY;
        using CGAL::ON_NEGATIVE_SIDE;
        const Side sa = plane.oriented_side(pa);
        const Side sb = plane.oriented_side(pb);
        const Side sc = plane.oriented_side(pc);
        return (sa == ON_NEGATIVE_SIDE && sb == ON_NEGATIVE_SIDE && sc == ON_NEGATIVE_SIDE);
      }

      void gl_draw_distortion(const typename K::Plane_3& plane,
                              const int star_id = -1) const
      {
        GLboolean was = (::glIsEnabled(GL_LIGHTING));
        if(was)
          ::glDisable(GL_LIGHTING);

        bool draw_all = (star_id < 0);
        std::size_t start = draw_all ? 0 : star_id;
        std::size_t N = draw_all ? m_stars.size() : (star_id + 1);

        std::set<Facet_ijk> done;
        for(std::size_t i = start; i < N; i++)
        {
          Star_handle star = m_stars[i];
          typename Star::Facet_set_iterator fit = star->begin_restricted_facets();
          typename Star::Facet_set_iterator fitend = star->end_restricted_facets();
          for(; fit != fitend; fit++)
          {
            Facet f = *fit;
            if(done.find(Facet_ijk(f)) != done.end())
              continue;
            done.insert(Facet_ijk(f));

            const Point_3& pa = transform_from_star_point(f.first->vertex((f.second+1)%4)->point(), star);
            const Point_3& pb = transform_from_star_point(f.first->vertex((f.second+2)%4)->point(), star);
            const Point_3& pc = transform_from_star_point(f.first->vertex((f.second+3)%4)->point(), star);
            if(!is_above_plane(plane, pa, pb, pc))
              continue;

            FT max_distortion = 0.;
            for (int i = 0; i < 3; i++) 
            {
              int index_1 = (f.second + i + 1) % 4;
              int index_2 = (f.second + (i + 1) % 3 + 1) % 4;
              FT distortion = m_stars[f.first->vertex(index_1)->info()]->metric().compute_distortion(
                 m_stars[f.first->vertex(index_2)->info()]->metric());
              max_distortion = (std::max)(distortion, max_distortion);
            }
            max_distortion = 120.*(max_distortion - 1.);
            float rgf = static_cast<float>((std::max)(0., 255. - max_distortion));
            gl_draw_triangle<K>(pa, pb, pc, FACES_ONLY, rgf, rgf, 255);
          }
        }
        if(was)
          ::glEnable(GL_LIGHTING);
      }

      void gl_draw_inconsistent_facets(const typename K::Plane_3& plane,
                                       const int star_id = -1) const
      {
        GLboolean was = (::glIsEnabled(GL_LIGHTING));
        if(!was)
          ::glEnable(GL_LIGHTING);

        ::glPolygonOffset(1.f, 0.1f);
        ::glEnable(GL_POLYGON_OFFSET_FILL);
        ::glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        bool draw_all = (star_id < 0);
        std::size_t start = draw_all ? 0 : star_id;
        std::size_t N = draw_all ? m_stars.size() : (star_id + 1);

        for(std::size_t i = start; i < N; i++)
        {
          Star_handle star = m_stars[i];
          typename Star::Facet_set_iterator fit = star->begin_restricted_facets();
          typename Star::Facet_set_iterator fitend = star->end_restricted_facets();
          for(; fit != fitend; fit++)
          {
            Facet f = *fit;
            if(is_consistent(f))
              continue;

            const Point_3& pa = transform_from_star_point(f.first->vertex((f.second+1)%4)->point(), star);
            const Point_3& pb = transform_from_star_point(f.first->vertex((f.second+2)%4)->point(), star);
            const Point_3& pc = transform_from_star_point(f.first->vertex((f.second+3)%4)->point(), star);

            if(is_above_plane(plane, pa, pb, pc))
              gl_draw_triangle<K>(pa,pb,pc,EDGES_AND_FACES, 227,27,27);
          }
        }
        ::glDisable(GL_POLYGON_OFFSET_FILL);
        if(!was)
          ::glDisable(GL_LIGHTING);
      }   
            
#ifdef USE_ANISO_TIMERS
      void reset_timers() const
      {
        m_pick_valid_timer = 0.;
        m_insert_timer = 0.;
        m_fill_queue_timer  = 0.;
        m_compute_steiner_timer = 0.;
        m_create_star_timer = 0.;
      }
      void report_timers() const
      {
        std::cout << "Timers : " << std::endl;
        std::cout << "\tfill_queue      : " << m_fill_queue_timer << " sec." << std::endl;
        std::cout << "\tinsert          : " << m_insert_timer << " sec." << std::endl;
        std::cout << "\tcompute_steiner : " << m_compute_steiner_timer << " sec." << std::endl;
        std::cout << "\tpick_valid      : " << m_pick_valid_timer << " sec." << std::endl;
        std::cout << "\tdual_intersection : " << Star::m_compute_dual_intersection_timer << " sec." << std::endl;
        m_aabb_tree.report_timers();
      }
#endif

    public:
      Surface_star_set_3(const Criteria* criteria_,
        const Metric_field* metric_field_,
        const Constrain_surface* const pconstrain_,
        const int nb_initial_points = 10,
        const RefinementCondition& rc_ = RefinementCondition())
        :
        m_pConstrain(pconstrain_),
        m_metric_field(metric_field_), 
        m_criteria(criteria_), 
        m_stars(), 
        m_refine_queue(),
        m_ch_triangulation(),
        m_refinement_condition(rc_),
        m_aabb_tree(100/*insertion buffer size*/),
        m_kd_tree(m_stars),
        vertex_with_smoothing_counter(0),
        vertex_without_smoothing_counter(0)
      {
        initialize_stars(nb_initial_points); // initialize poles + points on surface
        m_ch_triangulation.infinite_vertex()->info() = -10;

#ifndef NO_USE_AABB_TREE_OF_BBOXES 
        build_aabb_tree();  // update bboxes and rebuild
#endif
#ifdef USE_ANISO_TIMERS
        reset_timers();
#endif
      }

      void clear() 
      {
        m_stars.clear();
        m_aabb_tree.clear();
        m_kd_tree.clear();
        m_ch_triangulation.clear();
        m_refine_queue.clear();
      }
      
      ~Surface_star_set_3() 
      {
        delete m_metric_field;
        delete m_criteria;
        Star_iterator si = m_stars.begin();
        Star_iterator siend = m_stars.end();
        for (; si != siend; si++)
          delete (*si);
      }
    private:
      Surface_star_set_3(const Surface_star_set_3&) {}
      Surface_star_set_3& operator=(const Surface_star_set_3&) {}

    };
  }//namespace Anisotropic_mesh_3
}//namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_SURFACE_STAR_SET_3_H

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

#ifndef CGAL_ANISOTROPIC_MESH_3_STRETCHED_DELAUNAY_3
#define CGAL_ANISOTROPIC_MESH_3_STRETCHED_DELAUNAY_3

#include <iostream>
#include <fstream>
#include <utility>
#include <time.h>
#include <CGAL/Profile_counter.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Random.h>

#include <CGAL/bbox.h>
#include <CGAL/Metric.h>
#include <CGAL/Criteria.h>
#include <CGAL/Delaunay_traits_3.h>
#include <CGAL/Constrain_surface_3.h>
#include <CGAL/Triangulation_cell_base_with_domain_info_3.h>

#include <CGAL/helpers/combinatorics_helper.h>
#include <CGAL/gl_draw/drawing_helper.h>

namespace CGAL{
  namespace Anisotropic_mesh_3{


    template<typename K, typename KExact = K>
    class Stretched_Delaunay_3 : public CGAL::Delaunay_triangulation_3<
      Delaunay_traits_3<K, KExact>,
      CGAL::Triangulation_data_structure_3<
      CGAL::Triangulation_vertex_base_with_info_3<int/*std::size_t*/, Delaunay_traits_3<K, KExact> >,
      CGAL::Triangulation_cell_base_with_domain_info_3<Delaunay_traits_3<K, KExact>, Constrain_surface_3<K> > > >
    {
    public:
      typedef Delaunay_traits_3<K, KExact> Traits;
      typedef CGAL::Delaunay_triangulation_3<
        Traits,
        CGAL::Triangulation_data_structure_3<
        CGAL::Triangulation_vertex_base_with_info_3<int/*std::size_t*/, Traits>,
        CGAL::Triangulation_cell_base_with_domain_info_3<Traits, Constrain_surface_3<K> >
        >
      > Base;
      typedef Stretched_Delaunay_3<K, KExact> Self;
      typedef CGAL::Triangulation_data_structure_3<
        CGAL::Triangulation_vertex_base_with_info_3<int/*std::size_t*/, Traits>,
        CGAL::Triangulation_cell_base_with_domain_info_3<Traits, Constrain_surface_3<K> >
      > DS;

      typedef int Index;

      typedef Constrain_surface_3<K>          Constrain_surface;
      typedef Metric_base<K, KExact>          Metric;
      typedef typename K::FT                  FT;
      typedef typename K::Point_3             Point_3;
      typedef typename K::Point_3             TPoint_3; //Transformed_point_3
      typedef typename K::Sphere_3            Sphere;
      typedef typename K::Vector_3            Vector_3;
      typedef typename K::Line_3              Line_3;
      typedef typename K::Ray_3               Ray_3;
      typedef typename K::Plane_3             Plane_3;
      typedef typename Base::Segment            Segment;
      typedef typename Base::Triangle           Triangle;
      typedef typename Base::Tetrahedron        Tetrahedron;
      typedef typename Base::Vertex             Vertex;
      typedef typename Base::Cell               Cell;
      typedef typename Base::Edge               Edge;
      typedef typename Base::Facet              Facet;
      typedef typename DS::Vertex_handle      Vertex_handle;
      typedef typename DS::Cell_handle        Cell_handle;
      typedef ::CGAL::Anisotropic_mesh_3::Stretched_criteria<K, KExact> Stretched_criteria;
      typedef Criteria_base<K>          Criteria;
      typedef CGAL::Bbox<K>             Bbox;

      typedef std::vector<Facet>                        Facet_vector;
      typedef std::set<Facet>                           Facet_set;
      typedef std::vector<Cell_handle>                  Cell_handle_vector;
      typedef std::vector<Vertex_handle>                Vertex_handle_vector;
      typedef std::vector<Point_3>                      Point_vector;
      typedef typename Facet_set::iterator              Facet_set_iterator;
      typedef typename Facet_vector::iterator           Facet_handle;
      typedef typename Cell_handle_vector::iterator     Cell_handle_handle;
      typedef typename Vertex_handle_vector::iterator   Vertex_handle_handle;

    private:
      typedef typename KExact::Point_3                      Exact_Point_3;
      typedef typename KExact::Point_3                      Exact_TPoint_3;
      typedef typename KExact::Vector_3                     Exact_Vector_3;
      typedef CGAL::Cartesian_converter<K, KExact> To_exact;
      typedef CGAL::Cartesian_converter<KExact, K> Back_from_exact;
      To_exact to_exact;
      Back_from_exact back_from_exact;

    private:
      Point_3 m_center_point; //before transformation by m_metric.transform
      Vertex_handle m_center; //after transformation by m_metric.transform
      Traits *m_traits;
      Metric m_metric;
      const Constrain_surface* m_pConstrain;
      Stretched_criteria* m_criteria;
      bool m_is_surface_star; // allows to handle stars with center on medial axis

    private:
      mutable bool m_is_topological_disk;
      mutable bool m_is_valid_topo_disk;
      mutable Bbox m_bbox;
      mutable bool m_is_valid_bbox;

    public:
      int ellipsoid_color;

#ifdef USE_ANISO_TIMERS
    public:
      static double m_compute_dual_intersection_timer;
#endif

    private:
      static const Index index_of_infinite_vertex = -10;
    public:
      static Index infinite_vertex_index() { return index_of_infinite_vertex; }
      bool is_infinite() const
      {
        return (infinite_vertex_index() == index_in_star_set());
      }

    public:
      inline bool is_inside(const Cell_handle &cell) const
      {
        if(this->is_infinite(cell))
          return false;
#ifdef ANISO_DEBUG //added to see coordinates in debugger
        const Point_3& p0 = cell->vertex(0)->point();
        const Point_3& p1 = cell->vertex(1)->point();
        const Point_3& p2 = cell->vertex(2)->point();
        const Point_3& p3 = cell->vertex(3)->point();
#endif
        return cell->template is_inside<Self>(*this, *constrain(), *traits());
      }
#ifdef ANISO_USE_INSIDE_EXACT
      inline bool is_inside_exact(const Cell_handle &cell) const
      {
        if(this->is_infinite(cell))
          return false;
        return cell->template is_inside_exact<Self>(*this, *constrain(), *traits());
      }
#endif
      inline bool is_infinite_vertex(const Vertex_handle &v) const
      {
        return Base::is_infinite(v);
      }
      inline bool is_infinite(const Cell_handle& cell) const
      {
        return Base::is_infinite(cell);
      }

    public:
      Vertex_handle center() const      { return m_center; }
      Metric metric() const             { return m_metric; }
      Traits* traits() const            { return m_traits; }
      const Constrain_surface* constrain() const { return m_pConstrain; }
      const Stretched_criteria* criteria() const { return m_criteria; }
      void set_criteria(const Criteria* criteria_)
      {
        delete m_criteria;
        m_criteria = new Stretched_criteria(*m_traits, criteria_);
      }

    public:
      Point_3& center_point()             { return m_center_point; }
      const Point_3& center_point() const { return m_center_point; }

      Index& index_in_star_set()             { return m_center->info(); }
      const Index& index_in_star_set() const { return m_center->info(); }

      bool& is_surface_star()             { return m_is_surface_star; }
      const bool& is_surface_star() const { return m_is_surface_star; }

    public:
      const Bbox& bbox(const bool verbose = false) const
      {
        update_bbox(verbose);
        return m_bbox;
      }

    protected: // star configuration caching
      mutable bool is_cache_dirty;
      mutable Facet_set restricted_facets_cache; //restricted and incident to m_center
      mutable Cell_handle_vector neighboring_cells_cache;
      mutable Cell_handle_vector neighboring_finite_cells_cache;
      mutable Vertex_handle_vector neighboring_vertices_cache;

//      FT squared_bounding_radius;

      inline void update_star_caches() const
      {
        if (!is_cache_dirty)
          return;

        CGAL_PROFILER("[update_star_caches]");
        neighboring_cells_cache.clear();
        neighboring_finite_cells_cache.clear();
        restricted_facets_cache.clear();
        neighboring_vertices_cache.clear();
//        squared_bounding_radius = DBL_MAX;

        if(Base::dimension() == 2)
        {
          typename Base::Finite_facets_iterator fit = this->finite_facets_begin();
          typename Base::Finite_facets_iterator fend = this->finite_facets_end();
          for(; fit != fend; ++fit)
            if(is_in_star(*fit) && is_restricted(*fit))
              restricted_facets_cache.insert(this->make_canonical(*fit));
        }
        else if(Base::dimension() > 2)
        {
          // update neighboring vertices
          std::back_insert_iterator<Vertex_handle_vector> neighboring_vertices_insertor(neighboring_vertices_cache);
          Base::adjacent_vertices(m_center, neighboring_vertices_insertor);

          // update neighboring cells
          std::back_insert_iterator<Cell_handle_vector> neighboring_cells_insertor(neighboring_cells_cache);
          Base::incident_cells(m_center, neighboring_cells_insertor);

          std::back_insert_iterator<Cell_handle_vector>
            neighboring_finite_cells_insertor(neighboring_finite_cells_cache);
          Base::finite_incident_cells(m_center, neighboring_finite_cells_insertor);

          // update boundary facets
          if(is_surface_star())
          {
            Cell_handle_handle ci = neighboring_finite_cells_cache.begin();
            Cell_handle_handle cend = neighboring_finite_cells_cache.end();
            for (; ci != cend; ci++)
            {
              int center_index = (*ci)->index(m_center);
              for(int i = 0; i < 4; i++)
              {
                if(i == center_index) continue;
                Facet f = this->make_canonical(Facet(*ci, i));
                Point_3 p;
                if(is_restricted(f, p, true)) //updates surface delaunay center, if needed
                  restricted_facets_cache.insert(f);
              }
            }
          }
          // update bounding radius
          //squared_bounding_radius = -1.0;
          //ci = neighboring_cells_cache.begin();
          //cend = neighboring_cells_cache.end();
          //Traits::Compute_squared_distance_3 csd = m_traits->compute_squared_distance_3_object();
          //for (; ci != cend; ci++)
          //{
          //  if (is_infinite(*ci))
          //  {
          //    squared_bounding_radius = DBL_MAX;
          //    break;
          //  }
          //  Point_3 cc = (*ci)->circumcenter(*m_traits);
          //  FT squared_radius = csd(cc, m_center->point()) * 4.01; // a number bigger than 2^2
          //  if (squared_radius > squared_bounding_radius)
          //    squared_bounding_radius = squared_radius;
          //}
          //if (squared_bounding_radius < 0)
          //  squared_bounding_radius = DBL_MAX;
        }

        is_cache_dirty = false;
      }

      void invalidate_cache() const
      {
        is_cache_dirty = true;
        m_is_valid_bbox = false;
        m_is_valid_topo_disk = false;
      }

public:
      Facet make_canonical(const Facet& f) const
      {
        //returns always the same (c,i) for f and its mirror facet
        if(this->dimension() < 3)
          return f;

        Facet f2 = this->mirror_facet(f);
        Index i_f = f.first->vertex(f.second)->info();
        Index i_f2 = f2.first->vertex(f2.second)->info();
#ifdef ANISO_DEBUG
        if(i_f == -10)
          std::cerr << "Warning : index i_f is " << i_f << std::endl;
#endif
        return (i_f > i_f2) ? f : f2;
        // note : infinite vertex has index -10
        // we choose the one with the greatest index to get the facet on the finite side
      }

public:
      void update_bbox(const bool verbose = false) const
      {
        if(m_is_valid_bbox)
          return;

        if(verbose) std::cout << "Update Bbox...";
        CGAL_PROFILER("[update_bbox]");

        if(this->dimension() < 3)
          m_bbox = m_pConstrain->get_bbox();  //should be found by every request to aabb_tree
        else if(this->is_topological_disk())
          m_bbox = this->surface_bbox();
        else
          m_bbox = this->volume_bbox();
        // Note : {surface Delaunay balls} is included inside {volume Delaunay balls}
        m_is_valid_bbox = true;
        if(verbose) std::cout << "done.\n";
      }


      Bbox surface_bbox() const// compute bbox of incident surface Delaunay balls
      {
        typename Traits::Compute_squared_distance_3 csd
          = m_traits->compute_squared_distance_3_object();

        Bbox bb = m_center->point().bbox();
        Facet_set_iterator fi = begin_restricted_facets();
        Facet_set_iterator fend = end_restricted_facets();
        for(; fi != fend; fi++)
        {
          Point_3 p;
#ifdef ANISO_USE_EXACT
          compute_exact_dual_intersection(*fi, p);
#else
          compute_dual_intersection(*fi, p);
#endif
          TPoint_3 tp = m_metric.transform(p);
          FT squared_radius = csd(tp, m_center->point());
          Sphere s(tp, squared_radius);
          bb = bb + s.bbox();
        }
        return m_metric.inverse_transform(bb);
      }

      Bbox volume_bbox() const//// compute bbox of incident cells circumspheres
      {
        typename Traits::Compute_squared_distance_3 csd
          = m_traits->compute_squared_distance_3_object();

        Bbox bb = m_center->point().bbox();
        Cell_handle_handle ci = begin_finite_star_cells();
        Cell_handle_handle cend = end_finite_star_cells();
        for (; ci != cend; ci++)
        {
          TPoint_3 cc = (*ci)->circumcenter(*m_traits);
          FT squared_radius = csd(cc, m_center->point());

#ifdef ANISO_APPROXIMATE_SPHERE_BBOX
          double radius = sqrt(squared_radius);
          Bbox bbs(cc.x()-radius, cc.y()-radius, cc.z()-radius,
                  cc.x()+radius, cc.y()+radius, cc.z()+radius);
          bb = bb + bbs;
#else
          Sphere s(cc, squared_radius);
          bb = bb + s.bbox();
#endif
        }
        return m_metric.inverse_transform(bb);
      }

    public:
      inline FT compute_volume(const Cell_handle &cell)
      {
        return m_criteria->compute_volume(cell->vertex(0)->point(), cell->vertex(1)->point(),
          cell->vertex(2)->point(), cell->vertex(3)->point());
      }

      inline FT compute_volume(const Facet &facet)
      {
        return m_criteria->compute_volume(
          facet.first->vertex((facet.second + 1) % 4)->point(),
          facet.first->vertex((facet.second + 2) % 4)->point(),
          facet.first->vertex((facet.second + 3) % 4)->point());
      }

      inline FT compute_radius_edge_ratio_overflow(const Cell_handle &cell)
      {
        return m_criteria->radius_edge_ratio_overflow(
          cell->vertex(0)->point(), cell->vertex(1)->point(),
          cell->vertex(2)->point(), cell->vertex(3)->point());
      }

      inline FT compute_radius_edge_ratio_overflow(const Facet &facet)
      {
        return m_criteria->radius_edge_ratio_overflow(
          facet.first->vertex((facet.second + 1) % 4)->point(),
          facet.first->vertex((facet.second + 2) % 4)->point(),
          facet.first->vertex((facet.second + 3) % 4)->point());
      }

      inline FT compute_squared_radius_edge_ratio(const Facet& facet)
      {
        return m_criteria->compute_squared_radius_edge_ratio(
          facet.first->vertex((facet.second + 1) % 4)->point(),
          facet.first->vertex((facet.second + 2) % 4)->point(),
          facet.first->vertex((facet.second + 3) % 4)->point());
      }

      inline FT compute_circumradius_overflow(const Cell_handle &cell)
      {
        return m_criteria->circumradius_overflow(
          cell->vertex(0)->point(), cell->vertex(1)->point(),
          cell->vertex(2)->point(), cell->vertex(3)->point());
      }

      inline FT compute_circumradius_overflow(const Facet &facet)
      {
        // surface Delaunay ball radius
        Point_3 p;
        compute_dual_intersection(facet,p);
        TPoint_3 tp = m_metric.transform(p);
        TPoint_3 tp2 = facet.first->vertex((facet.second + 1) % 4)->point();
        FT sqr = CGAL::squared_distance(tp, tp2);
        return sqr - m_criteria->criteria->squared_circumradius;
      }

      inline FT compute_squared_circumradius(const Facet &facet)
      {
        return m_criteria->compute_squared_circumradius(
          facet.first->vertex((facet.second + 1) % 4)->point(),
          facet.first->vertex((facet.second + 2) % 4)->point(),
          facet.first->vertex((facet.second + 3) % 4)->point());
      }

      inline FT compute_squared_circumradius(const TPoint_3& p1,
                                             const TPoint_3& p2,
                                             const TPoint_3& p3)
      {
        return m_criteria->compute_squared_circumradius(p1, p2, p3);
      }

      inline FT compute_sliverity_overflow(const Cell_handle &cell)
      {
        return m_criteria->sliverity_overflow(
          cell->vertex(0)->point(), cell->vertex(1)->point(),
          cell->vertex(2)->point(), cell->vertex(3)->point());
      }

      inline FT compute_sliverity_overflow(const Point_3& p1,
                                           const Point_3& p2,
                                           const Point_3& p3,
                                           const Point_3& p4) const
      {
        return m_criteria->sliverity_overflow(
          m_metric.transform(p1), m_metric.transform(p2), 
          m_metric.transform(p3), m_metric.transform(p4));
      }
          
      inline bool is_sliver(const TPoint_3& p1,
                            const TPoint_3& p2,
                            const TPoint_3& p3,
                            const TPoint_3& p4) const
      {
        return (compute_sliverity_overflow(p1, p2, p3, p4) > 0.);
      }

    public:
      template<typename OutputIndexIterator>
      void finite_star_vertices(OutputIndexIterator oit)
      {
        std::vector<Vertex_handle> vertices;
        finite_adjacent_vertices(m_center, std::back_inserter(vertices));
        Vertex_handle_handle it = vertices.begin();
        Vertex_handle_handle itend = vertices.end();
        for(; it != itend; it++)
          *oit++ = (*it)->info();
      }

    public:
      // star configuration
      inline Facet_set_iterator begin_restricted_facets() const
      {
        update_star_caches();
        return restricted_facets_cache.begin();
      }
      inline Facet_set_iterator end_restricted_facets() const
      {
        return restricted_facets_cache.end();
      }
      inline Cell_handle_handle begin_star_cells() const
      {
        update_star_caches();
        return neighboring_cells_cache.begin();
      }
      inline Cell_handle_handle end_star_cells() const
      {
        return neighboring_cells_cache.end();
      }
      inline Cell_handle_handle begin_finite_star_cells() const
      {
        update_star_caches();
        return neighboring_finite_cells_cache.begin();
      }
      inline Cell_handle_handle end_finite_star_cells() const
      {
        return neighboring_finite_cells_cache.end();
      }
      inline Vertex_handle_handle begin_neighboring_vertices() const
      {
        update_star_caches();
        return neighboring_vertices_cache.begin();
      }
      inline Vertex_handle_handle end_neighboring_vertices() const
      {
        return neighboring_vertices_cache.end();
      }
      //inline FT get_squared_bounding_radius() {
      //  update_star_caches();
      //  return squared_bounding_radius;
      //}
      inline int get_finite_star_cell_count() const
      {
        update_star_caches();
        return (int)neighboring_finite_cells_cache.size();
      }
      inline bool is_boundary_star() const
      {
        update_star_caches();
        return (restricted_facets_cache.size() > 0);
      }

      unsigned int nb_restricted_facets() const
      {
        unsigned int count = 0;
        Facet_set_iterator fit = begin_restricted_facets();
        Facet_set_iterator fend = end_restricted_facets();
        for(; fit != fend; ++fit)
          count++;
        return count;
      }

      bool is_restricted(const Facet& f) const
      {
        Point_3 c;
        return is_restricted(f,c,false/*do not compute intersection*/);
      }
      bool is_restricted(const Facet& f,
                         Point_3& surface_delaunay_ball_center) const
      {
        return is_restricted(f, surface_delaunay_ball_center, true);
      }
      bool is_restricted(const Facet& f,
                         Point_3& surface_delaunay_ball_center,
                         const bool compute_intersection) const
      {
        if(Base::is_infinite(f))
          return false;
        else if(Base::dimension() == 2)
          return is_restricted_2_in_3(f, surface_delaunay_ball_center, true/*=exact*/);
        else if(Base::dimension() == 3)
        {
          Cell_handle c1 = f.first;
          Cell_handle c2 = f.first->neighbor(f.second);
#ifdef ANISO_USE_INSIDE_EXACT
          bool inside1 = is_inside_exact(c1);
          bool inside2 = is_inside_exact(c2);
#else
          bool inside1 = is_inside(c1);
          bool inside2 = is_inside(c2);
#endif
          if((inside1 && !inside2) || (!inside1 && inside2))
          {
            if(!compute_intersection)
              return true;
#ifdef ANISO_USE_EXACT
            if(!compute_exact_dual_intersection(f, surface_delaunay_ball_center))
              std::cerr << "Warning : is_restricted exact can't find an intersection point.\n";
#else
            if(!compute_dual_intersection(f, surface_delaunay_ball_center))
              std::cerr << "Warning : is_restricted can't find an intersection point.\n";
#endif
            return true;
          }
          else return false;
        }
        return false;
      }

      bool is_topological_disk() const
      {
        update_is_topological_disk();
        return m_is_topological_disk;
      }

      void update_is_topological_disk() const
      {
        if(m_is_valid_topo_disk)
          return;

        CGAL_PROFILER("[Update is_topological_disk]");
        if(!is_surface_star() || this->dimension() < 3)
        {
          m_is_topological_disk = false;
          return;
        }
        Facet_set_iterator fit = this->begin_restricted_facets();
        Facet_set_iterator fend = this->end_restricted_facets();
        if(fit == fend)
        {
          m_is_topological_disk = false;// no restricted facet --> false
          return;
        }

        std::set<Index> vertices_around;
        for(; fit != fend; fit++)
        {
          Facet f = *fit;
          for(unsigned int i = 1; i < 4; i++)
          {
            Vertex_handle v = f.first->vertex((f.second + i)%4);
            if(v != m_center)
            {
              Index j = v->info();
              if(vertices_around.find(j) == vertices_around.end())
                vertices_around.insert(j);
              else vertices_around.erase(j);
            }
          }
        }
        m_is_topological_disk = vertices_around.empty();
        m_is_valid_topo_disk = true;
      }

      inline bool is_same(const Cell_handle &c, const Cell_handle &d)
      {
        int cids[4], dids[4];
        for (int i = 0; i < 4; i++)
        {
          cids[i] = c->vertex(i)->info();
          dids[i] = d->vertex(i)->info();
        }
        return is_same_ids<4>(cids, dids);
      }
      inline bool is_same(const Cell_handle &c, int *vertices)
      {
        int cids[4];
        for (int i = 0; i < 4; i++)
          cids[i] = c->vertex(i)->info();
        return is_same_ids<4>(cids, vertices);
      }

      bool has_cell(const Cell_handle &cell)
      {
        int cids[4], dids[4];
        for (int i = 0; i < 4; i++)
          cids[i] = cell->vertex(i)->info();

        Cell_handle_handle ci = begin_star_cells();
        Cell_handle_handle cend = end_star_cells();
        for (; ci != cend; ci++)
        {
          for (int i = 0; i < 4; i++)
            dids[i] = (*ci)->vertex(i)->info();
          if (is_same_ids<4>(cids, dids))
            return true;
        }
        return false;
      }
      bool has_cell_ref(int *vertices, Cell_handle &cell)
      {
        int dids[4];
        Cell_handle_handle ci = begin_star_cells();
        Cell_handle_handle cend = end_star_cells();
        for (; ci != cend; ci++) {
          for (int i = 0; i < 4; i++)
            dids[i] = (*ci)->vertex(i)->info();
          if (is_same_ids<4>(vertices, dids))
          {
            cell = *ci;
            return true;
          }
        }
        return false;
      }

      bool has_facet_ref(int *vertices, Facet &facet) const
      {
        int dids[3];
        Facet_set_iterator fi = begin_restricted_facets();
        Facet_set_iterator fiend = end_restricted_facets();
        for (; fi != fiend; fi++)
        {
          for (int i = 1; i <= 3; i++)
            dids[i - 1] = fi->first->vertex((fi->second + i) % 4)->info();
          if (is_same_ids<3>(vertices, dids))
          {
            facet = *fi;
            return true;
          }
        }
        return false;
      }

      bool is_same(const Facet &c, const Facet &d) const
      {
        int cids[3], dids[3];
        int cidp = 0, didp = 0;
        for (int i = 0; i < 4; i++)
        {
          if (i != c.second)
            cids[cidp++] = c.first->vertex(i)->info();
          if (i != d.second)
            dids[didp++] = d.first->vertex(i)->info();
        }
        // sorting and compare
        for (int i = 2; i >= 0; i--)
        {
          for (int j = 1; j <= i; j++)
          {
            if (cids[j] < cids[j - 1])
              std::swap(cids[j], cids[j - 1]);
            if (dids[j] < dids[j - 1])
              std::swap(dids[j], dids[j - 1]);
          }
          if (cids[i] != dids[i])
            return false;
        }
        return true;
      }

      bool has_facet(const Facet &facet) const
      {
        Facet_set_iterator fi = begin_restricted_facets(); // does update caches
        Facet_set_iterator fend = end_restricted_facets();
        for (; fi != fend; fi++)
          if (is_same(*fi, facet))
            return true;
        return false;
      }

      bool has_vertex(const int i)
      {
        typename Base::Finite_vertices_iterator vit = this->finite_vertices_begin();
        typename Base::Finite_vertices_iterator vend = this->finite_vertices_end();
        for(; vit != vend; vit++)
          if(vit->info() == i)
            return true;
        return false;
      }

      void print_vertices(const bool print_points = false) const
      {
        std::cout << "\t Star_" << index_in_star_set();
        std::cout << " vertices ("<<(this->number_of_vertices())<<" v.)\t (";
        typename std::set<Vertex_handle> vertices;
        if(begin_restricted_facets() == end_restricted_facets())
          std::cout << "empty";
        else
        {
          for(Facet_set_iterator fit = begin_restricted_facets();
              fit != end_restricted_facets();
              fit++)
          {
            Facet f = *fit;
            vertices.insert(f.first->vertex((f.second + 1) % 4));
            vertices.insert(f.first->vertex((f.second + 2) % 4));
            vertices.insert(f.first->vertex((f.second + 3) % 4));
          }
          for(typename std::set<Vertex_handle>::iterator vit = vertices.begin();
              vit != vertices.end(); vit++)
          {
            if(print_points) std::cout << "\n\t";
            std::cout << (*vit)->info() << " ";
            if(print_points) std::cout << " : " << (*vit)->point();
          }
        }
        std::cout << ")\t(";
        typename Base::Finite_vertices_iterator vit = this->finite_vertices_begin();
        typename Base::Finite_vertices_iterator vend = this->finite_vertices_end();
        for(; vit != vend; vit++)
          std::cout << vit->info() << " ";
        std::cout << ")\n";
      }

      void print_faces() const
      {
        std::cout << "\t Star_" << index_in_star_set() << " faces :\n";
        Facet_set facets;
        Cell_handle_handle ci = begin_finite_star_cells();
        Cell_handle_handle cend = end_finite_star_cells();
        for (; ci != cend; ci++)
        {
          int center_index = (*ci)->index(m_center);
          for(int i = 0; i < 4; i++)
          {
            if(i == center_index) continue;
            facets.insert(this->make_canonical(Facet(*ci, i)));
          }
        }
        Facet_set_iterator it = facets.begin();
        for(; it != facets.end(); it++)
        {
          Facet f = *it;
          std::cout << "\t\t f(" << f.first->vertex((f.second+1)%4)->info()
                        << " " << f.first->vertex((f.second+2)%4)->info()
                        << " " << f.first->vertex((f.second+3)%4)->info() << ") ";
          if(is_restricted(f))
            std::cout << "restricted";
          std::cout << ",\n";
        }
      }

      bool is_inside_bbox(const Point_3& p) const
      {
        update_bbox();
        return p.x() >= m_bbox.xmin() && p.x() <= m_bbox.xmax()
            && p.y() >= m_bbox.ymin() && p.y() <= m_bbox.ymax()
            && p.z() >= m_bbox.zmin() && p.z() <= m_bbox.zmax();
      }

      inline bool is_conflicted(const TPoint_3 &tp,
                                Cell_handle &in_which_cell) const
      {
        if(Base::dimension() <= 1)
          return true;
        else if(!is_inside_bbox(m_metric.inverse_transform(tp)))
          return false;
        else if(is_topological_disk())
          return is_in_a_surface_delaunay_ball(tp, in_which_cell);
        else
          return is_in_a_volume_delaunay_ball(tp, in_which_cell);
      }

      inline bool is_in_a_surface_delaunay_ball(const TPoint_3& tp,
                                               Cell_handle& in_which_cell) const
      {
        CGAL_PROFILER("[is_in_a_surface_delaunay_ball]");
        Facet_set_iterator fi = begin_restricted_facets();
        Facet_set_iterator fend = end_restricted_facets();
        for(; fi != fend; fi++)
          if(is_in_surface_delaunay_ball(tp, *fi, in_which_cell))
            return true;
        return false;
      }

      inline bool is_in_a_volume_delaunay_ball(const TPoint_3& tp,
                                               Cell_handle& in_which_cell) const
      {
        CGAL_PROFILER("[is_in_a_volume_delaunay_ball]");
        if(this->dimension() < 3) return true;
        Cell_handle_handle ci = begin_star_cells();
        Cell_handle_handle cend = end_star_cells();
        for (; ci != cend; ci++)
        {
          if(is_in_delaunay_ball(tp, *ci))
          {
            in_which_cell = *ci;
            return true;
          }
        }
        return false;
      }

private :
      bool is_in_delaunay_ball(const TPoint_3& tp,
                               const Cell_handle& c) const
      {
        CGAL_PROFILER("[conflict volumeDB test]");
        return (Base::side_of_sphere(c, tp) == CGAL::ON_BOUNDED_SIDE);
        // i.e. ON_BOUNDARY or ON_BOUNDED_SIDE
      }

      bool is_in_surface_delaunay_ball(const TPoint_3& tp,
                                       const Facet& f,
                                       Cell_handle& in_which_cell) const
      {
        CGAL_PROFILER("[conflict surfaceDB test]");
        if(Base::dimension() < 3)
          return true;

        typename Traits::Compute_squared_distance_3 csd
          = m_traits->compute_squared_distance_3_object();
        Point_3 center;
#ifdef ANISO_USE_EXACT
        compute_exact_dual_intersection(f, center);
#else
        compute_dual_intersection(f, center);
#endif
        TPoint_3 tcenter = m_metric.transform(center);

        FT sq_radius = csd(tcenter, m_center->point()); // M_star
        if(csd(tcenter, tp) >= sq_radius)
          return false;

        Cell_handle c1 = f.first;
        Cell_handle c2 = f.first->neighbor(f.second);
        if(Base::side_of_sphere(c1, tp) == CGAL::ON_BOUNDED_SIDE)
        {
          in_which_cell = c1;
          return true;
        }
        else if(Base::side_of_sphere(c2, tp) == CGAL::ON_BOUNDED_SIDE)
        {
          in_which_cell = c2;
          return true;
        }
        return false;
      }

public:
      void remove(Vertex_handle v)
      {
        Base::remove(v);
        invalidate_cache();
      }
      void remove(const Index& i)
      {
        typename Base::Finite_vertices_iterator vit = this->finite_vertices_begin();
        typename Base::Finite_vertices_iterator vend = this->finite_vertices_end();
        for(; vit != vend; vit++)
        {
          if(vit->info() == i)
          {
            this->remove(vit);
            break;//warning : vit is broken
          }
        }
      }

  private:
      // warning : point should be transformed before this insert function is called
      Vertex_handle insert(const TPoint_3 &tp, const Cell_handle &cell)
      {
        Vertex_handle retval = Base::insert(tp, cell);
        invalidate_cache();
        return retval;
      }
      // warning : point should be transformed before this insert function is called
      Vertex_handle insert(const TPoint_3 &tp)
      {
        Vertex_handle retval = Base::insert(tp);
        invalidate_cache();
        return retval;
      }
      // warning : point should be transformed before this insert function is called
      Vertex_handle insert_to_star_(const TPoint_3 &tp,
                                    const bool conditional)
      {
        if(Base::dimension() < 3 || !conditional)
          return this->insert(tp);

        Cell_handle cell;
        if(is_conflicted(tp, cell))
          return this->insert(tp, cell);

        else // point is not in conflict with cells incident to center
          return Vertex_handle();
      }

  public:
      Vertex_handle insert_to_star(const Point_3 &p,
                                   const int id,
                                   const bool conditional)
      {
        Vertex_handle v = insert_to_star_(m_metric.transform(p), conditional);
        if(v != center() && v != Vertex_handle())
          v->info() = id;
        return v;
      }

      int simulate_insert_to_star(const Point_3& p, const int id)
      {
        CGAL_HISTOGRAM_PROFILER("V", this->number_of_vertices());
        TPoint_3 tp = m_metric.transform(p);
        bool found_vertex = false;
        Vertex_handle vh;

#ifndef ANISO_BRUTE_FORCE_SIMULATE_INSERT_TO_STAR
        int li, lj;
        typename Base::Locate_type lt;
        Cell_handle c = Base::locate(tp, lt, li, lj, m_center->cell());
        if(lt ==  Base::VERTEX){
          found_vertex = true;
          vh = c->vertex(li);
        }
#else
        for(typename Base::Finite_vertices_iterator it =  this->finite_vertices_begin();
            (! found_vertex) && (it != this->finite_vertices_end());
            ++it){
          if(p == it->point()){
            found_vertex = true;
            vh = it;
          }
        }
#endif
        if(found_vertex) // already in star
        {
          std::cout << "Warning : simulate_insert_to_star re-inserts same point.\n";
          return vh->info(); //
        }
        else             // in conflict
        {
          Cell_handle cell;
          if(is_conflicted(tp, cell))
            return id;
        }
        return -1; // no conflict
      }

      template<typename BoundaryFacetsOutputIterator>
      bool find_conflicts(const Point_3& p, BoundaryFacetsOutputIterator oit)
      {
        int dim = Base::dimension();
        if(dim <= 1)
          return false;
        else if(dim == 2)
        {
          Facet_set_iterator fit = begin_restricted_facets();
          Facet_set_iterator fend = end_restricted_facets();
          for(; fit != fend; fit++)
            if(is_restricted(*fit))
              *oit++ = *fit;
        }
        else
        {
          TPoint_3 tp = m_metric.transform(p);
          Cell_handle ch; //warning : dimension should be 3 to use find_conflicts!
          if(is_in_a_volume_delaunay_ball(tp, ch))
          {
            Base::find_conflicts(tp, ch, oit, Emptyset_iterator());
            return true;
          }
        }
        return false;
      }

      std::size_t clean() //remove useless vertices
      {
        typedef typename std::pair<TPoint_3, int> PPoint;
        std::vector<PPoint> backup;
        Vertex_handle_handle vit = begin_neighboring_vertices();
        Vertex_handle_handle vend = end_neighboring_vertices();
        // backup star vertices
        //center first
        backup.push_back(PPoint(m_center->point(), index_in_star_set()));
        //other vertices
        for(; vit != vend; vit++)
          if(!is_infinite_vertex(*vit))
            backup.push_back(PPoint((*vit)->point(), (*vit)->info()));

        std::size_t nbv = this->number_of_vertices();
        this->clear();

        //center first
        m_center = insert(backup[0].first);
        m_center->info() = backup[0].second;
        //other vertices
        for(std::size_t i = 1; i < backup.size(); i++)
          insert(backup[i].first)->info() = backup[i].second;
        this->infinite_vertex()->info() = index_of_infinite_vertex;

        invalidate_cache();
        return (nbv - backup.size());
      }

      // Note (JT) : this function is commented because it is not used anywhere
      // Note 2 : p should be transformed before anything
      //bool can_form_sliver(const Point_3 &p, const FT &squared_radius_bound)
      //{
      //  Cell_handle_handle ci = begin_star_cells();
      //  Cell_handle_handle cend = end_star_cells();
      //  for (; ci != cend; ci++)
      //  {
      //    if (Base::side_of_sphere((*ci), p) == CGAL::ON_UNBOUNDED_SIDE)
      //      continue;
      //    Cell_handle c = *ci;
      //    Point_3 points[3];
      //    Point_3 centerpoint = m_center->point();
      //    int k = 0;
      //    for (int i = 0; i < 4; i++)
      //    {
      //      if (c->vertex(i) != m_center)
      //        points[k++] = c->vertex(i)->point();
      //    }
      //    assert(k == 3);
      //    for (int i = 0; i < 3; i++)
      //    {
      //      FT sliver_overflow = m_criteria->sliverity_overflow(
      //        centerpoint, points[i], points[(i + 1) % 3], p);
      //      if (sliver_overflow <= 0.0)
      //        continue;
      //      if (m_criteria->compute_volume(
      //        centerpoint, points[i], points[(i + 1) % 3], p) < 1e-7)
      //        continue;
      //      // here, we have to check coplanar and colinear cases
      //      if (m_criteria->compute_squared_circumradius(centerpoint, points[i],
      //        points[(i + 1) % 3], p) >= squared_radius_bound)
      //        continue;
      //      return true;
      //    }
      //  }
      //  return false;
      //}

      inline void get_inverse_transformed_points(Point_3 *ps, const Facet &facet) const
      {
        int k = 0;
        for (int i = 0; i < 4; i++)
          if (i != facet.second)
            ps[k++] = m_metric.inverse_transform(facet.first->vertex(i)->point());
      }

      inline void get_points(TPoint_3 *ps, const Facet &facet) const
      {
        int k = 0;
        for (int i = 0; i < 4; i++)
          if (i != facet.second)
            ps[k++] = facet.first->vertex(i)->point();
      }
      inline void get_points(TPoint_3 *ps, Cell_handle cell) const
      {
        int k = 0;
        for (int i = 0; i < 4; i++)
          ps[k++] = cell->vertex(i)->point();
      }


      TPoint_3 compute_circumcenter(const Facet &facet) const
      {
#ifdef ANISO_USE_CC_EXACT
        return back_from_exact(compute_exact_circumcenter(facet));
#else
        return compute_exact_circumcenter(facet);
#endif
      }

      Exact_TPoint_3 compute_exact_circumcenter(const Facet& facet) const
      {
        TPoint_3 ps[3];
        get_points(ps, facet);
#ifdef ANISO_USE_CC_EXACT
        return m_traits->construct_circumcenter_3_object()(
                    to_exact(ps[0]), to_exact(ps[1]), to_exact(ps[2]));
#else
        return m_traits->construct_circumcenter_3_object()(ps[0], ps[1], ps[2]);
#endif
      }

      Exact_TPoint_3 compute_exact_circumcenter(Cell_handle cell) const
      {
#ifdef ANISO_USE_CC_EXACT
        TPoint_3 ps[4];
        get_points(ps, cell);
        return m_traits->construct_circumcenter_3_object()(
                 to_exact(ps[0]), to_exact(ps[1]), to_exact(ps[2]), to_exact(ps[3]));
#else
        return cell->circumcenter(*m_traits);
#endif
      }
      TPoint_3 compute_circumcenter(Cell_handle cell) const
      {
        return compute_exact_circumcenter(cell);
      }

      template<typename TPoint> //exact or not
      TPoint compute_circumcenter(const TPoint &p0, const TPoint &p1,
                                  const TPoint &p2, const TPoint &p3) const
      {
        return m_traits->construct_circumcenter_3_object()(p0, p1, p2, p3);
      }

      bool is_encroached(const Point_3 &testp, typename Constrain_surface::EdgeIterator &e)
      {
        typename Traits::Side_of_bounded_sphere_3 o = m_traits->side_of_bounded_sphere_3_object();
        typename Constrain_surface::EdgeIterator ei = m_pConstrain->edge_begin();
        typename Constrain_surface::EdgeIterator eiend = m_pConstrain->edge_end();
        for (; ei != eiend; ei++)
        {
          Point_3 p1 = ei->first;
          Point_3 p2 = ei->second;
          if (CGAL::abs(p1.x() - testp.x()) + CGAL::abs(p1.y() - testp.y())
            + CGAL::abs(p1.z() - testp.z()) < 1e-3)
            continue;
          if (CGAL::abs(p2.x() - testp.x()) + CGAL::abs(p2.y() - testp.y())
            + CGAL::abs(p2.z() - testp.z()) < 1e-3)
            continue;
          if (o(p1, p2, testp) != CGAL::ON_UNBOUNDED_SIDE)
          {
            e = ei;
            return true;
          }
        }
        return false;
      }

      typename K::Object_3 constrain_ray_intersection(const Point_3 &p1, //source
                                                      const Point_3 &p2) const //target
      {
        if(p1 == p2)
        {
#ifdef ANISO_VERBOSE
          std::cout << "CONSTRAIN_RAY_INTERSECTION : source cannot be == target ";
          std::cout << " ("<< p1 << ")" << std::endl;
#endif
          return typename K::Object_3();
        }
        Point_3 res;
        typename K::Ray_3 ray(p1, p2);
        if (m_pConstrain->intersection(make_object(ray)).assign(res))
          return make_object(res);
        else return typename K::Object_3();
      }

      typename K::Object_3 constrain_segment_intersection(const Point_3 &p1, //point1
                                                          const Point_3 &p2) const //point2
      {
        if(p1 == p2)
        {
#ifdef ANISO_VERBOSE
          std::cout << "CONSTRAIN_SEGMENT_INTERSECTION : source cannot be == target ";
          std::cout << " ("<< p1 << ")" << std::endl;
#endif
          return typename K::Object_3();
        }
        Point_3 res;
        typename K::Segment_3 seg(p1, p2);
        if (m_pConstrain->intersection(make_object(seg)).assign(res))
          return make_object(res);
        else return typename K::Object_3();
      }

      bool is_restricted_2_in_3(const Facet& f,
                                Point_3& surface_delaunay_ball_center,
                                const bool exact = true) const
      {
        if(Base::dimension() != 2)
          return false;
        CGAL_PROFILER("[is_restricted_2_in_3]");

        Triangle tr = m_metric.inverse_transform(Base::triangle(f));
        Vector_3 n = tr.supporting_plane().orthogonal_vector();
        if(exact)
        {
#ifdef ANISO_USE_CC_EXACT
          Point_3 fc = back_from_exact(m_metric.inverse_transform(compute_exact_circumcenter(f)));
#else
          Point_3 fc = m_metric.inverse_transform(compute_circumcenter(f));
#endif
          return (constrain_ray_intersection(fc, fc-n).assign(surface_delaunay_ball_center)
               || constrain_ray_intersection(fc, fc+n).assign(surface_delaunay_ball_center));
        }
        else
        {
          Point_3 fc = m_metric.inverse_transform(compute_circumcenter(f));
          return (constrain_ray_intersection(fc, fc-n).assign(surface_delaunay_ball_center)
               || constrain_ray_intersection(fc, fc+n).assign(surface_delaunay_ball_center));
        }
      }

      typename K::Object_3 dual(const Facet& facet) const
      {
        if(Base::dimension() == 2)
        {
          Point_3 fc = m_metric.inverse_transform(compute_circumcenter(facet));
          Triangle tr = m_metric.inverse_transform(Base::triangle(facet));
          return make_object(Line_3(fc, tr.supporting_plane().orthogonal_vector()));
        }

        Point_3 ps[3];
        get_inverse_transformed_points(ps, facet);

        Cell_handle c1 = facet.first;
        Cell_handle c2 = c1->neighbor(facet.second);
        bool f1 = !is_infinite(c1);
        bool f2 = !is_infinite(c2);

        Point_3 fc = m_metric.inverse_transform(compute_circumcenter(facet));
        Vector_3 t_n = this->dual_support(c1, facet.second).to_vector();
        Vector_3 n = m_metric.inverse_transform(t_n);
        n = std::sqrt(1./(n*n)) * n;

        if(f1)
        {
          if (f2)
          {
            Point_3 cp1 = m_metric.inverse_transform(c1->circumcenter(*(m_traits)));
            Point_3 cp2 = m_metric.inverse_transform(c2->circumcenter(*(m_traits)));
            return make_object(Segment(cp1, cp2));
          }
          else // !f2
          {
            Point_3 cp = m_metric.inverse_transform(c1->circumcenter(*(m_traits)));
            Point_3 ps3 = m_metric.inverse_transform(c1->vertex(facet.second)->point());

            CGAL::Orientation o1 = CGAL::orientation(ps[0], ps[1], ps[2], ps3);
            CGAL::Orientation o2 = CGAL::orientation(ps[0], ps[1], ps[2], fc + n);
            CGAL::Orientation o3 = CGAL::orientation(ps[0], ps[1], ps[2], cp);

            if(o3 == COPLANAR) //fc = cp
            {
              if(fc != cp)
              {
                std::cerr << "fc != cp with cp on the facet..." << std::endl;
                std::cerr << "fc : " << fc << std::endl;
                std::cerr << "cp : " << cp << std::endl;
              }
              if(o1 == o2) //n points towards ps3
                n = -n;
              return make_object(Ray_3(cp, cp + n));
            }
            else if(o1 == o3)
              return make_object(Ray_3(cp, fc));
            else if(o1 != o3)
              return make_object(Ray_3(cp, Point_3(cp + Vector_3(fc, cp))));
          }
        }
        else // !f1
        {
          if (f2)
          {
            Point_3 cp = m_metric.inverse_transform(c2->circumcenter(*(m_traits)));
            Point_3 c2_ps3 = m_metric.inverse_transform(c2->vertex(c2->index(c1))->point()); //4th pt of c2

            CGAL::Orientation o1 = CGAL::orientation(ps[0], ps[1], ps[2], c2_ps3);
            CGAL::Orientation o2 = CGAL::orientation(ps[0], ps[1], ps[2], fc + n);
            CGAL::Orientation o3 = CGAL::orientation(ps[0], ps[1], ps[2], cp);

            if(o3 == COPLANAR) //fc = cp
            {
              if(fc != cp)
              {
                std::cerr << "fc != cp with cp on the facet..." << std::endl;
                std::cerr << "fc : " << fc << std::endl;
                std::cerr << "cp : " << cp << std::endl;
              }
              if(o1 == o2) //n points towards c2_ps3
                n = -n;
              return make_object(Ray_3(cp, cp + n));
            }
            else if(o1 == o3)
              return make_object(Ray_3(cp, fc));
            else if(o1 != o3)
              return make_object(Ray_3(cp, Point_3(cp + Vector_3(fc, cp))));
          }
        }
        return typename K::Object_3();
      }

      bool compute_exact_dual_intersection(const Facet &facet,
                                           Point_3& p,
                                           const bool use_cache = false,
                                           const bool verbose = true,
                                           const bool super_verbose = false) const
      {
#ifdef USE_ANISO_TIMERS
        std::clock_t start = clock();
#endif
        if(!is_restricted(facet))// point is not computed here
          return false;

        if(super_verbose)
          facet_indices(facet);

        if(use_cache && facet.first->is_facet_visited(facet.second))
        {
          p = facet.first->get_facet_surface_center(facet.second);
          return true;
        }
        bool ret_val = false;
        CGAL_PROFILER("[compute_exact_dual_intersection]");
        if(Base::dimension() == 2)
          ret_val = compute_exact_dual_intersection_2(facet, p);
        else
        {
          Point_3 ps[3];
          get_inverse_transformed_points(ps, facet);

          Cell_handle c1 = facet.first;
          Cell_handle c2 = c1->neighbor(facet.second);
          bool f1 = !is_infinite(c1);
          bool f2 = !is_infinite(c2);

          if(super_verbose)
            std::cerr << "(case " << f1 << " " << f2 << ")";

          Exact_Point_3 fc = m_metric.inverse_transform(compute_exact_circumcenter(facet));
          Exact_Vector_3 t_n = to_exact(this->dual_support(c1, facet.second).to_vector());
          Exact_Vector_3 n = m_metric.inverse_transform(t_n);
          n = std::sqrt(1./(n*n)) * n;

          if(!f1 && f2)
          {
            Cell_handle tmp = c2;
            c2 = c1;    f2 = false;
            c1 = tmp;   f1 = true;
          }

          if(f1)
          {
            if (f2)
            {
              Exact_Point_3 cp1 = m_metric.inverse_transform(compute_exact_circumcenter(c1));
              Exact_Point_3 cp2 = m_metric.inverse_transform(compute_exact_circumcenter(c2));
              if(m_pConstrain->intersection(back_from_exact(cp1), back_from_exact(cp2)).assign(p))
                ret_val = true;
            }
            else // !f2
            {
              Exact_Point_3 cp = m_metric.inverse_transform(compute_exact_circumcenter(c1));
              Point_3 ps3 = m_metric.inverse_transform(c1->vertex(facet.second)->point());

              CGAL::Orientation o1 = CGAL::orientation(to_exact(ps[0]), to_exact(ps[1]),
                                                       to_exact(ps[2]), to_exact(ps3));
              CGAL::Orientation o2 = CGAL::orientation(to_exact(ps[0]), to_exact(ps[1]),
                                                       to_exact(ps[2]), fc + n);
              CGAL::Orientation o3 = CGAL::orientation(to_exact(ps[0]), to_exact(ps[1]),
                                                       to_exact(ps[2]), cp);

              if(o3 == COPLANAR) //fc = cp
              {
                if(fc != cp)
                {
                  std::cerr << "fc != cp with cp on the facet...[exact "<<f1<< " "<<f2 <<"]" << std::endl;
                  std::cerr << "fc : " << fc << std::endl;
                  std::cerr << "cp : " << cp << std::endl;
                }
                if(o1 == o2) //n points towards ps3
                  n = -n;
                Exact_Point_3 ep = cp + n;
                if(constrain_ray_intersection(back_from_exact(cp), back_from_exact(ep)).assign(p))
                  ret_val = true;
              }
              else if(o1 == o3 && constrain_ray_intersection(back_from_exact(cp), back_from_exact(fc)).assign(p))
                ret_val = true;
              else if(o1 != o3)
              {
                Exact_Point_3 target(cp + Exact_Vector_3(fc, cp));
                if(constrain_ray_intersection(back_from_exact(cp), back_from_exact(target)).assign(p))
                  ret_val = true;
              }
              if(super_verbose)
              {
                std::cerr << "\t(o " << o1 << " " << o2 << " " << o3 << ")\n";
                std::cerr << "\t(cp " << cp << ")\n";
                std::cerr << "\t(fc " << fc << ")\n";
                std::cerr << "\t(returns " << ret_val << ")\n";
              }
            }
          }

          if(ret_val)//do this in dimension 3 only (mirror_facet could crash, otherwise)
            set_facet_cache(facet, p);
          else if(verbose)
          {
            std::cerr.precision(15);
            std::cerr << "Oops! no intersection [exact]" << std::endl;
            std::cerr << "Case:  " << f1 << " " << f2 << std::endl;
            std::cerr << "Star:  " << m_center->info() << std::endl;
            std::cerr << "Facet 1: " << std::endl;
            std::cerr << c1->vertex((facet.second + 1) % 4)->info() << " ";
            std::cerr << to_exact(m_metric.inverse_transform(c1->vertex((facet.second + 1) % 4)->point())) << std::endl;
            std::cerr << c1->vertex((facet.second + 2) % 4)->info() << " ";
            std::cerr << to_exact(m_metric.inverse_transform(c1->vertex((facet.second + 2) % 4)->point())) << std::endl;
            std::cerr << c1->vertex((facet.second + 3) % 4)->info() << " ";
            std::cerr << to_exact(m_metric.inverse_transform(c1->vertex((facet.second + 3) % 4)->point())) << std::endl;
            std::cerr << c1->vertex(facet.second)->info() << " ";
            std::cerr << to_exact(m_metric.inverse_transform(c1->vertex(facet.second)->point()));
            std::cerr << " (.second)" << std::endl;

            std::cerr << "Cell 2: " << std::endl;
            std::cerr << c2->vertex(0)->info() << " ";
            std::cerr << to_exact(m_metric.inverse_transform(c2->vertex(0)->point())) << std::endl;
            std::cerr << c2->vertex(1)->info() << " ";
            std::cerr << to_exact(m_metric.inverse_transform(c2->vertex(1)->point())) << std::endl;
            std::cerr << c2->vertex(2)->info() << " ";
            std::cerr << to_exact(m_metric.inverse_transform(c2->vertex(2)->point())) << std::endl;
            std::cerr << c2->vertex(3)->info() << " ";
            std::cerr << to_exact(m_metric.inverse_transform(c2->vertex(3)->point())) << std::endl;
          }
          if(super_verbose)
            std::cerr << "Normal ("<< n <<")" << std::endl;
        }
#ifdef USE_ANISO_TIMERS
        m_compute_dual_intersection_timer += (clock()-start+0.) / ((double)CLOCKS_PER_SEC);
#endif
        return ret_val;
      }

      void set_facet_cache(const Facet& facet,
                           const Point_3& p) const
      {
        Facet mf = this->mirror_facet(facet);
        facet.first->set_facet_surface_center(facet.second, p);
        mf.first->set_facet_surface_center(mf.second, p);
        facet.first->set_facet_visited(facet.second);
        mf.first->set_facet_visited(mf.second);
      }

      //compute_exact_dual_intersection if triangulation is of dimension 2
      bool compute_exact_dual_intersection_2(const Facet &facet,
                                             Point_3& p) const
      {
        bool ret_val = false;
        Point_3 p1, p2;
        TPoint_3 fc = back_from_exact(m_metric.inverse_transform(compute_exact_circumcenter(facet)));
        Triangle tr = m_metric.inverse_transform(Base::triangle(facet));
        Vector_3 n = tr.supporting_plane().orthogonal_vector();

        bool b1 = false; bool b2 = false;
        if(constrain_ray_intersection(fc, fc-n).assign(p1)) b1 = true;
        if(constrain_ray_intersection(fc, fc+n).assign(p2)) b2 = true;
        if(b1 && b2)
        {
          if (CGAL::squared_distance(fc, p1) > CGAL::squared_distance(fc, p2))
            p = p2;
          else p = p1;
          ret_val = true;
        }
        else if(b1 && !b2) { p = p1 ; ret_val = true; }
        else if(!b1 && b2) { p = p2 ; ret_val = true; }
        return ret_val;
      }


      bool compute_dual_intersection(const Facet &facet,
                                     Point_3& p,
                                     const bool use_cache = true,
                                     const bool verbose = true,
                                     const bool super_verbose = false) const
      {
#ifdef USE_ANISO_TIMERS
        std::clock_t start = clock();
#endif
        if(!is_restricted(facet))// point is not computed here
          return false;

        if(use_cache && facet.first->is_facet_visited(facet.second))
        {
          p = facet.first->get_facet_surface_center(facet.second);
          return true;
        }
        bool ret_val = false;
        CGAL_PROFILER("[compute_dual_intersection]");

        if(Base::dimension() == 2)
        {
          Point_3 p1, p2;
          Point_3 fc = m_metric.inverse_transform(compute_circumcenter(facet));
          Triangle tr = m_metric.inverse_transform(Base::triangle(facet));
          Vector_3 n = tr.supporting_plane().orthogonal_vector();

          bool b1 = false;
          bool b2 = false;
          if(constrain_ray_intersection(fc, fc - n).assign(p1)) b1 = true;
          if(constrain_ray_intersection(fc, fc + n).assign(p2)) b2 = true;
          if(b1 && b2)
          {
            if (CGAL::squared_distance(fc, p1) > CGAL::squared_distance(fc, p2))
              p = p2;
            else p = p1;
            ret_val = true;
          }
          else if(b1 && !b2) { p = p1 ; ret_val = true; }
          else if(!b1 && b2) { p = p2 ; ret_val = true; }
        }
        else
        {
          Point_3 ps[3];
          get_inverse_transformed_points(ps, facet);

          Cell_handle c1 = facet.first;
          Cell_handle c2 = c1->neighbor(facet.second);
          bool f1 = !is_infinite(c1);
          bool f2 = !is_infinite(c2);

          if(super_verbose)
            std::cerr << "(case " << f1 << " " << f2 << ")";

          Point_3 fc = m_metric.inverse_transform(compute_circumcenter(facet));
          Vector_3 t_n = this->dual_support(c1, facet.second).to_vector();
          Vector_3 n = m_metric.inverse_transform(t_n);
          n = std::sqrt(1./(n*n)) * n;

          if(!f1 && f2)
          {
            Cell_handle tmp = c2;
            c2 = c1;    f2 = false;
            c1 = tmp;   f1 = true;
          }
          
          if(f1)
          {
            if (f2)
            {
              Point_3 cp1 = m_metric.inverse_transform(c1->circumcenter(*(m_traits)));
              Point_3 cp2 = m_metric.inverse_transform(c2->circumcenter(*(m_traits)));

              if(cp1 == cp2){
                p = cp1;
                ret_val = true;
              }
              else if(m_pConstrain->intersection(cp1, cp2).assign(p))
                ret_val = true;
            }
            else // !f2
            {
              Point_3 cp = m_metric.inverse_transform(c1->circumcenter(*(m_traits)));
              Point_3 ps3 = m_metric.inverse_transform(c1->vertex(facet.second)->point());

              CGAL::Orientation o1 = CGAL::orientation(ps[0], ps[1], ps[2], ps3);
              CGAL::Orientation o2 = CGAL::orientation(ps[0], ps[1], ps[2], fc + n);
              CGAL::Orientation o3 = CGAL::orientation(ps[0], ps[1], ps[2], cp);

              if(o3 == COPLANAR) //fc = cp
              {
                if(fc != cp)
                {
                  std::cerr << "fc != cp with cp on the facet..." << std::endl;
                  std::cerr << "fc : " << fc << std::endl;
                  std::cerr << "cp : " << cp << std::endl;
                }
                if(o1 == o2) //n points towards ps3
                  n = -n;
                if(constrain_ray_intersection(cp, cp + n).assign(p))
                  ret_val = true;
              }
              else if(o1 == o3 && constrain_ray_intersection(cp, fc).assign(p))
                ret_val = true;
              else if(o1 != o3 && constrain_ray_intersection(cp, Point_3(cp + Vector_3(fc, cp))).assign(p))
                ret_val = true;
              if(super_verbose)
              {
                std::cerr << "(o " << o1 << " " << o2 << " " << o3 << ")";
                std::cerr << "(cp " << cp << ")";
              }
            }
          }

          if(ret_val)//do this in dimension 3 only (mirror_facet could crash, otherwise)
            set_facet_cache(facet, p);
          else if(verbose)
          {
            std::cerr.precision(15);
            std::cerr << "Oops! no intersection!" << std::endl;
            std::cerr << "Case:  " << f1 << " " << f2 << std::endl;
            std::cerr << "Star:  " << m_center->info() << std::endl;
            std::cerr << "Facet 1: " << std::endl;
            std::cerr << c1->vertex((facet.second + 1) % 4)->info() << " ";
            std::cerr << m_metric.inverse_transform(c1->vertex((facet.second + 1) % 4)->point()) << std::endl;
            std::cerr << c1->vertex((facet.second + 2) % 4)->info() << " ";
            std::cerr << m_metric.inverse_transform(c1->vertex((facet.second + 2) % 4)->point()) << std::endl;
            std::cerr << c1->vertex((facet.second + 3) % 4)->info() << " ";
            std::cerr << m_metric.inverse_transform(c1->vertex((facet.second + 3) % 4)->point()) << std::endl;
            std::cerr << c1->vertex(facet.second)->info() << " ";
            std::cerr << m_metric.inverse_transform(c1->vertex(facet.second)->point()) << " (.second)" << std::endl;
            std::cerr << "Cell 2: " << std::endl;
            std::cerr << c2->vertex(0)->info() << " ";
            std::cerr << m_metric.inverse_transform(c2->vertex(0)->point()) << std::endl;
            std::cerr << c2->vertex(1)->info() << " ";
            std::cerr << m_metric.inverse_transform(c2->vertex(1)->point()) << std::endl;
            std::cerr << c2->vertex(2)->info() << " ";
            std::cerr << m_metric.inverse_transform(c2->vertex(2)->point()) << std::endl;
            std::cerr << c2->vertex(3)->info() << " ";
            std::cerr << m_metric.inverse_transform(c2->vertex(3)->point()) << std::endl;
          }
          if(super_verbose)
            std::cerr << "Normal ("<< n <<")" << std::endl;
        }
#ifdef USE_ANISO_TIMERS
        m_compute_dual_intersection_timer += (clock()-start+0.) / ((double)CLOCKS_PER_SEC);
#endif
        return ret_val;
      }

      bool compute_steiner_dual_intersection(Point_3 &p,
                                             const TPoint_3 &tcandidate_1, // first point within the ball
                                             const TPoint_3 &tcandidate_2, // second point within the ball
                                             const Facet &facet,
                                             const FT circumradius,
                                             const bool verbose = true) const
      {
        CGAL_PROFILER("[compute_steiner_dual_intersection]");

        Point_3 candidate_1 = m_metric.inverse_transform(tcandidate_1);
        Point_3 candidate_2 = m_metric.inverse_transform(tcandidate_2);
        Point_3 ccf;
        compute_dual_intersection(facet, ccf);

        if(Base::dimension() == 2) //need to check if correct todo
        {
          TPoint_3 tfacetp = compute_circumcenter(facet);
          Point_3 facetp = m_metric.inverse_transform(tfacetp);
          Triangle tr = m_metric.inverse_transform(Base::triangle(facet));
          Vector_3 n = tr.supporting_plane().orthogonal_vector();
          return (constrain_ray_intersection(facetp, facetp - n).assign(p)
                 || constrain_ray_intersection(facetp, facetp + n).assign(p));
        }
        else
        {
          if(candidate_1 == candidate_2) // just in case
            return false;

          if(constrain_segment_intersection(candidate_1, candidate_2).assign(p))
          {
            //check if not too far away from the intersection of the dual & surface (should never happen)
            TPoint_3 tccf = m_metric.transform(ccf);
            TPoint_3 tp = m_metric.transform(p);
            if(std::sqrt(CGAL::squared_distance(tccf, tp)) > circumradius)
            {
              std::cout << "intersection point outside of picking ball in compute_steiner_dual_intersection" << std::endl;
              std::cout << "that should not happen !" << std::endl;
              std::cout << "candidates : " << std::endl;
              std::cout << tcandidate_1 << std::endl << tcandidate_2 << std::endl;
              std::cout << "tccf : " << tccf << std::endl;
              std::cout << "p : " << p << std::endl;
              std::cout << "dist : " << std::sqrt(CGAL::squared_distance(tccf, tp)) << " and circumradius : " << circumradius << std::endl;
              return false;
            }
            else
              return true;
          }
          else
          {
//            double eval1 = candidate_1.x()*candidate_1.x()/100. + candidate_1.y()*candidate_1.y() + candidate_1.z()*candidate_1.z() - 1.0;
//            double eval2 = candidate_2.x()*candidate_2.x()/100. + candidate_2.y()*candidate_2.y() + candidate_2.z()*candidate_2.z() - 1.0;
//            double evalc = ccf.x()*ccf.x()/100. + ccf.y()*ccf.y() + ccf.z()*ccf.z() - 1.0;
//            std::cout << "no intersection..." << std::endl;
//            std::cout << "candidate_1 : " << candidate_1 << " " << eval1 << std::endl;
//            std::cout << "candidate_2 : " << candidate_2 << " " << eval2 << std::endl;
//            std::cout << "ccf : " << ccf << " " << evalc << std::endl;
//            std::cout << "circum : " << circumradius << std::endl;
            return false;
          }
        }
      }

      void facet_indices(const Facet& f) const
      {
        std::cout << "Facet[";
        Vertex_handle v1 = f.first->vertex((f.second+1)%4);
        Vertex_handle v2 = f.first->vertex((f.second+2)%4);
        Vertex_handle v3 = f.first->vertex((f.second+3)%4);
        std::cout << v1->info() << " " << v2->info() << " " << v3->info();
        std::cout << "]";
      }

      bool debug_steiner_point(const Point_3& steiner_point,
                               const Facet& f,
                               const bool more = false) const
      {
        //std::cout << "<-debug_steiner_point\n";
        bool bug = false;
        Point_3 p;
        compute_exact_dual_intersection(f,p, false, true, more);

        const TPoint_3& pf = f.first->vertex((f.second+1)&3)->point();
        const TPoint_3& center = m_metric.transform(p);
        const TPoint_3& steiner = m_metric.transform(steiner_point);

        Sphere s(center, CGAL::squared_distance(center, pf));
        if(s.has_on_unbounded_side(steiner))
        {
          facet_indices(f);
          std::cout << "\n\tSteiner point ("
            <<steiner_point<< ") outside exact surface Delaunay ball\n";
          std::cout << "\tp is " << p << std::endl;
          bug = true;
        }

        //Point_3 p2;
        //compute_dual_intersection(f,p2, false,true,more);
        //const TPoint_3& center2 = m_metric.transform(p);
        //Sphere s2(center2, CGAL::squared_distance(center2, pf));
        //if(s2.has_on_unbounded_side(steiner))
        //{
        //  facet_indices(f);
        //  std::cout << "\n\tSteiner point ("
        //    <<steiner_point<< ") outside inexact surface Delaunay ball\n";
        //  bug = true;
        //}
        //std::cout << "->\n";        
        return bug;
      }

      //not used yet
      //bool is_facet_encroached(const TPoint_3 &tp, const Facet &facet)
      //{
      //  TPoint_3 tc = m_metric.transform(compute_exact_dual_intersection(facet));
      //  TPoint_3 tq = facet.first->vertex((facet.second + 1) % 4)->point();
      //  Traits::Compute_squared_distance_3 o = m_traits->compute_squared_distance_3_object();
      //  return o(tc, tp) < o(tc, tq);
      //}

      ////not used yet
      //bool is_encroached(const Point_3 &p, Facet &facet)
      //{
      //  TPoint_3 tp = m_metric.transform(p);

      //  Facet_set_iterator fi = begin_restricted_facets();
      //  Facet_set_iterator fend = end_restricted_facets();
      //  for (; fi != fend; fi++)
      //  {
      //    if (is_facet_encroached(tp, *fi))
      //    {
      //      facet = *fi;
      //      return true;
      //    }
      //  }
      //  return false;
      //}

      void gl_draw_center() const
      {
        if(this->is_surface_star())
        {
          ::glPointSize(5.);
          ::glColor3f(0., 0., 0.);
        }
        else
        {
          ::glPointSize(3.);
          ::glColor3f(0.5f, 0.5f, 0.5f);
        }
        ::glBegin(GL_POINTS);
        ::glVertex3d(m_center_point.x(), m_center_point.y(), m_center_point.z());
        ::glEnd();
      }

      void gl_draw(const typename K::Plane_3& plane,
                   const bool draw_edges = true) const
      {
        if(! this->is_surface_star())
          return;
        gl_draw_center();
        Facet_set_iterator fit = begin_restricted_facets();
        Facet_set_iterator fend = end_restricted_facets();
        for(; fit != fend; fit++)
        {
          const Cell_handle& cell = (*fit).first;
          const int& i = (*fit).second;

          const Point_3& pa = m_metric.inverse_transform(cell->vertex((i+1)&3)->point());
          const Point_3& pb = m_metric.inverse_transform(cell->vertex((i+2)&3)->point());
          const Point_3& pc = m_metric.inverse_transform(cell->vertex((i+3)&3)->point());
          if(is_above_plane(plane, pa, pb, pc))
            gl_draw_triangle<K>(pa, pb, pc, EDGES_AND_FACES, 205., 175., 149);
        }
      }

      void gl_draw_cell(const typename K::Plane_3& plane) const
      {
        gl_draw_center();

        Cell_handle_handle ci = begin_finite_star_cells();
        Cell_handle_handle cend = end_finite_star_cells();
        for (; ci != cend; ci++)
        {
          std::vector<Point_3> cell_points(4);
          cell_points[0] = m_metric.inverse_transform((*ci)->vertex(0)->point());
          cell_points[1] = m_metric.inverse_transform((*ci)->vertex(1)->point());
          cell_points[2] = m_metric.inverse_transform((*ci)->vertex(2)->point());
          cell_points[3] = m_metric.inverse_transform((*ci)->vertex(3)->point());
          if(is_above_plane(plane, cell_points[0], cell_points[1], cell_points[2], cell_points[3]))
          {
            for(int i = 0; i < 4; i++)
            {
              gl_draw_triangle<K>(cell_points[i%4],
                                  cell_points[(i+1)%4],
                                  cell_points[(i+2)%4],EDGES_ONLY, 255, 117, 44);
            }
          }
        }
      }

      void gl_draw_metric(const typename K::Plane_3& plane,
                          double coeff) const
      {
        if(!this->is_surface_star())
          return;
        if(!is_above_plane(plane, this->center_point()))
          return;

        Point_3 p = this->center_point();
        Vector_3 vec = CGAL::NULL_VECTOR;
        double val = 0.;

        val = 1/m_metric.get_min_eigenvalue();
        m_metric.get_min_eigenvector(vec);
        ::glColor3f(0.,0.,1.f);
        ::gl_draw_arrow<K>(p, p+val*coeff*vec);

        val = 1/m_metric.get_max_eigenvalue();
        m_metric.get_max_eigenvector(vec);
        ::glColor3f(1.f,0.,0.);
        ::gl_draw_arrow<K>(p, p+val*coeff*vec);

        val = 1/m_metric.get_third_eigenvalue();
        m_metric.get_third_eigenvector(vec);
        ::glColor3f(0.,1.f,0.);
        ::gl_draw_arrow<K>(p, p+val*coeff*vec);
      }

      bool is_above_plane(const typename K::Plane_3& plane,
                          const typename K::Point_3& p) const
      {
        using CGAL::ON_NEGATIVE_SIDE;
        return (plane.oriented_side(p) == ON_NEGATIVE_SIDE);
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

      bool is_above_plane(const typename K::Plane_3& plane,
                          const typename K::Point_3& pa,
                          const typename K::Point_3& pb,
                          const typename K::Point_3& pc,
                          const typename K::Point_3& pd) const
      {
        typedef typename K::Oriented_side Side;
        using CGAL::ON_NEGATIVE_SIDE;
        const Side sd = plane.oriented_side(pd);
        return (sd == ON_NEGATIVE_SIDE && is_facet_above_plane<K>(plane,pa,pb,pc));
      }

      //void gl_draw_all(const typename K::Plane_3& plane) const
      //{
      //  gl_draw_center();
      //  if(! is_surface_star())
      //    return;
      //  typename Base::Finite_facets_iterator fit2 = this->finite_facets_begin();
      //  typename Base::Finite_facets_iterator fend2 = this->finite_facets_end();
      //  for(; fit2 != fend2; fit2++)
      //  {
      //    Facet f = *fit2;
      //    const Point_3& pa = m_metric.inverse_transform(f.first->vertex((f.second+1)&3)->point());
      //    const Point_3& pb = m_metric.inverse_transform(f.first->vertex((f.second+2)&3)->point());
      //    const Point_3& pc = m_metric.inverse_transform(f.first->vertex((f.second+3)&3)->point());
      //    if(is_above_plane(plane, pa, pb, pc))
      //    {
      //      if(is_in_star(f) && is_restricted(f))
      //        gl_draw_triangle<Kernel>(pa, pb, pc, EDGES_AND_FACES, 205., 175., 149);
      //      else
      //        gl_draw_triangle<Kernel>(pa, pb, pc, EDGES_ONLY/*, 139., 137., 137.*/);
      //    }
      //  }
      //}

      void gl_draw_dual(const typename K::Plane_3& plane) const
      {
        GLboolean was = (::glIsEnabled(GL_LIGHTING));
        if(was)
          ::glDisable(GL_LIGHTING);

        typename Base::Finite_facets_iterator fit = this->finite_facets_begin();
        typename Base::Finite_facets_iterator fend = this->finite_facets_end();
        for(; fit != fend; fit++)
        {
          Facet f = *fit;
          if(!is_in_star(f)) continue;
          const Point_3& pa = m_metric.inverse_transform(f.first->vertex((f.second+1)&3)->point());
          const Point_3& pb = m_metric.inverse_transform(f.first->vertex((f.second+2)&3)->point());
          const Point_3& pc = m_metric.inverse_transform(f.first->vertex((f.second+3)&3)->point());
          if(!is_above_plane(plane, pa, pb, pc))
            continue;

          if(is_restricted(f)) { ::glColor3f(0.,0.,1.f); ::glLineWidth(2.);}
          else                 { ::glColor3f(0.,0.,0.);  ::glLineWidth(1.);}

          typename K::Object_3 o = dual(f);
          typename K::Line_3 l;
          typename K::Ray_3 r;
          typename K::Segment_3 s;
          if(CGAL::assign(s, o))
            gl_draw_segment<K>(s.source(), s.target());
          else if(CGAL::assign(r, o))
            gl_draw_segment<K>(r.source(), r.source() + r.to_vector()
                * (m_pConstrain->get_bounding_radius() * 10.0 / std::sqrt(r.to_vector() * r.to_vector())));
          else if(CGAL::assign(l, o))
            std::cout << "gl_draw_dual : line dual\n";

          //experimental stuff
          Point_3 p;
          if(compute_dual_intersection(f, p))
          {
            ::glColor3f(1.f,0.,0.);
            ::glPointSize(5.);
            ::glBegin(GL_POINTS);
            ::glVertex3d(p.x(),p.y(),p.z());
            ::glEnd();

            //Point_3 centroid = CGAL::centroid(pa, pb, pc);
            //gl_draw_segment<K>(p, centroid);
          }
        }
        if(was)
          ::glEnable(GL_LIGHTING);
      }

      bool is_in_star(const Facet& f) const
      {
        for(unsigned int i = 0; i < 4; i++)
        {
          Vertex_handle v = f.first->vertex((f.second+i)%4);
          if(v == m_center)
          {
            int index = f.first->index(m_center);
            if(index != f.second)
              return true;
          }
        }
        return false;
      }

      void gl_draw_surface_delaunay_balls(const typename K::Plane_3& plane) const
      {
        if(!is_above_plane(plane, this->center_point()))
          return;

        Facet_set_iterator fit = begin_restricted_facets();
        Facet_set_iterator fend = end_restricted_facets();
        for(; fit != fend; fit++)
        {
          Facet f = *fit;
          Point_3 ce;
          if(!is_above_plane(plane, m_metric.inverse_transform(f.first->vertex((f.second+1) % 4)->point()),
                                    m_metric.inverse_transform(f.first->vertex((f.second+2) % 4)->point()),
                                    m_metric.inverse_transform(f.first->vertex((f.second+3) % 4)->point())) )
            continue;
#ifdef ANISO_USE_EXACT
          this->compute_exact_dual_intersection(f, ce);
#else
          this->compute_dual_intersection(f, ce);
#endif
          TPoint_3 tce = m_metric.transform(ce);
          TPoint_3 tp2 = f.first->vertex((f.second+1) % 4)->point();
          //Point_3 p2 = m_metric.inverse_transform(tp2);

          FT sqr = CGAL::squared_distance(tce, tp2);
          //gl_draw_sphere<K>(typename K::Sphere_3(ce, sqr));

            //ellipsoid's a, b & c
          FT a = std::sqrt(sqr)/m_metric.get_max_eigenvalue();
          FT b = std::sqrt(sqr)/m_metric.get_min_eigenvalue();
          FT c = std::sqrt(sqr)/m_metric.get_third_eigenvalue();

            //rotation to align on the eigenvectors
          Vector_3 v1,v2,vn;
          m_metric.get_max_eigenvector(v1);
          m_metric.get_min_eigenvector(v2);
          m_metric.get_third_eigenvector(vn);

          ::GLdouble rot_mat[16];
          rot_mat[0] = v1.x(); rot_mat[4] = v2.x(); rot_mat[8] = vn.x();  rot_mat[12] = ce.x();
          rot_mat[1] = v1.y(); rot_mat[5] = v2.y(); rot_mat[9] = vn.y();  rot_mat[13] = ce.y();
          rot_mat[2] = v1.z(); rot_mat[6] = v2.z(); rot_mat[10] = vn.z(); rot_mat[14] = ce.z();
          rot_mat[3] = 0.; rot_mat[7] = 0.; rot_mat[11] = 0.; rot_mat[15] = 1.;

          ::glMatrixMode (GL_MODELVIEW);
          ::glPushMatrix();
          ::glMultMatrixd(rot_mat);
          if(ellipsoid_color == 1)
            gl_draw_ellipsoid<K>(CGAL::ORIGIN, 10, 10, a, b, c, 232, 67, 183);
          else if(ellipsoid_color == 2)
            gl_draw_ellipsoid<K>(CGAL::ORIGIN, 10, 10, a, b, c, 20, 20, 233);
          else
          {
            //gl_draw_ellipsoid<K>(CGAL::ORIGIN, 10, 10, a, b, c, 245, 230, 93);
          }
          ::glPopMatrix();

            //a b & c visu
          if(ellipsoid_color == 1 || ellipsoid_color == 2)
          {
            ::glColor3f(1.f,0.,0.);
            gl_draw_segment<K>(ce, ce+a*v1);
            ::glColor3f(0.,0.,1.f);
            gl_draw_segment<K>(ce, ce+b*v2);
            ::glColor3f(0.,1.f,0.);
            gl_draw_segment<K>(ce, ce+c*vn);

              //center
            ::glColor3d(33,224,237);
            ::glPointSize(10.);
            ::glBegin(GL_POINTS);
            ::glVertex3d(ce.x(),ce.y(),ce.z());
            ::glEnd();
          }
        }
      }

      void invalidate()
      {
        this->clear();
        m_center_point = CGAL::ORIGIN;
        m_center = Vertex_handle();
        invalidate_cache();
      }

      bool is_valid(const bool verbose = false) const
      {
        return (m_center != Vertex_handle()) && Base::is_valid(verbose);
      }

      void reset(const Point_3 &centerpoint,
                 const int &index,
                 const Metric &metric_,
                 const bool is_surface_star = true)
      {
        this->clear();
        m_center_point = centerpoint;
        m_metric = metric_;
        m_center = Base::insert(m_metric.transform(centerpoint));
        m_center->info() = index;
        m_is_surface_star = is_surface_star;
        this->infinite_vertex()->info() = index_of_infinite_vertex;

        invalidate_cache();
      }

    public:
      Stretched_Delaunay_3(const Criteria* criteria_,
                           const Constrain_surface* pconstrain_surface,
                           const bool is_surface_star = true) :
        Base(*(m_traits = new Traits())),
        m_center_point(CGAL::ORIGIN),
        m_metric(),
        m_pConstrain(pconstrain_surface),
        m_criteria(new Stretched_criteria(*m_traits, criteria_)),
        m_is_surface_star(is_surface_star),
        m_is_topological_disk(false),
        m_is_valid_topo_disk(false),
        m_is_valid_bbox(false),
        is_cache_dirty(true),
        restricted_facets_cache(),
        neighboring_cells_cache(),
        neighboring_finite_cells_cache()
      {
        ellipsoid_color = 0;
        m_center = Vertex_handle();
        m_bbox = m_pConstrain->get_bbox(); // in M_euclidean
        this->infinite_vertex()->info() = index_of_infinite_vertex;
      }

      Stretched_Delaunay_3(const Criteria* criteria_,
                           const Point_3 &centerpoint,
                           const int &index,
                           const Metric &metric_,
                           const Constrain_surface* pconstrain_surface,
                           const bool is_surface_star = true) :
        Base(*(m_traits = new Traits())),
        m_center_point(centerpoint),
        m_metric(metric_),
        m_pConstrain(pconstrain_surface),
        m_criteria(new Stretched_criteria(*m_traits, criteria_)),
        m_is_surface_star(is_surface_star),
        m_is_topological_disk(false),
        m_is_valid_topo_disk(false),
        m_is_valid_bbox(false),
        is_cache_dirty(true),
        restricted_facets_cache(),
        neighboring_cells_cache(),
        neighboring_finite_cells_cache()
      {
        ellipsoid_color = 0;
        m_center = Base::insert(m_metric.transform(centerpoint));
        m_center->info() = index;
        m_bbox = m_pConstrain->get_bbox(); // in M_euclidean
        this->infinite_vertex()->info() = index_of_infinite_vertex;
      }

      ~Stretched_Delaunay_3()
      {
        delete m_traits;
        delete m_criteria;
        //delete m_pConstrain;
      }
    };

#ifdef USE_ANISO_TIMERS
    template<typename K, typename KExact>
    double Stretched_Delaunay_3<K, KExact>::m_compute_dual_intersection_timer = 0.;
#endif
  }
}

#endif //CGAL_ANISOTROPIC_MESH_3_STRETCHED_DELAUNAY_3

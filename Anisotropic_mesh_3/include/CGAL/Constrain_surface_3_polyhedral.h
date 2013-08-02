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

#ifndef CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_POLYHEDRAL_H
#define CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_POLYHEDRAL_H

#include <CGAL/Cartesian.h>
#include <CGAL/Monge_via_jet_fitting.h>

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <iterator>

#include "../../demo/Polyhedron_type.h"

#include <CGAL/Constrain_surface_3_ex.h>
#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <Eigen/Dense>
#include <CGAL/Metric.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>

// for Mesh_3
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>

#include <CGAL/gl_draw/drawing_helper.h>
#include <CGAL/helpers/tiling_helper.h>

namespace CGAL
{
  namespace Anisotropic_mesh_3
  {

// a refined facet with a tag
template <class Refs, class T>
class Enriched_facet : public CGAL::HalfedgeDS_face_base<Refs, T>
{
  // tag
  std::size_t m_tag;

public:
  // life cycle
  // no constructors to repeat, since only
  // default constructor mandatory
  Enriched_facet()
    : m_tag(0)  {}
  // tag
  const std::size_t& tag() const { return m_tag; }
  std::size_t& tag() { return m_tag; }
};

// a refined halfedge with a tag
template <class Refs, class Tprev, class Tvertex, class Tface>
class Enriched_halfedge : public CGAL::HalfedgeDS_halfedge_base<Refs,Tprev,Tvertex,Tface>
{
private:
  // general purpose tag
  std::size_t m_tag;

public:
  // life cycle
  Enriched_halfedge()
    : m_tag(0)  {}
  // tag
  const std::size_t& tag() const { return m_tag;  }
  std::size_t& tag()             { return m_tag;  }
};

//enum Point_geo_classification { PARABOLIC = 0, ELLIPTIC = 1, HYPERBOLIC = -1};

// a refined vertex with a tag
template <class Refs, class T, class P>
class Enriched_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P>
{
  // tag
  std::size_t m_tag;
  int m_classification; // point classification (elliptic, hyperbolic, etc.)
  Vector m_normal;
  FT m_aniso_ratio;

#if 1//def TILING_DRAW
  std::vector<Point> m_tile_points;
#endif

public:
  // life cycle
  Enriched_vertex()  {}
  // repeat mandatory constructors
  Enriched_vertex(const P& pt)
    : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt),
      m_tag(0), m_classification(0), m_normal(CGAL::NULL_VECTOR), m_aniso_ratio(0) {}

  Enriched_vertex(const P& pt, const int& classification, const Vector& n)
    : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt),
      m_tag(0), m_classification(classification), m_normal(n), m_aniso_ratio(0) {}

  // tag
  std::size_t& tag() {  return m_tag; }
  const std::size_t& tag() const {  return m_tag; }
  // classification
  int& classification() {  return m_classification; }
  const int& classification() const {  return m_classification; }
  // normal
  Vector& normal() {  return m_normal; }
  const Vector& normal() const {  return m_normal; }
  // aniso_ratio
  FT& aniso_ratio() {  return m_aniso_ratio; }
  const FT& aniso_ratio() const {  return m_aniso_ratio; }
  // tile_points
  std::vector<Point>& tile_points() {  return m_tile_points; }
  const std::vector<Point>& tile_points() const {  return m_tile_points; }
}; 

struct Enriched_items : public CGAL::Polyhedron_items_3
{
  // wrap vertex
  template<class Refs, class Traits>
  struct Vertex_wrapper
  {
    typedef typename Traits::Point_3 Point;
    typedef Enriched_vertex<Refs, CGAL::Tag_true, Point> Vertex;
  };

  // wrap face
  template<class Refs, class Traits>
  struct Face_wrapper
  {
    typedef Enriched_facet<Refs, CGAL::Tag_true> Face;
  };

  // wrap halfedge
  template<class Refs, class Traits>
  struct Halfedge_wrapper
  {
    typedef Enriched_halfedge<Refs,
                              CGAL::Tag_true,
                              CGAL::Tag_true,
                              CGAL::Tag_true> Halfedge;
  };
};

template<typename K, typename Polyhedron = CGAL::Polyhedron_3<K, Enriched_items> >
class Constrain_surface_3_polyhedral : 
    public Constrain_surface_3_ex<K, typename Constrain_surface_3<K>::Point_container> 
    {
    public:
      typedef Constrain_surface_3_ex<K, typename Constrain_surface_3<K>::Point_container> Base;

      typedef typename Base::Object_3      Object_3;
      typedef typename Base::FT                          FT;
      typedef typename Base::Vector_3                    Vector_3;
      typedef typename Base::Point_3                     Point_3;
      typedef typename Base::Segment_3                   Segment_3;
      typedef typename Base::Triangle_3                  Triangle_3;
      typedef typename Base::Plane_3                     Plane_3;
      typedef typename Base::Oriented_side               Oriented_side;
      typedef typename Base::Pointset                    Pointset;

      typedef typename CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;

      typedef typename CGAL::AABB_polyhedron_triangle_primitive<K, Polyhedron> Primitive;
      typedef typename CGAL::AABB_traits<K, Primitive> Traits;
      typedef typename CGAL::AABB_tree<Traits> Tree;
      typedef typename Tree::Object_and_primitive_id Object_and_primitive_id;
      typedef typename Tree::Point_and_primitive_id Point_and_primitive_id;

      typedef typename CGAL::Monge_via_jet_fitting<K> My_Monge_via_jet_fitting;
      typedef typename My_Monge_via_jet_fitting::Monge_form     My_Monge_form;

      typedef typename Constrain_surface_3<K>::Point_container Point_container;

      typedef typename Polyhedron::Vertex_handle    Vertex_handle;
      typedef typename Polyhedron::Facet_handle     Facet_handle;
      typedef typename Polyhedron::Halfedge_handle  Halfedge_handle;
      typedef typename Polyhedron::Halfedge_around_vertex_circulator HV_circulator;

      // for Mesh_3
      typedef typename CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
      typedef typename CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
      typedef typename CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

      // exact
      typedef typename CGAL::Exact_predicates_exact_constructions_kernel KExact;
      typedef typename KExact::FT                           Exact_FT;
      typedef typename KExact::Point_3                      Exact_Point_3;
      typedef typename KExact::Vector_3                     Exact_Vector_3;
      typedef typename KExact::Segment_3                    Exact_Segment_3;
      typedef typename KExact::Triangle_3                   Exact_Triangle_3;
      typedef typename KExact::Plane_3                      Exact_Plane_3;

      typedef CGAL::Cartesian_converter<K, KExact> To_exact;
      typedef CGAL::Cartesian_converter<KExact, K> Back_from_exact;
      To_exact to_exact;
      Back_from_exact back_from_exact;

    public:
      Mesh_domain* domain;
      Polyhedron   m_polyhedron;
      FT           bounding_radius;
      Tree         *tree;
      FT           minx, miny, minz, maxx, maxy, maxz;
      FT           nearest_start_try_radius;
      int          point_count;

      std::vector<Vertex_handle> m_vertices;
      std::vector<typename Eigen::Matrix3d> m_metrics;

#ifdef CGAL_DEBUG_OUTPUT_VECTOR_FIELD
      std::ofstream* min_vector_field_os;
      std::ofstream* max_vector_field_os;
      std::ofstream* normals_field_os;
#endif

    protected:
      mutable C3t3 m_c3t3;

      mutable double m_max_curvature;
      mutable double m_min_curvature;
      mutable bool m_cache_max_curvature;
      mutable bool m_cache_min_curvature;
      
    public:
      virtual FT get_bounding_radius() const
      {
        return bounding_radius;
      }

      virtual std::string name() const
      {
        return std::string("Polyhedron");
      }

      virtual typename CGAL::Bbox_3 get_bbox() const
      {
        //note : compute_bounding_box was run in constructor
        return typename CGAL::Bbox_3(minx, miny, minz, maxx, maxy, maxz);
      }

      virtual Oriented_side side_of_constraint(const Point_3 &p) const 
      {
        if ((bool)domain->is_in_domain_object()(p))
          return CGAL::ON_POSITIVE_SIDE;
        else
          return CGAL::ON_NEGATIVE_SIDE;
      }

    public:
      struct by_distance_to_p0
      {
        Point_3 p0;

        bool operator()(const Point_3 &pi, const Point_3& pj)
        {
          return CGAL::squared_distance(p0, pi) < CGAL::squared_distance(p0, pj);
        }

        by_distance_to_p0(const Point_3 &pini) : p0(pini) { }
      };

      virtual Object_3 intersection(const Point_3 &p0, const Point_3 &p1) const
      {
        typedef typename Tree::Object_and_primitive_id Object_and_primitive_id;
        std::list<Object_and_primitive_id> intersections;
        Segment_3 seg(p0, p1);
        std::vector<Point_3> intersection_points;

        tree->all_intersections(seg, std::back_inserter(intersections));

        typename std::list<Object_and_primitive_id>::iterator it;
        for(it = intersections.begin(); it != intersections.end(); it++)
        {
          Object_and_primitive_id ob = *it;
          Object_3 object = ob.first;
          Point_3 p;
          if(object.assign(p))
            intersection_points.push_back(p);
        }

        if(intersection_points.size()%2 == 0)
          return Object_3();

        std::sort(intersection_points.begin(), intersection_points.end(), by_distance_to_p0(p0));

#ifdef ANISO_DEBUG_INTERSECTION
        if(intersection_points.size() > 1)
        {
          std::cout << "intersection poly : " << std::endl << p0 << std::endl << p1 << std::endl;
          std::cout << intersection_points.size() << " intersections" << std::endl;
          for(typename std::vector<Point_3>::iterator vit = intersection_points.begin(); vit != intersection_points.end(); ++vit)
            std::cout << *vit << " || " << CGAL::squared_distance(*vit, p0) << std::endl;
        }
#endif
        return make_object(intersection_points.front());
      }

      FT compute_sq_approximation(const Point_3& p) const
      {
        return tree->squared_distance(p);
      }

      virtual void tensor_frame_on_input_point(Vertex_handle v, 
        Vector_3 &v0/*n*/, Vector_3 &v1/*vmax*/, Vector_3 &v2/*vmin*/,//unit vectors
        FT& c1/*cmax*/, FT& c2/*cmin*/,
        const FT& epsilon) 
      {
        std::vector<Vertex_handle> points;
        find_nearest_vertices(v, points, 36/*, max_sqd=0 by default*/);
//          nearest_start_try_radius*nearest_start_try_radius);

        std::vector<Point_3> points_geo;
        for(std::size_t i=0; i<points.size(); ++i)
          points_geo.push_back(points[i]->point());

        My_Monge_form monge_form;
        My_Monge_via_jet_fitting monge_fit;
        monge_form = monge_fit(points_geo.begin(), points_geo.end(), 4, 4);

        FT maxc = monge_form.principal_curvatures(0);
        FT minc = monge_form.principal_curvatures(1);
        v0 = monge_form.normal_direction();

        if(minc >= 0.)
        {
          c1 = (std::max)(epsilon, fabs(maxc));
          c2 = (std::max)(epsilon, fabs(minc));
          v1 = monge_form.maximal_principal_direction();
          v2 = monge_form.minimal_principal_direction();
        }
        else if(maxc <= 0.)
        {
          c2 = (std::max)(epsilon, fabs(maxc));
          c1 = (std::max)(epsilon, fabs(minc));
          v2 = monge_form.maximal_principal_direction();
          v1 = monge_form.minimal_principal_direction();
        }
        else //minc < 0 && maxc > 0
        {
          FT abs_min = fabs(minc);
          if(abs_min < maxc)
          {
            c1 = (std::max)(epsilon, fabs(maxc));
            c2 = (std::max)(epsilon, fabs(minc));
            v1 = monge_form.maximal_principal_direction();
            v2 = monge_form.minimal_principal_direction();
          }
          else
          {
            c2 = (std::max)(epsilon, fabs(maxc));
            c1 = (std::max)(epsilon, fabs(minc));
            v2 = monge_form.maximal_principal_direction();
            v1 = monge_form.minimal_principal_direction();
          }
        }
      }

      void tensor_frame(const Point_3 &p,
                        Vector_3 &v0,     //unit normal
                        Vector_3 &v1,     //unit eigenvector
                        Vector_3 &v2,     //unit eigenvector
                        double& e0,       //eigenvalue corresponding to v0
                        double& e1,       //eigenvalue corresponding to v1
                        double& e2) const //eigenvalue corresponding to v2
      {
        //choose r_b for tensor blending
        Facet_handle facet = find_nearest_facet(p);
        std::vector<Vertex_handle> neighbors;
        find_neighbors(facet, neighbors);
        
        FT sum_dist = 0.;
        std::vector<FT> sq_distances;
        for(std::size_t i = 0; i < neighbors.size(); ++i)
        {
          sq_distances.push_back(CGAL::squared_distance(p, neighbors[i]->point()));
          sum_dist += std::sqrt(sq_distances[i]);
        }
        FT rb = sum_dist / ((int)neighbors.size());
        FT wpp_inv = 1./(rb * rb);

        // blend
        Eigen::Matrix3d m = Eigen::Matrix3d::Zero();
        FT wsum = 0.;
        for (std::size_t i = 0; i < neighbors.size(); ++i)
        { 
          FT w = std::exp(-wpp_inv * sq_distances[i]);
          wsum += w;
          std::size_t index = neighbors[i]->tag();
          for(int j = 0; j < 3; j++)
            for(int k = 0; k < 3; k++)
              m(j,k) = m(j,k) + w * m_metrics[index](j,k);
        }

        FT wsum_inv = 1. / wsum;
        for(int j = 0; j < 3; j++)
          for(int k = 0; k < 3; k++)
            m(j,k) = m(j,k) * wsum_inv;
        
        Eigen::EigenSolver<Eigen::Matrix3d> es(m, true/*compute eigenvectors and values*/);
        const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvalueType& vals = es.eigenvalues();
        const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvectorsType& vecs = es.eigenvectors();

        e0 = (std::real(vals[0]) >= 0.) ? std::real(vals[0]) : 0.;
        e1 = (std::real(vals[1]) >= 0.) ? std::real(vals[1]) : 0.;
        e2 = (std::real(vals[2]) >= 0.) ? std::real(vals[2]) : 0.;
                
        v0 = get_eigenvector(vecs.col(0));
        v1 = get_eigenvector(vecs.col(1));
        v2 = get_eigenvector(vecs.col(2));
      }

      void find_neighbors(Facet_handle facet, 
                          std::vector<Vertex_handle>& neighbors) const
      {
        Halfedge_handle he1 = facet->halfedge();
        Halfedge_handle he2 = he1->next();
        Halfedge_handle he3 = he2->next();
        neighbors.push_back(he1->vertex());
        neighbors.push_back(he1->opposite()->vertex());
        neighbors.push_back(he2->vertex());
        neighbors.push_back(he2->opposite()->vertex());
        neighbors.push_back(he3->vertex());
        neighbors.push_back(he3->opposite()->vertex());
      }

      bool are_equal(const Vector_3& v1,        //unit vector
                     const Vector_3& v2) const  //unit vector
      {
        double cosine = std::abs(v1*v2);
        if(std::abs(cosine - 1.) < 0.2)
          return true;
        std::cout << "[" << cosine << "]";
        return false;
      }
      Vector_3 normalize(const Vector_3& v) const
      {
        return std::sqrt(1./(v*v)) * v;
      }

      Vector_3 get_vector(const Eigen::Vector3d& v) const
      {
       return Vector_3(v[0], v[1], v[2]);
      }
      Vector_3 get_eigenvector(const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvectorsType::ConstColXpr& v) const
      {
        return Vector_3(std::real(v[0]), std::real(v[1]), std::real(v[2]));
      }

      virtual Point_container initial_points(const int nb) const 
      {
        typedef typename std::vector<std::pair<
          typename Mesh_domain::Point_3, 
          typename Mesh_domain::Index> > PointIndex_container;
        PointIndex_container points_indices;
        std::back_insert_iterator<PointIndex_container> inserter(points_indices);
        domain->construct_initial_points_object()(inserter, nb);

        Point_container points;
        std::size_t i;
        std::size_t N = points_indices.size();
        for(i = 0; i < N; i++)
          points.push_back(points_indices[i].first);
        return points;
      }

      inline Vertex_handle find_nearest_vertex(const Point_3 &p) const
      {
        Point_and_primitive_id pp = tree->closest_point_and_primitive(p);
        typename Polyhedron::Facet_handle facet = pp.second;
        typename Polyhedron::Facet::Halfedge_handle ci = facet->halfedge();
        return ci->vertex();
      }

      inline Facet_handle find_nearest_facet(const Point_3 &p) const
      {
        Point_and_primitive_id pp = tree->closest_point_and_primitive(p);
        Facet_handle facet = pp.second;
        return facet;
      }

      void find_nearest_vertices(Vertex_handle v,
                                 std::vector<Vertex_handle>& points,
                                 int count,
                                 const double max_sqd = 0.)
      {
        std::set<Vertex_handle> visited;
        std::map<FT, Vertex_handle/*sqd to center : closest come first*/> ring;

        ring.insert(std::make_pair<FT,Vertex_handle>(0.,v));       
        visited.insert(v);

        while(points.size() < count && !ring.empty())
        {
          typename std::map<FT, Vertex_handle>::iterator it = ring.begin();
          Vertex_handle vi = (*it).second;
          ring.erase(it);

          points.push_back(vi);
          find_ring_vertices(v->point(), vi, visited, ring, max_sqd);
        }
      }

      void find_ring_vertices(const Point_3& center,
                              Vertex_handle v,
                              std::set<Vertex_handle>& visited,
                              std::map<FT,Vertex_handle>& ring,
                              const double max_sqd = 0.)/**/
      {
        HV_circulator h = v->vertex_begin();
        HV_circulator hend = h;
        do
        {
          Vertex_handle vv = h->opposite()->vertex();
          if(visited.find(vv) == visited.end()) //vertex not already visited
          {
            visited.insert(vv);
            FT sqd = CGAL::squared_distance(center, vv->point());
            if(max_sqd == 0. || sqd < max_sqd)
              ring.insert(std::make_pair<FT,Vertex_handle>(sqd,vv));
          }
          h++;
        }
        while(h != hend);
      }

      inline void find_nearest_vertices_kanle(Vertex_handle v, 
                                        std::vector<Vertex_handle>& points,
                                        int count) 
      {
        const Point_3& p = v->point();
        FT dist = nearest_start_try_radius;

        int iter = 0, max_iter = 100;
        double max_search_dist = 1e5;
        bool found = false;

        while (!found && iter < max_iter)
        {
          ++iter;
          points.clear();
          find_nearest_vertices(p, points, dist);
          if ((int)points.size() > count * 2) {
            dist = dist * 0.8;
          } else if ((int)points.size() < count / 2) {
            dist = dist * 1.2;
          } else {
            nearest_start_try_radius = nearest_start_try_radius * 0.6 + dist * 0.4;
            found = true;
          }
        }
        if(!found && iter == max_iter){
          while( (int)points.size() < count / 2 ){
            dist = dist * 1.01;
            points.clear();
            find_nearest_vertices(p, points, dist);
            if(dist > max_search_dist){
                std::cerr << "Not enough points close enough to be found in find_nearest_vertices" << std::endl;
                break;
            }
          }
          nearest_start_try_radius = nearest_start_try_radius * 0.6 + dist * 0.4;
          if((int)points.size() > count * 2){
            int diff = (int)points.size() - (count * 2);
            for(int i = 0; i<diff; ++i)
                points.pop_back();
          }
        }
      }

      inline void find_nearest_vertices(const Point_3 &p, std::vector<Vertex_handle> &points, const FT max_dist)
      {
        std::list<typename Polyhedron::Vertex_handle> visited;
        std::list<typename Polyhedron::Vertex_handle> queue;
        queue.push_back(find_nearest_vertex(p));

        FT squared_max_dist = max_dist * max_dist;
        while (!queue.empty()) {
          typename Polyhedron::Vertex_handle top = queue.front();
          queue.pop_front();

          typename std::list<typename Polyhedron::Vertex_handle>::iterator i = visited.begin();
          typename std::list<typename Polyhedron::Vertex_handle>::iterator iend = visited.end();
          bool has_visited = false;
          for (; i != iend; i++)
            if ((*i) == top) {
              has_visited = true;
              break;
            }
            if (has_visited)
              continue;

            visited.push_back(top);

            Point_3 q = top->point();
            FT dist = (p.x() - q.x()) * (p.x() - q.x()) + 
              (p.y() - q.y()) * (p.y() - q.y()) + 
              (p.z() - q.z()) * (p.z() - q.z());
            if (dist > squared_max_dist)
              continue;

            points.push_back(top);

            typename Polyhedron::Halfedge_handle end_e = top->halfedge();
            typename Polyhedron::Halfedge_handle e = end_e;
            do {
              e = e->opposite();
              queue.push_back(e->vertex());
              e = e->prev();
            } while (e != end_e);
        }
      }

      void set_domain_and_polyhedron(const char* filename)
      {
        std::cout << "Loading polyhedron..." << std::endl;
        std::ifstream input(filename);
        if(!input)
          std::cout << "\nWarning : file does not exist" << std::endl;
        Polyhedron p;
        input >> p;
        set_domain_and_polyhedron(p);
      }

      void set_domain_and_polyhedron(const Polyhedron& p)
      {
        m_polyhedron = p;
        std::cout << "[" << m_polyhedron.size_of_vertices() << "]" << std::endl;
        std::cout << "\nCreating domain..." << std::endl;
        domain = new Mesh_domain(m_polyhedron);
      }

      void set_aabb_tree()
      {
        std::cout << "\nCreating AABB tree..." << std::endl;
        tree = new Tree(m_polyhedron.facets_begin(), m_polyhedron.facets_end());
        std::cout << "\nAccelerating distance queries..." << std::endl;
        tree->accelerate_distance_queries();
      }

      void compute_bounding_box(/*std::vector<Vertex_handle>& points*/)
      {
        std::cout << "\nComputing bounding box..." << std::endl;

        /*typename Polyhedron::Vertex_iterator vi = m_polyhedron.vertices_begin();
        typename Polyhedron::Vertex_iterator viend = m_polyhedron.vertices_end();
        for (; vi != viend; vi++)
          points.push_back(vi);
        */  
        minx = miny = minz = DBL_MAX;
        maxx = maxy = maxz = DBL_MIN;
        //for (int i = 0; i < (int)points.size(); i++) 
        std::size_t index = 0;
        typename Polyhedron::Vertex_iterator vi = m_polyhedron.vertices_begin();
        typename Polyhedron::Vertex_iterator viend = m_polyhedron.vertices_end();
        for (; vi != viend; ++vi, ++index)
        {
          vi->tag() = index;
          m_vertices.push_back(vi);//[index] = vi;
          Point_3 pi = vi->point();

          FT x = pi.x(), y = pi.y(), z = pi.z();
          if (x < minx) minx = x;   if (x > maxx) maxx = x;
          if (y < miny) miny = y;   if (y > maxy) maxy = y;
          if (z < minz) minz = z;   if (z > maxz) maxz = z;
        }
        bounding_radius = sqrt(
          (maxx - minx) * (maxx - minx) + 
          (maxy - miny) * (maxy - miny) + 
          (maxz - minz) * (maxz - minz)) * 0.5 * 1.1;
        nearest_start_try_radius = 
          sqrt((4 * bounding_radius * bounding_radius) / ((FT)m_vertices.size() * 2.0));
      }

      void get_metrics_from_file(std::ifstream& input,
                                 const FT& epsilon,
                                 const FT& en_factor = 1.0)
      {
        std::cout << "reading metrics from file" << std::endl;
        size_t nv;
        input >> nv;
        if(nv != m_vertices.size())
        {
          std::cerr << "input metric file is not adapted to the polyhedron" << std::endl;
          return;
        }

        FT x,y,z;
        FT e1, e2, en;
        Vector_3 v1, v2, vn;
        Point_3 pi;

        for (std::size_t i = 0; i < m_vertices.size(); ++i)
        {
          input >> x >> y >> z;
          pi = Point_3(x,y,z);
          input >> e1 >> e2;
          input >> x >> y >> z;
          v1 = Vector_3(x,y,z);
          input >> x >> y >> z;
          v2 = Vector_3(x,y,z);
          input >> x >> y >> z;
          vn = Vector_3(x,y,z);

          e1 = std::sqrt(std::abs(e1));
          e2 = std::sqrt(std::abs(e2));
          en = en_factor*(std::max)(e1, e2);

          m_vertices[i]->normal() = vn;

          Metric_base<K, K> M(vn, v1, v2, en, e1, e2, epsilon);
          Eigen::Matrix3d transf = M.get_transformation();
          m_metrics.push_back(transf.transpose()*transf);

#ifdef CGAL_DEBUG_OUTPUT_VECTOR_FIELD
          FT f = get_bounding_radius()*0.02;
          (*normals_field_os) << "2 " << pi << " " << (pi+f*vn) << std::endl;
          if(e1>e2)
          {
            (*max_vector_field_os) << "2 " << pi << " " << (pi+f*v1) << std::endl;
            (*min_vector_field_os) << "2 " << pi << " " << (pi+f*v2) << std::endl;
          }
          else
          {
            (*min_vector_field_os) << "2 " << pi << " " << (pi+f*v1) << std::endl;
            (*max_vector_field_os) << "2 " << pi << " " << (pi+f*v2) << std::endl;
          }
#endif
        }
      }

      void compute_local_metric(const FT& epsilon,
                                const FT& en_factor,
                                const bool smooth_metric = false)
      {
        std::cout << "\nComputing local metric..." << std::endl;
        m_max_curvature = 0.;
        m_min_curvature = DBL_MAX;

        std::ofstream out("poly_curvatures.txt");
        out << m_vertices.size() << std::endl;

        for (std::size_t i = 0; i < m_vertices.size(); ++i)
        {
          if (i % 100 == 0)
            std::cout << ".";

          Vector_3 vn, v1, v2;
          FT e1, e2, en;
          tensor_frame_on_input_point(m_vertices[i], 
            vn/*normal*/, v1/*vmax*/, v2/*vmin*/,
            e1/*cmax*/, e2/*cmin*/, epsilon);

          en = en_factor*(std::max)(e1, e2);
          m_max_curvature = (std::max)(m_max_curvature, en);
          m_min_curvature = (std::min)(m_min_curvature, (std::min)(e1, e2));

          e1 = std::sqrt(e1);
          e2 = std::sqrt(e2);
          en = en_factor*(std::max)(e1, e2);

          const Point_3& pi = m_vertices[i]->point();
//          tensor_frame_on_point(pi, vn, v1, v2, e1, e2, epsilon);

          out << pi.x() << " " << pi.y() << " " << pi.z() << std::endl;
          out << e1 << " " << e2 << "     ";
          out << v1.x() << " " << v1.y() << " " << v1.z() << "     ";
          out << v2.x() << " " << v2.y() << " " << v2.z() << "     ";
          out << vn.x() << " " << vn.y() << " " << vn.z() << std::endl;

#ifdef CGAL_DEBUG_OUTPUT_VECTOR_FIELD
          if(!smooth_metric)
          {
            FT f = get_bounding_radius()*0.02;
            (*max_vector_field_os) << "2 " << pi << " " << (pi+f*v1) << std::endl;
            (*min_vector_field_os) << "2 " << pi << " " << (pi+f*v2) << std::endl;
            (*normals_field_os) << "2 " << pi << " " << (pi+f*vn) << std::endl;
          }
#endif

          std::size_t index = m_vertices[i]->tag();
          if(i != index) 
            std::cerr << "Error1 in indices in compute_local_metric!" << std::endl;

          m_vertices[i]->normal() = vn;

          Metric_base<K, K> m(vn, v1, v2, en, e1, e2, epsilon);
          Eigen::Matrix3d transf = m.get_transformation();
          m_metrics.push_back(transf.transpose()*transf);

          //if(m_metrics.size() != i+1 )
          //  std::cerr << "Error3 in indices in compute_local_metric!" << std::endl;
        }

        m_cache_max_curvature = true;
        m_cache_min_curvature = true;
        std::cout << "done." << std::endl;
      }

      void heat_kernel_smoothing(const FT& epsilon)
      {
        std::cout << "Smoothing the metric field" << std::endl;
        int smooth_step_n = 3;

        while(smooth_step_n--)
        {
          std::cout << "computing new metrics... steps left : " << smooth_step_n << std::endl;

          std::vector<typename Eigen::Matrix3d> temp_m_metrics;
          for (std::size_t i = 0; i < m_vertices.size(); ++i)
          {
            Vertex_handle vi = m_vertices[i];

            std::vector<Vertex_handle> points;
            std::set<Vertex_handle> visited;
            std::map<FT,Vertex_handle> ring;
            std::vector<FT> sq_distances;

            //find first ring neighbours
            visited.insert(vi);
            find_ring_vertices(vi->point(), vi, visited, ring);
            typename std::map<FT, Vertex_handle>::iterator it = ring.begin();
            while(it != ring.end())
              points.push_back((*it++).second);

            points.push_back(vi);

            FT dist_sum = 0.;
            for(std::size_t j=0; j<points.size(); ++j)
            {
              sq_distances.push_back(CGAL::squared_distance(vi->point(), points[j]->point()));
              dist_sum += std::sqrt(sq_distances[j]);
            }
            FT sigma = dist_sum/((FT) points.size());
            FT coeff = -1./(2*sigma*sigma);

            //coeff W_sigma
            std::vector<FT> wsigmas(points.size());
            FT w_sigmas_sum = 0.;

            for(std::size_t j=0; j<points.size(); ++j)
            {
              FT sq_dist = sq_distances[j];
              FT w = std::exp(coeff * sq_dist);
              w_sigmas_sum += w;
              wsigmas[j] = w;
            }

            for(std::size_t j=0; j<points.size(); ++j)
              wsigmas[j] /= w_sigmas_sum;

            //compute diffusion for each coeff of the tensor
            Eigen::Matrix3d new_tensor_at_vi = Eigen::Matrix3d::Zero();

            for(int j = 0; j < 3; j++)
              for(int k = 0; k < 3; k++)
                for(std::size_t l=0; l<points.size(); ++l)
                  new_tensor_at_vi(j,k) += wsigmas[l] * m_metrics[points[l]->tag()](j,k); //neighbours

            //compute new principal directions & curvature values
            double e0,e1,e2;
            Vector_3 v0,v1,v2;

            Eigen::EigenSolver<Eigen::Matrix3d> es(new_tensor_at_vi, true/*compute eigenvectors and values*/);
            const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvalueType& vals = es.eigenvalues();
            const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvectorsType& vecs = es.eigenvectors();

            e0 = (std::real(vals[0]) >= 0.) ? std::real(vals[0]) : 0.;
            e1 = (std::real(vals[1]) >= 0.) ? std::real(vals[1]) : 0.;
            e2 = (std::real(vals[2]) >= 0.) ? std::real(vals[2]) : 0.;

            v0 = get_eigenvector(vecs.col(0));
            v1 = get_eigenvector(vecs.col(1));
            v2 = get_eigenvector(vecs.col(2));

            Vector_3 ni = vi->normal();
            FT sp0 = std::abs(v0*ni);
            FT sp1 = std::abs(v1*ni);
            FT sp2 = std::abs(v2*ni);
            FT max_sp = (std::max)(sp0,(std::max)(sp1,sp2));

            FT temp_e, temp_sp;
            Vector_3 temp_v;

            if(sp1 == max_sp)
            {
              temp_e = e0;
              e0 = e1;
              e1 = e2;
              e2 = temp_e;

              temp_v = v0;
              v0 = v1;
              v1 = v2;
              v2 = temp_v;

              temp_sp = sp0;
              sp0 = sp1;
              sp1 = sp2;
              sp2 = temp_sp;
            }
            else if(sp2 == max_sp)
            {
              temp_e = e0;
              e0 = e2;
              e2 = e1;
              e1 = temp_e;

              temp_v = v0;
              v0 = v2;
              v2 = v1;
              v1 = temp_v;

              temp_sp = sp0;
              sp0 = sp2;
              sp2 = sp1;
              sp1 = temp_sp;
            }

            //projection on the tangent plane
            v1 = v1-(v1*ni)*ni;
            v2 = v2-(v2*ni)*ni;

            v1 = v1/CGAL::sqrt(v1*v1);
            v2 = v2/CGAL::sqrt(v2*v2);

            //orthogonalize v2 to v1...
            v2 = v2 - (v2*v1)*v1;
            v2 = v2/CGAL::sqrt(v2*v2);

#ifdef ANISO_DEBUG_HEAT_KERNEL
            std::cout << "point is : " << vi->tag() << " || " << vi->point() << std::endl;

            std::cout << points.size() << " neighbours" << std::endl;
            for(std::size_t j=0; j<points.size(); ++j)
                std::cout << " " << points[j]->tag() << " ";
            std::cout << std::endl;

            std::cout << "sigma sum : " << w_sigmas_sum << std::endl;
            for(std::size_t j=0; j<points.size(); ++j)
                std::cout << wsigmas[j] << " ";
            std::cout << std::endl;

            std::cout << "old tensor : " << std::endl << m_metrics[i] << std::endl;
            std::cout << "new tensor : " << std::endl << new_tensor_at_vi << std::endl;

            std::cout << "ni : " << ni << std::endl;
            std::cout << "evs & vs: " << std::endl;
            std::cout << e0 << " || " << v0 << " || " << sp0 << std::endl; //should be the normal in new basis
            std::cout << e1 << " || " << v1 << " || " << sp1 << std::endl;
            std::cout << e2 << " || " << v2 << " || " << sp2 << std::endl;
#endif

            //new metric
            Metric_base<K, K> M(ni, v1, v2, std::sqrt(e0), std::sqrt(e1), std::sqrt(e2), epsilon);
            Eigen::Matrix3d transf = M.get_transformation();
            temp_m_metrics.push_back(transf.transpose()*transf);

 #ifdef CGAL_DEBUG_OUTPUT_VECTOR_FIELD
            if(smooth_step_n == 0)
            {
              Point_3 pi = vi->point();
              FT f = get_bounding_radius()*0.02;
              (*normals_field_os) << "2 " << pi << " " << (pi+f*v0) << std::endl;
              if(e1 > e2)
              {
                (*max_vector_field_os) << "2 " << pi << " " << (pi+f*v1) << std::endl;
                (*min_vector_field_os) << "2 " << pi << " " << (pi+f*v2) << std::endl;
              }
              else
              {
                (*min_vector_field_os) << "2 " << pi << " " << (pi+f*v1) << std::endl;
                (*max_vector_field_os) << "2 " << pi << " " << (pi+f*v2) << std::endl;
              }
            }
#endif
          }

          //update all the metrics
          m_metrics = temp_m_metrics;
        }

        std::cout << "smoothed the metric field!" << std::endl;
      }

      void gl_draw_tiling(const Plane_3& plane,
                          const int star_id = -1) const
      {
        bool draw_all = (star_id < 0);
        std::size_t start = draw_all ? 0 : star_id;
        std::size_t N = draw_all ? m_vertices.size() : (star_id + 1);

        for(std::size_t i = start; i < N; i++)
        {
          Vertex_handle vi = m_vertices[i];
          /*
          if(!is_above_plane(plane, vi->point()))
            continue;
          */

          ::glDisable(GL_LIGHTING);
          ::glColor3ub( rand()%255, rand()%255, rand()%255 );

          ::glBegin(GL_POINTS);
          ::glVertex3d(vi->point().x(), vi->point().y(), vi->point().z());
          ::glEnd();

          std::vector<Point_3> tile_points = vi->tile_points();

          ::glBegin(GL_LINE_LOOP);
          for (std::size_t j = 0; j < tile_points.size(); ++j)
          {
            ::glVertex3d(tile_points[j].x(), tile_points[j].y(), tile_points[j].z());
          }
          ::glEnd();
          ::glEnable(GL_LIGHTING);
        }
      }

      void anisotropy_bounding_with_tiling(const FT& epsilon)
      {
        time_t start = time(NULL);

        Tiling<Polyhedron> Aniso_tiling;

        double tolerance_angle = 0.3;
        std::ifstream tiling_input("tiling_angle.txt");
        tiling_input >> tolerance_angle;

        if(tolerance_angle < 0)
          return;
        else
          std::cout << "tolerance is : " << tolerance_angle << std::endl;

        for (std::size_t i = 0; i < m_vertices.size(); ++i)
        {
          Vertex_handle vi = m_vertices[i];
          std::cout << "i : " << i << std::endl;
          //check on curvatures here, skip if it's ~ a plane

          //compute the aniso_ratio
          Aniso_tiling.compute(vi, tolerance_angle);

          FT vi_eps = vi->aniso_ratio();

          std::cout << "vi eps : " << vi_eps << std::endl;

          //compute new principal directions & curvature values
          double e0,e1,e2;
          Vector_3 v0,v1,v2;
          Eigen::Matrix3d m = m_metrics[i];

          Eigen::EigenSolver<Eigen::Matrix3d> es(m, true/*compute eigenvectors and values*/);
          const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvalueType& vals = es.eigenvalues();
          const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvectorsType& vecs = es.eigenvectors();

          e0 = std::real(vals[0]); //should all be > 0 here
          e1 = std::real(vals[1]);
          e2 = std::real(vals[2]);
          assert(e0 > 0 && e1 > 0 && e2 > 0);

          e0 = std::sqrt(e0);
          e1 = std::sqrt(e1);
          e2 = std::sqrt(e2);

          v0 = get_eigenvector(vecs.col(0));
          v1 = get_eigenvector(vecs.col(1));
          v2 = get_eigenvector(vecs.col(2));

          std::cout << "es : " << e0 << " " << e1 << " " << e2 << std::endl;

          FT emax = (std::max)((std::max)(e0, e1), e2);
          FT emin = (std::min)((std::min)(e0, e1), e2);

          if(emin > (emax/vi_eps)) // don't want to increase the anisotropy
          {
            std::cout << "increasing aniso ignored" << std::endl;
            continue;
          }

          if(emin == e0)
            e0 = emax / vi_eps;
          if(emin == e1)
            e1 = emax / vi_eps;
          if(emin == e2)
            e2 = emax / vi_eps;

          std::cout << "es : " << e0 << " " << e1 << " " << e2 << std::endl;
          std::cout << "v0 : " << v0 << std::endl;
          std::cout << "v1 : " << v1 << std::endl;
          std::cout << "v2 : " << v2 << std::endl;

          Metric_base<K, K> M(v0, v1, v2, e0, e1, e2, epsilon);
          Eigen::Matrix3d transf = M.get_transformation();
          m_metrics[i] = transf.transpose()*transf;
        }

        time_t end = time(NULL);
        double time_seconds = end - start;

        std::cout << "duration: " << time_seconds << " s." << std::endl;
      }

      struct Facet_Compare
      {
        private:
          Point_3 ref_point;
        public:
        bool operator()(const Facet_handle& fh1, const Facet_handle& fh2) const
        { //compare facets by comparing the furthest points from the ref point.
          Point_3 fh1_a = fh1->halfedge()->vertex()->point();
          FT sq_d_p1a = CGAL::squared_distance(ref_point, fh1_a);
          Point_3 fh1_b = fh1->halfedge()->next()->vertex()->point();
          FT sq_d_p1b = CGAL::squared_distance(ref_point, fh1_b);
          Point_3 fh1_c = fh1->halfedge()->next()->next()->vertex()->point();
          FT sq_d_p1c = CGAL::squared_distance(ref_point, fh1_c);
          FT sq_d_p1 = (std::max)((std::max)(sq_d_p1a, sq_d_p1b), sq_d_p1c);

          Point_3 fh2_a = fh2->halfedge()->vertex()->point();
          FT sq_d_p2a = CGAL::squared_distance(ref_point, fh2_a);
          Point_3 fh2_b = fh2->halfedge()->next()->vertex()->point();
          FT sq_d_p2b = CGAL::squared_distance(ref_point, fh2_b);
          Point_3 fh2_c = fh2->halfedge()->next()->next()->vertex()->point();
          FT sq_d_p2c = CGAL::squared_distance(ref_point, fh2_c);
          FT sq_d_p2 = (std::max)((std::max)(sq_d_p2a, sq_d_p2b), sq_d_p2c);

          //don't want to miss some facets due to same distance
          if(sq_d_p1 == sq_d_p2)
            return fh1 < fh2;

          //general result
          return sq_d_p1 < sq_d_p2;
        }

        Facet_Compare(const Point_3& p) : ref_point(p) {}
      };

      void find_facet_first_ring(const Vertex_handle& v,
                                 std::set<Facet_handle, Facet_Compare>& facet_first_ring)
      {
        facet_first_ring.clear();
        HV_circulator h = v->vertex_begin();
        HV_circulator hend = h;
        do
        {
          facet_first_ring.insert(h->facet());
          h++;
        }
        while(h != hend);
      }

      void find_facet_first_ring(const Facet_handle& f,
                                 std::set<Facet_handle, Facet_Compare>& facet_queue)
      {
        std::vector<Vertex_handle> v(3);
        v[0] = f->halfedge()->vertex();
        v[1] = f->halfedge()->next()->vertex();
        v[2] = f->halfedge()->next()->next()->vertex();

        for(int i=0; i<3; ++i)
        {
          HV_circulator h = v[i]->vertex_begin();
          HV_circulator hend = h;
          do
          {
            facet_queue.insert(h->facet());
            h++;
          }
          while(h != hend);
        }
      }

      bool project_on_plane_along_vector(const Exact_Point_3& p,
                                         const Exact_Plane_3& pl,
                                         const Exact_Vector_3& n,
                                         Exact_Point_3& res)
      {
        Exact_FT a = pl.a();
        Exact_FT b = pl.b();
        Exact_FT c = pl.c();
        Exact_FT d = pl.d();

        Exact_FT nx = n.x();
        Exact_FT ny = n.y();
        Exact_FT nz = n.z();

        Exact_FT lambda;
        Exact_FT den = a*nx + b*ny + c*nz;

        if(den == 0)
          return false;
        else
        {
          lambda = -(a*p.x() + b*p.y() + c*p.z() + d) / den;
          res = p + lambda * n;
          return true;
        }
      }

      bool approx_bound_reached_in_triangle(const Plane_3& pl_k1k2,
                                            const Point_3& pa,
                                            const Point_3& pb,
                                            const Point_3& pc,
                                            const FT& h)
      {
        Point_3 ppa = pl_k1k2.projection(pa);
        Point_3 ppb = pl_k1k2.projection(pb);
        Point_3 ppc = pl_k1k2.projection(pc);

        FT da = std::sqrt(CGAL::squared_distance(pa, ppa));
        FT db = std::sqrt(CGAL::squared_distance(pb, ppb));
        FT dc = std::sqrt(CGAL::squared_distance(pc, ppc));

        std::cout << "approximation at triangle points : " << std::endl;
        std::cout << pa << " || " << da << std::endl;
        std::cout << pb << " || " << db << std::endl;
        std::cout << pc << " || " << dc << std::endl;

        if((da > h) && (db > h) && (dc > h))
          std::cout << "triangle has all 3 points above approximation bound..." << std::endl;

        return ((da > h) || (db > h) || (dc > h));
      }

      bool facet_reaches_approx_bound(const Point_3& p,
                                      const Facet_handle& fh,
                                      const Plane_3& pl,
                                      const FT& h,
                                      const Segment_3& seg,
                                      Point_3& pi,
                                      std::set<Facet_handle, Facet_Compare>& facet_queue,
                                      FT& max_dist_from_p)
      {
        //project seg on the facet
        Exact_Point_3 pa = to_exact(fh->halfedge()->vertex()->point());
        Exact_Point_3 pb = to_exact(fh->halfedge()->next()->vertex()->point());
        Exact_Point_3 pc = to_exact(fh->halfedge()->next()->next()->vertex()->point());
        Exact_Triangle_3 tr_fh(pa, pb, pc);
        Exact_Plane_3 pl_fh = tr_fh.supporting_plane();
        Exact_Plane_3 pl_k1k2 = to_exact(pl);
        Exact_Point_3 pseg_ki_source, pseg_ki_target;
        Exact_Segment_3 seg_ki(to_exact(seg));
        bool skip_inter = false;

        if(!(project_on_plane_along_vector(seg_ki.source(), pl_fh, pl_k1k2.orthogonal_vector(), pseg_ki_source)))
        {
#ifdef ANISO_DEBUG_APPROX_ANISO
          std::cout << "Warning : plan of the facet is parallel to the normal of the plan made of eigen vectors" << std::endl;
#endif
          find_facet_first_ring(fh, facet_queue);
          skip_inter = true;
        }
        project_on_plane_along_vector(seg_ki.target(), pl_fh, pl_k1k2.orthogonal_vector(), pseg_ki_target);
        Exact_Segment_3 pseg_ki(pseg_ki_source, pseg_ki_target); //projected seg_ki on facet plane

#ifdef ANISO_DEBUG_APPROX_ANISO
        std::cout << std::endl << "entered facet reaches approx bound with : " << std::endl;
        std::cout << fh->halfedge()->vertex()->tag() << " " << pa << std::endl;
        std::cout << fh->halfedge()->next()->vertex()->tag() << " " << pb << std::endl;
        std::cout << fh->halfedge()->next()->next()->vertex()->tag() << " " << pc << std::endl;
        std::cout << pseg_ki_source << std::endl << pseg_ki_target << std::endl;
#endif

        if(skip_inter || CGAL::do_intersect(tr_fh, pseg_ki))
        {
          //if(approx_bound_reached_in_triangle(pl_k1k2, pa, pb, pc, h) ||
          //   (pa == p || pb == p || pc == p)) //if in the first ring around p, need to compute anyway
          {
            //compute intersection of triangle & segment
            typename KExact::Object_3 inter_obj = CGAL::intersection(tr_fh, pseg_ki);

            Exact_Point_3 pres;
            Exact_Segment_3 sres;
            if(inter_obj.assign(pres))
            {
#ifdef ANISO_DEBUG_APPROX_ANISO
              std::cout << "intersection is a point : " << pres << std::endl;
#endif

              if(pres == pa || pres == pb || pres == pc)
              {
                //it is equal to initial p, skipping this facet
                return false;
              }

              Exact_Point_3 ppres = pl_k1k2.projection(pres);
              FT sq_approx_pres = CGAL::squared_distance(back_from_exact(pres),
                                                         back_from_exact(ppres));

              if(sq_approx_pres > (h*h))
              {
                pi = back_from_exact(ppres);
                return true;
              }
            }
            else if(inter_obj.assign(sres))
            {
              Exact_Point_3 source = sres.source();
              Exact_Point_3 target = sres.target();
              Exact_Point_3 psource = pl_k1k2.projection(source);
              Exact_Point_3 ptarget = pl_k1k2.projection(target);
              FT source_h = std::sqrt(CGAL::squared_distance(back_from_exact(source),
                                                             back_from_exact(psource)));
              FT target_h = std::sqrt(CGAL::squared_distance(back_from_exact(target),
                                                             back_from_exact(ptarget)));

#ifdef ANISO_DEBUG_APPROX_ANISO
              std::cout << "intersection is a segment " << std::endl;
              std::cout << source << std::endl << target << std::endl;
#endif

              if((source_h > h) != (target_h > h)) //xor
              {
                Exact_Vector_3 svec(source, target);
                Exact_Point_3 p0 = source + ((h - source_h) / (target_h - source_h)) * svec;
                pi = back_from_exact(pl_k1k2.projection(p0));
#ifdef ANISO_DEBUG_APPROX_ANISO
                std::cout << "xor is 0-1 / 1-0 : " << h << " " << source_h << " " << target_h << std::endl;
                std::cout << "projected : " << p0 << " |-> " << pi << std::endl;
#endif
                return true;
              }
              else if((source_h > h) && (target_h > h))
              {
                //std::cout << "warning : xor 1 1!" << std::endl;
              }
              else //xor 0 0
              {//check if either intersection is furthest point from p
#ifdef ANISO_DEBUG_APPROX_ANISO
                std::cout << "xor 0 0 " << source_h << " " << target_h << std::endl;
#endif
                FT dist_psource_p = CGAL::sqrt(CGAL::squared_distance(back_from_exact(psource), p));
                FT dist_ptarget_p = CGAL::sqrt(CGAL::squared_distance(back_from_exact(ptarget), p));

                if(dist_psource_p > max_dist_from_p)
                {
                  max_dist_from_p = dist_psource_p;
                  pi = back_from_exact(psource);
#ifdef ANISO_DEBUG_APPROX_ANISO
                  std::cout << "source>maxdist : " << dist_psource_p << " new pi : " << pi << std::endl;
#endif
                }
                if(dist_ptarget_p > max_dist_from_p)
                {
                  max_dist_from_p = dist_ptarget_p;
                  pi = back_from_exact(ptarget);
#ifdef ANISO_DEBUG_APPROX_ANISO
                  std::cout << "target>maxdist : " << dist_ptarget_p << " new pi : " << pi << std::endl;
#endif
                }
              }
            }
            else
            {//no intersection but CGAL::do_intersect says there is an intersection...
              std::cout << "problem in intersection computations of projected segment and triangle" << std::endl;
            }
          }

          find_facet_first_ring(fh, facet_queue);
#ifdef ANISO_DEBUG_APPROX_ANISO
          std::cout << "intersect but does not reach approximation error : find neighbours of fh and carry on" << std::endl;
          std::cout << "added new first facet ring. Size of ring is : " << facet_queue.size() << std::endl;
#endif
        }
#ifdef ANISO_DEBUG_APPROX_ANISO
        else
          std::cout << "no intersection" << std::endl;
#endif

        return false;
      }

      void compute_aniso_ratio_at_point(const Vertex_handle& v,
                                        FT& ratio,
                                        std::vector<Point_3>& points,
                                        std::vector<Vector_3>& dirs, //k1+, k1-, k2+, k2-, n
                                        const FT& h)
      {
        Point_3 p = v->point();

#ifdef ANISO_DEBUG_APPROX_ANISO
        std::cout << "computing ratio (polyhedron) for : " << p << std::endl;
        std::cout << dirs[0] << std::endl << dirs[1] << std::endl;
        std::cout << dirs[2] << std::endl << dirs[3] << std::endl;
        std::cout << dirs[4] << std::endl;
#endif

        Plane_3 pl_k1k2(p, p + dirs[0], p + dirs[2]);

        points.clear();
        points.resize(4, p);

        //ring of facets around f
        Facet_Compare fc(p);
        std::set<Facet_handle, Facet_Compare> facet_queue(fc);

        for(int i=0; i<4; ++i) //k1+, k1-, k2+, k2-
        {
          std::set<Facet_handle> done; //facets already checked
          Vector_3 dir_i = dirs[i];
          Segment_3 seg_ki(p, p + dir_i*(get_bounding_radius() * 2.0 / std::sqrt(dir_i * dir_i)));
          FT max_distance_to_p = 0; // distance between pi and p

          find_facet_first_ring(v, facet_queue);

#ifdef ANISO_DEBUG_APPROX_ANISO
          std::cout << "direction : k" << i << " " << dirs[i] << std::endl;
          std::cout << "first ring around p has size : " << facet_queue.size() << std::endl;
#endif

          while(!(facet_queue.empty()))
          {
            Facet_handle fh = *(facet_queue.begin());
            facet_queue.erase(facet_queue.begin());

            if(done.find(fh) != done.end())
              continue;
            done.insert(fh);

            if(facet_reaches_approx_bound(p, fh, pl_k1k2, h, seg_ki, points[i], facet_queue, max_distance_to_p))
              break;
          }

          if(facet_queue.empty() && points[i] == p)
          {
            std::cout << "Warning : facet_queue empty but p[i]==p" << std::endl;
          }

          //to be safe : check that points[i] belongs to dirs[i];
          //todo
        }
        FT d0 = std::sqrt(CGAL::squared_distance(p, points[0]));
        FT d1 = std::sqrt(CGAL::squared_distance(p, points[1]));
        FT d2 = std::sqrt(CGAL::squared_distance(p, points[2]));
        FT d3 = std::sqrt(CGAL::squared_distance(p, points[3]));
        ratio = (std::max)((d0+d1)/(d2+d3), (d2+d3)/(d0+d1));

#ifdef ANISO_DEBUG_APPROX_ANISO
        std::cout << "final points : " << std::endl;
        std::cout << points[0] << std::endl << points[1] << std::endl << points[2] << std::endl << points[3] << std::endl;
        std::cout << "ds : " << d0 << " " << d1 << " " << d2 << " " << d3 << std::endl;
        std::cout << "long ratio : " << ratio << std::endl;
#endif
      }

      void compute_aniso_ratio_with_approximation(const FT& epsilon,
                                                  const FT& approximation)
      {
        std::cout << "compute aniso ratio w/ approx" << std::endl;
        for (std::size_t i = 0; i < m_vertices.size(); ++i)
        {
          Vertex_handle vi = m_vertices[i];
          std::cout << "i : " << i << std::endl;

          //check on curvatures here, skip if it's ~ a plane


          //compute principal directions & curvature values
          std::vector<double> e(3);
          std::vector<Vector_3> v(3);
          int indn = -1, ind1 = -1, ind2 = -1;
          FT maxsp = -1e308;

          Eigen::Matrix3d m = m_metrics[i];

          Eigen::EigenSolver<Eigen::Matrix3d> es(m, true/*compute eigenvectors and values*/);
          const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvalueType& vals = es.eigenvalues();
          const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvectorsType& vecs = es.eigenvectors();

          //should all be > 0 here
          e[0] = std::real(vals[0]);
          e[1] = std::real(vals[1]);
          e[2] = std::real(vals[2]);

          e[0] = std::sqrt(e[0]);
          e[1] = std::sqrt(e[1]);
          e[2] = std::sqrt(e[2]);

          v[0] = get_eigenvector(vecs.col(0));
          v[1] = get_eigenvector(vecs.col(1));
          v[2] = get_eigenvector(vecs.col(2));

#ifdef ANISO_DEBUG_APPROX_ANISO
          std::cout << "es : " << e[0] << " " << e[1] << " " << e[2] << std::endl;
          std::cout << "vs" << std::endl << v[0] << std::endl << v[1] << std::endl << v[2] << std::endl;
#endif

          FT emax = (std::max)((std::max)(e[0], e[1]), e[2]);
          FT emin = (std::min)((std::min)(e[0], e[1]), e[2]);

          Vector_3 vi_n = vi->normal();
          for(int it=0; it<3; it++)
          {
            //One of the eigenvectors should be vi_n, but Eigen doesn't always manage to find (exactly)
            //in the eigen vectors of a matrix M one of the eigenvectors M was created with...
            FT sp = std::abs(v[it]*vi_n);
            if(sp > maxsp)
            {
              maxsp = sp;
              indn = it;
            }
          }

          ind1 = (indn + 1)%3;
          ind2 = (indn + 2)%3;

#ifdef ANISO_DEBUG_APPROX_ANISO
          std::cout << "index : " << indn << " " << ind1 << " " << ind2 << std::endl;
          std::cout << "normal check : " << std::endl << v[indn] << std::endl << vi_n << std::endl;
          std::cout << "other two : " << std::endl << v[ind1] << std::endl << v[ind2] << std::endl;
#endif

          //compute local anisotropy
          FT vi_eps;
          std::vector<Vector_3> dirs(5);

          dirs[0] = v[ind1];
          dirs[1] = -v[ind1];
          dirs[2] = v[ind2];
          dirs[3] = -v[ind2];
          dirs[4] = v[indn];

          compute_aniso_ratio_at_point(vi, vi_eps, vi->tile_points(), dirs, approximation);

          std::cout << "exited vi_eps computations : " << vi->point();
          std::cout << " ratios : " << vi_eps << " " << emax/emin << std::endl;

          //compute new eigenvalues
          if(emin > (emax/vi_eps)) // don't want to increase the anisotropy
          {
            std::cout << "increasing aniso ignored"<< std::endl;
            continue;
          }

          if(emin == e[0])
            e[0] = emax / vi_eps;
          if(emin == e[1])
            e[1] = emax / vi_eps;
          if(emin == e[2])
            e[2] = emax / vi_eps;

          std::cout << "es : " << e[0] << " " << e[1] << " " << e[2] << std::endl;
          std::cout << "v0 : " << v[0] << std::endl;
          std::cout << "v1 : " << v[0] << std::endl;
          std::cout << "v2 : " << v[2] << std::endl;

          Metric_base<K, K> M(v[indn], v[ind1], v[ind2], e[indn], e[ind1], e[ind2], epsilon);
          Eigen::Matrix3d transf = M.get_transformation();
          m_metrics[i] = transf.transpose()*transf;
        }
      }

      virtual void compute_poles(std::set<Point_3>& poles) const
      {
        compute_triangulation_poles(m_c3t3, std::inserter(poles, poles.end()), get_bbox());
      }

      virtual Pointset get_surface_points(unsigned int nb,
                                          double facet_distance_coeff /*= 0.05*/) const
      {
        std::vector<Point_3> all_points;
        typename C3t3::Triangulation::Finite_vertices_iterator v;
        for(v = m_c3t3.triangulation().finite_vertices_begin();
            v != m_c3t3.triangulation().finite_vertices_end();
            ++v)
        {
          if(m_c3t3.in_dimension(v) == 1 || m_c3t3.in_dimension(v) == 2)
            all_points.push_back(v->point());
        }

        if(std::size_t(nb) > all_points.size())
        {
#ifdef ANISO_VERBOSE
          std::cout << "C3T3 mesh did not contain enough points (" << all_points.size();
          std::cout << "). Generating a new c3t3...";
          std::cout << "Facet distance is " << facet_distance_coeff << std::endl;
#endif
          FT r = this->get_bounding_radius();
          // run Mesh_3
          Mesh_criteria criteria(CGAL::parameters::facet_angle = 25.,
                                 CGAL::parameters::facet_size = r * 0.05,
                                 CGAL::parameters::facet_distance = r * facet_distance_coeff);
                                 // cell criteria are ignored
          m_c3t3 = CGAL::make_mesh_3<C3t3>(*domain, criteria,
                                           CGAL::parameters::no_perturb(),
                                           CGAL::parameters::no_exude());
#ifdef ANISO_OUTPUT_MESH_FOR_POLES
          std::ofstream out("mesh_3_temp.mesh");
          m_c3t3.output_to_medit(out);
#endif
#ifdef ANISO_OUTPUT_POINTS_FOR_MEDIAL_AXIS
          std::ofstream out_pts("mesh_3_points.pts");
          this->output_points(m_c3t3, out_pts);
#endif
          return get_surface_points(nb, facet_distance_coeff/3.);
        }

        std::random_shuffle(all_points.begin(), all_points.end());
        return Pointset(all_points.begin(), (all_points.begin() + nb));
      }

      virtual double global_max_curvature() const
      {
        if(m_cache_max_curvature)
          return m_max_curvature;

        m_max_curvature = 0.;
        m_min_curvature = DBL_MAX;
        typename C3t3::Triangulation::Finite_vertices_iterator v;
        for(v = m_c3t3.triangulation().finite_vertices_begin();
            v != m_c3t3.triangulation().finite_vertices_end();
            ++v)
        {
          double minc, maxc;
          min_max_curvatures(v->point(), minc, maxc);
          m_max_curvature = (std::max)(m_max_curvature, maxc);
          m_min_curvature = (std::min)(m_min_curvature, minc);
        }
        m_cache_max_curvature = true;
        m_cache_min_curvature = true;
        return m_max_curvature;
      }

      virtual double global_min_curvature() const
      {
        if(m_cache_min_curvature)
          return m_min_curvature;

        global_max_curvature(); //computes both max and min

        return m_min_curvature;
      }

      void min_max_curvatures(const Point_3& p,
                              double& minc, 
                              double& maxc) const
      {
        Vector_3 v0, v1, v2;
        double e0, e1, e2;
        tensor_frame(p, v0, v1, v2, e0, e1, e2);
        minc = (std::min)(e0, (std::min)(e1, e2));
        maxc = (std::max)(e0, (std::max)(e1, e2));
      }

      void initialize(const FT& epsilon,
                      const FT& en_factor = 1.0,
                      const FT& approximation = 1.0,
                      const bool smooth_metric = false)
      {
        set_aabb_tree();
        
        point_count = (int)m_polyhedron.size_of_vertices();
        m_vertices.reserve(point_count);
        m_metrics.reserve(point_count);

        compute_bounding_box();        

        std::ifstream metric_input("metrics.txt");
        if(metric_input)
          get_metrics_from_file(metric_input, epsilon, en_factor);
        else
          compute_local_metric(epsilon, en_factor, smooth_metric);

        std::ifstream useless("smoothing.txt");
        if(useless)
          heat_kernel_smoothing(epsilon);

        std::ifstream approx_aniso("approx_aniso.txt");
        std::ifstream tiling_input("tiling_angle.txt");
        if(approx_aniso)
          compute_aniso_ratio_with_approximation(epsilon, approximation);
        else if(tiling_input)
          anisotropy_bounding_with_tiling(epsilon);
      }

      void gl_draw_intermediate_mesh_3(const Plane_3& plane) const
      {
        gl_draw_c3t3<C3t3, Plane_3>(m_c3t3, plane);
      }

      Constrain_surface_3_polyhedral(const char *filename,
                                     const FT& epsilon,
                                     const FT& en_factor = 1.0,
                                     const FT& approximation = 1.0,
                                     const bool smooth_metric = false) 
        : m_vertices(),
          m_metrics(),
          m_cache_max_curvature(false), 
          m_cache_min_curvature(false)
        {
#ifdef CGAL_DEBUG_OUTPUT_VECTOR_FIELD
          min_vector_field_os = new std::ofstream("vector_field_min.polylines.cgal");
          max_vector_field_os = new std::ofstream("vector_field_max.polylines.cgal");
          normals_field_os = new std::ofstream("vector_field_normals.polylines.cgal");
#endif
          set_domain_and_polyhedron(filename);
          initialize(epsilon, en_factor, approximation, smooth_metric);
        }

      Constrain_surface_3_polyhedral(const Polyhedron& p, 
                                     const FT& epsilon,
                                     const FT& en_factor = 1.0,
                                     const FT& approximation = 1.0,
                                     const bool smooth_metric = false) 
        : m_vertices(),
          m_metrics(),
          m_cache_max_curvature(false), 
          m_cache_min_curvature(false)
        {
#ifdef CGAL_DEBUG_OUTPUT_VECTOR_FIELD
          min_vector_field_os = new std::ofstream("vector_field_min.polylines.cgal");
          max_vector_field_os = new std::ofstream("vector_field_max.polylines.cgal");
          normals_field_os = new std::ofstream("vector_field_normals.polylines.cgal");
#endif
          set_domain_and_polyhedron(p);
          initialize(epsilon, en_factor, approximation, smooth_metric);
        }

      Constrain_surface_3_polyhedral* clone() const
      {
        return new Constrain_surface_3_polyhedral(m_polyhedron);
      }

      ~Constrain_surface_3_polyhedral() 
      { 
        delete domain;
#ifdef CGAL_DEBUG_OUTPUT_VECTOR_FIELD
        (*min_vector_field_os).close();
        (*max_vector_field_os).close();
        (*normals_field_os).close();
        delete min_vector_field_os;
        delete max_vector_field_os;
        delete normals_field_os;
#endif
      }
    }; //Constrain_surface_3_polyhedral
  } //Anisotropic_mesh_3
} //CGAL

#endif //CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_POLYHEDRAL_H

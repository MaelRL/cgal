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

#include "../../demo/Polyhedron_type.h"

#include <CGAL/Constrain_surface_3.h>
#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>
#include <CGAL/Enriched_items.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <Eigen/Dense>
#include <CGAL/Metric.h>
#include <CGAL/helpers/metric_helper.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>

// for Mesh_3
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>

#include <CGAL/gl_draw/drawing_helper.h>
#include <CGAL/helpers/c3t3_polyhedron_builder.h>

#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K,
         typename Polyhedron = CGAL::Polyhedron_3<K, Enriched_items> >
class Constrain_surface_3_polyhedral :
    public Constrain_surface_3<K>
{
public:
  typedef Constrain_surface_3<K>                     Base;
  typedef Polyhedron                                 Polyhedron_type;

  typedef typename Base::Object_3                    Object_3;
  typedef typename Base::FT                          FT;
  typedef typename Base::Vector_3                    Vector_3;
  typedef typename Base::Point_3                     Point_3;
  typedef typename Base::Segment_3                   Segment_3;
  typedef typename Base::Triangle_3                  Triangle_3;
  typedef typename Base::Plane_3                     Plane_3;
  typedef typename Base::Oriented_side               Oriented_side;

  typedef typename Base::Pointset                    Point_container;
  typedef typename Base::Colored_poly                Colored_polyhedron;

  typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;

  typedef CGAL::AABB_polyhedron_triangle_primitive<K, Polyhedron> Primitive;
  typedef CGAL::AABB_traits<K, Primitive>                         Traits;
  typedef CGAL::AABB_tree<Traits>                                 Tree;
  typedef typename Tree::Object_and_primitive_id Object_and_primitive_id;
  typedef typename Tree::Point_and_primitive_id Point_and_primitive_id;

  typedef typename Polyhedron::Vertex_handle    Vertex_handle;
  typedef typename Polyhedron::Facet_handle     Facet_handle;
  typedef typename Polyhedron::Halfedge_handle  Halfedge_handle;
  typedef typename Polyhedron::Halfedge_around_vertex_circulator HV_circulator;

  // for Mesh_3
  typedef typename CGAL::Mesh_triangulation_3<Mesh_domain>::type    Tr;
  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr>               C3t3;
  typedef CGAL::Mesh_criteria_3<Tr>                                 Mesh_criteria;

  // exact
  typedef CGAL::Exact_predicates_exact_constructions_kernel KExact;
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

#ifdef CGAL_DEBUG_OUTPUT_VECTOR_FIELD
  std::ofstream* min_vector_field_os;
  std::ofstream* max_vector_field_os;
  std::ofstream* normals_field_os;
#endif

protected:
  mutable C3t3 m_c3t3;

public:
  C3t3 c3t3() const { return m_c3t3; }

  virtual FT get_bounding_radius() const { return bounding_radius; }
  virtual std::string name() const { return std::string("Polyhedron"); }
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
  struct by_distance_to_p
  {
    Point_3 p;

    bool operator()(const Point_3 &pi, const Point_3& pj)
    {
      return CGAL::squared_distance(p, pi) < CGAL::squared_distance(p, pj);
    }

    by_distance_to_p(const Point_3 &pini) : p(pini) { }
  };

  virtual Object_3 intersection(const Point_3 &p0, const Point_3 &p1, const Point_3 &ref) const
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

    std::sort(intersection_points.begin(), intersection_points.end(), by_distance_to_p(ref));

#ifdef ANISO_DEBUG_INTERSECTION
    if(intersection_points.size() > 1)
    {
      std::cout << "intersection poly : " << std::endl << p0 << std::endl << p1 << std::endl << ref << std::endl;
      std::cout << intersection_points.size() << " intersections" << std::endl;
      for(typename std::vector<Point_3>::iterator vit = intersection_points.begin(); vit != intersection_points.end(); ++vit)
        std::cout << *vit << " || " << CGAL::squared_distance(*vit, p0) << std::endl;
    }
#endif
    return make_object(intersection_points.back());;
  }

  FT compute_sq_approximation(const Point_3& p) const
  {
    return tree->squared_distance(p);
  }

  void find_neighbors(Facet_handle facet,
                      std::vector<Vertex_handle>& neighbors) const
  {
    Halfedge_handle he1 = facet->halfedge();
    Halfedge_handle he2 = he1->next();
    Halfedge_handle he3 = he2->next();
    neighbors.push_back(he1->vertex());
    neighbors.push_back(he1->opposite()->next()->vertex());
    neighbors.push_back(he2->vertex());
    neighbors.push_back(he2->opposite()->next()->vertex());
    neighbors.push_back(he3->vertex());
    neighbors.push_back(he3->opposite()->next()->vertex());
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

  virtual Point_container initial_points(const int nb) const
  {
    typedef std::vector<std::pair<typename Mesh_domain::Point_3,
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
                             std::size_t count,
                             const double max_sqd = 0.) const
  {
    std::set<Vertex_handle> visited;
    std::map<FT, Vertex_handle/*sqd to center : closest come first*/> ring;

    ring.insert(std::make_pair(0.,v));
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
                          const double max_sqd = 0.) const
  {
    HV_circulator h = v->vertex_begin();
    HV_circulator hend = h;
    do
    {
      Vertex_handle vv = h->opposite()->vertex();
      std::pair<typename std::set<Vertex_handle>::iterator, bool> is_insert_successful;
      is_insert_successful = visited.insert(vv);
      if(is_insert_successful.second)
      {
        FT sqd = CGAL::squared_distance(center, vv->point());
        if(max_sqd == 0. || sqd < max_sqd)
          ring.insert(std::make_pair(sqd,vv));
      }
      h++;
    }
    while(h != hend);
  }

  inline void find_nearest_vertices_kanle(Vertex_handle v,
                                          std::vector<Vertex_handle>& points,
                                          std::size_t count)
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

  inline void find_nearest_vertices(const Point_3 &p,
                                    std::vector<Vertex_handle> &points,
                                    const FT max_dist) const
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
    {
      std::cout << "\n Input file does not exist" << std::endl;
      exit(0);
    }
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

  virtual void compute_poles(std::set<Point_3>& poles) const
  {
    compute_triangulation_poles(m_c3t3, std::inserter(poles, poles.end()), get_bbox());
  }

  void build_colored_polyhedron(Colored_polyhedron& poly) const
  {
    std::cout << "trying to convert c3t3 to polyhedron and outputing it (P)...";
    Complex_3_in_triangulation_3_polyhedron_builder<C3t3, Colored_polyhedron> builder(m_c3t3);
    poly = Colored_polyhedron();
    poly.delegate(builder);
    std::ofstream out("colored_poly.off");
    out << poly;
    std::cout << "done" << std::endl;
  }

  virtual Point_container get_surface_points(unsigned int nb,
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
                             CGAL::parameters::facet_distance = r * facet_distance_coeff,
                             CGAL::parameters::facet_topology = MANIFOLD);
                             // cell criteria are ignored
      m_c3t3 = CGAL::make_mesh_3<C3t3>(*domain, criteria,
                                       CGAL::parameters::no_perturb(),
                                       CGAL::parameters::no_exude());
#if 1//def ANISO_OUTPUT_MESH_FOR_POLES
      std::ofstream out("mesh_3_temp.mesh");
      m_c3t3.output_to_medit(out);
#endif
#ifdef ANISO_OUTPUT_POINTS_FOR_MEDIAL_AXIS
      std::ofstream out_pts("mesh_3_points.pts");
      this->output_points(m_c3t3, out_pts);
#endif
      return get_surface_points(nb, facet_distance_coeff/3.);
    }
#ifdef ANISO_VERBOSE
    std::cout << "c3t3 has enough points : " << all_points.size() << " faces : " << m_c3t3.number_of_facets_in_complex() << std::endl;
#endif
    std::random_shuffle(all_points.begin(), all_points.end());
    return Point_container(all_points.begin(), (all_points.begin() + nb));
  }

  void initialize()
  {
    set_aabb_tree();

    point_count = (int)m_polyhedron.size_of_vertices();
    m_vertices.reserve(point_count);

    compute_bounding_box();
  }

  void gl_draw_intermediate_mesh_3(const Plane_3& plane) const
  {
    gl_draw_c3t3<C3t3, Plane_3>(m_c3t3, plane);
  }

  Constrain_surface_3_polyhedral(const char *filename)
    :
      m_vertices()
  {
#ifdef CGAL_DEBUG_OUTPUT_VECTOR_FIELD
    min_vector_field_os = new std::ofstream("vector_field_min.polylines.cgal");
    max_vector_field_os = new std::ofstream("vector_field_max.polylines.cgal");
    normals_field_os = new std::ofstream("vector_field_normals.polylines.cgal");
#endif
    set_domain_and_polyhedron(filename);
    initialize();
  }

  Constrain_surface_3_polyhedral(const Polyhedron& p)
    :
      m_vertices()
  {
#ifdef CGAL_DEBUG_OUTPUT_VECTOR_FIELD
    min_vector_field_os = new std::ofstream("vector_field_min.polylines.cgal");
    max_vector_field_os = new std::ofstream("vector_field_max.polylines.cgal");
    normals_field_os = new std::ofstream("vector_field_normals.polylines.cgal");
#endif
    set_domain_and_polyhedron(p);
    initialize();
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

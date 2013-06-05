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


// a refined vertex with a tag
template <class Refs, class T, class P>
class Enriched_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P>
{
  // tag
  std::size_t m_tag;

public:
  // life cycle
  Enriched_vertex()  {}
  // repeat mandatory constructors
  Enriched_vertex(const P& pt)
    : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt),
      m_tag(0)  {}

  // tag
  std::size_t& tag() {  return m_tag; }
  const std::size_t& tag() const {  return m_tag; }
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

      typedef typename Base::FT                          FT;
      typedef typename Base::Vector_3                    Vector_3;
      typedef typename Base::Point_3                     Point_3;
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
      std::vector<Vector_3> m_normals;

#ifdef CGAL_DEBUG_OUTPUT_VECTOR_FIELD
      std::ofstream min_vector_field_os;
      std::ofstream max_vector_field_os;
      std::ofstream normals_field_os;
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
      FT compute_sq_approximation(const Point_3& p) const
      {
        return tree->squared_distance(p);
      }

      virtual void tensor_frame_on_input_point(Vertex_handle v, 
        Vector_3 &v0/*n*/, Vector_3 &v1/*vmax*/, Vector_3 &v2/*vmin*/,//unit vectors
        FT& c1/*cmax*/, FT& c2/*cmin*/,
        const FT& epsilon) 
      {
        std::vector<Point_3> points;
        find_nearest_vertices(v, points, 36, 0.);
//          nearest_start_try_radius*nearest_start_try_radius);
        points.insert(points.begin(), v->point());

        My_Monge_form monge_form;
        My_Monge_via_jet_fitting monge_fit;
        monge_form = monge_fit(points.begin(), points.end(), 4, 4);

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
                                 std::vector<Point_3>& points,
                                 int count,
                                 const double max_sqd)
      {
        std::set<Vertex_handle> visited;
        std::list<Vertex_handle> ring;
        find_ring_vertices(v->point(), v, points, visited, ring);

        while(points.size() < count && !ring.empty())
        {
          Vertex_handle vi = ring.front();
          ring.pop_front();
          find_ring_vertices(v->point(), vi, points, visited, ring, max_sqd);
        }
      }

      void find_ring_vertices(const Point_3& center,
                              Vertex_handle v,
                              std::vector<Point_3>& points,
                              std::set<Vertex_handle>& visited,
                              std::list<Vertex_handle>& ring,
                              const double max_sqd = 0.)/**/
      {
        HV_circulator h = v->vertex_begin();
        HV_circulator hend = h;
        do
        {
          //the 1-ring should be complete, even if nbv > count
          Vertex_handle vv = h->opposite()->vertex();
          if(visited.find(vv) == visited.end())
          {
            visited.insert(vv);
            if(max_sqd == 0. || CGAL::squared_distance(center, vv->point()) < max_sqd)
            {
              ring.push_back(vv);
              points.push_back(vv->point());
            }
          }
          h++;
        }
        while(h != hend);
      }

      inline void find_nearest_vertices_kanle(Vertex_handle v, 
                                        std::vector<Point_3>& points, 
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

      inline void find_nearest_vertices(const Point_3 &p, std::vector<Point_3> &points, const FT max_dist) 
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

            points.push_back(q);

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

      void compute_bounding_box(/*std::vector<Point_3>& points*/)
      {
        std::cout << "\nComputing bounding box..." << std::endl;

        /*typename Polyhedron::Vertex_iterator vi = m_polyhedron.vertices_begin();
        typename Polyhedron::Vertex_iterator viend = m_polyhedron.vertices_end();
        for (; vi != viend; vi++)
          points.push_back(vi->point());
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

      void get_metrics_from_file(std::ifstream& input, const FT& epsilon)
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

        for(size_t i=0; i<nv; ++i)
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

          en = (std::max)(std::abs(e1), std::abs(e2));
          m_normals.push_back(vn);

          Metric_base<K, K> M(vn, v1, v2, en, e1, e2, epsilon);
          Eigen::Matrix3d transf = M.get_transformation();
          m_metrics.push_back(transf.transpose()*transf);
        }
      }

      void compute_local_metric(const FT& epsilon)
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

          const Point_3& pi = m_vertices[i]->point();
//          tensor_frame_on_point(pi, vn, v1, v2, e1, e2, epsilon);

          out << pi.x() << " " << pi.y() << " " << pi.z() << std::endl;
          out << e1 << " " << e2 << "     ";
          out << v1.x() << " " << v1.y() << " " << v1.z() << "     ";
          out << v2.x() << " " << v2.y() << " " << v2.z() << "     ";
          out << vn.x() << " " << vn.y() << " " << vn.z() << std::endl;

          en = (std::max)(e1, e2);

          m_max_curvature = (std::max)(m_max_curvature, en);
          m_min_curvature = (std::min)(m_min_curvature, (std::min)(e1, e2));

#ifdef CGAL_DEBUG_OUTPUT_VECTOR_FIELD
          FT f = get_bounding_radius()*0.05;
          max_vector_field_os << "2 " << pi << " " << (pi+f*v1) << std::endl;
          min_vector_field_os << "2 " << pi << " " << (pi+f*v2) << std::endl;
          normals_field_os << "2 " << pi << " " << (pi+f*vn) << std::endl;
#endif

          std::size_t index = m_vertices[i]->tag();
          if(i != index) 
            std::cerr << "Error1 in indices in compute_local_metric!" << std::endl;

          m_normals.push_back(vn);
          if(m_normals.size() != i+1 )
            std::cerr << "Error2 in indices in compute_local_metric!" << std::endl;

          Metric_base<K, K> M(vn, v1, v2, en, e1, e2, epsilon);
          Eigen::Matrix3d transf = M.get_transformation();
          m_metrics.push_back(transf.transpose()*transf);

          //if(m_metrics.size() != i+1 )
          //  std::cerr << "Error3 in indices in compute_local_metric!" << std::endl;
        }

        m_cache_max_curvature = true;
        m_cache_min_curvature = true;
        std::cout << "done." << std::endl;
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

      void initialize(const FT& epsilon, const bool smooth_metric = false)
      {
        set_aabb_tree();
        
        point_count = (int)m_polyhedron.size_of_vertices();
        m_vertices.reserve(point_count);
        m_metrics.reserve(point_count);
        m_normals.reserve(point_count);

        compute_bounding_box();        

        std::ifstream metric_input("metrics.txt");
        if(metric_input)
          get_metrics_from_file(metric_input, epsilon);
        else
          compute_local_metric(epsilon);

      }

      void gl_draw_intermediate_mesh_3(const Plane_3& plane) const
      {
        gl_draw_c3t3<C3t3, Plane_3>(m_c3t3, plane);
      }

      Constrain_surface_3_polyhedral(const char *filename, 
                                     const FT& epsilon,
                                     const bool smooth_metric = false) 
        : m_vertices(),
          m_metrics(),
          m_normals(),
          m_cache_max_curvature(false), 
          m_cache_min_curvature(false)
        {
#ifdef CGAL_DEBUG_OUTPUT_VECTOR_FIELD
          min_vector_field_os = std::ofstream("vector_field_min.polylines.cgal");
          max_vector_field_os = std::ofstream("vector_field_max.polylines.cgal");
          normals_field_os = std::ofstream("vector_field_normals.polylines.cgal");
#endif
          set_domain_and_polyhedron(filename);
          initialize(epsilon, smooth_metric);
        }

      Constrain_surface_3_polyhedral(const Polyhedron& p, 
                                     const FT& epsilon,
                                     const bool smooth_metric = false) 
        : m_vertices(),
          m_metrics(),
          m_normals(),
          m_cache_max_curvature(false), 
          m_cache_min_curvature(false)
        {
#ifdef CGAL_DEBUG_OUTPUT_VECTOR_FIELD
          min_vector_field_os = std::ofstream("vector_field_min.polylines.cgal");
          max_vector_field_os = std::ofstream("vector_field_max.polylines.cgal");
          normals_field_os = std::ofstream("vector_field_normals.polylines.cgal");
#endif
          set_domain_and_polyhedron(p);
          initialize(epsilon, smooth_metric);
        }

      Constrain_surface_3_polyhedral* clone() const
      {
        return new Constrain_surface_3_polyhedral(m_polyhedron);
      }

      ~Constrain_surface_3_polyhedral() 
      { 
        delete domain;
#ifdef CGAL_DEBUG_OUTPUT_VECTOR_FIELD
        min_vector_field_os.close();
        max_vector_field_os.close();
        normals_field_os.close();
#endif
      }
    }; //Constrain_surface_3_polyhedral
  } //Anisotropic_mesh_3
} //CGAL

#endif //CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_POLYHEDRAL_H

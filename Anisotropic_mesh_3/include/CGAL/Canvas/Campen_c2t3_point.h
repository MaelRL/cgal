#ifndef CGAL_ANISOTROPIC_MESH_3_CAMPEN_C2T3_POINT_H
#define CGAL_ANISOTROPIC_MESH_3_CAMPEN_C2T3_POINT_H

#include <CGAL/Canvas/canvas_config.h>
#include <CGAL/Canvas/canvas_enum.h>
#include <CGAL/Canvas/canvas_point.h>
#include <CGAL/Canvas/canvas_helper.h>

#include <CGAL/Metric.h>
#include <CGAL/helpers/metric_helper.h>

#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/assertions.h>

#include <boost/array.hpp>
#include <boost/unordered_map.hpp>
#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K, typename Metric_field, typename C3t3>
class Campen_canvas;

template<typename K, typename Metric_field, typename C3t3>
class Campen_canvas_point :
    public Canvas_point<K, Campen_canvas<K, Metric_field, C3t3> >
{
private:
  typedef Campen_canvas_point<K, Metric_field, C3t3>                  Self;

public:
  typedef std::vector<Campen_canvas_point>            Campen_canvas_point_vector;
  typedef Self*                                       Campen_canvas_point_handle;

  typedef int                          Vertex_Info; // index of the canvas point
  typedef int                          Cell_Info; // index of the subdomain

  typedef Campen_canvas<K, Metric_field, C3t3>                        Canvas;
  typedef Canvas_point<K, Canvas>                                     Base;

  typedef typename C3t3::Triangulation                                Tr;
  typedef typename Tr::Geom_traits                                    Gt;

  typedef typename K::FT                                              FT;
  typedef typename Gt::Weighted_point                                 Weighted_point_3;
  typedef typename Gt::Bare_point                                     Point_3;
  typedef typename Gt::Vector_3                                       Vector_3;
  typedef typename Base::Metric                                       Metric;

  typedef Eigen::Matrix<FT, 2, 1>                                     Vector2d;
  typedef typename Base::Vector3d                                     Vector3d;

  typedef typename Tr::Vertex_handle                                  Vertex_handle;
  typedef std::vector<Vertex_handle>                                  Vertex_handle_vector;
  typedef typename Vertex_handle_vector::iterator                     Vertex_handle_handle;
  typedef typename Tr::Edge                                           Edge;
  typedef typename Tr::Facet                                          Facet;
  typedef std::vector<Facet>                                          Facet_vector;
  typedef typename Facet_vector::iterator                             Facet_handle;
  typedef typename Tr::Cell_handle                                    Cell_handle;

  Vertex_handle m_v; // corresponding vertex in the c3t3

  bool ignore_children;

  mutable bool m_is_vertices_cache_dirty, m_is_facets_cache_dirty;
  mutable Vertex_handle_vector m_adjacent_vertices_cache; // on the complex boundary
  mutable Facet_vector m_incident_facets_cache; // on the complex boundary

  C3t3& c3t3() { return this->canvas()->m_c3t3; }
  const C3t3& c3t3() const { return this->canvas()->m_c3t3; }

  // cache related functions
  struct In_complex_filter
  {
    const C3t3& c;

    In_complex_filter(const C3t3& c_) : c(c_) { }
    bool operator() (const Facet& f) const
    {
      // filters facets that are not on the boundary of the complex
      return !(c.is_in_complex(f));
    }
  };

  template<class OutputIterator>
  OutputIterator
  incident_facets_on_complex_boundary(OutputIterator facets) const
  {
    return c3t3().triangulation().tds().incident_facets(m_v, facets,
                                                        In_complex_filter(c3t3()));
  }

  void set_adjacent_vertices_on_complex_boundary() const
  {
    // brute forcy...

    CGAL_assertion(m_v != Vertex_handle());
    CGAL_assertion(c3t3().triangulation().dimension() > 1);

    // copy
    update_facets_cache();
    Facet_vector facets = m_incident_facets_cache;
    CGAL_assertion(facets.size() > 2);

    // decide the orientation from the first facet
    const Facet& front = facets.front();
    Cell_handle c = front.first;
    std::size_t second = front.second;

    if(!c3t3().is_in_complex(c))
    {
      Cell_handle c_mirror = c->neighbor(second);
      second = c_mirror->index(c);
      c = c_mirror;
    }
    CGAL_assertion(!c3t3().triangulation().is_infinite(c));

    // get the two other vertices of the face
    // todo: there might be something prettier here
    Vertex_handle vh_a, vh_b;
    for(std::size_t i=0; i<4; ++i)
    {
      if(i == second || c->vertex(i) == m_v)
        continue;
      if(vh_a == Vertex_handle())
        vh_a = c->vertex(i);
      else
        vh_b = c->vertex(i);
    }
    CGAL_postcondition(vh_a != Vertex_handle() && vh_b != Vertex_handle());

    Vertex_handle curr, last;
    // Now look at the orientation of the tet (a, b, m_v, second)
    CGAL::Orientation o = CGAL::orientation(vh_a->point().point(),
                                            vh_b->point().point(),
                                            m_v->point().point(),
                                            c->vertex(second)->point().point());
    if(o == CGAL::POSITIVE) // ((ab^ac),ad) > 0
    {
      // since we want the normal to point outside
      curr = vh_a;
      last = vh_b;
    }
    else
    {
      CGAL_assertion(o == CGAL::NEGATIVE); // ZERO case should not happen
      curr = vh_b;
      last = vh_a;
    }

/*
    std::cout << "----------------------" << std::endl;
    std::cout << "all available: " << m_v->info() << std::endl;
    for(std::size_t i=0; i<facets.size(); ++i)
    {
      const Facet& f = facets[i];
      const Cell_handle c = f.first;

      for(int j=0; j<4; ++j)
        std::cout << c->vertex(j)->info() << " ";

      std::cout << "( " << f.second << " )" << std::endl;
    }
*/

    facets.erase(facets.begin());
    m_adjacent_vertices_cache.push_back(curr);

    // ugly part
    Vertex_handle prev_curr = curr;
    while(!facets.empty())
    {
//      std::cout << "looking for " << curr->info() << std::endl;

      typename Facet_vector::iterator it = facets.begin();
      typename Facet_vector::iterator end = facets.end();
      for(; it!=end; ++it)
      {
        c = it->first;
        second = it->second;

        // get the two other vertices
        Vertex_handle v1, v2;
        for(std::size_t i=0; i<4; ++i)
        {
          if(i == second || c->vertex(i) == m_v)
            continue;
          if(v1 == Vertex_handle())
            v1 = c->vertex(i);
          else
            v2 = c->vertex(i);
        }

//        std::cout << "available: " << v1->info() << " and " << v2->info() << std::endl;

        // look for curr
        if(v1 == curr)
        {
          curr = v2;
          facets.erase(it);
          break;
        }
        else if(v2 == curr)
        {
          curr = v1;
          facets.erase(it);
          break;
        }
      }

      CGAL_assertion(prev_curr != curr);

      prev_curr = curr;
      m_adjacent_vertices_cache.push_back(curr);
    }

    CGAL_assertion(m_adjacent_vertices_cache.back() == last);
    CGAL_assertion(m_adjacent_vertices_cache.size() == m_incident_facets_cache.size());
  }

  void set_adjacent_vertices_on_complex_boundary_old() const
  {
    // doesn't work if the dihedral angles are too small or too large.

    // We're going to need to have them oriented to compute the relative angles
    // that are needed to flatten a path along the surface, so we might as well
    // do it here !

    // We orient them so the normal points outward from the surface

    CGAL_assertion(m_v != Vertex_handle());
    CGAL_assertion(c3t3().triangulation().dimension() > 1);

    // Build a map<Vh, Vh> (each entry is one edge of the link)
    // and we order the vertices by walking the link (the map) at the end
    boost::unordered_multimap<Vertex_handle, Vertex_handle> oriented_edges;

    Facet_handle fit = incident_facets_in_complex_begin();
    Facet_handle end = incident_facets_in_complex_end();
    std::size_t nb_of_inc_facets = end - fit;

    for(; fit!=end; ++fit)
    {
      Cell_handle c = fit->first;
      std::size_t second = fit->second;

      // we want to have the cell that is inside the complex to avoid the inf vertex
      if(!c3t3().is_in_complex(c))
      {
        Cell_handle c_mirror = c->neighbor(second);
        second = c_mirror->index(c);
        c = c_mirror;
      }
      CGAL_assertion(!c3t3().triangulation().is_infinite(c));

      // get the two other vertices of the face
      // todo: there might be something prettier here
      Vertex_handle vh_a, vh_b;
      for(std::size_t i=0; i<4; ++i)
      {
        if(i == second || c->vertex(i) == m_v)
          continue;


        if(vh_a == Vertex_handle())
          vh_a = c->vertex(i);
        else
          vh_b = c->vertex(i);
      }
      CGAL_postcondition(vh_a != Vertex_handle() && vh_b != Vertex_handle());

      // Now look at the orientation of the tet (a, b, m_v, second)
      CGAL::Orientation o = CGAL::orientation(vh_a->point().point(),
                                              vh_b->point().point(),
                                              m_v->point().point(),
                                              c->vertex(second)->point().point());
      if(o == CGAL::POSITIVE) // ((ab^ac),ad) > 0
      {

        oriented_edges[vh_b] = vh_a; // since we want the normal to point outside
      }
      else
      {
        CGAL_assertion(o == CGAL::NEGATIVE); // ZERO case should not happen
        oriented_edges[vh_a] = vh_b;
      }
    }
    CGAL_postcondition(nb_of_inc_facets == oriented_edges.size());

    // We can now put the (oriented) vertices in cache
    typename boost::unordered_multimap<Vertex_handle, Vertex_handle>::iterator
                                                    it = oriented_edges.begin();
    Vertex_handle init = it->first;

    while(true)
    {
      Vertex_handle edge_beginning = it->first;

      m_adjacent_vertices_cache.push_back(edge_beginning);
      it = oriented_edges.find(it->second);
      CGAL_postcondition(it != oriented_edges.end());

      if(it->first == init)
        break;
    }
    CGAL_postcondition(m_adjacent_vertices_cache.size() == oriented_edges.size());
  }

  void update_vertices_cache() const
  {
    if(!m_is_vertices_cache_dirty)
      return;

    m_adjacent_vertices_cache.clear();
    set_adjacent_vertices_on_complex_boundary();

    m_is_vertices_cache_dirty = false;
  }

  void update_facets_cache() const
  {
    if(!m_is_facets_cache_dirty)
      return;

    m_incident_facets_cache.clear();
    incident_facets_on_complex_boundary(std::back_inserter(m_incident_facets_cache));

    m_is_facets_cache_dirty = false;
  }

  void update_canvas_point_caches() const
  {
    update_facets_cache();
    update_vertices_cache();
  }

  void invalidate_caches() const
  {
    m_is_vertices_cache_dirty = true;
    m_is_facets_cache_dirty = true;
  }

  inline Vertex_handle_handle adjacent_vertices_in_complex_begin() const
  {
    update_vertices_cache();
    return m_adjacent_vertices_cache.begin();
  }

  inline Vertex_handle_handle adjacent_vertices_in_complex_end() const
  {
    CGAL_assertion(!m_is_vertices_cache_dirty);
    return m_adjacent_vertices_cache.end();
  }

  inline Facet_handle incident_facets_in_complex_begin() const
  {
    update_facets_cache();
    return m_incident_facets_cache.begin();
  }

  inline Facet_handle incident_facets_in_complex_end() const
  {
    CGAL_assertion(!m_is_facets_cache_dirty);
    return m_incident_facets_cache.end();
  }

  FT compute_flattening_angle(std::size_t in, std::size_t center, std::size_t out) const
  {
    // query the canvas' angle map to know the angle between [in-center] and [center-out]
    typename Canvas::Angle_map::const_iterator mit =
                      this->canvas()->m_angles.find(std::make_pair(center, in));
    CGAL_postcondition(mit != this->canvas()->m_angles.end());
    FT angle_in = mit->second;

    mit = this->canvas()->m_angles.find(std::make_pair(center, out));
    CGAL_postcondition(mit != this->canvas()->m_angles.end());
    FT angle_out = mit->second;

    return std::fmod(angle_out - angle_in + 2.*CGAL_PI, 2.*CGAL_PI);
  }

  Point_3 compute_unfolded_coordinates(const Point_3& in,
                                       const Point_3& center,
                                       const Vector3d& normal,
                                       const FT length, const FT angle) const
  {
    // 'in' and 'center' are unfolded point coordinates
    // knowing the length and the angle, we compute the new unfolded point

#if (VERBOSITY > 20)
    std::cout << "Computing unfolding of new point..." << std::endl;
    std::cout << "Previous two: " << in << " || " << center << std::endl;
    std::cout << "normal: " << normal.transpose() << std::endl;
    std::cout << "point at length: " << length << " and angle: " << angle << std::endl;
#endif

    Vector3d v;
    v(0) = in.x() - center.x();
    v(1) = in.y() - center.y();
    v(2) = in.z() - center.z();
    v = v / v.norm();

    Eigen::Matrix3d p, q, r;
    p = normal * normal.transpose();
    q(0,0) = 0; q(0,1) = -normal.z(); q(0,2) = normal.y();
    q(1,0) = normal.z(); q(1,1) = 0; q(1,2) = -normal.x();
    q(2,0) = -normal.y(); q(2,1) = normal.x(); q(2,2) = 0;

    FT c = std::cos(angle);
    FT s = std::sin(angle);

    r = p + c*(Eigen::Matrix3d::Identity() - p) + s*q;
    Vector3d out = r * v;

    Point_3 new_point(center.x() + length * out(0),
                      center.y() + length * out(1),
                      center.z() + length * out(2));

    // debug stuff
    CGAL_precondition(CGAL::abs(v.norm()) > 1e-5);
    Eigen::Matrix3d sq_r = r*r.transpose();
    CGAL_postcondition(CGAL::abs(sq_r.norm() - CGAL::sqrt(3.)) < 1e-5);
    CGAL_postcondition(CGAL::abs(r.determinant() - 1.) < 1e-5);

    Vector3d new_vec(new_point.x() - center.x(),
                     new_point.y() - center.y(),
                     new_point.z() - center.z());
//    FT new_dot = CGAL::abs(new_vec.dot(normal));
//    CGAL_postcondition(new_dot < 1e-5);
    FT length_check = CGAL::abs(new_vec.norm() - length);
    CGAL_postcondition(length_check < 1e-5);
    new_vec = new_vec / new_vec.norm();
    FT new_cos = new_vec.dot(v);

    Vector3d new_normal = v.cross(new_vec);
    FT new_sin = new_normal.norm();
    new_sin *= (new_normal.dot(normal) > 0) ? 1 : -1;

    CGAL_postcondition(CGAL::abs(c - new_cos) < 1e-5);
    CGAL_postcondition(CGAL::abs(s - new_sin) < 1e-5);

    return new_point;
  }

  Eigen::Matrix3d rotate_metric(const Eigen::Matrix3d& f,
                                const Vector3d& edge_normal,
                                const Vector3d& unfolding_normal) const
  {
    // rotate the metric from one plane to another

    // todo
    // I'm sure there's something way smarter than expliciting the rotation,
    // desconstructing the metric, applying the rotation and rebuilding the metric
    // (and it's costly)

    CGAL_precondition(std::abs(edge_normal.norm() - 1.) < 1e-10 &&
                      std::abs(unfolding_normal.norm() - 1.) < 1e-10);

    // get the rotation axis (note that the directions do not matter)
    Vector3d rot_axis = edge_normal.cross(unfolding_normal);
    FT ra_norm = rot_axis.norm();

    if(std::abs(ra_norm < 1e-5)) // parallel planes
      return f;

    rot_axis = rot_axis / ra_norm;

    FT c = edge_normal.dot(unfolding_normal);
    if(c > FT(1.)) c = 1.;
    if(c < FT(-1.)) c = -1.;

    FT theta = std::acos(c); // dihedral angle between the planes... Should atan2 be used instead ?
    FT s = std::sin(theta);

#if (VERBOSITY > 25)
    std::cout << "dihedral: " << theta << std::endl;
#endif

    const Vector3d& v = rot_axis;

    // Rodrigues formula to get the rotation matrix
    Eigen::Matrix3d sm = Eigen::Matrix3d::Zero(); // skew_matrix
    sm(0,0) = 0; sm(0,1) = -v(2); sm(0,2) = v(1);
    sm(1,0) = v(2); sm(1,1) = 0; sm(1,2) = -v(0);
    sm(2,0) = -v(1); sm(2,1) = v(0); sm(2,2) = 0;

    Eigen::Matrix3d rot_m = Eigen::Matrix3d::Identity(); // rot_m
    rot_m = rot_m + s*sm + (1-c)*sm*sm;

#if (VERBOSITY > 25)
    std::cout << "rotation matrix : " << std::endl << rot_m << std::endl;
#endif

    // some debug
    Eigen::Matrix3d m_check = rot_m*rot_m.transpose();
    CGAL_postcondition((m_check - Eigen::Matrix3d::Identity()).norm() < 1e-5);
    Vector3d rotated_edge_normal = rot_m * edge_normal;
    CGAL_postcondition((rotated_edge_normal.cross(unfolding_normal)).norm() < 1e-5 &&
                        rotated_edge_normal.dot(unfolding_normal) > 0);
    // end debug

    // the metric if M = F^TF with F = V^T D V
    // the rotation is given by F_rot = (RV)^T D (RV), R being the rotation...
    typename Gt::Vector_3 v0, v1, v2;
    FT e0, e1, e2;
    get_eigen_vecs_and_vals<K>(f, v0, v1, v2, e0, e1, e2);

    v0 = transform<typename Gt::Vector_3>(rot_m, v0);
    v1 = transform<typename Gt::Vector_3>(rot_m, v1);
    v2 = transform<typename Gt::Vector_3>(rot_m, v2);

    return build_UDUt<K>(v0, v1, v2, e0, e1, e2);
  }

  template<typename Array_k, typename Array_kp1>
  void output_ancestor_edge(std::size_t n,
                            const Array_kp1& folded_points,
                            const Array_kp1& unfolded_points,
                            const Array_k& metrics) const
  {
    std::ofstream out("unfolding_visualization.mesh");

    CGAL_precondition(n < folded_points.size());

    // outputs :
    // - the folded edge (the path on the initial mesh)
    // - the unfolded edge (edges moved to a common plane while preserving angles
    //   and lengths)
    // - the rotated metrics

    out << "MeshVersionFormatted 1" << std::endl;
    out << "Dimension 3" << std::endl;
    out << "Vertices" << std::endl;
    out << 2.*n /*edges*/ + 4*(n-1) /*metrics*/ << std::endl;

    for(std::size_t i=0; i<n; ++i)
      out << folded_points[i] << " " << i << std::endl;

    for(std::size_t i=0; i<n; ++i)
      out << unfolded_points[i] << " " << i << std::endl;

    for(std::size_t i=0; i<n-1; ++i)
    {
      // metric on the edge unfolded[i] unfolded[i+1]
      const Eigen::Matrix3d& m = metrics[i];
      const Point_3& mid_point = CGAL::midpoint(unfolded_points[i], unfolded_points[i+1]);

      FT e0, e1, e2;
      Vector_3 v0, v1, v2;
      get_eigen_vecs_and_vals<K>(m, v0, v1, v2, e0, e1, e2);

      // note that there's no sqrt() since this is F (the square root of M) and not M
      e1 = 1./std::abs(e1);
      e2 = 1./std::abs(e2);
      e0 = 1./std::abs(e0);

      FT scaling = 0.1; // something pretty todo

      out << mid_point << " " << i << std::endl;
      out << mid_point + scaling * e0 * v0 << " " << i << std::endl;
      out << mid_point + scaling * e1 * v1 << " " << i << std::endl;
      out << mid_point + scaling * e2 * v2 << " " << i << std::endl;
    }

    out << "Edges" << std::endl;
    out << 2.*(n-1) /*edges*/ + 3*(n-1) /*metrics*/ << std::endl;

    for(std::size_t i=1; i<n; ++i)
    {
      out << i << " " << i+1 << " 1" << std::endl; // folded
      out << n+i << " " << n+i+1 << " 2" << std::endl; // unfolded
    }

    for(std::size_t i=2*n+1; i<=2*n+4*(n-1);)
    {
      out << i << " " << i + 1 << " 3" << std::endl;
      out << i << " " << i + 2 << " 4" << std::endl;
      out << i << " " << i + 3 << " 5" << std::endl;
      i += 4;
    }

    out << "End" << std::endl;
  }

  // this function is the heart of the painter
  bool compute_closest_seed(const Self& anc)
  {
    // returns true if we improved the distance

#if (VERBOSITY > 20)
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "compute closest seed for : " << this->index();
    std::cout << " (" << this->point  () << ") ";
    std::cout << "curr. dist: " << this->distance_to_closest_seed();
    std::cout << " anc: " << anc.index() << " & ancdist: ";
    std::cout << anc.distance_to_closest_seed() << std::endl;
#endif
    CGAL_assertion(anc.state() == KNOWN);

    const int k = 8; // depth of the ancestor edge
    FT d = FT_inf;

    // the path is 'this' to ancestor1, ancestor1 to ancestor2, etc.
    // stored as 'this', ancestor1, ancestor2, etc.
    boost::array<std::size_t, k+1> ancestor_path;
    for(int i=2; i<k+1; ++i)
      ancestor_path[i] = -1;
    ancestor_path[0] = this->index();
    ancestor_path[1] = anc.index();

    // compute the metric at the edges
    boost::array<Eigen::Matrix3d, k> path_metrics; // THESE ARE 'F', NOT 'M' !!!
    for(int i=0; i<k; ++i)
    {
      if(i >= 1)
      {
        if(this->canvas()->get_point(ancestor_path[i]).ancestor() ==
                                                   static_cast<std::size_t>(-1))
          break;

        ancestor_path[i+1] = this->canvas()->get_point(ancestor_path[i]).ancestor();
      }

      std::size_t e0 = ancestor_path[i];
      std::size_t e1 = ancestor_path[i+1];

      const Metric& m0 = this->canvas()->get_point(e0).metric();
      const Metric& m1 = this->canvas()->get_point(e1).metric();

      path_metrics[i] = get_interpolated_transformation(m0, m1);
    }

    Vector3d unfolded_edge = Vector3d::Zero(); // between 'this' and the i-th unfolded ancestor
    Vector3d unfolding_plane_normal;
    boost::array<Point_3, k+1> folded_points; // only needed for debug
    boost::array<Point_3, k+1> unfolded_points;
    boost::array<Vector3d, k> edge_segments;
    boost::array<Vector3d, k> unfolded_edge_segments;
    boost::array<Eigen::Matrix3d, k> rotated_path_metrics;
    std::size_t next_id = this->index(); // the descendant (opposite of ancestor)
    std::size_t curr_id = this->index();
    std::size_t prev_id = anc.index();

    unfolded_points[0] = this->point();
    folded_points[0] =  this->point();
    unfolded_points[1] = anc.point();

    // set up the unfolding plane's normal
    bool found = false;
    Facet_handle fh = this->incident_facets_in_complex_begin();
    Facet_handle fend = this->incident_facets_in_complex_end();
    for(; fh!=fend; ++fh)
    {
      Cell_handle ch = fh->first;
      std::size_t second = fh->second;

      if(!c3t3().is_in_complex(ch))
      {
        ch = ch->neighbor(fh->second);
        second = ch->index(fh->first);
      }

#if (VERBOSITY > 30)
      std::cout << "consider facet :" << std::endl;
      std::cout << ch->vertex((second+1)%4)->info() << std::endl;
      std::cout << ch->vertex((second+2)%4)->info() << std::endl;
      std::cout << ch->vertex((second+3)%4)->info() << std::endl;
      std::cout << "trying to find : " << m_v->info() << " and " << anc.m_v->info() << std::endl;
#endif

      Facet f(ch, second);
      CGAL_postcondition(c3t3().is_in_complex(f.first));

      int useless;
      if(c3t3().triangulation().has_vertex(f, m_v, useless) &&
         c3t3().triangulation().has_vertex(f, anc.m_v, useless))
      {
        found = true;
        CGAL_assertion(this->canvas()->m_facet_normals.find(f) !=
                       this->canvas()->m_facet_normals.end());
        unfolding_plane_normal = this->canvas()->m_facet_normals[f];
        CGAL_postcondition(std::abs(unfolding_plane_normal.norm() - 1.) < 1e-5);
        break;
      }
    }
    CGAL_postcondition(found);

    std::size_t best_i = -1;
    for(int i=1; i<=k; ++i)
    {
#if (VERBOSITY > 20)
      std::cout << "~~~~~~ depth i: " << i << std::endl;
#endif
      const Base& next_p = this->canvas()->get_point(next_id);
      const Base& curr_p = this->canvas()->get_point(curr_id);
      const Base& prev_p = this->canvas()->get_point(prev_id);

      folded_points[i] = prev_p.point();

#if (VERBOSITY > 25)
      std::cout << "current triplet on the ancestor path : "
                << next_id << " " << curr_id << " " << prev_id << std::endl;
      std::cout << next_p.point() << std::endl;
      std::cout << curr_p.point() << std::endl;
      std::cout << prev_p.point() << std::endl;
#endif

      Vector3d edge_segment;
      edge_segment(0) = prev_p.point().x() - curr_p.point().x();
      edge_segment(1) = prev_p.point().y() - curr_p.point().y();
      edge_segment(2) = prev_p.point().z() - curr_p.point().z();
      FT edge_segment_length = edge_segment.norm();
      CGAL_assertion(edge_segment_length > 1e-16);
      edge_segments[i-1] = edge_segment;

      // the next point to be unfolded is at distance edge_segment_length from
      // the previous unfolded point and with an angle
      if(i > 1)
      {
        FT angle = compute_flattening_angle(next_id, curr_id, prev_id);

        unfolded_points[i] = compute_unfolded_coordinates(unfolded_points[i-2],
                                                          unfolded_points[i-1],
                                                          unfolding_plane_normal,
                                                          edge_segment_length,
                                                          angle);
      }

      unfolded_edge_segments[i-1](0) = unfolded_points[i].x() - unfolded_points[i-1].x();
      unfolded_edge_segments[i-1](1) = unfolded_points[i].y() - unfolded_points[i-1].y();
      unfolded_edge_segments[i-1](2) = unfolded_points[i].z() - unfolded_points[i-1].z();
      unfolded_edge += unfolded_edge_segments[i-1];

#if (VERBOSITY > 25)
      std::cout << "unfolding normal : " << unfolding_plane_normal.transpose() << std::endl;
      std::cout << "unfolded " << prev_id << " (" << folded_points[i] << ") to : "
                << unfolded_points[i] << std::endl;
      std::cout << "compare distance from " << this->index() << " to " << prev_id << " "
                << CGAL::sqrt(CGAL::squared_distance(this->point(), folded_points[i])) << " "
                << CGAL::sqrt(CGAL::squared_distance(this->point(), unfolded_points[i])) << std::endl;

      std::cout << "compare distance from " << curr_id << " to " << prev_id << " "
                << CGAL::sqrt(CGAL::squared_distance(folded_points[i-1], folded_points[i])) << " "
                << CGAL::sqrt(CGAL::squared_distance(unfolded_points[i-1], unfolded_points[i])) << std::endl;

      if(i > 1)
      {
        std::cout << "check the angles (folded unfolded)" << std::endl;
        Vector3d folded_v_in, folded_v_out, unfolded_v_in, unfolded_v_out;

        folded_v_in(0) = folded_points[i-2].x() - folded_points[i-1].x();
        folded_v_in(1) = folded_points[i-2].y() - folded_points[i-1].y();
        folded_v_in(2) = folded_points[i-2].z() - folded_points[i-1].z();
        folded_v_in = folded_v_in / folded_v_in.norm();

        folded_v_out(0) = folded_points[i].x() - folded_points[i-1].x();
        folded_v_out(1) = folded_points[i].y() - folded_points[i-1].y();
        folded_v_out(2) = folded_points[i].z() - folded_points[i-1].z();
        folded_v_out = folded_v_out / folded_v_out.norm();

        unfolded_v_in(0) = unfolded_points[i-2].x() - unfolded_points[i-1].x();
        unfolded_v_in(1) = unfolded_points[i-2].y() - unfolded_points[i-1].y();
        unfolded_v_in(2) = unfolded_points[i-2].z() - unfolded_points[i-1].z();
        unfolded_v_in = unfolded_v_in / unfolded_v_in.norm();

        unfolded_v_out(0) = unfolded_points[i].x() - unfolded_points[i-1].x();
        unfolded_v_out(1) = unfolded_points[i].y() - unfolded_points[i-1].y();
        unfolded_v_out(2) = unfolded_points[i].z() - unfolded_points[i-1].z();
        unfolded_v_out = unfolded_v_out / unfolded_v_out.norm();

        std::cout << folded_v_in.transpose() << std::endl;
        std::cout << folded_v_out.transpose() << std::endl;
        std::cout << unfolded_v_in.transpose() << std::endl;
        std::cout << unfolded_v_out.transpose() << std::endl;

        std::cout << folded_v_in.dot(folded_v_out) << " " << unfolded_v_in.dot(unfolded_v_out) << std::endl; // that's the cos
      }
#endif

      FT ancestor_edge_length = unfolded_edge.norm();
      Vector3d normalized_unfolded_edge = unfolded_edge / ancestor_edge_length;

      // we must rotate the metric of the folded edge so that it is appropriate
      // for the unfolded edge
      Vector3d folded_edge_segment_normal;

      std::size_t i1, i2;
      if(prev_id < curr_id)
      {
        i1 = prev_id;
        i2 = curr_id;
      }
      else
      {
        i1 = curr_id;
        i2 = prev_id;
      }

      std::pair<std::size_t, std::size_t> e(i1, i2);
      CGAL_assertion(this->canvas()->m_edge_normals.find(e) !=
                     this->canvas()->m_edge_normals.end());

      folded_edge_segment_normal = this->canvas()->m_edge_normals[e];

      rotated_path_metrics[i-1] = rotate_metric(path_metrics[i-1],
                                                folded_edge_segment_normal,
                                                unfolding_plane_normal);

/*
 * Below is wrong: while this rotation will send the folded edge on the unfolded_edge,
 * it does not provide the correct rotation for the metric.
 * Leaving it here commented so the future generations can learn from the mistakes
 * of the past ones.
      rotated_path_metrics[i-1] = rotate_metric(path_metrics[i-1],
                                                edge_segments[i-1],
                                                unfolded_edge_segments[i-1]);
*/

#if (VERBOSITY > 25)
      std::cout << "folded edge segment: " << edge_segment.transpose() << std::endl;
      std::cout << "normalized unfolded (full) edge: " << normalized_unfolded_edge.transpose() << std::endl;
      std::cout << "norm of the unfolded (full) edge: " << ancestor_edge_length << std::endl;
#endif

      // compute the distance for the current depth (i) by splitting the unfolded
      // edge in segments.
      // The metric for each segment is drawn from the metric of the folded edge
      FT dist_to_ancestor = 0.;
      for(int j=0; j<i; ++j)
      {
#if (VERBOSITY > 20)
        std::cout << "adding part of depth: " << j+1 << " out of " << i << std::endl;
#endif
        const Vector3d& unfolded_edge_segment = unfolded_edge_segments[j];
        const Eigen::Matrix3d& f = rotated_path_metrics[j];

        Vector3d transformed_unfolded_edge = f * normalized_unfolded_edge;
        FT sp = unfolded_edge_segment.dot(normalized_unfolded_edge);

        // length of the normalized unfolded edge in the metric of the edge segment
        FT l = transformed_unfolded_edge.norm();
        dist_to_ancestor += sp * l;

#if (VERBOSITY > 30)
        std::cout << "unfolded edge segment: " << unfolded_edge_segment.transpose() << std::endl;
        std::cout << "rotated transformation :" << std::endl << f << std::endl;
        std::cout << "current transformed (unfolded full) edge: " << transformed_unfolded_edge.transpose() << std::endl;
        std::cout << "dist_to_anc: " << dist_to_ancestor << " sp: " << sp << " l: " << l << std::endl << std::endl;
#endif
      }

      dist_to_ancestor = (std::max)(dist_to_ancestor, 0.);

      // add ancestor edge length to the distance at that ancestor
      FT dist_at_anc = prev_p.distance_to_closest_seed();
      FT new_d = dist_at_anc + dist_to_ancestor;

#if (VERBOSITY > 25)
      std::cout << "potential update: " << new_d << " (old: " << dist_at_anc
                << " +new: " << dist_to_ancestor << ")" << std::endl;
#endif

      if(new_d < d)
      {
#if (VERBOSITY > 25)
        std::cout << "better length at depth " << i << std::endl;
#endif
        d = new_d;
        best_i = i;
      }

      // checks if we can go farther up in the ancestor path
      if(prev_p.ancestor() == static_cast<std::size_t>(-1) || i == k)
      {
        break;
      }

      next_id = curr_id;
      curr_id = prev_id;
      prev_id = prev_p.ancestor();
    }

#if (VERBOSITY > 20)
    std::cout << "distance with that anc: " << d << std::endl;
#endif

//    if(d < this->distance_to_closest_seed())
    if(this->distance_to_closest_seed() - d > 1e-8) // to avoid useless updates
    {
#if (VERBOSITY > 20)
      std::cout << "improving distance at " << Base::index() << " "
                << " from: " << this->distance_to_closest_seed()
                << " to " << d << std::endl;
      if(this->ancestor())
        std::cout << "prev anc: " << this->ancestor() << std::endl;
      else
        std::cout << "no prev anc" << std::endl;
      std::cout << "new anc: " << anc.index() << std::endl;
#endif

      // remove 'index' from the previous ancestor's children (if needed)
      if(this->ancestor() != static_cast<std::size_t>(-1))
        this->canvas()->get_point(this->ancestor()).remove_from_children(this->index());

      this->depth() = best_i;
      this->ancestor() = anc.index();
      this->distance_to_closest_seed() = d;
      this->closest_seed_id() = anc.closest_seed_id();

      // add 'index' to the new ancestor's children
      anc.m_children.insert(this->index());

//      if(!this->children().empty())
//      {
//        std::cout << "dealing with the " << this->children().size()
//                  << " descendant(s) of " << this->index() << std::endl;
//        std::cout << "position : " << this->point() << std::endl;
//      }

      if(ignore_children)
        return true;

      while(!this->children().empty())
      {
        Self& cp = this->canvas()->get_point(*(this->children().begin()));
        cp.reset_descendants();
        CGAL_postcondition(cp.ancestor_path_length()); // checks for circular ancestry
      }

      return true;
    }
    return false;
  }

  PQ_state update_neighbors_distances(std::vector<Campen_canvas_point*>& trial_pq)
  {
    // consider all the neighbors of a KNOWN point and compute their distance to 'this'
#if (VERBOSITY > 15)
    std::cout << "update neighbors of " << Base::index() << std::endl;
#endif
    CGAL_assertion(this->state() == KNOWN);

    PQ_state pqs_ret = NOTHING_TO_DO;

    Vertex_handle_handle it = adjacent_vertices_in_complex_begin();
    Vertex_handle_handle end = adjacent_vertices_in_complex_end();
    for(; it!=end; ++it)
    {
      Vertex_handle v = *it;
      Self& cp = this->canvas()->get_point(v->info());
      if(cp.state() == KNOWN)
        continue;
      else if(cp.state() == TRIAL)
      {
        // note that we don't insert in trial_pq since it's already in
        if(cp.compute_closest_seed(*this))
          pqs_ret = REBUILD_TRIAL;
      }
      else // cp.state == FAR
      {
        CGAL_assertion(cp.state() == FAR);

        // note that cp.distance_to_closest_seed is not necessarily FT_inf here :
        // if we're refining, we've assigned FAR to all points after inserting a new
        // seed, therefore we must verify that compute_closest_seed is an update
        // before inserting it in the trial_queue
        if(cp.compute_closest_seed(*this))
        {
          CGAL_assertion(cp.distance_to_closest_seed() != FT_inf);
          cp.state() = TRIAL;
          trial_pq.push_back(&cp);
          std::push_heap(trial_pq.begin(), trial_pq.end(),
                         Canvas_point_comparer<Self>());
        }
      }
      CGAL_assertion(cp.distance_to_closest_seed() != FT_inf);
    }
    return pqs_ret;
  }

  Campen_canvas_point() : Base(), m_v(NULL) { }
  Campen_canvas_point(const Weighted_point_3& p, const std::size_t index,
                      Vertex_handle v, Canvas* canvas)
    :
      Base(p.point(), index, canvas),
      m_v(v),
      ignore_children(false),
      m_is_vertices_cache_dirty(true),
      m_is_facets_cache_dirty(true),
      m_adjacent_vertices_cache(),
      m_incident_facets_cache()
  { }
};

} // namespace Anisotropic_mesh_3
} // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_CAMPEN_C2T3_POINT_H

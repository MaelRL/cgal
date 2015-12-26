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
  typedef typename Base::Metric                                       Metric;

  typedef Eigen::Matrix<FT, 2, 1>                                     Vector2d;
  typedef typename Base::Vector3d                                     Vector3d;

  typedef typename Tr::Vertex_handle                                  Vertex_handle;
  typedef std::vector<Vertex_handle>                                  Vertex_handle_vector;
  typedef typename Vertex_handle_vector::iterator                     Vertex_handle_handle;
  typedef typename Tr::Facet                                          Facet;
  typedef std::vector<Facet>                                          Facet_vector;
  typedef typename Facet_vector::iterator                             Facet_handle;
  typedef typename Tr::Cell_handle                                    Cell_handle;

  Vertex_handle m_v; // corresponding vertex in the c3t3

  C3t3& c3t3() { return this->canvas()->m_c3t3; }
  const C3t3& c3t3() const { return this->canvas()->m_c3t3; }

  mutable bool m_is_vertices_cache_dirty, m_is_facets_cache_dirty;

  mutable Vertex_handle_vector m_adjacent_vertices_cache; // on the complex boundary
  mutable Facet_vector m_incident_facets_cache; // on the complex boundary

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
    return c3t3().triangulation().tds().incident_facets(m_v, facets, In_complex_filter(c3t3()));
  }

  void set_adjacent_vertices_on_complex_boundary() const
  {
    // We're going to need to have them oriented to compute the relative angles
    // that are needed to flatten a path along the surface, so we might as well
    // do it here !

    // We orient them so the normal points outward from the surface

    CGAL_assertion(m_v != Vertex_handle());
    CGAL_assertion(c3t3().triangulation().dimension() > 1);

    // Build a map<Vh, Vh> (each entry is one oriented edge) and we link them together
    // by going through the map at the end to build the oriented ring
    boost::unordered_map<Vertex_handle, Vertex_handle> oriented_edges;

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

      // find the position of the center of the star in the cell
      std::size_t center_id = 0;
      for(; center_id<4; ++center_id)
      {
        if(m_v == c->vertex(center_id))
          break;
      }

      // todo: there might be something prettier here
      Vertex_handle vh_a, vh_b;
      for(std::size_t i=0; i<4; ++i)
      {
        if(i == second || i == center_id)
          continue;
        if(vh_a == Vertex_handle())
          vh_a = c->vertex(i);
        else
          vh_b = c->vertex(i);
      }
      CGAL_postcondition(vh_a != Vertex_handle() && vh_b != Vertex_handle());

      // Now look at the orientation of the tet (a, b, center_id, second)
      CGAL::Orientation o = CGAL::orientation(vh_a->point().point(),
                                              vh_b->point().point(),
                                              m_v->point().point(),
                                              c->vertex(second)->point().point());
      if(o == CGAL::POSITIVE) // ((ab^ac),ad) > 0
        oriented_edges[vh_b] = vh_a; // since we want the normal to point outside !
      else
      {
        CGAL_assertion(o == CGAL::NEGATIVE); // ZERO case should not happen
        oriented_edges[vh_a] = vh_b;
      }
    }
    CGAL_postcondition(nb_of_inc_facets == oriented_edges.size());

    // We can now put the (oriented) vertices in cache
    typename boost::unordered_map<Vertex_handle, Vertex_handle>::iterator it = oriented_edges.begin();
    Vertex_handle init = it->first;

//    std::cout << "adjacency of : " << this->index() << std::endl;

    while(true)
    {
      Vertex_handle edge_beginning = it->first;
//      std::cout << it->first->info() << " ";

      m_adjacent_vertices_cache.push_back(edge_beginning);
      it = oriented_edges.find(it->second);
      CGAL_postcondition(it != oriented_edges.end());

      if(it->first == init)
        break;
    }
//    std::cout << std::endl;
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

    std::cout << "Computing unfolding of new point..." << std::endl;
    std::cout << "Previous two: " << in << " || " << center << std::endl;
    std::cout << "normal: " << normal.transpose() << std::endl;
    std::cout << "point at length: " << length << " and angle: " << angle << std::endl;

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
    CGAL_precondition(CGAL::abs(v.norm()) > 1e-10);
    Eigen::Matrix3d asd = r*r.transpose();
    CGAL_postcondition(CGAL::abs(asd.norm() - CGAL::sqrt(3.)) < 1e-10);
    CGAL_postcondition(CGAL::abs(r.determinant() - 1.) < 1e-10);

    Vector3d new_vec(new_point.x() - center.x(),
                     new_point.y() - center.y(),
                     new_point.z() - center.z());
    CGAL_postcondition(CGAL::abs(new_vec.dot(normal)) < 1e-10);
    CGAL_postcondition(CGAL::abs(new_vec.norm() - length) < 1e-10);
    new_vec = new_vec / new_vec.norm();
    FT new_cos = new_vec.dot(v);

    Vector3d new_normal = v.cross(new_vec);
    FT new_sin = new_normal.norm();
    new_sin *= (new_normal.dot(normal) > 0) ? 1 : -1;

    CGAL_postcondition(CGAL::abs(c - new_cos) < 1e-10);
    CGAL_postcondition(CGAL::abs(s - new_sin) < 1e-10);

    return new_point;
  }

  Eigen::Matrix3d rotate_metric(const Eigen::Matrix3d& f,
                                const Vector3d& orig,
                                const Vector3d& end) const
  {
    // normalize the vectors
    Vector3d n_orig = orig / orig.norm();
    Vector3d n_end = end / end.norm();

    // check some degenerate (and easy) cases
    if(n_orig == n_end || n_orig == -n_end)
      return f;

    // Rodrigues formula to get the rotation matrix in a standard case
    Vector3d v = n_orig.cross(n_end);
    FT s = v.norm(); // that's the sin of the angle since the vectors are normalized
    FT c = n_orig.dot(n_end); // that's the cos of the angle since the vectors are normalized

    Eigen::Matrix3d sm = Eigen::Matrix3d::Zero(); // skew_matrix
    sm(0,0) = 0; sm(0,1) = -v(2); sm(0,2) = v(1);
    sm(1,0) = v(2); sm(1,1) = 0; sm(1,2) = -v(0);
    sm(2,0) = -v(1); sm(2,1) = v(0); sm(2,2) = 0;

    Eigen::Matrix3d rot_m = Eigen::Matrix3d::Identity(); // rot_m
    CGAL_precondition(s != 0.);
    rot_m = rot_m + sm + sm*sm*(1-c)/(s*s);

    CGAL_postcondition(CGAL::abs((rot_m*rot_m.transpose()).norm() - CGAL::sqrt(3.)) < 1e-10);
    CGAL_postcondition(CGAL::abs((rot_m * orig).dot(end) - end.dot(end)) < 1e-10);

    // the metric if F = V^T D V
    // the rotation is F_rot = (RV)^T D (RV), R being the rotation...
    typename Gt::Vector_3 v0, v1, v2;
    FT e0, e1, e2;
    get_eigen_vecs_and_vals<K>(f, v0, v1, v2, e0, e1, e2);

    v0 = transform<typename Gt::Vector_3>(rot_m, v0);
    v1 = transform<typename Gt::Vector_3>(rot_m, v1);
    v2 = transform<typename Gt::Vector_3>(rot_m, v2);

    return build_UDUt<K>(v0, v1, v2, e0, e1, e2);
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
    boost::array<Eigen::Matrix3d, k> path_metrics;
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
    std::vector<Point_3> unfolded_points(k+1);
    std::vector<Vector3d> edge_segments(k);
    std::vector<Vector3d> unfolded_edge_segments(k);
    std::size_t next_id = this->index(); // the descendant (opposite of ancestor)
    std::size_t curr_id = this->index();
    std::size_t prev_id = anc.index();

    unfolded_points[0] = this->point();
    unfolded_points[1] = anc.point();

    for(int i=1; i<=k; ++i)
    {
      const Base& next_p = this->canvas()->get_point(next_id);
      const Base& curr_p = this->canvas()->get_point(curr_id);
      const Base& prev_p = this->canvas()->get_point(prev_id);

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

      if(i == 1)
      {
        // the unfolding plane could be taken as the normal
        // of the surface at the center, but is there a point ?

        // edge_segment = (a,b,c), (-c,-c,a+b) is orthogonal but is (0,0,0) if c=0
        // and b=-a so we need to test out that case. In that degenerate case,
        // edge_semgent = (a,-a,0) and we can take (a,a,0) as orthogonal vector
        if(edge_segment(2) == 0. && edge_segment(0) == -edge_segment(1))
        {
          unfolding_plane_normal(0) = edge_segment(0);
          unfolding_plane_normal(1) = edge_segment(0);
          unfolding_plane_normal(2) = 0;
        }
        else
        {
          unfolding_plane_normal(0) = -edge_segment(2);
          unfolding_plane_normal(1) = -edge_segment(2);
          unfolding_plane_normal(2) = edge_segment(0) + edge_segment(1);
        }

        CGAL_postcondition(CGAL::abs(unfolding_plane_normal.dot(edge_segment)) < 1e-10);
        unfolding_plane_normal = unfolding_plane_normal / unfolding_plane_normal.norm();
      }
      else // i > 1
      {
      // the next point to be unfolded is at distance edge_segment_length from
      // the previous unfolded point and with an angle
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

      FT ancestor_edge_length = unfolded_edge.norm();
      Vector3d normalized_unfolded_edge = unfolded_edge / ancestor_edge_length;

      // compute the distance for the current depth (i) by splitting the unfolded
      // edge in segments.
      // The metric for each segment is drawn from the metric of the folded edge
      // segments (but this metric needs to also be unfolded)
      FT dist_to_ancestor = 0.;
      for(int j=0; j<i; ++j)
      {
        const Vector3d& unfolded_edge_segment = unfolded_edge_segments[j];

        // interpolate between both metric and transform the normalized edge
        // then we have (transformed_edge).norm() = || e ||_M = sqrt(e^t M e).
        // Since we have unfolded the edge, the metric must be rotated
        // by the rotation from the folded edge segment to the unfolded edge segment
        const Eigen::Matrix3d& f = rotate_metric(path_metrics[j],
                                                 edge_segments[j], // folded
                                                 unfolded_edge_segment);

        Vector3d transformed_unfolded_edge = f * normalized_unfolded_edge;
        FT sp = unfolded_edge_segment.dot(normalized_unfolded_edge);
        // length of the normalized unfolded edge in the metric of the edge segment
        FT l = transformed_unfolded_edge.norm();
        dist_to_ancestor += sp * l;

#if (VERBOSITY > 30)
        std::cout << "unfolded to : " << unfolded_points[i] << std::endl;
        std::cout << "unfolding normal : " << unfolding_plane_normal.transpose() << std::endl;
        std::cout << "folded edge segment: " << edge_segment.transpose() << std::endl;
        std::cout << "unfolded edge segment: " << unfolded_edge_segment.transpose() << std::endl;
        std::cout << "normalized unfolded (full) edge: " << normalized_unfolded_edge.transpose() << std::endl;
        std::cout << "rotated tranfsormation metric:" << std::endl << f << std::endl << std::endl;
        std::cout << "transformed edge: " << transformed_unfolded_edge.transpose() << std::endl;
        std::cout << "dist_to_anc: " << dist_to_ancestor << " sp: " << sp << " l: " << l << std::endl;
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
        std::cout << "better length with a depth " << i << std::endl;
#endif
        d = new_d;
      }

      // checks if we can go farther up in the ancestor path
      if(prev_p.ancestor() == static_cast<std::size_t>(-1))
        break;

      next_id = curr_id;
      curr_id = prev_id;
      prev_id = prev_p.ancestor();
    }

#if (VERBOSITY > 20)
    std::cout << "distance with that anc: " << d << std::endl;
#endif

    if(d < this->distance_to_closest_seed())
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
      m_is_vertices_cache_dirty(true),
      m_is_facets_cache_dirty(true),
      m_adjacent_vertices_cache(),
      m_incident_facets_cache()
  { }
};

} // namespace Anisotropic_mesh_3
} // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_CAMPEN_C2T3_POINT_H

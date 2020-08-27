// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Nico Kruithof <Nico@nghk.nl>
//                 Mael Rouxel-Labbé

#ifndef CGAL_P2T2_DELAUNAY_TRIANGULATION_ON_SQUARE_FLAT_TORUS_2_H
#define CGAL_P2T2_DELAUNAY_TRIANGULATION_ON_SQUARE_FLAT_TORUS_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/Periodic_2_triangulation_2/Square_flat_torus_2/Triangulation_on_square_flat_torus_2.h>
#include <CGAL/Periodic_2_triangulation_2/Square_flat_torus_2/Triangulation_face_base_on_square_flat_torus_2.h>
#include <CGAL/Periodic_2_triangulation_vertex_base_2.h>
#include <CGAL/Periodic_2_triangulation_2/internal/Periodic_2_Delaunay_triangulation_remove_traits_2.h>

#include <CGAL/algorithm.h>
#include <CGAL/Default.h>
#include <CGAL/iterator.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_utils_2.h>

#include <array>
#include <iostream>
#include <iterator>
#include <set>
#include <utility>
#include <vector>

namespace CGAL {

template < class Gt_,
           class Tds_ = CGAL::Default >
class Delaunay_triangulation_on_square_flat_torus_2
  : public Triangulation_on_square_flat_torus_2<
             Gt_,
             typename Default::Get<Tds_,
               Triangulation_data_structure_2<
                 Periodic_2_triangulation_vertex_base_2<Gt_>,
                 Triangulation_face_base_on_square_flat_torus_2<Gt_> > >::type>
{
public:
  typedef Gt_                                                     Geom_traits;

  typedef typename Default::Get<Tds_,
                                Triangulation_data_structure_2<
                                  Periodic_2_triangulation_vertex_base_2<Gt_>,
                                  Triangulation_face_base_on_square_flat_torus_2<Gt_> > >::type
                                                                  Triangulation_data_structure;

private:
  typedef Geom_traits                                             Gt;
  typedef Triangulation_data_structure                            Tds;
  typedef Delaunay_triangulation_on_square_flat_torus_2<Gt, Tds>  Self;

public:
  typedef Triangulation_on_square_flat_torus_2<Gt, Tds>           Base;

  typedef typename Gt::Periodic_2_offset_2                        Offset;
  typedef typename Gt::Domain                                     Domain;
  typedef std::array<int, 2>                                      Covering_sheets;

  typedef typename Gt::FT                                         FT;
  typedef typename Gt::Point_2                                    Point;
  typedef typename Gt::Segment_2                                  Segment;
  typedef typename Gt::Triangle_2                                 Triangle;
  typedef typename Gt::Iso_rectangle_2                            Iso_rectangle;

  typedef std::pair<Point, Offset>                                Periodic_point;
  typedef std::array<std::pair<Point, Offset>, 2>                 Periodic_segment;
  typedef std::array<std::pair<Point, Offset>, 3>                 Periodic_triangle;
  typedef std::array<std::pair<Point, Offset>, 4>                 Periodic_tetrahedron;

  typedef typename Base::size_type                                size_type;
  typedef typename Base::Locate_type                              Locate_type;

  typedef typename Base::Vertex_handle                            Vertex_handle;
  typedef typename Base::Vertex_circulator                        Vertex_circulator;
  typedef typename Base::Vertex_iterator                          Vertex_iterator;
  typedef typename Base::Edge                                     Edge;
  typedef typename Base::Edge_circulator                          Edge_circulator;
  typedef typename Base::Edge_iterator                            Edge_iterator;
  typedef typename Base::Face_handle                              Face_handle;
  typedef typename Base::Face_circulator                          Face_circulator;
  typedef typename Base::Face_iterator                            Face_iterator;
  typedef typename Base::Finite_vertices_iterator                 Finite_vertices_iterator;
  typedef typename Base::Finite_edges_iterator                    Finite_edges_iterator;
  typedef typename Base::Finite_faces_iterator                    Finite_faces_iterator;
  typedef typename Base::All_faces_iterator                       All_faces_iterator;

  typedef typename Base::Periodic_segment_iterator                Periodic_segment_iterator;
  typedef typename Base::Periodic_triangle_iterator               Periodic_triangle_iterator;

  //Tag to distinguish Delaunay from regular triangulations
  typedef Tag_false                                               Weighted_tag;

  // Tag to distinguish periodic triangulations from others
  typedef Tag_true                                                Periodic_tag;

private:
  class Point_hider
  {
  public:
    template <class InputIterator>
    inline void set_vertices(InputIterator, InputIterator) const { }

    inline void hide_point(Face_handle, const Point &) { CGAL_triangulation_assertion(false); }
    inline void hide(Point&, Face_handle) const { CGAL_triangulation_assertion(false); }
    inline void do_hide(const Point &, Face_handle) const { CGAL_triangulation_assertion(false); }

    template <class Conflict_tester>
    inline void hide_points(Vertex_handle, const Conflict_tester &) { CGAL_triangulation_assertion(false); }

    template <class Tester>
    inline bool replace_vertex(const Point&, Vertex_handle, const Tester&) const { return true; }
    inline Vertex_handle replace_vertex(Face_handle f, int index, const Point &) { return f->vertex(index); }
    inline void reinsert_vertices(Vertex_handle) { }
  };

  class Conflict_tester
  {
    const Self* tr_ptr;
    Point p;
    mutable Offset o;

  public:
    /// Constructor
    Conflict_tester(const Self* tr_ptr) : tr_ptr(tr_ptr), p(Point()) { }
    Conflict_tester(const Point &pt, const Self* tr_ptr) : tr_ptr(tr_ptr), p(pt) { }

    /// returns true if the circumcircle of 'f' contains 'p'
    bool operator()(const Face_handle f, const Offset& off) const
    {
      return (tr_ptr->side_of_circle(f, p, tr_ptr->combine_offsets(o, off), true) == ON_BOUNDED_SIDE);
    }

//    bool operator()(const Face_handle f, const Point& pt, const Offset& off) const
//    {
//      return (tr_ptr->side_of_circle(f, pt, o + off, true) == ON_BOUNDED_SIDE);
//    }

    int compare_weight(Point, Point) const { return 0; }

    bool test_initial_face(Face_handle f, const Offset& off) const
    {
      if(!(operator()(f, off)))
        CGAL_triangulation_assertion(false);
      return true;
    }

    void set_point(const Point &_p) { p = _p; }
    const Point& point() const { return p; }
    void set_offset(const Offset& off) const { o = off; }
    const Offset& get_offset() const { return o; }
  };

  class Cover_manager
  {
    Self& tr;

  public:
    Cover_manager(Self& tr) : tr(tr) { }

    void create_initial_triangulation() { tr.create_initial_triangulation(); }

    template <class FaceIt>
    void insert_unsatisfying_elements(Vertex_handle v, const FaceIt begin, const FaceIt end) {
      tr.insert_faces_with_too_big_circumradius(v, begin, end);
    }

    template <class FaceIt>
    void delete_unsatisfying_elements(const FaceIt begin, const FaceIt end) {
      tr.delete_faces_with_too_big_circumradius(begin, end);
    }

    bool can_be_converted_to_1_sheet() const { return tr.can_be_converted_to_1_sheet(); }

    bool update_cover_data_during_management(Face_handle new_f,
                                             const std::vector<Face_handle>& new_faces,
                                             const bool abort_if_cover_change)
    {
      return tr.update_cover_data_during_management(new_f, new_faces, abort_if_cover_change);
    }
  };

  template <class TriangulationR2>
  struct Vertex_remover
  {
    typedef TriangulationR2                                          Triangulation_R2;
    typedef typename std::vector<Point>::iterator                    Hidden_points_iterator;

    typedef std::pair<Vertex_handle, Vertex_handle>                  Vertex_pair;
    typedef std::map<Vertex_pair, Edge>                              Vertex_pair_Edge_map;

    typedef typename Triangulation_R2::Triangulation_data_structure  TDSE;
    typedef typename Triangulation_R2::Face_handle                   FaceE_handle;
    typedef typename Triangulation_R2::Vertex_handle                 VertexE_handle;
    typedef typename Triangulation_R2::Edge                          EdgeE;
    typedef typename Triangulation_R2::Finite_faces_iterator         Finite_facesE_iterator;

    typedef std::pair<VertexE_handle, VertexE_handle>                VertexE_pair;
    typedef std::map<Vertex_pair, EdgeE>                             Vertex_pair_EdgeE_map;

    typedef typename Vertex_pair_EdgeE_map::iterator                 Vertex_pair_EdgeE_map_it;

    Vertex_remover(const Self* t, Triangulation_R2& tmp_) : _t(t), tmp(tmp_) { }

    const Self* _t;
    Triangulation_R2& tmp;

    void add_hidden_points(Face_handle) { }
    Hidden_points_iterator hidden_points_begin() { return hidden.begin(); }
    Hidden_points_iterator hidden_points_end() { return hidden.end(); }

    std::vector<Point> hidden;
  };

private:
  struct Face_handle_hash
    : public CGAL::cpp98::unary_function<Face_handle, std::size_t>
  {
    std::size_t operator()(const Face_handle f) const
    {
      return boost::hash<typename Face_handle::pointer>()(&*f);
    }
  };

  typedef std::unordered_set<Face_handle, Face_handle_hash>       Too_big_circumdisks_set;
  typedef typename Too_big_circumdisks_set::const_iterator        Too_big_circumdisks_set_it;

public:
#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using Base::cw;
  using Base::ccw;

  using Base::combine_offsets;
  using Base::get_offset;
  using Base::get_neighbor_offset;
  using Base::extract_offset;

  using Base::construct_point;
  using Base::periodic_point;
  using Base::orientation;
  using Base::point;
  using Base::construct_segment;

  using Base::dimension;
  using Base::domain;
  using Base::geom_traits;
  using Base::tds;
  using Base::number_of_vertices;
  using Base::edges_begin;
  using Base::edges_end;
  using Base::faces_begin;
  using Base::faces_end;

  using Base::locate;
  using Base::is_1_cover;
#endif

private:
  /// This threshold should be chosen such that if all Delaunay balls have a squared radius smaller than this,
  /// we can be sure that there are no self-edges anymore.
  FT squared_circumradius_threshold;

  /// This container stores all the faces whose circumdisk squared radius is larger
  /// than the treshold `squared_circumradius_threshold`.
  Too_big_circumdisks_set faces_with_too_big_circumradius;

public:
  /// \name Constructors
  Delaunay_triangulation_on_square_flat_torus_2(const Gt& gt)
    : Base(gt)
  {
    update_cover_data_after_setting_domain(gt.get_domain());
  }

  Delaunay_triangulation_on_square_flat_torus_2(const Domain& domain = Domain())
    : Delaunay_triangulation_on_square_flat_torus_2(Gt(domain))
  { }

  /// Constructor with a range of points
  template <class InputIterator>
  Delaunay_triangulation_on_square_flat_torus_2(InputIterator first, InputIterator last,
                                                const Gt& gt = Gt())
    : Delaunay_triangulation_on_square_flat_torus_2(gt)
  {
    insert(first, last);
  }

  /// Copy
  // @todo (can't be "= default" because some members are pointers)
  Delaunay_triangulation_on_square_flat_torus_2(const Delaunay_triangulation_on_square_flat_torus_2<Gt, Tds>& tr) = delete;
  Delaunay_triangulation_on_square_flat_torus_2& operator=(const Delaunay_triangulation_on_square_flat_torus_2&) = delete;
  void swap(Delaunay_triangulation_on_square_flat_torus_2& tr) = delete;

  void clear()
  {
    Base::clear();
    clear_covering_data();
  }

  /// Domain setting
public:
  template <typename D> // that's domain, but need an abstract level to compile
  void update_cover_data_after_setting_domain(const D& d)
  {
    squared_circumradius_threshold = 0.0625 * d.systole_sq_length(); // case of Lattice
  }

  void update_cover_data_after_setting_domain(const Iso_rectangle& d)
  {
    // the criterion is that the largest circumdisk must have a diameter smaller than c/2
    // (c being the square side), thus we need a squared circumdisk radius smaller than c^2/16
    squared_circumradius_threshold = 0.0625 * square(d.xmax() - d.xmin());
  }

  void set_domain(const Domain& domain) override
  {
    clear();
    Base::set_domain(domain);
    update_cover_data_after_setting_domain(geom_traits().get_domain());
  }

public:
  /// Determines whether the point p lies on the (un-)bounded side of
  /// the circle through the points p0, p1 and p2
  Oriented_side side_of_oriented_circle(const Point& p0, const Point& p1, const Point& p2, const Point& p,
                                        bool perturb) const
  {
    Oriented_side os = geom_traits().side_of_oriented_circle_2_object()(p0, p1, p2, p);
    if((os != ON_ORIENTED_BOUNDARY) || (!perturb))
      return os;

    // We are now in a degenerate case => we do a symbolic perturbation.

    // We sort the points lexicographically.
    const Point * points[4] = { &p0, &p1, &p2, &p };
    std::sort(points, points + 4, typename Base::Perturbation_order(this));

    // We successively look whether the leading monomial, then 2nd monomial
    // of the determinant has non null coefficient.
    // 2 iterations are enough (cf paper)
    for(int i=3; i>0; --i)
    {
      if(points[i] == &p)
        return ON_NEGATIVE_SIDE; // since p0 p1 p2 are non collinear and positively oriented

      Orientation o;
      if(points[i] == &p2 && (o = orientation(p0, p1, p)) != COLLINEAR)
        return Oriented_side(o);
      if(points[i] == &p1 && (o = orientation(p0, p, p2)) != COLLINEAR)
        return Oriented_side(o);
      if(points[i] == &p0 && (o = orientation(p, p1, p2)) != COLLINEAR)
        return Oriented_side(o);
    }

    CGAL_triangulation_assertion(false);
    return ON_NEGATIVE_SIDE;
  }

  /// Determines whether the point (p,o) lies on the (un-)bounded side of
  /// the circle through the points (p0,o0), (p1,o1) and (p2,o2)
  Oriented_side side_of_oriented_circle(const Point& p0, const Point& p1, const Point& p2, const Point& p,
                                        const Offset& o0, const Offset& o1, const Offset& o2, const Offset& o,
                                        const bool perturb) const
  {
    Oriented_side os = geom_traits().side_of_oriented_circle_2_object()(p0, p1, p2, p, o0, o1, o2, o);
    if((os != ON_ORIENTED_BOUNDARY) || (!perturb))
      return os;

    // We are now in a degenerate case => we do a symbolic perturbation.
    // We sort the points lexicographically.
    Periodic_point pts[4] = { std::make_pair(p0, o0), std::make_pair(p1, o1),
                              std::make_pair(p2, o2), std::make_pair(p, o) };
    const Periodic_point *points[4] = { &pts[0], &pts[1], &pts[2], &pts[3] };

    std::sort(points, points + 4, typename Base::Perturbation_order(this));

    // We successively look whether the leading monomial, then 2nd monomial
    // of the determinant has non null coefficient.
    // 2 iterations are enough (cf paper)
    for(int i = 3; i > 0; --i)
    {
      if(points[i] == &pts[3])
        return ON_NEGATIVE_SIDE; // since p0 p1 p2 are non collinear and positively oriented

      Orientation orient;
      if((points[i] == &pts[2]) && ((orient = orientation(p0, p1, p, o0, o1, o)) != COLLINEAR))
        return Oriented_side(orient);
      if((points[i] == &pts[1]) && ((orient = orientation(p0, p, p2, o0, o, o2)) != COLLINEAR))
        return Oriented_side(orient);
      if((points[i] == &pts[0]) && ((orient = orientation(p, p1, p2, o, o1, o2)) != COLLINEAR))
        return Oriented_side(orient);
    }

    CGAL_triangulation_assertion(false);
    return ON_NEGATIVE_SIDE;
  }

  /// Determines whether the point p lies on the (un-)bounded side of
  /// the circle through the vertices of f
  Oriented_side side_of_oriented_circle(const Face_handle f, const Point& p,
                                        const bool perturb = false) const
  {
    Oriented_side os = ON_NEGATIVE_SIDE;

    int i = 0;
    // TODO: optimize which copies to check depending on the offsets in
    // the face.
    while(os == ON_NEGATIVE_SIDE && i < 4)
    {
      os = side_of_oriented_circle(f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point(), p,
                                   get_offset(f, 0), get_offset(f, 1), get_offset(f, 2), combine_offsets(Offset(), extract_offset(i)),
                                   perturb);
      ++i;
    }

    return os;
  }

  Bounded_side side_of_circle(const Face_handle f,
                              const Point& p, const Offset& o = Offset(),
                              const bool perturb = false) const
  {
#ifdef CGAL_DEBUG_P2T2
    std::cout << "side_of_circle" << std::endl;
    std::cout << "face: " << f->vertex(0)->point() << " (" << &*(f->vertex(0)) << ") | "
                                  << f->vertex(1)->point() << " (" << &*(f->vertex(0)) << ") | "
                                  << f->vertex(2)->point() << " (" << &*(f->vertex(0)) << ")" << std::endl;
    std::cout << "vertex offs: " << get_offset(f->vertex(0)) << " | "
                                 << get_offset(f->vertex(1)) << " | "
                                 << get_offset(f->vertex(2)) << std::endl;
    std::cout << "face offs: " << extract_offset(f->offset(0)) << " | "
                               << extract_offset(f->offset(1)) << " | "
                               << extract_offset(f->offset(2)) << std::endl;
    std::cout << "face points: " << point(f, 0) << " 0 | "
                                 << point(f, 1) << " 0 | "
                                 << point(f, 2) << " 0 | "
                                 << point(f, 0) << " 0" << std::endl;
    std::cout << "p: " << p << std::endl;
    std::cout << "o: " << o << std::endl;
    std::cout << "cp: " << construct_point(p, o) << std::endl;
#endif

    return enum_cast<Bounded_side>(
             side_of_oriented_circle(
               f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point(), p,
               this->get_offset(f, 0), this->get_offset(f, 1), this->get_offset(f, 2), o,
               perturb));
  }

  bool incircle(const int x, const int j, const int k, const int l,
                std::vector<Face_handle>&,
                std::vector<Vertex_handle>& w,
                std::vector<int>&)
  {
    return side_of_oriented_circle(w[j]->point(), w[k]->point(), w[l]->point(), w[x]->point(), true) ==  ON_POSITIVE_SIDE;
  }

  bool incircle(const int x, const int j, const int k, const int l,
                std::vector<Face_handle> &,
                std::vector<Vertex_handle> &w,
                std::vector<Offset> &o,
                std::vector<int> &)
  {
    return side_of_oriented_circle(w[j]->point(), w[k]->point(), w[l]->point(), w[x]->point(),
                                   o[j], o[k], o[l], o[x], true) ==  ON_POSITIVE_SIDE;
  }

  /// Checks whether f->vertex(i) lies outside the circumcircle of the face nb
  inline bool locally_Delaunay(Face_handle f, int i, Face_handle nb)
  {
    CGAL_BRANCH_PROFILER("locally_Delaunay(), simplicity check failures", tmp);

    bool simplicity_criterion = is_1_cover() && f->has_zero_offsets() && nb->has_zero_offsets();

    const Point *p[4];
    for(int index = 0; index < 3; ++index)
      p[index] = &nb->vertex(index)->point();

    p[3] = &f->vertex(i)->point();

    Oriented_side os;
    if(simplicity_criterion)
    {
      // No periodic offsets
      os = side_of_oriented_circle(*p[0], *p[1], *p[2], *p[3], true);
    }
    else
    {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);

      Offset off[4];

      for(int index=0; index<3; ++index)
        off[index] = get_offset(nb, index);

      off[3] = combine_offsets(get_offset(f, i), get_neighbor_offset(f, i));

      os = side_of_oriented_circle(*p[0], *p[1], *p[2], *p[3],
                                   off[0], off[1], off[2], off[3], true);
    }

    return (ON_POSITIVE_SIDE != os);
  }

  Comparison_result compare_squared_circumradius_to_threshold(const Periodic_point& p0,
                                                              const Periodic_point& p1,
                                                              const Periodic_point& p2,
                                                              const FT threshold) const
  {
//    std::cout << "threshold: " << squared_circumradius_threshold << std::endl;
//    std::cout << p0.first << " " << p0.second << " " << point(p0) << std::endl;
//    std::cout << p1.first << " " << p1.second << " " << point(p1) << std::endl;
//    std::cout << p2.first << " " << p2.second << " " << point(p2) << std::endl;

    return geom_traits().compare_squared_radius_2_object()(p0.first, p1.first, p2.first,
                                                           p0.second, p1.second, p2.second,
                                                           threshold);
  }

  Comparison_result compare_squared_circumradius_to_threshold(const Face_handle face, const FT threshold) const
  {
    Periodic_point p0 = periodic_point(face, 0);
    Periodic_point p1 = periodic_point(face, 1);
    Periodic_point p2 = periodic_point(face, 2);

    return compare_squared_circumradius_to_threshold(p0, p1, p2, threshold);
  }

  /// Constructs the circumcenter of the face f, respects the offset
  Point circumcenter(const Face_handle f) const
  {
    return construct_circumcenter(f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point(),
                                  get_offset(f, 0), get_offset(f, 1), get_offset(f, 2));
  }

  Point construct_circumcenter(const Point& p1, const Point& p2, const Point& p3,
                               const Offset& o1, const Offset& o2, const Offset& o3) const
  {
    return geom_traits().construct_circumcenter_2_object()(p1, p2, p3, o1, o2, o3);
  }

public:
  /// \name Conflict checking

  template <class OutputIteratorBoundaryEdges,
            class OutputIteratorFaces,
            class OutputIteratorInternalEdges>
  Triple<OutputIteratorBoundaryEdges, OutputIteratorFaces, OutputIteratorInternalEdges>
  find_conflicts(const Point& p,
                 Face_handle f,
                 OutputIteratorBoundaryEdges beit,
                 OutputIteratorFaces fit,
                 OutputIteratorInternalEdges ieit) const
  {
    CGAL_triangulation_precondition(number_of_vertices() != 0);

    std::vector<Edge> edges;
    edges.reserve(16);
    std::vector<Face_handle> faces;
    faces.reserve(16);

    Conflict_tester tester(p, this);

    Triple<typename std::back_insert_iterator<std::vector<Edge> >,
           typename std::back_insert_iterator<std::vector<Face_handle> >,
           OutputIteratorInternalEdges> tit =
        Base::find_conflicts(f, tester, make_triple(std::back_inserter(edges),
                                                    std::back_inserter(faces),
                                                    ieit));
    ieit = tit.third;

    // Reset the conflict flag on the boundary.
    for(const Edge& e : edges)
    {
      e.first->neighbor(e.second)->tds_data().clear();
      *beit++ = e;
    }

    // Reset the conflict flag in the conflict faces.
    for(Face_handle f : faces)
    {
      f->tds_data().clear();
      *fit++ = f;
    }

    for(Vertex_handle vc : this->v_offsets)
      vc->clear_offset();
    this->v_offsets.clear();

    return make_triple(beit, fit, ieit);
  }

  template <class OutputIteratorBoundaryEdges,
            class OutputIteratorFaces>
  std::pair<OutputIteratorBoundaryEdges, OutputIteratorFaces>
  find_conflicts(const Point &p,
                 Face_handle f,
                 OutputIteratorBoundaryEdges beit,
                 OutputIteratorFaces fit) const
  {
    Triple<OutputIteratorBoundaryEdges, OutputIteratorFaces, Emptyset_iterator> t =
      find_conflicts(p, f, beit, fit, Emptyset_iterator());

    return std::make_pair(t.first, t.second);
  }

  template <class OutputIteratorFaces>
  OutputIteratorFaces
  find_conflicts(const Point &p,
                 Face_handle f,
                 OutputIteratorFaces fit) const
  {
    Triple<Emptyset_iterator, OutputIteratorFaces, Emptyset_iterator> t =
      find_conflicts(p, f, Emptyset_iterator(), fit, Emptyset_iterator());

    return t.second;
  }

public:
  /// \name Insertion-Removal
  Vertex_handle insert(const Point& p,
                       Face_handle start = Face_handle())
  {
    Conflict_tester tester(p, this);
    Point_hider hider;
    Cover_manager cover_manager(*this);

    Vertex_handle v = Base::insert_in_conflict(p, start, tester, hider, cover_manager);
    CGAL_assertion(v != Vertex_handle());

    return v;
  }

  Vertex_handle insert(const Point& p,
                       Locate_type lt,
                       Face_handle f,
                       int li)
  {
    Conflict_tester tester(p, this);
    Point_hider hider;
    Cover_manager cover_manager(*this);

    Vertex_handle v = Base::insert_in_conflict(p, lt, f, li, tester, hider, cover_manager);
    CGAL_assertion(v != Vertex_handle());

    return v;
  }

  Vertex_handle push_back(const Point& p) { return insert(p); }

  /// Insertion with info
  template <class InputIterator>
  std::ptrdiff_t insert(InputIterator first, InputIterator last,
                        bool is_large_point_set = true)
  {
    if(first == last)
      return 0;

    size_type n = number_of_vertices();

    // The heuristic discards the existing triangulation so it can only be
    // applied to empty triangulations.
    if(n != 0)
      is_large_point_set = false;

    std::vector<Point> points(first, last);
    std::vector<Vertex_handle> dummy_points;
    typename std::vector<Point>::iterator pbegin = points.begin();

    if(is_large_point_set)
    {
      std::vector<Vertex_handle> dummy_points = this->insert_dummy_points();
    }
    else if(!is_1_cover())
    {
      CGAL::cpp98::random_shuffle(points.begin(), points.end());
      pbegin = points.begin();

      for(;;)
      {
        if(pbegin == points.end())
          return number_of_vertices() - n;

        insert(*pbegin);
        ++pbegin;

        if(is_1_cover())
          break;
      }
    }

    CGAL_triangulation_assertion(is_1_cover());

    // Organize the points to improve runtime
    spatial_sort(pbegin, points.end(), geom_traits());

    Face_handle hint;
    Conflict_tester tester(*pbegin, this);
    Point_hider hider;
    Cover_manager cover_manager(*this);

    // Actual insertion
    std::vector<Vertex_handle> double_vertices =
      Base::insert_in_conflict(points.begin(), points.end(), hint, tester, hider, cover_manager);

    CGAL_assertion_code(for(Vertex_handle v : double_vertices))
    CGAL_assertion(v != Vertex_handle());

    if(is_large_point_set)
    {
      typedef CGAL::Periodic_2_Delaunay_triangulation_remove_traits_2<Gt> P2removeT;
      typedef CGAL::Delaunay_triangulation_2<P2removeT> DT;
      typedef Vertex_remover<DT> Remover;

      P2removeT remove_traits(domain());
      DT dt(remove_traits);
      Remover remover(this, dt);
      Conflict_tester t(this);

      for(const Vertex_handle dummy_v : dummy_points)
      {
        if(std::find(double_vertices.begin(), double_vertices.end(), dummy_v) == double_vertices.end())
          Base::remove(dummy_v, remover, t, cover_manager);
      }
    }

    return number_of_vertices() - n;
  }

  void remove(Vertex_handle v)
  {
    typedef CGAL::Periodic_2_Delaunay_triangulation_remove_traits_2<Gt> P2removeT;
    typedef CGAL::Delaunay_triangulation_2<P2removeT> Euclidean_triangulation;
    typedef Vertex_remover<Euclidean_triangulation> Remover;

    P2removeT remove_traits(geom_traits());
    Euclidean_triangulation tmp(remove_traits);
    Remover remover(this, tmp);
    Conflict_tester ct(this);
    Cover_manager cover_manager(*this);

    Base::remove(v, remover, ct, cover_manager);
    CGAL_triangulation_expensive_assertion(is_valid());
  }

  template <typename InputIterator>
  std::ptrdiff_t remove(InputIterator first, InputIterator beyond)
  {
    std::size_t n = number_of_vertices();
    while(first != beyond)
      remove(*first++);

    return n - number_of_vertices();
  }

  /// \name Displacement

  Vertex_handle move_if_no_collision(Vertex_handle v, const Point& p)
  {
    Locate_type lt;
    int li;
    Vertex_handle inserted;
    Face_handle loc = locate(p, lt, li, v->face());

    if(lt == Base::VERTEX)
      return v;
    else
      /// This can be optimized by checking whether we can move v->point() to p
      return insert(p, lt, loc, li);
  }

  Vertex_handle move_point(Vertex_handle v, const Point& p)
  {
    if(v->point() == p)
      return v;

    Vertex_handle w = move_if_no_collision(v, p);
    if(w != v)
    {
      remove(v);
      return w;
    }
    return v;
  }

  /// \name Check - Query

  Vertex_handle nearest_vertex_2D(const Point& p, Face_handle f) const
  {
    CGAL_triangulation_precondition(dimension() == 2);
    f = locate(p, f);

    typename Geom_traits::Compare_distance_2 compare_distance = geom_traits().compare_distance_2_object();

    Vertex_handle nn = f->vertex(0);
    if(compare_distance(p, f->vertex(1)->point(), nn->point()) == SMALLER)
      nn = f->vertex(1);
    if(compare_distance(p, f->vertex(2)->point(), nn->point()) == SMALLER)
      nn = f->vertex(2);

    look_nearest_neighbor(p, f, 0, nn);
    look_nearest_neighbor(p, f, 1, nn);
    look_nearest_neighbor(p, f, 2, nn);

    return nn;
  }

  void look_nearest_neighbor(const Point& p,
                             Face_handle f,
                             int i,
                             Vertex_handle& nn) const
  {
    Face_handle  ni = f->neighbor(i);
    if(this->side_of_oriented_circle(ni, p, true) != ON_POSITIVE_SIDE)
      return;

    typename Geom_traits::Compare_distance_2 compare_distance = geom_traits().compare_distance_2_object();
    i = ni->index(f);
    if(compare_distance(p, ni->vertex(i)->point(), nn->point()) == SMALLER)
      nn = ni->vertex(i);

    // recursive exploration of triangles whose circumcircle contains p
    look_nearest_neighbor(p, ni, ccw(i), nn);
    look_nearest_neighbor(p, ni, cw(i), nn);
  }

  /// Returns the vertex closest to p, the point location will start from f
  Vertex_handle nearest_vertex(const Point& p, Face_handle f = Face_handle()) const
  {
    switch(dimension())
    {
      case -2:
        return Vertex_handle();
      case 2:
        return nearest_vertex_2D(p, f);
      default:
        CGAL_triangulation_assertion(false);
        break;
    }

    return Vertex_handle();
  }

public:
  /// \name Dual

  /// Returns the dual of f, which is the circumcenter of f.
  Point dual(Face_handle f) const
  {
    CGAL_triangulation_precondition(dimension() == 2);
    return circumcenter(f);
  }

  /// Returns the dual of e, which is always a segment in the periodic triangulation.
  Segment dual(const Edge& e) const
  {
    // dimension==2
    Face_handle nb = e.first->neighbor(e.second);
    Point p0 = dual(e.first);
    Point p1 = dual(nb);
    Offset o = combine_offsets(Offset(), get_neighbor_offset(e.first, e.second));
    Segment s = construct_segment(p0, p1, o, Offset());

    return s;
  }

  /// Returns the dual of the edge pointed to by ec.
  Segment dual(const Edge_circulator& ec) const { return dual(*ec); }
  /// Returns the dual of the edge pointed to by ei.
  Segment dual(const Edge_iterator& ei) const { return dual(*ei); }

  template <class Stream>
  Stream& draw_dual(Stream& ps)
  {
    Edge_iterator eit = edges_begin(), eend = edges_end();
    for(; eit!=eend; ++eit)
      ps << dual(eit);

    return ps;
  }

  /// \name Methods regarding the covering
  /// \{

protected:
  void create_initial_triangulation()
  {
    CGAL_triangulation_assertion(faces_with_too_big_circumradius.empty());

    for(Face_iterator iter=faces_begin(), end_iter=faces_end(); iter!=end_iter; ++iter)
    {
      if(compare_squared_circumradius_to_threshold(iter, squared_circumradius_threshold) != CGAL::SMALLER)
        faces_with_too_big_circumradius.insert(iter);
    }
  }

  template <class FaceIt>
  void insert_faces_with_too_big_circumradius(Vertex_handle /*v*/, FaceIt begin, const FaceIt end)
  {
    for(; begin != end; ++begin)
      if(compare_squared_circumradius_to_threshold(*begin, squared_circumradius_threshold) != CGAL::SMALLER)
        faces_with_too_big_circumradius.insert(*begin);
  }

  void insert_faces_with_too_big_circumradius(Face_iterator begin, Face_iterator end)
  {
    for(; begin != end; ++begin)
      if(compare_squared_circumradius_to_threshold(begin, squared_circumradius_threshold) != CGAL::SMALLER)
        faces_with_too_big_circumradius.insert(begin);
  }

  template <class FaceIt>
  void delete_faces_with_too_big_circumradius(FaceIt begin, const FaceIt end)
  {
    for(; begin != end; ++begin)
      faces_with_too_big_circumradius.erase(*begin);
  }

  // returns 'true/false' depending on whether the cover would (or has, if 'abort_if_cover_change'
  // is set to 'false') change.
  bool update_cover_data_during_management(Face_handle new_f,
                                           const std::vector<Face_handle>& new_faces,
                                           const bool abort_if_cover_change)
  {
    if(compare_squared_circumradius_to_threshold(new_f, squared_circumradius_threshold) != CGAL::SMALLER)
    {
      if(is_1_cover())
      {
        // Whether we are changing the cover or simply aborting, we need to get rid of the new faces
        tds().delete_faces(new_faces.begin(), new_faces.end());

        if(!abort_if_cover_change)
          this->convert_to_9_sheeted_covering();

        return true;
      }
      else
      {
        faces_with_too_big_circumradius.insert(new_f);
      }
    }

    return false;
  }

  virtual void update_cover_data_after_converting_to_9_sheeted_covering() override
  {
    for(Face_iterator iter = faces_begin(), end_iter = faces_end(); iter != end_iter; ++iter)
      if(compare_squared_circumradius_to_threshold(iter, squared_circumradius_threshold) != CGAL::SMALLER)
        faces_with_too_big_circumradius.insert(iter);
  }

  virtual void clear_covering_data() override
  {
    faces_with_too_big_circumradius.clear();
  }

public:
  bool can_be_converted_to_1_sheet() const { return faces_with_too_big_circumradius.empty(); }

  /// Checks whether the triangulation is a valid simplicial complex in the one cover.
  /// Uses an edge-length-criterion.
  bool is_extensible_triangulation_in_1_sheet_h1() const
  {
    if(!is_1_cover())
      return can_be_converted_to_1_sheet();

    return is_extensible_triangulation_in_1_sheet_h2();
  }

  /// Checks whether the triangulation is a valid simplicial complex in the one cover.
  /// Uses a criterion based on the maximal radius of the circumscribing circle.
  bool is_extensible_triangulation_in_1_sheet_h2() const
  {
    for(Periodic_triangle_iterator tit = this->periodic_triangles_begin(Base::UNIQUE);
        tit != this->periodic_triangles_end(Base::UNIQUE); ++tit)
    {
      if(compare_squared_circumradius_to_threshold(tit->at(0), tit->at(1), tit->at(2),
                                                   squared_circumradius_threshold) != CGAL::SMALLER)
        return false;
    }

    return true;
  }

public:
  /// \name Checking
  bool is_valid(bool verbose = false, int level = 0) const
  {
    // Check the parent
    bool result = Base::is_valid(verbose, level);

    // Check in_sphere:
    if(dimension() == 2)
    {
      const Point *p[4];
      Offset off[4];
      for(Face_iterator fit = faces_begin(); fit != this->faces_end(); ++fit)
      {
        for(int i=0; i<3; ++i)
        {
          p[i] = &fit->vertex(i)->point();
          off[i] = get_offset(fit, i);
        }

        /// Check whether the vertices of the neighbor lie outside the circumcircle of the face
        for(int i=0; i<3; ++i)
        {
          p[3] = &fit->vertex(i)->point();
          off[3] = combine_offsets(get_offset(fit, i), get_neighbor_offset(fit, i));

          result &= ON_POSITIVE_SIDE !=
              side_of_oriented_circle(*p[0], *p[1], *p[2], *p[3],
                                      off[0], off[1], off[2], off[3],
                                      false);
          CGAL_triangulation_assertion(result);
        }
      }
    }

    return result;
  }
};

} // namespace CGAL

#endif // CGAL_P2T2_DELAUNAY_TRIANGULATION_ON_SQUARE_FLAT_TORUS_2_H

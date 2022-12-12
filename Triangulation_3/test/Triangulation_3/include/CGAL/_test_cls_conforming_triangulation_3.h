#include "_test_cls_iterator.h"
#include "_test_cls_circulator.h"
#include <CGAL/Testsuite/Triangulation_23/test_move_semantic.h>
#include <CGAL/Testsuite/use.h>

#include <CGAL/Random.h>
#include <CGAL/use.h>

// - constraint iterator?
// - naming: insert() or insert_constraint() for constraints
// - clear_constraints() exists?
// - incident_constraints() accessor?
// - API to insert ranges of points? API to insert ranges of points & ranges of constraints?

// - @todo: move is_conforming() into is_valid()

enum Intersection_type
{
  NO_INTERSECTION = 0,
  INTERSECTION_WITHOUT_CONSTRUCTION,
  INTERSECTION
};

template <class Triang, class Pt>
void
_test_cdt_throwing(const Pt& p0, const Pt& p1, const Pt& p2, const Pt& p3,
                   const Intersection_type& intersection_type)
{
  std::cout << "test_cdt_throwing [" << p0 << "] - [" << p1 << "] || [" << p2 << "] - [" << p3 << "]" << std::endl;

  try
  {
    Triang tr;
    tr.insert_constraint(p0, p1);
    tr.insert_constraint(p2, p3);
  }
  catch (typename Triang::Intersection_of_constraints_exception& /*e*/)
  {
    std::cout << "threw, expected: " << intersection_type << std::endl;

    // There must have been an intersection
    assert(intersection_type != NO_INTERSECTION);

    // If the intersection requires no construction, then only 'No_constraint_intersection_tag' throws
    if(intersection_type == INTERSECTION_WITHOUT_CONSTRUCTION)
    {
      assert((std::is_same<typename Triang::Itag, CGAL::No_constraint_intersection_tag>::value));
    }
    else // threw and it's not a construction-less intersection ---> real intersection
    {
      assert(intersection_type == INTERSECTION);
#if !defined(CGAL_NO_DEPRECATED_CODE) && defined(CGAL_NO_DEPRECATION_WARNINGS)
      assert((std::is_same<typename Triang::Itag, CGAL::No_constraint_intersection_tag>::value) ||
             (std::is_same<typename Triang::Itag, CGAL::No_intersection_tag>::value) ||
             (std::is_same<typename Triang::Itag, CGAL::No_constraint_intersection_requiring_constructions_tag>::value));
#else
      assert((std::is_same<typename Triang::Itag, CGAL::No_constraint_intersection_tag>::value) ||
             (std::is_same<typename Triang::Itag, CGAL::No_constraint_intersection_requiring_constructions_tag>::value));
#endif
    }

    return;
  }

  if(intersection_type == INTERSECTION_WITHOUT_CONSTRUCTION)
  {
    // Even with an intersection without construction, 'No_constraint_intersection_tag' should have thrown
    assert(!(std::is_same<typename Triang::Itag, CGAL::No_constraint_intersection_tag>::value));
  }
  else if(intersection_type == INTERSECTION)
  {
    assert(!(std::is_same<typename Triang::Itag, CGAL::No_constraint_intersection_tag>::value) &&
           !(std::is_same<typename Triang::Itag, CGAL::No_constraint_intersection_requiring_constructions_tag>::value));
#if !defined(CGAL_NO_DEPRECATED_CODE) && defined(CGAL_NO_DEPRECATION_WARNINGS)
    assert(!(std::is_same<typename Triang::Itag, CGAL::No_intersection_tag>::value));
#endif
  }
}

template <class Triangulation>
void
_test_cls_conforming_delaunay_3(const Triangulation&)
{
  typedef Triangulation                      Cls;

  static_assert(std::is_nothrow_move_constructible<Cls>::value, "move cstr is missing");
  static_assert(std::is_nothrow_move_assignable<Cls>::value, "move assignment is missing");

  using Location_policy = typename Test_location_policy<Cls>::Location_policy;

  using Geom_traits = typename Cls::Geom_traits;
  using Triangulation_data_structure = typename Cls::Triangulation_data_structure;

  using Point = typename Cls::Point;
  using Segment = typename Cls::Segment;
  using Triangle = typename Cls::Triangle;
  using Tetrahedron = typename Cls::Tetrahedron;

  using Vertex = typename Cls::Vertex;
  using Cell = typename Cls::Cell;
  using Facet = typename Cls::Facet;
  using Edge = typename Cls::Edge;

  using size_type = typename Cls::size_type;

  using Vertex_handle = typename Cls::Vertex_handle;
  using Cell_handle = typename Cls::Cell_handle;
  using Vertex_iterator = typename Cls::Vertex_iterator;
  using Cell_iterator = typename Cls::Cell_iterator;
  using Locate_type = typename Cls::Locate_type;

  using Finite_vertices_iterator = typename Cls::Finite_vertices_iterator;
  using Finite_cells_iterator = typename Cls::Finite_cells_iterator;

  CGAL_USE_TYPE(Location_policy);
  CGAL_USE_TYPE(Segment);
  CGAL_USE_TYPE(Triangle);
  CGAL_USE_TYPE(Tetrahedron);
  CGAL_USE_TYPE(Vertex);
  CGAL_USE_TYPE(Cell);
  CGAL_USE_TYPE(Vertex_iterator);
  CGAL_USE_TYPE(Cell_iterator);

  CGAL_USE_TYPE(typename Cls::Constraint_hierarchy_tag);
  CGAL_USE_TYPE(typename Cls::Periodic_tag);
  CGAL_USE_TYPE(typename Cls::Weighted_tag);

  // +++ We define now some points and constraints for building triangulations +++ //
  int n = 100;
  CGAL::Random rnd;
  std::cout << "Seed = " << rnd.get_seed() << std::endl;

  using Point_vector = std::vector<Point>;

  using Point_pair_list = std::list<std::pair<Point, Point> >;
  using Segment_list = std::list<Segment>;

  // points
  Point_vector v1d;
  for(int i=0; i<10; ++i)
    v1d.emplace_back(0,0,i);

  Point_vector v2d;
  for(int a=0;a!=10;++a)
    for(int b=0;b!=10;++b)
      for(int d=0;d!=5;++d)
        v2d.emplace_back(a*b - d*a + (a-b)*10 + a,
                         a   - b+d + 5*b,
                        0);

  Point_vector v3d_1;
  for(int a=0;a!=10;++a)
    for(int b=0;b!=10;++b)
      for(int d=0;d!=5;++d)
        v3d_1.emplace_back(a*b - d*a + (a-b)*10 + a,
                           a   - b+d + 5*b,
                           a*a - d*d + b);

  Point_vector v3d_2;
  for(int a=0;a!=4;++a)
    for(int b=0;b!=4;++b)
      for(int d=0;d!=4;++d)
        v3d_2.emplace_back(a*b - d*a)*10 + a,
                           (a - b + d + 5*b)*100,
                           a*a - d*d - b);

  // constraints
  Point_pair_list empty_pp_constraints;
  Segment_list empty_seg_constraints;

  Point_pair_list pp_constraints {{ Point(0,0,0), Point(1,0,0) }};
  Segment_list seg_constraints_2 { {Point(0,0,0), Point(1,0,0)}, {Point(1,0,0), Point(1,0,1)} };
  Segment_list seg_constraints_3 { {Point(0,0,0), Point(1,0,0)},
                                   {Point(1,0,0), Point(1,0,1)},
                                   {Point(1,0,1), Point(1,1,1)} };

  Geom_traits gt;

  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  /// CONSTRUCTORS

  // Do not test basic constructors with points, those are already tested by _test_cls_delaunay_triangulation_3()

  // 0D
  // ---
  Cls T0(std::cbegin(empty_pp_constraints), std::cend(empty_pp_constraints));
  assert(T0.dimension() == -1);
  assert(T0.number_of_vertices() == 0);
  assert(T0.number_of_finite_cells() == 0);
  assert(T0.is_valid());

  // ---
  Cls T0b(std::cbegin(empty_pp_constraints), std::cend(empty_pp_constraints), gt);
  assert(T0b.dimension() == -1);
  assert(T0b.number_of_vertices() == 0);
  assert(T0b.number_of_finite_cells() == 0);
  assert(std::distance(T0b.conforming_edges_begin(), T0b.conforming_edges_end()) == 0);
  assert(T0b.is_valid());

  // ---
  Cls T0c(std::cbegin(empty_seg_constraints), std::cend(empty_seg_constraints));
  assert(T0c.dimension() == -1);
  assert(T0c.number_of_vertices() == 0);
  assert(T0c.number_of_finite_cells() == 0);
  assert(std::distance(T0c.conforming_edges_begin(), T0c.conforming_edges_end()) == 0);
  assert(T0c.is_valid());

  // ---
  Cls T0d(std::cbegin(empty_seg_constraints), std::cend(empty_seg_constraints), gt);
  assert(T0d.dimension() == -1);
  assert(T0d.number_of_vertices() == 0);
  assert(T0d.number_of_finite_cells() == 0);
  assert(std::distance(T0d.conforming_edges_begin(), T0d.conforming_edges_end()) == 0);
  assert(T0d.is_valid());

  // 1D
  // ---
  Cls T1(std::cbegin(pp_constraints), std::cend(pp_constraints));
  assert(T1.dimension() == 1);
  assert(T1.number_of_vertices() == 2);
  assert(T1.number_of_finite_cells() == 0);
  assert(std::distance(T1.conforming_edges_begin(), T1.conforming_edges_end()) == pp_constraints.size());
  assert(T1.is_valid());

  // 2D
  // ---
  Cls T2(std::cbegin(seg_constraints_2), std::cend(seg_constraints_2));
  assert(T2.dimension() == 2);
  assert(T2.number_of_vertices() == 3);
  assert(T2.number_of_finite_cells() == 0);
  assert(std::distance(T2.conforming_edges_begin(), T2.conforming_edges_end()) == seg_constraints_2.size());
  assert(T2.is_valid());

  // 3D
  // ---
  Cls T3(std::cbegin(seg_constraints_3), std::cend(seg_constraints_3));
  assert(T3.dimension() == 3);
  assert(T3.number_of_vertices() == 4);
  assert(T3.number_of_finite_cells() == 1);
  assert(std::distance(T3.conforming_edges_begin(), T3.conforming_edges_end()) == seg_constraints_3.size());
  assert(T3.is_valid());

  // --- Copy
  Cls T0e(T0d);
  assert(T0e.dimension() == -1);
  assert(T0e.number_of_vertices() == 0);
  assert(T0e.number_of_finite_cells() == 0);
  assert(T0e.is_valid());
  assert(T0e == T0d);

  Cls T1b(T1);
  assert(T1b.dimension() == 1);
  assert(T1b.number_of_vertices() == 2);
  assert(T1b.number_of_finite_cells() == 0);
  assert(std::distance(T1b.conforming_edges_begin(), T1b.conforming_edges_end()) == pp_constraints.size());
  assert(T1b.is_valid());
  assert(T1b == T1);

  Cls T2b(T2);
  assert(T2b.dimension() == 2);
  assert(T2b.number_of_vertices() == 3);
  assert(T2b.number_of_finite_cells() == 0);
  assert(std::distance(T2b.conforming_edges_begin(), T2b.conforming_edges_end()) == seg_constraints_2.size());
  assert(T2b.is_valid());
  assert(T2b == T2);

  Cls T3b(T3);
  assert(T3b.dimension() == 3);
  assert(T3b.number_of_vertices() == 4);
  assert(T3b.number_of_finite_cells() == 1);
  assert(std::distance(T3b.conforming_edges_begin(), T3b.conforming_edges_end()) == seg_constraints_3.size());
  assert(T3b.is_valid());
  assert(T3b == T3);

  // --- Assignment
  Cls T0f = T0d;
  assert(T0f.dimension() == -1);
  assert(T0f.number_of_vertices() == 0);
  assert(T0f.number_of_finite_cells() == 0);
  assert(T0f.is_valid());
  assert(T0f == T0d);

  Cls T1c = T1;
  assert(T1c.dimension() == 1);
  assert(T1c.number_of_vertices() == 2);
  assert(T1c.number_of_finite_cells() == 0);
  assert(std::distance(T1c.conforming_edges_begin(), T1c.conforming_edges_end()) == pp_constraints.size());
  assert(T1c.is_valid());
  assert(T1c == T1);

  Cls T2c = T2;
  assert(T2c.dimension() == 2);
  assert(T2c.number_of_vertices() == 3);
  assert(T2c.number_of_finite_cells() == 0);
  assert(std::distance(T2c.conforming_edges_begin(), T2c.conforming_edges_end()) == seg_constraints_2.size());
  assert(T2c.is_valid());
  assert(T2c == T2);

  Cls T3c = T3;
  assert(T3c.dimension() == 3);
  assert(T3c.number_of_vertices() == 4);
  assert(T3c.number_of_finite_cells() == 1);
  assert(std::distance(T3c.conforming_edges_begin(), T3c.conforming_edges_end()) == seg_constraints_3.size());
  assert(T3c.is_valid());
  assert(T3c == T3);

  // --- Move
  namespace test_tr_23 = CGAL::Testsuite::Triangulation_23;
  test_tr_23::test_move_semantic(T0);
  test_tr_23::test_move_semantic(T1);
  test_tr_23::test_move_semantic(T2);
  test_tr_23::test_move_semantic(T3);

  // --- Swap
  Cls T0g;
  T0.swap(T0g);
  assert(T0g.dimension() == -1);
  assert(T0g.number_of_vertices() == 0);
  assert(T0g.number_of_finite_cells() == 0);
  assert(T0g.is_valid());
  assert(T0g == T0);
  T0g.swap(T0);

  Cls T3d;
  T3.swap(T3d);
  assert(T3d.dimension() == 2);
  assert(T3d.number_of_vertices() == 3);
  assert(T3d.number_of_finite_cells() == 1);
  assert(std::distance(T3d.conforming_edges_begin(), T3d.conforming_edges_end()) == seg_constraints_3.size());
  assert(T3d.is_valid());
  T3.swap(T3d);

  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  /// Insertion of constraints

  // --- 0D
  T0.insert_constraint(v1d[0], v1d[0]); // should give an assertion

  // --- 1D
  T0.insert_constraint(v1d[0], v1d[2]);
  assert(T1.dimension() == 1);
  assert(T1.number_of_vertices() == 2);
  assert(T1.number_of_finite_cells() == 0);
  assert(T1.is_valid());
  assert(std::distance(T0.conforming_edges_begin(), T0.conforming_edges_end()) == 1);
  assert(T0.is_constrained(*(T0.finite_edges_begin())));

  T0.insert_constraint(Segment{v1d[0], v1d[2]});
  assert(T1.dimension() == 1);
  assert(T1.number_of_vertices() == 2);
  assert(T1.number_of_finite_cells() == 0);
  assert(T1.is_valid());
  assert(std::distance(T0.conforming_edges_begin(), T0.conforming_edges_end()) == 1);
  assert(T0.is_constrained(*(T0.finite_edges_begin())));

  // --- exactly the same
  T0.insert_constraint(v1d[0], v1d[2]);
  assert(T1.dimension() == 1);
  assert(T1.number_of_vertices() == 2);
  assert(T1.number_of_finite_cells() == 0);
  assert(T1.is_valid());
  assert(std::distance(T0.conforming_edges_begin(), T0.conforming_edges_end()) == 1);
  assert(T0.is_constrained(*(T0.finite_edges_begin())));

  // --- same constraint, in the other direction
  T0.insert_constraint(v1d[2], v1d[0]);
  assert(T1.dimension() == 1);
  assert(T1.number_of_vertices() == 2);
  assert(T1.number_of_finite_cells() == 0);
  assert(std::distance(T0.conforming_edges_begin(), T0.conforming_edges_end()) == 1);
  assert(T0.incident_constraints())
  assert(T1.is_valid());

  // --- constraint sharing a vertex with a previous constraint
  T0.insert_constraint(v1d[2], v1d[4]);
  assert(T1.dimension() == 1);
  assert(T1.number_of_vertices() == 4);
  assert(T1.number_of_finite_cells() == 0);
  assert(std::distance(T0.conforming_edges_begin(), T0.conforming_edges_end()) == 3);
  assert(T1.is_valid());

  // --- disconnected constraint
  T0.insert_constraint(v1d[6], v1d[7]);
  assert(T1.dimension() == 1);
  assert(T1.number_of_vertices() == 5);
  assert(T1.number_of_finite_cells() == 0);
  assert(std::distance(T0.conforming_edges_begin(), T0.conforming_edges_end()) == 4);
  assert(T1.is_valid());

  // --- overlapping constraint
  T0.insert_constraint(v1d[6], v1d[8]);
  assert(T1.dimension() == 1);
  assert(T1.number_of_vertices() == 6);
  assert(T1.number_of_finite_cells() == 0);
  assert(std::distance(T0.conforming_edges_begin(), T0.conforming_edges_end()) == 5);
  assert(T1.is_valid());

  // --- another overlapping constraint
  T0.insert_constraint(v1d[1], v1d[3]);
  assert(T1.dimension() == 1);
  assert(T1.number_of_vertices() == 8);
  assert(T1.number_of_finite_cells() == 0);
  assert(std::distance(T0.conforming_edges_begin(), T0.conforming_edges_end()) == 7);
  assert(T1.is_valid());

  // --- existing constraint (overlapping)
  T0.insert_constraint(v1d[1], v1d[4]);
  assert(T1.dimension() == 1);
  assert(T1.number_of_vertices() == 8);
  assert(T1.number_of_finite_cells() == 0);
  assert(std::distance(T0.conforming_edges_begin(), T0.conforming_edges_end()) == 7);
  assert(T1.is_valid());

  // --- Random
  for(int i=0; i<10; ++i)
  {
    Cls T1d;
    T1d.insert(std::cbegin(v1d), std::cend(v1d));
    for(int i=0; i<10; ++i)
    {
      const int i0 = rnd.get_int(0,10);
      const int i1 = rnd.get_int(0,10);
      T1d.insert_constraint(v1d[i0], v1d[i1]);
      assert(T1d.is_valid())
    }
  }

  // determinism
  Cls Ta(std::cbegin(v1d), std::cend(v1d)), Tb(std::cbegin(v1d), std::cend(v1d));
  assert(Ta == Tb);

  for (Finite_vertices_iterator ita = Ta.finite_vertices_begin(),
                                itb = Tb.finite_vertices_begin(),
                                end = Ta.finite_vertices_end(); ita != end; ++ita, ++itb)
  {
    assert(ita->point() == itb->point());
  }

  for (Finite_cells_iterator ita = Ta.finite_cells_begin(),
                             itb = Tb.finite_cells_begin(),
                             end = Ta.finite_cells_end(); ita != end; ++ita, ++itb)
  {
    assert(ita->vertex(0)->point() == itb->vertex(0)->point());
    assert(ita->vertex(1)->point() == itb->vertex(1)->point());
    assert(ita->vertex(2)->point() == itb->vertex(2)->point());
    assert(ita->vertex(3)->point() == itb->vertex(3)->point());
  }

  // ---
  // 2D

  Cls T2d;
  T2d.insert_constraint({Point{0,0,0}, Point{1,0,0}});
  T2d.insert_constraint({Point{1,0,0}, Point{1,1,0}});
  assert(T2d.dimension() == 2);
  assert(T2d.number_of_vertices() == 3);
  assert(T2d.number_of_finite_cells() == 0);
  assert(std::distance(T2d.conforming_edges_begin(), T2d.conforming_edges_end()) == 2);
  assert(T2d.is_valid());

  // --- Same constraint
  T2d.insert_constraint({Point{1,0,0}, Point{1,1,0}});
  assert(T2d.dimension() == 2);
  assert(T2d.number_of_vertices() == 3);
  assert(T2d.number_of_finite_cells() == 0);
  assert(std::distance(T2d.conforming_edges_begin(), T2d.conforming_edges_end()) == 2);
  assert(T2d.is_valid());

  // --- Same constraint, opposite direction
  T2d.insert_constraint({Point{1,1,0}, Point{1,0,0}});
  assert(T2d.dimension() == 2);
  assert(T2d.number_of_vertices() == 3);
  assert(T2d.number_of_finite_cells() == 0);
  assert(std::distance(T2d.conforming_edges_begin(), T2d.conforming_edges_end()) == 2);
  assert(T2d.is_valid());

  // --- Disconnected constraint
  T2d.insert_constraint({Point{2,2,0}, Point{3,3,0}});
  assert(T2d.dimension() == 2);
  assert(T2d.number_of_vertices() == 5);
  assert(T2d.number_of_finite_cells() == 0);
  assert(std::distance(T2d.conforming_edges_begin(), T2d.conforming_edges_end()) == 3);
  assert(T2d.is_valid());

  // --- Sharing a vertex
  T2d.insert_constraint({Point{1,1,0}, Point{1,2,0}});
  assert(T2d.dimension() == 2);
  assert(T2d.number_of_vertices() == 6);
  assert(T2d.number_of_finite_cells() == 0);
  assert(std::distance(T2d.conforming_edges_begin(), T2d.conforming_edges_end()) == 4);
  assert(T2d.is_valid());

  // --- Sharing two vertices
  T2d.insert_constraint({Point{1,1,0}, Point{2,2,0}});
  assert(T2d.dimension() == 2);
  assert(T2d.number_of_vertices() == 6);
  assert(T2d.number_of_finite_cells() == 0);
  assert(std::distance(T2d.conforming_edges_begin(), T2d.conforming_edges_end()) == 5);
  assert(T2d.is_valid());

  // --- Overlapping
  T2d.insert_constraint({Point{-2,0,0}, Point{2,0,0}});
  assert(T2d.dimension() == 2);
  assert(T2d.number_of_vertices() == 7);
  assert(T2d.number_of_finite_cells() == 0);
  assert(std::distance(T2d.conforming_edges_begin(), T2d.conforming_edges_end()) == 7);
  assert(T2d.is_valid());

  // --- Overlapping
  T2d.insert_constraint({Point{-1,0,0}, Point{1,0,0}});
  assert(T2d.dimension() == 2);
  assert(T2d.number_of_vertices() == 8);
  assert(T2d.number_of_finite_cells() == 0);
  assert(std::distance(T2d.conforming_edges_begin(), T2d.conforming_edges_end()) == 8);
  assert(T2d.is_valid());

  // --- Crossing in interior
  T2d.insert_constraint({Point{2,3,0}, Point{3,2,0}});
  assert(T2d.dimension() == 2);
  assert(T2d.number_of_vertices() == 10);
  assert(T2d.number_of_finite_cells() == 0);
  assert(std::distance(T2d.conforming_edges_begin(), T2d.conforming_edges_end()) == 11);
  assert(T2d.is_valid());

  // --- Crossing at vertex
  T2d.insert_constraint({Point{3,4,0}, Point{3,2,0}});
  assert(T2d.dimension() == 2);
  assert(T2d.number_of_vertices() == 12);
  assert(T2d.number_of_finite_cells() == 0);
  assert(std::distance(T2d.conforming_edges_begin(), T2d.conforming_edges_end()) == 12);
  assert(T2d.is_valid());

  // --- Random
  for(int i=0; i<10; ++i)
  {
    Cls T2e;
    T2e.insert(std::cbegin(v2d), std::cend(v2d));
    for(int i=0; i<100; ++i)
    {
      const int i0 = rnd.get_int(0,10);
      const int i1 = rnd.get_int(0,10);
      T2e.insert_constraint(v2d[i0], v2d[i1]);
      assert(T2e.dimension() == 2);
      assert(T2e.is_valid());
    }
  }

  // ---
  // 3D

  Cls T3e;
  T3e.insert_constraint({Point{0,0,0}, Point{1,0,0}});
  T3e.insert_constraint({Point{0,0,0}, Point{0,1,0}});
  T3e.insert_constraint({Point{0,0,0}, Point{0,0,1}});
  assert(T3e.dimension() == 3);
  assert(T3e.number_of_vertices() == 4);
  assert(T3e.number_of_finite_cells() == 1);
  assert(std::distance(T3e.conforming_edges_begin(), T3e.conforming_edges_end()) == 3);
  assert(T3e.is_valid());

  // --- Same constraint
  T3e.insert_constraint({Point{0,0,0}, Point{1,0,0}});
  assert(T3e.dimension() == 3);
  assert(T3e.number_of_vertices() == 4);
  assert(T3e.number_of_finite_cells() == 1);
  assert(std::distance(T3e.conforming_edges_begin(), T3e.conforming_edges_end()) == 3);
  assert(T3e.is_valid());

  // --- Same constraint, opposite direction
  T3e.insert_constraint({Point{1,0,0}, Point{0,0,0}});
  assert(T3e.dimension() == 3);
  assert(T3e.number_of_vertices() == 4);
  assert(T3e.number_of_finite_cells() == 1);
  assert(std::distance(T3e.conforming_edges_begin(), T3e.conforming_edges_end()) == 3);
  assert(T3e.is_valid());

  // --- Disconnected constraint
  T3e.insert_constraint({Point{1,1,1}, Point{2,2,2}});
  assert(T3e.dimension() == 3);
  assert(T3e.number_of_vertices() == 6);
  assert(std::distance(T3e.conforming_edges_begin(), T3e.conforming_edges_end()) == 4);
  assert(T3e.is_valid());

  // --- Sharing a vertex
  T3e.insert_constraint({Point{1,1,1}, Point{2,2,1}});
  assert(T3e.dimension() == 3);
  assert(T3e.number_of_vertices() == 7);
  assert(std::distance(T3e.conforming_edges_begin(), T3e.conforming_edges_end()) == 5);
  assert(T3e.is_valid());

  // --- Sharing two vertices
  T3e.insert_constraint({Point{2,2,1}, Point{2,2,2}});
  assert(T3e.dimension() == 3);
  assert(T3e.number_of_vertices() == 7);
  assert(std::distance(T3e.conforming_edges_begin(), T3e.conforming_edges_end()) == 6);
  assert(T3e.is_valid());

  // --- Overlapping
  T3e.insert_constraint({Point{2,2,0}, Point{2,2,3}});
  assert(T3e.dimension() == 3);
  assert(T3e.number_of_vertices() == 9);
  assert(std::distance(T3e.conforming_edges_begin(), T3e.conforming_edges_end()) == 8);
  assert(T3e.is_valid());

  // --- Crossing
  T3e.insert_constraint({Point{2,1,4}, Point{2,3,3}});
  assert(T3e.dimension() == 3);
  assert(T3e.number_of_vertices() == 11);
  assert(std::distance(T3e.conforming_edges_begin(), T3e.conforming_edges_end()) == 11);
  assert(T3e.is_valid());

  // --- Random
  for(int i=0; i<10; ++i)
  {
    Cls T3f;
    T3f.insert(std::cbegin(v3d_1), std::cend(v3d_1));
    for(int i=0; i<100; ++i)
    {
      const int i0 = rnd.get_int(0,10);
      const int i1 = rnd.get_int(0,10);
      T2e.insert_constraint(v3d_1[i0], v3d_1[i1]);
      assert(T2f.dimension() == 2);
      assert(T2f.is_valid());
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  /// Throwing (constraint tags)

  // test throwing/not throwing on intersecting constraints
  _test_cdt_throwing<Triang>(Point(0, 0, 0), Point(1, 1, 0), Point(2, 2, 0), Point(2, 3, 0), NO_INTERSECTION);
  _test_cdt_throwing<Triang>(Point(0, 0, 0), Point(1, 1, 0), Point(1, 1, 0), Point(2, 2, 0), NO_INTERSECTION); // common point
  _test_cdt_throwing<Triang>(Point(2, 2, 0), Point(0, 0, 0), Point(2, 2, 0), Point(3, 3, 0), NO_INTERSECTION); // ^
  _test_cdt_throwing<Triang>(Point(0, 0, 0), Point(2, 2, 0), Point(1, 1, 0), Point(3, 3, 0), INTERSECTION_WITHOUT_CONSTRUCTION); // overlapping
  _test_cdt_throwing<Triang>(Point(2, 2, 0), Point(0, 0, 0), Point(1, 1, 0), Point(3, 3, 0), INTERSECTION_WITHOUT_CONSTRUCTION); // ^
  _test_cdt_throwing<Triang>(Point(2, 2, 0), Point(0, 0, 0), Point(0, 0, 0), Point(3, 3, 0), INTERSECTION_WITHOUT_CONSTRUCTION); // ^
  _test_cdt_throwing<Triang>(Point(2, 2, 0), Point(0, 0, 0), Point(3, 3, 0), Point(0, 0, 0), INTERSECTION_WITHOUT_CONSTRUCTION); // ^
  _test_cdt_throwing<Triang>(Point(0, 0, 0), Point(3, 3, 0), Point(1, 1, 0), Point(2, 2, 0), INTERSECTION_WITHOUT_CONSTRUCTION); // contains
  _test_cdt_throwing<Triang>(Point(3, 3, 0), Point(0, 0, 0), Point(1, 1, 0), Point(2, 2, 0), INTERSECTION_WITHOUT_CONSTRUCTION); // ^
  _test_cdt_throwing<Triang>(Point(1, 1, 0), Point(2, 2, 0), Point(3, 3, 0), Point(0, 0, 0), INTERSECTION_WITHOUT_CONSTRUCTION); // ^
  _test_cdt_throwing<Triang>(Point(2, 2, 0), Point(1, 1, 0), Point(3, 3, 0), Point(0, 0, 0), INTERSECTION_WITHOUT_CONSTRUCTION); // ^
  _test_cdt_throwing<Triang>(Point(3, 3, 0), Point(0, 0, 0), Point(0, 0, 0), Point(3, 3, 0), INTERSECTION_WITHOUT_CONSTRUCTION); // same constraint
  _test_cdt_throwing<Triang>(Point(3, 3, 0), Point(0, 0, 0), Point(1, 1, 0), Point(1, 1, 0), INTERSECTION_WITHOUT_CONSTRUCTION); // degenerate entry
  _test_cdt_throwing<Triang>(Point(0, 0, 0), Point(0, 0, 0), Point(0, 0, 0), Point(0, 0, 0), NO_INTERSECTION); // degenerate same entry

  // extremity on the interior of another segment
  _test_cdt_throwing<Triang>(Point(0, 0, 0), Point(2, 0, 0), Point(1, 0, 0), Point(1, 4, 0), INTERSECTION_WITHOUT_CONSTRUCTION);

  // non-aligned
  _test_cdt_throwing<Triang>(Point(0, 0, 0), Point(3, 3, 0), Point(1, 0, 0), Point(4, 0, 0), NO_INTERSECTION);
  _test_cdt_throwing<Triang>(Point(0, 2, 0), Point(2, 2, 0), Point(1, 0, 0), Point(1, 3, 0), INTERSECTION); // generic intersection
  _test_cdt_throwing<Triang>(Point(0, 0, 0), Point(2, 2, 0), Point(1, 3, 0), Point(1, 0, 0), INTERSECTION);

  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  /// Duality

  // @todo

  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  /// Traversal

  typename Cls::Conforming_edges_iterator ceit = T2.conforming_edges_begin(),
                                          cend = T2.conforming_edges_end();


  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  /// Input/Output



}
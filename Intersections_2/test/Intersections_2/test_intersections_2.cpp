// 2D intersection tests.

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/internal/Exact_type_selector.h>

#include <CGAL/Intersections_2/Triangle_2_Triangle_2.h>
#include <CGAL/Intersections_2/Iso_rectangle_2_Line_2.h>
#include <CGAL/Intersections_2/Iso_rectangle_2_Ray_2.h>

#include <CGAL/Intersection_traits_2.h>
#include <CGAL/Line_2_Iso_rectangle_2_intersection.h>
#include <CGAL/Intersections_2/Bbox_2_Circle_2.h>
#include <CGAL/Intersections_2/Bbox_2_Point_2.h>
#include <CGAL/Intersections_2/Circle_2_Iso_rectangle_2.h>
#include <CGAL/Intersections_2/Circle_2_Point_2.h>
#include <CGAL/Intersections_2/Point_2_Point_2.h>

#include <vector>
#include <iostream>
#include <cassert>

const double epsilon = 0.001;

struct randomint {
  randomint() ;
  int	get() const { return sequence[cur]; }
  int next() { cur = (cur+1)%11; return get();}
private:
  int sequence[11];
  int cur;
};

inline randomint::randomint()
{
  cur = 0;
  sequence[0] = 19;
  sequence[1] = 5;
  sequence[2] = 17;
  sequence[3] = 13;
  sequence[4] = 29;
  sequence[5] = 2;
  sequence[6] = 23;
  sequence[7] = 31;
  sequence[8] = 3;
  sequence[9] = 37;
  sequence[10] = 11;
}

randomint ri;

inline double to_nt(int d)
{
    return double(d);
}

template < typename K >
struct Test {

  typedef typename K::Point_2               P;
  typedef typename K::Line_2                L;
  typedef typename K::Segment_2             S;
  typedef typename K::Ray_2                 R;
  typedef typename K::Triangle_2            T;
  typedef typename K::Iso_rectangle_2       Rec;
  typedef typename K::Circle_2              C;
  typedef std::vector<P>              Pol;


  template < typename Type >
  bool approx_equal_nt(const Type &t1, const Type &t2)
  {
	if (t1 == t2)
		return true;
	if (CGAL::abs(t1 - t2) / (CGAL::max)(CGAL::abs(t1), CGAL::abs(t2)) < epsilon)
		return true;
	std::cout << " Approximate comparison failed between : " << t1 << "  and  " << t2 << "\n";
	return false;
  }

  template < typename Type >
  bool approx_equal(const Type&t1, const Type&t2)
  {
	return t1 == t2;
	// we need approx equal to check approx kernels, but maybe we should only test with exact kernels
	// (approx kernels were useful before, when the text output was checked by diff ?)
	// idea : test containment with intervals ?  or use some "epsilon double"?
	// I need to convert the text output to exact rationals in the source...
	// Well, for now the current scheme works.
  }

  bool approx_equal(const P & p, const P & q)
  {
	return approx_equal_nt(p.x(), q.x()) && approx_equal_nt(p.y(), q.y());
  }

  bool approx_equal(const Pol & p, const Pol & q)
  {
	if (p.size() != q.size())
		return false;
	for(typename Pol::const_iterator itp = p.begin(), itq = q.begin(); itp != p.end(); ++itp, ++itq)
		if (!approx_equal(*itp, *itq))
			return false;
	return true;
  }

  template < typename O1, typename O2>
  void check_no_intersection(const O1& o1, const O2& o2)
  {
	assert(!CGAL::do_intersect(o1, o2));
	assert(!CGAL::do_intersect(o2, o1));
        assert(!CGAL::intersection(o2, o1));
    
	//check with the functors
	typename CGAL::Kernel_traits<O1>::Kernel::Do_intersect_2 do_2;
	typename CGAL::Kernel_traits<O1>::Kernel::Intersect_2 i_2;
	assert(!do_2(o1, o2));
        assert(!i_2(o1, o2));
	assert(!do_2(o2, o1));
        assert(!i_2(o2, o1));
  }

  template < typename Res, typename O1, typename O2 >
  void check_intersection(const O1& o1, const O2& o2)
  {
	Res tmp;
	assert(CGAL::do_intersect(o1, o2));
	assert(CGAL::assign(tmp, CGAL::intersection(o1, o2)));
	assert(CGAL::do_intersect(o2, o1));
	assert(CGAL::assign(tmp, CGAL::intersection(o2, o1)));
  }

  template < typename Res, typename O1, typename O2 >
  void check_intersection(const O1& o1, const O2& o2, const Res& result, bool do_opposite=true)
  {
	Res tmp;
	assert(CGAL::do_intersect(o1, o2));
	assert(CGAL::assign(tmp, CGAL::intersection(o1, o2)));
	assert(approx_equal(tmp, result));
	if (do_opposite) {
	  assert(CGAL::do_intersect(o2, o1));
	  assert(CGAL::assign(tmp, CGAL::intersection(o2, o1)));
	  assert(approx_equal(tmp, result));
	}
  }

  template < typename O >
  void check_intersection(const O& o)
  {
    return check_intersection(o, o, o);
  }

  P p(int x, int y)
  {
    int w = ri.next();
    return P(to_nt(x*w), to_nt(y*w), to_nt(w));
  }

  void B_P()
  {
    CGAL::Bbox_2 bb(0,0,10,10);
    P p(1,0), bl(0,0), tr(10,10);
    C c(bl,1);
    Rec r(bl,tr);
    check_intersection(bb,p,p,true);
    check_intersection(c,p,p,true);
    assert(do_intersect(r,c));
  }
  
  void L_L()
  {
    std::cout << "Line - Line\n";

    // no intersection
    check_no_intersection  (L(p(0, 0), p(10,10)), L(p(8,7), p(1, 0)));

    // point intersection
    check_intersection     (L(p(0, 0), p(10, 0)), L(p( 1, 7), p(  1,  -2)), P(1,0));
    check_intersection     (L(p(0,-1), p(10, 0)), L(p( 2, 1), p(  8,  -6)), P(3.42105,-0.657895));
    check_intersection<P>  (L(p(0, 0), p( 0, 4)), L(p(-1,-3), p(-10, -10)));

    // line intersection
    check_intersection     (L(p( 1, 1), p( 5, 5)));
    check_intersection<L>  (L(p( 1, 1), p( 5, 5)), L(p(3, 3), p( 7,  7))); // ps0 < ps1 < pt0 < pt1
    check_intersection<L>  (L(p( 5, 5), p( 1, 1)), L(p(3, 3), p( 7,  7))); // L0 & L1 have opposite directions
    check_intersection<L>  (L(p( 0, 0), p(10, 0)), L(p(0, 0), p(10,  0))); // ps0 == ps1 < pt0 == pt1
    check_intersection<L>  (L(p(10, 0), p( 0, 0)), L(p(0, 0), p(10,  0))); // ps1 == pt0 < ps0 == pt1
    check_intersection<L>  (L(p( 0, 0), p(10, 0)), L(p(1, 0), p( 8,  0))); // ps0 < ps1 < pt1 < pt0
    check_intersection<L>  (L(p( 0, 0), p(10, 0)), L(p(8, 0), p( 1,  0))); // L0 & L1 have opposite directions
  }

  void S_S()
  {
    std::cout << "Segment - Segment\n";

    // no intersection
    check_no_intersection  (S(p(29,  16), p( 28,   9)), S(p( 30,  12), p( 29,   6)));
    check_no_intersection  (S(p( 0,   0), p(103567,9826)), S(p(10000,3782), p(76250,83736)));

    // point intersection
    check_intersection     (S(p( 0,   0), p( 10,   0)), S(p(  1,   6), p(  1,  -3)), P(1, 0));
    check_intersection     (S(p( 0,  -1), p( 10,   0)), S(p(  2,   1), p(  8,  -6)), P(3.42105, -0.657895));
    check_intersection     (S(p( 0,   0), p( 10,   0)), S(p(  0,   0), p(  4,  -7)), P(0, 0)); // meeting at an extremity

    // segment intersection
    check_intersection     (S(p( 0,   0), p( 10,   0)));
    check_intersection     (S(p( 0,   0), p( 10,   0)), S(p(  1,   0), p(  8,   0)), S(P(  1,   0), P(  8,   0)));
    check_intersection     (S(p(68, 106), p(192, 106)), S(p(150, 106), p(255, 106)), S(P(150, 106), P(192, 106)));
    check_intersection     (S(p( 1,  10), p(  1,   2)), S(p(  1,   7), p(  1,   3)), S(P(  1,   3), P(  1,   7)));
  }

  void R_R()
  {
    std::cout << "Ray - Ray\n";

    // no intersection
    check_no_intersection  (R(p( 3,   4), p(  5,   7)), R(p(  2,   0), p(  2,   4)));

    // point intersection
    check_intersection     (R(p( 2,  -1), p(  2,   1)), R(p(  1,   3), p(  2,   3)), P(2, 3));
    check_intersection     (R(p( 0,  -1), p( 10,   0)), R(p(  2,   1), p(  8,  -6)), P(3.42105, -0.657895));
    check_intersection     (R(p( 0,   0), p( 10,   0)), R(p(  1,   6), p(  1,  -3)), P(1, 0));
    check_intersection     (R(p( 0,   0), p( 10,   0)), R(p(  15,  0), p( 18,  30)), P(15, 0)); // R1's source on R0's path
    check_intersection     (R(p( 0,   0), p( 10,   0)), R(p(  0,   0), p( -3,   4)), P(0, 0)); // same source but different directions

    // segment intersection (same supporting line, different directions)
    check_intersection<S>  (R(p( 2,   4), p(  6,   1)), R(p(  14, -5), p(10, -2)));
    check_intersection<S>  (R(p( 2,   4), p( 10,  -2)), R(p(   6,  1), p(-2,  7)));

    // ray intersection
    check_intersection     (R(p( 0,   0), p( 10,   0)));
    check_intersection<R>  (R(p( 0,   0), p( 10,   0)), R(p(  -1,  0), p(0,   0))); // R0 'runs' into R1's source
  }

  void S_R()
  {
    std::cout << "Segment - Ray\n";

    // no intersection
    check_no_intersection  (S(p( 2,  -1), p(  2,   1)), R(p(  1,   3), p(  2,   3)));

    // point intersection
    check_intersection     (S(p( 0,  -1), p( 10,   0)), R(p(  2,   1), p(  8,  -6)), P(3.42105, -0.657895));
    check_intersection     (S(p( 0,   0), p( 10,   0)), R(p(  1,   6), p(  1,  -3)), P(1, 0));
    check_intersection     (S(p( 0,   0), p( 10,   0)), R(p(  1,   6), p(  1,  -3)), P(1, 0));
    check_intersection     (S(p( 0,   0), p( 10,   0)), R(p(  0,   0), p(-10,   4)), P(0, 0)); // start of ray is exactly on the segment
    check_intersection     (S(p( 0,   0), p( 10,   0)), R(p(  4,   0), p(-10,   4)), P(4, 0)); // start of ray is a segment extremity

    // segment intersection
    check_intersection     (S(p(  0,   0), p(  1,  0)), R(p(  0,   0), p( 10,   0)), S(P(0, 0), P(1,0)));
    check_intersection<S>  (S(p(  3,   7), p(  2,  5)), R(p(  1,   3), p(  4,   9)));
  }

  void L_R()
  {
    std::cout << "Line - Ray\n";

    // no intersection
    check_no_intersection  (L(p( 2,  -1), p(  2,   1)), R(p(  1,   -3), p(  1,   3)));

    // point intersection
    check_intersection     (L(p( 2,  -1), p(  2,   1)), R(p(  1,   3), p(  2,   3)), P(2, 3));
    check_intersection     (L(p( 0,  -1), p( 10,   0)), R(p(  2,   1), p(  8,  -6)), P(3.42105, -0.657895));
    check_intersection     (L(p( 0,   0), p( 10,   0)), R(p(  1,   6), p(  1,  -3)), P(1, 0));
    check_intersection     (L(p( 0,   0), p(  4,   5)), R(p(  1,   8), p(  7,   2)), P(4, 5)); // ray starts on the line

    // ray intersection
    check_intersection<R>  (L(p( 2,  -1), p(  2,   1)), R(p(  2,   -3), p(  2,   3)));
    check_intersection<R>  (L(p( 2,  -1), p(  2,   1)), R(p(  2,    3), p(  2,  -3))); // opposite direction
  }

  void S_L()
  {
    std::cout << "Segment - Line\n";

    // no intersection
    check_no_intersection  (S(p( 2,  -1), p(  2,   1)), L(p(  1,   3), p(  2,   3)));

    // point intersection
    check_intersection     (S(p( 0,  -1), p( 10,   0)), L(p(  2,   1), p(  8,  -6)), P(3.42105, -0.657895));
    check_intersection     (S(p( 0,   0), p( 10,   0)), L(p(  1,   6), p(  1,  -3)), P(1, 0));
    check_intersection     (S(p( 0,   0), p( 10,   0)), L(p(  1,   6), p(  2,  12)), P(0, 0));

    // segment intersection
    check_intersection<S>  (S(p(-3,   5), p( 12,   1)), L(p( 12,   1), p( 27,  -3)));
    check_intersection<S>  (S(p(-3,   5), p( 12,   1)), L(p(-18,   9), p( 27,  -3)));
    check_intersection<S>  (S(p(-3,   5), p( 12,   1)), L(p( 27,  -3), p(-18,   9)));
  }

  void T_T()
  {
    std::cout << "Triangle - Triangle\n";

    // no intersection
    check_no_intersection  (T(p( -10,-10), p(  0,  10), p( 20, -5)), T(p(   90, -10), p(100,  10), p(120, -5)));

    // point intersection
    check_intersection     (T(p( -10,  0), p( 10,   0), p(0,  3)), T(p( -12,   3), p( 12,   3), p(1, 5)), P(0, 3)); // intersection on an edge
    check_intersection     (T(p( -25, -4), p( 13, -18), p(0,  3)), T(p( -12,   5), p( 12,   5), p(0, 3)), P(0, 3)); // intersection at a vertex

    // segment intersection
    check_intersection<S>  (T(p( -10,  0), p( 10,   0), p(0, -3)), T(p(  -8,   0), p( 12,   0), p(1, 5)));
    check_intersection<S>  (T(p( -10,  0), p( 10,   0), p(0, -3)), T(p(  -8,   0), p(  8,   0), p(1, 5)));
    check_intersection<S>  (T(p( -10,  0), p( 10,   0), p(0, -3)), T(p( -10,   0), p( 10,   0), p(1, 5)));
    check_intersection<S>  (T(p( -10,  0), p( 10,   0), p(0, -3)), T(p(  10,   0), p(-10,   0), p(1, 5)));
    check_intersection<S>  (T(p( -10,  0), p( 10,   0), p(0, -3)), T(p( -10,   0), p( 10,   0), p(1, 5)));
    check_intersection<S>  (T(p( -10,  0), p( 10,   0), p(0, -3)), T(p( -12,   0), p( 12,   0), p(1, 5)));

    // triangle intersection
    check_intersection<T>  (T(p(   0, 10), p(-10, -10), p( 20,  -5)), T(p( -10, -10), p(  0,  10), p( 20, -5)));
    check_intersection<T>  (T(p( -12,  1), p(  5,   3), p( -7, -15)), T(p(  29,  -2), p(  0, -13), p(1, 21)));

    // polygon intersection
    Pol pol0;
    pol0.push_back(P(-6, -4));
    pol0.push_back(P( -5.11111, -0.222222 ));
    pol0.push_back(P( 0, 10 ));
    pol0.push_back(P( 8, 4 ));
    check_intersection     (T(p(   0, 10), p(-10, -10), p( 20, -5)), T(p(   2,  30), p( -6,  -4), p(15, 8)), pol0, false);

    Pol pol1;
    pol1.push_back(P( 8,  4));
    pol1.push_back(P( 0, 10 ));
    pol1.push_back(P( -5.11111, -0.222222 ));
    pol1.push_back(P(-6, -4));
    check_intersection     (T(p( -10,-10), p(  0,  10), p( 20, -5)), T(p(   2,  30), p( -6,  -4), p(15, 8)), pol1, false);

    Pol pol2;
    pol2.push_back(P( 10.2222, 2.33333 ));
    pol2.push_back(P( 1.96923, 8.52308 ));
    pol2.push_back(P( -0.680851, 8.6383 ));
    pol2.push_back(P( -5.94872, -1.89744 ));
    pol2.push_back(P( -3.96178, -8.99363 ));
    pol2.push_back(P( 3.5, -7.75 ));
    check_intersection     (T(p( -10,-10), p(  0,  10), p( 20, -5)), T(p(   -9,   9), p( 14,   8), p(-2, -16)), pol2, false);
  }

  void L_T()
  {
    std::cout << "Line - Triangle\n";

    // no intersection
    check_no_intersection  (L(p(-10,   0), p( 10,   0)), T(p(  -12,   3), p(  12,   3), p(1,  5)));

    // point intersection
    check_intersection<P>  (L(p( -8,  30), p(  8,  30)), T(p(    2,  30), p(  14,   2), p(-7, -2)));
    check_intersection<P>  (L(p( -8,  30), p( -7,  30)), T(p(    2,  30), p(  14,   2), p(-7, -2)));

    // segment intersection
    check_intersection<S>  (L(p( -1,  -1), p(  0,  -1)), T(p(  -10,   0), p(  10,   0), p(0, -4)));
    check_intersection<S>  (L(p( -1,  -1), p(  0,  -1)), T(p(  -10,   0), p(  15,   2), p(0, -4)));
    check_intersection<S>  (L(p(-10,   0), p( 10,   0)), T(p(   -8,   0), p(   8,   0), p(1,  5)));
    check_intersection<S>  (L(p(-10,   0), p( 10,   0)), T(p(  -12,   0), p(  12,   0), p(1,  5)));
    check_intersection<S>  (L(p(  0,  10), p(-10, -10)), T(p(    2,  30), p(  -6,  -4), p(15, 8)));
    check_intersection<S>  (L(p(-12,   1), p(  5,   3)), T(p(   29,  -2), p(   0, -13), p( 1, 21)));
    check_intersection<S>  (L(p(-10, -10), p(  0,  10)), T(p(    2,  30), p(  -6,  -4), p(15,  8)));
    check_intersection<S>  (L(p(-10, -10), p(  0,  10)), T(p(   -9,   9), p(  14,   8), p(-2,-16)));
  }

  void R_T()
  {
    std::cout << "Ray - Triangle\n";

    // no intersection
    check_no_intersection  (R(p(-10,   0), p( 10,   0)), T(p(  -12,   3), p(  12,   3), p(1,  5)));

    // point intersection
    check_intersection<P>  (R(p(-2, -16), p(  4,  -20)), T(p(   -9,   9), p(  14,   8), p(-2, -16)));
    check_intersection<P>  (R(p(-8, -1),  p(  -8, -12)), T(p(  -12,   2), p(  10,   3), p(-4,  -4)));
    check_intersection<P>  (R(p(-8, 30),  p(   4,  30)), T(p(    2,  30), p(  14,   2), p(-7,  -2)));
    check_intersection<P>  (R(p(-8, 30),  p(  -7,  30)), T(p(    2,  30), p(  14,   2), p(-7,  -2)));

    // segment intersection
    check_intersection<S>  (R(p( -1,  -1), p(  0,  -1)), T(p(  -10,   0), p(  10,   0), p(0, -4)));
    check_intersection<S>  (R(p( -1,  -1), p(  0,  -1)), T(p(  -10,   0), p(  15,   2), p(0, -4)));
    check_intersection<S>  (R(p(-10,   0), p( 10,   0)), T(p(   -8,   0), p(   8,   0), p(1,  5)));
    check_intersection<S>  (R(p(-10,   0), p( 10,   0)), T(p(  -12,   0), p(  12,   0), p(1,  5)));
    check_intersection<S>  (R(p(  0,  10), p(-10, -10)), T(p(    2,  30), p(  -6,  -4), p(15, 8)));
    check_intersection<S>  (R(p(-12,   1), p(  5,   3)), T(p(   29,  -2), p(   0, -13), p( 1, 21)));
    check_intersection<S>  (R(p(-10, -10), p(  0,  10)), T(p(    2,  30), p(  -6,  -4), p(15,  8)));
    check_intersection<S>  (R(p(-10, -10), p(  0,  10)), T(p(   -9,   9), p(  14,   8), p(-2,-16)));
  }

  void S_T()
  {
    std::cout << "Segment - Triangle\n";

    // no intersection
    check_no_intersection  (S(p(-10, -10), p(  0,  10)), T(p(   90, -10), p( 100,  10), p(120, -5)));
    check_no_intersection  (S(p(-10,   0), p( 10,   0)), T(p(  -12,   3), p(  12,   3), p(1,  5)));

    // point intersection
    check_intersection<P>  (S(p(-2, -16), p(  4,  -20)), T(p(   -9,   9), p(  14,   8), p(-2, -16)));
    check_intersection<P>  (S(p(-8,  -1), p(  -8, -12)), T(p(  -12,   2), p(  10,   3), p(-4,  -4)));
    check_intersection<P>  (S(p(-8, -12), p(  -8,  -1)), T(p(  -12,   2), p(  10,   3), p(-4,  -4)));
    check_intersection<P>  (S(p(-8,  30), p(   1,  30)), T(p(   -2,  30), p(  14,   2), p(-7,  -2)));

    // segment intersection
    check_intersection<S>  (S(p( -1,  -1), p(  0,  -1)), T(p(  -10,   0), p(  10,   0), p(0, -4)));
    check_intersection<S>  (S(p( -1,  -1), p(  0,  -1)), T(p(  -10,   0), p(  15,   2), p(0, -4)));
    check_intersection<S>  (S(p(-10,   0), p( 10,   0)), T(p(   -8,   0), p(   8,   0), p(1,  5)));
    check_intersection<S>  (S(p(-10,   0), p( 10,   0)), T(p(  -12,   0), p(  12,   0), p(1,  5)));
    check_intersection<S>  (S(p(  0,  10), p(-10, -10)), T(p(    2,  30), p(  -6,  -4), p(15, 8)));
    check_intersection<S>  (S(p(-12,   1), p(  5,   3)), T(p(   29,  -2), p(   0, -13), p( 1, 21)));
    check_intersection<S>  (S(p(-10, -10), p(  0,  10)), T(p(    2,  30), p(  -6,  -4), p(15,  8)));
    check_intersection<S>  (S(p(-10, -10), p(  0,  10)), T(p(   -9,   9), p(  14,   8), p(-2,-16)));
  }

  void P_P()
  {
    std::cout << "Point - Point\n";
    check_no_intersection<P>  (p(  8, 4), p(-4,  8));
    check_intersection<P>     (p(  8, 4), p( 8,  4));
  }

  void P_T()
  {
    std::cout << "Point - Triangle\n";

    // no intersection
    check_no_intersection  (p(  8,   6), T(p(    4,   0), p(  12,   4), p(-4,  8)));

    // point intersection
    check_intersection<P>  (p(  8,   4), T(p(    4,   0), p(  12,   4), p(-4,  8)));
    check_intersection<P>  (p(  8,   5), T(p(    4,   0), p(  12,   4), p(-4,  8)));
    check_intersection<P>  (p(  4,   0), T(p(    4,   0), p(  12,   4), p(-4,  8)));
    check_intersection<P>  (p(  12,  4), T(p(    4,   0), p(  12,   4), p(-4,  8)));
    check_intersection<P>  (p(  -4,  8), T(p(    4,   0), p(  12,   4), p(-4,  8)));
  }

  void L_Rec()
  {
    std::cout << "Line - Iso_rectangle\n";

    // no intersection
    check_no_intersection  (L(p( 18,  6), p( 16,  4)), Rec(p( 2,  0), p(6,  3)));

    // point intersection
    check_intersection     (L(p( -1,  0), p( 4,   5)), Rec(p( 0,  0), p(1,  1)), P(0, 1));
    check_intersection<P>  (L(p( -5, 10), p(-1, -12)), Rec(p(-3, -1), p(2,  14)));

    // segment intersection
    check_intersection<S>  (L(p( 18,  6), p( 0,   0)), Rec(p( 2,  0), p(6,  3)));
    check_intersection<S>  (L(p( 18,  6), p( 0,   0)), Rec(p( 2,  0), p(6,  3)));
    check_intersection<S>  (L(p( 2,  14), p( 2, -14)), Rec(p( 2,  0), p(6,  3)));
    check_intersection<S>  (L(p( 6,   1), p( 6,   2)), Rec(p( 2,  0), p(6,  3)));
    check_intersection<S>  (L(p(-1,   3), p(-2,   3)), Rec(p( 2,  0), p(6,  3)));
  }

  void R_Rec()
  {
    std::cout << "Ray - Iso_rectangle\n";

    // no intersection
    check_no_intersection  (R(p( 18,  6), p( 16,  4)), Rec(p( 2,  0), p(6,  3)));

    // point intersection
    check_intersection     (R(p( -1,  0), p( 4,   5)), Rec(p( 0,  0), p(1,  1)), P(0, 1));
    check_intersection     (R(p(  0,  2), p(-4,  14)), Rec(p( 0,  0), p(5,  2)), P(0, 2));
    check_intersection<P>  (R(p( -5, 10), p(-1, -12)), Rec(p(-3, -1), p(2,  14)));

    // segment intersection
    check_intersection<S>  (R(p( 18,  6), p( 0,   0)), Rec(p( 2,  0), p(6,  3)));
    check_intersection<S>  (R(p( 2,  14), p( 2, -14)), Rec(p( 2,  0), p(6,  3)));
    check_intersection<S>  (R(p( 2,   1), p( 2,   4)), Rec(p( 2,  0), p(6,  3)));
    check_intersection<S>  (R(p(-2,   3), p(-1,   3)), Rec(p( 2,  0), p(6,  3)));
  }

  void S_Rec()
  {
    std::cout << "Segment - Iso_rectangle\n";

    // no intersection
    check_no_intersection  (S(p( 18,  6), p( 16,  4)), Rec(p( 2,  0), p(6,  3)));

    // point intersection
    check_intersection     (S(p( -1,  0), p( 4,   5)), Rec(p( 0,  0), p(1,  1)), P(0, 1));
    check_intersection     (S(p(  0,  2), p(-4,  14)), Rec(p( 0,  0), p(5,  2)), P(0, 2));
    check_intersection<P>  (S(p( -5, 10), p(-1, -12)), Rec(p(-3, -1), p(2,  14)));

    // segment intersection
    check_intersection<S>  (S(p( 18,  6), p( 0,   0)), Rec(p( 2,  0), p(6,  3)));
    check_intersection<S>  (S(p( 2,  14), p( 2, -14)), Rec(p( 2,  0), p(6,  3)));
    check_intersection<S>  (S(p( 6,   1), p( 6,   2)), Rec(p( 2,  0), p(6,  3)));
    check_intersection<S>  (S(p(-2,   3), p( 6,   3)), Rec(p( 2,  0), p(6,  3)));
  }

  void Rec_Rec()
  {
    std::cout << "Iso_rectangle - Iso_rectangle\n";

    // no intersection
    check_no_intersection  (Rec(p( -4, -12), p(12, 23)), Rec(p( -4,  24), p(  5, 26)));
    check_no_intersection  (Rec(p( -4, -12), p(12, 23)), Rec(p( 13, -24), p( 14, 15)));

    // point intersection
    check_intersection<Rec> (Rec(p( 10,  12), p(30, 40)), Rec(p(  30,   40), p( 31, 42))/*, p(30, 40)*/);
    check_intersection<Rec> (Rec(p( 10,  12), p(30, 40)), Rec(p(  30,  -13), p( 33, 12))/*, p(30, 12)*/);

    // segment intersection
    check_intersection<Rec>  (Rec(p( 3,  5), p(4, 6)), Rec(p( 2, 4), p( 6, 5)));
    check_intersection<Rec>  (Rec(p( 3,  5), p(4, 6)), Rec(p( 1, 1), p( 3, 8)));
    check_intersection<Rec>  (Rec(p( 3,  5), p(9, 9)), Rec(p( 1, 4), p( 3, 8)));

    // Iso rectangle intersection
    check_intersection     (Rec(p( 10,  12), p(30, 40)));
    check_intersection     (Rec(p( 10,  12), p(30, 40)), Rec(p(  25,   40), p( 26,  103)), Rec(P(25, 40), P(26, 40)));
  }

  void T_Rec()
  {
    std::cout << "Triangle - Iso_rectangle\n";

    // no intersection
    check_no_intersection  (Rec(p( 10,  12), p(30, 40)), T(p(    4,   0), p(  12,   4), p(-4,  8)));

    // point intersection
    check_intersection<P>  (Rec(p( 0,  0), p(1, 1)), T(p(  -1,  0), p(  0,  0), p( 0,  -1))); // intersection at a vertex
    check_intersection<P>  (Rec(p( 0,  0), p(1, 1)), T(p(   0,  0), p( -1,  0), p( 0,  -1))); // inversed triangle
    check_intersection<P>  (Rec(p( 0,  0), p(1, 1)), T(p(  -1,  0), p(  2,  3), p(-4,   6))); // intersection on an edge of the triangle
    check_intersection     (Rec(p( 0,  0), p(3, 3)), T(p( -10,  0), p(  0,  2), p(-1,   4)), p(0, 2)); // intersection on an edge of the iso rectangle

    // segment intersection
    check_intersection<S>  (Rec(p( 0,  0), p(3, 3)), T(p( -10,  0), p(  0,  0), p(  0, 3)));
    check_intersection<S>  (Rec(p(-2,  2), p(3, 3)), T(p( -15, 12), p( -1,  3), p(  2, 3)));
    check_intersection<S>  (Rec(p(-2,  2), p(3, 3)), T(p( -15, 12), p( -4,  3), p(  2, 3)));
    check_intersection<S>  (Rec(p(-2,  2), p(3, 3)), T(p( -15, 12), p( -4,  3), p( 15, 3)));

    // triangle intersection
    check_intersection<T>  (Rec(p( 0,  0), p(1, 1)), T(p(  -1,  0), p( -1,  2), p(  2,  2)));
    check_intersection<T>  (Rec(p( 0,  0), p(1, 1)), T(p(  -1,  0), p(  2,  2), p( -1,  2)));
    check_intersection<T>  (Rec(p( 0,  0), p(2, 2)), T(p(   0,  0), p(  1,  0), p(  0,  1)));
    check_intersection<T>  (Rec(p( 0,  0), p(3, 3)), T(p(   1,  1), p(  2,  1), p(  1,  2)));

    // polygon intersection
    check_intersection<Pol>(Rec(p(   0,   0), p(  1,   1)), T(p( -1, -2), p( -1,   2), p( 5,   2)));
    check_intersection<Pol>(Rec(p( 100, 100), p(200, 200)), T(p(150, 50), p(250, 170), p(50, 170)));
  }
  
  void run()
  {
    std::cout << "2D Intersection tests with Kernel: " << typeid(K).name() << std::endl;
    B_P();
    L_L();
    S_S();
    R_R();
    S_R();
    L_R();
    S_L();
    T_T();
    L_T();
    R_T();
    S_T();
    P_T();
    P_P();
    L_Rec();
    R_Rec();
    S_Rec();
    Rec_Rec();
    T_Rec();
  }
};

int main()
{
  Test< CGAL::Simple_cartesian<typename CGAL::internal::Exact_field_selector<void*>::Type > >().run();
  Test< CGAL::Cartesian<double>   >().run();
  Test< CGAL::Homogeneous<double> >().run();
  Test< CGAL::Exact_predicates_inexact_constructions_kernel >().run();
  Test< CGAL::Exact_predicates_exact_constructions_kernel >().run();
}

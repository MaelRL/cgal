// This is based on the paper :
// Konukoglu et al., A recursive anisotropic fast marching approach to reaction diffusion equation

#include <CGAL/Starset.h>
#include <CGAL/Stretched_Delaunay_2.h>

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <Eigen/Dense>
#include <omp.h>
#include <boost/array.hpp>

#include <iostream>
#include <limits>
#include <set>
#include <vector>

using namespace CGAL::Anisotropic_mesh_2;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef typename K::FT                                       FT;
typedef Stretched_Delaunay_2<K>                              Star;
typedef typename Star::Point_2                               Point_2;
typedef typename Star::Metric                                Metric;
typedef typename Star::Traits                                Traits;

typedef std::set<std::size_t>                                Simplex;
typedef boost::array<std::size_t, 2>                         Edge;
typedef boost::array<std::size_t, 3>                         Tri;

typedef Eigen::Matrix<double, 2, 1>                          Vector2d;

//typedef CGAL::Anisotropic_mesh_2::Euclidean_metric_field<K>* MF;
typedef CGAL::Anisotropic_mesh_2::Custom_metric_field<K>* MF;

typedef CGAL::Exact_predicates_exact_constructions_kernel    KExact;
typedef typename KExact::Point_2                             EPoint;
typedef typename KExact::Segment_2                           ESegment;
typedef typename KExact::Line_2                              ELine;
typedef CGAL::Cartesian_converter<K, KExact>                 To_exact;
typedef CGAL::Cartesian_converter<KExact, K>                 Back_from_exact;

To_exact to_exact;
Back_from_exact back_from_exact;

#define FILTER_SEEDS_OUTSIDE_GRID
#define USE_RECURSIVE_UPDATES
//#define TMP_REFINEMENT_UGLY_HACK

#define verbose 10
const FT FT_inf = std::numeric_limits<FT>::infinity();

// using a lot of global variables, it's ugly but passing them by function is tedious
std::size_t vertices_nv = 1;

// grid related stuff
Point_2 center(0.,0.);
const FT grid_side = 5.0;
FT offset_x = center.x() - grid_side/2.; // offset is the bottom left point
FT offset_y = center.y() - grid_side/2.;
FT n = 100.; // number of points per side of the grid
FT step = grid_side / (n-1);

// the metric field and the seeds
MF mf;
std::vector<Point_2> seeds;
std::vector<Metric> seeds_m;

// how big of an update the new value needs to be
FT recursive_tolerance = 1e-15;

int n_refine = 0;
#ifdef TMP_REFINEMENT_UGLY_HACK
FT best_ref_x = FT_inf, best_ref_y = FT_inf;
#endif

std::clock_t start;

enum FMM_state
{
  CHANGED = 0, // changed implies known. It is used to not add elements that are already in the queue
  KNOWN,
  TRIAL,
  FAR
};

enum PQ_state
{
  NOTHING_TO_DO = 0,
  REBUILD_TRIAL,
  REBUILD_CHANGED,
  REBUILD_BOTH
};

template<typename Gp>
struct Grid_point_comparer
{
  bool operator()(Gp const * const gp1, Gp const * const gp2)
  {
    return gp1->distance_to_closest_seed > gp2->distance_to_closest_seed;
  }
};

int insert_new_seed(const FT x, const FT y)
{
#ifdef FILTER_SEEDS_OUTSIDE_GRID
    if(x < offset_x || x > center.x() + grid_side/2. ||
       y < offset_y || y > center.y() + grid_side/2. )
    {
#if (verbose > 0)
      std::cout << "filtered : " << x << " " << y << std::endl;
#endif
      return seeds.size();
    }
#endif
#if (verbose > 0)
  std::cout << "new seed: " << x << " " << y << " added" << std::endl;
#endif
  seeds.push_back(Point_2(x, y));
  seeds_m.push_back(mf->compute_metric(seeds.back()));
  return seeds.size();
}

int build_seeds()
{
  std::ifstream in("bambimboum.mesh");
  std::string word;
  std::size_t useless, nv, dim;
  FT r_x, r_y;

  in >> word >> useless; //MeshVersionFormatted i
  in >> word >> dim; //Dimension d
  in >> word >> nv;
  std::cout << "nv: " << nv << std::endl;
  assert(dim == 2);

  std::size_t min_nv = (std::min)(nv, vertices_nv);

  seeds.reserve(min_nv);
  seeds_m.reserve(min_nv);

  for(std::size_t i=0; i<nv; ++i)
  {
    in >> r_x >> r_y >> useless;
    insert_new_seed(r_x, r_y);
    if(seeds.size() >= vertices_nv)
      break;
  }
#if (verbose > 0)
  std::cout << "seeds: " << seeds.size() << std::endl;
#endif
  return seeds.size();
}

struct Grid_point
{
  typedef boost::array<Grid_point*, 4> Neighbors;
  typedef boost::array<const Grid_point*, 2> Optimal_neighbors;

  Point_2 point;
  std::size_t index;
  FT distance_to_closest_seed;
  std::size_t closest_seed_id;
  FMM_state state;
  Neighbors neighbors; // 4 neighbors ORDERED (left, up, right, down)
                       // if no neighbor (borders for example) then NULL
  Metric metric;

  const Grid_point* ancestor; // used to debug & draw nice things
  // the neighbors of this for which we have reached the lowest distance_to_closest_seed
  // this is a 2-array with at least 1 non null pointer in [0]
  Optimal_neighbors opt_neighbors;

  bool compute_closest_seed_1D(const Grid_point* gp)
  {
#if (verbose > 15)
    std::cout << "compute closest 1D " << index << " to " << gp->index;
    std::cout << " ds " << gp->distance_to_closest_seed;
    std::cout << " resp. " << gp->closest_seed_id << std::endl;
#endif

    CGAL_assertion(gp->state == KNOWN || gp->state == CHANGED);
    bool changed = false;

    Vector2d v;
    v(0) = gp->point.x() - point.x();
    v(1) = gp->point.y() - point.y();
    const Eigen::Matrix2d& m = metric.get_mat();
    FT neighbor_d = std::sqrt(v.transpose() * m * v);
    FT dcs_at_gp = gp->distance_to_closest_seed;
    FT d = dcs_at_gp + neighbor_d;

#if (verbose > 10)
    std::cout << "min 1D; min : " << d;
    std::cout << " vs current best: " << distance_to_closest_seed << std::endl;
#endif

    if(distance_to_closest_seed - d > recursive_tolerance * d)
    {
      changed = true;
      distance_to_closest_seed = d;

      opt_neighbors[0] = gp;
      opt_neighbors[1] = NULL; // just in case a _1D becomes better than a _2D

      // below is useless fixme
      closest_seed_id = gp->closest_seed_id;
      ancestor = gp;

#if (verbose > 10)
      std::cout << "1D new best for " << index << " : " << distance_to_closest_seed << std::endl;
#endif
    }
    return changed;
  }

  FT eval_at_p(FT p, const FT a, const FT b, const FT c, const FT d, const FT e) const
  {
    if(std::abs(c*p*p + d*p + e) < 1e-16) // avoids stuff like -1e-17
      return p*a+(1.-p)*b;
    CGAL_assertion(c*p*p + d*p + e >= 0);
    return p*a+(1.-p)*b + std::sqrt(c*p*p + d*p + e);
  }

  FT compute_min_distance_2D(const Grid_point* gp, const Grid_point* gq,
                             FT& p) const
  {
    CGAL_assertion((gp->state == KNOWN || gp->state == CHANGED) &&
                   (gq->state == KNOWN || gq->state == CHANGED) );

    // in 2D, we can solve the minimization analytically.
    // We're solving :
    // min_{p\in [0,1]} f(p) = F(x)*p + F(y)*(1-p) + sqrt(v(p)^t * M * v(p))
    // vp = p*gp+(1-p)*gq
    // local minimum is given by f'(p_0) = 0 _but_ we have to compare f(p_0)
    // with f(0) and f(1) as we're only interested in solutions that live in [0,1]
    // If there's no solution in [0,1], then there's no (local) minimum in [0,1]
    // and the minimum is obtained at 0 or at 1.

    const FT dcs_at_gp = gp->distance_to_closest_seed;
    const FT dcs_at_gq = gq->distance_to_closest_seed;
    Vector2d vp, vq;
    vp(0) = gp->point.x() - point.x();
    vp(1) = gp->point.y() - point.y();
    vq(0) = gq->point.x() - point.x();
    vq(1) = gq->point.y() - point.y();

    const Eigen::Matrix2d& m = metric.get_mat();
    const FT vpmvp = vp.transpose() * m * vp;
    const FT vpmvq = vp.transpose() * m * vq;
    const FT vqmvq = vq.transpose() * m * vq;

    // expressing f(p) as ap+(1-p)b+sqrt(cp^2+dp+e)
    // the tiny trick is to write f'(p) = 0 and then square the expression to get rid of sqrt

    const FT a = dcs_at_gp;
    const FT b = dcs_at_gq;
    const FT c = vpmvp + vqmvq - 2.*vpmvq;
    const FT d = 2.*(vpmvq - vqmvq);
    const FT e = vqmvq;

    const FT A = 4*(a-b)*(a-b)*c-4*c*c;
    const FT B = 4*(a-b)*(a-b)*d-4*c*d;
    const FT C = 4*(a-b)*(a-b)*e-d*d;
    FT delta = B*B-4*A*C;

    if(delta >= -1e17)
      delta = (std::max)(0., delta);

#if (verbose > 20)
    std::cout << "a: " << a << std::endl;
    std::cout << "b: " << b << std::endl;
    std::cout << "c: " << c << std::endl;
    std::cout << "d: " << d << std::endl;
    std::cout << "e: " << e << std::endl;
    std::cout << "A: " << A << std::endl;
    std::cout << "B: " << B << std::endl;
    std::cout << "C: " << C << std::endl;
    std::cout << "delta: " << delta << std::endl;
#endif

    FT min_id = 0.;
    FT min = eval_at_p(0., a, b, c, d, e);

    FT d1 = eval_at_p(1., a, b, c, d, e);
    if(min > d1)
    {
      min_id = 1.;
      min = d1;
    }

    // get solutions
    if(A == 0 && B == 0)
    {
      std::cout << "Warning: constant derivative in compute_min_distance_2D : " << C << std::endl;
      // min is found at {0} or {1} (computed above)
    }
    else if(A == 0.)
    {
      p = -C/B;
      FT dp = eval_at_p(p, a, b, c, d, e);
      if(min > d && p >= 0 && p <= 1.)
      {
        min_id = p;
        min = dp;
      }
    }
    else if(delta >= 0)
    {
      FT sqrt_delta = std::sqrt(delta);
      FT denom = 1. / (2.*A);
      double p1 = (-B + sqrt_delta) * denom;
      double p2 = (-B - sqrt_delta) * denom;

      FT dp1 = eval_at_p(p1, a, b, c, d, e);
      FT dp2 = eval_at_p(p2, a, b, c, d, e);

//      std::cout << "p1 : " << p1 << " dp1: " << dp1 << std::endl;
//      std::cout << "p2 : " << p2 << " dp2: " << dp2 << std::endl;

      if(min > dp1 && p1 >= 0 && p1 <= 1.)
      {
        min_id = p1;
        min = dp1;
      }
      if(min > dp2 && p2 >= 0 && p2 <= 1.)
      {
        min_id = p2;
        min = dp2;
      }
    }

    p = min_id;
//    std::cout << "final min at : " << p << " d: " << min << std::endl;
    return min;
  }

  void determine_ancestor(FT lambda = -1.)
  {
#if (verbose > 10)
    std::cout << "determine ancestor for " << index << " d: ";
    std::cout << distance_to_closest_seed << std::endl;
#endif

    if(this->closest_seed_id != static_cast<std::size_t>(-1)) // already has a color
      return;

    CGAL_assertion(opt_neighbors[0]);
    const Grid_point* gp = opt_neighbors[0];
    const Grid_point* gq = opt_neighbors[1];

    bool is_init_p = gp->closest_seed_id != static_cast<std::size_t>(-1);
    bool is_init_q = gq?(gq->closest_seed_id != static_cast<std::size_t>(-1)):false;
    CGAL_assertion((is_init_p || is_init_q) && "gp & gq have uninitialized seeds...");

    if(gp && gq)
    compute_min_distance_2D(gp, gq, lambda);

    if(lambda == 0) // not sure this is correct
    {
      closest_seed_id = gq->closest_seed_id;
      ancestor = gq;
      return;
    }
    else if(lambda == 1)
    {
      closest_seed_id = gp->closest_seed_id;
      ancestor = gp;
      return;
    }

    // quick filters
    if(is_init_p && !is_init_q)
    {
      closest_seed_id = gp->closest_seed_id;
      ancestor = gp;
      return;
    }
    else if(!is_init_p && is_init_q)
    {
      closest_seed_id = gq->closest_seed_id;
      ancestor = gq;
      return;
    }

    bool same_colors = (gp->closest_seed_id == gq->closest_seed_id);

#if (verbose > 20)
    std::cout << "colors are the same ? " << same_colors << std::endl;
#endif

    if(same_colors)
    {
      closest_seed_id = gp->closest_seed_id;
      ancestor = gp;
      return;
    }

    FT gp_d = gp->distance_to_closest_seed;
    FT gq_d = gq->distance_to_closest_seed;
    const Grid_point* best_g = NULL;
#if 0
    // compute p_mid, the limit value to decide which seed id to give depending on p
    // p_mid defined as the point of [pq] on the dual..

    FT Cp, Cq;
    if(gq->index == index + n || gq->index == index - n)
    {
      fixme you have to use the proper Cp & Cq, also verify the computations for
      a non uniform metric field
      Cp = std::sqrt(5.); // m_1 + 2*m_2 + m_3
      Cq = std::sqrt(5.);
    }
    if(gq->index == index + 1 || gq->index == index - 1)
    {
      Cp = std::sqrt(5.); // m_1 - 2*m_2 + m_3
      Cq = std::sqrt(5.);
    }

    FT lambda_mid = (std::sqrt(2.)*step*Cq + gq_d - gp_d)/(Cp + Cq);
    lambda_mid /= step; // to get it in [0,1] if it's in [0, step]
    std::cout << "lambda_mid: " << lambda_mid << std::endl;

    CGAL_assertion(lambda_mid >= -1e-10 && lambda_mid <= 1+1e-10);

    // the point on [pq] is given by p + lambda_mid*(q-p)
    // if we find another point on the dual close, we can approximate the dual
    // by a line and compute the intersection of the dual with 'this p' or
    // 'this q' and deduce in which voronoi cell 'this' is...

    // consider something like below (scaling of the figure is... anisotropic AHAH)

    // THIS ----p'--- p
    //  |     /      /
    //  |    /      /
    //  |   / <----/----- some lambda_mid' living on [p'q'] --> r'
    //  |  /      /
    //  | /      /
    //  q'      /
    //  |      /  lambda_mid around here gives a point r
    //  |     /
    //  |    /
    //  |   /
    //  |  /
    //  | /
    //  |/
    //  q
    // and the dual is given by

    FT rx = gp->point.x() + lambda_mid * (gq->point.x() - gp->point.x());
    FT ry = gp->point.y() + lambda_mid * (gq->point.y() - gp->point.y());
    Point_2 r(rx, ry);

    // let's compute that r' fellow... we need p' and q' first...
    FT rppx = 0.5*(point.x() + gp->point.x());
    FT rppy = 0.5*(point.y() + gp->point.y());
    FT rqpx = 0.5*(point.x() + gq->point.x());
    FT rqpy = 0.5*(point.y() + gq->point.y());
    Point_2 pp(rppx, rppy), qp(rqpx, rqpy);

    // finding lambda_mid' is solving the same thing than finding lambda_prime
    // the value at 'this' is d, the values at p and q are gp_d, gq_d
    FT Cpp, Cqp;

    // can use the index of gp or gq since the order doesn't change
    if(gq->index == index + n || gq->index == index - n)
    {
      Cpp = std::sqrt(5.); // m_1 + 2*m_2 + m_3
      Cqp = std::sqrt(5.);
    }
    if(gq->index == index + 1 || gq->index == index - 1)
    {
      Cpp = std::sqrt(5.); // m_1 - 2*m_2 + m_3
      Cqp = std::sqrt(5.);
    }

    // interpolate the values for value at pp and qp:
    FT gpp_d = 0.5*(gp_d + d);
    FT gqp_d = 0.5*(gq_d + d);

    FT semi_step = 0.5*step;
    FT lambda_midp = (std::sqrt(2.)*semi_step*Cqp + gqp_d - gpp_d)/(Cpp + Cqp);
    lambda_midp /= semi_step; // to get it in [0,1] if it's in [0, step]

    // create r'
    FT rpx = pp.x() + lambda_midp * (qp.x() - pp.x());
    FT rpy = pp.y() + lambda_midp * (qp.y() - pp.y());
    Point_2 rp(rpx, rpy);

    // compute the intersection of the line [r,rp] with [this, gp] and [this, gq]
    const typename K::Segment_2 thisp(point, gp->point);
    const typename K::Segment_2 thisq(point, gq->point);
    const typename K::Segment_2 rrp(r, rp);
    ESegment ethisp = to_exact(thisp);
    ESegment ethisq = to_exact(thisq);
    ELine errp_l = (to_exact(rrp)).supporting_line();

    CGAL::cpp11::result_of<KExact::Intersect_2(ESegment, ELine)>::type
        resultp = CGAL::intersection(ethisp, errp_l);
    CGAL::cpp11::result_of<KExact::Intersect_2(ESegment, ELine)>::type
        resultq = CGAL::intersection(ethisq, errp_l);

    if(lambda_mid < 1e-10 || lambda_mid > 1-1e-10)
      std::cout << "WARNING: CLOSE TO DEGEN CASE IN LAMBDA_MID'" << std::endl;

    const EPoint* res;
    CGAL_assertion((resultp && (res = boost::get<EPoint>(&*resultp))) ||
                   (resultq && (res = boost::get<EPoint>(&*resultq))));

    // handle the degen case of intersection being 'this'... todo

    if (resultp)
    {
      std::cout << "intersected on [this gp]" << std::endl;
      best_g = gq;
    }
    if (resultq)
    {
      std::cout << "intersected on 'this gq'" << std::endl;
      best_g = gp;
    }
#else
    // Version n°2 : we compute the intersection of the dual with [this,p] and
    // [this,q] (we know that the dual intersects p-q).
    // depending on which side is intersected, we can deduce the color of 'this'

    // fixme, rename max & min in metric.h and everywhere...
    // we never care who's max and who's min, just their position !

    // we're solving d(this) + d_this(this-mp) = d(p) + d_p(p-mp)

    const Eigen::Matrix2d& mp = gp->metric.get_mat();
    const Eigen::Matrix2d& mq = gq->metric.get_mat();
    const Eigen::Matrix2d& m = metric.get_mat();

    FT pe1 = std::sqrt(mp(0,0)), pe2 = std::sqrt(mp(1,1));
    FT qe1 = std::sqrt(mq(0,0)), qe2 = std::sqrt(mq(1,1));
    FT me1 = std::sqrt(m(0,0)), me2 = std::sqrt(m(1,1));

    FT lambda_p = -1., lambda_q = -1.;
    if(gp->index == index + 1 || gp->index == index - 1)
    {
      // [this gp] is horizontal and [this,gq] is vertical
      // compute the intersection on [this,gp] & [this gq]
      lambda_p = (gp_d - distance_to_closest_seed + pe1*step)/(me1 + pe1);
      lambda_q = (gq_d - distance_to_closest_seed + qe2*step)/(me2 + qe2);
    }
    else if(gp->index == index + n || gp->index == index - n)
    {
      // 'this gp' is vertical and [this,gq] is horizontal
      // compute the intersection on [this,gp] & [this gq]
      lambda_p = (gp_d - distance_to_closest_seed + pe2*step)/(me2 + pe2);
      lambda_q = (gq_d - distance_to_closest_seed + qe1*step)/(me1 + qe1);
    }

#ifdef DEBUG_LAMBDAS
    // check that we're computing the correct lambdas :
    FT mpx = point.x() + lambda_p * (gp->point.x() - point.x()) / step;
    FT mpy = point.y() + lambda_p * (gp->point.y() - point.y()) / step;
    FT mqx = point.x() + lambda_q * (gq->point.x() - point.x()) / step;
    FT mqy = point.y() + lambda_q * (gq->point.y() - point.y()) / step;

    Vector2d p_to_mp, this_to_mp, q_to_mq, this_to_mq;
    //these are flipped, but it doesn't matter since v^t M v = -v^t M -v
    p_to_mp(0) = gp->point.x() - mpx;
    p_to_mp(1) = gp->point.y() - mpy;
    q_to_mq(0) = gq->point.x() - mqx;
    q_to_mq(1) = gq->point.y() - mqy;
    this_to_mp(0) = point.x() - mpx;
    this_to_mp(1) = point.y() - mpy;
    this_to_mq(0) = point.x() - mqx;
    this_to_mq(1) = point.y() - mqy;

    FT d_p_to_mp = std::sqrt(p_to_mp.transpose()*mp*p_to_mp);
    FT d_q_to_mq = std::sqrt(q_to_mq.transpose()*mq*q_to_mq);
    FT d_this_to_mp = std::sqrt(this_to_mp.transpose()*m*this_to_mp);
    FT d_this_to_mq = std::sqrt(this_to_mq.transpose()*m*this_to_mq);

    FT diffp = distance_to_closest_seed + d_this_to_mp - (gp_d + d_p_to_mp);
    FT diffq = distance_to_closest_seed + d_this_to_mq - (gq_d + d_q_to_mq);

    // we simplified the absolute values of the system because we need lamba in [0, step]
    // the assertion will be false if the solution we found is outside of this interval
    if(lambda_p >= -1e-10 && lambda_p <= step+1e-10)
      CGAL_assertion(std::abs(diffp) < 1e-10 );
    if(lambda_q >= -1e-10 && lambda_q <= step+1e-10)
      CGAL_assertion(std::abs(diffq) < 1e-10 );
#endif

#if (verbose > 20)
    std::cout << "lambdas: " << lambda_p << " " << lambda_q << std::endl;
    std::cout << "step: " << step << " gpd/q: " << gp_d << " " << gq_d << std::endl;
#endif

    // if the colors are different, then the dual intersects the segment [p,q]
    // if we find which of the segment [this p] or [this q] intersects the dual
    // we know which colors to associate to 'this'!

    bool intersect_p = lambda_p >= -1e-2 && lambda_p <= step+1e-10;
    bool intersect_q = lambda_q >= -1e-2 && lambda_q <= step+1e-10;

    if(intersect_p)
    {
      if(intersect_q)
      {
        // we intersect both [this p] and [this q]...
#if (verbose > 5)
        std::cerr << "mega warning [this p] and [this q]" << std::endl;
#endif
        best_g = (lambda_q > lambda_p)?gp:gq; // this is shady... fixme...?
      }
      else
        best_g = gq; // we intersect [this,p] therefore 'this' is on the side of q
    }
    else if(intersect_q)
      best_g = gp; // we intersect [this q] therefore 'this' is on the side of p
    else
      CGAL_assertion(false && "we intersect nothing... ?");
#endif

    closest_seed_id = best_g->closest_seed_id;
    ancestor = best_g;

#if (verbose > 5)
    std::cout << "ancestor shen @ " << index << " color is : " << closest_seed_id << std::endl;
#endif
  }

  bool compute_closest_seed_2D(const Grid_point* gp, const Grid_point* gq)
  {
    FT gp_d = gp->distance_to_closest_seed;
    FT gq_d = gq->distance_to_closest_seed;

#if (verbose > 15)
    std::cout << "compute closest 2D " << index << " to ";
    std::cout << gp->index << " & " << gq->index;
    std::cout << " ds: " << gp_d << " " << gq_d << " ";
    std::cout << " resp. " << gp->closest_seed_id << " " << gq->closest_seed_id << std::endl;
#endif

    CGAL_assertion((gp->state == KNOWN || gp->state == CHANGED) &&
                   (gq->state == KNOWN || gq->state == CHANGED) );
    bool changed = false;

    FT lambda = 0.;
    const FT d = compute_min_distance_2D(gp, gq, lambda);
    CGAL_assertion(lambda >= -1e-10 && lambda <= 1+1e-10);

#if (verbose > 10)
    std::cout << "min 2D " << d << " vs current best: " << distance_to_closest_seed << std::endl;
#endif

#ifdef BRUTE_FORCE_CHECK_OPTIMAL_P
    // this is pretty much debug code to verify that compute_min_distance_2D() is correct
    std::size_t query_inters = 1e5;
    FT query_step = 1./static_cast<FT>(query_inters);
    FT min_lambda;
    FT min_d = FT_inf;

    for(std::size_t i=0; i<=query_inters; ++i)
    {
      // segment PQ
      const FT x = query_step*(i*gp->point.x() + (query_inters-i)*gq->point.x());
      const FT y = query_step*(i*gp->point.y() + (query_inters-i)*gq->point.y());

      // arc PQ (just to check what would be the value if we interpolated on a cicle arc
      // rather than a segment)

      // pq has been cleverly (not...) chosen to be in counter clockwise angle
      // thus if we do min_angle + add_angle, we must start from gq and do
      // a pi/2 increase...
//      FT min_angle;
//      if(gq->index == index-1)
//        min_angle = M_PI;
//      else if(gq->index == index + n)
//        min_angle = M_PI_2;
//      else if(gq->index == index + 1)
//        min_angle = 0;
//      else // (gq->index == index - n)
//        min_angle = -M_PI_2;
//      const FT add_angle = query_step*i*M_PI_2;
//      const FT x = point.x() + step * std::cos(min_angle + add_angle);
//      const FT y = point.y() + step * std::sin(min_angle + add_angle);

      const FT d_at_que = query_step*(i*gp->distance_to_closest_seed +
                                     (query_inters-i)*gq->distance_to_closest_seed);
      Vector2d v;
      v(0) = point.x() - x;
      v(1) = point.y() - y;
      const Eigen::Matrix2d& m = metric.get_mat();
      FT neighbor_d = std::sqrt(v.transpose() * m * v);
      FT d_que = d_at_que + neighbor_d;
      if(d_que < min_d)
      {
        min_lambda = i*query_step;
        min_d = d_que;
      }
    }
    std::cout << query_step << std::endl;
    std::cout << "lambda, min_lambda: " << lambda << " " << min_lambda << std::endl;
    std::cout << ":min_d d " << min_d << " " << d << std::endl;
    CGAL_assertion(std::abs(min_d - d) < 1e-2*d);
    CGAL_assertion(std::abs(min_lambda - lambda) < query_step);
#endif
    // -------------------------------------------------------------------------

    // tolerance to ignore unsignificant changes
    if(distance_to_closest_seed - d > recursive_tolerance * d)
    {
      changed = true;
      distance_to_closest_seed = d;

      opt_neighbors[0] = gp;
      opt_neighbors[1] = gq;

      // below is useless fixme
      closest_seed_id = gp->closest_seed_id;
      ancestor = gp;

#if (verbose > 10)
      std::cout << "2D new best for " << index << " : " << distance_to_closest_seed << std::endl;
#endif
    }

    return changed;
  }

  bool compute_closest_seed(const Grid_point* anc,
                            const bool ignore_ancestor = false)
  {
#if (verbose > 10)
    std::cout << "compute closest high level " << index;
    std::cout << " state : " << state;
    std::cout << " dist : " << distance_to_closest_seed;
    std::cout << " color : " << closest_seed_id;
    std::cout << " ancestor: " << anc->index << std::endl;
#endif

    bool changed = false;
    for(std::size_t i=0; i<neighbors.size(); ++i)
    {
      const Grid_point* gp = neighbors[i];
      const Grid_point* gq = neighbors[(i+1)%(neighbors.size())];

      if(gp && (gp->state == KNOWN || gp->state == CHANGED))
      {
        if(gq && (gq->state == KNOWN || gq->state == CHANGED))
        {
          if(ignore_ancestor || gp == anc || gq == anc)
            if(compute_closest_seed_2D(gp, gq))
              changed = true;
        }
        else // we could do it++ here since we know the next 'it' has gp->state != KNOWN
        {
          if(ignore_ancestor || gp == anc)
            if(compute_closest_seed_1D(gp))
              changed = true;
        }
      }
    }

#if (verbose > 10)
    std::cout << "color and value: " << closest_seed_id << " " << distance_to_closest_seed << std::endl;
#endif
    return changed;
  }

  PQ_state update_neighbors_distances(std::vector<Grid_point*>& changed_pq,
                                      std::vector<Grid_point*>& trial_pq) const
  {
#if (verbose > 10)
    std::cout << "update neighbors of " << index << std::endl;
#endif

    PQ_state pqs_ret = NOTHING_TO_DO;

    // consider all the known neighbors and compute the value to that
    CGAL_assertion(state == KNOWN || state == CHANGED);

    for(std::size_t i=0; i<neighbors.size(); ++i)
    {
      Grid_point* gp = neighbors[i]; // a neighbor of 'this'
      if(!gp)
        continue;

#ifdef USE_RECURSIVE_UPDATES
      if(gp->state == KNOWN || gp->state == CHANGED) // that's the recursive part
      {
        // recompute the distance in the direction gp-'this'
        FT mem = gp->distance_to_closest_seed;
        std::size_t ancestor_mem = (gp->ancestor)?gp->ancestor->index:-1;
        if(gp->compute_closest_seed(this))
        {
#if (verbose > 15)
          std::cout << "new value at " << gp->index << " : " << gp->distance_to_closest_seed;
          std::cout << " with ancestor: " << gp->ancestor->index;
          std::cout << " (mem: " << mem << " ancestor: " << ancestor_mem << ") ";
          std::cout << " diff: " << gp->distance_to_closest_seed-mem << std::endl;
#endif
          if(gp->state == KNOWN)
          {
            // if it's already CHANGED we only need to reorder the CHANGED queue
            gp->state = CHANGED;
            changed_pq.push_back(gp);
            std::push_heap(changed_pq.begin(), changed_pq.end(), Grid_point_comparer<Grid_point>());
          }
          else
          {
            if(pqs_ret == NOTHING_TO_DO || pqs_ret == REBUILD_CHANGED)
              pqs_ret = REBUILD_CHANGED;
            else // pqs_ret == REBUILD_TRIAL || pqs_ret == REBUILD_BOTH
              pqs_ret = REBUILD_BOTH;
          }
        }
      }
      else
#endif
      if(gp->state == TRIAL)
      {
        // note that we don't insert in trial_pq since it's already in
        if(gp->compute_closest_seed(this))
        {
          if(pqs_ret == NOTHING_TO_DO || pqs_ret == REBUILD_TRIAL)
            pqs_ret = REBUILD_TRIAL;
          else // pqs_ret == REBUILD_CHANGED || pqs_ret == REBUILD_BOTH
            pqs_ret = REBUILD_BOTH;
        }
      }
      else // gp->state == FAR
      {
        CGAL_assertion(gp->state == FAR);
        if(gp->compute_closest_seed(this))
        {
          gp->state = TRIAL;
          trial_pq.push_back(gp);
          std::push_heap(trial_pq.begin(), trial_pq.end(),
                                             Grid_point_comparer<Grid_point>());
        }
      }
    }

    return pqs_ret;
  }

  bool is_contained(const Point_2& p,
                    const Grid_point* next) const
  {
    //    n1
    //    |
    // n0-gp-n2
    //    |
    //    n3

    Grid_point *n0=neighbors[0], *n1=neighbors[1],
               *n2=neighbors[2], *n3=neighbors[3];

    bool is_p_leq_n0_x = (n0) ? p.x() <= n0->point.x() : false;
    bool is_p_geq_n1_y = (n1) ? p.y() >= n1->point.y() : false;
    bool is_p_geq_n2_x = (n2) ? p.x() >= n2->point.x() : false;
    bool is_p_leq_n3_y = (n3) ? p.y() <= n3->point.y() : false;

    if(is_p_leq_n0_x)
      next = n0;
    else if(is_p_geq_n1_y)
      next = n1;
    else if(is_p_geq_n2_x)
      next = n2;
    else if (is_p_leq_n3_y)
      next = n3;
    else
      next = this; // we're inside !
    return (!is_p_leq_n0_x && !is_p_geq_n1_y && !is_p_geq_n2_x && !is_p_leq_n3_y);
  }

  void initialize_from_point(const Point_2& p,
                             const std::size_t seed_id)
  {
    Vector2d v;
    v(0) = p.x() - point.x();
    v(1) = p.y() - point.y();
    const Eigen::Matrix2d& m = metric.get_mat();
    FT d = std::sqrt(v.transpose() * m * v);
    closest_seed_id = seed_id;
    distance_to_closest_seed = d;
    state = TRIAL;
    ancestor = NULL;
  }

  Grid_point(const Point_2& point_, const std::size_t index_)
    : point(point_), index(index_), distance_to_closest_seed(FT_inf),
      closest_seed_id(-1), state(FAR), neighbors(),
      metric(mf->compute_metric(point)), ancestor(NULL), opt_neighbors()
  {
    for(std::size_t i=0; i<neighbors.size(); ++i)
      neighbors[i] = NULL;
  }
};

struct Geo_grid
{
  typedef std::vector<Grid_point*> Grid_point_vector;

  std::vector<Grid_point> points;
  Grid_point_vector changed_points;
  Grid_point_vector trial_points;

  const Grid_point* locate_point(const Point_2& p) const
  {
    bool found = false;
    const Grid_point& next = points[0];
    while(!found)
      found = (&next)->is_contained(p, &next);
    return &next;
  }

  void locate_and_initialize(const Point_2& p,
                             const std::size_t seed_id)
  {
    // only usable for a full grid
    int index_x = std::floor((p.x()-offset_x)/step);
    int index_y = std::floor((p.y()-offset_y)/step);

    Grid_point* gp = &(points[index_y*n + index_x]);
    if(gp->closest_seed_id != static_cast<std::size_t>(-1))
    {
      std::cerr << "WARNING: a new seed is overwriting the closest seed id";
      std::cerr << " of a grid point! previous index is: " << gp->closest_seed_id << std::endl;
    }

    if(gp->state == TRIAL) // it's already a trial point, we can't accept two seeds for one grid pt
      CGAL_assertion(false && "the grid is not dense enough for the input...");
    gp->initialize_from_point(p, seed_id);

    trial_points.push_back(gp);
    std::push_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());

#if (verbose > 15)
    std::cout << "looking for p: " << p.x() << " " << p.y() << std::endl;
    std::cout << "found: " << gp->index << " [" << gp->point.x() << ", " << gp->point.y() << "]" << std::endl;
    std::cout << "step reminder: " << step << std::endl;
    std::cout << "looking for p: " << p.x() << " " << p.y() << " " << std::endl;
    std::cout << "found gp: " << gp->index << " [" << gp->point.x() << ", " << gp->point.y() << "]" << std::endl;
    std::cout << "index: " << index_x << " " << index_y << " " << std::endl;
    std::cout << "offset: " << offset_x << " " << offset_y << " " << std::endl;
#endif

    // if we want to initialize a bit more, we initialize all the vertices of the quad
    // that contains p
    // is this really useful... ?

/* //fixme these neighbors don't always exist
    CGAL_assertion(gp->neighbors[1] && gp->neighbors[2] && (gp->neighbors[2])->neighbors[1]);
    Grid_point* gq = gp->neighbors[2];
    Grid_point* gr = gq->neighbors[1];
    Grid_point* gs = gp->neighbors[1]; // which is also gr->nei[0]

    gq->initialize_from_point(p, seed_id);
    trial_points.push_back(gq);
    std::push_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());

    gr->initialize_from_point(p, seed_id);
    trial_points.push_back(gr);
    std::push_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());

    gs->initialize_from_point(p, seed_id);
    trial_points.push_back(gs);
    std::push_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());

    std::cout << "neighbors are : " << gq->index << " " << gq->point.x() << " " << gq->point.y() << std::endl;
    std::cout << "neighbors are : " << gr->index << " " << gr->point.x() << " " << gr->point.y() << std::endl;
    std::cout << "neighbors are : " << gs->index << " " << gs->point.x() << " " << gs->point.y() << std::endl;
*/
  }

  void initialize_geo_grid()
  {
    std::cout << "grid initialization" << std::endl;

    // create the grid points (either a full grid or a smart grid (todo) )
    for(unsigned int j=0; j<n; ++j)
    {
      for(unsigned int i=0; i<n; ++i)
      {
        Point_2 p(offset_x+i*step, offset_y+j*step); // fill from bot left to top right
        Grid_point gp(p, i + j*n);
        points.push_back(gp);
      }
    }

    // assign the neighbors
    for(unsigned int j=0; j<n; ++j)
    {
      for(unsigned int i=0; i<n; ++i)
      {
        std::size_t curr_id = i + j*n;
        if(i != 0) // there is a neighbor left
          points[curr_id].neighbors[0] = &(points[i-1 + j*n]);
        if(j != n-1) // there is a neighbor above
          points[curr_id].neighbors[1] = &(points[i + (j+1)*n]);
        if(i != n-1) // there is a neighbor right
          points[curr_id].neighbors[2] = &(points[i+1 + j*n]);
        if(j != 0) // there is a neighbor below
          points[curr_id].neighbors[3] = &(points[i + (j-1)*n]);
      }
    }
    std::cout << "neighbors assigned" << std::endl;

    // seeds belong to a quad, find it
    for(std::size_t i=0; i<seeds.size(); ++i)
    {
      const Point_2& p = seeds[i];
//      const Grid_point* gp = locate_point(p); // todo
//      gp->initialize_from_point(p, i);
      locate_and_initialize(p, i);
    }

    std::cout << "grid initialized" << std::endl;
  }

  bool print_states() const
  {
    int changed_count=0, known_count=0, trial_count=0, far_count=0;

    for(std::size_t i=0; i<points.size(); ++i)
    {
      const Grid_point& gp = points[i];
      if(gp.state == CHANGED)
        changed_count++;
      if(gp.state == KNOWN)
        known_count++;
      if(gp.state == TRIAL)
        trial_count++;
      if(gp.state == FAR)
        far_count++;
    }
    std::cout << "changed: " << changed_count << " known: " << known_count;
    std::cout << " trial: " << trial_count << " far: " << far_count << std::endl;

    return false;
  }

  void determine_ancestors()
  {
#if (verbose > 5)
    std::cout << "Determine ancestors" << std::endl;
#endif

    // TODO do something smarter when spreading a new color in an already colored
    // grid... One can probably just start from the new seed WITHOUT RESETING EVERYTHING
    // and then spread and only consider stuff if one of the possible ancestors
    // has the new color something like that...

    Grid_point_vector known_points;
    for(std::size_t i=0; i<points.size(); ++i)
    {
      // reset the colors (but not the values!)
      Grid_point* gp = &(points[i]);
      gp->closest_seed_id = -1;
      gp->ancestor = NULL;

      known_points.push_back(gp);
    }
    std::make_heap(known_points.begin(), known_points.end(), Grid_point_comparer<Grid_point>());

    // only re initialize the closest_seed_id for each point-seed
    for(std::size_t i=0; i<seeds.size(); ++i)
    {
      // if you use 4 pts initilization  for a seed, you need to modify here too...
      const Point_2& p = seeds[i];
      int index_x = std::floor((p.x()-offset_x)/step);
      int index_y = std::floor((p.y()-offset_y)/step);
      Grid_point* gp = &(points[index_y*n + index_x]);
      gp->closest_seed_id = i;
    }

    bool is_kp_empty = known_points.empty();
    while(!is_kp_empty)
    {
      Grid_point* gp = known_points.front();
      std::pop_heap(known_points.begin(), known_points.end(), Grid_point_comparer<Grid_point>());
      known_points.pop_back();

      gp->determine_ancestor();
      is_kp_empty = known_points.empty();
    }

    std::cerr << "End of determine ancestors. time: ";
    std::cerr << ( std::clock() - start ) / (double) CLOCKS_PER_SEC << std::endl;
  }

  void geo_grid_loop()
  {
    std::cout << "main loop" << std::endl;

    bool is_cp_empty = changed_points.empty();
    bool is_t_empty = trial_points.empty();

    Grid_point* gp;
    while(!is_cp_empty || !is_t_empty)
    {
#if (verbose > 5)
      std::cout << "Queue sizes. Trial: " << trial_points.size() << " Changed: " << changed_points.size() << std::endl;
#endif
#if (verbose > 15)
      std::cout << "changed heap: " << std::endl;
      for (std::vector<Grid_point*>::iterator it = changed_points.begin();
           it != changed_points.end(); ++it)
        std::cout << (*it)->index << " " << (*it)->distance_to_closest_seed << std::endl;
      std::cout << std::endl;

      std::cout << "trial heap: " << std::endl;
      for (std::vector<Grid_point*>::iterator it = trial_points.begin();
           it != trial_points.end(); ++it)
        std::cout << (*it)->index << " " << (*it)->distance_to_closest_seed << std::endl;
      std::cout << std::endl;
#endif

      if(!is_cp_empty)
      {
        gp = changed_points.front();
        CGAL_assertion(gp && gp->state == CHANGED);
        std::pop_heap(changed_points.begin(), changed_points.end(), Grid_point_comparer<Grid_point>());
        changed_points.pop_back();
      }
      else // !is_t_empty
      {
        gp = trial_points.front();
        CGAL_assertion(gp && gp->state == TRIAL);
        std::pop_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());
        trial_points.pop_back();
      }

//      if(print_states())
//        break;

#if (verbose > 10)
      std::cout << "picked n° " << gp->index << " (" << gp->point.x() << ", " << gp->point.y() << ") ";
      std::cout << "at distance : " << gp->distance_to_closest_seed << " from " << gp->closest_seed_id << std::endl;
#endif

      gp->state = KNOWN;
      PQ_state pqs = gp->update_neighbors_distances(changed_points, trial_points);

      if(pqs == REBUILD_TRIAL || pqs == REBUILD_BOTH)
        std::make_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());
      if(pqs == REBUILD_CHANGED || pqs == REBUILD_BOTH)
        std::make_heap(changed_points.begin(), changed_points.end(), Grid_point_comparer<Grid_point>());

      is_cp_empty = changed_points.empty();
      is_t_empty = trial_points.empty();
    }

    std::cerr << "End of geo_grid_loop. time: ";
    std::cerr << ( std::clock() - start ) / (double) CLOCKS_PER_SEC << std::endl;

    determine_ancestors();
    //debug();
  }

  void debug()
  {
#ifdef USE_RECURSIVE_UPDATES
    CGAL_assertion(trial_points.empty() && changed_points.empty());

    std::cout << "BRUTE FORCE CHECK THAT WE ARE REALLY FINISHED : " << std::endl;
    for(std::size_t i=0; i<points.size(); ++i)
    {
      Grid_point* gp = &(points[i]);
      std::cout << "point " << i << " min distance is supposedly: ";
      std::cout << gp->distance_to_closest_seed << std::endl;
      boost::array<Grid_point*, 4>::const_iterator it = gp->neighbors.begin(),
                                                   iend = gp->neighbors.end();
      for(; it!=iend; ++it)
      {
        const Grid_point* gq = *it;
        if(gq)
          CGAL_assertion(!gp->compute_closest_seed(gq));
      }
    }
#endif

    std::cerr << "End of debug. time: ";
    std::cerr << ( std::clock() - start ) / (double) CLOCKS_PER_SEC << std::endl;
  }

  void build_grid()
  {
    initialize_geo_grid();
    geo_grid_loop();
  }

  Point_2 compute_refinement_point()
  {
    // todo
    CGAL_assertion(false);
    return Point_2();
  }

  void refresh_grid_point_states()
  {
    CGAL_assertion(trial_points.empty());
    for(std::size_t i=0; i<points.size(); ++i)
      points[i].state = FAR;
  }

  void refresh_grid_after_seed_creation()
  {
    // todo change this function so that we can insert multiple points at once

    refresh_grid_point_states();

    std::size_t seed_id = seeds.size() - 1;
    const Point_2& p = seeds.back();
    locate_and_initialize(p, seed_id);

    geo_grid_loop();
  }

  void refine_grid()
  {
#ifdef TMP_REFINEMENT_UGLY_HACK
    std::cout << "insert: " << best_ref_x << " " << best_ref_y << std::endl;
    Point_2 p(best_ref_x, best_ref_y);
    CGAL_assertion(best_ref_x != FT_inf && best_ref_y != FT_inf);
#else
    Point_2 p = compute_refinement_point();
#endif
    vertices_nv = insert_new_seed(p.x(), p.y());
    refresh_grid_after_seed_creation();

#ifdef TMP_REFINEMENT_UGLY_HACK
    best_ref_x = FT_inf; best_ref_y = FT_inf;
#endif
  }

  void output_grid(const std::string str_base) const
  {
    std::ofstream out((str_base + ".mesh").c_str());
    std::ofstream out_bb((str_base + ".bb").c_str());

    out << "MeshVersionFormatted 1" << std::endl;
    out << "Dimension 2" << std::endl;
    out << "Vertices" << std::endl;
    out << points.size() << std::endl;
    out_bb << "2 1 " << points.size() << " 2" << std::endl;

    std::set<Tri> simplices;
    for(std::size_t i=0; i<points.size(); ++i)
    {
      const Grid_point& gp = points[i];
      out << gp.point.x() << " " << gp.point.y() << " " << i+1 << std::endl;

//      out_bb << gp.closest_seed_id << std::endl;
      out_bb << gp.distance_to_closest_seed << std::endl;

      const Grid_point* n0 = gp.neighbors[0];
      const Grid_point* n1 = gp.neighbors[1];
      const Grid_point* n2 = gp.neighbors[2];
      const Grid_point* n3 = gp.neighbors[3];

//      CGAL_assertion((!n0 || n0->point.x() == gp.point.x() - step) &&
//                     (!n1 || n1->point.y() == gp.point.y() + step) &&
//                     (!n2 || n2->point.x() == gp.point.x() + step) &&
//                     (!n3 || n3->point.y() == gp.point.y() - step) );

      if(n1 && n2)
      {
        Tri simplex;
        simplex[0] = gp.index;
        simplex[1] = n1->index;
        simplex[2] = n2->index;
        simplices.insert(simplex);
      }
      if(n0 && n3)
      {
        Tri simplex;
        simplex[0] = gp.index;
        simplex[1] = n0->index;
        simplex[2] = n3->index;
        simplices.insert(simplex);
      }

//      // we want to draw a botleft-topright diagonal so we only draw
//      // gp-n0-n1 and gp-n2-n3 (if possible obviously)
//      typename Neighbors::const_iterator it = gp.neighbors.begin(),
//                                         iend = gp.neighbors.end();
//      for(; it!=iend; ++it)
//      {
//        const Grid_point* gq = *it;
//        const Grid_point* gr = *((++it));
//        if(gq && gr)
//        {
//          Simplex simplex;
//          simplex.insert(gp.index);
//          simplex.insert(gq->index);
//          simplex.insert(gr->index);
//          simplices.insert(simplex);
//        }
//      }
    }

    out << "Triangles" << std::endl;
    out << simplices.size() << std::endl;
    for(typename std::set<Tri>::iterator it = simplices.begin();
                                         it != simplices.end(); ++it)
    {
      const Tri& t = *it;
      std::set<std::size_t> materials;

      for(std::size_t i=0; i<t.size(); ++i)
      {
        out << t[i] + 1 << " ";
        materials.insert(points[t[i]].closest_seed_id);
      }
      std::size_t mat = (materials.size()==1)?(*(materials.begin())):(seeds.size());
      out << mat << std::endl;
      //if only one value at all points, we're okay
      // if more than one, print a unique frontier material
    }

    out_bb << "End" << std::endl;
    out << "End" << std::endl;
  }

  void add_triangles_edges_to_dual_edges(const Tri& tr,
                                         std::set<Edge>& dual_edges) const
  {
    Edge e;
    e[0] = tr[0]; e[1] = tr[1]; dual_edges.insert(e);
    e[0] = tr[0]; e[1] = tr[2]; dual_edges.insert(e);
    e[0] = tr[1]; e[1] = tr[2]; dual_edges.insert(e);
  }

  void compute_dual(std::set<Edge>& dual_edges,
                    std::set<Tri>& dual_triangles) const
  {
#if (verbose > 5)
    std::cout << "dual computations" << std::endl;
#endif

#ifdef TMP_REFINEMENT_UGLY_HACK
    FT farthest_dual_distance = 0.;
    const Grid_point* ref;
    FT back_up_farthest = 0.;
    const Grid_point* backup; // separate from ref because it could be farther
#endif

    for(std::size_t i=0; i<points.size(); ++i)
    {
      const Grid_point& gp = points[i];
      Simplex dual_simplex;
      dual_simplex.insert(gp.closest_seed_id);

      boost::array<Grid_point*, 4>::const_iterator it = gp.neighbors.begin(),
          iend = gp.neighbors.end();
      for(; it!=iend; ++it)
      {
        const Grid_point* gq = *it;
        if(!gq)
          continue;

        dual_simplex.insert(gq->closest_seed_id);
      }

      // don't care if dual_simplex.size() == 1
      if(dual_simplex.size() == 2) // an edge!
      {
        typename Simplex::const_iterator it = dual_simplex.begin();
        Edge e; e[0] = *it; e[1] = (*++it);
        dual_edges.insert(e);

#ifdef TMP_REFINEMENT_UGLY_HACK
        if(dual_triangles.empty() && gp.distance_to_closest_seed > back_up_farthest)
        {
          // this is just in case we don't find any triangle...
          back_up_farthest = gp.distance_to_closest_seed;
          backup = &gp;
        }
#endif
      }
      else if(dual_simplex.size() >= 3)
      {
        if(dual_simplex.size() == 3) // a triangle!
        {
          typename Simplex::const_iterator it = dual_simplex.begin();
          Tri tr; tr[0] = *it; tr[1] = (*++it); tr[2] = (*++it);
          dual_triangles.insert(tr);

          add_triangles_edges_to_dual_edges(tr, dual_edges);
        }
        else
          std::cerr << "WARNING: COSPHERICITY..." << std::endl;

#ifdef TMP_REFINEMENT_UGLY_HACK
        if(gp.distance_to_closest_seed > farthest_dual_distance)
        {
          farthest_dual_distance = gp.distance_to_closest_seed;
          ref = &gp;
        }
#endif
      }
    }

#ifdef TMP_REFINEMENT_UGLY_HACK
    if(dual_triangles.empty())
    {
      std::cerr << "WARNING: no triangle captured, using back-up " << std::endl;
      best_ref_x = backup->point.x();
      best_ref_y = backup->point.y();
    }
    else
    {
      best_ref_x = ref->point.x();
      best_ref_y = ref->point.y();
    }
#endif
  }

  void output_dual(const std::string str_base) const
  {
    std::set<Edge> dual_edges;
    std::set<Tri> dual_triangles;
    compute_dual(dual_edges, dual_triangles);

#if (verbose > 5)
    std::cout << "captured: ";
    std::cout << dual_edges.size() << " edges, ";
    std::cout << dual_triangles.size() << " triangles" << std::endl;
#endif

    std::ofstream out((str_base + "_dual.mesh").c_str());
    out << "MeshVersionFormatted 1" << std::endl;
    out << "Dimension 2" << std::endl;
    out << "Vertices" << std::endl;
    out << seeds.size() << std::endl;
    for(std::size_t i=0; i<seeds.size(); ++i)
      out << seeds[i].x() << " " << seeds[i].y() << " " << i+1 << std::endl;

    out << "Triangles" << std::endl;
    out << dual_triangles.size() << std::endl;
    for(typename std::set<Tri>::iterator it = dual_triangles.begin();
                                         it != dual_triangles.end(); ++it)
    {
      const Tri& tr = *it;
      for(std::size_t i=0; i<tr.size(); ++i)
        out << tr[i] + 1 << " ";
      out << "1" << std::endl; // fixme add 'is_intersected' like in grid_gen.cpp
    }
    out << "End" << std::endl;
  }

  Geo_grid() : points(), changed_points(), trial_points() { }
};

void initialize()
{
//  mf = new Euclidean_metric_field<K>(10,1);
  mf = new Custom_metric_field<K>();

  vertices_nv = build_seeds();
  CGAL_assertion(vertices_nv > 0 && "No seed in domain..." );

#ifdef TMP_REFINEMENT_UGLY_HACK
  CGAL_assertion(n_refine > 0 && "don't enable the macro if you're not refining...");
  if(vertices_nv < 2)
    CGAL_assertion(false && "can't refine from a 1 seed diagram...");
#endif
}

int main(int, char**)
{
  std::cout.precision(17);
  std::freopen("log.txt", "w", stdout);

  double duration;
  start = std::clock();

  std::srand(0);
  initialize();

  Geo_grid gg;
  gg.build_grid();
  gg.output_grid("geo_grid_pre");
  gg.output_dual("geo_grid_pre");

  for(int i=0; i<n_refine; ++i)
  {
#ifndef TMP_REFINEMENT_UGLY_HACK
    std::cerr << "you forgot to enable the macro, birdhead" << std::endl;
#endif

    gg.refine_grid();
    std::ostringstream out;
    out << "geo_grid_ref" << i;
    gg.output_grid(out.str());
    gg.output_dual(out.str());
  }

  gg.output_grid("geo_grid");
  gg.output_dual("geo_grid");

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cerr << "dur: " << duration << std::endl;
}

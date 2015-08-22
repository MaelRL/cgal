#include <CGAL/Starset.h>
#include <CGAL/Stretched_Delaunay_3.h>

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <boost/array.hpp>
#include <boost/chrono.hpp>
#include <Eigen/Dense>
#include <omp.h>

#include <ctime>
#include <iostream>
#include <limits>
#include <set>
#include <vector>

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef typename K::FT                                       FT;
typedef typename K::Segment_3                                Segment;
typedef typename K::Vector_3                                 Vector;
typedef typename K::Triangle_3                               Triangle;
typedef typename K::Line_3                                   Line;

typedef Stretched_Delaunay_3<K>                              Star;
typedef typename Star::Point_3                               Point_3;
typedef typename Star::Metric                                Metric;
typedef typename Star::Traits                                Traits;

typedef typename Eigen::Matrix<double, 3, 1>                 Vector3d;

//typedef typename CGAL::Anisotropic_mesh_3::Euclidean_metric_field<K>* MF;
typedef typename CGAL::Anisotropic_mesh_3::Custom_metric_field<K>* MF;

typedef std::set<std::size_t>                                Simplex;
typedef boost::array<std::size_t, 2>                         Edge;
typedef boost::array<std::size_t, 3>                         Tri;
typedef boost::array<std::size_t, 4>                         Tet;

typedef CGAL::Exact_predicates_exact_constructions_kernel    KExact;
typedef typename KExact::Point_3                             EPoint;
typedef typename KExact::Segment_3                           ESegment;
typedef typename KExact::Line_3                              ELine;
typedef typename KExact::Triangle_3                          ETriangle;
typedef CGAL::Cartesian_converter<K, KExact>                 To_exact;
typedef CGAL::Cartesian_converter<KExact, K>                 Back_from_exact;

To_exact to_exact;
Back_from_exact back_from_exact;

#define FILTER_SEEDS_OUTSIDE_GRID
#define USE_RECURSIVE_UPDATES

#define TMP_REFINEMENT_UGLY_HACK

#define verbose 6
const FT FT_inf = std::numeric_limits<FT>::infinity();

// using a lot of global variables, it's ugly but passing them by function is tedious
std::size_t vertices_nv = 10;

// grid related stuff
Point_3 center(0.,0.,0.);
const FT grid_side = 4.0;
FT offset_x = center.x() - grid_side/2.; // offset is the bottom left point
FT offset_y = center.y() - grid_side/2.;
FT offset_z = center.z() - grid_side/2.;
FT n = 50.; // number of points per side of the grid
FT sq_n = n*n;
FT step = grid_side / (n-1);

// the metric field and the seeds
MF mf;
std::vector<Point_3> seeds;
std::vector<Metric> seeds_m;

// how big of an update the new value needs to be
FT recursive_tolerance = 1e-5;

int n_refine = 20;
#ifdef TMP_REFINEMENT_UGLY_HACK
FT best_ref_x = FT_inf, best_ref_y = FT_inf, best_ref_z = FT_inf;
#endif

std::clock_t start;

enum FMM_state
{
  CHANGED = 0, // CHANGED implies KNOWN. It is used to not add the same element multiple times to the CHANGED priority queue
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

int insert_new_seed(const FT x, const FT y, const FT z)
{
#ifdef FILTER_SEEDS_OUTSIDE_GRID
    if(x < offset_x || x > center.x() + grid_side/2. ||
       y < offset_y || y > center.y() + grid_side/2. ||
       z < offset_z || z > center.z() + grid_side/2.)
    {
      std::cout << "filtered : " << x << " " << y << " " << z << std::endl;
      return seeds.size();
    }
#endif
  std::cout << "new seed: " << x << " " << y << " " << z << " added" << std::endl;

  seeds.push_back(Point_3(x, y, z));
  seeds_m.push_back(mf->compute_metric(seeds.back()));
  return seeds.size();
}

int build_seeds()
{
  std::ifstream in("bambimboum.mesh");
  std::string word;
  std::size_t useless, nv, dim;
  FT r_x, r_y, r_z;

  in >> word >> useless; //MeshVersionFormatted i
  in >> word >> dim; //Dimension d
  in >> word >> nv;
  std::cout << "file has : " << nv << " pts" << std::endl;
  assert(dim == 3);

  std::size_t min_nv = (std::min)(nv, vertices_nv);

  seeds.reserve(min_nv);
  seeds_m.reserve(min_nv);

  for(std::size_t i=0; i<nv; ++i)
  {
    in >> r_x >> r_y >> r_z >> useless;
//    insert_new_seed(offset_x, offset_y, offset_z);
//    insert_new_seed(1.6, offset_y, offset_z);
//    insert_new_seed(1.6, 1.6, 1.6);
//    insert_new_seed(offset_x, offset_y+grid_side/2., offset_z);
//    insert_new_seed(offset_x, offset_y, offset_z+grid_side/2.);
    insert_new_seed(r_x, r_y, r_z);

    if(seeds.size() >= vertices_nv)
      break;
  }
  std::cout << "seeds: " << seeds.size() << std::endl;
  return seeds.size();
}

struct Grid_point
{
  typedef boost::array<Grid_point*, 6> Neighbors;
  typedef boost::array<const Grid_point*, 3> Optimal_neighbors;

  Point_3 point;
  std::size_t index;
  mutable FT distance_to_closest_seed;
  mutable std::size_t closest_seed_id;
  mutable  FMM_state state;
  Neighbors neighbors; // 6 neighbors ORDERED (above, left, back, right, front, below)
                       // if no neighbor (borders for example) then NULL
  Metric metric;

  const Grid_point* ancestor; // used to debug & draw nice things

  // the neighbors of this for which we have reached the lowest distance_to_closest_seed
  // this is a 3-array with at least 1 non null pointer in [0]
  Optimal_neighbors opt_neighbors;

  bool compute_closest_seed_1D(const Grid_point* gp)
  {
    std::cout << "compute closest 1D " << index << " to " << gp->index;
    std::cout << " ds " << gp->distance_to_closest_seed;
    std::cout << " resp. " << gp->closest_seed_id << std::endl;

    CGAL_assertion(gp->state == KNOWN || gp->state == CHANGED);
    bool changed = false;

    Vector3d v;
    v(0) = gp->point.x() - point.x();
    v(1) = gp->point.y() - point.y();
    v(2) = gp->point.z() - point.z();

    const Eigen::Matrix3d& m = metric.get_mat();
    FT neighbor_d = std::sqrt(v.transpose() * m * v);
    FT dcs_at_gp = gp->distance_to_closest_seed;
    FT d = dcs_at_gp + neighbor_d;

    std::cout << "min 1D; min : " << d;
    std::cout << " vs current best: " << distance_to_closest_seed << std::endl;

    if(distance_to_closest_seed - d > recursive_tolerance * d)
    {
      changed = true;
      distance_to_closest_seed = d;

      opt_neighbors[0] = gp;
      opt_neighbors[1] = NULL; // just in case a _1D becomes better than a _2D or _3D
      opt_neighbors[2] = NULL;

      // below is useless fixme
      closest_seed_id = static_cast<std::size_t>(-1);
      ancestor = gp;

      std::cout << "1D new best for " << index << " : " << distance_to_closest_seed << std::endl;
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
    Vector3d vp, vq;
    vp(0) = gp->point.x() - point.x();
    vp(1) = gp->point.y() - point.y();
    vp(2) = gp->point.z() - point.z();
    vq(0) = gq->point.x() - point.x();
    vq(1) = gq->point.y() - point.y();
    vq(2) = gq->point.z() - point.z();

    const Eigen::Matrix3d& m = metric.get_mat();
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

    // Write f'(p) as Ap^2 + Bp + C
    const FT A = 4*(a-b)*(a-b)*c-4*c*c;
    const FT B = 4*(a-b)*(a-b)*d-4*c*d;
    const FT C = 4*(a-b)*(a-b)*e-d*d;
    FT delta = B*B-4*A*C;

    if(delta >= -1e17)
      delta = (std::max)(0., delta);

//    std::cout << "a: " << a << std::endl;
//    std::cout << "b: " << b << std::endl;
//    std::cout << "c: " << c << std::endl;
//    std::cout << "d: " << d << std::endl;
//    std::cout << "e: " << e << std::endl;
//    std::cout << "A: " << A << std::endl;
//    std::cout << "B: " << B << std::endl;
//    std::cout << "C: " << C << std::endl;
//    std::cout << "delta: " << delta << std::endl;

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

  bool compute_closest_seed_2D(const Grid_point* gp, const Grid_point* gq)
  {
    FT gp_d = gp->distance_to_closest_seed;
    FT gq_d = gq->distance_to_closest_seed;

    std::cout << "compute closest 2D " << index << " to ";
    std::cout << gp->index << " & " << gq->index;
    std::cout << " ds: " << gp_d << " " << gq_d << " ";
    std::cout << " resp. " << gp->closest_seed_id << " " << gq->closest_seed_id << std::endl;

    CGAL_assertion((gp->state == KNOWN || gp->state == CHANGED) &&
                   (gq->state == KNOWN || gq->state == CHANGED) );
    bool changed = false;

    FT p = 0.;
    const FT d = compute_min_distance_2D(gp, gq, p);

    std::cout << "min 2D; min : " << d << " at " << p;
    std::cout << " vs current best: " << distance_to_closest_seed << std::endl;

    // tolerance to ignore unsignificant changes
    if(distance_to_closest_seed - d > recursive_tolerance * d)
    {
      changed = true;
      distance_to_closest_seed = d;

      opt_neighbors[0] = gp;
      opt_neighbors[1] = gq;
      opt_neighbors[2] = NULL;

      // below is useless fixme
      closest_seed_id = static_cast<std::size_t>(-1);
      ancestor = gp;

      std::cout << "2D new best for " << index << " : " << distance_to_closest_seed << std::endl;
    }
    return changed;
  }

  bool check_causality_at_root(const Grid_point* gp, const Grid_point* gq,
                               const Grid_point* gr, FT x)
  {
    if(x < 0)
      return false;

    const Eigen::Matrix3d& m = metric.get_mat();
    FT gp_d = gp->distance_to_closest_seed;
    FT gq_d = gq->distance_to_closest_seed;
    FT gr_d = gr->distance_to_closest_seed;

    FT l = step; // all edge lengths are the same
    FT lden = 1./l;

    // gradient:
    const FT gx = lden * (x - gp_d);
    const FT gy = lden * (x - gq_d);
    const FT gz = lden * (x - gr_d);

    // caracteristic direction
    const FT dx = m(0,0)*gx + m(0,1)*gy + m(0,2)*gz;
    const FT dy = m(1,0)*gx + m(1,1)*gy + m(1,2)*gz;
    const FT dz = m(2,0)*gx + m(2,1)*gy + m(2,2)*gz;
    const Vector v(dx, dy, dz);
    const Line char_line(this->point, v);

#if (verbose > 20)
    std::cout << "gradient: " << gx << " " << gy << " " << gz << std::endl;
    std::cout << "and char direction: " << dx << " " << dy << " " << dz << std::endl;
#endif

    // intersect with the triangle PQR...
    const Triangle triangle(gp->point, gq->point, gr->point);

    const ETriangle etriangle = to_exact(triangle);
    const ELine exact_char_line = to_exact(char_line);

    CGAL::cpp11::result_of<KExact::Intersect_3(ELine, ETriangle)>::type
        result = CGAL::intersection(exact_char_line, etriangle);

    if (result)
    {
      if (const EPoint* p = boost::get<EPoint>(&*result))
      {
        // compute the barycentric coordinates to check if it's within the triangle ?
        std::cout << "pt intersection: " << p->x() << " " << p->y() << " " << p->z() << std::endl;
        return true; //fixme...
      }
      else if(boost::get<ESegment>(&*result))
      {
        // the intersection cannot be a segment...
        CGAL_assertion(false);
      }
    }
    else // no intersection...
    {
      std::cout << "no intersection" << std::endl;
      return false;
    }
  }

  bool compute_min_distance_3D(const Grid_point* gp, const Grid_point* gq,
                               const Grid_point* gr, FT& d)
  {
    // See theory for what's happening in there...
    FT gp_d = gp->distance_to_closest_seed;
    FT gq_d = gq->distance_to_closest_seed;
    FT gr_d = gr->distance_to_closest_seed;

    FT l = step; // all edge lengths are the same at the moment
    FT lden = 1./l;

    Vector3d vp, vq, vr; // note that these "points" towards 'this'
    vp(0) = point.x() - gp->point.x();
    vp(1) = point.y() - gp->point.y();
    vp(2) = point.z() - gp->point.z();
    vq(0) = point.x() - gq->point.x();
    vq(1) = point.y() - gq->point.y();
    vq(2) = point.z() - gq->point.z();
    vr(0) = point.x() - gr->point.x();
    vr(1) = point.y() - gr->point.y();
    vr(2) = point.z() - gr->point.z();

    Eigen::Matrix3d p_matrix;
    p_matrix(0,0) = lden*vp(0); p_matrix(0,1) = lden*vp(1); p_matrix(0,2) = lden*vp(2);
    p_matrix(1,0) = lden*vq(0); p_matrix(1,1) = lden*vq(1); p_matrix(1,2) = lden*vq(2);
    p_matrix(2,0) = lden*vr(0); p_matrix(2,1) = lden*vr(1); p_matrix(2,2) = lden*vr(2);

    // P^-1 = P^t since the edges 'this-gp', 'this-gq', 'this-gr' are orthogonal to each other
    Eigen::Matrix3d p_matrix_m1 = p_matrix.transpose();

    const Eigen::Matrix3d& m = metric.get_mat();

    FT g1 = lden * ( p_matrix_m1(0,0) + p_matrix_m1(0,1) + p_matrix_m1(0,2));
    FT g2 = - lden * ( p_matrix_m1(0,0)*gp_d + p_matrix_m1(0,1)*gq_d + p_matrix_m1(0,2)*gr_d);

    FT g3 = lden * ( p_matrix_m1(1,0) + p_matrix_m1(1,1) + p_matrix_m1(1,2));
    FT g4 = - lden * ( p_matrix_m1(1,0)*gp_d + p_matrix_m1(1,1)*gq_d + p_matrix_m1(1,2)*gr_d);

    FT g5 = lden * ( p_matrix_m1(2,0) + p_matrix_m1(2,1) + p_matrix_m1(2,2));
    FT g6 = - lden * ( p_matrix_m1(2,0)*gp_d + p_matrix_m1(2,1)*gq_d + p_matrix_m1(2,2)*gr_d);

    FT w1 = m(0,0)*g1*g1 + m(1,1)*g3*g3 + m(2,2)*g5*g5 + 2*m(0,1)*g1*g3 + 2*m(2,0)*g1*g5 + 2*m(1,2)*g3*g5;
    FT w2 = 2*(m(0,0)*g1*g2 + m(1,1)*g3*g4 + m(2,2)*g5*g6
               + m(0,1)*(g1*g4+g2*g3) + m(2,0)*(g1*g6+g2*g5) + m(1,2)*(g3*g6+g4*g5));
    FT w3 = m(0,0)*g2*g2 + m(1,1)*g4*g4 + m(2,2)*g6*g6
            + 2*(m(0,1)*g2*g4 + m(2,0)*g2*g6 + m(1,2)*g4*g5);

    FT delta = w2*w2 - 4*w1*(w3-1.);

    if(delta >= -1e17)
      delta = (std::max)(0., delta);

    if(delta < 0.)
    {
      // no real positive roots...
      return false;
    }

    CGAL_assertion(delta >= 0.);

    FT sqrt_delta = std::sqrt(delta);
    FT w1_den = 1./(2.*w1);
    FT x1 = (-w2 - sqrt_delta) * w1_den;
    FT x2 = (-w2 + sqrt_delta) * w1_den;

#if (verbose > 20)
// debug stuff
    std::cout << "debug of compute_min_distance_3D" << std::endl;
    std::cout << "gpgqgr " << gp->index << " " << gq->index << " " << gr->index << std::endl;
    std::cout << "resp dist " << gp_d << " " << gq_d << " " << gr_d << std::endl;
    std::cout << "p_matrix: " << std::endl << p_matrix << std::endl;
    std::cout << "inversed p_matrix: " << std::endl << p_matrix_m1 << std::endl;
    std::cout << "g12: " << g1 << " " << g2 << std::endl;
    std::cout << "g34: " << g3 << " " << g4 << std::endl;
    std::cout << "g56: " << g5 << " " << g6 << std::endl;
    std::cout << "ws: " << w1 << " " << w2 << " " << w3 << std::endl;
    std::cout << "delta: " << delta << std::endl;
    std::cout << "roots: " << x1 << " " << x2 << std::endl;
#endif

    // check both roots for positivity and causality
    bool is_x1_acceptable = check_causality_at_root(gp, gq, gr, x1);
    bool is_x2_acceptable = check_causality_at_root(gp, gq, gr, x2);

    std::cout << "acceptable: " << is_x1_acceptable << " " << is_x2_acceptable << std::endl;

    d = FT_inf;
    if(is_x1_acceptable)
      d = x1;
    if(is_x2_acceptable)
      d = (std::min)(d, x2);

    return (is_x1_acceptable || is_x2_acceptable);
  }

  FT compute_lambda_intersection(const Grid_point* gp)
  {
    FT lambda = -1.;
    if(!gp || gp->closest_seed_id == static_cast<std::size_t>(-1))
      return lambda;

    FT gp_d = gp->distance_to_closest_seed;
    const Eigen::Matrix3d& mp = gp->metric.get_mat();
    const Eigen::Matrix3d& m = metric.get_mat();

    const FT me1 = std::sqrt(m(0,0)), me2 = std::sqrt(m(1,1)), me3 = std::sqrt(m(2,2));
    const FT pe1 = std::sqrt(mp(0,0)), pe2 = std::sqrt(mp(1,1)), pe3 = std::sqrt(mp(2,2));

    if(gp->closest_seed_id != static_cast<std::size_t>(-1))
    {
      if(gp->index == index + 1 || gp->index == index - 1) // p pight/left of 'this'
        lambda = (gp_d - distance_to_closest_seed + pe1*step)/(me1 + pe1);
      else if(gp->index == index + n || gp->index == index - n) // p back/fpont of 'this'
        lambda = (gp_d - distance_to_closest_seed + pe2*step)/(me2 + pe2);
      else if(gp->index == index + sq_n || gp->index == index - sq_n) // p above/below 'this'
        lambda = (gp_d - distance_to_closest_seed + pe3*step)/(me3 + pe3);
      else
        CGAL_assertion(false);
    }

    // debug & vepifications
    Vector3d p_to_mp, this_to_mp;
    FT mpx = point.x() + lambda * (gp->point.x() - point.x()) / step;
    FT mpy = point.y() + lambda * (gp->point.y() - point.y()) / step;
    FT mpz = point.z() + lambda * (gp->point.z() - point.z()) / step;
    p_to_mp(0) = gp->point.x() - mpx;
    p_to_mp(1) = gp->point.y() - mpy;
    p_to_mp(2) = gp->point.z() - mpz;
    this_to_mp(0) = point.x() - mpx;
    this_to_mp(1) = point.y() - mpy;
    this_to_mp(2) = point.z() - mpz;
    FT d_p_to_mp = std::sqrt(p_to_mp.transpose() * mp * p_to_mp);
    FT d_this_to_mp = std::sqrt(this_to_mp.transpose() * m * this_to_mp);
    FT diffp = distance_to_closest_seed + d_this_to_mp - (gp_d + d_p_to_mp);

    if(lambda >= -1e-10 && lambda <= step+1e-10)
      CGAL_assertion(std::abs(diffp) < 1e-10);

    return lambda;
  }

  void determine_ancestor_2D(const Grid_point* gp, const Grid_point* gq,
                             const FT lambda_p, const FT lambda_q)
  {
    CGAL_assertion(gp && gp->closest_seed_id != static_cast<std::size_t>(-1) &&
                   gq && gq->closest_seed_id != static_cast<std::size_t>(-1));
    std::cout << "determine ancestor 2D " << gp->index << " " << gq->index;
    std::cout << " colors: " << gp->closest_seed_id << " " << gq->closest_seed_id << std::endl;

    if(gp->closest_seed_id == gq->closest_seed_id)
    {
      closest_seed_id = gp->closest_seed_id; // this probably creates problems in the ancestor tree...
      ancestor = gp;
    }

    // we have to use the lambdas to decide which one to use as ancestor!
    bool intersect_p = lambda_p >= -1e-10 && lambda_p <= step+1e-10;
    bool intersect_q = lambda_q >= -1e-10 && lambda_q <= step+1e-10;

    const Grid_point* best_g;
    if(!intersect_p)
    {
      if(!intersect_q) // !p && !q
        CGAL_assertion(false && "we intersect nothing... ?");
      else
        best_g = gp; // we intersect [this q] therefore 'this' is on the side of p
    }
    else
    {
      if(!intersect_q)
        best_g = gq; // we intersect [this p] therefore 'this' is on the side of p
      else
      {
        std::cout << "mega warning: intersect [this p] [this q]" << std::endl;
        best_g = (lambda_q < lambda_p)?gq:gp; // this is shady... fixme...?
      }
    }

    closest_seed_id = best_g->closest_seed_id;
    ancestor = best_g;

    std::cout << "ancestor shen 2D @ " << index << " color is : " << closest_seed_id << std::endl;
  }

  void determine_ancestor()
  {
    std::cout << "determine ancestor for " << index << " d: " << distance_to_closest_seed << std::endl;

    if(this->closest_seed_id != static_cast<std::size_t>(-1))
      return;  // might be better to use a boolean like 'is_the_grid_point_closest_to_a_seed' so that
               // we can still spread a color on top of another while refining

    CGAL_assertion(opt_neighbors[0]);
    const Grid_point* gp = opt_neighbors[0];
    const Grid_point* gq = opt_neighbors[1];
    const Grid_point* gr = opt_neighbors[2];

    std::cout << " " << (gq?(gr?"3":"2"):"1") << "D" << std::endl;

    // gp         non NULL => best is achieved from 1D
    // gp, gq     non NULL => best is achieved from 2D
    // gp, gq, gr non NULL => best is achieved from 3D

    bool is_init_p = gp->closest_seed_id != static_cast<std::size_t>(-1);
    bool is_init_q = gq?(gq->closest_seed_id != static_cast<std::size_t>(-1)):false;
    bool is_init_r = gr?(gr->closest_seed_id != static_cast<std::size_t>(-1)):false;

    std::cout << "ini bools: " << is_init_p << " " << is_init_q << " " << is_init_r;
    std::cout << " with colors: " << gp->closest_seed_id;
    if(gq)
      std::cout << " " << gq->closest_seed_id << " ";
    else
      std::cout << " -1 ";

    if(gr)
      std::cout << gr->closest_seed_id << std::endl;
    else
      std::cout << "-1" << std::endl;

    CGAL_assertion((is_init_p || is_init_q || is_init_r) &&
                   "gp & gq & gr have uninitialized seeds...");

    // quick filters
    if(is_init_p && !is_init_q && !is_init_r)
    {
      closest_seed_id = gp->closest_seed_id;
      ancestor = gp;
      return;
    }
    else if(!is_init_p && is_init_q && !is_init_r)
    {
      closest_seed_id = gq->closest_seed_id;
      ancestor = gq;
      return;
    }
    else if(!is_init_p && !is_init_q && is_init_r)
    {
      closest_seed_id = gr->closest_seed_id;
      ancestor = gr;
      return;
    }
    else if(is_init_p && is_init_q && is_init_r &&
            gp->closest_seed_id == gq->closest_seed_id &&
            gp->closest_seed_id == gr->closest_seed_id)
    {
      closest_seed_id = gp->closest_seed_id;
      ancestor = gp;
      return;
    }

    // we compute the intersection of the dual with [this,p], [this,q] and [this,q]
    // (we know that the dual intersects the face p-q-r).
    // depending on which side is intersected, we can deduce the color of 'this'

    // we're solving d(this) + d_this(this-mp) = d(p) + d_p(p-mp) on three edges
    FT lambda_p = compute_lambda_intersection(gp);
    FT lambda_q = compute_lambda_intersection(gq);
    FT lambda_r = compute_lambda_intersection(gr);

    // we've dealt with the case 'only 1 neighbor initialized', let's do 2 now...
    if(is_init_p && is_init_q && !is_init_r)
      return determine_ancestor_2D(gp, gq, lambda_p, lambda_q);
    else if(is_init_p && !is_init_q && is_init_r)
      return determine_ancestor_2D(gp, gr, lambda_p, lambda_r);
    else if(!is_init_p && is_init_q && is_init_r)
      return determine_ancestor_2D(gq, gr, lambda_q, lambda_r);

    CGAL_assertion(gp && gq && gr);

    bool intersect_p = lambda_p >= -1e-10 && lambda_p <= step+1e-10;
    bool intersect_q = lambda_q >= -1e-10 && lambda_q <= step+1e-10;
    bool intersect_r = lambda_r >= -1e-10 && lambda_r <= step+1e-10;

    std::cout << "lambdas: " << lambda_p << " " << lambda_q << " " << lambda_r << std::endl;
    std::cout << "intersect bools " << intersect_p << " " << intersect_q << " " << intersect_r << std::endl;

    const Grid_point* best_g = NULL; // this will be the ancestor of 'this'
    if(!intersect_p)
    {
      if(!intersect_q) // we do not intersect either[this p] nor [this q]... what to do ?
      {
        if(!intersect_r) // !p && !q && !r (we intersect nothing... !?)
        {
          CGAL_assertion(false && "we intersect nothing... ?");
        }
        else // !p && !q && r
        {
          // gotta choose between p & q...
          if(gp->closest_seed_id != gq->closest_seed_id)
            std::cout << "mega warning [this p] and [this q]" << std::endl;
          best_g = (lambda_q < lambda_p)?gq:gp; // this is shady... fixme...?
        }
      }
      else // q
      {
        if(!intersect_r) // !p && q && !r
        {
          // gotta choose between p & r...
          if(gp->closest_seed_id != gr->closest_seed_id)
            std::cout << "mega warning [this p] and [this r]" << std::endl;
          best_g = (lambda_r < lambda_p)?gr:gp; // this is shady... fixme...?
        }
        else // !p && q && r
        {
          best_g = gp; // we intersect [this,q] [this,r] therefore 'this' is on the side of p
        }
      }
    }
    else // p
    {
      if(!intersect_q)
      {
        if(!intersect_r) // p && !q && !r
        {
          // gotta choose between q & r...
          if(gq->closest_seed_id != gr->closest_seed_id)
            std::cout << "mega warning [this q] and [this r]" << std::endl;
          best_g = (lambda_r < lambda_q)?gr:gq; // this is shady... fixme...?
        }
        else // p && !q && r
        {
          best_g = gq;
        }
      }
      else // q
      {
        if(!intersect_r) // p && q && !r
        {
          best_g = gr;
        }
        else // p && q && r
        {
          // gotta choose between p & q & r...
          if(gp->closest_seed_id != gq->closest_seed_id ||
             (gr && gp->closest_seed_id != gr->closest_seed_id) ||
             (gr && gq->closest_seed_id != gr->closest_seed_id))
            std::cout << "ultra mega warning [this p], [this q] and [this r]" << std::endl;
          best_g = (lambda_q < lambda_p)?((lambda_r < lambda_q)?gr:gq):((lambda_r < lambda_p)?gr:gp);
        }
      }
    }

    CGAL_assertion(best_g);
    closest_seed_id = best_g->closest_seed_id;
    ancestor = best_g;
  }

  FT compute_closest_seed_3D(const Grid_point* gp,
                             const Grid_point* gq,
                             const Grid_point* gr)
  {
    FT gp_d = gp->distance_to_closest_seed;
    FT gq_d = gq->distance_to_closest_seed;
    FT gr_d = gr->distance_to_closest_seed;

    std::cout << "compute closest 3D " << index;
    std::cout << " to " << gp->index << ", " << gq->index << " & " << gr->index << " ";
    std::cout << " ds: " << gp_d << " " << gq_d << " " << gr_d;
    std::cout << " resp. " << gp->closest_seed_id << " ";
    std::cout << gq->closest_seed_id << " " << gr->closest_seed_id << std::endl;

    CGAL_assertion((gp->state == KNOWN || gp->state == CHANGED) &&
                   (gq->state == KNOWN || gq->state == CHANGED) &&
                   (gr->state == KNOWN || gr->state == CHANGED));
    bool changed = false;

    FT d = 0.;
#if 0
    bool is_3D_optimal = compute_min_distance_3D(gp, gq, gr, d);

    // if false, the solution is actually on one of the faces of the tetrahedron
    // so we call the 2D...
    if(!is_3D_optimal)
    {
      FT lambda;
      FT d_f1 = compute_min_distance_2D(gp, gq, lambda);
      FT d_f2 = compute_min_distance_2D(gp, gr, lambda);
      FT d_f3 = compute_min_distance_2D(gq, gr, lambda);

      d = (std::min)((std::min)(d_f1, d_f2), d_f3);
      std::cout << "3D wasn't optimal so 2D was called and found: " << d << std::endl;
    }
#endif

#ifdef BRUTE_FORCE_CHECK_OPTIMAL_P
    // brute force check
    Vector3d v1, v2;
    v1(0) = gq->point.x() - gp->point.x();
    v1(1) = gq->point.y() - gp->point.y();
    v1(2) = gq->point.z() - gp->point.z();
    v2(0) = gr->point.x() - gp->point.x();
    v2(1) = gr->point.y() - gp->point.y();
    v2(2) = gr->point.z() - gp->point.z();

    // loops originally are :
    //for(std::size_t ni=0; ni<=pt_n; ni++)
    //  for(std::size_t nj=0; nj<=pt_n-ni; nj++)

    // To parallelize these two loops, we create a single loop and to make it
    // less tedious to count, we create a grid of the parallelogram and
    // then ignore the pts that are outside the triangle (when ni+nj > pt_n)

    std::size_t pt_n = 15, tot = pt_n*(pt_n+1);
    FT increment = 1./static_cast<FT>(pt_n);
    FT min_d = FT_inf;
    std::size_t min_id = 0;

    // actually using the parallel version very much slows the code atm.... false sharing ...?
    //#pragma omp parallel for (don't forget the pragma critical if you uncomment here)
    for(std::size_t nij = 0; nij <= tot; ++nij)
    {
      std::size_t ni = nij / (pt_n+1);
      std::size_t nj = nij % (pt_n+1);

      if(ni+nj > pt_n)
        continue;

      FT i = ni * increment;
      FT j = nj * increment;

      Point_3 p(gp->point.x() + i*v1(0) + j*v2(0),
                gp->point.y() + i*v1(1) + j*v2(1),
                gp->point.z() + i*v1(2) + j*v2(2));

      // compute the bary weights
      FT lambda_p, lambda_q, lambda_r;
      Vector3d v3;
      v3(0) = p.x() - gp->point.x();
      v3(1) = p.y() - gp->point.y();
      v3(2) = p.z() - gp->point.z();

      FT d11 = v1.dot(v1);
      FT d12 = v1.dot(v2);
      FT d22 = v2.dot(v2);
      FT d31 = v3.dot(v1);
      FT d32 = v3.dot(v2);
      FT den = 1. / (d11 * d22 - d12 * d12);
      lambda_q = (d22 * d31 - d12 * d32) * den;
      lambda_r = (d11 * d32 - d12 * d31) * den;
      lambda_p = 1. - lambda_q - lambda_r;

      // interpolated value at P
      FT dcs_at_p = lambda_p * gp_d + lambda_q * gq_d + lambda_r * gr_d;

      CGAL_assertion(lambda_p >= -1e-10 && lambda_p <= 1.+1e-10 &&
                     lambda_q >= -1e-10 && lambda_q <= 1.+1e-10 &&
                     lambda_r >= -1e-10 && lambda_r <= 1.+1e-10);

      Vector3d v;
      v(0) = p.x() - point.x();
      v(1) = p.y() - point.y();
      v(2) = p.z() - point.z();
      FT neighbor_d = std::sqrt(v.transpose()*metric.get_mat()*v);
      FT d = dcs_at_p + neighbor_d;

//#pragma omp critical
//{
//      std::cout << "Debug compute_closest_seed_3D" << std::endl;
//      std::cout << "the three points : " << std::endl;
//      std::cout << gp->point.x() << " " << gp->point.x() << " " << gp->point.z() << std::endl;
//      std::cout << gq->point.y() << " " << gq->point.y() << " " << gq->point.z() << std::endl;
//      std::cout << gr->point.z() << " " << gr->point.z() << " " << gr->point.z() << std::endl;
//      std::cout << "ij:" << i << " " << j << " and p : " << p.x() << " " << p.y() << " " << p.z() << std::endl;
//      std::cout << "barywhite: " << lambda_p << " " << lambda_q << " " << lambda_r << std::endl;
//      std::cout << "check: ";
//      std::cout << lambda_p*gp->point.x()+lambda_q*gq->point.x()+lambda_r*gr->point.x() << " ";
//      std::cout << lambda_p*gp->point.y()+lambda_q*gq->point.y()+lambda_r*gr->point.y() << " ";
//      std::cout << lambda_p*gp->point.z()+lambda_q*gq->point.z()+lambda_r*gr->point.z() << std::endl;
//}
//#pragma omp critical
//{
      if(d < min_d)
      {
        min_id = nij;
        min_d = dp;
      }
    }

#endif
    d = min_d;

    std::cout << "min 3D; min : " << d;
    std::cout << " vs current best: " << distance_to_closest_seed << std::endl;

    // tolerance to ignore unsignificant changes
    if(distance_to_closest_seed - d > recursive_tolerance * d)
    {
      changed = true;
      distance_to_closest_seed = d;
      opt_neighbors[0] = gp;
      opt_neighbors[1] = gq;
      opt_neighbors[2] = gr;

      // below is useless fixme
      closest_seed_id = static_cast<std::size_t>(-1);
      ancestor = gp;

      std::cout << "3D new best for " << index << " : " << distance_to_closest_seed << std::endl;
    }

    return changed;
  }

  boost::array<const Grid_point*, 4> get_adjacent_neighbors(const Grid_point* gp,
                                                            const Grid_point* gq)
  {
    // Find, in the neighbors of p, the points adjacent to q
    // for example, if q = right, the adjacent neighbors are back, above, front, below

    CGAL_assertion(gp && gq && gp != gq);
    boost::array<const Grid_point*, 4> adj_n;

    // an annoying case is the border of the grid: not all neighbors are != NULL
    // so if we want to find the adjacent neighbors, we have to know what is the
    // neighbor opposite of gq (AND THAT OPPOSITE NEIGHBOR VERY WELL COULD BE 'NULL'!!)
    // we therefore use simply the index of the neighbors: we know the opposite couples:
    // 0-5, 1-3, and 2-4

    const boost::array<Grid_point*, 6>& ns = gp->neighbors; // just to get a shorter name...

    // first of all, we need to find the index of gq in gp->neighbors
    std::size_t pos_gq = 6; // 6 so that it will fail an assert later if gq not found
    for(std::size_t i=0; i<ns.size(); ++i)
      if(ns[i] == gq)
        pos_gq = i;

    // hardcoding all the cases, because it's simpler to code (and to understand)
    // considering pairs change the order, which will change 'nij' (but it's not a big deal)
    if(pos_gq == 0 || pos_gq == 5) // we want 1 2 3 4
    {
      adj_n[0] = ns[1]; adj_n[1] = ns[2]; adj_n[2] = ns[3]; adj_n[3] = ns[4];
    }
    else if(pos_gq == 1 || pos_gq == 3) // we want 0 2 5 4
    {
      adj_n[0] = ns[0]; adj_n[1] = ns[2]; adj_n[2] = ns[5]; adj_n[3] = ns[4];
    }
    else if(pos_gq == 2 || pos_gq == 4) // we want 0 1 5 3
    {
      adj_n[0] = ns[0]; adj_n[1] = ns[1]; adj_n[2] = ns[5]; adj_n[3] = ns[3];
    }

    return adj_n;
  }

  bool compute_closest_seed(const Grid_point* anc)
  {
    // returns true if we improved the distance

    std::cout << "compute closest high level " << index;
    std::cout << " with ancestor: " << anc->index;
    bool changed = false;

    // get the neighbors of 'this' that are adjacent of ancestor
    // then 4 possible 3D faces: this-ancestor-2 from adj_n
    // if a tet fails, we have two possible 2D faces: this-ancestor-1 from tet
    // if no face is available, then compute 1D: this-ancestor
    boost::array<const Grid_point*, 4> adj_n = get_adjacent_neighbors(this, anc);

    std::cout << " neighbors: ";
    for(std::size_t i=0; i<4; ++i)
      if(adj_n[i])
        std::cout << adj_n[i]->index << " ";
      else
        std::cout << "NULL ";
    std::cout << std::endl;

    for(std::size_t i=0; i<adj_n.size(); ++i)
    {
      const Grid_point* gp = adj_n[i];
      const Grid_point* gq = adj_n[(i+1)%(adj_n.size())];

      if(gp && (gp->state == KNOWN || gp->state == CHANGED))
      {
        if(gq && (gq->state == KNOWN || gq->state == CHANGED))
        {
          // the tet 'this - ancestor - gp - gq'
          if(compute_closest_seed_3D(anc, gp, gq))
            changed = true;
        }
        else
        { // the face 'this - ancestor - gp'
          if(compute_closest_seed_2D(anc, gp))
              changed = true;
        }
      }
      else if(gq && (gq->state == KNOWN || gq->state == CHANGED))
      { // the face 'this - ancestor - gq'
        if(compute_closest_seed_2D(anc, gq))
          changed = true;
      }
      else
      {
        if(compute_closest_seed_1D(anc))
          changed = true;
      }
    }

    std::cout << "End high level. color and value: " << closest_seed_id << " ";
    std::cout << distance_to_closest_seed << std::endl;

    return changed;
  }

  PQ_state update_neighbors_distances(std::vector<Grid_point*>& changed_pq,
                                      std::vector<Grid_point*>& trial_pq) const
  {
    std::cout << "update neighbors of " << index << std::endl;
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
          std::cout << "new value at " << gp->index << " : " << gp->distance_to_closest_seed;
          std::cout << " with ancestor: " << gp->ancestor->index;
          std::cout << " (mem: " << mem << " ancestor: " << ancestor_mem << ") ";
          std::cout << " diff: " << gp->distance_to_closest_seed-mem << std::endl;
          if(gp->state == KNOWN)
          {
            // if it's already 'changed' we only need to reorder the 'changed' queue
            // but we don't need to insert anything
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
        // no need to insert in the trial queue since it already is in it
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
        gp->compute_closest_seed(this); // always returns true here so no need to test it
        gp->state = TRIAL;
        trial_pq.push_back(gp);
        std::push_heap(trial_pq.begin(), trial_pq.end(), Grid_point_comparer<Grid_point>());
      }
    }
    return pqs_ret;
  }

  bool is_contained(const Point_3& p,
                    const Grid_point* next) const
  {
    // above : n0
    //    n2
    //    |
    // n1-gp-n3
    //    |
    //    n4
    // below : n5

    Grid_point *n0 = neighbors[0], *n1 = neighbors[1], *n2 = neighbors[2],
               *n3 = neighbors[3], *n4 = neighbors[4], *n5 = neighbors[5];

    bool is_p_geq_n0_z = (n1) ? p.z() <= n0->point.z() : false;
    bool is_p_leq_n1_x = (n1) ? p.x() <= n1->point.x() : false;
    bool is_p_geq_n2_y = (n2) ? p.y() >= n2->point.y() : false;
    bool is_p_geq_n3_x = (n3) ? p.x() >= n3->point.x() : false;
    bool is_p_leq_n4_y = (n4) ? p.y() <= n4->point.y() : false;
    bool is_p_leq_n5_z = (n5) ? p.z() <= n5->point.z() : false;

    if(is_p_geq_n0_z) { next = n0; return false; }
    if(is_p_leq_n1_x) { next = n1; return false; }
    if(is_p_geq_n2_y) { next = n2; return false; }
    if(is_p_geq_n3_x) { next = n3; return false; }
    if(is_p_leq_n4_y) { next = n4; return false; }
    if(is_p_leq_n5_z) { next = n5; return false; }
    next = this; // we're inside !
    return true;
  }

  void initialize_from_point(const Point_3& p,
                             const std::size_t seed_id) const
  {
    Vector3d v;
    v(0) = p.x() - point.x();
    v(1) = p.y() - point.y();
    v(2) = p.z() - point.z();
    const Eigen::Matrix3d& m = metric.get_mat();
    FT d = std::sqrt(v.transpose() * m * v);
    closest_seed_id = seed_id;
    distance_to_closest_seed = d;
    state = TRIAL;
  }

  Grid_point(const Point_3& point_, const std::size_t index_)
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

  const Grid_point* locate_point(const Point_3& p) const
  {
    bool found = false;
    const Grid_point& next = points[0];
    while(!found)
      found = (&next)->is_contained(p, &next);
    return &next;
  }

  void locate_and_initialize(const Point_3& p,
                             const std::size_t seed_id)
  {
    // WARNING: only usable for a full grid since we assume the step size is constant
    int index_x = std::floor((p.x()-offset_x)/step);
    int index_y = std::floor((p.y()-offset_y)/step);
    int index_z = std::floor((p.z()-offset_z)/step);

    Grid_point* gp = &(points[index_z*sq_n + index_y*n + index_x]);

    if(gp->closest_seed_id != static_cast<std::size_t>(-1))
    {
      std::cerr << "WARNING: a new seed is overwriting the closest seed id";
      std::cerr << " of a grid point! previous index is: " << gp->closest_seed_id << std::endl;
    }

    if(gp->state == TRIAL) // it's already a trial point, we can't accept two seeds for one grid pt
      CGAL_assertion(false && "the grid is not dense enough for the input...");

    gp->initialize_from_point(p, seed_id);
    trial_points.push_back(gp);
    std::make_heap(trial_points.begin(), trial_points.end(), Grid_point_comparer<Grid_point>());

    std::cout << "step reminder: " << step << std::endl;
    std::cout << "looking for p: " << p.x() << " " << p.y() << " " << p.z() << std::endl;
    std::cout << "found gp: " << gp->index << " [" << gp->point.x() << ", " << gp->point.y() << " " << gp->point.z() << "]" << std::endl;
    std::cout << "index: " << index_x << " " << index_y << " " << index_z << std::endl;
    std::cout << "offset: " << offset_x << " " << offset_y << " " << offset_z << std::endl;
  }

  void initialize_geo_grid()
  {
    std::cout << "grid initialization" << std::endl;

    // create the grid points (either a full grid or a smart grid (todo) )
    for(unsigned int k=0; k<n; ++k)
    {
      for(unsigned int j=0; j<n; ++j)
      {
        for(unsigned int i=0; i<n; ++i)
        {
          Point_3 p(offset_x+i*step, offset_y+j*step, offset_z+k*step); // fill from bot left to top right
          Grid_point gp(p, i + j*n + k*sq_n);
          points.push_back(gp);
        }
      }
    }

    // assign the neighbors
    for(unsigned int k=0; k<n; ++k)
    {
      for(unsigned int j=0; j<n; ++j)
      {
        for(unsigned int i=0; i<n; ++i)
        {
          std::size_t curr_id = i + j*n + k*sq_n;
          if(k != n-1) // there is a neighbor above
            points[curr_id].neighbors[0] = &(points[i + j*n + (k+1)*sq_n]);
          if(i != 0) // there is a neighbor left
            points[curr_id].neighbors[1] = &(points[i-1 + j*n + k*sq_n]);
          if(j != n-1) // there is a neighbor back
            points[curr_id].neighbors[2] = &(points[i + (j+1)*n + k*sq_n]);
          if(i != n-1) // there is a neighbor right
            points[curr_id].neighbors[3] = &(points[i+1 + j*n + k*sq_n]);
          if(j != 0) // there is a neighbor front
            points[curr_id].neighbors[4] = &(points[i + (j-1)*n + k*sq_n]);
          if(k != 0) // there is a neighbor below
            points[curr_id].neighbors[5] = &(points[i + j*n + (k-1)*sq_n]);
        }
      }
    }
    std::cout << "neighbors assigned" << std::endl;

    // seeds belong to a cube, find it
    for(std::size_t i=0; i<seeds.size(); ++i)
    {
      const Point_3& p = seeds[i];
//      const Grid_point* gp = locate_point(p);
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
    std::cout << "Determine ancestors" << std::endl;
    CGAL_assertion(changed_points.empty() && trial_points.empty());

    // TODO do something smarter when spreading a new color in an already colored
    // grid... One can probably just start from the new seed WITHOUT RESETING EVERYTHING
    // and then spread and only consider stuff if one of the possible ancestors
    // has the new color something like that...

    Grid_point_vector known_points;
    for(std::size_t i=0; i<points.size(); ++i)
    {
      // reset the colors (but not the values!)

      // fixme, just leave the colors uninitialized in the first pass and
      // there's no need to reset them here...

      Grid_point* gp = &(points[i]);
      gp->closest_seed_id = -1;
      gp->ancestor = NULL;

      CGAL_assertion(gp->state == KNOWN);

      known_points.push_back(gp);
    }
    std::make_heap(known_points.begin(), known_points.end(), Grid_point_comparer<Grid_point>());

    // only re initialize the closest_seed_id for each point-seed
    for(std::size_t i=0; i<seeds.size(); ++i)
    {
      // if you use an 8 pts initilization for a seed, you need to modify here too...
      const Point_3& p = seeds[i];
      int index_x = std::floor((p.x()-offset_x)/step);
      int index_y = std::floor((p.y()-offset_y)/step);
      int index_z = std::floor((p.z()-offset_z)/step);
      Grid_point* gp = &(points[index_z*sq_n + index_y*n + index_x]);
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
      if(print_states())
        break;
      std::cout << " ------------------------------------------------------- " << std::endl;
      std::cout << "Queue sizes. Trial: " << trial_points.size() << " Changed: " << changed_points.size() << std::endl;

//      std::cout << "changed heap: " << std::endl;
//      for (std::vector<Grid_point*>::iterator it = changed_points.begin();
//           it != changed_points.end(); ++it)
//        std::cout << (*it)->index << " " << (*it)->distance_to_closest_seed << std::endl;
//      std::cout << std::endl;

//      std::cout << "trial heap: " << std::endl;
//      for (std::vector<Grid_point*>::iterator it = trial_points.begin();
//           it != trial_points.end(); ++it)
//        std::cout << (*it)->index << " " << (*it)->distance_to_closest_seed << std::endl;
//      std::cout << std::endl;

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

      std::cout << "picked n° " << gp->index << " (" << gp->point.x() << " " << gp->point.y() << " " << gp->point.z() << ")";
      std::cout << "at distance : " << gp->distance_to_closest_seed << " from " << gp->closest_seed_id << std::endl;

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

    // it would probably be best to consider all the possible 2D cases here...
    // but more generally, 3D being weaker than 2D is a problem to be fixed fixme

    // detail of the problem :
    // 2D can be solved analytically while 3D cannot. Therefore we might sometimes
    // lose the exact result when we 'unlock' the 3D if the min is achieved on a face

    for(std::size_t i=0; i<points.size(); ++i)
    {
      Grid_point* gp = &(points[i]);
      std::cout << "point " << i << " min distance is supposedly: ";
      std::cout << gp->distance_to_closest_seed << std::endl;
      boost::array<Grid_point*, 6>::const_iterator it = gp->neighbors.begin(),
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

  Point_3 compute_refinement_point()
  {
    // todo
    return Point_3();
  }

  void refresh_grid_after_new_seed_creation()
  {
    CGAL_assertion(changed_points.empty() && trial_points.empty());
    std::size_t seed_id = seeds.size() - 1;
    const Point_3& p = seeds.back();
    locate_and_initialize(p, seed_id);
    geo_grid_loop();
  }

  void refine_grid()
  {
#if (verbose > 0)
    std::cout << "Refine grid !" << std::endl;
#endif

#ifdef TMP_REFINEMENT_UGLY_HACK
    std::cerr << "insert: " << best_ref_x << " " << best_ref_y << " " << best_ref_z << std::endl;
    Point_3 p(best_ref_x, best_ref_y, best_ref_z);
    CGAL_assertion(best_ref_x != FT_inf && best_ref_y != FT_inf && best_ref_z != FT_inf);
#else
    Point_3 p = compute_refinement_point();
#endif
    vertices_nv = insert_new_seed(p.x(), p.y(), p.z());
    refresh_grid_after_new_seed_creation();

#ifdef TMP_REFINEMENT_UGLY_HACK
    best_ref_x = FT_inf; best_ref_y = FT_inf; best_ref_z = FT_inf;
#endif
  }

  void add_tet_to_simplices(std::set<Tet>& simplices, const Grid_point& gp,
                            const Grid_point* n1, const Grid_point* n2, const Grid_point* n3) const
  {
    Tet t;
    if(n1 && n2 && n3)
    {
      t[0] = gp.index;
      t[1] = n1->index;
      t[2] = n2->index;
      t[3] = n3->index;
      std::sort(t.begin(), t.end());
      simplices.insert(t);
    }
  }

  void add_tet_to_simplices(std::set<Tet>& simplices, const Grid_point& gp,
                            std::size_t i1, std::size_t i2, std::size_t i3) const
  {
    return add_tet_to_simplices(simplices, gp, gp.neighbors[i1], gp.neighbors[i2], gp.neighbors[i3]);
  }

  void output_grid(const std::string str_base) const
  {
    std::ofstream out((str_base + ".mesh").c_str());
    std::ofstream out_bb((str_base + ".bb").c_str());
    out << "MeshVersionFormatted 1" << std::endl;
    out << "Dimension 3" << std::endl;
    out << "Vertices" << std::endl;
    out << points.size() << std::endl;
    out_bb << "3 1 " << points.size() << " 2" << std::endl;

    std::set<Tet> simplices;
    for(std::size_t i=0; i<points.size(); ++i)
    {
      const Grid_point& gp = points[i];
      out << gp.point.x() << " " << gp.point.y() << " " << gp.point.z() << " " << i+1 << std::endl;

//      out_bb << gp.closest_seed_id << std::endl;
      out_bb << gp.distance_to_closest_seed << std::endl;

      // compute the tets...
      const boost::array<Grid_point*, 6>& ns = gp.neighbors; // just to get a shorter name...
      for(std::size_t i=0; i<ns.size(); ++i)
      {
#if 0
        // lazily : the eight corner tetrahedra (creates self intersections,
        // but it doesn't really matter except visually)
        add_tet_to_simplices(simplices, gp, 0, 1, 2);
        add_tet_to_simplices(simplices, gp, 0, 2, 3);
        add_tet_to_simplices(simplices, gp, 0, 3, 4);
        add_tet_to_simplices(simplices, gp, 0, 4, 1);
        add_tet_to_simplices(simplices, gp, 5, 1, 2);
        add_tet_to_simplices(simplices, gp, 5, 2, 3);
        add_tet_to_simplices(simplices, gp, 5, 3, 4);
        add_tet_to_simplices(simplices, gp, 5, 4, 1);
#else
        // a bit smarter : a cube is decomposed in 5 tets and we pick them
        // correctly so we have no overlap
        add_tet_to_simplices(simplices, gp, 0, 2, 3);
        add_tet_to_simplices(simplices, gp, 0, 1, 4);
        add_tet_to_simplices(simplices, gp, 3, 4, 5);
        add_tet_to_simplices(simplices, gp, 1, 2, 5);

        // the last one is a bit annoying since we have to grab neighbors of neighbors
        const Grid_point* n4 = gp.neighbors[4];
        const Grid_point* n3 = gp.neighbors[3];
        if(n4 && n3)
          add_tet_to_simplices(simplices, gp, n4->neighbors[3],
                               n4->neighbors[0], n3->neighbors[0]);
#endif
      }
    }

    out << "Tetrahedra" << std::endl;
    out << simplices.size() << std::endl;
    for(typename std::set<Tet>::iterator it = simplices.begin();
                                             it != simplices.end(); ++it)
    {
      const Tet& tet = *it;
      std::set<std::size_t> materials;
      for(std::size_t i=0; i<tet.size(); ++i)
      {
        out << tet[i] + 1 << " ";
        materials.insert(points[tet[i]].closest_seed_id);
      }
      std::size_t mat = (materials.size()==1)?(*(materials.begin())):(seeds.size());
      out << mat << std::endl;
    }

    std::set<Tri> triangles;
    for(typename std::set<Tet>::iterator it = simplices.begin();
                                         it != simplices.end(); ++it)
    {
      const Tet& tet = *it;
      Tri tr; // could use a permutation to make it look nicer, I suppose...
      tr[0] = tet[0]; tr[1] = tet[1]; tr[2] = tet[2];
      std::sort(tr.begin(), tr.end()); triangles.insert(tr);
      tr[0] = tet[0]; tr[1] = tet[1]; tr[2] = tet[3];
      std::sort(tr.begin(), tr.end()); triangles.insert(tr);
      tr[0] = tet[0]; tr[1] = tet[2]; tr[2] = tet[3];
      std::sort(tr.begin(), tr.end()); triangles.insert(tr);
      tr[0] = tet[1]; tr[1] = tet[2]; tr[2] = tet[3];
      std::sort(tr.begin(), tr.end()); triangles.insert(tr);
    }

    out << "Triangles" << std::endl;
    out << triangles.size() << std::endl;
    for(typename std::set<Tri>::iterator it = triangles.begin();
                                         it != triangles.end(); ++it)
    {
      const Tri& tr = *it;
      std::set<std::size_t> materials;
      for(std::size_t i=0; i<tr.size(); ++i)
      {
        out << tr[i] + 1 << " ";
        materials.insert(points[tr[i]].closest_seed_id);
      }
      std::size_t mat = (materials.size()==1)?(*(materials.begin())):(seeds.size());
      out << mat << std::endl;
    }
  }

  void add_triangles_edges_to_dual_edges(const Tri& tr,
                                         std::set<Edge>& dual_edges) const
  {
    Edge e;
    e[0] = tr[0]; e[1] = tr[1]; dual_edges.insert(e);
    e[0] = tr[0]; e[1] = tr[2]; dual_edges.insert(e);
    e[0] = tr[1]; e[1] = tr[2]; dual_edges.insert(e);
  }

  void add_tet_triangles_to_dual_triangles(const Tet& tet,
                                           std::set<Tri>& dual_triangles,
                                           std::set<Edge>& dual_edges) const
  {
    Edge e;
    e[0] = tet[0]; e[1] = tet[1]; dual_edges.insert(e);
    e[0] = tet[0]; e[1] = tet[2]; dual_edges.insert(e);
    e[0] = tet[0]; e[1] = tet[3]; dual_edges.insert(e);
    e[0] = tet[1]; e[1] = tet[2]; dual_edges.insert(e);
    e[0] = tet[1]; e[1] = tet[3]; dual_edges.insert(e);
    e[0] = tet[2]; e[1] = tet[3]; dual_edges.insert(e);

    Tri tr;
    tr[0] = tet[0]; tr[1] = tet[1]; tr[2] = tet[2]; dual_triangles.insert(tr);
    tr[0] = tet[0]; tr[1] = tet[1]; tr[2] = tet[3]; dual_triangles.insert(tr);
    tr[0] = tet[0]; tr[1] = tet[2]; tr[2] = tet[3]; dual_triangles.insert(tr);
    tr[0] = tet[1]; tr[1] = tet[2]; tr[2] = tet[3]; dual_triangles.insert(tr);
  }

  void compute_dual(std::set<Edge>& dual_edges,
                    std::set<Tri>& dual_triangles,
                    std::set<Tet>& dual_tets) const
  {
    std::cout << "dual computations" << std::endl;

#ifdef TMP_REFINEMENT_UGLY_HACK
    FT farthest_dual_distance = 0.;
    const Grid_point* ref;
    FT back_up_farthest = 0.;
    const Grid_point* backup; // separate from ref because it could be farther
#endif

//#pragma omp parallel for (slows things down at the moment... false sharing ...?)
// need to put pragma criticals if parallel is used...

    for(std::size_t i=0; i<points.size(); ++i)
    {
      const Grid_point& gp = points[i];
      Simplex dual_simplex;
      dual_simplex.insert(gp.closest_seed_id);

      boost::array<Grid_point*, 6>::const_iterator it = gp.neighbors.begin(),
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
      }
      else if(dual_simplex.size() == 3) // a triangle!
      {
        typename Simplex::const_iterator it = dual_simplex.begin();
        Tri tr; tr[0] = *it; tr[1] = (*++it); tr[2] = (*++it);
        dual_triangles.insert(tr);

        add_triangles_edges_to_dual_edges(tr, dual_edges);

        if(dual_tets.empty() && gp.distance_to_closest_seed > back_up_farthest)
        {
          // this is just in case we don't find any tet...
          back_up_farthest = gp.distance_to_closest_seed;
          backup = &gp;
        }
      }
      else if(dual_simplex.size() >= 4)
      {
        if(dual_simplex.size() == 4) // a tetrahedron!
        {
          typename Simplex::const_iterator it = dual_simplex.begin();
          Tet tet; tet[0] = *it; tet[1] = (*++it); tet[2] = (*++it); tet[3] = (*++it);
          dual_tets.insert(tet);

          add_tet_triangles_to_dual_triangles(tet, dual_triangles, dual_edges);
        }
        else
          std::cerr << "WARNING: COSPHERICITY..." << std::endl;

#ifdef TMP_REFINEMENT_UGLY_HACK
        //todo : this could be gp + its neighbors instead of gp alone
        if(gp.distance_to_closest_seed > farthest_dual_distance)
        {
          farthest_dual_distance = gp.distance_to_closest_seed;
          ref = &gp;
        }
#endif
      }
    }
#ifdef TMP_REFINEMENT_UGLY_HACK
    if(dual_tets.empty()) // if no tet-dual exists, take the farthest triangle-dual
    {
      if(dual_triangles.empty())
      {
        std::cerr << "Couldn't find a ref point, need more initial points" << std::endl;
        output_grid("failed_geo_grid");
        std::exit(EXIT_FAILURE);
      }

      std::cerr << "WARNING: no tet captured, using back-up" << std::endl;
      best_ref_x = backup->point.x();
      best_ref_y = backup->point.y();
      best_ref_z = backup->point.z();
    }
    else
    {
      best_ref_x = ref->point.x();
      best_ref_y = ref->point.y();
      best_ref_z = ref->point.z();
    }
#endif
  }

  bool is_triangle_intersected(const Tri& tri, const std::set<Edge>& edges) const
  {
    bool is_intersected = false;

    // need to switch to epeck for the correct intersections here...
    const std::size_t i0 = tri[0];
    const std::size_t i1 = tri[1];
    const std::size_t i2 = tri[2];
    const Point_3& p0 = seeds[i0];
    const Point_3& p1 = seeds[i1];
    const Point_3& p2 = seeds[i2];
    const Triangle triangle(p0, p1, p2);

#pragma omp parallel shared(is_intersected, p0, p1, p2) // hardly noticeable increase of speed...
{
    for(std::set<Edge>::const_iterator it = edges.begin(); it!=edges.end(); ++it)
    {
#pragma omp single nowait // hack to parallelize the std::set
{
#pragma omp flush (is_intersected)
      if(!is_intersected) // hack because we're not allowed to 'break' (or 'continue')
                          // inside a pragma omp for
      {
        const Edge& edge = *it;
        bool local_intersected = false;

        const Segment segment(seeds[edge[0]], seeds[edge[1]]);

        ESegment esegment = to_exact(segment);
        ETriangle etriangle = to_exact(triangle);
        CGAL::cpp11::result_of<KExact::Intersect_3(ESegment, ETriangle)>::type
            result = CGAL::intersection(esegment, etriangle);

        if (result)
        {
          if (const EPoint* p = boost::get<EPoint>(&*result))
          {
            const EPoint& ep = *p;
            local_intersected = (ep!=to_exact(p0) && ep!=to_exact(p1) && ep!=to_exact(p2));
          }
          else if(const ESegment* s = boost::get<ESegment>(&*result))
          {
            const EPoint& ep0 = to_exact(p0);
            const EPoint& ep1 = to_exact(p1);
            const EPoint& ep2 = to_exact(p2);
            local_intersected = ((s->source()!=ep0 && s->source()!=ep1 && s->source()!=ep2) ||
                                 (s->target()!=ep0 && s->target()!=ep1 && s->target()!=ep2));
          }
        }
#pragma omp critical
{
        if(local_intersected)
          is_intersected = true;
}
#pragma omp flush (is_intersected)
      } // end !is_intersected
} // end single nowait
    } // end for
} // end parallel
    return is_intersected;
  }

  void output_dual(const std::string str_base) const
  {
    std::set<Edge> dual_edges;
    std::set<Tri> dual_triangles;
    std::set<Tet> dual_tets;
    compute_dual(dual_edges, dual_triangles, dual_tets);

    std::cout << "captured: ";
    std::cout << dual_edges.size() << " edges, ";
    std::cout << dual_triangles.size() << " triangles, ";
    std::cout << dual_tets.size() << " tets" << std::endl;

    std::ofstream out((str_base + "_dual.mesh").c_str());
    std::ofstream outbb((str_base + "_dual.bb").c_str());
    out << "MeshVersionFormatted 1" << std::endl;
    out << "Dimension 3" << std::endl;
    out << "Vertices" << std::endl;
    out << seeds.size() << std::endl;
    for(std::size_t i=0; i<seeds.size(); ++i)
      out << seeds[i].x() << " " << seeds[i].y() << " " << seeds[i].z() << " " << i+1 << std::endl;

    outbb << "3 1 " << dual_tets.size() + dual_triangles.size() << " 1" << std::endl;

    // TETRAHEDRA
    out << "Tetrahedra" << std::endl;
    out << dual_tets.size() << std::endl;
    for(typename std::set<Tet>::iterator it = dual_tets.begin(); it != dual_tets.end(); ++it)
    {
      const Tet& tet = *it;
      for(std::size_t i=0; i<tet.size(); ++i)
        out << tet[i] + 1 << " ";
      out << "1" << std::endl;

      // very inefficient to recompute them every time but whatever
      const Metric& m0 = mf->compute_metric(seeds[tet[0]]);
      const Metric& m1 = mf->compute_metric(seeds[tet[1]]);
      const Metric& m2 = mf->compute_metric(seeds[tet[2]]);
      const Metric& m3 = mf->compute_metric(seeds[tet[3]]);

      FT gamma01 = m0.compute_distortion(m1);
      FT gamma02 = m0.compute_distortion(m2);
      FT gamma03 = m0.compute_distortion(m3);
      FT gamma12 = m1.compute_distortion(m2);
      FT gamma13 = m1.compute_distortion(m3);
      FT gamma23 = m2.compute_distortion(m3);

      FT max_dist = (std::max)((std::max)(gamma01, gamma02), gamma03);
      max_dist = (std::max)((std::max)(max_dist, gamma12), gamma13);
      max_dist = (std::max)(max_dist, gamma23);

      outbb << max_dist << std::endl;
    }

     // TRIANGLES
    out << "Triangles" << std::endl;
    out << dual_triangles.size() << std::endl;
    for(typename std::set<Tri>::iterator it = dual_triangles.begin();
                                         it != dual_triangles.end(); ++it)
    {
      const Tri& tr = *it;
      for(std::size_t i=0; i<tr.size(); ++i)
        out << tr[i] + 1 << " ";
      out << is_triangle_intersected(tr, dual_edges) << std::endl;

      // very inefficient to recompute them every time but whatever
      const Metric& m0 = mf->compute_metric(seeds[tr[0]]);
      const Metric& m1 = mf->compute_metric(seeds[tr[1]]);
      const Metric& m2 = mf->compute_metric(seeds[tr[2]]);

      FT gamma01 = m0.compute_distortion(m1);
      FT gamma02 = m0.compute_distortion(m2);
      FT gamma12 = m1.compute_distortion(m2);

      FT max_dist = (std::max)((std::max)(gamma01, gamma02), gamma12);
      outbb << max_dist << std::endl;
    }
    outbb << "End" << std::endl;

    // edges are not printed by medit, so we don't bother...

    out << "End" << std::endl;
  }

  Geo_grid() : points(), changed_points(), trial_points() { }
};

void initialize()
{
//  mf = new Euclidean_metric_field<K>(1., 1., 10.);
  mf = new Custom_metric_field<K>();

  vertices_nv = build_seeds();
  CGAL_assertion(vertices_nv > 0 && "No seed in domain..." );
}

int main(int, char**)
{
  std::cout.precision(17);
  std::freopen("geo_grid_log.txt", "w", stdout);

  start = std::clock();
  double duration;
  boost::chrono::system_clock::time_point t_start = boost::chrono::system_clock::now();

  std::srand(0);
  initialize();

  Geo_grid gg;
  gg.build_grid();
//  gg.output_grid("geo_grid_2_pre");
  gg.output_dual("geo_grid_2_pre");

  for(int i=0; i<n_refine; ++i)
  {
#ifndef TMP_REFINEMENT_UGLY_HACK
    std::cerr << "you forgot to enable the macro, birdhead" << std::endl;
#endif

    gg.refine_grid();
    std::ostringstream out;
    out << "geo_grid_2_ref_" << i;
//    gg.output_grid(out.str());
    gg.output_dual(out.str());
  }

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cerr << "End refinement: " << duration << std::endl;

  gg.output_grid("geo_grid_2");
  gg.output_dual("geo_grid_2");

  boost::chrono::system_clock::time_point t_end = boost::chrono::system_clock::now();
  boost::chrono::milliseconds t = boost::chrono::duration_cast<boost::chrono::milliseconds> (t_end-t_start);
//  std::cerr << "Wall: " << t.count() << std::endl;

  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cerr << "Total: " << duration << std::endl;
}

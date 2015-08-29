// This is based on the paper :
// Konukoglu et al., A recursive anisotropic fast marching approach to reaction diffusion equation

// The canvas used is an orthogonal grid

#include <CGAL/Starset.h>
#include <CGAL/Stretched_Delaunay_3.h>

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <CGAL/Voronoi_painter.h>

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

#define USE_RECURSIVE_UPDATES

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K, typename KExact, typename Canvas_point, typename Metric_field>
class Konukoglu_canvas;

template<typename K, typename KExact, typename Metric_field>
class Konukoglu_canvas_point :
    public Canvas_point<K>
{
private:
  typedef Konukoglu_canvas_point<K, KExact, Metric_field>     Self;

public:
  typedef boost::array<Konukoglu_canvas_point*, 6>            Neighbors;
  typedef Canvas_point<K>                                     Base;
  typedef Konukoglu_canvas<K, KExact, Self, Metric_field>     Canvas;


  typedef typename Base::FT                                   FT;
  typedef typename Base::Point_3                              Point_3;
  typedef typename Base::Vector3d                             Vector3d;

  typedef boost::array<const Konukoglu_canvas_point*, 3>      Optimal_neighbors;

// Neighbors
  Neighbors neighbors;

  // the neighbors of this for which we have reached the lowest distance_to_closest_seed
  // this is a 3-array with at least 1 non null pointer in [0]
  Optimal_neighbors opt_neighbors;

  // we need a lot of info from the grid with this algorithm...
  Canvas* canvas; // forced to take a pointer and not a ref since we copy points

  bool compute_closest_seed_1D(const Konukoglu_canvas_point* cp)
  {
#if (verbose > 15)
    std::cout << "compute closest 1D " << this->index << " to " << cp->index;
    std::cout << " ds " << cp->distance_to_closest_seed;
    std::cout << " resp. " << cp->closest_seed_id << std::endl;
#endif
    CGAL_assertion(cp->state == KNOWN || cp->state == CHANGED);
    bool changed = false;

    Vector3d v;
    v(0) = cp->point.x() - this->point.x();
    v(1) = cp->point.y() - this->point.y();
    v(2) = cp->point.z() - this->point.z();

    const Eigen::Matrix3d& m = this->metric.get_mat();
    FT neighbor_d = std::sqrt(v.transpose() * m * v);
    FT dcs_at_cp = cp->distance_to_closest_seed;
    FT d = dcs_at_cp + neighbor_d;

#if (verbose > 15)
    std::cout << "min 1D; min : " << d;
    std::cout << " vs current best: " << this->distance_to_closest_seed << std::endl;
#endif

    if(this->distance_to_closest_seed - d > canvas->recursive_tolerance * d)
    {
      changed = true;
      this->distance_to_closest_seed = d;

      opt_neighbors[0] = cp;
      opt_neighbors[1] = NULL; // just in case a _1D becomes better than a _2D or _3D
      opt_neighbors[2] = NULL;

      // below is useless fixme
      this->closest_seed_id = static_cast<std::size_t>(-1);
      this->ancestor = cp;

#if (verbose > 15)
      std::cout << "1D new best for " << this->index << " : "
                << this->distance_to_closest_seed << std::endl;
#endif
    }
    return changed;
  }

  FT eval_at_p(FT p, const FT a, const FT b,
               const FT c, const FT d, const FT e) const
  {
    if(std::abs(c*p*p + d*p + e) < 1e-16) // avoids stuff like -1e-17
      return p*a+(1.-p)*b;
    CGAL_assertion(c*p*p + d*p + e >= 0);
    return p*a+(1.-p)*b + std::sqrt(c*p*p + d*p + e);
  }

  FT compute_min_distance_2D(const Konukoglu_canvas_point* cp,
                             const Konukoglu_canvas_point* cq,
                             FT& p) const
  {
    CGAL_assertion((cp->state == KNOWN || cp->state == CHANGED) &&
                   (cq->state == KNOWN || cq->state == CHANGED) );

    // in 2D, we can solve the minimization analytically.
    // We're solving :
    // min_{p\in [0,1]} f(p) = F(x)*p + F(y)*(1-p) + sqrt(v(p)^t * M * v(p))
    // vp = p*cp+(1-p)*cq
    // local minimum is given by f'(p_0) = 0 _but_ we have to compare f(p_0)
    // with f(0) and f(1) as we're only interested in solutions that live in [0,1]
    // If there's no solution in [0,1], then there's no (local) minimum in [0,1]
    // and the minimum is obtained at 0 or at 1.

    const FT dcs_at_cp = cp->distance_to_closest_seed;
    const FT dcs_at_cq = cq->distance_to_closest_seed;
    Vector3d vp, vq;
    vp(0) = cp->point.x() - this->point.x();
    vp(1) = cp->point.y() - this->point.y();
    vp(2) = cp->point.z() - this->point.z();
    vq(0) = cq->point.x() - this->point.x();
    vq(1) = cq->point.y() - this->point.y();
    vq(2) = cq->point.z() - this->point.z();

    const Eigen::Matrix3d& m = this->metric.get_mat();
    const FT vpmvp = vp.transpose() * m * vp;
    const FT vpmvq = vp.transpose() * m * vq;
    const FT vqmvq = vq.transpose() * m * vq;

    // expressing f(p) as ap+(1-p)b+sqrt(cp^2+dp+e)
    // the tiny trick is to write f'(p) = 0 and then square the expression to get rid of sqrt

    const FT a = dcs_at_cp;
    const FT b = dcs_at_cq;
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

  bool compute_closest_seed_2D(const Konukoglu_canvas_point* cp,
                               const Konukoglu_canvas_point* cq)
  {
#if (verbose > 15)
    FT cp_d = cp->distance_to_closest_seed;
    FT cq_d = cq->distance_to_closest_seed;
    std::cout << "compute closest 2D " << this->index << " to ";
    std::cout << cp->index << " & " << cq->index;
    std::cout << " ds: " << cp_d << " " << cq_d << " ";
    std::cout << " resp. " << cp->closest_seed_id << " " << cq->closest_seed_id << std::endl;
#endif

    CGAL_assertion((cp->state == KNOWN || cp->state == CHANGED) &&
                   (cq->state == KNOWN || cq->state == CHANGED) );
    bool changed = false;

    FT p = 0.;
    const FT d = compute_min_distance_2D(cp, cq, p);

#if (verbose > 15)
    std::cout << "min 2D; min : " << d << " at " << p;
    std::cout << " vs current best: " << this->distance_to_closest_seed << std::endl;
#endif

    // tolerance to ignore unsignificant changes
    if(this->distance_to_closest_seed - d > canvas->recursive_tolerance * d)
    {
      changed = true;
      this->distance_to_closest_seed = d;

      opt_neighbors[0] = cp;
      opt_neighbors[1] = cq;
      opt_neighbors[2] = NULL;

      // below is useless fixme
      this->closest_seed_id = static_cast<std::size_t>(-1);
      this->ancestor = cp;

#if (verbose > 15)
      std::cout << "2D new best for " << this->index << " : "
                << this->distance_to_closest_seed << std::endl;
#endif
    }
    return changed;
  }

  bool check_causality_at_root(const Konukoglu_canvas_point* cp,
                               const Konukoglu_canvas_point* cq,
                               const Konukoglu_canvas_point* cr,
                               FT x)
  {
    if(x < 0)
      return false;

    typedef typename K::Vector_3                        Vector;
    typedef typename K::Line_3                          Line;
    typedef typename K::Triangle_3                      Triangle;

    typedef typename KExact::Point_3                    EPoint;
    typedef typename KExact::Segment_3                  ESegment;
    typedef typename KExact::Line_3                     ELine;
    typedef typename KExact::Triangle_3                 ETriangle;


    const Eigen::Matrix3d& m = this->metric.get_mat();
    FT cp_d = cp->distance_to_closest_seed;
    FT cq_d = cq->distance_to_closest_seed;
    FT cr_d = cr->distance_to_closest_seed;

    FT l = canvas->step; // all edge lengths are the same
    FT lden = 1./l;

    // gradient:
    const FT gx = lden * (x - cp_d);
    const FT gy = lden * (x - cq_d);
    const FT gz = lden * (x - cr_d);

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
    const Triangle triangle(cp->point, cq->point, cr->point);

    const ETriangle etriangle = to_exact(triangle);
    const ELine exact_char_line = to_exact(char_line);

    typename CGAL::cpp11::result_of<typename KExact::Intersect_3(ELine, ETriangle)>::type
                        result = CGAL::intersection(exact_char_line, etriangle);

    if(result)
    {
      if(const EPoint* p = boost::get<EPoint>(&*result))
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

  bool compute_min_distance_3D(const Konukoglu_canvas_point* cp,
                               const Konukoglu_canvas_point* cq,
                               const Konukoglu_canvas_point* cr,
                               FT& d)
  {
    // See theory for what's happening in there...
    FT cp_d = cp->distance_to_closest_seed;
    FT cq_d = cq->distance_to_closest_seed;
    FT cr_d = cr->distance_to_closest_seed;

    FT l = canvas->step; // all edge lengths are the same at the moment
    FT lden = 1./l;

    Vector3d vp, vq, vr; // note that these "points" towards 'this'
    vp(0) = this->point.x() - cp->point.x();
    vp(1) = this->point.y() - cp->point.y();
    vp(2) = this->point.z() - cp->point.z();
    vq(0) = this->point.x() - cq->point.x();
    vq(1) = this->point.y() - cq->point.y();
    vq(2) = this->point.z() - cq->point.z();
    vr(0) = this->point.x() - cr->point.x();
    vr(1) = this->point.y() - cr->point.y();
    vr(2) = this->point.z() - cr->point.z();

    Eigen::Matrix3d p_matrix;
    p_matrix(0,0) = lden*vp(0); p_matrix(0,1) = lden*vp(1); p_matrix(0,2) = lden*vp(2);
    p_matrix(1,0) = lden*vq(0); p_matrix(1,1) = lden*vq(1); p_matrix(1,2) = lden*vq(2);
    p_matrix(2,0) = lden*vr(0); p_matrix(2,1) = lden*vr(1); p_matrix(2,2) = lden*vr(2);

    // P^-1 = P^t since the edges 'this-cp', 'this-cq', 'this-cr' are orthogonal to each other
    Eigen::Matrix3d p_matrix_m1 = p_matrix.transpose();

    const Eigen::Matrix3d& m = this->metric.get_mat();

    FT g1 = lden * ( p_matrix_m1(0,0) + p_matrix_m1(0,1) + p_matrix_m1(0,2));
    FT g2 = - lden * ( p_matrix_m1(0,0)*cp_d + p_matrix_m1(0,1)*cq_d + p_matrix_m1(0,2)*cr_d);

    FT g3 = lden * ( p_matrix_m1(1,0) + p_matrix_m1(1,1) + p_matrix_m1(1,2));
    FT g4 = - lden * ( p_matrix_m1(1,0)*cp_d + p_matrix_m1(1,1)*cq_d + p_matrix_m1(1,2)*cr_d);

    FT g5 = lden * ( p_matrix_m1(2,0) + p_matrix_m1(2,1) + p_matrix_m1(2,2));
    FT g6 = - lden * ( p_matrix_m1(2,0)*cp_d + p_matrix_m1(2,1)*cq_d + p_matrix_m1(2,2)*cr_d);

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
    std::cout << "cpcqcr " << cp->index << " " << cq->index << " " << cr->index << std::endl;
    std::cout << "resp dist " << cp_d << " " << cq_d << " " << cr_d << std::endl;
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
    bool is_x1_acceptable = check_causality_at_root(cp, cq, cr, x1);
    bool is_x2_acceptable = check_causality_at_root(cp, cq, cr, x2);

    std::cout << "acceptable: " << is_x1_acceptable << " " << is_x2_acceptable << std::endl;

    d = FT_inf;
    if(is_x1_acceptable)
      d = x1;
    if(is_x2_acceptable)
      d = (std::min)(d, x2);

    return (is_x1_acceptable || is_x2_acceptable);
  }

  FT compute_lambda_intersection(const Konukoglu_canvas_point* cp)
  {
    // compute the intersection point P between this & cp such that
    // d_this + d(this,P) = d(P,cp) + d_p

    FT lambda = -1.;
    if(!cp || cp->closest_seed_id == static_cast<std::size_t>(-1))
      return lambda;

    FT cp_d = cp->distance_to_closest_seed;
    const Eigen::Matrix3d& mp = cp->metric.get_mat();
    const Eigen::Matrix3d& m = this->metric.get_mat();

    const FT me1 = std::sqrt(m(0,0)), me2 = std::sqrt(m(1,1)), me3 = std::sqrt(m(2,2));
    const FT pe1 = std::sqrt(mp(0,0)), pe2 = std::sqrt(mp(1,1)), pe3 = std::sqrt(mp(2,2));

    if(cp->closest_seed_id != static_cast<std::size_t>(-1))
    {
      if(cp->index == this->index + 1 || cp->index == this->index - 1) // p pight/left of 'this'
        lambda = (cp_d - this->distance_to_closest_seed + pe1 * canvas->step) / (me1 + pe1);
      else if(cp->index == this->index + canvas->n || cp->index == this->index - canvas->n) // p back/fpont of 'this'
        lambda = (cp_d - this->distance_to_closest_seed + pe2 * canvas->step) / (me2 + pe2);
      else if(cp->index == this->index + canvas->sq_n || cp->index == this->index - canvas->sq_n) // p above/below 'this'
        lambda = (cp_d - this->distance_to_closest_seed + pe3 * canvas->step) / (me3 + pe3);
      else
        CGAL_assertion(false);
    }

    // debug & vepifications
    Vector3d p_to_mp, this_to_mp;
    FT mpx = this->point.x() + lambda * (cp->point.x() - this->point.x()) / canvas->step;
    FT mpy = this->point.y() + lambda * (cp->point.y() - this->point.y()) / canvas->step;
    FT mpz = this->point.z() + lambda * (cp->point.z() - this->point.z()) / canvas->step;
    p_to_mp(0) = cp->point.x() - mpx;
    p_to_mp(1) = cp->point.y() - mpy;
    p_to_mp(2) = cp->point.z() - mpz;
    this_to_mp(0) = this->point.x() - mpx;
    this_to_mp(1) = this->point.y() - mpy;
    this_to_mp(2) = this->point.z() - mpz;
    FT d_p_to_mp = std::sqrt(p_to_mp.transpose() * mp * p_to_mp);
    FT d_this_to_mp = std::sqrt(this_to_mp.transpose() * m * this_to_mp);
    FT diffp = this->distance_to_closest_seed + d_this_to_mp - (cp_d + d_p_to_mp);

    if(lambda >= -1e-10 && lambda <= canvas->step + 1e-10)
      CGAL_assertion(std::abs(diffp) < 1e-10);

    return lambda;
  }

  void determine_ancestor_2D(const Konukoglu_canvas_point* cp,
                             const Konukoglu_canvas_point* cq,
                             const FT lambda_p, const FT lambda_q)
  {
    CGAL_assertion(cp && cp->closest_seed_id != static_cast<std::size_t>(-1) &&
                   cq && cq->closest_seed_id != static_cast<std::size_t>(-1));

#if (verbose > 15)
    std::cout << "determine ancestor 2D " << cp->index << " " << cq->index;
    std::cout << " colors: " << cp->closest_seed_id << " " << cq->closest_seed_id << std::endl;
#endif

    if(cp->closest_seed_id == cq->closest_seed_id)
    {
      this->closest_seed_id = cp->closest_seed_id; // this probably creates problems in the ancestor tree...
      this->ancestor = cp;
    }

    // we have to use the lambdas to decide which one to use as ancestor!
    bool intersect_p = lambda_p >= -1e-2 && lambda_p <= canvas->step + 1e-10;
    bool intersect_q = lambda_q >= -1e-2 && lambda_q <= canvas->step + 1e-10;
    //fixme -1e-2 is absurdly big for a tolerance...

    const Konukoglu_canvas_point* best_g;
    if(!intersect_p)
    {
      if(!intersect_q) // !p && !q
        CGAL_assertion(false && "we intersect nothing... ?");
      else
        best_g = cp; // we intersect [this q] therefore 'this' is on the side of p
    }
    else
    {
      if(!intersect_q)
        best_g = cq; // we intersect [this p] therefore 'this' is on the side of p
      else
      {
        std::cerr << "mega warning: intersect [this p] [this q]" << std::endl;
        best_g = (lambda_q < lambda_p)?cq:cp; // this is shady... fixme...?
      }
    }

    this->closest_seed_id = best_g->closest_seed_id;
    this->ancestor = best_g;

#if (verbose > 15)
    std::cout << "ancestor shen 2D @ " << this->index
              << " color is : " << this->closest_seed_id << std::endl;
#endif
  }

  void determine_ancestor()
  {
#if (verbose > 10)
    std::cout << "determine ancestor for " << this->index
              << " d: " << this->distance_to_closest_seed << std::endl;
#endif

    if(this->closest_seed_id != static_cast<std::size_t>(-1))
      return;
    // might be better to use a boolean like 'is_the_canvas_point_closest_to_a_seed' so that
    // we can still spread a color on top of another while refining

    CGAL_assertion(opt_neighbors[0]);
    const Konukoglu_canvas_point* cp = opt_neighbors[0];
    const Konukoglu_canvas_point* cq = opt_neighbors[1];
    const Konukoglu_canvas_point* cr = opt_neighbors[2];

#if (verbose > 10)
    std::cout << " " << (cq?(cr?"3":"2"):"1") << "D" << std::endl;
#endif

    // cp         non NULL => best is achieved from 1D
    // cp, cq     non NULL => best is achieved from 2D
    // cp, cq, cr non NULL => best is achieved from 3D

    bool is_init_p = cp->closest_seed_id != static_cast<std::size_t>(-1);
    bool is_init_q = cq?(cq->closest_seed_id != static_cast<std::size_t>(-1)):false;
    bool is_init_r = cr?(cr->closest_seed_id != static_cast<std::size_t>(-1)):false;

#if (verbose > 20)
    std::cout << "ini bools: " << is_init_p << " " << is_init_q << " " << is_init_r;
    std::cout << " with colors: " << cp->closest_seed_id;

    if(cq)
      std::cout << " " << cq->closest_seed_id << " ";
    else
      std::cout << " -1 ";

    if(cr)
      std::cout << cr->closest_seed_id << std::endl;
    else
      std::cout << "-1" << std::endl;
#endif

    CGAL_assertion((is_init_p || is_init_q || is_init_r) &&
                   "cp & cq & cr have uninitialized seeds...");

    // quick filters
    if(is_init_p && !is_init_q && !is_init_r)
    {
      this->closest_seed_id = cp->closest_seed_id;
      this->ancestor = cp;
      return;
    }
    else if(!is_init_p && is_init_q && !is_init_r)
    {
      this->closest_seed_id = cq->closest_seed_id;
      this->ancestor = cq;
      return;
    }
    else if(!is_init_p && !is_init_q && is_init_r)
    {
      this->closest_seed_id = cr->closest_seed_id;
      this->ancestor = cr;
      return;
    }
    else if(is_init_p && is_init_q && is_init_r &&
            cp->closest_seed_id == cq->closest_seed_id &&
            cp->closest_seed_id == cr->closest_seed_id)
    {
      this->closest_seed_id = cp->closest_seed_id;
      this->ancestor = cp;
      return;
    }

    // we compute the intersection of the dual with [this,p], [this,q] and [this,q]
    // (we know that the dual intersects the face p-q-r).
    // depending on which side is intersected, we can deduce the color of 'this'

    // we're solving d(this) + d_this(this-mp) = d(p) + d_p(p-mp) on three edges
    FT lambda_p = compute_lambda_intersection(cp);
    FT lambda_q = compute_lambda_intersection(cq);
    FT lambda_r = compute_lambda_intersection(cr);

    // we've dealt with the case 'only 1 neighbor initialized', let's do 2 now...
    if(is_init_p && is_init_q && !is_init_r)
      return determine_ancestor_2D(cp, cq, lambda_p, lambda_q);
    else if(is_init_p && !is_init_q && is_init_r)
      return determine_ancestor_2D(cp, cr, lambda_p, lambda_r);
    else if(!is_init_p && is_init_q && is_init_r)
      return determine_ancestor_2D(cq, cr, lambda_q, lambda_r);

    CGAL_assertion(cp && cq && cr);

    bool intersect_p = lambda_p >= -1e-2 && lambda_p <= canvas->step + 1e-10;
    bool intersect_q = lambda_q >= -1e-2 && lambda_q <= canvas->step + 1e-10;
    bool intersect_r = lambda_r >= -1e-2 && lambda_r <= canvas->step + 1e-10;
    //fixme tolerance too big...

#if (verbose > 20)
    std::cout << "lambdas: " << lambda_p << " " << lambda_q << " " << lambda_r << std::endl;
    std::cout << "intersect bools " << intersect_p << " " << intersect_q << " " << intersect_r << std::endl;
#endif

    // this will be the ancestor of 'this'
    const Konukoglu_canvas_point* best_g = NULL;
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
#if (verbose > 10)
          if(cp->closest_seed_id != cq->closest_seed_id)
            std::cerr << "mega warning [this p] and [this q]" << std::endl;
#endif
          best_g = (lambda_q < lambda_p)?cq:cp; // this is shady... fixme...?
        }
      }
      else // q
      {
        if(!intersect_r) // !p && q && !r
        {
          // gotta choose between p & r...
#if (verbose > 10)
          if(cp->closest_seed_id != cr->closest_seed_id)
            std::cerr << "mega warning [this p] and [this r]" << std::endl;
#endif
          best_g = (lambda_r < lambda_p)?cr:cp; // this is shady... fixme...?
        }
        else // !p && q && r
        {
          best_g = cp; // we intersect [this,q] [this,r] therefore 'this' is on the side of p
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
#if (verbose > 10)
          if(cq->closest_seed_id != cr->closest_seed_id)
            std::cerr << "mega warning [this q] and [this r]" << std::endl;
#endif
          best_g = (lambda_r < lambda_q)?cr:cq; // this is shady... fixme...?
        }
        else // p && !q && r
        {
          best_g = cq;
        }
      }
      else // q
      {
        if(!intersect_r) // p && q && !r
        {
          best_g = cr;
        }
        else // p && q && r
        {
          // gotta choose between p & q & r...
#if (verbose > 10)
          if(cp->closest_seed_id != cq->closest_seed_id ||
             (cr && cp->closest_seed_id != cr->closest_seed_id) ||
             (cr && cq->closest_seed_id != cr->closest_seed_id))
            std::cerr << "ultra mega warning [this p], [this q] and [this r]" << std::endl;
#endif
          best_g = (lambda_q < lambda_p)?((lambda_r < lambda_q)?cr:cq):((lambda_r < lambda_p)?cr:cp);
        }
      }
    }

    CGAL_assertion(best_g);
    this->closest_seed_id = best_g->closest_seed_id;
    this->ancestor = best_g;
  }

  FT compute_closest_seed_3D(const Konukoglu_canvas_point* cp,
                             const Konukoglu_canvas_point* cq,
                             const Konukoglu_canvas_point* cr)
  {
    FT cp_d = cp->distance_to_closest_seed;
    FT cq_d = cq->distance_to_closest_seed;
    FT cr_d = cr->distance_to_closest_seed;

#if (verbose > 15)
    std::cout << "compute closest 3D " << this->index;
    std::cout << " to " << cp->index << ", " << cq->index << " & " << cr->index << " ";
    std::cout << " ds: " << cp_d << " " << cq_d << " " << cr_d;
    std::cout << " resp. " << cp->closest_seed_id << " ";
    std::cout << cq->closest_seed_id << " " << cr->closest_seed_id << std::endl;
#endif

    CGAL_assertion((cp->state == KNOWN || cp->state == CHANGED) &&
                   (cq->state == KNOWN || cq->state == CHANGED) &&
                   (cr->state == KNOWN || cr->state == CHANGED));
    bool changed = false;

    FT d = 0.;
#if 0
    bool is_3D_optimal = compute_min_distance_3D(cp, cq, cr, d);

    // if false, the solution is actually on one of the faces of the tetrahedron
    // so we call the 2D...
    if(!is_3D_optimal)
    {
      FT lambda;
      FT d_f1 = compute_min_distance_2D(cp, cq, lambda);
      FT d_f2 = compute_min_distance_2D(cp, cr, lambda);
      FT d_f3 = compute_min_distance_2D(cq, cr, lambda);

      d = (std::min)((std::min)(d_f1, d_f2), d_f3);
      std::cout << "3D wasn't optimal so 2D was called and found: " << d << std::endl;
    }
#endif

#if 1//def BRUTE_FORCE_CHECK_OPTIMAL_P
    // brute force check
    Vector3d v1, v2;
    v1(0) = cq->point.x() - cp->point.x();
    v1(1) = cq->point.y() - cp->point.y();
    v1(2) = cq->point.z() - cp->point.z();
    v2(0) = cr->point.x() - cp->point.x();
    v2(1) = cr->point.y() - cp->point.y();
    v2(2) = cr->point.z() - cp->point.z();

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
    //#pragma omp parallel for (don't forget the pragma critical if you uncomment here
    for(std::size_t nij = 0; nij <= tot; ++nij)
    {
      std::size_t ni = nij / (pt_n+1);
      std::size_t nj = nij % (pt_n+1);

      if(ni+nj > pt_n)
        continue;

      FT i = ni * increment;
      FT j = nj * increment;

      Point_3 p(cp->point.x() + i*v1(0) + j*v2(0),
                cp->point.y() + i*v1(1) + j*v2(1),
                cp->point.z() + i*v1(2) + j*v2(2));

      // compute the bary weights
      FT lambda_p, lambda_q, lambda_r;
      Vector3d v3;
      v3(0) = p.x() - cp->point.x();
      v3(1) = p.y() - cp->point.y();
      v3(2) = p.z() - cp->point.z();

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
      FT dcs_at_p = lambda_p * cp_d + lambda_q * cq_d + lambda_r * cr_d;

      CGAL_assertion(lambda_p >= -1e-10 && lambda_p <= 1.+1e-10 &&
                     lambda_q >= -1e-10 && lambda_q <= 1.+1e-10 &&
                     lambda_r >= -1e-10 && lambda_r <= 1.+1e-10);

      Vector3d v;
      v(0) = p.x() - this->point.x();
      v(1) = p.y() - this->point.y();
      v(2) = p.z() - this->point.z();

      const Eigen::Matrix3d& m = this->metric.get_mat();
      FT neighbor_d = std::sqrt(v.transpose() * m * v);
      FT dp = dcs_at_p + neighbor_d;

#if (verbose > 25)
      std::cout << "Debug compute_closest_seed_3D" << std::endl;
      std::cout << "the three points : " << std::endl;
      std::cout << cp->point.x() << " " << cp->point.x() << " " << cp->point.z() << std::endl;
      std::cout << cq->point.y() << " " << cq->point.y() << " " << cq->point.z() << std::endl;
      std::cout << cr->point.z() << " " << cr->point.z() << " " << cr->point.z() << std::endl;
      std::cout << "ij:" << i << " " << j << " and p : " << p.x() << " " << p.y() << " " << p.z() << std::endl;
      std::cout << "barywhite: " << lambda_p << " " << lambda_q << " " << lambda_r << std::endl;
      std::cout << "check: ";
      std::cout << lambda_p*cp->point.x()+lambda_q*cq->point.x()+lambda_r*cr->point.x() << " ";
      std::cout << lambda_p*cp->point.y()+lambda_q*cq->point.y()+lambda_r*cr->point.y() << " ";
      std::cout << lambda_p*cp->point.z()+lambda_q*cq->point.z()+lambda_r*cr->point.z() << std::endl;
#endif

      if(dp < min_d)
      {
        min_id = nij;
        min_d = dp;
      }
    }

#if (verbose > 20)
    std::cout << "brute force calls min at : " << min_id << " min_d: " << min_d << std::endl;
#endif
#endif
    d = min_d;

#if (verbose > 15)
    std::cout << "min 3D; min : " << d;
    std::cout << " vs current best: " << this->distance_to_closest_seed << std::endl;
#endif

    // tolerance to ignore unsignificant changes
    if(this->distance_to_closest_seed - d > canvas->recursive_tolerance * d)
    {
      changed = true;
      this->distance_to_closest_seed = d;
      opt_neighbors[0] = cp;
      opt_neighbors[1] = cq;
      opt_neighbors[2] = cr;

      // below is useless fixme
      this->closest_seed_id = static_cast<std::size_t>(-1);
      this->ancestor = cp;

#if (verbose > 15)
      std::cout << "3D new best for " << this->index
                << " : " << this->distance_to_closest_seed << std::endl;
#endif
    }
    return changed;
  }

  boost::array<const Konukoglu_canvas_point*, 4>
  get_adjacent_neighbors(const Konukoglu_canvas_point* cp,
                         const Konukoglu_canvas_point* cq)
  {
    // Find, in the neighbors of p, the points adjacent to q
    // for example, if q = right, the adjacent neighbors are back, above, front, below

    CGAL_assertion(cp && cq && cp != cq);
    boost::array<const Konukoglu_canvas_point*, 4> adj_n;

    // an annoying case is the border of the grid: not all neighbors are != NULL
    // so if we want to find the adjacent neighbors, we have to know what is the
    // neighbor opposite of cq (AND THAT OPPOSITE NEIGHBOR VERY WELL COULD BE 'NULL'!!)
    // we therefore use simply the index of the neighbors: we know the opposite couples:
    // 0-5, 1-3, and 2-4

    const Neighbors& ns = cp->neighbors; // just to get a shorter name...

    // first of all, we need to find the index of cq in cp->neighbors
    std::size_t pos_cq = 6; // 6 so that it will fail an assert later if cq not found
    for(std::size_t i=0; i<ns.size(); ++i)
      if(ns[i] == cq)
        pos_cq = i;

    // hardcoding all the cases, because it's simpler to code (and to understand)
    // considering pairs change the order, which will change 'nij' (but it's not a big deal)
    if(pos_cq == 0 || pos_cq == 5) // we want 1 2 3 4
    {
      adj_n[0] = ns[1]; adj_n[1] = ns[2]; adj_n[2] = ns[3]; adj_n[3] = ns[4];
    }
    else if(pos_cq == 1 || pos_cq == 3) // we want 0 2 5 4
    {
      adj_n[0] = ns[0]; adj_n[1] = ns[2]; adj_n[2] = ns[5]; adj_n[3] = ns[4];
    }
    else if(pos_cq == 2 || pos_cq == 4) // we want 0 1 5 3
    {
      adj_n[0] = ns[0]; adj_n[1] = ns[1]; adj_n[2] = ns[5]; adj_n[3] = ns[3];
    }

    return adj_n;
  }

  bool compute_closest_seed(const Konukoglu_canvas_point* anc)
  {
    // returns true if we improved the distance

#if (verbose > 10)
    std::cout << "compute closest high level " << this->index;
    std::cout << " with ancestor: " << anc->index;
#endif
    bool changed = false;

    // get the neighbors of 'this' that are adjacent of ancestor
    // then 4 possible 3D faces: this-ancestor-2 from adj_n
    // if a tet fails, we have two possible 2D faces: this-ancestor-1 from tet
    // if no face is available, then compute 1D: this-ancestor
    boost::array<const Konukoglu_canvas_point*, 4> adj_n =
                                              get_adjacent_neighbors(this, anc);

#if (verbose > 10)
    std::cout << " neighbors: ";
    for(std::size_t i=0; i<4; ++i)
      if(adj_n[i])
        std::cout << adj_n[i]->index << " ";
      else
        std::cout << "NULL ";
    std::cout << std::endl;
#endif

    for(std::size_t i=0; i<adj_n.size(); ++i)
    {
      const Konukoglu_canvas_point* cp = adj_n[i];
      const Konukoglu_canvas_point* cq = adj_n[(i+1)%(adj_n.size())];

      if(cp && (cp->state == KNOWN || cp->state == CHANGED))
      {
        if(cq && (cq->state == KNOWN || cq->state == CHANGED))
        {
          // the tet 'this - ancestor - cp - cq'
          if(compute_closest_seed_3D(anc, cp, cq))
            changed = true;
        }
        else
        { // the face 'this - ancestor - cp'
          if(compute_closest_seed_2D(anc, cp))
              changed = true;
        }
      }
      else if(cq && (cq->state == KNOWN || cq->state == CHANGED))
      { // the face 'this - ancestor - cq'
        if(compute_closest_seed_2D(anc, cq))
          changed = true;
      }
      else
      {
        if(compute_closest_seed_1D(anc))
          changed = true;
      }
    }

#if (verbose > 10)
    std::cout << "End high level. color and value: " << this->closest_seed_id << " ";
    std::cout << this->distance_to_closest_seed << std::endl;
#endif

    return changed;
  }

  PQ_state update_neighbors_distances() const
  {
//    std::cout << "update neighbors of " << index << std::endl;
    PQ_state pqs_ret = NOTHING_TO_DO;

    // consider all the known neighbors and compute the value to that
    CGAL_assertion(this->state == KNOWN || this->state == CHANGED);

    for(std::size_t i=0; i<neighbors.size(); ++i)
    {
      Konukoglu_canvas_point* cp = neighbors[i]; // a neighbor of 'this'
      if(!cp)
        continue;

#ifdef USE_RECURSIVE_UPDATES
      if(cp->state == KNOWN || cp->state == CHANGED) // that's the recursive part
      {
#if (verbose > 10)
        FT mem = cp->distance_to_closest_seed;
        std::size_t ancestor_mem = (cp->ancestor)?cp->ancestor->index:-1;
#endif
        // recompute the distance in the direction cp-'this'
        if(cp->compute_closest_seed(this))
        {
#if (verbose > 10)
          std::cout << "new value at " << cp->index << " : " << cp->distance_to_closest_seed;
          std::cout << " with ancestor: " << cp->ancestor->index;
          std::cout << " (mem: " << mem << " ancestor: " << ancestor_mem << ") ";
          std::cout << " diff: " << cp->distance_to_closest_seed-mem << std::endl;
#endif
          if(cp->state == KNOWN)
          {
            // if it's already 'changed' we only need to reorder the 'changed' queue
            // but we don't need to insert anything
            cp->state = CHANGED;
            canvas->changed_points.push_back(cp);
            std::push_heap(canvas->changed_points.begin(), canvas->changed_points.end(),
                           Canvas_point_comparer<Konukoglu_canvas_point>());
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
      if(cp->state == TRIAL)
      {
        // no need to insert in the trial queue since it already is in it
        if(cp->compute_closest_seed(this))
        {
          if(pqs_ret == NOTHING_TO_DO || pqs_ret == REBUILD_TRIAL)
            pqs_ret = REBUILD_TRIAL;
          else // pqs_ret == REBUILD_CHANGED || pqs_ret == REBUILD_BOTH
            pqs_ret = REBUILD_BOTH;
        }
      }
      else // cp->state == FAR
      {
        CGAL_assertion(cp->state == FAR);
        cp->compute_closest_seed(this); // always returns true here so no need to test it
        cp->state = TRIAL;
        canvas->trial_points.push_back(cp);
        std::push_heap(canvas->trial_points.begin(), canvas->trial_points.end(),
                       Canvas_point_comparer<Konukoglu_canvas_point>());
      }
    }
    return pqs_ret;
  }


  // really ugly hack to get the correct neighbor update function...
  // A proper fix would be to always have a const ref to the canvas in all the
  // point types but it's not worth the work (you need to change many typedefs
  // in all the classes & such)
  PQ_state update_neighbors_distances(std::vector<Konukoglu_canvas_point*>& /*trial*/) const
  {
    return update_neighbors_distances(); // this will call the function above with changed & trial
  }

  template<typename MF>
  Konukoglu_canvas_point(const Point_3&, const std::size_t, MF)
  {
    // this is a dummy because the initialize() of the base canvas wants to
    // construct a (derived) point with a metric argument, but in the Konu case,
    // our point is constructed with the canvas (canvas includes metric field)
    // We also don't need this constructor because initialize() is overloaded to
    // use the point constructor with the canvas, but the compiler doesn't see that...
  }

  Konukoglu_canvas_point(const Point_3& _point,
                         const std::size_t _index,
                         Canvas* _canvas)
    :
      Base(_point, _index, _canvas->mf),
      neighbors(),
      opt_neighbors(),
      canvas(_canvas)
  {
    for(std::size_t i=0; i<neighbors.size(); ++i)
      neighbors[i] = NULL;
  }
};

template<typename K, typename KExact, typename Canvas_point, typename Metric_field>
class Konukoglu_canvas :
    public Grid_canvas<K, KExact, Canvas_point, Metric_field>
{
public:
  typedef Grid_canvas<K, KExact, Canvas_point, Metric_field> Base;

  typedef typename Base::FT                                  FT;
  typedef typename Base::Point_3                             Point_3;

  typedef typename Base::Canvas_point_vector                 Canvas_point_vector;

  // a priority queue specific to this algorithm
  Canvas_point_vector changed_points;

  // tolerance on how large of an update the new value needs to be
  FT recursive_tolerance;

  // the exact same function exists in the base class of the canvas, but since
  // the points use the canvas in the constructor and need to know the fully
  // derived canvas, it's put here again so *this is a Konukoglu_canvas
  // Could do CRTP and used derived() etc... but it's a hassle.
  void initialize()
  {
#if (verbose > 5)
    std::cout << "canvas initialization" << std::endl;
#endif

    // create the canvas points
    for(unsigned int k=0; k<this->n; ++k)
    {
      for(unsigned int j=0; j<this->n; ++j)
      {
        for(unsigned int i=0; i<this->n; ++i)
        {
          Point_3 p(this->offset_x + i*this->step,
                    this->offset_y + j*this->step,
                    this->offset_z + k*this->step); // fill from bot left to top right
          Canvas_point cp(p, i + j*this->n + k*this->sq_n, this);
          this->points.push_back(cp);
        }
      }
    }

    // assign the neighbors
    for(unsigned int k=0; k<this->n; ++k)
    {
      for(unsigned int j=0; j<this->n; ++j)
      {
        for(unsigned int i=0; i<this->n; ++i)
        {
          std::size_t curr_id = i + j*this->n + k*this->sq_n;
          if(k != this->n-1) // there is a neighbor above
            this->points[curr_id].neighbors[0] = &(this->points[i + j*this->n + (k+1)*this->sq_n]);
          if(i != 0) // there is a neighbor left
            this->points[curr_id].neighbors[1] = &(this->points[i-1 + j*this->n + k*this->sq_n]);
          if(j != this->n-1) // there is a neighbor back
            this->points[curr_id].neighbors[2] = &(this->points[i + (j+1)*this->n + k*this->sq_n]);
          if(i != this->n-1) // there is a neighbor right
            this->points[curr_id].neighbors[3] = &(this->points[i+1 + j*this->n + k*this->sq_n]);
          if(j != 0) // there is a neighbor front
            this->points[curr_id].neighbors[4] = &(this->points[i + (j-1)*this->n + k*this->sq_n]);
          if(k != 0) // there is a neighbor below
            this->points[curr_id].neighbors[5] = &(this->points[i + j*this->n + (k-1)*this->sq_n]);
        }
      }
    }
#if (verbose > 5)
    std::cout << "neighbors assigned" << std::endl;
#endif

    Base::compute_canvas_bbox();
    this->seeds.initialize_seeds();
    Base::locate_seeds_on_canvas();

#if (verbose > 5)
    std::cout << "canvas initialized" << std::endl;
#endif
  }

  void determine_ancestors()
  {
    std::clock_t start = std::clock();

#if (verbose > 0)
    std::cout << "Determine ancestors" << std::endl;
#endif
    CGAL_assertion(changed_points.empty() && this->trial_points.empty());

    // TODO do something smarter when spreading a new color in an already colored
    // grid... One can probably just start from the new seed WITHOUT RESETING EVERYTHING
    // and then spread and only consider stuff if one of the possible ancestors
    // has the new color something like that...

    Canvas_point_vector known_points;
    for(std::size_t i=0; i<this->points.size(); ++i)
    {
      // reset the colors (but not the values!)

      // fixme, just leave the colors uninitialized in the first pass and
      // there's no need to reset them here...

      Canvas_point* cp = &(this->points[i]);
      cp->closest_seed_id = -1;
      cp->ancestor = NULL;

      CGAL_assertion(cp->state == KNOWN);

      known_points.push_back(cp);
    }
    std::make_heap(known_points.begin(), known_points.end(),
                   Canvas_point_comparer<Canvas_point>());

    // only re initialize the closest_seed_id for each point-seed
    for(std::size_t i=0; i<this->seeds.size(); ++i)
    {
      // if you use an 8 pts initilization for a seed, you need to modify here too...
      const Point_3& p = this->seeds[i];
      int index_x = std::floor((p.x() - this->offset_x) / this->step);
      int index_y = std::floor((p.y() - this->offset_y) / this->step);
      int index_z = std::floor((p.z() - this->offset_z) / this->step);
      Canvas_point* cp = &(this->points[index_z * this->sq_n + index_y * this->n + index_x]);
      cp->closest_seed_id = i;
    }

    bool is_kp_empty = known_points.empty();
    while(!is_kp_empty)
    {
      Canvas_point* cp = known_points.front();
      std::pop_heap(known_points.begin(), known_points.end(),
                    Canvas_point_comparer<Canvas_point>());
      known_points.pop_back();

      cp->determine_ancestor();
      is_kp_empty = known_points.empty();
    }

    std::cerr << "End of determine ancestors. time: ";
    std::cerr << ( std::clock() - start ) / (double) CLOCKS_PER_SEC << std::endl;
  }

  void paint()
  {
    std::clock_t start = std::clock();
#if (verbose > 0)
    std::cout << "main loop" << std::endl;
#endif

    bool is_cp_empty = changed_points.empty();
    bool is_t_empty = this->trial_points.empty();

    Canvas_point* cp;
    while(!is_cp_empty || !is_t_empty)
    {
#if (verbose > 5)
      std::cout << "Queue sizes. Trial: " << this->trial_points.size() << " Changed: " << changed_points.size() << std::endl;
#endif

#if (verbose > 20)
      std::cout << "changed heap: " << std::endl;
      typename std::vector<Canvas_point*>::iterator it = changed_points.begin();
      for (; it != changed_points.end(); ++it)
        std::cout << (*it)->index << " " << (*it)->distance_to_closest_seed << std::endl;
      std::cout << std::endl;

      std::cout << "trial heap: " << std::endl;
      it = this->trial_points.begin();
      for (; it != this->trial_points.end(); ++it)
        std::cout << (*it)->index << " " << (*it)->distance_to_closest_seed << std::endl;
      std::cout << std::endl;
#endif
      if(!is_cp_empty)
      {
        cp = changed_points.front();
        CGAL_assertion(cp && cp->state == CHANGED);
        std::pop_heap(changed_points.begin(), changed_points.end(),
                      Canvas_point_comparer<Canvas_point>());
        changed_points.pop_back();
      }
      else // !is_t_empty
      {
        cp = this->trial_points.front();
        CGAL_assertion(cp && cp->state == TRIAL);
        std::pop_heap(this->trial_points.begin(), this->trial_points.end(),
                      Canvas_point_comparer<Canvas_point>());
        this->trial_points.pop_back();
      }

#if (verbose > 5)
      std::cout << "picked n° " << cp->index << " (" << cp->point.x() << " " << cp->point.y() << " " << cp->point.z() << ")";
      std::cout << "at distance : " << cp->distance_to_closest_seed << " from " << cp->closest_seed_id << std::endl;
#endif

      cp->state = KNOWN;
      PQ_state pqs = cp->update_neighbors_distances();

      if(pqs == REBUILD_TRIAL || pqs == REBUILD_BOTH)
        std::make_heap(this->trial_points.begin(), this->trial_points.end(),
                       Canvas_point_comparer<Canvas_point>());
      if(pqs == REBUILD_CHANGED || pqs == REBUILD_BOTH)
        std::make_heap(changed_points.begin(), changed_points.end(),
                       Canvas_point_comparer<Canvas_point>());

      is_cp_empty = changed_points.empty();
      is_t_empty = this->trial_points.empty();
    }

    std::cerr << "End of geo_grid_loop. time: ";
    std::cerr << ( std::clock() - start ) / (double) CLOCKS_PER_SEC << std::endl;

    determine_ancestors();
  }

  void debug()
  {
#ifdef USE_RECURSIVE_UPDATES
    CGAL_assertion(this->trial_points.empty() && changed_points.empty());

    std::cout << "BRUTE FORCE CHECK THAT WE ARE REALLY FINISHED : " << std::endl;

    // it would probably be best to consider all the possible 2D cases here...
    // but more generally, 3D being weaker than 2D is a problem to be fixed fixme

    // detail of the problem :
    // 2D can be solved analytically while 3D cannot. Therefore we might sometimes
    // lose the exact result when we 'unlock' the 3D if the min is achieved on a face

    for(std::size_t i=0; i<this->points.size(); ++i)
    {
      Canvas_point* cp = &(this->points[i]);
      std::cout << "point " << i << " min distance is supposedly: ";
      std::cout << cp->distance_to_closest_seed << std::endl;
      typename Base::Neighbors::const_iterator it = cp->neighbors.begin(),
                                               iend = cp->neighbors.end();
      for(; it!=iend; ++it)
      {
        const Canvas_point* cq = *it;
        if(cq)
          CGAL_assertion(!cp->compute_closest_seed(cq));
      }
    }
#endif
  }

  Konukoglu_canvas(const std::string& _canvas_str,
                   const std::string& _seeds_str,
                   const Point_3& _center,
                   const FT _side,
                   const std::size_t points_per_side,
                   const std::size_t _max_seeds_n,
                   const std::size_t _n_refine,
                   const Metric_field _mf,
                   const FT _tolerance)
    :
      Base(_canvas_str, _seeds_str,
           _center, _side, points_per_side,
           _max_seeds_n, _n_refine,
           _mf),
      changed_points(),
      recursive_tolerance(_tolerance)
  { }
};

} // namespace Anisotropic_mesh_3
} // namespace CGAL

using namespace CGAL::Anisotropic_mesh_3;

int main(int, char**)
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel      K;
  typedef CGAL::Exact_predicates_exact_constructions_kernel        KExact;

  typedef typename K::FT                                           FT;
  typedef typename K::Point_3                                      Point_3;

  typedef typename CGAL::Anisotropic_mesh_3::Euclidean_metric_field<K>* MF;
  //typedef typename CGAL::Anisotropic_mesh_3::Custom_metric_field<K>* MF;

  typedef Konukoglu_canvas_point<K, KExact, MF>                    Konukoglu_canvas_point;
  typedef Konukoglu_canvas<K, KExact, Konukoglu_canvas_point, MF>  Canvas;

  std::cout.precision(17);
//  std::freopen("log.txt", "w", stdout);

  std::clock_t start = std::clock();
  double duration;

  MF mf = new Euclidean_metric_field<K>(1., 1., 3.);
//  MF mf = new Custom_metric_field<K>();

  const std::string canvas_str = "grid";
  const std::string seeds_str = "base_mesh.mesh";
  std::size_t max_seeds_n = 1;
  std::size_t n_refine = 0;

  // canvas geometry
  Point_3 center(1., 1., 1.);
  const FT canvas_side = 2.;
  FT points_per_side = 50.; // number of points per side of the canvas

  // how big of an update the new value needs to be
  FT recursive_tolerance = 1e-5;

  Canvas canvas(canvas_str, seeds_str,
                center, canvas_side, points_per_side,
                max_seeds_n, n_refine,
                mf,
                recursive_tolerance);

  canvas.initialize();
  canvas.paint();
  canvas.refine();
  canvas.output_canvas_data_and_dual(canvas_str + "_tr");
//  canvas.cosphericity_sweeper();


  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cerr << "Total: " << duration << std::endl;
}

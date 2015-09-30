#ifndef CGAL_ANISOTROPIC_MESH_3_KONUKOGLU_POINT_H
#define CGAL_ANISOTROPIC_MESH_3_KONUKOGLU_POINT_H

#include <CGAL/Canvas/canvas_point.h>

#include <boost/array.hpp>
#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>

#define USE_RECURSIVE_UPDATES

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K, typename KExact, typename Metric_field>
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
  typedef Konukoglu_canvas<K, KExact, Metric_field>           Canvas;

  typedef typename Base::FT                                   FT;
  typedef typename Base::Point_3                              Point_3;
  typedef typename Base::Vector3d                             Vector3d;

  typedef boost::array<const Konukoglu_canvas_point*, 3>      Optimal_neighbors;

// Neighbors
  Neighbors neighbors;

  // the neighbors of 'this' for which we have reached the lowest distance_to_closest_seed
  // this is a 3-array with at least 1 non-null pointer at opt_neighbors[0]
  Optimal_neighbors opt_neighbors;

  // we need a lot of info from the grid with this algorithm...
  Canvas* canvas; // forced to take a pointer and not a ref since we copy points

  bool compute_closest_seed_1D(const Konukoglu_canvas_point* cp)
  {
#if (verbosity > 15)
    std::cout << "compute closest 1D " << this->index() << " to " << cp->index();
    std::cout << " ds " << cp->distance_to_closest_seed();
    std::cout << " resp. " << cp->closest_seed_id() << std::endl;
#endif
    CGAL_assertion(cp->state() == KNOWN || cp->state() == CHANGED);
    bool changed = false;

    Vector3d v;
    v(0) = cp->point().x() - this->point().x();
    v(1) = cp->point().y() - this->point().y();
    v(2) = cp->point().z() - this->point().z();

    const Eigen::Matrix3d& m = this->metric().get_mat();
    FT neighbor_d = std::sqrt(v.transpose() * m * v);
    FT dcs_at_cp = cp->distance_to_closest_seed();
    FT d = dcs_at_cp + neighbor_d;

#if (verbosity > 15)
    std::cout << "min 1D; min : " << d;
    std::cout << " vs current best: " << this->distance_to_closest_seed() << std::endl;
#endif

    if(this->distance_to_closest_seed() - d > canvas->recursive_tolerance * d)
    {
      changed = true;
      this->distance_to_closest_seed() = d;

      opt_neighbors[0] = cp;
      opt_neighbors[1] = NULL; // a _1D might become better than a _2D or _3D
      opt_neighbors[2] = NULL; // so we have to reset the others

      // below is useless fixme
      this->closest_seed_id() = static_cast<std::size_t>(-1);
      this->m_ancestor = cp;

#if (verbosity > 15)
      std::cout << "1D new best for " << this->index() << " : "
                << this->distance_to_closest_seed() << std::endl;
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
    CGAL_assertion((cp->state() == KNOWN || cp->state() == CHANGED) &&
                   (cq->state() == KNOWN || cq->state() == CHANGED) );

    // in 2D, we can solve the minimization analytically.
    // We're solving :
    // min_{p\in [0,1]} f(p) = F(x)*p + F(y)*(1-p) + sqrt(v(p)^t * M * v(p))
    // vp = p*cp+(1-p)*cq
    // local minimum is given by f'(p_0) = 0 _but_ we have to compare f(p_0)
    // with f(0) and f(1) as we're only interested in solutions that live in [0,1]
    // If there's no solution in [0,1], then there's no (local) minimum in [0,1]
    // and the minimum is obtained at 0 or at 1.

    const FT dcs_at_cp = cp->distance_to_closest_seed();
    const FT dcs_at_cq = cq->distance_to_closest_seed();
    Vector3d vp, vq;
    vp(0) = cp->point().x() - this->point().x();
    vp(1) = cp->point().y() - this->point().y();
    vp(2) = cp->point().z() - this->point().z();
    vq(0) = cq->point().x() - this->point().x();
    vq(1) = cq->point().y() - this->point().y();
    vq(2) = cq->point().z() - this->point().z();

    const Eigen::Matrix3d& m = this->metric().get_mat();
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
#if (verbosity > 15)
    FT cp_d = cp->distance_to_closest_seed();
    FT cq_d = cq->distance_to_closest_seed();
    std::cout << "compute closest 2D " << this->index() << " to ";
    std::cout << cp->index() << " & " << cq->index();
    std::cout << " ds: " << cp_d << " " << cq_d << " ";
    std::cout << " resp. " << cp->closest_seed_id() << " " << cq->closest_seed_id() << std::endl;
#endif

    CGAL_assertion((cp->state() == KNOWN || cp->state() == CHANGED) &&
                   (cq->state() == KNOWN || cq->state() == CHANGED) );
    bool changed = false;

    FT p = 0.;
    const FT d = compute_min_distance_2D(cp, cq, p);

#if (verbosity > 15)
    std::cout << "min 2D; min : " << d << " at " << p;
    std::cout << " vs current best: " << this->distance_to_closest_seed() << std::endl;
#endif

    // tolerance to ignore unsignificant changes
    if(this->distance_to_closest_seed() - d > canvas->recursive_tolerance * d)
    {
      changed = true;
      this->distance_to_closest_seed() = d;

      opt_neighbors[0] = cp;
      opt_neighbors[1] = cq;
      opt_neighbors[2] = NULL;

      // below is useless fixme
      this->closest_seed_id() = static_cast<std::size_t>(-1);
      this->m_ancestor = cp;

#if (verbosity > 15)
      std::cout << "2D new best for " << this->index() << " : "
                << this->distance_to_closest_seed() << std::endl;
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


    const Eigen::Matrix3d& m = this->metric().get_mat();
    FT cp_d = cp->distance_to_closest_seed();
    FT cq_d = cq->distance_to_closest_seed();
    FT cr_d = cr->distance_to_closest_seed();

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

#if (verbosity > 20)
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
    FT cp_d = cp->distance_to_closest_seed();
    FT cq_d = cq->distance_to_closest_seed();
    FT cr_d = cr->distance_to_closest_seed();

    FT l = canvas->step; // all edge lengths are the same at the moment
    FT lden = 1./l;

    Vector3d vp, vq, vr; // note that these vectors "point" towards 'this'
    vp(0) = this->point().x() - cp->point().x();
    vp(1) = this->point().y() - cp->point().y();
    vp(2) = this->point().z() - cp->point().z();
    vq(0) = this->point().x() - cq->point().x();
    vq(1) = this->point().y() - cq->point().y();
    vq(2) = this->point().z() - cq->point().z();
    vr(0) = this->point().x() - cr->point().x();
    vr(1) = this->point().y() - cr->point().y();
    vr(2) = this->point().z() - cr->point().z();

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

#if (verbosity > 20)
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
    if(!cp || cp->closest_seed_id() == static_cast<std::size_t>(-1))
      return lambda;

    FT cp_d = cp->distance_to_closest_seed();
    const Eigen::Matrix3d& mp = cp->metric().get_mat();
    const Eigen::Matrix3d& m = this->metric().get_mat();

    const FT me1 = std::sqrt(m(0,0)), me2 = std::sqrt(m(1,1)), me3 = std::sqrt(m(2,2));
    const FT pe1 = std::sqrt(mp(0,0)), pe2 = std::sqrt(mp(1,1)), pe3 = std::sqrt(mp(2,2));

    if(cp->closest_seed_id() != static_cast<std::size_t>(-1))
    {
      if(cp->index() == this->index() + 1 || cp->index() == this->index() - 1) // p pight/left of 'this'
        lambda = (cp_d - this->distance_to_closest_seed() + pe1 * canvas->step) / (me1 + pe1);
      else if(cp->index() == this->index() + canvas->n || cp->index() == this->index() - canvas->n) // p back/fpont of 'this'
        lambda = (cp_d - this->distance_to_closest_seed() + pe2 * canvas->step) / (me2 + pe2);
      else if(cp->index() == this->index() + canvas->sq_n || cp->index() == this->index() - canvas->sq_n) // p above/below 'this'
        lambda = (cp_d - this->distance_to_closest_seed() + pe3 * canvas->step) / (me3 + pe3);
      else
        CGAL_assertion(false);
    }

    // debug & verifications
    Vector3d p_to_mp, this_to_mp;
    FT mpx = this->point().x() + lambda * (cp->point().x() - this->point().x()) / canvas->step;
    FT mpy = this->point().y() + lambda * (cp->point().y() - this->point().y()) / canvas->step;
    FT mpz = this->point().z() + lambda * (cp->point().z() - this->point().z()) / canvas->step;
    p_to_mp(0) = cp->point().x() - mpx;
    p_to_mp(1) = cp->point().y() - mpy;
    p_to_mp(2) = cp->point().z() - mpz;
    this_to_mp(0) = this->point().x() - mpx;
    this_to_mp(1) = this->point().y() - mpy;
    this_to_mp(2) = this->point().z() - mpz;
    FT d_p_to_mp = std::sqrt(p_to_mp.transpose() * mp * p_to_mp);
    FT d_this_to_mp = std::sqrt(this_to_mp.transpose() * m * this_to_mp);
    FT diffp = this->distance_to_closest_seed() + d_this_to_mp - (cp_d + d_p_to_mp);

    if(lambda >= -1e-10 && lambda <= canvas->step + 1e-10)
      CGAL_assertion(std::abs(diffp) < 1e-10);

    return lambda;
  }

  void determine_ancestor_2D(const Konukoglu_canvas_point* cp,
                             const Konukoglu_canvas_point* cq,
                             const FT lambda_p, const FT lambda_q)
  {
    CGAL_assertion(cp && cp->closest_seed_id() != static_cast<std::size_t>(-1) &&
                   cq && cq->closest_seed_id() != static_cast<std::size_t>(-1));

#if (verbosity > 15)
    std::cout << "determine ancestor 2D " << cp->index() << " " << cq->index();
    std::cout << " colors: " << cp->closest_seed_id() << " " << cq->closest_seed_id() << std::endl;
#endif

    if(cp->closest_seed_id() == cq->closest_seed_id())
    {
      this->closest_seed_id() = cp->closest_seed_id(); // this probably creates problems in the ancestor tree...
      this->m_ancestor = cp;
    }

    // we have to use the lambdas to decide which one to use as ancestor!
    bool intersect_p = lambda_p >= -1e-2 && lambda_p <= canvas->step + 1e-10;
    bool intersect_q = lambda_q >= -1e-2 && lambda_q <= canvas->step + 1e-10;
    //fixme -1e-2 is absurdly big for a tolerance...

    const Konukoglu_canvas_point* best_g = NULL;
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

    this->closest_seed_id() = best_g->closest_seed_id();
    this->m_ancestor = best_g;

#if (verbosity > 15)
    std::cout << "ancestor shen 2D @ " << this->index()
              << " color is : " << this->closest_seed_id() << std::endl;
#endif
  }

  void determine_ancestor()
  {
#if (verbosity > 10)
    std::cout << "determine ancestor for " << this->index()
              << " d: " << this->distance_to_closest_seed() << std::endl;
#endif

    if(this->closest_seed_id() != static_cast<std::size_t>(-1))
      return;
    // might be better to use a boolean like 'is_the_canvas_point_closest_to_a_seed' so that
    // we can still spread a color on top of another while refining

    CGAL_assertion(opt_neighbors[0]);
    const Konukoglu_canvas_point* cp = opt_neighbors[0];
    const Konukoglu_canvas_point* cq = opt_neighbors[1];
    const Konukoglu_canvas_point* cr = opt_neighbors[2];

#if (verbosity > 10)
    std::cout << " " << (cq?(cr?"3":"2"):"1") << "D" << std::endl;
#endif

    // cp         non NULL => best is achieved from 1D
    // cp, cq     non NULL => best is achieved from 2D
    // cp, cq, cr non NULL => best is achieved from 3D

    bool is_init_p = cp->closest_seed_id() != static_cast<std::size_t>(-1);
    bool is_init_q = cq?(cq->closest_seed_id() != static_cast<std::size_t>(-1)):false;
    bool is_init_r = cr?(cr->closest_seed_id() != static_cast<std::size_t>(-1)):false;

#if (verbosity > 20)
    std::cout << "ini bools: " << is_init_p << " " << is_init_q << " " << is_init_r;
    std::cout << " with colors: " << cp->closest_seed_id();

    if(cq)
      std::cout << " " << cq->closest_seed_id() << " ";
    else
      std::cout << " -1 ";

    if(cr)
      std::cout << cr->closest_seed_id() << std::endl;
    else
      std::cout << "-1" << std::endl;
#endif

    CGAL_assertion((is_init_p || is_init_q || is_init_r) &&
                   "cp & cq & cr have uninitialized seeds...");

    // quick filters
    if(is_init_p && !is_init_q && !is_init_r)
    {
      this->closest_seed_id() = cp->closest_seed_id();
      this->m_ancestor = cp;
      return;
    }
    else if(!is_init_p && is_init_q && !is_init_r)
    {
      this->closest_seed_id() = cq->closest_seed_id();
      this->m_ancestor = cq;
      return;
    }
    else if(!is_init_p && !is_init_q && is_init_r)
    {
      this->closest_seed_id() = cr->closest_seed_id();
      this->m_ancestor = cr;
      return;
    }
    else if(is_init_p && is_init_q && is_init_r &&
            cp->closest_seed_id() == cq->closest_seed_id() &&
            cp->closest_seed_id() == cr->closest_seed_id())
    {
      this->closest_seed_id() = cp->closest_seed_id();
      this->m_ancestor = cp;
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

#if (verbosity > 20)
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
#if (verbosity > 10)
          if(cp->closest_seed_id() != cq->closest_seed_id())
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
#if (verbosity > 10)
          if(cp->closest_seed_id() != cr->closest_seed_id())
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
#if (verbosity > 10)
          if(cq->closest_seed_id() != cr->closest_seed_id())
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
#if (verbosity > 10)
          if(cp->closest_seed_id() != cq->closest_seed_id() ||
             (cr && cp->closest_seed_id() != cr->closest_seed_id()) ||
             (cr && cq->closest_seed_id() != cr->closest_seed_id()))
            std::cerr << "ultra mega warning [this p], [this q] and [this r]" << std::endl;
#endif
          best_g = (lambda_q < lambda_p)?((lambda_r < lambda_q)?cr:cq):((lambda_r < lambda_p)?cr:cp);
        }
      }
    }

    CGAL_assertion(best_g);
    this->closest_seed_id() = best_g->closest_seed_id();
    this->m_ancestor = best_g;
  }

  FT compute_closest_seed_3D(const Konukoglu_canvas_point* cp,
                             const Konukoglu_canvas_point* cq,
                             const Konukoglu_canvas_point* cr)
  {
    FT cp_d = cp->distance_to_closest_seed();
    FT cq_d = cq->distance_to_closest_seed();
    FT cr_d = cr->distance_to_closest_seed();

#if (verbosity > 15)
    std::cout << "compute closest 3D " << this->index();
    std::cout << " to " << cp->index() << ", " << cq->index() << " & " << cr->index() << " ";
    std::cout << " ds: " << cp_d << " " << cq_d << " " << cr_d;
    std::cout << " resp. " << cp->closest_seed_id() << " ";
    std::cout << cq->closest_seed_id() << " " << cr->closest_seed_id() << std::endl;
#endif

    CGAL_assertion((cp->state() == KNOWN || cp->state() == CHANGED) &&
                   (cq->state() == KNOWN || cq->state() == CHANGED) &&
                   (cr->state() == KNOWN || cr->state() == CHANGED));
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
    v1(0) = cq->point().x() - cp->point().x();
    v1(1) = cq->point().y() - cp->point().y();
    v1(2) = cq->point().z() - cp->point().z();
    v2(0) = cr->point().x() - cp->point().x();
    v2(1) = cr->point().y() - cp->point().y();
    v2(2) = cr->point().z() - cp->point().z();

    // loops originally are :
    //for(std::size_t ni=0; ni<=pt_n; ni++)
    //  for(std::size_t nj=0; nj<=pt_n-ni; nj++)

    // To parallelize these two loops, we create a single loop and to make it
    // less tedious to count, we create a grid of the parallelogram and
    // ignore the pts that are outside the triangle (that is, when ni+nj > pt_n)

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

      Point_3 p(cp->point().x() + i*v1(0) + j*v2(0),
                cp->point().y() + i*v1(1) + j*v2(1),
                cp->point().z() + i*v1(2) + j*v2(2));

      // compute the bary weights
      FT lambda_p, lambda_q, lambda_r;
      Vector3d v3;
      v3(0) = p.x() - cp->point().x();
      v3(1) = p.y() - cp->point().y();
      v3(2) = p.z() - cp->point().z();

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
      v(0) = p.x() - this->point().x();
      v(1) = p.y() - this->point().y();
      v(2) = p.z() - this->point().z();

      const Eigen::Matrix3d& m = this->metric().get_mat();
      FT neighbor_d = std::sqrt(v.transpose() * m * v);
      FT dp = dcs_at_p + neighbor_d;

#if (verbosity > 25)
      std::cout << "Debug compute_closest_seed_3D" << std::endl;
      std::cout << "the three points : " << std::endl;
      std::cout << cp->point() << std::endl;
      std::cout << cq->point() << std::endl;
      std::cout << cr->point() << std::endl;
      std::cout << "ij:" << i << " " << j << " and p : " << p << std::endl;
      std::cout << "barywhite: " << lambda_p << " " << lambda_q << " " << lambda_r << std::endl;
      std::cout << "check: ";
      std::cout << lambda_p*cp->point().x()+lambda_q*cq->point().x()+lambda_r*cr->point().x() << " ";
      std::cout << lambda_p*cp->point().y()+lambda_q*cq->point().y()+lambda_r*cr->point().y() << " ";
      std::cout << lambda_p*cp->point().z()+lambda_q*cq->point().z()+lambda_r*cr->point().z() << std::endl;
#endif

      if(dp < min_d)
      {
        min_id = nij;
        min_d = dp;
      }
    }

#if (verbosity > 20)
    std::cout << "brute force calls min at : " << min_id << " min_d: " << min_d << std::endl;
#endif
#endif
    d = min_d;

#if (verbosity > 15)
    std::cout << "min 3D; min : " << d;
    std::cout << " vs current best: " << this->distance_to_closest_seed() << std::endl;
#endif

    // tolerance to ignore unsignificant changes
    if(this->distance_to_closest_seed() - d > canvas->recursive_tolerance * d)
    {
      changed = true;
      this->distance_to_closest_seed() = d;
      opt_neighbors[0] = cp;
      opt_neighbors[1] = cq;
      opt_neighbors[2] = cr;

      // below is useless fixme
      this->closest_seed_id() = static_cast<std::size_t>(-1);
      this->m_ancestor = cp;

#if (verbosity > 15)
      std::cout << "3D new best for " << this->index()
                << " : " << this->distance_to_closest_seed() << std::endl;
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

#if (verbosity > 10)
    std::cout << "compute closest high level " << this->index();
    std::cout << " with ancestor: " << anc->index();
#endif
    bool changed = false;

    // get the neighbors of 'this' that are adjacent of ancestor
    // then 4 possible 3D faces: this-ancestor-2 from adj_n
    // if a tet fails, we have two possible 2D faces: this-ancestor-1 from tet
    // if no face is available, then compute 1D: this-ancestor
    boost::array<const Konukoglu_canvas_point*, 4> adj_n =
                                              get_adjacent_neighbors(this, anc);

#if (verbosity > 10)
    std::cout << " neighbors: ";
    for(std::size_t i=0; i<4; ++i)
      if(adj_n[i])
        std::cout << adj_n[i]->index() << " ";
      else
        std::cout << "NULL ";
    std::cout << std::endl;
#endif

    for(std::size_t i=0; i<adj_n.size(); ++i)
    {
      const Konukoglu_canvas_point* cp = adj_n[i];
      const Konukoglu_canvas_point* cq = adj_n[(i+1)%(adj_n.size())];

      if(cp && (cp->state() == KNOWN || cp->state() == CHANGED))
      {
        if(cq && (cq->state() == KNOWN || cq->state() == CHANGED))
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
      else if(cq && (cq->state() == KNOWN || cq->state() == CHANGED))
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

#if (verbosity > 10)
    std::cout << "End high level. color and value: " << this->closest_seed_id() << " ";
    std::cout << this->distance_to_closest_seed() << std::endl;
#endif

    return changed;
  }

  PQ_state update_neighbors_distances() const
  {
//    std::cout << "update neighbors of " << index << std::endl;
    PQ_state pqs_ret = NOTHING_TO_DO;

    // consider all the known neighbors and compute the value to that
    CGAL_assertion(this->state() == KNOWN || this->state() == CHANGED);

    for(std::size_t i=0; i<neighbors.size(); ++i)
    {
      Konukoglu_canvas_point* cp = neighbors[i]; // a neighbor of 'this'
      if(!cp)
        continue;

#ifdef USE_RECURSIVE_UPDATES
      if(cp->state() == KNOWN || cp->state() == CHANGED) // that's the recursive part
      {
#if (verbosity > 10)
        FT mem = cp->distance_to_closest_seed();
        std::size_t ancestor_mem = (cp->ancestor())?cp->ancestor()->index():-1;
#endif
        // recompute the distance in the direction cp-'this'
        if(cp->compute_closest_seed(this))
        {
#if (verbosity > 10)
          std::cout << "new value at " << cp->index() << " : " << cp->distance_to_closest_seed();
          std::cout << " with ancestor: " << cp->ancestor()->index();
          std::cout << " (mem: " << mem << " ancestor: " << ancestor_mem << ") ";
          std::cout << " diff: " << cp->distance_to_closest_seed() - mem << std::endl;
#endif
          if(cp->state() == KNOWN)
          {
            // if it's already 'changed' we only need to reorder the 'changed' queue
            // but we don't need to insert anything
            cp->state() = CHANGED;
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
      if(cp->state() == TRIAL)
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
        CGAL_assertion(cp->state() == FAR);

        // note that cp->distance_to_closest_seed is not necessarily FT_inf here :
        // if we're refining, we've assigned FAR to all points after inserting a new
        // seed, therefore we must verify that compute_closest_seed is an update
        // before inserting it in the trial_queue
        if(cp->compute_closest_seed(this))
        {
          cp->state() = TRIAL;
          canvas->trial_points.push_back(cp);
          std::push_heap(canvas->trial_points.begin(), canvas->trial_points.end(),
                         Canvas_point_comparer<Konukoglu_canvas_point>());
        }
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

} // namespace Anisotropic_mesh_3
} // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_KONUKOGLU_POINT_H

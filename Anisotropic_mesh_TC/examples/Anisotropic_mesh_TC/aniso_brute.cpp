#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Regular_triangulation.h>
#include <CGAL/Regular_triangulation_euclidean_traits.h>

#include <CGAL/helpers/combinatorics_helper.h>

#include <CGAL/Combination_enumerator.h>

#include <CGAL/algorithm.h>
#include <CGAL/Timer.h>
#include <CGAL/assertions.h>

#include <Eigen/Dense>

#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>

#define dim2
const int d = 2;
const int D = 5;

typedef CGAL::Epick_d< CGAL::Dimension_tag<D> >                  K;
typedef K::Point_d                                               Point_d;
typedef typename K::FT                                           FT;
typedef CGAL::Regular_triangulation_euclidean_traits<K>          Traits;
typedef typename CGAL::Triangulation_data_structure<typename K::Dimension,
CGAL::Triangulation_vertex<Traits, std::size_t> >                TDS;
typedef CGAL::Regular_triangulation<Traits, TDS>                 RTriangulation;

typedef typename K::Vector_d                                     Vector;
typedef typename Eigen::Matrix<double, d, 1>                     Pointd;
typedef typename Eigen::Matrix<double, d, 1>                     Vectord;
typedef typename Eigen::Matrix<double, D, 1>                     VectorD;
typedef typename Eigen::Matrix<double, d, d>                     Matrixd;
typedef typename Eigen::Matrix<double, D, D>                     MatrixD;

typedef Ordered_simplex_base<d+1>                                Simplex;
typedef typename Simplex_unordered_set<d+1>::type                Simplex_set;

typedef typename TDS::Vertex                                     Vertex;
typedef typename TDS::Vertex_handle                              Vertex_handle;
typedef typename TDS::Face                                       Face;
typedef typename TDS::Full_cell_handle                           Full_cell_handle;
typedef typename TDS::Full_cell::Vertex_handle_iterator          Vertex_h_iterator;
typedef typename std::vector<Full_cell_handle>                   Full_cell_handle_vector;
typedef typename Full_cell_handle_vector::iterator               Full_cell_handle_iterator;

typedef typename RTriangulation::Point                           Point;
typedef typename RTriangulation::Bare_point                      Bare_point;
typedef typename RTriangulation::Weighted_point                  Weighted_point;
typedef typename RTriangulation::Finite_vertex_iterator          Finite_vertex_iterator;
typedef typename RTriangulation::Finite_full_cell_iterator       Finite_full_cell_iterator;
typedef typename RTriangulation::Finite_full_cell_const_iterator Finite_full_cell_const_iterator;
typedef typename RTriangulation::Finite_vertex_const_iterator    Finite_vertex_const_iterator;
typedef typename RTriangulation::Vertex_const_handle             Vertex_const_handle;

std::ostream& operator<<(std::ostream& out, const Bare_point& p)
{
  for(int i=0; i<D; ++i)
    out << p[i] << " ";
  out << std::endl;
  return out;
}

FT volume(const Simplex& s,
          const std::vector<Pointd>& points)
{
  assert(s.size() == d+1);
  FT den = 1./(FT) fact(d);
  Matrixd m;
  for(int i=0; i<d; ++i)
  {
    Pointd v0 = points[s[0]];
    for(int j=0; j<d; ++j)
    {
      Pointd vj = points[s[j+1]];
      m(i,j) = vj(i)-v0(i);
    }
  }
  return std::abs(den * m.determinant());
}

// check if a simplex in R^D is degenerated in R^d
bool is_simplex_degenerated_in_Rd(const Simplex& s,
                                  const std::vector<Pointd>& points)
{
  Matrixd m;

  for(int i=0; i<d; ++i)
    for(int j=0; j<d; ++j)
      m(i,j) = points[s[i+1]](j)-points[s[0]](j);

  return (std::abs(m.determinant()) < 1e-10);
}


FT norm(const Point_d& p) // squared norm, actually
{
  FT n = 0.;
  for(int i=0; i<D; ++i)
    n += p[i]*p[i];
  return n;
}

Bare_point construct_bp(const VectorD& v)
{
#ifdef dim1
  return K().construct_point_d_object()(v(0), v(1));
#else
#ifdef dim2
  return K().construct_point_d_object()(v(0), v(1), v(2), v(3), v(4));
#else
  return K().construct_point_d_object()(v(0), v(1), v(2), v(3), v(4), v(5),
                                        v(6), v(7), v(8));
#endif
#endif
}

Pointd get_pointd_from_RD(const Bare_point& p)
{
  Pointd pp;
  for(int i=0; i<d; ++i)
    pp(i) = p[i];
  return pp;
}

Matrixd get_metric(const Pointd& p)
{
  // fixme Use the proper metric files, this function is a bit pathetic...
  Matrixd m;

#ifdef dim1
  m(0,0) = p(0);
#elif defined(dim2)

#ifdef HYPER
  FT x = p[0];
  FT y = p[1];
  FT delta = 0.1;
  FT tanhder = tanh((2.0 * x - sin(5.0 * y)) / delta);
  tanhder = (1.0 - tanhder * tanhder) / delta;
  FT x1 = 2. * tanhder + 3.0 * x * x + y * y;
  FT y1 = -tanhder * cos(5.0 * y) * 5.0 + 2.0 * x * y;

  FT r = sqrt(x1 * x1 + y1 * y1);
  FT x2 = -y1 / r;
  FT y2 = x1 / r;
  FT l1 = std::sqrt(x1*x1 + y1*y1);
  FT l2 = std::sqrt(x2*x2 + y2*y2);

  Vectord v1 = (1./l1) * Vectord(x1, y1);
  Vectord v2 = (1./l2) * Vectord(x2, y2);

  Eigen::Matrix2d eigen_m;
  eigen_m(0,0) = v1.x();  eigen_m(0,1) = v2.x();
  eigen_m(1,0) = v1.y();  eigen_m(1,1) = v2.y();

  Eigen::Matrix2d eigen_diag = Eigen::Matrix2d::Zero();
  eigen_diag(0,0) = l1;
  eigen_diag(1,1) = l2;

  Eigen::Matrix2d eigen_mtransp = eigen_m.transpose();
  m = eigen_m * eigen_diag * eigen_mtransp;
  m = m.transpose() * m;
#else
  double h = 0.3;

  double x = p[0];
  double y = p[1];
  double r = std::sqrt(x*x + y*y);

  if(r<0.1) // avoids problem near 0 (where v1 & v2 are undefined)
  {
    m(0,0) = 1./(h*h); m(0,1) = 0.;
    m(1,0) = 0.; m(1,1) = 1./(h*h);
    return m;
  }

  x /= r;
  y /= r;

  double lambda = std::exp(-0.5*std::abs(r*r-1));
  double h1 = 0.1 + (1-lambda); // artificial smoothing
  double h2 = 1.;

  Vectord v1(x, y);
  Vectord v2(-y, x);

  FT l1 = 1./(h*h1);
  FT l2 = 1./(h*h2);
#endif

  Eigen::Matrix2d eigen_m;
  eigen_m(0,0) = v1.x();  eigen_m(0,1) = v2.x();
  eigen_m(1,0) = v1.y();  eigen_m(1,1) = v2.y();

  Eigen::Matrix2d eigen_diag = Eigen::Matrix2d::Zero();
  eigen_diag(0,0) = l1;
  eigen_diag(1,1) = l2;

  Eigen::Matrix2d eigen_mtransp = eigen_m.transpose();
  m = eigen_m * eigen_diag * eigen_mtransp;
  m = m.transpose() * m;

#else
  m(0,0) = 1.; m(0,1) = 0.; m(0,2) = 0.;
  m(1,0) = 0.; m(1,1) = 1.; m(1,2) = 0.;
  m(2,0) = 0.; m(2,1) = 0.; m(2,2) = 1.;
#endif
  return m;
}

VectorD to_Q(const Vectord& p) // from a point in R^d to R^D on Q
{
  VectorD p_on_Q;
  for(int i=0; i<d; ++i)
    p_on_Q(i) = p(i);

  int ind = d;
  for(int i=0; i<d; ++i)
    for(int j=i; j<d; ++j)
      p_on_Q(ind++) = p(i)*p(j);

  return p_on_Q;
}

Vectord from_Q(const VectorD& p_on_Q) // from R^D on Q to R^d
{
  Vectord p;
  for(int i=0; i<d; ++i)
    p(i) = p_on_Q(i);
  return p;
}

VectorD to_S(const Vectord& p) // from R^d to R^D on the metric surface
{
  Matrixd m = get_metric(p);

  std::cout << "p : " << p.transpose() << std::endl;
  Vectord p_bar = m*p;
  std::cout << "pbar : " << p_bar.transpose() << std::endl;
  VectorD p_on_S;
  for(int i=0;i<d; ++i)
    p_on_S(i) = p_bar(i);

  int ind = d;
  for(int i=0; i<d; ++i)
  {
    for(int j=i; j<d; ++j)
    {
      if(j==i)
        p_on_S(ind++) = -0.5*m(i,i);
      else
        p_on_S(ind++) = -m(i,j);
    }
  }

  return p_on_S;
}

Vectord from_S(const VectorD& p_on_S) // from R^D on the metric surface to R^d
{
  Vectord p;
  Matrixd m;
  int i=0, j=0;
  for(int k=d; k<D; k++)
  {
    m(i,j) = p_on_S(k);
    m(j,i) = p_on_S(k);
    j++;
    if(j==d)
    {
      i++;
      j=0;
    }
  }

  Vectord tmp;
  for(int i=0; i<d; ++i)
    tmp(i) = p_on_S(i);

  p = m.inverse()*tmp;
  return p;
}

Weighted_point compute_wpoint(const RTriangulation& rt,
                              const Pointd& p)
{
  Matrixd m = get_metric(p);
  VectorD p_hat = to_S(p);
  std::cout << "hatp :" << p_hat.transpose() << std::endl;

  FT w = p_hat.norm()*p_hat.norm() - p.transpose()*m*p;
  std::cout << "weight: " << w << std::endl;

  const Traits& traits = rt.geom_traits();
  Bare_point ret = construct_bp(p_hat);
  return traits.construct_weighted_point_d_object()(ret, w);
}

void generate_wpoints(const RTriangulation& rt,
                      std::vector<Pointd>& points,
                      std::vector<Weighted_point>& wpoints,
                      const unsigned int N)
{
  srand(0);
//  for(unsigned int i=0; i<N; ++i)
//  {
//    double x = 1.+5.*((double) rand() / (RAND_MAX));
//    double y = 1.+5.*((double) rand() / (RAND_MAX));
//    double z = 1.+5.*((double) rand() / (RAND_MAX));
//-----------------------------------------------------------------------------
//  for(int i=1; i<3; ++i)
//    for(int j=1; j<3; ++j)
//    {
//      double x=i;
//      double y=j;
//-----------------------------------------------------------------------------
//  Eigen::Vector3d xs; xs(0) = -1.; xs(1) = 1.; xs(2) = -1.;
//  Eigen::Vector3d ys; ys(0) = -1.; ys(1) = -1.; ys(2) = 1.;
//  for(int i=0; i<3; ++i)
//  {
//    double x = xs(i), y = ys(i);
//-----------------------------------------------------------------------------

  std::ifstream in("../../../../Anisotropic_mesh_2/examples/Anisotropic_mesh_2/build/bambimboum.mesh");
  std::string word;
  int useless, nv, dim;
  FT x;

  in >> word >> useless; //MeshVersionFormatted i
  in >> word >> dim; //Dimension d
  in >> word >> nv;
  CGAL_assertion(d == dim);

  std::cout << "input files has " << nv << " vertices" << std::endl;
  for(int i=0; i<nv; ++i)
  {
    Pointd p;
    for(int j=0; j<d; ++j)
    {
      in >> x;
      p(j) = x;
    }
    in >> useless; // don't care about the vertex's reference

    points.push_back(p);
    Weighted_point wp = compute_wpoint(rt, p);
    wpoints.push_back(wp);

    if(wpoints.size() >= N)
      break;
  }
}

bool is_support_intersecting(const Simplex& fs,
                             const std::vector<Pointd>& points,
                             const std::vector<Weighted_point>& wpoints,
                             Bare_point& sol)
{
  // intersection with the tangent plane at the corresponding point
  // on the paraboloid of a point of a simplex

  // system is given by
  //   d(p,p0) - w0 = d(p,pi) - wi
  //   p \in T_{P_0}

  MatrixD m = MatrixD::Zero();
  VectorD l;
  Weighted_point p0 = wpoints[fs[0]];
  FT d0 = norm(p0.point());
  VectorD p0_on_Q = to_Q(points[fs[0]]);

  std::cout << "p0 on Q: " << p0_on_Q.transpose() << std::endl;

  //first part of the matrix (lines 0-d) (dual equation)
  for(int i=0; i<d; ++i)
  {
    Weighted_point pi = wpoints[fs[i+1]];
    FT di = norm(pi.point());
    l(i) = d0 - p0.weight() - di + pi.weight();
    for(int j=0; j<D; ++j)
      m(i,j) = 2*(p0.point()[j] - pi.point()[j]);
  }

  //second part of the matrix (lines d-D) (tangent plane equation)
  for(int i=d; i<D; ++i)
  {
    l(i) = p0_on_Q(i);
    m(i,i) = -1;
  }

  int pos = d;
  for(int i=0; i<d; ++i)
  {
    for(int j=i; j<d; ++j)
    {
      if(i == j)
        m(pos, i) = 2*points[fs[0]](i);
      else
      {
        m(pos, i) = points[fs[0]](j);
        m(pos, j) = points[fs[0]](i);
      }
      pos++;
    }
  }

  std::cout << "m @ supportintersection:" << std::endl << m << std::endl;
  std::cout << "LHS " << l.transpose() << std::endl;
  FT det = m.determinant();
  std::cout << "The determinant of m is " << det << std::endl;

  if(std::abs(det) > 1e-10)
  {
    MatrixD mm1 = m.inverse();
    std::cout << "The inverse of m is:\n" << mm1 << std::endl;
    VectorD vsol = mm1*l;
    sol = construct_bp(vsol);
    std::cout << "sol: " << vsol.transpose() << std::endl;

    Traits t;
    typename K::Squared_distance_d csd = t.squared_distance_d_object();

    for(int i=0; i<(d+1); ++i)
    {
      std::size_t id = fs[i];
      const Weighted_point& wi = wpoints[id];
      FT di = csd(sol, wi.point()) - wi.weight();
      std::cout << "dist² sol " << id << " : "  << di << std::endl;
    }

    return true;
  }
  else
  {
    std::cout << "null det" << std::endl;
    return false;
  }
}

struct H // hyperplanes describing the dual and error computation
{
  const Simplex& fs;
  const std::vector<Weighted_point>& wpoints;

  Vectord operator()(const Vectord& s) const
  {
    VectorD s_on_Q = to_Q(s);

//std::cout << "eval H. Point: " << s_on_Q.transpose() << std::endl;
//std::cout << "from init: " << s.transpose() << std::endl;

//    Traits t;
//    typename K::Squared_distance_d csd = t.squared_distance_d_object();
//    Bare_point sQp = construct_bp(s_on_Q);
//    for(int i=0; i<(d+1); ++i)
//    {
//      std::size_t id = fs[i];
//      const Weighted_point& wi = wpoints[id];
//      FT di = csd(sQp, wi.point()) - wi.weight();
//      std::cout << "dist² sQp " << id << " : "  << di << std::endl;
//    }

    Vectord ret = Vectord::Zero();
    for(int i=0; i<d; ++i)
    {
      const Weighted_point& wp0 = wpoints[fs[0]];
      const Weighted_point& wi = wpoints[fs[i+1]];
      for(int j=0; j<D; ++j)
      {
        ret(i) += 2*(wi.point()[j]-wp0.point()[j])*s_on_Q(j);
        ret(i) += wp0.point()[j]*wp0.point()[j] - wi.point()[j]*wi.point()[j];
      }
      ret(i) += wi.weight() - wp0.weight();
//      std::cout << "ret: " << i << " " << ret(i) << std::endl;
    }

    return ret;
  }

  H(const Simplex& fs_, const std::vector<Weighted_point>& wpoints_)
    : fs(fs_), wpoints(wpoints_)
  { }
};

struct J // Jacobian
{
  const Simplex& fs;
  const std::vector<Weighted_point>& wpoints;

  // the matrix is the jacobian of the system :
  //  d(p,p0) - w0 = d(p,pi) - wi
  //  p \in Q
  // that is : J(i,j) = \frac{\partial \sum_k 2*(p_k-pi_k)x_j}{\partial p_j}
  // where the x_j are x x² (if d=1), x y x² xy y² (if d=2), etc.

  Matrixd operator()(const Vectord& s) const
  {
//std::cout << "computing J for: " << s.transpose() << std::endl;

    // everything below is complicated and ugly...
    // surely something smarter exists (but this is correct)
    // maybe using CGAL::combi_enum ? todo

    Matrixd m;
    Weighted_point wp0 = wpoints[fs[0]];
//std::cout << "wp0: " << fs.vertices[0] << " || " << wp0.point();

    //the easy coefficients : those in the part x y z of \tilde{P}
    for(int i=0; i<d; ++i)
    {
      Weighted_point wpi = wpoints[fs[i+1]];
//std::cout << "wpi: " << fs.vertices[i+1] << " || " << wpi.point();
      for(int j=0; j<d; ++j) // x y z...
        m(i,j) = 2.*(wpi.point()[j] - wp0.point()[j]);
    }

//std::cout << "after ez, J: " << std::endl << m << std::endl;

    // and now for the hard part: x² xy xz y² yz z² etc...
    for(int l=0; l<d; ++l) // lines of m (the function)
    {
      const Weighted_point wpl = wpoints[fs[l+1]];
      for(int c=0; c<d; ++c) // columns of m (the derivative)
      {
        //ij are the indices of the Q matrix
        //max is the number of entries at line i (max=d-i)
        //count is the number of entries visited (lexico order)
        int i=0, j=0, max=d, count=0;
        while(max > 0)
        {
          assert(i<d && j<d);
          FT coeff = 0.;
          if(i == c) // the line Q(i,.) is x_i*x_j
          {
            if(j == c) // the coeff Q(i,j) is x_i*x_i
              coeff = 2*s(i);
            else
              coeff = s(j); // the coeff Q(i,j) is x_i*x_j
          }
          else if(j == c) // the coeff Q(i,j) is x_j*x_i
            coeff = s(i);

          assert(d+count < D);
          // "+d" since we've added the easy part already (the d first)
//std::cout << "lc: " << l << " " << c << " " << 2.*(wpl.point()[d+count]-wp0.point()[d+count])*coeff << std::endl;
          m(l, c) += 2.*(wpl.point()[d+count]-wp0.point()[d+count])*coeff;

          j++; count++; // switch from x_i*x_j to x_i*x_{j+1}
          if(j >= max) // reached the end of a line of Q
          {
            max--; i++; j=i; // switch from x_i*x_d to x_{i+1}*x_{i+1}
          }
        }
      }
    }
//std::cout << "new J: " << std::endl << m << std::endl;

#ifdef nothing //debug stuff with the formulas for d=1, d=2
    VectorD p,q,r;
    //---------------d=1
    FT x = s(0);
    p(0) = wpoints[fs.vertices[0]].point()[0];
    p(1) = wpoints[fs.vertices[0]].point()[1];
    q(0) = wpoints[fs.vertices[1]].point()[0];
    q(1) = wpoints[fs.vertices[1]].point()[1];
    m(0,0) = 2.*(p(0)-q(0) + 2*(p(1)-q(1))*x);

    //---------------d=2
    FT y = s(1);
    p(2) = wpoints[fs.vertices[0]].point()[2];
    p(3) = wpoints[fs.vertices[0]].point()[3];
    p(4) = wpoints[fs.vertices[0]].point()[4];

    q(2) = wpoints[fs.vertices[1]].point()[2];
    q(3) = wpoints[fs.vertices[1]].point()[3];
    q(4) = wpoints[fs.vertices[1]].point()[4];

    r(0) = wpoints[fs.vertices[2]].point()[0];
    r(1) = wpoints[fs.vertices[2]].point()[1];
    r(2) = wpoints[fs.vertices[2]].point()[2];
    r(3) = wpoints[fs.vertices[2]].point()[3];
    r(4) = wpoints[fs.vertices[2]].point()[4];

    std::cout << "p " << p.transpose() << std::endl;
    std::cout << "q" << q.transpose() << std::endl;
    std::cout << "r " << r.transpose() << std::endl;
    std::cout << "azsd: " << x << " " << y << std::endl;

    m(0,0) = -2.*(p(0)-q(0) + 2*(p(2)-q(2))*x + (p(3)-q(3))*y);
    m(0,1) = -2.*(p(1)-q(1) + (p(3)-q(3))*x + 2*(p(4)-q(4))*y);
    m(1,0) = -2.*(p(0)-r(0) + 2*(p(2)-r(2))*x + (p(3)-r(3))*y);
    m(1,1) = -2.*(p(1)-r(1) + (p(3)-r(3))*x + 2*(p(4)-r(4))*y);

    std::cout << "old J: " << std::endl << m << std::endl;
#endif
    return m;
  }

  J(const Simplex& fs_, const std::vector<Weighted_point>& wpoints_)
    : fs(fs_), wpoints(wpoints_)
  { }
};

void tinker_jacobian(Matrixd& m, const Vectord& h)
{
  //we add a diagonal matrix of the form delta*f_i(x) to J to avoid |J| ~ 0
  Matrixd m2 = Matrixd::Zero();
  for(int i=0; i<d; ++i)
  {
    m2(i,i) = m(i,i)/std::abs(m(i,i)); // sign
    m2(i,i) *= h(i);
  }

  m += m2;
}

void buse_shenanigans(const Simplex& fs,
                      const std::vector<Pointd>& points,
                      const std::vector<Weighted_point>& wpoints)
{
  std::cout << "let's do some shenanigans!" << std::endl;
  CGAL_assertion(d == 2);

  // We solve the system analytically. We recall the system :
  // a0 + a1 x + a2 y + a3 xy + a4 x^2 + a5 y^2 = 0
  // b0 + b1 x + b2 y + b3 xy + b4 x^2 + b5 y^2 = 0

  // and ai = (P1-P0)(i) (except a0 = sum (P1-P0)(i)*Ci
  // and bi = (P2-P0)(i) (except b0 = sum (P2-P0)(i)*Ci

  // first of all, we need the coefficients
  FT a0 = 68, a1 = -29, a2 = 5, a3 = -51, a4 = -26, a5 = 88;
  FT b0 = 97, b1 = -67, b2 = 58, b3 = 37, b4 = 29, b5 = 5;

  // then we can hardcode some matrices
  Eigen::Matrix<FT,6,6> A, B;

  A(0,0) = 0.; A(0,1) = 0.; A(0,2) = 1.; A(0,3) = 0.; A(0,4) = 0.; A(0,5) = 0.;
  A(1,0) = 0.; A(1,1) = 0.; A(1,2) = 0.; A(1,3) = 1.; A(1,4) = 0.; A(1,5) = 0.;
  A(2,0) = 0.; A(2,1) = 0.; A(2,2) = 0.; A(2,3) = 0.; A(2,4) = 1.; A(2,5) = 0.;
  A(3,0) = 0.; A(3,1) = 0.; A(3,2) = 0.; A(3,3) = 0.; A(3,4) = 0.; A(3,5) = 1.;
  A(4,0) = -a2*b5+a5*b2; A(4,1) = -a0*b5+a5*b0; A(4,2) = -a3*b5+a5*b3; A(4,3) = -a1*b5+a5*b1; A(4,4) = 0.; A(4,5) = -a4*b5+a5*b4;
  A(5,0) = -a0*b5+a5*b0; A(5,1) = -a0*b2+a2*b0; A(5,2) = -a1*b5+a5*b1; A(5,3) = -a0*b3-a1*b2+a2*b1+a3*b0; A(5,4) = -a4*b5+a5*b4; A(5,5) = -a1*b3+a2*b4+a3*b1-a4*b2;

  B(0,0) = 1.; B(0,1) = 0.; B(0,2) = 0.; B(0,3) = 0.; B(0,4) = 0.; B(0,5) = 0.;
  B(1,0) = 0.; B(1,1) = 1.; B(1,2) = 0.; B(1,3) = 0.; B(1,4) = 0.; B(1,5) = 0.;
  B(2,0) = 0.; B(2,1) = 0.; B(2,2) = 1.; B(2,3) = 0.; B(2,4) = 0.; B(2,5) = 0.;
  B(3,0) = 0.; B(3,1) = 0.; B(3,2) = 0.; B(3,3) = 1.; B(3,4) = 0.; B(3,5) = 0.;
  B(4,0) = 0.; B(4,1) = 0.; B(4,2) = 0.; B(4,3) = 0.; B(4,4) = 0.; B(4,5) = 0.;
  B(5,0) = 0.; B(5,1) = 0.; B(5,2) = 0.; B(5,3) = 0.; B(5,4) = 0.; B(5,5) = -a3*b4+a4*b3;

  // we take the generalized eigen decomposition
  Eigen::GeneralizedEigenSolver<Eigen::Matrix<FT,6,6> > ges;
  ges.compute(A, B, true /*compute vectors*/);

  // and we look at the solutions
  std::cout << "The (complex) numerators of the generalized eigenvalues are: " << std::endl
            << ges.alphas() << std::endl;
  std::cout << "The (real) denominators of the generalized eigenvalues are: " << std::endl
            << ges.betas() << std::endl;
  std::cout << "The (complex) generalized eigenvalues are (alphas./beta): " << std::endl
            << ges.eigenvalues() << std::endl;
//  std::cout << "The (complex) generalized eigenvectors are : " << std::endl
//            << ges.eigenvectors() << std::endl; // not implemented yet...
}

bool newton(const Simplex& fs,
            const std::vector<Pointd>& points,
            const std::vector<Weighted_point>& wpoints,
            Bare_point& sol)
{
  std::cout << "NEWTON START ------" << std::endl;

  if(!is_support_intersecting(fs, points, wpoints, sol)) // get the approximation
    return false;

  Vectord vsol;
  for(int i=0; i<d; ++i)
    vsol(i) = sol[i];  // ugly 'projection' on the paraboloid
                       // of the intersection with the tangent plane

  FT err = 1e-5;
  H h(fs, wpoints); // compute the error
  J j(fs, wpoints); // compute the jacobian

  int max = 100, count = 0;
  Vectord h_err = h(vsol);
  Matrixd jac = j(vsol);
  tinker_jacobian(jac, h_err);

//std::cout << "init error : " << h_err.norm() << std::endl;
//std::cout << "init jac: " << jac << " det: " << jac.determinant() << std::endl;
  while(h_err.norm() > err && ++count < max)
  {
//std::cout << "removing: " << (jac.inverse()*h_err).transpose() << std::endl;
    vsol = vsol - jac.inverse()*h(vsol);
//std::cout << "new point : " << vsol.transpose() << std::endl;

    jac = j(vsol);
    h_err = h(vsol);
    tinker_jacobian(jac, h_err);
//std::cout << "new norm : " << h(vsol).norm() << std::endl;
//std::cout << "new jacobian determinant : " << jac.determinant() << std::endl;
  }

  if(count == max)
  {
    std::cout << "didn't converge in " << max << " iterations... ";
    std::cout << "h: " << h_err.transpose() << std::endl;
    return false;
  }

  VectorD vsolD = to_Q(vsol);
  sol = construct_bp(vsolD);
  std::cout << "We converged in " << count << " iterations ! ";
  std::cout << "final solution : " << vsolD.transpose() << std::endl;
  std::cout << "final jac/h: " << jac.determinant() << " " << h_err.norm() << std::endl;

  buse_shenanigans(fs, points, wpoints);
  CGAL_assertion(false);

  return true;
}

bool is_intersection_point_in_power_cell(const Bare_point& sol,
                                         const Simplex& fs,
                                         const std::vector<Weighted_point>& wpoints,
                                         const RTriangulation& rt)
{
  std::cout << "Check if the point is in the Pcell of the simplex" << std::endl;

  typename K::Squared_distance_d csd =
      rt.geom_traits().squared_distance_d_object();

  const Weighted_point& wp0 = wpoints[fs[0]];
  FT sq_d = csd(sol, wp0.point()) - wp0.weight();
  for(int i=0; i<d; i++)
  {
    std::size_t id = fs[i+1];
    const Weighted_point& wp = wpoints[id];
    FT sq_dd = csd(sol, wp.point()) - wp.weight();
    if(std::abs(sq_d-sq_dd)/sq_d > 1e-5) // sq_d might be null, I guess... fixme (and the other /sq_d too)
    {
      std::cout.precision(20);
      std::cout << "WARNING: dds don't agree (sol isn't in the dual) " << sq_d << " " << sq_dd << std::endl;
    }
  }

  Finite_vertex_const_iterator fvit = rt.finite_vertices_begin();
  for(; fvit!=rt.finite_vertices_end(); ++fvit)
  {
    std::size_t n = fvit->data();
    if(fs.has(n))
      continue;

    FT sq_dv = csd(sol, fvit->point().point()) - fvit->point().weight();
    if(std::abs(sq_dv-sq_d)/sq_d>1e-5 && sq_dv < sq_d) // not in the power cell of Pi
    {
      std::cout << "closer : " << n << " " << sq_dv << " ref: " << sq_d << std::endl;
      return false;
    }
  }
  return true;
}

bool is_face_restricted(const Simplex& fs,
                        const std::vector<Pointd>& points,
                        const std::vector<Weighted_point>& wpoints,
                        const RTriangulation& rt)
{
  Bare_point sol;
  //  if(is_support_intersecting(tr, points, wpoints, sol))
  if(newton(fs, points, wpoints, sol))
    return is_intersection_point_in_power_cell(sol, fs, wpoints, rt);
  return false;
}

void get_faces_in_cell(Finite_full_cell_iterator ch,
                       std::vector<Simplex>& faces,
                       const RTriangulation& rt)
{
  std::cout << "get faces in cell" << std::endl;
  std::cout << "dim rt: " << rt.current_dimension() << std::endl;
  std::cout << "dim " <<  ch->maximal_dimension() << std::endl;

  const int max_d = ch->maximal_dimension();
  int actual = 0;
  while(actual <= max_d &&
        ch->vertex(actual) != Vertex_handle())
  {
    std::cout << actual << " ";
    for(int j=0; j<D; ++j)
      std::cout << ch->vertex(actual)->point().point()[j] << " ";
    std::cout << ch->vertex(actual)->data() << std::endl;
    actual++;
  }

  std::cout << "GALACTICA ACTUAL : " << actual << std::endl;
  if(actual < d+1)
    return;

  CGAL::Combination_enumerator<int> combi(d+1, 0, actual);

  for(; !combi.finished(); combi++)
  {
    std::cout << "combi ";
    std::vector<int> ids;
    for(int j=0; j<(d+1); ++j)
    {
      std::cout << ch->vertex(combi[j])->data() << " ";
      ids.push_back(ch->vertex(combi[j])->data());
    }
    Simplex fs(ids.begin(), ids.end());
    std::cout << std::endl;
    faces.push_back(fs);
  }
}

void draw_stuff(const std::map<Simplex, int>& simplices_to_draw,
                const std::vector<Pointd>& points)
{
  int inconsistencies_count = 0;
  int tri_per_simplex = (d+1)*d*(d-1)/6.; // 3 from d+1
  int tn = tri_per_simplex * simplices_to_draw.size();

  // draw stuff
#if 1
  std::ofstream out("aniso_regular.mesh");
  out << "MeshVersionFormatted 1" << std::endl;
  out << "Dimension " << d << std::endl;
  out << "Vertices" << std::endl << points.size() << std::endl;
  for(std::size_t i=0; i<points.size(); ++i)
  {
    for(int j=0; j<d; ++j)
      out << points[i](j) << " ";
    out << i+1 << std::endl;
  }

  out << "Triangles" << std::endl << tn << std::endl;
  std::map<Simplex, int>::const_iterator it = simplices_to_draw.begin();
  std::map<Simplex, int>::const_iterator itend = simplices_to_draw.end();
  for(; it!=itend; ++it)
  {
    const Simplex& s = it->first;
    if(it->second != d+1)
      inconsistencies_count++;

    CGAL::Combination_enumerator<int> combi(3, 0, d+1);
    for(; !combi.finished(); combi++)
      out << s[combi[0]]+1 << " " << s[combi[1]]+1 << " " << s[combi[2]]+1 << " 1" << std::endl;
  }

  out << "End" << std::endl;
#else
  std::ofstream out("aniso_regular.off");
  out << "OFF" << std::endl;
  out << points.size() << " " << tn << " 0" << std::endl;
  for(std::size_t i=0; i<points.size(); ++i)
  {
    for(int j=0; j<d; ++j)
      out << points[i](j) << " ";
    for(int j=0; j<(3-d); ++j) // ugly trick: fill with 0
      out << "0 ";
    out << std::endl;
  }

  std::map<Simplex, int>::const_iterator it = simplices_to_draw.begin();
  std::map<Simplex, int>::const_iterator itend = simplices_to_draw.end();
  for(; it!=itend; ++it)
  {
    const Simplex& s = it->first;

    std::cout << "itsecond: " << it->second << " (" << it->first;
    if(it->second != d+1)
      inconsistencies_count++;

    CGAL::Combination_enumerator<int> combi(3, 0, d+1);
    for(; !combi.finished(); combi++)
      out << "3 " << s[combi[0]] << " " << s[combi[1]] << " " << s[combi[2]] << std::endl;
  }
#endif

  std::cout << inconsistencies_count << " inconsistencies" << std::endl;
}

void compute_restricted_facets(const RTriangulation& rt,
                               const std::vector<Pointd>& points,
                               const std::vector<Weighted_point>& wpoints)
{
  std::cout << "computing restricting facets" << std::endl;

  std::map<Simplex, int> simplices_to_draw;
  std::vector<Simplex> faces;

  //loop over vertices and loop over incident full cells of said vertices
  //to count faces d+1 times (if they are consistent)

  // note that this is always the case since we're in brute force mode...

  Finite_vertex_const_iterator vhit = rt.finite_vertices_begin();
  Finite_vertex_const_iterator vhitend = rt.finite_vertices_end();
  for(; vhit!=vhitend; ++vhit)
  {
    std::cout << "Adding the faces of vertex: " << vhit->data() << std::endl;

    // if the metric field is uniform rt.current_dim() = d
    // and we can't use incident_faces()
    if(rt.current_dimension() == d)
    {
      Full_cell_handle_vector cells;
      std::back_insert_iterator<Full_cell_handle_vector>
          incident_cells_insertor(cells);
      rt.incident_full_cells(vhit.base(), incident_cells_insertor);

      std::cout << "degen rt w cellsn " << cells.size() << std::endl;
      Full_cell_handle_iterator cit = cells.begin();
      Full_cell_handle_iterator cend = cells.end();
      for(; cit!=cend; ++cit)
      {
        if(rt.is_infinite(*cit))
          continue;

        std::vector<int> ids;
        for(int j=0; j<(d+1); ++j)
          ids.push_back((*cit)->vertex(j)->data());
        Simplex fs(ids.begin(), ids.end());
        std::cout << "add: " << fs << std::endl;
        faces.push_back(fs);
      }
    }
    else // get the faces of dimension d
    {
      std::vector<Face> ffaces;
      std::back_insert_iterator<std::vector<Face> > out(ffaces);
      rt.incident_faces(vhit.base(), d, out);

      typename std::vector<Face>::iterator fit = ffaces.begin();
      typename std::vector<Face>::iterator fend = ffaces.end();
      for(; fit<fend; ++fit)
      {
        if(rt.is_infinite(*fit))
          continue;

        std::vector<int> ids;
        for(int j=0; j<(d+1); ++j)
          ids.push_back(fit->vertex(j)->data());
          Simplex fs(ids.begin(), ids.end());
        faces.push_back(fs);
      }
    }
  }

  std::cout << faces.size() << " incident faces to consider" << std::endl;
  typename std::vector<Simplex>::iterator it = faces.begin();
  typename std::vector<Simplex>::iterator end = faces.end();
  for(; it!=end; ++it)
  {
    const Simplex& s = *it;

    //avoiding the degenerate case of a d-face in R^D being flat in R^d
    if(is_simplex_degenerated_in_Rd(s, points))
      continue;

    if(is_face_restricted(s, points, wpoints, rt))
      add_to_map(s, simplices_to_draw);
    else
    {
      std::cout << "non degenerate, non restricted: " ;
      std::cout << s[0] << " " << s[1] << " "  << s[2] << std::endl;
    }
  }

  std::cout << "total faces (restricted) " << faces.size() << " (" << simplices_to_draw.size() << ")" << std::endl;

  std::cout << "drawing stuff: " << simplices_to_draw.size() << std::endl;
  draw_stuff(simplices_to_draw, points);
}

void refine_biggest_face(RTriangulation& rt,
                         std::vector<Pointd>& points,
                         std::vector<Weighted_point>& wpoints)
{
  std::cout << "refine biggest" << std::endl;

  std::vector<Simplex> faces;

  //find biggest face
  Finite_full_cell_iterator ffcit = rt.finite_full_cells_begin();
  Finite_full_cell_iterator ffcend = rt.finite_full_cells_end();
  for(; ffcit!=ffcend; ++ffcit)
  {
    assert(ffcit->is_valid());
    get_faces_in_cell(ffcit, faces, rt);
  }

  std::cout << faces.size() << " faces" << std::endl;

  Bare_point sol;
  bool found = false;
  while(!found && !faces.empty())
  {
    int id = -1;
    FT max = -1e30;
    for(std::size_t i=0; i<faces.size(); ++i)
    {
      FT v = volume(faces[i], points);
      std::cout << "new v: " << v << std::endl;
      if(v > max)
      {
        id = i;
        max = v;
      }
    }

    std::cout << "max: " << id << " max: " << max << " || ";
    std::cout << faces[id][0] << " ";
    std::cout << faces[id][1] << " ";
    std::cout << faces[id][2] << std::endl;

    //refine biggest face
    //  if(is_support_intersecting(tr, points, wpoints, sol))
    if(newton(faces[id], points, wpoints, sol))
    {
      if(is_intersection_point_in_power_cell(sol, faces[id], wpoints, rt))
      {

        if(sol[0] < 1 || sol[0] > 2 || sol[1] < 1 || sol[1] > 2)
        {
          std::cout << "outside" << std::endl;
          faces.erase(faces.begin()+id);
          continue;
        }

        std::cout << "all good: point is computed" << std::endl;
        found = true;
      }
      else
      {
        faces.erase(faces.begin()+id);
        std::cout << "point not in power cell! " << faces.size() << " left" << std::endl;
      }
    }
    else
      std::cout << "no newton..." << std::endl;
  }

  if(faces.empty())
  {
    std::cout << "nothing to refine!" << std::endl;
    return;
  }

  Pointd p = get_pointd_from_RD(sol);
  Weighted_point wp = compute_wpoint(rt, p);

  points.push_back(p);
  wpoints.push_back(wp);
  RTriangulation::Vertex_handle v = rt.insert(wp);
  v->data() = rt.number_of_vertices()-1;

}

int main(int, char **)
{
  std::freopen("log.txt", "w", stdout); // redirect std::cout to log.txt

  const unsigned int N = 5;
  CGAL::Timer cost;

  std::vector<Pointd> points;
  std::vector<RTriangulation::Weighted_point> wpoints;

  RTriangulation rt(D);
  rt.infinite_vertex()->data() = -10;
  generate_wpoints(rt, points, wpoints, N);
  CGAL_assertion(rt.empty());

  // insert the wpoints in the RTriangulation
  cost.reset();
  cost.start();
  std::cout << "  Regular RTriangulation of " << N << " points in dim " << D << std::endl;

  for(std::size_t i=0; i<wpoints.size(); ++i)
  {
    std::cout << "i:" << i << std::endl;
    std::cout << "insert: ";
    for(int j=0; j<D; ++j)
      std::cout << wpoints[i].point()[j] << " ";
    std::cout << " || " << wpoints[i].weight() << std::endl;
    std::cout << "rt n: " <<  rt.number_of_vertices() << std::endl;
    RTriangulation::Vertex_handle v = rt.insert(wpoints[i]);
    v->data() = i;
  }

  std::cout << " done in " << cost.time() << " seconds." << std::endl;
  std::cout << "finite cells: " << rt.number_of_finite_full_cells() << std::endl;
  std::cout << "dim rt: " << rt.current_dimension() << std::endl;

/*
  refine_biggest_face(rt, points, wpoints);
  std::cout << "refine biggest done in " << cost.time() << " seconds." << std::endl;
  std::cout << "finite cells: " << rt.number_of_finite_full_cells() << std::endl;
  std::cout << "dim rt: " << rt.current_dimension() << std::endl;
*/

  compute_restricted_facets(rt, points, wpoints);

  return 0;
}

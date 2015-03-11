#include <CGAL/Starset.h>
#include <CGAL/Stretched_Delaunay_2.h>

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include <fstream>
#include <iostream>
#include <vector>

#define verb
#define debug

using namespace CGAL::Anisotropic_mesh_2;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef typename K::FT                                       FT;
typedef Stretched_Delaunay_2<K>                              Star;
typedef typename Star::Metric                                Metric;
typedef typename Star::Traits                                Traits;

typedef typename Eigen::Vector2d                             Vectord;
typedef Vectord                                              Point;
typedef Vectord                                              IPoint;

const int d = 2; // immersion space dimension
const int n = 2; // dimension of the manifold
const int D = 2*d+3; // embedding space dimension
typedef Eigen::Matrix<FT, d, 1>                             Vectord;
typedef Eigen::Matrix<FT, d, d>                             Matrixdd;
typedef Eigen::Matrix<FT, D, 1>                             VectorD;

const std::size_t N = 100;
const FT rho = 0.2; // covering criterion
const FT omega = 100; // weight
const FT delta = 1e-3; // approximation for the tangent vector computation (algo 2)

VectorD random_nu, random_mu; // orthogonal vectors to compute Psi_ij

class Position // used to create a memory of Psi_ij at t
{
public:
  const IPoint t; // the coordinates of the point
  const int i,j; // for Psi_ij

  Position(const IPoint& t_, const int i_, const int j_)
    : t(t_), i(i_), j(j_)
  { }

  bool operator<(const Position& p) const
  {
    if(this->i != p.i) return this->i < p.i;
    if(this->j != p.j) return this->j < p.j;
    for(int k=0; k<d; ++k)
      if(this->t(k) != p.t(k))
        return this->t(k) < p.t(k);
    return false;
  }

  void print() const
  {
    std::cout << "t: " << t.transpose() << " ij: " << i << " " << j << std::endl;
  }
};

void print_memory(const std::map<Position, VectorD>& m)
{
  std::cout << "PRINT MEM: " << m.size() << std::endl;
  std::map<Position, VectorD>::const_iterator it = m.begin();
  std::map<Position, VectorD>::const_iterator iend = m.end();
  for(; it!=iend; ++it)
  {
    it->first.print();
    std::cout << it->second.transpose() << std::endl;
  }
}

void generate_points(std::vector<Point>& points)
{
  double r_x, r_y;
#if 0
  for(std::size_t i=0; i<N; ++i)
  {
    r_x = ((double) rand() / (RAND_MAX));
    r_y = ((double) rand() / (RAND_MAX));
    Vectord v; v(0) = r_x; v(1) = r_y;
    points.push_back(v);
  }
#else
  for(int i=0; i<10; ++i)
  {
    for(int j=0; j<10; ++j)
    {
      r_x = i*0.1;
      r_y = j*0.1;
      Point v;
      v(0) = r_x; v(1) = r_y;
      points.push_back(v);
    }
  }
//  std::random_shuffle(points.begin(), points.end());
#endif
  assert(points.size() == N);
  std::cout << "generated initial points" << std::endl;
}

Matrixdd immersion_matrix(const Point& p) // needs to be contractive wrt the metric
{
  Matrixdd phi;
  phi(0,0) = 0.1; phi(0,1) = 0.;
  phi(1,0) = 0.; phi(1,1) = 0.1;

  return phi;
}

Matrixdd transformation(const Point& p) // F_M (sqrt of the metric)
{
  Matrixdd m;
  m(0,0) = 1.; m(0,1) = 0.;
  m(1,0) = 0.; m(1,1) = 1.;

  return m;
}

Eigen::Matrix2d metric(const Point& p) // full metric
{
  Eigen::Matrix2d m = transformation(p);
  return m.transpose()*m;
}

FT met_dist(const Vectord& v1, const Vectord& v2) // distance in the metric
{
  return (std::max)((transformation(v1)*(v2-v1)).norm(), (transformation(v2)*(v2-v1)).norm());
}

FT dist(const std::vector<FT>& v1, const std::vector<FT>& v2) // euclidean distance
{
  assert(v1.size() == v2.size());
  FT n = 0.;
  for(std::size_t i=0; i<v1.size(); ++i)
    n += (v2[i] - v1[i]) * (v2[i] - v1[i]);

  return std::sqrt(n);
}

void compute_eigen_matrices(const Matrixdd& m,
                            Matrixdd& sigma,
                            Matrixdd& u)
{
  Eigen::EigenSolver<Matrixdd> es(m, true);
  const Eigen::EigenSolver<Matrixdd>::EigenvectorsType& vecs = es.eigenvectors();
  const Eigen::EigenSolver<Matrixdd>::EigenvalueType& vals = es.eigenvalues();

  sigma(0,0) = std::real(vals[0]);
  sigma(1,1) = std::real(vals[1]);

  u(0,0) = std::real(vecs(0,0)); u(0,1) = std::real(vecs(0,1));
  u(1,0) = std::real(vecs(1,0)); u(1,1) = std::real(vecs(1,1));
}

Vectord immerse(const Vectord t)
{
  return immersion_matrix(t)*t;
}

bool find_acceptable_bucket(const int v_id,
                            std::set<int>& bucket,
                            const std::vector<Point>& points)
{
  //find a bucket where the distance to all the other points is > 2*rho
  bool acceptable =  true;
  Point v = points[v_id];
  for(std::set<int>::iterator it=bucket.begin(); it!=bucket.end(); ++it)
  {
    Point w = points[*it];
    FT d = (v-w).norm();
//    std::cout << "ijd: " << v_id << " " << *it << " " << d << std::endl;
    if(d < 2. * rho)
      acceptable = false;
  }
  if(acceptable)
    bucket.insert(v_id);
  return acceptable;
}

int build_partition(std::vector<std::set<int> >& buckets,
                    const std::vector<Vectord>& points)
{
  for(std::size_t i=0; i<points.size(); ++i)
  {
    std::size_t j=0;
    bool found = false;
    while(!found && j<buckets.size())
    {
      if(find_acceptable_bucket(i, buckets[j++], points))
        found = true;
    }

    if(!found)
    {
      std::set<int> new_bucket;
      new_bucket.insert(i);
      buckets.push_back(new_bucket);
    }
  }

  return buckets.size();
}

FT bump_function(const Point x, // the center of the bump
                 const IPoint t) // the query
{
  IPoint ix = immerse(x); // very unefficient to compute them here fixme
  FT norm = (t-ix).norm();

  if(norm >= rho)
    return 0.;

  FT frac = norm*norm / (rho*rho);
  FT ret = std::exp(-1./(1.-frac));
  return ret;
}

FT full_bump_function(const Point x,
                      const IPoint t,
                      const std::vector<Vectord>& points)
{
  // full bump = bump(x, t) / sum_i bump(i, t)

  FT sum = 0.;

  for(std::size_t i=0; i<points.size(); ++i)
  {
    Point v = points[i];
    FT v_bump = bump_function(v, t);
    sum += v_bump;
  }

  return bump_function(x, t) / sum;
}

Matrixdd compute_correction_matrix(const Point p)
{
  Matrixdd m = metric(p);
  Matrixdd phi = immersion_matrix(p);
#if 0 // using the eigen decomp of the Phi * M
  // C = (Sigma^-2 - Id)^1/2 U^T
  Matrixdd sigma = Matrixdd::Zero();
  Matrixdd u = Matrixdd::Zero();
  compute_eigen_matrices(phi*m, sigma, u);

  if(sigma(0,0) > 1. || sigma(0,0) <= 0 || sigma(1,1) > 1. || sigma(1,1) <= 0 )
  {
    std::cerr << "es: " << sigma(0,0) << " " << sigma(1,1) << std::endl;
    assert(0);
  }

  FT e0 = std::sqrt(1./(sigma(0,0)*sigma(0,0))-1.); // (Sigma^-2 - Id)^1/2
  FT e1 = std::sqrt(1./(sigma(1,1)*sigma(1,1))-1.);

  Matrixdd c = Eigen::Matrix2d::Zero();
  c(0,0) = e0; c(1,1) = e1;

  c = c*u.transpose();
  return c;
#else
  // C = ((Sigma^-1*F_M*Sigma^-1) - Id)^1/2 U^T
  Matrixdd sigma = Matrixdd::Zero();
  Matrixdd u = Matrixdd::Zero();
  compute_eigen_matrices(phi, sigma, u);

  Matrixdd sigma_inv = Matrixdd::Zero();
  for(int i=0; i<d; ++i)
    sigma_inv(i, i) = 1./sigma(i,i);

  Matrixdd c = (sigma_inv * m * sigma_inv - Matrixdd::Identity());
  if(c == Matrixdd::Zero())
    return c;

  c = c.pow(0.5);
  c *= u.transpose();
  return c;
#endif
}

std::vector<FT> embed_1(const Vectord p,
                        const std::vector<std::set<int > >& buckets,
                        const std::vector<Point>& points)
{
  // First version of the embedding, brutal & high dimensional embedding space
  IPoint t = immerse(p);

  std::size_t K = buckets.size();
  std::vector<FT> ret(d+2*n*K, 0.);

  for(int i=0; i<d; ++i)
    ret[i] = t(i);

  std::size_t pos = d;
  for(std::size_t i=0; i<K; ++i)
  {
    const std::set<int>& bucket = buckets[i];
    for(std::set<int>::iterator bit=bucket.begin(); bit!=bucket.end(); ++bit)
    {
      Vectord x = points[*bit];
      Matrixdd m = compute_correction_matrix(x);

      for(int j=0; j<n; ++j)
      {
        FT L = full_bump_function(x, t, points);
        FT coeff = std::sqrt(L)/omega;

//        std::cout << "lcoeff: " << L << " " << coeff << std::endl;

        ret[pos+j] += coeff * std::sin(omega * (m*t)(j));
        ret[pos+n+j] += coeff * std::cos(omega * (m*t)(j));
      }
    }
    pos += 2*n;
  }
  assert(pos == d+2*n*K);
  return ret;
}

void tests_1(const std::vector<Point>& points,
             const std::vector<std::set<int > >& buckets)
{
  // embed stuff
  Vectord p0 = points[804];
  Vectord p1 = points[4];
  std::vector<FT> vec0 = embed_1(p0, buckets, points);
  std::vector<FT> vec1 = embed_1(p1, buckets, points);

  std::cout << "embeds to the dimension: " << vec0.size() << " : " << std::endl;
  for(std::size_t i=0; i<vec0.size(); ++i)
    std::cout << vec0[i] << " ";
  std::cout << std::endl;

  std::cout << "embeds to the dimension: " << vec1.size() << " : " << std::endl;
  for(std::size_t i=0; i<vec1.size(); ++i)
    std::cout << vec1[i] << " ";
  std::cout << std::endl;

  std::cout << "point0: " << p0.transpose() << std::endl;
  std::cout << "point1: " << p1.transpose() << std::endl;
  std::cout << "dist: " << (p0-p1).norm() << " in the metric: " << met_dist(p0, p1) << std::endl;
  std::cout << "dist in R^D: " << dist(vec0, vec1) << std::endl;
}

void algo_1(const std::vector<Point>& points)
{
  std::vector<std::set<int> > buckets;
  build_partition(buckets, points);

#ifdef verb
  std::cout << points.size() << " points" << std::endl;
  std::cout << buckets.size() << " buckets" << std::endl;

  std::size_t total = 0;
  for(std::size_t i=0; i<buckets.size(); ++i)
  {
    total += buckets[i].size();
    std::cout << "bucket : " << i << " has " << buckets[i].size() << " pts" << std::endl;
  }
  std::cout << "total : "  << total << std::endl;
  assert(total == points.size());

  tests_1(points, buckets);
#endif
}

VectorD compute_step_ij_embed_2(const IPoint& t,
                                const std::vector<Point>& points,
                                std::map<Position, VectorD>& memory,
                                int Ci, int Cj);

void compute_orthogonal_vectors(const IPoint& t, VectorD& nu, VectorD& mu,
                                const std::vector<Point>& points,
                                std::map<Position, VectorD>& memory,
                                int Ci, int Cj)
{
  std::cout << "compute_orthogonal_vectors +++++++++++++++";
  std::cout <<  t.transpose() << " @@ " << Ci << " " << Cj << std::endl;

  if(Cj == 0) { Cj=n-1; Ci--; }
  else Cj--;

  std::vector<VectorD> tangent_vectors(d);
  nu = random_nu, mu = random_mu;

  VectorD v2 = compute_step_ij_embed_2(t, points, memory, Ci, Cj);
  for(int k=0; k<d; ++k)
  {
    Vectord ek = Vectord::Zero(); ek[k] = 1.;
    VectorD v1 = compute_step_ij_embed_2(t+delta*ek, points, memory, Ci, Cj);
    tangent_vectors[k] = (v1 - v2) / delta;
//    std::cout << "v1 @ k: " << k << " " << v1.transpose() << std::endl;
//    std::cout << "Tk @ k: " << k << " " << tangent_vectors[k].transpose() << std::endl;
  }

  tangent_vectors.push_back(nu);
  tangent_vectors.push_back(mu);

//  std::cout << "computed TKs, now orthogonormalizing" << std::endl;
  for(int i=0; i<d+2; ++i)
  {
    VectorD& tvi = tangent_vectors[i];
//    std::cout << i << " ||| " << tangent_vectors[i].transpose() << std::endl;
    for(int j=0; j<i; ++j)
    {
      VectorD& tvj = tangent_vectors[j];
      FT alpha = tvi.dot(tvj) / tvj.dot(tvj);
      tvi -= alpha*tvj;
    }

    tvi /= tvi.norm();
  }

  nu = tangent_vectors[d];
  mu = tangent_vectors[d+1];

#ifdef debug
//  std::cout << "compute nu mu @ " << Ci << " " << Cj << " for " << t.transpose() << std::endl;
//  for(int i=0; i<d+2; ++i)
//    std::cout << i << " ||| " << tangent_vectors[i].transpose() << std::endl;

  for(int k=0; k<d+2; ++k) // all Tk + nu + mu
  {
    VectorD v1 = tangent_vectors[k];
    for(int kk=0; kk<d+2; ++kk)
    {
      VectorD v2 = tangent_vectors[kk];
//      std::cout << "v1v2 debug: " << k << " " << kk << " " << v1.dot(v2) << std::endl;
      if(k != kk)
        assert(std::abs(v1.dot(v2)) < 1e-8);
      else
        assert(std::abs(v1.norm()-1.) < 1e-8);
    }
  }
#endif
}

VectorD compute_step_ij_embed_2(const IPoint& t,
                                const std::vector<Point>& points,
                                std::map<Position, VectorD>& memory,
                                int Ci, int Cj)
{
  // compute Psi_ij(t)

  std::cout << "compute_step_ij_embed_2 ----------------- ";
  std::cout << t.transpose() << " @@ " << Ci << " " << Cj << std::endl;

  std::map<Position, VectorD>::iterator it = memory.find(Position(t, Ci, Cj));
  if(it != memory.end())
  {
//    print_memory(memory);
    std::cout << "already computed " << t.transpose() << " " << Ci << " " << Cj << std::endl;
    std::cout << "val: " << it->second.transpose() << std::endl;
    return it->second;
  }

  VectorD ret = VectorD::Zero();
  for(int i=0; i<d; ++i)
    ret[i] = t(i);

  if(Ci < 0)
  {
    assert(Cj == n-1);
    std::cout << "Add to memory: " << t.transpose() << " " << Ci << " " << Cj << std::endl;
    std::cout << "val: " << ret.transpose() << std::endl;
    memory[Position(t, Ci, Cj)] = ret;
    return ret;
  }

  assert(Ci <= int(points.size()-1) && Cj <= n-1);

  // the loop stops at (i,j) = (Ci, Cj+1) [this is possibly (i,j) = (Ci+1,0)]
  int i=0, j=0;
  while(i<Ci || (i==Ci && j<=Cj))
  {
    std::cout << "stepij: ~~~~~~~~~~~~~~~~~~~~ ";
    std::cout << t.transpose() << " " << i << " " << j << std::endl;

    it = memory.find(Position(t, i, j));
    if(it != memory.end())
    {
//      print_memory(memory);
      std::cout << "already computed (iw) ";
      std::cout << t.transpose() << " " << i << " " << j << std::endl;
      std::cout << "val: " << it->second.transpose() << std::endl;
      ret += it->second;
    }
    else
    {
      Point x = points[i]; // fixme (m is computed j too many times)
      Matrixdd m = compute_correction_matrix(x);

      VectorD nu, mu;
      compute_orthogonal_vectors(t, nu, mu, points, memory, i, j);

      FT L = full_bump_function(x, t, points);
      FT coeff = std::sqrt(L)/omega;
      ret += coeff*(std::sin(omega*(m*t)(j))*nu + std::cos(omega*(m*t)(j))*mu);

      std::cout << "Add to memory (iw) : " << t.transpose() << " " << i << " " << j << std::endl;
      std::cout << "val: " << ret.transpose() << std::endl;
      memory[Position(t, i, j)] = ret;
    }

//    std::cout << "ret is now: " << ret.transpose() << std::endl;

    j++;
    if(j == n) { i++; j=0; }
  }
  std::cout << "end of while with : " << Ci << " " << Cj << std::endl;

  return ret;
}

void initialize_random_vectors()
{
  random_nu = VectorD::Random(); //FIXME ("drawn from the surface of the unit sphere")
  random_mu = VectorD::Random();
  for(int i=0; i<d; ++i) // orthogonal to any vertex of the unit sphere of R^d
  {
    random_nu(i) = 0;
    random_mu(i) = 0;
  }

  std::cout << "initialized numu to : " << std::endl;
  std::cout << random_nu.transpose() << std::endl << random_mu.transpose() << std::endl;
}

VectorD embed_2(const Point& p,
                const std::vector<Point>& points)
{
  std::cout << "embed call for p : " << p.transpose() << std::endl;
  IPoint t = immerse(p);
  std::map<Position, VectorD> memory;
  return compute_step_ij_embed_2(t, points, memory, points.size()-1, n-1);
}

void tests_2(const std::vector<Point>& points)
{
  initialize_random_vectors();

  // embed stuff
  Vectord p0 = points[804];
  Vectord p1 = points[4];
  VectorD vec0 = embed_2(p0, points);
  VectorD vec1 = embed_2(p1, points);

  std::cout << "embeds to the dimension: " << vec0.size() << " : " << std::endl;
  for(std::size_t i=0; i<D; ++i)
    std::cout << vec0(i) << " ";
  std::cout << std::endl;

  std::cout << "embeds to the dimension: " << vec1.size() << " : " << std::endl;
  for(std::size_t i=0; i<D; ++i)
    std::cout << vec1(i) << " ";
  std::cout << std::endl;

  std::cout << "point0: " << p0.transpose() << std::endl;
  std::cout << "point1: " << p1.transpose() << std::endl;
  std::cout << "dist: " << (p0-p1).norm() << " in the metric: " << met_dist(p0, p1) << std::endl;
  std::cout << "dist in R^D: " << (vec0-vec1).norm() << std::endl;
}

int main(int, char**)
{
  std::cout.precision(10);
//  std::freopen("log.txt", "w", stdout);
  std::srand(0);

  std::cout << "constants: " << std::endl;
  std::cout << "d: " << d << " n:" << n << " D:" << D << " N: " << N << " ";
  std::cout << "delta: " << delta << " rho: " << rho << " omega: " << omega << std::endl;

  std::vector<Point> points;
  generate_points(points);

//  algo_1(points);
  tests_2(points);

  return 0;
}


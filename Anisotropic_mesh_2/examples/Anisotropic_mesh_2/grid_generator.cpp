#include <CGAL/Starset.h>
#include <CGAL/Stretched_Delaunay_2.h>

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <Eigen/Dense>

#include <fstream>
#include <iostream>
#include <vector>

using namespace CGAL::Anisotropic_mesh_2;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef typename K::FT                                       FT;
typedef Stretched_Delaunay_2<K>                              Star;
typedef typename Star::Point_2                               Point_2;
typedef typename Star::TPoint_2                              TPoint_2;
typedef typename Star::Metric                                Metric;
typedef typename Star::Traits                                Traits;

typedef typename Eigen::Matrix<double, 5, 1>                 Vector5d;
typedef typename Eigen::Matrix<double, 2, 1>                 Vector2d;

typedef std::set<std::size_t> Simplex;

// using a lot of global variables, it's ugly but passing them by function is tedious

std::size_t vertices_nv = 2000;

Point_2 center(0.,0.);
//#define R2
#ifdef R2
const FT grid_side = 4;
FT offset_x = -grid_side/2.; // offset is the bottom left point
FT offset_y = -grid_side/2.; // todo normalize this with aniso_mesh_2's rectangle whose offset is the center of the rectangle...
#else
const FT grid_side = 4;
FT offset_x = -grid_side/2.;//-2.75; // offset is the bottom left point
FT offset_y = -grid_side/2.;//0.9;
#endif

std::vector<Vector2d> R2seeds;
std::vector<FT> r2ws;
std::vector<Point_2> seeds;
std::vector<Vector5d> R5seeds;
std::vector<Metric> seeds_m;
std::vector<FT> ws;

std::set<Simplex> simplices;
std::vector<int> random_colors;

std::vector<std::vector<std::size_t> > witness_grid; // grid_size * 3
std::vector<std::vector<std::size_t> > first_position; // could be a vector of booleans...
std::vector<std::vector<std::size_t> > second_position;

Vector5d compute_hat(const Point_2& p,
                     const Metric& m)
{
  Eigen::Vector2d p_base;
  Vector5d p_hat;
  Eigen::Vector2d p_circle;
  Eigen::Vector3d p_wave;
  Eigen::Matrix2d mat = m.get_mat();
#if 1
  p_base[0] = p.x(); p_base[1] = p.y();
  p_circle = mat*p_base;
  p_wave[0] = -0.5*mat(0,0); p_wave[1] = -mat(0,1); p_wave[2] = -0.5*mat(1,1);
  p_hat[0] = p_circle[0]; p_hat[1] = p_circle[1];
  p_hat[2] = p_wave[0]; p_hat[3] = p_wave[1]; p_hat[4] = p_wave[2];
#elif 1
  p_base[0] = p.x(); p_base[1] = p.y();
  p_circle = mat*p_base;
  p_wave[0] = mat(0,0); p_wave[1] = mat(0,1); p_wave[2] = mat(1,1);
  p_hat[0] = p_circle[0]; p_hat[1] = p_circle[1];
  p_hat[2] = p_wave[0]; p_hat[3] = p_wave[1]; p_hat[4] = p_wave[2];
#else
  p_hat[0] = p.x(); p_hat[1] = p.y();
  p_hat[2] = mat(0,0); p_hat[3] = mat(0,1); p_hat[4] = mat(1,1);
#endif

  return p_hat;
}

template<typename Metric_field>
int build_seeds(const Metric_field& mf)
{
#if 1
  std::ifstream in("bambimboum.mesh");
//  std::ifstream in("/home/mrouxel/cgal/Anisotropic_mesh_TC/examples/Anisotropic_mesh_TC/build/aniso_TC.mesh");
  std::string word;
  std::size_t useless, nv, dim;
  FT r_x, r_y;

  in >> word >> useless; //MeshVersionFormatted i
  in >> word >> dim; //Dimension d
  in >> word >> nv;
  std::cout << "nv: " << nv << std::endl;
  assert(dim == 2);
  if(vertices_nv < nv)
    nv = vertices_nv;

  seeds.resize(nv);
  R5seeds.resize(nv);
  seeds_m.resize(nv);
  ws.resize(nv);

  for(std::size_t i=0; i<nv; ++i)
  {
    in >> r_x >> r_y >> useless;
#else
  for(std::size_t i=0; i<seeds.size(); ++i)
  {
    double r_x = offset_x + grid_side * ((double) rand() / (RAND_MAX));
    double r_y = offset_y + grid_side * ((double) rand() / (RAND_MAX));
#endif
    seeds[i] = Point_2(r_x, r_y);
    std::cout << "Seeds i: " << i << " " << seeds[i].x() << " " << seeds[i].y() << std::endl;

    seeds_m[i] = mf->compute_metric(seeds[i]);
    R5seeds[i] = compute_hat(seeds[i], seeds_m[i]);
    ws[i] = R5seeds[i].norm()*R5seeds[i].norm() - seeds[i].x()*R5seeds[i](0)
                                                  - seeds[i].y()*R5seeds[i](1);

//    std::cout << "check " << std::endl;
//    std::cout << seeds[i].x() << " " << seeds[i].y() << std::endl;
//    std::cout << seeds_m[i].get_mat() << std::endl;
//    std::cout << R5seeds[i](0) << " " << R5seeds[i](1) << std::endl;
//    std::cout << ws[i] << std::endl;
  }
  return nv;
}

int read_tangent_plane_seeds()
{
//  std::ifstream in("/home/mrouxel/cgal/Anisotropic_mesh_TC/examples/Anisotropic_mesh_TC/build/projected_points.txt");
  std::ifstream in("/home/mrouxel/cgal/Tangential_complex/test/Tangential_complex/build/projected_points_120.txt");
  std::size_t nv, dim;
  FT x, y, w;

  in >> dim >> nv;
  assert(dim == 2); //tmp
  if(vertices_nv < nv)
    nv = vertices_nv;
  std::cout << "nv: " << nv << std::endl;
  R2seeds.resize(nv);
  r2ws.resize(nv);

  for(std::size_t i=0; i<nv; ++i)
  {
    in >> x >> y >> w;
    R2seeds[i] = (Vector2d() << x, y ).finished();
    r2ws[i] = w;
  }
  return nv;
}

// comparator needs to be outside for c++03...
struct Comp
{
  const std::vector<FT>& m_sqds;
  bool operator()(std::size_t i, std::size_t j) const
  {
    return m_sqds[i] < m_sqds[j];
  }

  Comp(const std::vector<FT>& sqds_) : m_sqds(sqds_) { }
};

template<typename MF>
FT value_at_point(const Point_2& p, std::size_t p_id,
                  const MF* mf,
                  const Traits* traits)
{
  Metric m_p = mf->compute_metric(p);
  Vector2d pv = (Vector2d() << p.x(), p.y()).finished();
  typename Star::Traits::Compute_squared_distance_2 csd =
      traits->compute_squared_distance_2_object();

  // FIRST MODE : COMPUTE THE DISTANCE IN THE METRIC BETWEEN THE POINTS IN R^2
  TPoint_2 tp(m_p.transform(p));
  TPoint_2 tc(m_p.transform(center));
  //return csd(tp, tc);

  Vector5d p_hat(compute_hat(p, m_p));
  // SECOND MODE : COMPUTE THE EUCLIDEAN DISTANCE BETWEEN THE POINTS IN R^5
  //return (p_hat - m_hat).norm();

  // THIRD MODE : BUILD A VORONOI
  FT min = 1e30;
  unsigned int min_id = 0;
  std::vector<FT> sqds(vertices_nv, 1e30);
  std::vector<std::size_t> ids(vertices_nv);

  for(int l=0; l<1; ++l) // set l = 2 to draw borders
  {
    for(std::size_t k=0; k<vertices_nv; ++k)
    {
#ifndef R2
      // EUCLIDEAN DISTANCE IN R^5
//    double sq_d = (p_hat - R5seeds[k]).norm();

      // WEIGHTED DISTANCE IN R^5
//    double sq_d = std::pow((p_hat - R5seeds[k]).norm(),2.) - ws[k];

      // DU WANG : dx(p,x) < dx(q,x)
    double sq_d = csd(tp, m_p.transform(seeds[k]));

      // LABELLE SHEWCHUK : dp(p,x) < dq(q,x)
      TPoint_2 tp2 = (seeds_m[k]).transform(p);
      TPoint_2 ts = (seeds_m[k]).transform(seeds[k]);
//      double sq_d = csd(tp2, ts);
#else
      // WEIGHTED DISTANCE IN THE TANGENT PLANE (seeds are the projected p_hat)
      FT sq_d = std::pow((pv - R2seeds[k]).norm(), 2.) - r2ws[k];
#endif
      sqds[k] = sq_d;
      ids[k] = k;

      if(l==0 && sq_d<min) // first pass (compute the distance)
      {
        min_id = random_colors[k];
        min = sq_d;
      }
      else if(l==1 && ((random_colors[k]!=min_id && std::abs(sq_d-min)<0.1*min) ||
                       sq_d < 0.01*min)) // second pass
      {
        min_id = 1.01*vertices_nv; // drawing borders and center points
        break;
      }
    }
  }

  Comp c(sqds);
  std::sort(ids.begin(), ids.end(), c);

  //for the witness
  witness_grid.push_back(std::vector<std::size_t>(ids.begin(), ids.begin()+3));
  first_position[ids[0]].push_back(p_id);
  second_position[ids[1]].push_back(p_id);

  // check if we are in the dual of a triangle (somewhat, not exactly a witness)
  for(int i=0; i<vertices_nv-1; ++i)
    assert(sqds[ids[i]] < sqds[ids[i+1]]);

  FT alpha = 0.1;
  if(std::abs(sqds[ids[0]] - sqds[ids[1]]) < alpha*sqds[ids[0]] &&
     std::abs(sqds[ids[0]] - sqds[ids[2]]) < alpha*sqds[ids[0]]) // in the dual
  {
    std::set<std::size_t> simplex;
    simplex.insert(ids[0]); simplex.insert(ids[1]); simplex.insert(ids[2]);
    simplices.insert(simplex);
  }

  return min_id;
}

template<typename MF>
void full_grid(const MF* metric_field, const Traits* traits)
{
  double n = 1000.; // number of points per side
  double a = grid_side; // length of the side
  double step = a/n;
  double sq_n = n*n;
  double tri_n = 2*(n-1)*(n-1);
  std::size_t counter = 0;

  std::cout << "full grid with: " << sq_n << " vertices" << std::endl;

  std::ofstream out("grid.mesh");
  out << "MeshVersionFormatted 1" << std::endl;
  out << "Dimension 2" << std::endl;
  out << "Vertices" << std::endl;
  out << (int) sq_n << std::endl;

  std::ofstream out_bb("grid.bb");
  out_bb << "2 1 " << (int) sq_n << " 2" << std::endl;

  for(unsigned int i=0; i<n; ++i)
  {
    for(unsigned int j=0; j<n; ++j)
    {
      //filling from bot left to top right
      Point_2 p(offset_x+j*step, offset_y+i*step);
      out << p.x() << " " << p.y() << " " << ++counter << std::endl;

      FT ret = value_at_point(p, metric_field, traits);
      out_bb << ret << std::endl;
    }
  }
  out << "Triangles" << std::endl;
  out << (int) tri_n << std::endl;
  for(int i=1; i<n; ++i)
  {
    for(int j=0; j<(n-1); ++j)
    {
      int id1 = j*n + i;
      int id2 = j*n + (i+1);
      int id3 = (j+1)*n + i;
      int id4 = (j+1)*n + (i+1);
      out << id1 << " " << id2 << " " << id4 << " 1" << std::endl;
      out << id1 << " " << id3 << " " << id4 << " 1" << std::endl;
    }
  }
  out << "End" << std::endl;
}

struct Quad
{
  // p1---p0
  // |    |
  // p2---p3

  const std::size_t i0, i1, i2, i3;

  bool is_too_small(FT min_vol, const std::vector<Point_2>& points) const
  {
    FT dx = points[i0].x() - points[i1].x();
    FT dy = points[i1].y() - points[i2].y();
    return dx*dy < min_vol;
  }

  bool is_too_big(FT max_vol, const std::vector<Point_2>& points) const
  {
    FT dx = points[i0].x() - points[i1].x();
    FT dy = points[i1].y() - points[i2].y();
    return dx*dy > max_vol;
  }

  bool has_same_colors(const std::vector<FT>& values) const
  {
    CGAL_assertion(values.find(i0) != values.end() &&
                   values.find(i1) != values.end() &&
                   values.find(i2) != values.end() &&
                   values.find(i3) != values.end());
    FT v0 = values[i0];
    if(v0 != values[i1])
      return false;
    if(v0 != values[i2])
      return false;
    if(v0 != values[i3])
      return false;
    return true;
  }

  std::size_t number_of_colors(const std::vector<FT>& values) const
  {
    std::set<int> colors;
    colors.insert(values[i0]); colors.insert(values[i1]);
    colors.insert(values[i2]); colors.insert(values[i3]);
    return colors.size();
  }

  Quad(const std::size_t i0_, const std::size_t i1_,
       const std::size_t i2_, const std::size_t i3_)
    : i0(i0_), i1(i1_), i2(i2_), i3(i3_)
  { }
};

bool is_simplex_in_witness_complex(std::size_t l0,
                                   std::size_t l1,
                                   std::size_t l2)
{
  // we already know that there exists w for 012 so 0 and 01 are taken care of
  // we must find 1, 2, 12, and 02 in witness_grid

  // there might be a more subtle approach but we're going for clarity !
  bool found_1 = false, found_2 = false, found_12 = false, found_02 = false;

  //quick filters (need to have l1 and l2 in a 0th position somewhere, etc.)
  if(first_position[l1].empty())
    return false;
  if(first_position[l2].empty())
    return false;
  if((first_position[l0].empty() || second_position[l2].empty()) &&
     (first_position[l2].empty() || second_position[l0].empty()))
    return false;
  if((first_position[l1].empty() || second_position[l2].empty()) &&
     (first_position[l2].empty() || second_position[l1].empty()))
    return false;

  //proper checks
  for(std::size_t w=0; w<witness_grid.size(); ++w)
  {
    if(witness_grid[w][0] == l0 && witness_grid[w][1] == l2)
      found_02 = true;
    else if(witness_grid[w][0] == l1)
    {
      found_1 = true;
      if(witness_grid[w][1] == l2)
        found_12 = true;
    }
    else if(witness_grid[w][0] == l2)
    {
      found_2 = true;
      if(witness_grid[w][1] == l1)
        found_12 = true;
      else if(witness_grid[w][1] == l0)
        found_02 = true;
    }
  }

  return found_1 && found_2 && found_12 && found_02;
}

template<typename MF>
void split_quad(const Quad q,
                std::list<Quad>& quads_to_test,
                std::vector<FT>& values,
                std::vector<Point_2>& points,
                const MF* mf, const Traits* traits)
{
  // compute the new 5 pts & their respective color
  // p1---p4---p0
  // |    |    |
  // p5---p6---p7
  // |    |    |
  // p2---p8---p3

  FT x0 = points[q.i1].x();
  FT x1 = points[q.i0].x();
  FT y0 = points[q.i2].y();
  FT y1 = points[q.i1].y();

  FT xmid = (x0+x1)/2.;
  FT ymid = (y0+y1)/2.;

  std::size_t i4 = points.size();
  Point_2 p4(xmid, y1);
  points.push_back(p4);
  values.push_back(value_at_point(p4, i4, mf, traits));

  std::size_t i5 = points.size();
  Point_2 p5(x0, ymid);
  points.push_back(p5);
  values.push_back(value_at_point(p5, i5, mf, traits));

  std::size_t i6 = points.size();
  Point_2 p6(xmid, ymid);
  points.push_back(p6);
  values.push_back(value_at_point(p6, i6, mf, traits));

  std::size_t i7 = points.size();
  Point_2 p7(x1, ymid);
  points.push_back(p7);
  values.push_back(value_at_point(p7, i7, mf, traits));

  std::size_t i8 = points.size();
  Point_2 p8(xmid, y0);
  points.push_back(p8);
  values.push_back(value_at_point(p8, i8, mf, traits));

  // create the new four quads and add them to the queue
  Quad q1(q.i0, i4, i6, i7);
  Quad q2(i4, q.i1, i5, i6);
  Quad q3(i6, i5, q.i2, i8);
  Quad q4(i7, i6, i8, q.i3);

  quads_to_test.push_back(q1);
  quads_to_test.push_back(q2);
  quads_to_test.push_back(q3);
  quads_to_test.push_back(q4);
}

template<typename MF>
void smart_grid(const MF* mf, const Traits* traits)
{
  // idea is to create some kind of quadtree and refine a square
  // if its four corners don't have all the same colors (and it's not too small)
  FT min_vol = grid_side*grid_side*1e-7;
  FT max_vol = grid_side*grid_side*1e-3;

  std::list<Quad> quads_to_test;
  std::list<Quad> final_quads;
  std::vector<Point_2> points;
  std::vector<FT> values;

  // create the first quad
  FT l = grid_side / 2.;
  Point_2 p0(l,l), p1(-l,l), p2(-l,-l), p3(l,-l);
  points.push_back(p0);
  points.push_back(p1);
  points.push_back(p2);
  points.push_back(p3);
  values.push_back(value_at_point(p0, 0, mf, traits));
  values.push_back(value_at_point(p1, 1, mf, traits));
  values.push_back(value_at_point(p2, 2, mf, traits));
  values.push_back(value_at_point(p3, 3, mf, traits));

  Quad q0(0, 1, 2, 3);
  quads_to_test.push_back(q0);

  // refine the quad set
  while(!quads_to_test.empty())
  {
    const Quad& q = quads_to_test.front();

//    std::cout << quads_to_test.size() << std::endl;

    //experiment : refine more near 3pts dual
    if(false && !q.is_too_small(0.01*min_vol, points) && q.number_of_colors(values) >= 3)
      split_quad(q, quads_to_test, values, points, mf, traits);

    else if(q.is_too_small(min_vol, points) || q.has_same_colors(values))
      final_quads.push_back(q);

    else if(/*q.is_too_big(max_vol, points) ||*/ !q.has_same_colors(values))
      split_quad(q, quads_to_test, values, points, mf, traits);

    quads_to_test.pop_front();
  }

  // output the quad set (easily transformed into a triangulation) -------------
  CGAL_assertion(points.size() == values.size());
  std::size_t grid_side = values.size();
  std::cout << "check: " << points.size() << " " << values.size() << " vertices" << std::endl;

  std::ofstream out("smart_grid.mesh");
  out << "MeshVersionFormatted 1" << std::endl;
  out << "Dimension 2" << std::endl;
  out << "Vertices" << std::endl;
  out << grid_side << std::endl;

  std::ofstream out_bb("smart_grid.bb");
  out_bb << "2 1 " << grid_side << " 2" << std::endl;

  int counter = 0;
  for(int i=0; i!=grid_side; ++i)
  {
    const Point_2& p = points[i];
    out << p.x() << " " << p.y() << " " << ++counter << std::endl;
    out_bb << values[i] << std::endl;
  }
  out_bb << "End" << std::endl;

  out << "Triangles" << std::endl;
  out << 2*final_quads.size() << std::endl;
  std::list<Quad>::iterator it = final_quads.begin(), iend = final_quads.end();
  for(; it!=iend; ++it)
  {
    const Quad& q = *it;
    out << q.i0+1 << " " << q.i2+1 << " " << q.i3+1 << " 1" << std::endl;
    out << q.i0+1 << " " << q.i1+1 << " " << q.i2+1 << " 2" << std::endl;
  }
  out << "End" << std::endl;

  //smart dual : if a quad has 3 different values, the dual exists -------------
#if 0//def SMART_DUAL
  simplices.clear();
  for(it=final_quads.begin(); it!=iend; ++it)
  {
    const Quad& q = *it;
    std::set<std::size_t> simplex;
    simplex.insert(values[q.i0]); simplex.insert(values[q.i1]);
    simplex.insert(values[q.i2]);
    if(simplex.size() == 3)
      simplices.insert(simplex);
    else
    {
      simplex.insert(values[q.i3]);
      if(simplex.size() == 3)
        simplices.insert(simplex);
    }
  }
#else
  //witness --------------------------------------------------------------------
  CGAL_assertion(witness_grid.size() == grid_side);
  simplices.clear();
  for(std::size_t w=0; w<witness_grid.size(); ++w)
  {
    std::cout << w << " out of " << witness_grid.size() << std::endl;
    std::size_t l0 = witness_grid[w][0];
    std::size_t l1 = witness_grid[w][1];
    std::size_t l2 = witness_grid[w][2];
    std::set<std::size_t> simplex;
    simplex.insert(l0); simplex.insert(l1); simplex.insert(l2);

    if(simplices.find(simplex) != simplices.end()) // it's already in the complex
      continue;

    if(is_simplex_in_witness_complex(l0, l1, l2))
    {
      std::cout << "added a new simple: " << simplices.size() << std::endl;
      simplices.insert(simplex);
    }
  }
#endif
}

void output_simplices()
{
  std::cout << "captured: " << simplices.size() << " simplices" << std::endl;

  std::ofstream outd("grid_dual.mesh");
  outd << "MeshVersionFormatted 1" << std::endl;
  outd << "Dimension 2" << std::endl;
  outd << "Vertices" << std::endl;
  outd << vertices_nv << std::endl;
  for(int i=0; i<vertices_nv; ++i)
#ifdef R2
    outd << R2seeds[i].x() << " " << R2seeds[i].y() << " " << i+1 << std::endl;
#else
    outd << seeds[i].x() << " " << seeds[i].y() << " " << i+1 << std::endl;
#endif

  outd << "Triangles" << std::endl;
  outd << simplices.size() << std::endl;
  for(typename std::set<Simplex>::iterator it = simplices.begin();
      it != simplices.end(); ++it)
  {
    const Simplex& s = *it;
    typename Simplex::iterator sit = s.begin();
    typename Simplex::iterator siend = s.end();
    for(; sit!=siend; ++sit)
    {
//      std::cout << *sit+1 << " ";
      outd << *sit+1 << " ";
    }
    outd << "1" << std::endl;
//    std::cout << std::endl;
  }
  outd << "End" << std::endl;
}

template<typename MF>
void draw(const MF* metric_field)
{
  Traits* traits = new Traits();

#ifdef R2
  vertices_nv = read_tangent_plane_seeds();
#else
  vertices_nv = build_seeds(metric_field);
#endif

  // for random colors
  random_colors.resize(vertices_nv);
  for(int i=0; i<vertices_nv; ++i)
    random_colors[i] = i;
//  std::random_shuffle(random_colors.begin(), random_colors.end());

  first_position.resize(vertices_nv);
  second_position.resize(vertices_nv);

//  full_grid(metric_field, traits);
  smart_grid(metric_field, traits);

  output_simplices();
}

int main(int, char**)
{
  std::freopen("grid_log.txt", "w", stdout);
  std::srand(0);

  // to build a ~dual
  simplices.clear();
  witness_grid.clear();
  first_position.clear();
  second_position.clear();

  //Euclidean_metric_field<K>* metric_field = new Euclidean_metric_field<K>(16,1);
  Custom_metric_field<K>* metric_field = new Custom_metric_field<K>();

  draw(metric_field);
}

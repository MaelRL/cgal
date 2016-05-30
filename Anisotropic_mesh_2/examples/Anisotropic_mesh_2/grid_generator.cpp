// Voronoi painters for simple anisotropic distances (Du/Wang, Labelle/Shewchuk, etc.)

#include <CGAL/Starset.h>
#include <CGAL/Stretched_Delaunay_2.h>

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <Eigen/Dense>
#include <boost/array.hpp>
#include <omp.h>

#include <fstream>
#include <iostream>
#include <vector>

using namespace CGAL::Anisotropic_mesh_2;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef typename K::FT                                       FT;
typedef typename K::Segment_2                                Segment;
typedef typename K::Triangle_2                               Triangle;

typedef Stretched_Delaunay_2<K>                              Star;
typedef typename Star::Point_2                               Point_2;
typedef typename Star::TPoint_2                              TPoint_2;
typedef typename Star::Metric                                Metric;
typedef typename Star::Traits                                Traits;

typedef Eigen::Matrix<double, 5, 1>                          Vector5d;
typedef Eigen::Matrix<double, 2, 1>                          Vector2d;

typedef boost::array<std::size_t, 2>                         Edge;
typedef boost::array<std::size_t, 3>                         Tri;

//typedef CGAL::Anisotropic_mesh_2::Euclidean_metric_field<K>* MF;
typedef CGAL::Anisotropic_mesh_2::Custom_metric_field<K>* MF;

typedef CGAL::Exact_predicates_exact_constructions_kernel    KExact;
typedef typename KExact::Point_2                             EPoint;
typedef typename KExact::Segment_2                           ESegment;
typedef typename KExact::Triangle_2                          ETriangle;
typedef CGAL::Cartesian_converter<K, KExact>                 To_exact;
typedef CGAL::Cartesian_converter<KExact, K>                 Back_from_exact;

To_exact to_exact;
Back_from_exact back_from_exact;

// -----------------------------------------------------------------------------

// define (only) one of these for the dual
//#define APPROXIMATE_DUAL
#define SMART_DUAL
//#define WITNESS_DUAL

#define FILTER_SEEDS_OUTSIDE_GRID
//#define R2 // tangent plane mode

//#define TMP_REFINEMENT_UGLY_HACK

// -----------------------------------------------------------------------------

// using a lot of global variables, it's ugly but passing them by function is tedious
std::size_t vertices_nv = 100;

Point_2 center(0.,0.);
#ifdef R2
const FT grid_side = 100;
FT offset_x = -grid_side/2.; // offset is the bottom left point
FT offset_y = -grid_side/2.; // todo normalize this with aniso_mesh_2's rectangle whose offset is the center of the rectangle...
#else
const FT grid_side = 4.0; // length of a side
FT offset_x = center.x() - grid_side/2.; // offset is the bottom left point
FT offset_y = center.y() - grid_side/2.;
#endif
std::size_t max_grid_size = 1e6;

// metric & traits
MF mf;
Traits* traits;

// seeds
std::vector<Vector2d> R2seeds;
std::vector<FT> r2ws;
std::vector<Point_2> seeds;
std::vector<Vector5d> R5seeds;
std::vector<Metric> seeds_m;
std::vector<FT> ws;

std::set<Tri> simplices; // THESE ARE THE SIMPLICES OF THE DUAL OF THE GRID

// witness
std::vector<std::vector<std::size_t> > witness_grid;
std::vector<std::vector<std::size_t> > first_position; // could be a vector of booleans...
std::vector<std::vector<std::size_t> > second_position;

// refine
int n_refine = 25;
#ifdef TMP_REFINEMENT_UGLY_HACK
// farthest point memory
FT farthest_x = 1e30, farthest_y = 1e30;
FT farthest_d = 0.;
#endif

// ----------------- FUNCTIONS START -----------------------

std::set<Edge> build_edges()
{
  std::set<Edge> edges;
  typename std::set<Tri>::const_iterator it = simplices.begin(), iend = simplices.end();
  for(; it!=iend; ++it)
  {
    //the tris are sorted, so the edges are sorted as long as the indices are sorted!
    const Tri& tri = *it;
    Edge e;
    e[0] = tri[0]; e[1] = tri[1]; edges.insert(e);
    e[0] = tri[0]; e[1] = tri[2]; edges.insert(e);
    e[0] = tri[1]; e[1] = tri[2]; edges.insert(e);
  }

  return edges;
}

bool is_triangle_intersected(const Tri& tri, const std::set<Edge>& edges)
{
  bool is_intersected = false;

  // need to switch to epeck for the correct intersections here...
  const std::size_t i0 = tri[0];
  const std::size_t i1 = tri[1];
  const std::size_t i2 = tri[2];
  const Point_2& p0 = seeds[i0];
  const Point_2& p1 = seeds[i1];
  const Point_2& p2 = seeds[i2];
  const Triangle triangle(p0, p1, p2);

#pragma omp parallel shared(is_intersected, edges, p0, p1, p2)
{
  for(std::set<Edge>::const_iterator it = edges.begin(); it!=edges.end(); ++it)
  {
#pragma omp single nowait // hack to parallelize the std::set loop
{
#pragma omp flush (is_intersected)
    if(!is_intersected) // hack because we're not allowed to 'break' (nor to 'continue')
                        // inside a pragma omp for
    {
      const Edge& edge = *it;
      bool local_intersected = false;

      const Segment segment(seeds[edge[0]], seeds[edge[1]]);

      ESegment esegment = to_exact(segment);
      ETriangle etriangle = to_exact(triangle);
      CGAL::cpp11::result_of<KExact::Intersect_2(ESegment, ETriangle)>::type
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
{      if(local_intersected)
        is_intersected = true;
}
#pragma omp flush (is_intersected)
    } // end !is_intersected
} // end single nowait
  } // end for
} // end parallel
  return is_intersected;
}

// -----------------------------------------------------------------------------
// WITNESS FUNCTIONS
// -----------------------------------------------------------------------------

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

// -----------------------------------------------------------------------------
// SEED FUNCTIONS
// -----------------------------------------------------------------------------

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

int insert_new_seed(const FT x, const FT y)
{
#ifdef FILTER_SEEDS_OUTSIDE_GRID
    if(x < offset_x || x > center.x() + grid_side/2. ||
       y < offset_y || y > center.y() + grid_side/2. )
    {
      std::cout << "filtered : " << x << " " << y << std::endl;
      return seeds.size();
    }
#endif

    seeds.push_back(Point_2(x, y));
    const Point_2& p = seeds.back();
    seeds_m.push_back(mf->compute_metric(p));
    R5seeds.push_back(compute_hat(p, seeds_m.back()));
    const Vector5d& v = R5seeds.back();
    ws.push_back(v.norm() * v.norm() - p.x()*v(0) - p.y()*v(1));

//    std::cout << "Seeds i: " << seeds.size()-1 << " " << p.x() << " " << p.y() << std::endl;
//    std::cout << "check " << std::endl;
//    std::cout << seeds_m.back().get_mat() << std::endl;
//    std::cout << v(0) << " " << v(1) << std::endl;
//    std::cout << ws.back() << std::endl;
  return seeds.size();
}

int build_seeds()
{
  std::ifstream in(str_seed.c_str());
  CGAL_precondition(in);
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
  R5seeds.reserve(min_nv);
  seeds_m.reserve(min_nv);
  ws.reserve(min_nv);

  for(std::size_t i=0; i<nv; ++i)
  {
    in >> r_x >> r_y >> useless;
    insert_new_seed(r_x, r_y);

    if(seeds.size() >= vertices_nv)
      break;
  }
  std::cout << "seeds: " << seeds.size() << std::endl;
  return seeds.size();
}

int read_tangent_plane_seeds()
{
  std::ifstream in("/home/mrouxell/cgal/Anisotropic_mesh_TC/examples/Anisotropic_mesh_TC/build/projected_points.txt");
  //  std::ifstream in("/home/mrouxell/cgal/Anisotropic_mesh_TC/examples/Anisotropic_mesh_TC/build/projected_points_after.txt");
  //  std::ifstream in("/home/mrouxell/cgal/Tangential_complex/test/Tangential_complex/build/projected_points_984.txt");
  std::size_t nv, dim;
  FT x, y, w;

  in >> dim >> nv;
  assert(dim == 2);
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

// -----------------------------------------------------------------------------
// GRID FUNCTIONS
// -----------------------------------------------------------------------------

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

std::size_t value_at_point(const Point_2& p, std::size_t p_id)
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
  std::size_t min_id = 0;
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
        min_id = k;
        min = sq_d;
      }
      else if(l==1 && ((k!=min_id && std::abs(sq_d-min)<0.1*min) ||
                       sq_d < 0.01*min)) // second pass
      {
        min_id = 1.01*vertices_nv; // drawing borders and center points
        break;
      }
    }
  }

  Comp c(sqds);
  std::sort(ids.begin(), ids.end(), c);

#ifdef WITNESS_DUAL
  //for the witness
  witness_grid.push_back(std::vector<std::size_t>(ids.begin(), ids.begin()+3));
  first_position[ids[0]].push_back(p_id);
  second_position[ids[1]].push_back(p_id);
#elif defined(APPROXIMATE_DUAL)
  for(int i=0; i<vertices_nv-1; ++i)
  {
//    std::cout.precision(17);
//    std::cout << sqds[ids[i+1]] - sqds[ids[i]] << std::endl;
    assert(sqds[ids[i+1]] >= sqds[ids[i]]);
  }

  // check if we are in the dual of a triangle (somewhat!)
  // that is : the three closest points are more or less at the same distance
  FT alpha = 1e-1;
  if(std::abs(sqds[ids[0]] - sqds[ids[1]]) < alpha*sqds[ids[0]] &&
     std::abs(sqds[ids[0]] - sqds[ids[2]]) < alpha*sqds[ids[0]])
  {
    Tri t; t[0] = ids[0]; t[1] = ids[1]; t[2] = ids[2];
    simplices.insert(t);
    std::cout << "added one ! " << std::abs(sqds[ids[0]] - sqds[ids[2]]) << std::endl;
  }
#endif

  return min_id;
}

void full_grid()
{
  double n = 300.; // number of points per side
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

      FT ret = value_at_point(p, i*n+j);
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

// -----------------------------------------------------------------------------
// SMART GRID
// -----------------------------------------------------------------------------

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

  bool is_too_distorted(FT max_distortion, const std::vector<Point_2>& points) const
  {
    Metric m0 = mf->compute_metric(points[i0]);
    Metric m1 = mf->compute_metric(points[i1]);

    FT gamma01 = m0.compute_distortion(m1);
      if(gamma01 > max_distortion)
        return true;

    Metric m2 = mf->compute_metric(points[i2]);
    FT gamma02 = m0.compute_distortion(m2);
    if(gamma02 > max_distortion)
      return true;

    Metric m3 = mf->compute_metric(points[i3]);
    FT gamma03 = m0.compute_distortion(m3);
    if(gamma03 > max_distortion)
      return true;

    FT gamma12 = m1.compute_distortion(m2);
    if(gamma12 > max_distortion)
      return true;
    FT gamma13 = m1.compute_distortion(m3);
    if(gamma13 > max_distortion)
      return true;

    FT gamma23 = m2.compute_distortion(m3);
    if(gamma23 > max_distortion)
      return true;

    return false;
  }

  bool has_same_colors(const std::vector<std::size_t>& values) const
  {
    CGAL_assertion(i0 < values.size() && i1 < values.size() &&
                   i2 < values.size() && i3 < values.size());
    FT v0 = values[i0];
    if(v0 != values[i1])
      return false;
    if(v0 != values[i2])
      return false;
    if(v0 != values[i3])
      return false;
    return true;
  }

  std::size_t number_of_colors(const std::vector<std::size_t>& values) const
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

void split_quad(const Quad q,
                std::list<Quad>& quads_to_test,
                std::vector<std::size_t>& values,
                std::vector<Point_2>& points)
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
  points.push_back(p4); values.push_back(value_at_point(p4, i4));

  std::size_t i5 = points.size();
  Point_2 p5(x0, ymid);
  points.push_back(p5); values.push_back(value_at_point(p5, i5));

  std::size_t i6 = points.size();
  Point_2 p6(xmid, ymid);
  points.push_back(p6);  values.push_back(value_at_point(p6, i6));

  std::size_t i7 = points.size();
  Point_2 p7(x1, ymid);
  points.push_back(p7);  values.push_back(value_at_point(p7, i7));

  std::size_t i8 = points.size();
  Point_2 p8(xmid, y0);
  points.push_back(p8);  values.push_back(value_at_point(p8, i8));

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

#ifdef TMP_REFINEMENT_UGLY_HACK
void ugly_farthest_computation(const Tri& tri,
                               const std::vector<Point_2>& points,
                               const std::vector<std::size_t>& values)
{
  //tri is a triangle with different colors at each vertex
  typename Star::Traits::Compute_squared_distance_2 csd =
      traits->compute_squared_distance_2_object();

  for(std::size_t i=0; i<tri.size(); ++i)
  {
    const Point_2& p = points[tri[i]];
    std::size_t k = values[tri[i]];
    CGAL_assertion(k < seeds.size());

    // DU WANG : dx(p,x) < dx(q,x)
    const Metric& m_p = mf->compute_metric(p);
    TPoint_2 tp(m_p.transform(p));
    double sq_d = csd(tp, m_p.transform(seeds[k]));

    // LABELLE SHEWCHUK : dp(p,x) < dq(q,x)
    TPoint_2 tp2 = (seeds_m[k]).transform(p);
    TPoint_2 ts = (seeds_m[k]).transform(seeds[k]);
//    double sq_d = csd(tp2, ts);

    if(sq_d > farthest_d)
    {
#pragma omp critical
{
      farthest_d = sq_d;
      farthest_x = p.x();
      farthest_y = p.y();
}
    }
  }
}
#endif

void insert_simplex_if_colored(const std::size_t i0, const std::size_t i1,
                               const std::size_t i2,
                               const std::vector<std::size_t>& values,
                               const std::vector<Point_2>& points)
{
  if(values[i0] != values[i1] && values[i0] != values[i2] &&
     values[i1] != values[i2])
  {
    Tri t;
    t[0] = values[i0]; t[1] = values[i1]; t[2] = values[i2];
#pragma omp critical
    simplices.insert(t);
#ifdef TMP_REFINEMENT_UGLY_HACK
    t[0] = i0; t[1] = i1; t[2] = i2;
    ugly_farthest_computation(t, points, values);
#endif
  }
}

void output_smart_grid(const std::list<Quad>& final_quads,
                       const std::vector<Point_2>& points,
                       const std::vector<std::size_t>& values)
{
  // OUTPUT the quad set (easily transformed into a triangulation by adding the diagonal)
  CGAL_assertion(points.size() == values.size());
  std::size_t values_n = values.size();
  std::cout << "check: " << points.size() << " " << values.size() << " vertices" << std::endl;

  std::ofstream out("smart_grid_LS.mesh");
  out << "MeshVersionFormatted 1" << '\n';
  out << "Dimension 2" << '\n';
  out << "Vertices" << '\n';
  out << values_n << '\n';

  std::ofstream out_bb("smart_grid_LS.bb");
  out_bb << "2 1 " << values_n << " 2" << '\n';

  int counter = 0;
  for(std::size_t i=0; i!=values_n; ++i)
  {
    const Point_2& p = points[i];
    out << p.x() << " " << p.y() << " " << ++counter << '\n';
    out_bb << values[i] << '\n';
  }
  out_bb << "End" << '\n';

  out << "Triangles" << '\n';
  out << 2*final_quads.size() << '\n';
  std::list<Quad>::const_iterator it = final_quads.begin(), iend = final_quads.end();
  for(; it!=iend; ++it)
  {
    const Quad& q = *it;
    out << q.i0+1 << " " << q.i2+1 << " " << q.i3+1 << " 1" << '\n';
    out << q.i0+1 << " " << q.i1+1 << " " << q.i2+1 << " 2" << '\n';
  }
  out << "End" << '\n';
}

// todo parallelize this like it was done for grid_gen_3, be careful with the
// pragma omp critical in other functions!
void adapted_grid(const bool refine,
                  const bool output = true)
{
  std::cout << "smart grid !" << std::endl;
#ifdef TMP_REFINEMENT_UGLY_HACK
  if(refine)
  {
    CGAL_precondition(farthest_x != 1e30 && farthest_y != 1e30);
    vertices_nv = insert_new_seed(farthest_x, farthest_y);
    farthest_d = 0.; farthest_x = 1e30; farthest_y = 1e30;
  }
#endif

  // idea is to create some kind of quadtree and refine a quad following
  // a criterion based on the value at its vertices
  FT min_vol = grid_side*grid_side*1e-6;
  FT max_vol = grid_side*grid_side*1e-3;

  std::list<Quad> quads_to_test;
  std::list<Quad> final_quads; // quads that won't be refined anymore
  std::vector<Point_2> points;
  std::vector<std::size_t> values;

  // create the initial quad
  FT l = grid_side / 2.;
  FT x0 = center.x() - l, x1 = center.x() + l;
  FT y0 = center.y() - l, y1 = center.y() + l;

  Point_2 p0(x1,y1), p1(x0,y1), p2(x0,y0), p3(x1,y0);
  points.push_back(p0);
  points.push_back(p1);
  points.push_back(p2);
  points.push_back(p3);
  values.push_back(value_at_point(p0, 0));
  values.push_back(value_at_point(p1, 1));
  values.push_back(value_at_point(p2, 2));
  values.push_back(value_at_point(p3, 3));

  Quad q0(0, 1, 2, 3);
  quads_to_test.push_back(q0);

  // refine the quad set
  while(!quads_to_test.empty())
  {
    const Quad& q = quads_to_test.front();

    ///experiment : refine more near 3pts dual
//    if(!q.is_too_small(0.01*min_vol, points) && q.number_of_colors(values) >= 3)
//      split_quad(q, quads_to_test, values, points, mf, traits);
    ///experiment : refine for the distortion
//if(points.size() < 1e6 && (q.is_too_distorted(1.1, points) || points.size() < 100)

    if(q.is_too_small(min_vol, points) || points.size() > max_grid_size)
    {
      final_quads.push_back(q);
    }
    else if(q.is_too_big(max_vol, points))
    {
       split_quad(q, quads_to_test, values, points);
    }
    else if(!q.has_same_colors(values))
    {
      split_quad(q, quads_to_test, values, points);
    }
    else
      final_quads.push_back(q);

    quads_to_test.pop_front();
  }

  if(output)
    output_smart_grid(final_quads, points, values);

  //smart dual : if a quad has 3 different values, the dual exists
#ifdef SMART_DUAL
  simplices.clear();
  for(std::list<Quad>::iterator it=final_quads.begin(); it!=final_quads.end(); ++it)
  {
    const Quad& q = *it;
    insert_simplex_if_colored(q.i0, q.i1, q.i2, values, points);
    insert_simplex_if_colored(q.i0, q.i2, q.i3, values, points);
  }

#elif defined(WITNESS_DUAL)
  //witness
  CGAL_assertion(witness_grid.size() == points.size());
  simplices.clear();
  for(std::size_t w=0; w<witness_grid.size(); ++w)
  {
    std::cout << w << " out of " << witness_grid.size() << std::endl;
    std::size_t l0 = witness_grid[w][0];
    std::size_t l1 = witness_grid[w][1];
    std::size_t l2 = witness_grid[w][2];
    Tri t;
    t[0] = l0; t[1] = l1; t[2] = l2;

    if(simplices.find(t) != simplices.end()) // it's already in the complex
      continue;

    if(is_simplex_in_witness_complex(l0, l1, l2))
    {
      std::cout << "added a new simple: " << simplices.size() << std::endl;
      simplices.insert(t);
    }
  }
#endif
}

// -----------------------------------------------------------------------------
// OUTPUT FUNCTIONS
// -----------------------------------------------------------------------------
void output_simplices()
{
  std::cout << "captured: " << simplices.size() << " simplices" << std::endl;

  std::ofstream outd("grid_dual_LS.mesh");
  outd << "MeshVersionFormatted 1" << '\n';
  outd << "Dimension 2" << '\n';
  outd << "Vertices" << '\n';
  outd << seeds.size() << '\n';
  for(std::size_t i=0; i<seeds.size(); ++i)
#ifdef R2
    outd << R2seeds[i].x() << " " << R2seeds[i].y() << " " << i+1 << '\n';
#else
    outd << seeds[i].x() << " " << seeds[i].y() << " " << i+1 << '\n';
#endif

  const std::set<Edge>& edges = build_edges();

  outd << "Triangles" << '\n';
  outd << simplices.size() << '\n';
  for(typename std::set<Tri>::iterator it = simplices.begin();
      it != simplices.end(); ++it)
  {
    const Tri& tri = *it;
    for(std::size_t i=0; i<tri.size(); ++i)
      outd << tri[i]+1 << " ";
    outd << is_triangle_intersected(tri, edges) << '\n';
  }
  outd << "End" << '\n';
}

// -----------------------------------------------------------------------------
// MAIN FUNCTIONS
// -----------------------------------------------------------------------------
void initialize()
{
  traits = new Traits();
  //  mf = new Euclidean_metric_field<K>(16,1);
  mf = new Custom_metric_field<K>();

#ifdef R2
  vertices_nv = read_tangent_plane_seeds();
#else
  vertices_nv = build_seeds();
#endif
}

void build_grid(const bool refine = false,
                const bool output = true)
{
  simplices.clear();
#ifdef WITNESS_DUAL
  witness_grid.clear();
  first_position.resize(vertices_nv, 0.);
  second_position.resize(vertices_nv, 0.);
#endif

  adapted_grid(refine);
  if(output)
    output_simplices();
}

int main(int, char**)
{
//  std::freopen("log.txt", "w", stdout);
  std::srand(0);

  initialize();
  build_grid();

  for(int i=0; i<n_refine; ++i)
  {
    std::cout << "refine: " << i << std::endl;
    build_grid(true /*refine*/, (i%100==0) /*output*/);
  }

  std::cout << "end of program" << std::endl;
}

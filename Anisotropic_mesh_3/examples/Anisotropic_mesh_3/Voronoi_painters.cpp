// Voronoi painters for simple anisotropic distances (Du/Wang, Labelle/Shewchuk, etc.)

#include <CGAL/Starset.h>
#include <CGAL/Stretched_Delaunay_3.h>

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <boost/array.hpp>
#include <Eigen/Dense>
#include <omp.h>

#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef typename K::FT                                       FT;
typedef typename K::Segment_3                                Segment;
typedef typename K::Triangle_3                               Triangle;

typedef Stretched_Delaunay_3<K>                              Star;
typedef typename Star::Point_3                               Point_3;
typedef typename Star::TPoint_3                              TPoint_3;
typedef typename Star::Metric                                Metric;
typedef typename Star::Traits                                Traits;

//typedef typename CGAL::Anisotropic_mesh_3::Euclidean_metric_field<K>* MF;
typedef typename CGAL::Anisotropic_mesh_3::Custom_metric_field<K>* MF;

typedef std::set<std::size_t>                                Simplex;
typedef boost::array<std::size_t, 2>                         Edge;
typedef boost::array<std::size_t, 3>                         Tri;
typedef boost::array<std::size_t, 4>                         Tet;

typedef CGAL::Exact_predicates_exact_constructions_kernel    KExact;
typedef typename KExact::Point_3                             EPoint;
typedef typename KExact::Segment_3                           ESegment;
typedef typename KExact::Triangle_3                          ETriangle;
typedef CGAL::Cartesian_converter<K, KExact>                 To_exact;
typedef CGAL::Cartesian_converter<KExact, K>                 Back_from_exact;

To_exact to_exact;
Back_from_exact back_from_exact;

// -----------------------------------------------------------------------------

// pick (only) one to compute the dual
//#define APPROXIMATE_DUAL
#define SMART_DUAL

//#define REFINE_NEAR_CENTERS_ONLY
#define TMP_REFINEMENT_UGLY_HACK
#define FILTER_SEEDS_OUTSIDE_GRID

// -----------------------------------------------------------------------------

// using a lot of global variables, it's ugly but passing them by function is tedious
std::size_t vertices_nv = 10;
std::string str_seeds = "bambimboum_wip.mesh";

Point_3 center(0.6, 0.6, 0.6);
const FT grid_side = 1.0;
FT offset_x = center.x() - grid_side/2.; // offset is the bottom left point
FT offset_y = center.y() - grid_side/2.;
FT offset_z = center.z() - grid_side/2.;
std::size_t max_local_grid_nv = 1e5;

// metric & traits
MF mf;
Traits* traits;

// seeds
std::vector<Point_3> seeds;
std::vector<Metric> seeds_m;

std::set<Tet> simplices; // THESE ARE THE SIMPLICES OF THE DUAL OF THE GRID

// refine stuff
std::size_t n_refine = 0;
#ifdef TMP_REFINEMENT_UGLY_HACK
// farthest point memory
FT farthest_x = 1e30, farthest_y = 1e30, farthest_z = 1e30;
FT farthest_d = 0.;
#endif

// -----------------------------------------------------------------------------

std::set<Edge> compute_edges()
{
  std::set<Edge> edges;
  for(typename std::set<Tet>::const_iterator it = simplices.begin();
                                             it != simplices.end(); ++it)
  {
    //the tets are sorted, so the edges are sorted as long as the indices are sorted!
    const Tet& tet = *it;
    Edge e;
    e[0] = tet[0]; e[1] = tet[1]; edges.insert(e);
    e[0] = tet[0]; e[1] = tet[2]; edges.insert(e);
    e[0] = tet[0]; e[1] = tet[3]; edges.insert(e);
    e[0] = tet[1]; e[1] = tet[2]; edges.insert(e);
    e[0] = tet[1]; e[1] = tet[3]; edges.insert(e);
    e[0] = tet[2]; e[1] = tet[3]; edges.insert(e);
  }

  return edges;
}

std::set<Tri> compute_triangles()
{
  std::set<Tri> triangles;
  for(typename std::set<Tet>::const_iterator it = simplices.begin();
                                             it != simplices.end(); ++it)
  {
    //the tets are sorted, so the edges are sorted as long as the indices are sorted!
    const Tet& tet = *it;
    Tri t;
    t[0] = tet[0]; t[1] = tet[1]; t[2] = tet[2]; triangles.insert(t);
    t[0] = tet[0]; t[1] = tet[1]; t[2] = tet[3]; triangles.insert(t);
    t[0] = tet[0]; t[1] = tet[2]; t[2] = tet[3]; triangles.insert(t);
    t[0] = tet[1]; t[1] = tet[2]; t[2] = tet[3]; triangles.insert(t);
  }

  return triangles;
}

bool has_simplex_n_colors(const std::size_t i0, const std::size_t i1,
                          const std::size_t i2, const std::size_t i3,
                          const std::size_t n, const std::vector<std::size_t>& values)
{
  std::set<std::size_t> colors;
  colors.insert(values[i0]);
  colors.insert(values[i1]);
  colors.insert(values[i2]);
  colors.insert(values[i3]);

  return (colors.size() >= n);
}

bool is_triangle_intersected(const Tri& tri, const std::set<Edge>& edges)
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

#pragma omp parallel shared(is_intersected, p0, p1, p2)
{
  for(std::set<Edge>::const_iterator it = edges.begin(); it!=edges.end(); ++it)
  {
#pragma omp single nowait // hack to parallelize the std::set
{
#pragma omp flush (is_intersected)
    if(!is_intersected) // hack because we're not allowed to break (or continue)
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

int insert_new_seed(const FT x, const FT y, const FT z)
{
#ifdef FILTER_SEEDS_OUTSIDE_GRID
  if(x < offset_x || x > center.x() + grid_side/2. ||
     y < offset_y || y > center.y() + grid_side/2. ||
     z < offset_z || z > center.y() + grid_side/2.)
  {
    std::cout << "filtered : " << x << " " << y << " " << z << std::endl;
    return seeds.size();
  }
#endif
  seeds.push_back(Point_3(x, y, z));
  seeds_m.push_back(mf->compute_metric(seeds.back()));
  std::cout << "inserted " << x << " " << y << " " << z << std::endl;

  return seeds.size();
}

int build_seeds()
{
  std::ifstream in(str_seeds.c_str());
  std::string word;
  std::size_t useless, nv, dim;
  FT r_x, r_y, r_z;

  in >> word >> useless; //MeshVersionFormatted i
  in >> word >> dim; //Dimension d
  in >> word >> nv;
  std::cout << "nv: " << nv << std::endl;
  assert(dim == 3);

  std::size_t min_nv = (std::min)(nv, vertices_nv);
  seeds.reserve(min_nv);
  seeds_m.reserve(min_nv);

  for(std::size_t i=0; i<nv; ++i)
  {
    in >> r_x >> r_y >> r_z >> useless;
    insert_new_seed(r_x, r_y, r_z);

    if(seeds.size() >= vertices_nv)
      break;
  }
  std::cout << seeds.size() << " seeds" << std::endl;
  return seeds.size();
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

FT value_at_point(const Point_3& p)
{
  Metric m_p = mf->compute_metric(p);
  TPoint_3 tp(m_p.transform(p));
  typename Star::Traits::Compute_squared_distance_3 csd =
      traits->compute_squared_distance_3_object();

  FT min = 1e30;
  std::size_t min_id = 0;
  std::vector<FT> sqds(vertices_nv, 1e30);
  std::vector<std::size_t> ids(vertices_nv);

  for(int l=0; l<1; ++l) // set l = 2 to draw borders
  {
    for(std::size_t k=0; k<vertices_nv; ++k)
    {
      // DU WANG : dx(p,x) < dx(q,x)
      double sq_d = csd(tp, m_p.transform(seeds[k]));

      // LABELLE SHEWCHUK : dp(p,x) < dq(q,x)
      TPoint_3 tp2 = (seeds_m[k]).transform(p);
      TPoint_3 ts = (seeds_m[k]).transform(seeds[k]);
//      double sq_d = csd(tp2, ts);

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
        min_id = 1.01 * vertices_nv; // drawing borders and center points
        break;
      }
    }
  }

#ifdef APPROXIMATE_DUAL
  if(ids.size() < 4)
    return min_id;

  Comp c(sqds);
  std::sort(ids.begin(), ids.end(), c);

  for(int i=0; i<vertices_nv-1; ++i)
  {
//    std::cout << sqds[ids[i]] << " " << sqds[ids[i+1]] << std::endl;
    CGAL_assertion(sqds[ids[i+1]] >= sqds[ids[i]]);
  }

  // check if we are in the dual of a tet (somewhat, not exactly a witness)
  FT alpha = 0.1;
  if(std::abs(sqds[ids[0]] - sqds[ids[1]]) < alpha*sqds[ids[0]] &&
     std::abs(sqds[ids[0]] - sqds[ids[2]]) < alpha*sqds[ids[0]] &&
     std::abs(sqds[ids[0]] - sqds[ids[3]]) < alpha*sqds[ids[0]]) // in the dual of the tet 0123
  {
    if(ids.size() > 4)
    {
      if(std::abs(sqds[ids[0]] - sqds[ids[4]]) < alpha*sqds[ids[0]]) // 4 is close too!
      {
#pragma omp critical
{
        std::cout << "close tet" << std::endl;
        std::cout << sqds[ids[0]] << " " << sqds[ids[1]] << " " << sqds[ids[2]] << " " << sqds[ids[3]] << " ";
        std::cout << sqds[ids[4]] << std::endl;
        std::cout << "YOU VE GOT A CASE OF COSPHERITY at : " << p.x() << " " << p.y() << " " << p.z() << " ";
        std::cout << "(r: " << std::sqrt(p.x()*p.x() + p.y()*p.y() + p.z()*p.z()) << ")" << std::endl;
      }
}
    }

    Tet t; t[0] = ids[0]; t[1] = ids[1]; t[2] = ids[2]; t[3] = ids[3];
    std::sort(t.begin(), t.end());
    simplices.insert(simplex);
  }
#endif

  return min_id;
}

void full_grid()
{
  double n = 100.; // number of points per side
  double a = grid_side; // length of the side
  double step = a/n;
  double n_cube = n*n*n;
  double tet_n = -1;
  std::size_t counter = 0;

  std::cout << "full grid with: " << n_cube << " vertices" << std::endl;

  std::ofstream out("grid.mesh");
  out << "MeshVersionFormatted 1" << std::endl;
  out << "Dimension 3" << std::endl;
  out << "Vertices" << std::endl;
  out << (int) n_cube << std::endl;

  std::ofstream out_bb("grid.bb");
  out_bb << "3 1 " << (int) n_cube << " 2" << std::endl;

  for(unsigned int k=0; k<n; ++k)
  {
    for(unsigned int j=0; j<n; ++j)
    {
      for(unsigned int i=0; i<n; ++i)
      {
        Point_3 p(offset_x + i*step, offset_y + j*step, offset_z + k*step);
        out << p.x() << " " << p.y() << " " << p.z() << " " << ++counter << std::endl;

        FT ret = value_at_point(p);
        out_bb << ret << std::endl;
      }
    }
  }
  out << "Tetrahedra" << std::endl;
  out << 5*(n-1)*(n-1)*(n-1) << std::endl;
  for(std::size_t i=0; i<(n-1); ++i)
  {
    for(std::size_t j=0; j<(n-1); ++j)
    {
      for(std::size_t j=0; j<(n-1); ++j)
      {
        // define the cube 01234567

        // + 1 at the end for medit
        std::size_t id_0 = i + j*n + k*n*n + 1;
        std::size_t id_1 = i+1 + j*n + k*n*n + 1;
        std::size_t id_2 = i+1 + (j+1)*n + k*n*n + 1;
        std::size_t id_3 = i + (j+1)*n + k*n*n + 1;
        std::size_t id_4 = i + j*n + (k+1)*n*n + 1;
        std::size_t id_5 = i+1 + (j+1)*n + (k+1)*n*n + 1;
        std::size_t id_6 = i+1 + (j+1)*n + (k+1)*n*n + 1;
        std::size_t id_7 = i + j*n + (k+1)*n*n + 1;

        out << id_0 << " " << id_1 << " " << id_2 << " " << id_5 << " " << i+j*n+k*n*n << std::endl;
        out << id_0 << " " << id_2 << " " << id_5 << " " << id_7 << " " << i+j*n+k*n*n << std::endl;
        out << id_0 << " " << id_2 << " " << id_3 << " " << id_7 << " " << i+j*n+k*n*n << std::endl;
        out << id_0 << " " << id_4 << " " << id_5 << " " << id_7 << " " << i+j*n+k*n*n << std::endl;
        out << id_2 << " " << id_5 << " " << id_6 << " " << id_7 << " " << i+j*n+k*n*n << std::endl;
      }
    }
  }
  out << "End" << std::endl;
}

struct Cube
{
  // ^ y
  // |
  // 0--> x

  // top level
  // p1---p0
  // |    |
  // p2---p3

  //low level
  // p5---p4
  // |    |
  // p6---p7

  boost::array<std::size_t, 8> ids;

  bool is_too_small(FT min_vol, const std::vector<Point_3>& points) const
  {
    FT dx = points[0].x() - points[1].x();
    FT dy = points[0].y() - points[3].y();
    FT dz = points[0].z() - points[4].z();
//    if(dx*dy*dz < min_vol)
//      std::cout << "too small " << std::endl;
    return dx*dy*dz < min_vol;
  }

  bool is_too_big(FT max_vol, const std::vector<Point_3>& points) const
  {
    FT dx = points[0].x() - points[1].x();
    FT dy = points[0].y() - points[3].y();
    FT dz = points[0].z() - points[4].z();
    return dx*dy*dz > max_vol;
  }

  std::size_t number_of_colors(const std::vector<std::size_t>& values) const
  {
    std::set<FT> vals;
    boost::array<std::size_t, 8>::const_iterator it = ids.begin(), iend = ids.end();
    for(; it!=iend; ++it)
      vals.insert(values[*it]);
    return vals.size();
  }

  bool has_same_colors(const std::vector<std::size_t>& values) const
  {
    FT v0 = values[ids.front()];
    boost::array<std::size_t, 8>::const_iterator it = ids.begin(), iend = ids.end();
    ++it; // ignore the first one since it's v0
    for(; it!=iend; ++it)
      if(v0 != values[*it])
        return false;
    return true;
  }

  Cube(const boost::array<std::size_t, 8>& ids_) : ids(ids_) { }
  Cube(const std::size_t i0, const std::size_t i1, const std::size_t i2,
       const std::size_t i3, const std::size_t i4 ,const std::size_t i5,
       const std::size_t i6, const std::size_t i7) : ids()
  {
    ids[0] = i0; ids[1] = i1; ids[2] = i2; ids[3] = i3;
    ids[4] = i4; ids[5] = i5; ids[6] = i6; ids[7] = i7;
  }
};

void split_cube(const Cube q,
                std::list<Cube>& cubes_to_test,
                std::vector<std::size_t>& values,
                std::vector<Point_3>& points)
{
  // compute the new 5 pts & their respective color
  //top level (order defined from Oz)
  // p1---p8---p0
  // |     |    |
  // p11---p10---p9
  // |     |    |
  // p2---p12---p3

  //mid level
  // p15---p14---p13
  // |      |    |
  // p18---p17---p16
  // |      |    |
  // p21---p20---p19

  //lowest level
  // p5---p22---p4
  // |      |    |
  // p25---p24---p23
  // |      |    |
  // p6---p26---p7

  //gather the interesting values (2 of each)
  FT x0 = points[q.ids[1]].x(); FT x1 = points[q.ids[0]].x();
  FT y0 = points[q.ids[2]].y(); FT y1 = points[q.ids[1]].y();
  FT z0 = points[q.ids[4]].z(); FT z1 = points[q.ids[0]].z();

  FT xmid = (x0 + x1) / 2.;
  FT ymid = (y0 + y1) / 2.;
  FT zmid = (z0 + z1) / 2.;

  //new pts coordinates (19!):
  std::size_t i8 = points.size();
  Point_3 p8(xmid, y1, z1);
  points.push_back(p8); values.push_back(value_at_point(p8));

  std::size_t i9 = points.size();
  Point_3 p9(x1, ymid, z1);
  points.push_back(p9); values.push_back(value_at_point(p9));

  std::size_t i10 = points.size();
  Point_3 p10(xmid, ymid, z1);
  points.push_back(p10); values.push_back(value_at_point(p10));

  std::size_t i11 = points.size();
  Point_3 p11(x0, ymid, z1);
  points.push_back(p11); values.push_back(value_at_point(p11));

  std::size_t i12 = points.size();
  Point_3 p12(xmid, y0, z1);
  points.push_back(p12); values.push_back(value_at_point(p12));

  std::size_t i13 = points.size();
  Point_3 p13(x1, y1, zmid);
  points.push_back(p13); values.push_back(value_at_point(p13));

  std::size_t i14 = points.size();
  Point_3 p14(xmid, y1, zmid);
  points.push_back(p14); values.push_back(value_at_point(p14));

  std::size_t i15 = points.size();
  Point_3 p15(x0, y1, zmid);
  points.push_back(p15); values.push_back(value_at_point(p15));

  std::size_t i16 = points.size();
  Point_3 p16(x1, ymid, zmid);
  points.push_back(p16); values.push_back(value_at_point(p16));

  std::size_t i17 = points.size();
  Point_3 p17(xmid, ymid, zmid);
  points.push_back(p17); values.push_back(value_at_point(p17));

  std::size_t i18 = points.size();
  Point_3 p18(x0, ymid, zmid);
  points.push_back(p18); values.push_back(value_at_point(p18));

  std::size_t i19 = points.size();
  Point_3 p19(x1, y0, zmid);
  points.push_back(p19); values.push_back(value_at_point(p19));

  std::size_t i20 = points.size();
  Point_3 p20(xmid, y0, zmid);
  points.push_back(p20); values.push_back(value_at_point(p20));

  std::size_t i21 = points.size();
  Point_3 p21(x0, y0, zmid);
  points.push_back(p21); values.push_back(value_at_point(p21));

  std::size_t i22 = points.size();
  Point_3 p22(xmid, y1, z0);
  points.push_back(p22); values.push_back(value_at_point(p22));

  std::size_t i23 = points.size();
  Point_3 p23(x1, ymid, z0);
  points.push_back(p23); values.push_back(value_at_point(p23));

  std::size_t i24 = points.size();
  Point_3 p24(xmid, ymid, z0);
  points.push_back(p24); values.push_back(value_at_point(p24));

  std::size_t i25 = points.size();
  Point_3 p25(x0, ymid, z0);
  points.push_back(p25); values.push_back(value_at_point(p25));

  std::size_t i26 = points.size();
  Point_3 p26(xmid, y0, z0);
  points.push_back(p26); values.push_back(value_at_point(p26));

  // one cube gives 8 cubes :
  Cube c0(q.ids[0], i8, i10, i9, i13, i14, i17, i16);
  Cube c1(i8, q.ids[1], i11, i10, i14, i15, i18, i17);
  Cube c2(i10, i11, q.ids[2], i12, i17, i18, i21, i20);
  Cube c3(i9, i10, i12, q.ids[3], i16, i17, i20, i19);
  Cube c4(i13, i14, i17, i16, q.ids[4], i22, i24, i23);
  Cube c5(i14, i15, i18, i17, i22, q.ids[5], i25, i24);
  Cube c6(i17, i18, i21, i20, i24, i25, q.ids[6], i26);
  Cube c7(i16, i17, i20, i19, i23, i24, i26, q.ids[7]);

  cubes_to_test.push_back(c0);
  cubes_to_test.push_back(c1);
  cubes_to_test.push_back(c2);
  cubes_to_test.push_back(c3);
  cubes_to_test.push_back(c4);
  cubes_to_test.push_back(c5);
  cubes_to_test.push_back(c6);
  cubes_to_test.push_back(c7);
}

#ifdef TMP_REFINEMENT_UGLY_HACK
template<typename Simplex>
void ugly_farthest_computation(const Simplex& s,
                               const std::vector<Point_3>& points,
                               const std::vector<std::size_t>& values)
{
  //s is a _GRID_ simplex with different colors (values) at each vertex
  typename Star::Traits::Compute_squared_distance_2 csd =
      traits->compute_squared_distance_2_object();

  for(std::size_t i=0, ss=s.size(); i<ss; ++i)
  {
    const Point_3& p = points[s[i]];
    std::size_t k = values[s[i]];

    // DU WANG : dx(p,x) < dx(q,x)
    Metric m_p = mf->compute_metric(p);
    TPoint_3 tp(m_p.transform(p));
    double sq_d = csd(tp, m_p.transform(seeds[k]));

    // LABELLE SHEWCHUK : dp(p,x) < dq(q,x)
    TPoint_3 tp2 = (seeds_m[k]).transform(p);
    TPoint_3 ts = (seeds_m[k]).transform(seeds[k]);
//    double sq_d = csd(tp2, ts);

    if(sq_d > farthest_d)
    {
#pragma omp critical
{
      farthest_d = sq_d;
      farthest_x = p.x();
      farthest_y = p.y();
      farthest_z = p.z();
}
    }
  }
}

void get_ref_pt_from_triangle(const std::size_t i0, const std::size_t i1,
                              const std::size_t i2,
                              const std::vector<std::size_t>& values,
                              const std::vector<Point_3>& points)
{
  if(values[i0] == values[i1] || values[i1] == values[i2] || values[i0] == values[i2])
    return; // only keep faces

  Tri tr;
  tr[0] = values[i0]; tr[1] = values[i0]; tr[2] = values[i2];
  ugly_farthest_computation(tr, points, values);
}
#endif


void insert_simplex_if_colored(const std::size_t i0, const std::size_t i1,
                               const std::size_t i2, const std::size_t i3,
                               const std::vector<std::size_t>& values,
                               const std::vector<Point_3>& points)
{
  // quick filter
  if(values[i0] == values[i1] && values[i1] == values[i2] &&
     values[i2] == values[i3])
    return;

  if(values[i0] != values[i1] && values[i0] != values[i2] && values[i0] != values[i3] &&
     values[i1] != values[i2] && values[i1] != values[i3] &&
     values[i2] != values[i3])
  {
    Tet t;
    t[0] = values[i0]; t[1] = values[i1]; t[2] = values[i2]; t[3] = values[i3];
    std::sort(t.begin(), t.end());

    std::cout << "In the dual at : " << i0 << " " << i1 << " " << i2 << " " << i3 << std::endl;
    std::cout << "vals: " << t[0] << " " << t[1] << " " << t[2] << " " << t[3] << std::endl;
#pragma omp critical
    simplices.insert(t);
#ifdef TMP_REFINEMENT_UGLY_HACK
    t[0] = i0; t[1] = i1; t[2] = i2; t[3] = i3;
    ugly_farthest_computation(t, points, values);
  }
  else
  {
    get_ref_pt_from_triangle(i0, i1, i2, values, points);
    get_ref_pt_from_triangle(i0, i1, i3, values, points);
    get_ref_pt_from_triangle(i0, i2, i3, values, points);
    get_ref_pt_from_triangle(i1, i2, i3, values, points);
#endif
  }
}

void output_adapted_grid_points(std::ofstream& out_pts,
                                std::ofstream& out_bb,
                                std::size_t offset,
                                const std::vector<Point_3>& points,
                                const std::vector<std::size_t>& values)
{
  // we're adding new pts to existing points thus we add an offset
  CGAL_assertion(points.size() == values.size());

  std::size_t grid_nv = values.size();
  int counter = ++offset; // '++' because medit likes fortran
  for(std::size_t i=0; i!=grid_nv; ++i)
  {
    const Point_3& p = points[i];
    out_pts << p.x() << " " << p.y() << " " << p.z() << " " << ++counter << '\n';
    out_bb << values[i] << '\n';
  }
}

void output_adapted_grid_tets(std::ofstream& out_tets,
                              std::size_t offset,
                              const std::list<Cube>& final_cubes)
{
  // output the cube set (transformed into a triangulation)
  // the tets have LOCAL coordinates, so we need to add an offset

  ++offset; // '++' because medit likes fortran
  std::list<Cube>::const_iterator it = final_cubes.begin(), iend = final_cubes.end();
  for(; it!=iend; ++it)
  {
    const Cube& q = *it;
    out_tets << q.ids[0]+offset << " " << q.ids[1]+offset << " " << q.ids[2]+offset << " " << q.ids[5]+offset << " ";
    out_tets << "1" << '\n'; //has_simplex_n_colors(q.ids[0], q.ids[1], q.ids[2], q.ids[5], 4, values) << '\n';
    out_tets << q.ids[0]+offset << " " << q.ids[2]+offset << " " << q.ids[5]+offset << " " << q.ids[7]+offset << " ";
    out_tets << "1" << '\n'; //has_simplex_n_colors(q.ids[0], q.ids[2], q.ids[5], q.ids[7], 4, values) << '\n';
    out_tets << q.ids[0]+offset << " " << q.ids[2]+offset << " " << q.ids[3]+offset << " " << q.ids[7]+offset << " ";
    out_tets << "1" << '\n'; //has_simplex_n_colors(q.ids[0], q.ids[2], q.ids[3], q.ids[7], 4, values) << '\n';
    out_tets << q.ids[0]+offset << " " << q.ids[4]+offset << " " << q.ids[5]+offset << " " << q.ids[7]+offset << " ";
    out_tets << "1" << '\n'; //has_simplex_n_colors(q.ids[0], q.ids[4], q.ids[5], q.ids[7], 4, values) << '\n';
    out_tets << q.ids[2]+offset << " " << q.ids[5]+offset << " " << q.ids[6]+offset << " " << q.ids[7]+offset << " ";
    out_tets << "1" << '\n'; //has_simplex_n_colors(q.ids[2], q.ids[5], q.ids[6], q.ids[7], 4, values) << '\n';
  }
}

void output_adapted_grid_hex(std::ofstream& out_hex,
                             std::size_t offset,
                             const std::list<Cube>& final_cubes)
{
  // output the cube set

  ++offset; // '++' because medit likes fortran
  std::list<Cube>::const_iterator it = final_cubes.begin(), iend = final_cubes.end();
  for(; it!=iend; ++it)
  {
    const Cube& q = *it;
    out_hex << q.ids[0]+offset << " " << q.ids[1]+offset << " "
            << q.ids[2]+offset << " " << q.ids[3]+offset << " "
            << q.ids[4]+offset << " " << q.ids[5]+offset << " "
            << q.ids[6]+offset << " " << q.ids[7]+offset << " ";
    out_hex << "1" << '\n';
  }
}

void adapted_grid(const bool refine,
                  const bool output_grid = false)
{
#ifdef TMP_REFINEMENT_UGLY_HACK
  if(refine)
  {
    CGAL_precondition(farthest_x != 1e30 && farthest_y != 1e30 && farthest_z != 1e30);
    std::cout << "inserting: " << farthest_x << " " << farthest_y << " " << farthest_z << std::endl;
    vertices_nv = insert_new_seed(farthest_x, farthest_y, farthest_z);
    std::cout << "now: " << vertices_nv << " seeds" << std::endl;
    farthest_d = 0.; farthest_x = 1e30; farthest_y = 1e30; farthest_z = 1e30;
  }
#endif

  std::ofstream out("adapted_grid.mesh");
  std::ofstream out_bb("adapted_grid.bb");

  std::size_t grid_nv = 0, grid_ntet = 0, glob_offset = 0;

  // idea is to create some kind of octree and refine a cube
  // if its eight corners don't have all the same colors (and it's not too small)
  FT min_vol = grid_side*grid_side*grid_side*1e-7;
  FT max_vol = grid_side*grid_side*grid_side*1e-3;

  std::list<Cube> initial_cubes;
  std::vector<Point_3> points;
  std::vector<std::size_t> values;

  // create the first cube
  FT l = grid_side / 2.;
  FT x0 = center.x() - l, x1 = center.x() + l;
  FT y0 = center.x() - l, y1 = center.x() + l;
  FT z0 = center.x() - l, z1 = center.x() + l;

  Point_3 p0(x1,y1,z1), p1(x0,y1,z1), p2(x0,y0,z1), p3(x1,y0,z1);
  Point_3 p4(x1,y1,z0), p5(x0,y1,z0), p6(x0,y0,z0), p7(x1,y0,z0);
  points.push_back(p0); values.push_back(value_at_point(p0));
  points.push_back(p1); values.push_back(value_at_point(p1));
  points.push_back(p2); values.push_back(value_at_point(p2));
  points.push_back(p3); values.push_back(value_at_point(p3));
  points.push_back(p4); values.push_back(value_at_point(p4));
  points.push_back(p5); values.push_back(value_at_point(p5));
  points.push_back(p6); values.push_back(value_at_point(p6));
  points.push_back(p7); values.push_back(value_at_point(p7));

  Cube c0(0, 1, 2, 3, 4, 5, 6, 7);

// SPLIT
  split_cube(c0, initial_cubes, values, points);

  // remark: there is an obvious (small) redundancy in the points of the grid
  // coming from the borders of each region but it does not create any issue
  // (and it's simpler to do it that way...)

  omp_set_num_threads(8);
#pragma omp parallel shared(grid_nv, grid_ntet)
{
    int thread_ID = omp_get_thread_num();
    typename std::list<Cube>::iterator lcit = initial_cubes.begin();
    std::advance(lcit, thread_ID);

    std::list<Cube> local_final_cubes;
    std::vector<Point_3> local_points;
    std::vector<std::size_t> local_values;
    std::size_t local_offset; // 'local' ids need to be shifted for them to be 'global' in output

    for(std::size_t i=0; i<lcit->ids.size(); ++i)
    {
      local_points.push_back(points[lcit->ids[i]]);
      local_values.push_back(values[lcit->ids[i]]);
    }

    const Cube ci(0, 1, 2, 3, 4, 5, 6, 7); // cube n°i in the local coordinates
    std::list<Cube> local_cubes_to_test;
    local_cubes_to_test.push_back(ci);

    // refine the cube set
    while(!local_cubes_to_test.empty())
    {
      const Cube& q = local_cubes_to_test.front();

      if(local_points.size() > max_local_grid_nv || q.is_too_small(min_vol, local_points))
      {
        local_final_cubes.push_back(q);
      }
      else if(q.is_too_big(max_vol, points))
      {
        split_cube(q, local_cubes_to_test, local_values, local_points);
      }
#ifdef REFINE_NEAR_CENTERS_ONLY
      else if(q.number_of_colors(local_values) >= 4)
#else
      else if(!q.has_same_colors(local_values))
#endif
      {
        split_cube(q, local_cubes_to_test, local_values, local_points);
        if(local_points.size()%1000 == 0)
        {
          std::cout << local_points.size() << " points for thread: " << thread_ID
                    << " (" << ((FT) local_points.size()/(FT) max_local_grid_nv)*100. << "%)" << std::endl;
        }
      }
      else
        local_final_cubes.push_back(q);

      local_cubes_to_test.pop_front();
    }

    std::cout << "end of refinement for thread " << thread_ID << std::endl;

#ifdef SMART_DUAL
  for(std::list<Cube>::iterator it = local_final_cubes.begin();
      it != local_final_cubes.end(); ++it)
  {
    const Cube& q = *it;
    // split the cube in five tetrahedra and insert a tet in simplices if colored
    // THESE (CAN) COMPUTE THE POINT IN A DUAL THAT IS FARTHEST FROM ITS SEED
    insert_simplex_if_colored(q.ids[0], q.ids[1], q.ids[2], q.ids[5], local_values, local_points);
    insert_simplex_if_colored(q.ids[0], q.ids[2], q.ids[5], q.ids[7], local_values, local_points);
    insert_simplex_if_colored(q.ids[0], q.ids[2], q.ids[3], q.ids[7], local_values, local_points);
    insert_simplex_if_colored(q.ids[0], q.ids[4], q.ids[5], q.ids[7], local_values, local_points);
    insert_simplex_if_colored(q.ids[2], q.ids[5], q.ids[6], q.ids[7], local_values, local_points);
  }

  std::cout << "End of simplices computation for thread " << thread_ID << std::endl;
#endif

  // OUTPUT all the regions to a single file. From now on, it's more or less
  // sequential but we don't want to lose local data yet !
  if(output_grid)
  {
#pragma omp critical
{
    std::cout << "adding " << local_values.size() << " from thread " << thread_ID << std::endl;
    grid_nv += local_values.size();
    grid_ntet += 5*local_final_cubes.size();
}

#pragma omp barrier
#pragma omp single
{
    out << "MeshVersionFormatted 1" << '\n';
    out << "Dimension 3" << '\n';
    out << "Vertices" << '\n';
    out << grid_nv << '\n';

    out_bb << "3 1 " << grid_nv << " 2" << '\n';
}

#pragma omp critical // print the points
{
    local_offset = glob_offset;
//    std::cout << "loc off " << glob_offset << " @ " << thread_ID << std::endl;
    output_adapted_grid_points(out, out_bb, local_offset, local_points, local_values);

    //increment the offset for the next thread
    glob_offset += local_values.size();

    //clears to free memory asap (shouldn't really be needed)
    local_points.clear(); local_values.clear();
//    std::cout << "pts output done for thread " << thread_ID << std::endl;
}

#pragma omp barrier // waiting for all the threads to have written their points
#pragma omp single
{
    out << "Tetrahedra" << '\n';
    out << grid_ntet << '\n';
}

#pragma omp critical // print the tetrahedra
{
    output_adapted_grid_tets(out, local_offset, local_final_cubes);

    // locals shouldn't be needed anymore...
    local_points.clear(); local_final_cubes.clear(); local_values.clear();
//    std::cout << "Output done for thread " << thread_ID << std::endl;
}

/*
// medit is bugged so it segfaults when you use hexs...
// if you want to use it anyway, comment the tets (right above)
#pragma omp barrier // waiting for all the threads to have written their points
#pragma omp single
{
    out << "Hexahedra" << '\n';
    out << grid_ntet/5 << '\n';
}

#pragma omp critical // print the tetrahedra
{
    output_adapted_grid_hex(out, local_offset, local_final_cubes);

    // locals shouldn't be needed anymore...
    local_points.clear(); local_final_cubes.clear(); local_values.clear();
    std::cout << "Output done for thread " << thread_ID << std::endl;
}
*/

#pragma omp barrier // waiting for all the threads to have written their tets
#pragma omp single
{
    out << "End" << '\n';
    out_bb << "End" << '\n';
}
  }
} // end of parallel region

  std::cout << "after parallel : " << grid_nv << " pts" << std::endl;
}

void output_dual(bool compute_intersections = false)
{
  std::cout << "captured: " << simplices.size() << " simplices (seeds: "
            << seeds.size() << ")" << std::endl;

  std::ofstream outd("grid_dual.mesh");
  outd << "MeshVersionFormatted 1" << '\n';
  outd << "Dimension 3" << '\n';
  outd << "Vertices" << '\n';
  outd << vertices_nv << '\n';
  for(std::size_t i=0; i<vertices_nv; ++i)
    outd << seeds[i].x() << " " << seeds[i].y() << " " << seeds[i].z() << " " << i+1 << '\n';

  outd << "Tetrahedra" << '\n';
  outd << simplices.size() << '\n';
  for(typename std::set<Tet>::iterator it = simplices.begin();
                                       it != simplices.end(); ++it)
  {
    const Tet& t = *it;
    for(std::size_t i=0; i<t.size(); ++i)
      outd << t[i] + 1 << " ";
    outd << "1" << '\n';
  }

  // kinda bad but easier to compute self intersections
//  std::set<Edge> edges;
//  if(compute_intersections)
//    edges = compute_edges();

//  const std::set<Tri>& triangles = compute_triangles();
//  outd << "Triangles" << '\n';
//  outd << triangles.size() << '\n';
//  for(typename std::set<Tri>::iterator it = triangles.begin();
//                                       it != triangles.end(); ++it)
//  {
//    const Tri& tri = *it;
//    for(std::size_t i=0; i<tri.size(); ++i)
//      outd << tri[i]+1 << " ";
//    if(compute_intersections)
//      outd << is_triangle_intersected(tri, edges) << '\n';
//    else
//      outd << "1" << '\n';
//  }

  outd << "End" << '\n';
}

// -----------------------------------------------------------------------------
// MAIN FUNCTIONS
// -----------------------------------------------------------------------------
void initialize()
{
  mf = new Euclidean_metric_field<K>();
//  mf = new Custom_metric_field<K>();

  traits = new Traits();
  vertices_nv = build_seeds();
}

void build_grid(const bool refine = false)
{
  simplices.clear();
  adapted_grid(refine, true/*output grid*/);
  output_dual(true/*output dual*/);
}

int main(int, char**)
{
//  std::freopen("log.txt", "w", stdout);
  std::srand(0);

  initialize();
  build_grid();

  for(std::size_t i=0; i<n_refine; ++i)
  {
    std::cout << "refine: " << i << std::endl;
    simplices.clear();
    bool output = (i%100)==0 || i==(n_refine-1);

    adapted_grid(true/*refine*/, output);
    if(output)
      output_dual(true);
  }

  std::cout << "EoP" << std::endl;
}

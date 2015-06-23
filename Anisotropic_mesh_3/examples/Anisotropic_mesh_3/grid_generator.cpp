#include <CGAL/Starset.h>
#include <CGAL/Stretched_Delaunay_3.h>

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <Eigen/Dense>

#include <fstream>
#include <iostream>
#include <vector>

//#define SMART_DUAL
//#define REFINE_NEAR_CENTERS_ONLY

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef typename K::FT                                       FT;
typedef Stretched_Delaunay_3<K>                              Star;
typedef typename Star::Point_3                               Point_3;
typedef typename Star::TPoint_3                              TPoint_3;
typedef typename Star::Metric                                Metric;
typedef typename Star::Traits                                Traits;

typedef std::set<std::size_t> Simplex;

// pick (only) one to compute the dual
//#define APPROXIMATE_DUAL
#define SMART_DUAL

#define FILTER_SEEDS_OUTSIDE_GRID

// using a lot of global variables, it's ugly but passing them by function is tedious
std::size_t vertices_nv = 200;

Point_3 center(0.75,0.75,0.75);

const FT grid_side = 0.5;
FT offset_x = center.x() - grid_side/2.; // offset is the bottom left point
FT offset_y = center.y() - grid_side/2.;
FT offset_z = center.z() - grid_side/2.;

std::vector<Point_3> seeds;
std::vector<Metric> seeds_m;

std::set<Simplex> simplices;

template<typename Metric_field>
int build_seeds(const Metric_field& mf)
{
  std::ifstream in("bambimboum.mesh");
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

#ifdef FILTER_SEEDS_OUTSIDE_GRID
    if(r_x <= offset_x || r_x >= center.x() + grid_side/2 ||
       r_y <= offset_y || r_y >= center.y() + grid_side/2 ||
       r_z <= offset_z || r_z >= center.y() + grid_side/2)
    {
      std::cout << "filtered : " << r_x << " " << r_y << " " << r_z << std::endl;
      continue;
    }
#endif

    seeds.push_back(Point_3(r_x, r_y, r_z));
//    std::cout << "Seeds i: " << i << " " << seeds[i].x() << " ";
//    std::cout << seeds[i].y() << " " << seeds[i].z() << std::endl;

    seeds_m.push_back(mf->compute_metric(seeds[i]));

    if(seeds.size() == vertices_nv)
      break;
  }
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

template<typename MF>
FT value_at_point(const Point_3& p, const MF* mf, const Traits* traits)
{
  Metric m_p = mf->compute_metric(p);
  TPoint_3 tp(m_p.transform(p));

  typename Star::Traits::Compute_squared_distance_3 csd =
      traits->compute_squared_distance_3_object();

  FT min = 1e30;
  unsigned int min_id = 0;
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
        min_id = 1.01*vertices_nv; // drawing borders and center points
        break;
      }
    }
  }

  Comp c(sqds);
  std::sort(ids.begin(), ids.end(), c);

  // check if we are in the dual of a triangle (somewhat, not exactly a witness)
  for(int i=0; i<vertices_nv-1; ++i)
  {
//    std::cout << sqds[ids[i]] << " " << sqds[ids[i+1]] << std::endl;
    assert(sqds[ids[i+1]] >= sqds[ids[i]]);
  }

#ifdef APPROXIMATE_DUAL
  FT alpha = 0.1;
  if(std::abs(sqds[ids[0]] - sqds[ids[1]]) < alpha*sqds[ids[0]] &&
     std::abs(sqds[ids[0]] - sqds[ids[2]]) < alpha*sqds[ids[0]] &&
     std::abs(sqds[ids[0]] - sqds[ids[3]]) < alpha*sqds[ids[0]]) // in the dual
  {
    std::set<std::size_t> simplex;
    simplex.insert(ids[0]); simplex.insert(ids[1]);
    simplex.insert(ids[2]); simplex.insert(ids[3]);
    simplices.insert(simplex);
  }
#endif

  return min_id;
}

template<typename MF>
void full_grid(const MF* metric_field, const Traits* traits)
{
  double n = 1000.; // number of points per side
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

  for(unsigned int i=0; i<n; ++i)
  {
    for(unsigned int j=0; j<n; ++j)
    {
      for(unsigned int k=0; k<n; ++k)
      {
        Point_3 p(offset_x+j*step, offset_y+i*step, offset_z+k*step);
        out << p.x() << " " << p.y() << " " << p.z() << " " << ++counter << std::endl;

        FT ret = value_at_point(p, metric_field, traits);
        out_bb << ret << std::endl;
      }
    }
  }
  out << "Tetrahedra" << std::endl;
  out << (int) tet_n << std::endl;
  for(int i=1; i<(n-1); ++i)
  {
    for(int j=0; j<(n-1); ++j)
    {
      std::cout << "fixme" << std::endl;
      CGAL_assertion(false);
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

  std::vector<std::size_t> ids;

  bool is_too_small(FT min_vol, const std::vector<Point_3>& points) const
  {
    FT dx = points[0].x() - points[1].x();
    FT dy = points[0].y() - points[3].y();
    FT dz = points[0].z() - points[4].z();
//    if(dx*dy*dz < min_vol)
//      std::cout << "toosmall " << std::endl;
    return dx*dy*dz < min_vol;
  }

  bool is_too_big(FT max_vol, const std::vector<Point_3>& points) const
  {
    FT dx = points[0].x() - points[1].x();
    FT dy = points[0].y() - points[3].y();
    FT dz = points[0].z() - points[4].z();
    return dx*dy*dz > max_vol;
  }

  std::size_t number_of_colors(const std::vector<FT>& values) const
  {
    std::set<FT> vals;
    std::vector<std::size_t>::const_iterator it = ids.begin(), iend = ids.end();
    for(; it!=iend; ++it)
      vals.insert(values[*it]);
    return vals.size();
  }

  bool has_same_colors(const std::vector<FT>& values) const
  {
    FT v0 = values[ids.front()];
    std::vector<std::size_t>::const_iterator it = ids.begin()++, iend = ids.end();
    for(; it!=iend; ++it)
      if(v0 != values[*it])
        return false;
    return true;
  }

  Cube(const std::vector<std::size_t>& ids_) : ids(ids_) { }
  Cube(const std::size_t i0, const std::size_t i1, const std::size_t i2,
       const std::size_t i3, const std::size_t i4 ,const std::size_t i5,
       const std::size_t i6, const std::size_t i7) : ids()
  {
    ids.push_back(i0); ids.push_back(i1); ids.push_back(i2); ids.push_back(i3);
    ids.push_back(i4); ids.push_back(i5); ids.push_back(i6); ids.push_back(i7);
  }
};

template<typename MF>
void split_cube(const Cube q,
                std::list<Cube>& cubes_to_test,
                std::vector<FT>& values,
                std::vector<Point_3>& points,
                const MF* mf, const Traits* traits)
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
  FT x0 = points[q.ids[1]].x();
  FT x1 = points[q.ids[0]].x();
  FT y0 = points[q.ids[2]].y();
  FT y1 = points[q.ids[1]].y();
  FT z0 = points[q.ids[4]].z();
  FT z1 = points[q.ids[0]].z();

  FT xmid = (x0 + x1) / 2.;
  FT ymid = (y0 + y1) / 2.;
  FT zmid = (z0 + z1) / 2.;

  //new pts coordinates (19!):
  std::size_t i8 = points.size();
  Point_3 p8(xmid,y1,z1);
  points.push_back(p8);
  values.push_back(value_at_point(p8, mf, traits));

  std::size_t i9 = points.size();
  Point_3 p9(x1,ymid,z1);
  points.push_back(p9);
  values.push_back(value_at_point(p9, mf, traits));

  std::size_t i10 = points.size();
  Point_3 p10(xmid,ymid,z1);
  points.push_back(p10);
  values.push_back(value_at_point(p10, mf, traits));

  std::size_t i11 = points.size();
  Point_3 p11(x0,ymid,z1);
  points.push_back(p11);
  values.push_back(value_at_point(p11, mf, traits));

  std::size_t i12 = points.size();
  Point_3 p12(xmid,y0,z1);
  points.push_back(p12);
  values.push_back(value_at_point(p12, mf, traits));

  std::size_t i13 = points.size();
  Point_3 p13(x1,y1,zmid);
  points.push_back(p13);
  values.push_back(value_at_point(p13, mf, traits));

  std::size_t i14 = points.size();
  Point_3 p14(xmid,y1,zmid);
  points.push_back(p14);
  values.push_back(value_at_point(p14, mf, traits));

  std::size_t i15 = points.size();
  Point_3 p15(x0,y1,zmid);
  points.push_back(p15);
  values.push_back(value_at_point(p15, mf, traits));

  std::size_t i16 = points.size();
  Point_3 p16(x1,ymid,zmid);
  points.push_back(p16);
  values.push_back(value_at_point(p16, mf, traits));

  std::size_t i17 = points.size();
  Point_3 p17(xmid,ymid,zmid);
  points.push_back(p17);
  values.push_back(value_at_point(p17, mf, traits));

  std::size_t i18 = points.size();
  Point_3 p18(x0,ymid,zmid);
  points.push_back(p18);
  values.push_back(value_at_point(p18, mf, traits));

  std::size_t i19 = points.size();
  Point_3 p19(x1,y0,zmid);
  points.push_back(p19);
  values.push_back(value_at_point(p19, mf, traits));

  std::size_t i20 = points.size();
  Point_3 p20(xmid,y0,zmid);
  points.push_back(p20);
  values.push_back(value_at_point(p20, mf, traits));

  std::size_t i21 = points.size();
  Point_3 p21(x0,y0,zmid);
  points.push_back(p21);
  values.push_back(value_at_point(p21, mf, traits));

  std::size_t i22 = points.size();
  Point_3 p22(xmid,y1,z0);
  points.push_back(p22);
  values.push_back(value_at_point(p22, mf, traits));

  std::size_t i23 = points.size();
  Point_3 p23(x1,ymid,z0);
  points.push_back(p23);
  values.push_back(value_at_point(p23, mf, traits));

  std::size_t i24 = points.size();
  Point_3 p24(xmid,ymid,z0);
  points.push_back(p24);
  values.push_back(value_at_point(p24, mf, traits));

  std::size_t i25 = points.size();
  Point_3 p25(x0,ymid,z0);
  points.push_back(p25);
  values.push_back(value_at_point(p25, mf, traits));

  std::size_t i26 = points.size();
  Point_3 p26(xmid,y0,z0);
  points.push_back(p26);
  values.push_back(value_at_point(p26, mf, traits));

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

template<typename MF>
void smart_grid(const MF* mf, const Traits* traits)
{
  // idea is to create some kind of octree and refine a square
  // if its four corners don't have all the same colors (and it's not too small)
  FT min_vol = grid_side*grid_side*grid_side*1e-7;
  FT max_vol = grid_side*grid_side*grid_side*1e-3;

  std::list<Cube> cubes_to_test;
  std::list<Cube> final_cubes;
  std::vector<Point_3> points;
  std::vector<FT> values;

  // create the first cube
  FT l = grid_side / 2.;
  FT x0 = offset_x - l, x1 = offset_x + l;
  FT y0 = offset_y - l, y1 = offset_y + l;
  FT z0 = offset_z - l, z1 = offset_z + l;

  Point_3 p0(x1,y1,z1), p1(x0,y1,z1), p2(x0,y0,z1), p3(x1,y0,z1);
  Point_3 p4(x1,y1,z0), p5(x0,y1,z0), p6(x0,y0,z0), p7(x1,y0,z0);
  points.push_back(p0); points.push_back(p1); points.push_back(p2); points.push_back(p3);
  points.push_back(p4); points.push_back(p5); points.push_back(p6); points.push_back(p7);
  values.push_back(value_at_point(p0, mf, traits));
  values.push_back(value_at_point(p1, mf, traits));
  values.push_back(value_at_point(p2, mf, traits));
  values.push_back(value_at_point(p3, mf, traits));
  values.push_back(value_at_point(p4, mf, traits));
  values.push_back(value_at_point(p5, mf, traits));
  values.push_back(value_at_point(p6, mf, traits));
  values.push_back(value_at_point(p7, mf, traits));

  Cube c0(0,1,2,3,4,5,6,7);
  cubes_to_test.push_back(c0);

  // refine the cube set
  while(!cubes_to_test.empty())
  {
    const Cube& q = cubes_to_test.front();

    if(q.is_too_small(min_vol, points) || q.has_same_colors(values) || points.size() > 1e6)
      final_cubes.push_back(q);
#ifdef REFINE_NEAR_CENTERS_ONLY
    else if(q.number_of_colors(values) >= 4)
#else
    else if(!q.has_same_colors(values))
#endif
    {
      split_cube(q, cubes_to_test, values, points, mf, traits);
      std::cout << "Split. Now: " << points.size() << " points" << std::endl;
    }
    else
      final_cubes.push_back(q);

    cubes_to_test.pop_front();
  }

  // output the cube set (easily transformed into a triangulation) -------------
  CGAL_assertion(points.size() == values.size());
  std::size_t grid_side = values.size();
  std::cout << "check: " << points.size() << " " << values.size() << " vertices" << std::endl;

  std::ofstream out("smart_grid.mesh");
  out << "MeshVersionFormatted 1" << std::endl;
  out << "Dimension 3" << std::endl;
  out << "Vertices" << std::endl;
  out << grid_side << std::endl;

  std::ofstream out_bb("smart_grid.bb");
  out_bb << "3 1 " << grid_side << " 2" << std::endl;

  int counter = 0;
  for(int i=0; i!=grid_side; ++i)
  {
    const Point_3& p = points[i];
    out << p.x() << " " << p.y() << " " << p.z() << " " << ++counter << std::endl;
    out_bb << values[i] << std::endl;
  }
  out_bb << "End" << std::endl;

  out << "Tetrahedra" << std::endl;
  out << 5*final_cubes.size() << std::endl;
  std::list<Cube>::iterator it = final_cubes.begin(), iend = final_cubes.end();
  for(; it!=iend; ++it)
  {
    const Cube& q = *it;
    out << q.ids[0]+1 << " " << q.ids[1]+1 << " " << q.ids[2]+1 << " " << q.ids[5]+1 << " 1" << std::endl;
    out << q.ids[0]+1 << " " << q.ids[2]+1 << " " << q.ids[5]+1 << " " << q.ids[7]+1 << " 2" << std::endl;
    out << q.ids[0]+1 << " " << q.ids[2]+1 << " " << q.ids[3]+1 << " " << q.ids[7]+1 << " 3" << std::endl;
    out << q.ids[0]+1 << " " << q.ids[4]+1 << " " << q.ids[5]+1 << " " << q.ids[7]+1 << " 4" << std::endl;
    out << q.ids[2]+1 << " " << q.ids[5]+1 << " " << q.ids[6]+1 << " " << q.ids[7]+1 << " 5" << std::endl;
  }
  out << "End" << std::endl;

#ifdef SMART_DUAL
  // split the cube into tets. If a tet has 4 different values, the dual exists
  simplices.clear();
  for(it=final_cubes.begin(); it!=iend; ++it)
  {
    const Cube& q = *it;
// test all five tetrahedra
    std::set<std::size_t> vals;
    vals.insert(values[q.ids[0]]); vals.insert(values[q.ids[1]]);
    vals.insert(values[q.ids[2]]); vals.insert(values[q.ids[5]]);
    if(vals.size() == 4)
      simplices.insert(vals);
//------------
    vals.clear();
    vals.insert(values[q.ids[0]]); vals.insert(values[q.ids[2]]);
    vals.insert(values[q.ids[5]]); vals.insert(values[q.ids[7]]);
    if(vals.size() == 4)
      simplices.insert(vals);
//------------
    vals.clear();
    vals.insert(values[q.ids[0]]); vals.insert(values[q.ids[2]]);
    vals.insert(values[q.ids[3]]); vals.insert(values[q.ids[7]]);
    if(vals.size() == 4)
      simplices.insert(vals);
//------------
    vals.clear();
    vals.insert(values[q.ids[0]]); vals.insert(values[q.ids[4]]);
    vals.insert(values[q.ids[5]]); vals.insert(values[q.ids[7]]);
    if(vals.size() == 4)
      simplices.insert(vals);
//------------
    vals.clear();
    vals.insert(values[q.ids[2]]); vals.insert(values[q.ids[5]]);
    vals.insert(values[q.ids[6]]); vals.insert(values[q.ids[7]]);
    if(vals.size() == 4)
      simplices.insert(vals);
  }
#endif
}

void output_simplices()
{
  std::cout << "captured: " << simplices.size() << " simplices" << std::endl;

  std::ofstream outd("grid_dual.mesh");
  outd << "MeshVersionFormatted 1" << std::endl;
  outd << "Dimension 3" << std::endl;
  outd << "Vertices" << std::endl;
  outd << vertices_nv << std::endl;
  for(int i=0; i<vertices_nv; ++i)
    outd << seeds[i].x() << " " << seeds[i].y() << " " << seeds[i].z() << " " << i+1 << std::endl;

  outd << "Tetrahedra" << std::endl;
  outd << simplices.size() << std::endl;
  for(typename std::set<Simplex>::iterator it = simplices.begin();
                                           it != simplices.end(); ++it)
  {
    const Simplex& s = *it;
    typename Simplex::iterator sit = s.begin();
    typename Simplex::iterator siend = s.end();
    for(; sit!=siend; ++sit)
      outd << *sit+1 << " ";

    outd << "1" << std::endl;
  }
  outd << "End" << std::endl;
}

template<typename MF>
void draw(const MF* metric_field)
{
  Traits* traits = new Traits();
  vertices_nv = build_seeds(metric_field);

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

//  Euclidean_metric_field<K>* metric_field = new Euclidean_metric_field<K>();
  Custom_metric_field<K>* metric_field = new Custom_metric_field<K>();

  draw(metric_field);
}

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_sphere_2.h>
#include <CGAL/Projection_sphere_traits_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>
#include <CGAL/surface_neighbor_coordinates_3.h>

#include <boost/iterator/transform_iterator.hpp>

#include <chrono>
#include <deque>
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;

typedef K::FT                                                     FT;
typedef K::Point_3                                                Point_3;

typedef CGAL::Projection_sphere_traits_3<K>                       Projection_traits;
typedef CGAL::Delaunay_triangulation_sphere_2<Projection_traits>  Projected_DToS2;
typedef CGAL::Delaunay_triangulation_3<K>                         Delaunay_triangulation_3;

typedef std::pair<Point_3, FT>                                    Point_coordinate;
typedef std::vector<Point_coordinate>                             Point_coordinate_vector;

typedef std::deque<Point_3>                                       Point_container;

template <typename DTOS>
void output_triangulation_to_off(const DTOS& dtos,
                                 const char* filename)
{
  std::map<typename DTOS::Vertex_handle, int> ids;

  std::ofstream out(filename);
  out.precision(17);

  out << "OFF\n";
  out << dtos.number_of_vertices() << " " << dtos.number_of_faces() << " 0\n";

  int counter = 0;
  typename DTOS::All_vertices_iterator vit = dtos.vertices_begin(),
                                       vend = dtos.vertices_end();
  for(; vit!=vend; ++vit)
  {
    typename DTOS::Vertex_handle vh = vit;
    ids[vh] = counter++;
    out << vh->point() << "\n";
  }

  typename DTOS::All_faces_iterator fit = dtos.all_faces_begin(),
                                    fend = dtos.all_faces_end();
  for(; fit!=fend; ++fit)
  {
    out << "3 " << ids[fit->vertex(0)] << " "
                << ids[fit->vertex(1)] << " "
                << ids[fit->vertex(2)] << "\n" << std::endl;
  }
}

void test_interpolation(const Point_container& points,
                        const Point_container& queries)
{
  std::cout << " -- test interpolation with " << queries.size() << " queries and a range of points" << std::endl;

  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  for(const Point_3& query : queries)
  {
    K::Vector_3 normal(query - CGAL::ORIGIN);

    Point_coordinate_vector coords;
    CGAL::Triple<std::back_insert_iterator<Point_coordinate_vector>, FT, bool> result =
      CGAL::surface_neighbor_coordinates_3(points.begin(), points.end(), query, normal, std::back_inserter(coords), K());
  }

  std::chrono::steady_clock::time_point interpolation_time = std::chrono::steady_clock::now();
  std::cout << "Time elapsed (interpolation without triangulation): "
            << std::chrono::duration_cast<std::chrono::milliseconds>(interpolation_time - start_time).count()
            << "ms" << std::endl;
}

template <typename DT>
void test_interpolation(const DT& dt,
                        const Point_container& queries)
{
  std::cout << " -- test interpolation with " << queries.size() << " queries (" << typeid(DT).name() << ")" << std::endl;

  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  for(const Point_3& query : queries)
  {
    K::Vector_3 normal(query - CGAL::ORIGIN);

    Point_coordinate_vector coords;
    CGAL::Triple<std::back_insert_iterator<Point_coordinate_vector>, FT, bool> result =
      CGAL::surface_neighbor_coordinates_3(dt, query, normal, std::back_inserter(coords));

    std::set<Point_coordinate> sorted_coords(coords.begin(), coords.end());

//    std::cout << "coordinates of: " << query << std::endl;
//    for(const Point_coordinate& pc : sorted_coords)
//      std::cout << "pt: " << pc.first << " c: " << pc.second << std::endl;
  }

  std::chrono::steady_clock::time_point interpolation_time = std::chrono::steady_clock::now();
  std::cout << "Time elapsed (interpolation with triangulation): "
            << std::chrono::duration_cast<std::chrono::milliseconds>(interpolation_time - start_time).count()
            << "ms" << std::endl;
}

template <typename DT1, typename DT2>
void test_interpolation(const DT1& dt3, const DT2& dtos,
                        const Point_container& queries)
{
  std::cout << " -- test interpolation with " << queries.size() << " queries (both triangulations)" << std::endl;

  for(const Point_3& query : queries)
  {
    K::Vector_3 normal(query - CGAL::ORIGIN);

    Point_coordinate_vector coords_dt3;
    Point_coordinate_vector coords_dtos;

    CGAL::Triple<std::back_insert_iterator<Point_coordinate_vector>, FT, bool> result_dt3 =
      CGAL::surface_neighbor_coordinates_3(dt3, query, normal, std::back_inserter(coords_dt3));
    CGAL::Triple<std::back_insert_iterator<Point_coordinate_vector>, FT, bool> result_dtos =
      CGAL::surface_neighbor_coordinates_3(dtos, query, normal, std::back_inserter(coords_dtos));

    std::set<Point_coordinate> sorted_coords_dt3(coords_dt3.begin(), coords_dt3.end());
    std::set<Point_coordinate> sorted_coords_dtos(coords_dtos.begin(), coords_dtos.end());

    if(!result_dtos.third || !result_dt3.third)
    {
      std::cout << "baaaah" << std::endl;
      std::exit(0);
    }

    bool same_result = true;

    if(sorted_coords_dt3.size() != sorted_coords_dtos.size())
      same_result = false;

    if(same_result)
    {
      const FT eps = std::numeric_limits<FT>::epsilon();
      typename std::set<Point_coordinate>::iterator sorted_coords_dt3_it = sorted_coords_dt3.begin();
      typename std::set<Point_coordinate>::iterator sorted_coords_dtos_it = sorted_coords_dtos.begin();

      for(std::size_t i=0, n=sorted_coords_dt3.size(); i<n; ++i)
      {
        if(CGAL::abs(sorted_coords_dt3_it->first.x() - sorted_coords_dtos_it->first.x()) > eps ||
           CGAL::abs(sorted_coords_dt3_it->first.y() - sorted_coords_dtos_it->first.y()) > eps ||
           CGAL::abs(sorted_coords_dt3_it->first.z() - sorted_coords_dtos_it->first.z()) > eps ||
           CGAL::abs(sorted_coords_dt3_it->second - sorted_coords_dtos_it->second) > eps)
        {
          same_result = false;
          break;
        }

        ++sorted_coords_dt3_it;
        ++sorted_coords_dtos_it;
      }
    }

    if(!same_result)
    {
      std::cerr << "DIFFERENT RESULT: " << std::endl;
      std::cout << "coordinates of: " << query << std::endl;

      std::cout << "DT3:" << std::endl;
      for(const Point_coordinate& pc : sorted_coords_dt3)
        std::cout << "pt: " << pc.first << " c: " << pc.second << std::endl;

      std::cout << "DTOS:" << std::endl;
      for(const Point_coordinate& pc : sorted_coords_dtos)
        std::cout << "pt: " << pc.first << " c: " << pc.second << std::endl;

      std::exit(0);
    }
  }

  std::cout << "Clear!" << std::endl;
}

Projected_DToS2 construct_DTOS(const Point_container& points)
{
  Projection_traits traits(Point_3(0,0,0));
  Projected_DToS2 dtos(traits);

  Projection_traits::Construct_projected_point_3 cst =
      traits.construct_projected_point_3_object();

  std::cout << " -- Construct DTOS" << std::endl;
  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  dtos.insert(boost::make_transform_iterator(points.begin(), cst),
              boost::make_transform_iterator(points.end(), cst));

  std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();
  std::cout << "Time elapsed (dtos): "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()
            << "ms" << std::endl;

  std::cout << dtos.number_of_vertices() << " vertices" << std::endl;
  std::cout << dtos.number_of_faces() << " faces" << std::endl;

  output_triangulation_to_off(dtos, "triangulation.off");

  return dtos;
}

Delaunay_triangulation_3 construct_DT3(const Point_container& points)
{
  std::cout << " -- Construct DT3" << std::endl;
  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  Delaunay_triangulation_3 dt3(points.begin(), points.end());

  std::cout << dt3.number_of_vertices() << " vertices" << std::endl;
  std::cout << dt3.number_of_finite_cells() << " cells" << std::endl;

  std::chrono::steady_clock::time_point construction_time = std::chrono::steady_clock::now();
  std::cout << "Time elapsed (dt3): "
            << std::chrono::duration_cast<std::chrono::milliseconds>(construction_time - start_time).count()
            << "ms" << std::endl;

  return dt3;
}

int main()
{
  std::cout << std::setprecision(17) << std::fixed;

// #define USE_DT3
#define USE_DTOS


  int n = 1000000;
  Point_container points;

  std::cout << " -- Generate " << n << " random points on a sphere." << std::endl;
  CGAL::Random rnd = CGAL::get_default_random();
  CGAL::Random_points_on_sphere_3<Point_3> g(1, rnd);
  CGAL::cpp11::copy_n(g, n, std::back_inserter(points));

#ifdef USE_DT3
  const Delaunay_triangulation_3& dt3 = construct_DT3(points);
#endif

#ifdef USE_DTOS
  const Projected_DToS2& dtos = construct_DTOS(points);
#endif

  int qn = 100000;
  Point_container queries;
  CGAL::cpp11::copy_n(g, qn, std::back_inserter(queries));

//  std::cout << "Timings..." << std::endl;
//  test_interpolation(points, queries);

#ifdef USE_DT3
  test_interpolation(dt3, queries);
#endif

#ifdef USE_DTOS
  test_interpolation(dtos, queries);
#endif

#if defined(USE_DT3) && defined(USE_DTOS)
  std::cout << "Comparing both..." << std::endl;
  test_interpolation(dt3, dtos, queries);
#endif

  assert(dtos.is_valid());
}

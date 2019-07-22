#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/point_generators_2.h>

#include <chrono>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::FT                                               Coord_type;
typedef K::Point_2                                          Point;
typedef CGAL::Delaunay_triangulation_2<K>                   Delaunay_triangulation;

// Resulting points-coordinates pairs will be stored in an object of this type
typedef std::vector<std::pair<Point, Coord_type> >          Point_coordinate_vector;

int main()
{
  int n = 1000000;
  std::vector<Point> points;

  points.push_back(Point(-2, -2));
  points.push_back(Point(-2,  2));
  points.push_back(Point( 2, -2));
  points.push_back(Point( 2,  2));

  CGAL::Random rnd = CGAL::get_default_random();
  CGAL::Random_points_in_square_2<Point> g(1, rnd);
  std::copy_n(g, n, std::back_inserter(points));

  std::chrono::steady_clock::time_point construction_start_time = std::chrono::steady_clock::now();

  Delaunay_triangulation dt(points.begin(), points.end());

  std::chrono::steady_clock::time_point construction_end_time = std::chrono::steady_clock::now();
  std::cout << "Time elapsed (dt2): "
            << std::chrono::duration_cast<std::chrono::milliseconds>(construction_end_time - construction_start_time).count()
            << "ms" << std::endl;

  int qn = 100000;
  std::vector<Point> queries;
  std::copy_n(g, qn, std::back_inserter(queries));

  std::chrono::steady_clock::time_point nn_start_time = std::chrono::steady_clock::now();

  for(const Point& query : queries)
  {
    Point_coordinate_vector coords;
    CGAL::Triple<std::back_insert_iterator<Point_coordinate_vector>, K::FT, bool> result =
        CGAL::natural_neighbor_coordinates_2(dt, query, std::back_inserter(coords));

//    if(!result.third)
//    {
//      std::cout << "The coordinate computation was not successful." << std::endl;
//      std::cout << "The point (" << p << ") lies outside the convex hull." << std::endl;
//    }

//    K::FT norm = result.second;
//    std::cout << "Coordinate computation successful." << std::endl;
//    std::cout << "Normalization factor: " << norm << std::endl;

//    std::cout << "Coordinates for point: (" << p << ") are the following: " << std::endl;
//    for(std::size_t i=0; i<coords.size(); ++i)
//      std::cout << "  Point: (" << coords[i].first << ") coeff: " << coords[i].second << std::endl;
  }

  std::chrono::steady_clock::time_point nn_end_time = std::chrono::steady_clock::now();
  std::cout << "Time elapsed (nn): "
            << std::chrono::duration_cast<std::chrono::milliseconds>(nn_end_time - nn_start_time).count()
            << "ms" << std::endl;

  return EXIT_SUCCESS;
}

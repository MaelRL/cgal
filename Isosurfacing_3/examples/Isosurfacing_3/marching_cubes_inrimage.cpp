#include <CGAL/Simple_cartesian.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/Explicit_cartesian_grid_domain.h>
#include <CGAL/Isosurfacing_3/Marching_cubes_3.h>

#include <CGAL/boost/graph/IO/OFF.h>

#include <memory>
#include <vector>

using Kernel = CGAL::Simple_cartesian<double>;
using Point = typename Kernel::Point_3;
using Grid = CGAL::Cartesian_grid_3<Kernel>;

using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;

int main(int, char**)
{
  const std::string fname = CGAL::data_file_path("images/skull_2.9.inr");

  // load volumetric image from a file
  CGAL::Image_3 image;
  if(!image.read(fname))
  {
    std::cerr << "Error: Cannot read image file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  // convert image to a Cartesian grid
  std::shared_ptr<Grid> grid = std::make_shared<Grid>(image);

  // create a domain from the grid
  auto domain = CGAL::Isosurfacing::create_explicit_cartesian_grid_domain<Kernel>(grid);

  // prepare collections for the output indexed mesh
  Point_range points;
  Polygon_range polygons;

  // execute marching cubes with an isovalue of 2.9
  CGAL::Isosurfacing::marching_cubes(domain, 2.9, points, polygons);

  // save output indexed mesh to a file, in the OFF format
  CGAL::IO::write_OFF("result.off", points, polygons);

  return EXIT_SUCCESS;
}

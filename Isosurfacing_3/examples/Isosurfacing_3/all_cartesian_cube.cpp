#include <CGAL/Cartesian_grid_3.h>
#include <CGAL/Cartesian_grid_domain.h>
#include <CGAL/Dual_contouring_3.h>
#include <CGAL/Marching_cubes_3.h>
#include <CGAL/TC_marching_cubes_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/IO/OFF.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef typename Kernel::FT FT;
typedef typename Kernel::Point_3 Point;
typedef typename Kernel::Vector_3 Vector;

typedef CGAL::Cartesian_grid_3<Kernel> Grid;

typedef std::vector<Point> Point_range;
typedef std::vector<std::vector<std::size_t>> Polygon_range;

int main() {
    // create a cartesian grid with 100^3 grid points and the bounding box [-1, 1]^3
    Grid grid(100, 100, 100, {-1, -1, -1, 1, 1, 1});

    // calculate the value at all grid points
    for (std::size_t x = 0; x < grid.xdim(); x++) {
        for (std::size_t y = 0; y < grid.ydim(); y++) {
            for (std::size_t z = 0; z < grid.zdim(); z++) {

                const FT pos_x = x * grid.get_spacing()[0] + grid.get_bbox().xmin();
                const FT pos_y = y * grid.get_spacing()[1] + grid.get_bbox().ymin();
                const FT pos_z = z * grid.get_spacing()[2] + grid.get_bbox().zmin();

                // manhattan distance to the origin
                grid.value(x, y, z) = std::max({std::abs(pos_x), std::abs(pos_y), std::abs(pos_z)});

                // the normal depends on the side of the cube
                if (grid.value(x, y, z) == std::abs(pos_x)) {
                    grid.gradient(x, y, z) = Vector(std::copysign(1.0, pos_x), 0, 0);
                } else if (grid.value(x, y, z) == std::abs(pos_y)) {
                    grid.gradient(x, y, z) = Vector(0, std::copysign(1.0, pos_y), 0);
                } else if (grid.value(x, y, z) == std::abs(pos_z)) {
                    grid.gradient(x, y, z) = Vector(0, 0, std::copysign(1.0, pos_z));
                }
            }
        }
    }

    // create a domain from the grid
    CGAL::Isosurfacing::Cartesian_grid_domain<Kernel> domain(grid);

    // prepare collections for the results
    Point_range points_mc, points_tmc, points_dc;
    Polygon_range polygons_mc, polygons_tmc, polygons_dc;

    // execute marching cubes, topologically correct marching cubes and dual contouring with an isovalue of 0.8
    CGAL::Isosurfacing::make_triangle_mesh_using_marching_cubes(domain, 0.8, points_mc, polygons_mc);
    CGAL::Isosurfacing::make_triangle_mesh_using_tmc(domain, 0.8, points_tmc, polygons_tmc);
    CGAL::Isosurfacing::make_quad_mesh_using_dual_contouring(domain, 0.8, points_dc, polygons_dc);

    // save the results in the OFF format
    CGAL::IO::write_OFF("result_mc.off", points_mc, polygons_mc);
    CGAL::IO::write_OFF("result_tmc.off", points_tmc, polygons_tmc);
    CGAL::IO::write_OFF("result_dc.off", points_dc, polygons_dc);
}

#include "CLI11.hpp"
#include "Cartesian_grid_3.h"
#include "Cartesian_grid_oracle.h"
#include "Dual_contouring_3.h"
#include "Marching_cubes_3.h"

#include "types.h"
#include <CGAL/Timer.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <fstream>
#include <iostream>
#include <math.h>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef Kernel::Triangle_3 Triangle;
typedef std::vector<Triangle>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<Kernel, Iterator> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

inline Kernel::FT distance_to_mesh(const Tree& tree, const Point_3& p) {
	const Point_3& x = tree.closest_point(p);
	return std::sqrt((p - x).squared_length());
}

int main(int argc, char** argv) {
	CLI::App app;
	app.description("Computes the offset to an input mesh using Marching Cubes, Dual Contouring and Alpha Wrapping. To see the options use --help.");

	FT offset_value = 0.1;
	int n_voxel_points = 30;
	std::string input_name;
	bool use_bbox = false;
	app.add_option("--offset,-v", offset_value, "Offset value")->required()->check(CLI::NonNegativeNumber);
	app.add_option("--np,-n", n_voxel_points, "Number of grid points in every dimension")->required()->check(CLI::PositiveNumber);
	app.add_option("--input,-i", input_name, "Input mesh, OFF file format")->required()->check(CLI::ExistingFile);
	app.add_flag("--bbox", use_bbox, "Ensures that vertices stay inside their voxel during dual contouring")->default_val(true);

	CLI11_PARSE(app, argc, argv);

	// Mesh mesh_input;
	Point_range points_in;
	Polygon_range polygons_in;
	std::cout << "Reading input mesh...";
	if (!CGAL::IO::read_OFF(input_name, points_in, polygons_in)) {
		std::cout << "failed" << std::endl;
		exit(-1);
	}
	std::cout << "done (" << points_in.size() << " vertices, " << polygons_in.size() << " triangles)" << std::endl;

	// fill triangle soup
	std::vector<Triangle> triangles;
	triangles.reserve(polygons_in.size());
	for (const auto& polygon : polygons_in) {
		if (polygon.size() < 3) {
			continue;
		}
		for (int i = 1; i < polygon.size() - 1; ++i) {
			triangles.push_back(Triangle(points_in[polygon[0]], points_in[polygon[i]], points_in[polygon[i + 1]]));
		}
	}

	// compute bounding box
	std::cout << "Compute loose bounding box...";
	CGAL::Bbox_3 aabb_grid;
	for (const auto& t : triangles) {
		aabb_grid += t.bbox();
	}

	// compute min span
	const FT dx = aabb_grid.x_span();
	const FT dy = aabb_grid.y_span();
	const FT dz = aabb_grid.z_span();
	const FT min_span = std::min(dx, std::min(dy, dz));

	Vector_3 aabb_increase_vec = Vector_3(offset_value * 1.1, offset_value * 1.1, offset_value * 1.1);
	aabb_grid += (Point_3(aabb_grid.xmax(), aabb_grid.ymax(), aabb_grid.zmax()) + aabb_increase_vec).bbox();
	aabb_grid += (Point_3(aabb_grid.xmin(), aabb_grid.ymin(), aabb_grid.zmin()) - aabb_increase_vec).bbox();
	std::cout << "done" << std::endl;

	// construct AABB tree from the soup
	std::cout << "Build AABB tree...";
	Tree tree(triangles.begin(), triangles.end());
	std::cout << "done" << std::endl;

	Grid grid(n_voxel_points, n_voxel_points, n_voxel_points, aabb_grid);

	CGAL::Cartesian_grid_oracle<Kernel> grid_oracle(grid);

	std::cout << "Init grid...";
	const std::size_t size_k = grid.xdim();
	const std::size_t size_j = grid.ydim();
	const std::size_t size_i = grid.zdim();
#pragma omp parallel for
	for (int z = 0; z < grid.zdim(); z++) {
		for (int y = 0; y < grid.ydim(); y++) {
			for (int x = 0; x < grid.xdim(); x++) {
				const auto& q = grid_oracle.position(x, y, z);
				grid.value(x, y, z) = distance_to_mesh(tree, q);
			}
		}
	}
	std::cout << "done" << std::endl;


	// Marching Cubes
	std::string mc_output_name("out_mc.off");
	Point_range mc_points_out;
	Polygon_range mc_polygons_out;

	std::cout << "Run Marching Cubes...";
	CGAL::make_triangle_mesh_using_marching_cubes(grid_oracle, offset_value, mc_points_out, mc_polygons_out);
	std::cout << "done (" << mc_points_out.size() << " vertices)" << std::endl;

	std::cout << "Save output mesh...";
	CGAL::IO::write_OFF(mc_output_name, mc_points_out, mc_polygons_out);
	std::cout << "done" << std::endl;

	// Dual contouring
	std::string dc_output_name("out_dc.off");
	Point_range dc_points_out;
	Polygon_range dc_polygons_out;

	std::cout << "Run Dual Contouring...";
	CGAL::make_quad_mesh_using_dual_contouring(grid_oracle, offset_value, dc_points_out, dc_polygons_out, use_bbox);
	std::cout << "done (" << dc_points_out.size() << " vertices)" << std::endl;

	std::cout << "Save output mesh...";
	CGAL::IO::write_OFF(dc_output_name, dc_points_out, dc_polygons_out);
	std::cout << "done" << std::endl;

	// Alpha wrap
	// compute alpha value as twice min_span / #grid points
	const FT alpha_value = 2.0 * min_span / (FT)n_voxel_points;
	std::cout << "Alpha: " << alpha_value << std::endl;
	std::cout << "Offset: " << offset_value << std::endl;

	std::cout << "Run alpha wrap...";
	Mesh out_wrap;
	CGAL::Timer timer;
	timer.start();
	CGAL::alpha_wrap_3(points_in, polygons_in, alpha_value, offset_value, out_wrap);
	timer.stop();
	std::cout << "done (" << num_vertices(out_wrap) << " vertices, " << timer.time() << " s)" << std::endl;

	std::cout << "Save output mesh...";
	std::string wrap_output_name("out_wrap.off");
	CGAL::IO::write_polygon_mesh(wrap_output_name, out_wrap, CGAL::parameters::stream_precision(17));
	std::cout << "done" << std::endl;
}

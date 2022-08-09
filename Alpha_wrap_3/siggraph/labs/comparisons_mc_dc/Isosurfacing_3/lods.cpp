#include "CLI11.hpp"

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


int main(int argc, char** argv) {
	CLI::App app;
	app.description("Computes levels of details via Alpha Wrapping. To see the options use --help.");

	std::string input_name;
	app.add_option("--input,-i", input_name, "Input mesh, OFF file format")->required()->check(CLI::ExistingFile);

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
	CGAL::Bbox_3 bbox;
	for (const auto& t : triangles) {
		bbox += t.bbox();
	}

	// compute min span
	const FT dx = bbox.x_span();
	const FT dy = bbox.y_span();
	const FT dz = bbox.z_span();
	const FT min_span = std::min(dx, std::min(dy, dz));
	std::cout << "done" << std::endl;

	// Alpha wrap
	std::vector<double> alphas =  {5.0,  10.0, 20.0,  50.0, 100.0, 200.0, 400.0};
	std::vector<double> offsets = {10.0, 20.0, 40.0, 100.0, 200.0, 400.0, 800.0};

	for (int i = 0; i < alphas.size(); i++)
		{
			const FT alpha_value = min_span / (FT)alphas[i];
			const FT offset_value = min_span / (FT)offsets[i];

			std::cout << "Alpha wrapping (" << alphas[i] << ", " << offsets[i] << ")...";
			Mesh out_wrap;
			CGAL::Timer timer;
			timer.start();
			CGAL::alpha_wrap_3(points_in, polygons_in, alpha_value, offset_value, out_wrap);
			timer.stop();
			std::cout << "done (" << num_vertices(out_wrap) << " vertices, " << timer.time() << " s)" << std::endl;

			std::cout << "Save output mesh...";
			std::string name(input_name);
			name = name.substr(name.find_last_of("/") + 1, name.length() - 1);
			name = name.substr(0, name.find_last_of("."));
			std::string output_name = name
				+ "_" + std::to_string(static_cast<int>(alphas[i]))
				+ "_" + std::to_string(static_cast<int>(offsets[i])) + ".off";
			CGAL::IO::write_polygon_mesh(output_name, out_wrap, CGAL::parameters::stream_precision(17));
			std::cout << "done" << std::endl;

		}
}

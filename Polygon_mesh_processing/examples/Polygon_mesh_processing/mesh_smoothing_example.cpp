#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/mesh_smoothing.h>
#include <CGAL/Polygon_mesh_processing/shape_smoothing.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[]){

    namespace PMP = CGAL::Polygon_mesh_processing;
    const char* filename = argc > 1 ? argv[1] : "data/eight.off";
    std::ifstream input(filename);

    Mesh mesh;
    if (!input || !(input >> mesh) || mesh.is_empty()) {
      std::cerr << "Not a valid .off file." << std::endl;
      return 1;
    }

    const unsigned int repeat = 2;
    const unsigned int nb_iterations = 5;
    const double gradient_descent_precision = 1e-6;

    for(unsigned int t = 0 ; t < repeat; ++t)
    {
      PMP::smooth_angles(mesh,
        PMP::parameters::number_of_iterations(nb_iterations));

      PMP::smooth_areas(mesh,
        PMP::parameters::gradient_descent_precision(gradient_descent_precision));

      PMP::smooth_angles(mesh,
        PMP::parameters::number_of_iterations(nb_iterations));
    }

    std::ofstream output("data/eight_smoothed.off");
    output << mesh;
    output.close();


    return 0;
}

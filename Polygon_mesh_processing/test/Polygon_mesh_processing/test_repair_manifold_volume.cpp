#define CGAL_PMP_REPAIR_MANIFOLD_VOLUME_DEBUG

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/repair_manifold_volume.h>

#include <CGAL/IO/polygon_soup_io.h>

#include <fstream>

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using FT = K::FT;
using Point_3 = K::Point_3;

using Points = std::vector<Point_3>;
using Face = std::array<std::size_t, 3>;
using Faces = std::vector<Face>;

using Mesh = CGAL::Surface_mesh<Point_3>;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/blobby.off");

  Points points;
  Faces faces;
  if(!CGAL::IO::read_polygon_soup(filename, points, faces))
  {
    std::cerr << "Input mesh is not a valid file." << std::endl;
    return EXIT_FAILURE;
  }

  Mesh mesh;
  PMP::experimental::repair_manifold_volume_mini(points, faces, mesh,
                                                 CGAL::parameters::verbose(true),
                                                 CGAL::parameters::verbose(true));

  std::cout << num_vertices(mesh) << " nv " << num_faces(mesh) << " nf" << std::endl;

  CGAL::IO::write_polygon_mesh("mesh_fixed.off", mesh, CGAL::parameters::stream_precision(17));

  return EXIT_SUCCESS;
}

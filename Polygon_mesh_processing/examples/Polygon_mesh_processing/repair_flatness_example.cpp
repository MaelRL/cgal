#define CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE
#define CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
#define CGAL_PMP_STITCHING_DEBUG
#define CGAL_PMP_SMOOTHING_VERBOSE

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/IO/STL_reader.h>
#include <CGAL/IO/PLY_reader.h>

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/repair_flatness.h>
#include <CGAL/Polygon_mesh_processing/smooth_mesh.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                          Kernel;
typedef Kernel::FT                                              FT;
typedef Kernel::Point_3                                         Point_3;
typedef CGAL::Surface_mesh<Point_3>                             Surface_mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

template <typename K, typename Mesh>
void read_mesh(const char* filename, Mesh& sm)
{
  typedef typename K::Point_3                                   Point;

  std::ifstream in(filename, std::ios::binary);
  if(!in.good())
  {
    std::cerr << "Error: can't read file: " << filename << std::endl;
    std::exit(1);
  }

  std::string fn(filename);
  if(fn.substr(fn.find_last_of(".") + 1) == "stl")
  {
    std::vector<Point> points;
    std::vector<std::vector<int> > faces;
    if(!CGAL::read_STL(in, points, faces))
    {
      std::cerr << "Error: cannot open STL mesh\n";
      return;
    }

    std::cout << "Cleaning polygon soup..." << std::endl;
    PMP::repair_polygon_soup(points, faces);

    if(!CGAL::Polygon_mesh_processing::orient_polygon_soup(points, faces))
      std::cerr << "W: File does not describe a polygon mesh" << std::endl;

    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, faces, sm);
  }
  else if(fn.substr(fn.find_last_of(".") + 1) == "off")
  {
    if(!in || !(in >> sm))
    {
      std::cerr << "Error: cannot OFF open mesh\n";
      return;
    }
  }
  else if(fn.substr(fn.find_last_of(".") + 1) == "ply")
  {
    std::vector<Point> points;
    std::vector<std::vector<std::size_t> > polygons;
    std::vector<CGAL::Color> fcolors;
    std::vector<CGAL::Color> vcolors;

    if(!(CGAL::read_PLY(in, points, polygons, fcolors, vcolors)))
    {
      std::cerr << "Error: cannot open PLY mesh\n";
      return;
    }

    std::cout << "Cleaning polygon soup..." << std::endl;
    PMP::repair_polygon_soup(points, polygons);

    if(!CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons))
      std::cerr << "W: File does not describe a polygon mesh" << std::endl;

    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, sm);
  }
  else
  {
    std::cerr << "Unknown file type" << std::endl;
    return;
  }
}

template <typename TriangleMesh>
bool repair_local_self_intersection_with_smoothing(TriangleMesh& tmesh)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor       edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor       face_descriptor;
  typedef std::pair<face_descriptor, face_descriptor>                       Face_pair;

  std::vector<Face_pair> self_inter;

  PMP::self_intersections(tmesh, std::back_inserter(self_inter));

  // Here would be a good place to check that those are indeed local self-intersections...

  std::set<face_descriptor> faces;
  for(const Face_pair& f : self_inter)
  {
    faces.insert(f.first);
    faces.insert(f.second);
  }

  typedef typename boost::property_map<TriangleMesh, CGAL::edge_is_feature_t>::type  EIFMap;
  EIFMap eif = get(CGAL::edge_is_feature, tmesh);

  std::map<edge_descriptor, bool> is_border_of_selection;
  for(face_descriptor f : faces)
  {
    for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(f, tmesh), tmesh))
    {
      // meet it once, switch to 'true'; meet it twice, switch back to 'false'
      edge_descriptor e = edge(h, tmesh);
      is_border_of_selection[e] = !(is_border_of_selection[e]);
    }
  }

  FT cos_angle(std::cos(60. * CGAL_PI / 180.));

  for(const auto& ep : is_border_of_selection)
  {
    if(ep.second /*|| PMP::internal::is_sharp<TriangleMesh, Kernel>(tmesh, halfedge(ep.first, tmesh), cos_angle)*/)
      put(eif, ep.first, true);
  }

  std::ofstream out_sel("results/constrained_edges.selection.txt");
  out_sel << std::endl << std::endl; // edge selection is on the third line

  int counter = 0;
  std::map<vertex_descriptor, int> vid;
  for(vertex_descriptor v : vertices(tmesh))
    vid[v] = counter++;

  counter = 0;
  for(edge_descriptor e : edges(tmesh))
  {
    if(get(eif, e))
    {
      std::cout << e << " is constrained" << std::endl;
      ++counter;
      out_sel << vid[source(e, tmesh)] << " " << vid[target(e, tmesh)] << " ";
    }
  }
  std::cout << "total : " << counter << " out of " << num_edges(tmesh) << std::endl;

  PMP::smooth_mesh(faces, tmesh, CGAL::parameters::edge_is_constrained_map(eif)
                                                  .number_of_iterations(100)
                                                  .use_safety_constraints(false));

  return !PMP::does_self_intersect(tmesh);
}

int main(int argc, char** argv)
{
  Surface_mesh surface_mesh;

  const char* filename = argv[1];
  read_mesh<Kernel>(filename, surface_mesh);

  if(!CGAL::is_triangle_mesh(surface_mesh))
  {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Input mesh: " << num_vertices(surface_mesh) << " nv "
                              << num_edges(surface_mesh) << " ne "
                              << num_faces(surface_mesh) << " nf" << std::endl;


//  PMP::remove_degenerate_faces(surface_mesh);
//  PMP::repair_flatness(surface_mesh);

//  std::cout << "self inters (pre): " << PMP::does_self_intersect(surface_mesh) << std::endl;
//  PMP::remove_self_intersections(surface_mesh);
//  std::cout << "self inters (post): " << PMP::does_self_intersect(surface_mesh) << std::endl;

  repair_local_self_intersection_with_smoothing(surface_mesh);

  std::cout << "Done!" << std::endl;

  std::cout << "self inters (post post): " << PMP::does_self_intersect(surface_mesh) << std::endl;
  std::ofstream("results/final.off") << std::setprecision(17) << surface_mesh;

  return EXIT_SUCCESS;
}

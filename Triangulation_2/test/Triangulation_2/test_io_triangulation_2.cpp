// #define CGAL_T2_BUILD_FROM_RANGES_DEBUG

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>

#include <CGAL/Real_timer.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/point_generators_2.h>

#include <fstream>
#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel               K;
typedef K::Point_2                                                        Point_2;

//typedef CGAL::Triangulation_2<K>                                          Tr;
typedef CGAL::Delaunay_triangulation_2<K>                                 Tr;

typedef CGAL::Triangulation_vertex_base_2<K>                              Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K>                    Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                       Tds;
typedef CGAL::Exact_predicates_tag                                        Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds, Itag>          CDT;

typedef CGAL::Polygon_2<K>                                                Polygon_2;

typedef CGAL::Random_points_in_square_2<Point_2>                          Generator;

void constrained_random_test(const std::size_t n)
{
  std::ifstream in("/home/mrouxell/git/CGAL/cgal/GraphicsView/demo/Triangulation_2/data/norway.edg");
  std::size_t ne;
  Point_2 s0, s, t;
  Polygon_2 polygon;

  in >> ne;
  if(ne < 2)
  {
    std::cerr << "Error: expected a polygon in input (at least 2 edges)" << std::endl;
    std::exit(1);
  }

  in >> s0 >> t;
  polygon.push_back(s0);
  polygon.push_back(t);

  while(in >> s >> t)
    polygon.push_back(t);

  std::cout << "Input polygon of size: " << polygon.size() << std::endl;

  CGAL::Real_timer timer;
  timer.start();

  CDT cdt;
  cdt.insert_constraint(polygon.vertices_begin(), polygon.vertices_end(), true);

  timer.stop();
  std::cout << "Built classically: " << timer.time() << std::endl;
  std::cout << "Built classically: " << timer.time() << std::endl;
  std::cout << std::distance(cdt.finite_vertices_begin(), cdt.finite_vertices_end()) << " nv" << std::endl;
  std::cout << std::distance(cdt.finite_edges_begin(), cdt.finite_edges_end()) << " ne" << std::endl;
  std::cout << std::distance(cdt.finite_faces_begin(), cdt.finite_faces_end()) << " nf" << std::endl;
  std::cout << std::distance(cdt.all_faces_begin(), cdt.all_faces_end()) << " anf" << std::endl;

  // Extract ranges
  std::vector<Point_2> points(cdt.number_of_vertices());
  std::unordered_map<CDT::Vertex_handle, std::size_t> vs_to_ids;
  std::size_t id = 0;
  for(auto vit=cdt.finite_vertices_begin(); vit!=cdt.finite_vertices_end(); ++vit)
  {
    points[id] = cdt.point(vit);
    vs_to_ids[vit] = id++;
  }

  std::vector<std::pair<std::size_t, std::size_t> > edges;
  for(auto eit=cdt.finite_edges_begin(); eit!=cdt.finite_edges_end(); ++eit)
  {
    edges.emplace_back(vs_to_ids[eit->first->vertex(Tr::ccw(eit->second))],
                       vs_to_ids[eit->first->vertex(Tr::cw(eit->second))]);
#ifdef CGAL_T2_BUILD_FROM_RANGES_DEBUG
    std::cout << "Edge " << edges.back().first << " " << edges.back().second << std::endl;
#endif
  }

  timer.reset();
  timer.start();

  CDT cdt_b;
  cdt_b.range_input_v2(points, edges);

  timer.stop();
  std::cout << "Built with V2 range_input: " << timer.time() << std::endl;
  std::cout << std::distance(cdt_b.finite_vertices_begin(), cdt_b.finite_vertices_end()) << " nv" << std::endl;
  std::cout << std::distance(cdt_b.finite_edges_begin(), cdt_b.finite_edges_end()) << " ne" << std::endl;
  std::cout << std::distance(cdt_b.finite_faces_begin(), cdt_b.finite_faces_end()) << " nf" << std::endl;
  std::cout << std::distance(cdt_b.all_faces_begin(), cdt_b.all_faces_end()) << " anf" << std::endl;

  std::cout << "Validity: " << cdt_b.is_valid(true) << std::endl;
}

void random_test(const std::size_t n)
{
  CGAL::Random r;
  std::cout << "Random seed = " << r.get_seed() << std::endl;

  std::vector<Point_2> points;
  CGAL::Random_points_in_square_2<Point_2> rpg(1.0, r);
  for(std::size_t i=0; i<n; ++i)
    points.push_back(*rpg++);

  CGAL::Real_timer timer;
  timer.start();

  // -----------------------------------------------------------------------------------------------
  // BIG WARNING: If you insert a range into a DT2, it is Hilbert-sorted so the vertex handles
  // when doing t2.finite_vertices()... will NOT be in the same order as points were in the range
  // passed to the constructor, so you can't use that same range for the reader
  // -----------------------------------------------------------------------------------------------

  Tr T2a;
  CGAL::spatial_sort(points.begin(), points.end(), T2a.geom_traits());

#ifdef CGAL_T2_BUILD_FROM_RANGES_DEBUG
  for(std::size_t i=0; i<n; ++i)
    std::cout << "Point " << i << " = " << points[i] << std::endl;
#endif

  for(const auto& p : points)
    T2a.insert(p);

  timer.stop();

  std::cout << "Built classically: " << timer.time() << std::endl;
  std::cout << std::distance(T2a.finite_vertices_begin(), T2a.finite_vertices_end()) << " nv" << std::endl;
  std::cout << std::distance(T2a.finite_edges_begin(), T2a.finite_edges_end()) << " ne" << std::endl;
  std::cout << std::distance(T2a.finite_faces_begin(), T2a.finite_faces_end()) << " nf" << std::endl;
  std::cout << std::distance(T2a.all_faces_begin(), T2a.all_faces_end()) << " anf" << std::endl;

#ifdef CGAL_T2_BUILD_FROM_RANGES_DEBUG
  CGAL::draw(T2a);
#endif

  // Extract edges
  std::unordered_map<Tr::Vertex_handle, std::size_t> vs_to_ids;
  std::size_t id = 0;
  for(auto vit=T2a.finite_vertices_begin(); vit!=T2a.finite_vertices_end(); ++vit)
    vs_to_ids[vit] = id++;

  std::vector<std::pair<std::size_t, std::size_t> > edges;
  for(auto eit=T2a.finite_edges_begin(); eit!=T2a.finite_edges_end(); ++eit)
  {
    edges.emplace_back(vs_to_ids[eit->first->vertex(Tr::ccw(eit->second))],
                       vs_to_ids[eit->first->vertex(Tr::cw(eit->second))]);
#ifdef CGAL_T2_BUILD_FROM_RANGES_DEBUG
    std::cout << "Edge " << edges.back().first << " " << edges.back().second << std::endl;
#endif
  }

  std::cout << " === " << std::endl;
  std::cout << points.size() << " points and " << edges.size() << " edges" << std::endl;

  std::cout << " === " << std::endl;

  timer.reset();
  timer.start();

  Tr T2bn;
  T2bn.range_input_naive(points, edges);

  timer.stop();
  std::cout << "Built with naive range_input: " << timer.time() << std::endl;
  std::cout << std::distance(T2bn.finite_vertices_begin(), T2bn.finite_vertices_end()) << " nv" << std::endl;
  std::cout << std::distance(T2bn.finite_edges_begin(), T2bn.finite_edges_end()) << " ne" << std::endl;
  std::cout << std::distance(T2bn.finite_faces_begin(), T2bn.finite_faces_end()) << " nf" << std::endl;
  std::cout << std::distance(T2bn.all_faces_begin(), T2bn.all_faces_end()) << " anf" << std::endl;

  std::cout << "Validity: " << T2bn.is_valid(true) << std::endl;

  std::cout << " === " << std::endl;

  timer.reset();
  timer.start();

  Tr T2b0;
  T2b0.range_input_v0(points, edges);

  timer.stop();
  std::cout << "Built with V0 range_input: " << timer.time() << std::endl;
  std::cout << std::distance(T2b0.finite_vertices_begin(), T2b0.finite_vertices_end()) << " nv" << std::endl;
  std::cout << std::distance(T2b0.finite_edges_begin(), T2b0.finite_edges_end()) << " ne" << std::endl;
  std::cout << std::distance(T2b0.finite_faces_begin(), T2b0.finite_faces_end()) << " nf" << std::endl;
  std::cout << std::distance(T2b0.all_faces_begin(), T2b0.all_faces_end()) << " anf" << std::endl;

  std::cout << "Validity: " << T2b0.is_valid(true) << std::endl;

  std::cout << " === " << std::endl;

  timer.reset();
  timer.start();

  Tr T2b1;
  T2b1.range_input_v1(points, edges);

  timer.stop();
  std::cout << "Built with V1 range_input: " << timer.time() << std::endl;
  std::cout << std::distance(T2b1.finite_vertices_begin(), T2b1.finite_vertices_end()) << " nv" << std::endl;
  std::cout << std::distance(T2b1.finite_edges_begin(), T2b1.finite_edges_end()) << " ne" << std::endl;
  std::cout << std::distance(T2b1.finite_faces_begin(), T2b1.finite_faces_end()) << " nf" << std::endl;
  std::cout << std::distance(T2b1.all_faces_begin(), T2b1.all_faces_end()) << " anf" << std::endl;

  std::cout << "Validity: " << T2b1.is_valid(true) << std::endl;

  std::cout << " === " << std::endl;

  timer.reset();
  timer.start();

  Tr T2b2;
  T2b2.range_input_v2(points, edges);
  timer.stop();

  std::cout << "Built with V2 range_input: " << timer.time() << std::endl;
  std::cout << std::distance(T2b2.finite_vertices_begin(), T2b2.finite_vertices_end()) << " nv" << std::endl;
  std::cout << std::distance(T2b2.finite_edges_begin(), T2b2.finite_edges_end()) << " ne" << std::endl;
  std::cout << std::distance(T2b2.finite_faces_begin(), T2b2.finite_faces_end()) << " nf" << std::endl;
  std::cout << std::distance(T2b2.all_faces_begin(), T2b2.all_faces_end()) << " anf" << std::endl;

  std::cout << "Validity: " << T2b2.is_valid(true) << std::endl;

//  std::cout << " === " << std::endl;

//  timer.reset();
//  timer.start();

//  Tr T2b3;
//  T2b3.range_input_v3(points, edges);
//  timer.stop();

//  std::cout << "Built with V3 range_input: " << timer.time() << std::endl;
//  std::cout << std::distance(T2b3.finite_vertices_begin(), T2b3.finite_vertices_end()) << " nv" << std::endl;
//  std::cout << std::distance(T2b3.finite_edges_begin(), T2b3.finite_edges_end()) << " ne" << std::endl;
//  std::cout << std::distance(T2b3.finite_faces_begin(), T2b3.finite_faces_end()) << " nf" << std::endl;
//  std::cout << std::distance(T2b3.all_faces_begin(), T2b3.all_faces_end()) << " anf" << std::endl;

//  std::cout << "Validity: " << T2b3.is_valid(true) << std::endl;
}

int main(int argc, char** argv)
{
  int n = (argc > 1) ? std::atoi(argv[1]) : 100000;

//  constrained_random_test(n);
  random_test(n);

  std::cout << "done" << std::endl;
  return EXIT_SUCCESS;
}

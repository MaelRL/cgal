#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/internal/Generic_P2T2/Periodic_2_Delaunay_triangulation_2_generic.h>

#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K                                                   GT;
typedef GT::Vector_2                                        Vector;

typedef CGAL::Periodic_2_triangulation_vertex_base_2_generic<GT>    Vb;
typedef CGAL::Periodic_2_triangulation_face_base_2_generic<GT>      Fb;

typedef CGAL::Triangulation_data_structure_2<Vb, Fb>        Tds;
typedef CGAL::Periodic_2_Delaunay_triangulation_2_generic<GT, Tds>  PDT;

typedef PDT::Vertex_handle                                  Vertex_handle;
typedef PDT::Point                                          Point;

std::vector<Point> generate_lattice()
{
  std::vector<Point> pts;

  for(int i=0; i<10; ++i)
    for(int j=0; j<10; ++j)
      pts.emplace_back(i * 0.14, j * 0.32);

  return pts;
}

int main()
{
  // Remark: the basis
  // basis = CGAL::make_array(Vector(4, 1), Vector(-2.5, -1));
  // with the point set
  // std::vector<Point> pts { Point(0, 0), Point(-0.4, -0.6) };
  // are not in general position. If using inexact constructions,
  // this yields an inconsistent triangulation (i.e. the are 4
  // points forming a rectangle, but different periodic copies
  // of this rectangle are triangulated differently), making
  // the triangulation break down.
  // --> Use exact constructions for now, in the future hopefully
  // we replace these constructions with predicates.

  CGAL::cpp11::array<Vector, 2> basis;
  // the reduced version of this basis should be equivalent to
  // CGAL::make_array(Vector(-0.5, 1), Vector(1.5, 0));
  basis = CGAL::make_array(Vector(4, 1), Vector(-2.5, -1));

#if 1
  std::vector<Point> pts { Point(0, 0), Point(-0.21, -0.6) };
#else
  std::vector<Point> pts = generate_lattice();
#endif

  PDT T(pts.begin(), pts.end(), basis);
  PDT::Vertex_handle vh = T.insert(Point(12345.6, 5432.1));
  PDT::Vertex_handle vh2 = T.insert(Point(0, 0.1));
  PDT::Vertex_handle vh3 = T.insert(Point(0.5, -0.4));
  PDT::Vertex_handle vh4 = T.insert(Point(0.2, -0.2));
  PDT::Vertex_handle vh5 = T.insert(Point(0, -0.52));
  PDT::Vertex_handle vh6 = T.insert(Point(0.6, 0));
  std::cout << "Is simplicial complex (vh6): " << T.is_simplicial_complex() << std::endl;
  PDT::Vertex_handle vh7 = T.insert(Point(0.8, 0.55));
  std::cout << "Is simplicial complex (vh7): " << T.is_simplicial_complex() << std::endl;

  std::cout << "Number of vertices: " << T.number_of_vertices() << std::endl;
  std::cout << "Number of edges: " << T.number_of_edges() << std::endl;
  std::cout << "Number of faces: " << T.number_of_faces() << std::endl;

  // Draw everything surrounding the vertex in red.
  PDT::Face_handle some_face = PDT::Face_handle();
  std::set<PDT::Face_handle> incident_faces2 = T.incident_faces(vh7);
  for (PDT::Face_handle fh : incident_faces2) {
    fh->has_color = true;
    fh->color = CGAL::Color(255, 0, 0);
    some_face = fh;
  }
  some_face->has_color = true;
  some_face->color = CGAL::Color(255, 0, 192);

  // Make the neighbours of the cell blueish.
  for (int i = 0; i < 3; ++i) {
    PDT::Face_handle nbfh = T.neighbor(some_face, i);

    if (nbfh->has_color) {
      CGAL::Color c = nbfh->color;
      nbfh->color = CGAL::Color(c.red(), c.green(), 128);
    } else {
      nbfh->has_color = true;
      nbfh->color = CGAL::Color(255, 255, 128);
    }
  }

//  draw(T.dt2, "Post-Insertion");
  T.draw_dt2("dt2.off");

  T.convert_to_1_cover();
  T.draw_p2t2("before.off");

  std::cout << "T.p2t2.nv: " << T.p2t2.number_of_vertices() << std::endl;
  std::cout << "T.p2t2.nf: " << T.p2t2.number_of_faces() << std::endl;

  T.insert(Point(0.3, 0.12));

  std::cout << "T.p2t2.nv: " << T.p2t2.number_of_vertices() << std::endl;
  std::cout << "T.p2t2.nf: " << T.p2t2.number_of_faces() << std::endl;

  CGAL_assertion(T.p2t2.tds().is_valid());
  CGAL_assertion(T.p2t2.is_valid());
  std::cout << "All good!" << std::endl;

  T.draw_p2t2("after.off");

  return EXIT_SUCCESS;
}

#ifndef CGAL_ANISOTROPIC_MESH_3_CANVAS_TETRAHEDRON_INTERSECTION_H
#define CGAL_ANISOTROPIC_MESH_3_CANVAS_TETRAHEDRON_INTERSECTION_H

#include <CGAL/Bbox_3.h>
#include <CGAL/box_intersection_d.h>

#include <CGAL/assertions.h>

#include <boost/array.hpp>
#include <boost/functional/hash.hpp>

#include <algorithm>
#include <cstddef>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename Canvas>
struct Intersect_tetrahedra
{
  // In this struct, we are only interested in 'real' intersections, that is
  // if the intersection of the two tets is not simply a common vertex, edge or
  // face of the tets.

  typedef typename Canvas::Kernel                                Kernel;

  typedef typename Kernel::FT                                    FT;
  typedef typename Kernel::Point_3                               Point_3;
  typedef typename Kernel::Segment_3                             Segment_3;
  typedef typename Kernel::Triangle_3                            Triangle_3;
  typedef typename Kernel::Tetrahedron_3                         Tetrahedron_3;

  typedef typename Canvas::Primal_edge                           Primal_edge;
  typedef typename Canvas::Primal_triangle                       Primal_triangle;
  typedef typename Canvas::Primal_tetrahedron                    Primal_tetrahedron;
  typedef typename Primal_edge::BSimplex                         ESimplex;
  typedef typename Primal_triangle::BSimplex                     TrSimplex;
  typedef typename Primal_tetrahedron::BSimplex                  TSimplex;
  typedef typename Canvas::Primal_tetrahedra_container           Primal_tetrahedra_container;
  typedef typename Canvas::Primal_triangles_container            Primal_triangles_container;
  typedef typename Primal_triangles_container::iterator          PTrC_iterator;
  typedef typename Primal_tetrahedra_container::iterator         PTC_iterator;

  typedef CGAL::Bbox_3                                           Bbox;
  typedef CGAL::Box_intersection_d::Box_with_handle_d<FT, 3, PTC_iterator>
                                                                 Box;

  Canvas* canvas;

  enum Point_loc
  {
    OUTSIDE = -1,
    FIRST_VERTEX, SECOND_VERTEX, THIRD_VERTEX, FOURTH_VERTEX,
    INSIDE
  };

  Point_loc is_vertex_contained_in_tet(const Point_3& t,
                                       const Point_3& p, const Point_3& q,
                                       const Point_3& r, const Point_3& s) const
  {
    CGAL::Sign o1 = CGAL::orientation(p, q, r, t);
    CGAL::Sign o2 = CGAL::orientation(p, r, s, t);
    CGAL::Sign o3 = CGAL::orientation(p, s, q, t);
    CGAL::Sign o4 = CGAL::orientation(q, s, r, t);

#if (VERBOSITY > 30)
    std::cout << "is: " << t << " contained in : " << std::endl;
    std::cout << p << " || " << q << " || " << r << " || " << s << std::endl;
    std::cout << "orientations: " << o1 << " " << o2 << " " << o3 << " " << o4 << std::endl;
#endif

    if((o1 >= 0 && o2 >=0 && o3 >= 0 && o4 >= 0) ||
       (o1 <= 0 && o2 <=0 && o3 <= 0 && o4 <= 0))
    {
      if(o1 == 0 && o2 == 0 && o3 == 0)
        return FIRST_VERTEX;
      else if(o1 == 0 && o3 == 0 && o4 == 0)
        return SECOND_VERTEX;
      else if(o1 == 0 && o2 == 0 && o4 == 0)
        return THIRD_VERTEX;
      else if(o2 == 0 && o3 == 0 && o4 == 0)
        return FOURTH_VERTEX;
      else
        return INSIDE;
    }
    return OUTSIDE;
  }

  Point_loc is_vertex_contained_in_tet(const Point_3& t, const TSimplex& tet) const
  {
    Point_3& p = canvas->seeds[tet[0]];
    Point_3& q = canvas->seeds[tet[1]];
    Point_3& r = canvas->seeds[tet[2]];
    Point_3& s = canvas->seeds[tet[3]];
    return is_vertex_contained_in_tet(t, p, q, r, s);
  }

  Point_loc is_vertex_contained_in_tet(const std::size_t point_id,
                                       const TSimplex& tet) const
  {
    Point_3& t = canvas->seeds[point_id];
    return is_vertex_contained_in_tet(t, tet);
  }

  ESimplex get_opposite_edge(const TSimplex& tet,
                             const std::size_t id_1, const std::size_t id_2) const
  {
    ESimplex e;
    std::size_t pos = 0;
    for(std::size_t i=0; i<tet.size(); ++i)
    {
      if(tet[i] != id_1 && tet[i] != id_2)
      {
        CGAL_assertion(pos < e.size() && "id1/id2 was not contained in tet...");
        e[pos++] = tet[i];
      }
    }
    return e;
  }

  bool check_intersection_edge_tet(const ESimplex& e,
                                   const TSimplex& tet) const
  {
    // At this point we haven't detected any degenerate intersection cases that
    // wouldn't really be 'strict' intersections
    Segment_3 edge(canvas->seeds[e[0]], canvas->seeds[e[1]]);
    Tetrahedron_3 tetra(canvas->seeds[tet[0]], canvas->seeds[tet[1]],
                        canvas->seeds[tet[2]], canvas->seeds[tet[3]]);
    return CGAL::do_intersect(edge, tetra);
  }

  TrSimplex get_opposite_face(const TSimplex& tet,
                              const std::size_t id) const
  {
    TrSimplex tr;
    std::size_t pos = 0;
    for(std::size_t i=0; i<tet.size(); ++i)
    {
      if(tet[i] != id)
      {
        CGAL_assertion(pos < tr.size() && "id was not contained in tet...");
        tr[pos++] = tet[i];
      }
    }
    return tr;
  }

  bool check_intersection_triangle_tet(const TrSimplex& tr,
                                       const TSimplex& tet) const
  {
    // At this point we haven't detected any degenerate intersection cases that
    // wouldn't really be 'strict' intersections
    Triangle_3 tri(canvas->seeds[tr[0]], canvas->seeds[tr[1]], canvas->seeds[tr[2]]);
    Tetrahedron_3 tetra(canvas->seeds[tet[0]], canvas->seeds[tet[1]],
                        canvas->seeds[tet[2]], canvas->seeds[tet[3]]);
    return CGAL::do_intersect(tri, tetra);
  }

  void mark_primal_triangles_as_interected(PTC_iterator p_tet_1)
  {
    // very unefficient but it doesn't really matter...
    // a possible fix todo is to include pointers to subfaces in the primal simplex class
    for(int i=0; i<4; ++i)
    {
      boost::array<std::size_t, 3> face;
      const Primal_tetrahedron& tet = *p_tet_1;
      face[0] = tet[i]; face[1] = tet[(i+1)%4]; face[2] = tet[(i+2)%4];
      std::sort(face.begin(), face.end());
      Primal_triangle t(face);
      PTrC_iterator it = canvas->primal_triangles.find(t);
      CGAL_assertion(it != canvas->primal_triangles.end());
      it->m_is_intersected = true;
    }
  }

  void mark_primal_tetrahedra_as_intersected(PTC_iterator p_tet_1,
                                             PTC_iterator p_tet_2)
  {
    p_tet_1->m_is_intersected = true;
    p_tet_2->m_is_intersected = true;

    // only needed for visualization
#ifdef COMPUTE_PRIMAL_ALL_DIMENSIONS
    // if the macro isn't defined, we only care about tets
    mark_primal_triangles_as_interected(p_tet_1);
    mark_primal_triangles_as_interected(p_tet_2);
#endif
  }

  void operator() (const Box* box_1, const Box* box_2)
  {
    PTC_iterator p_tet_1 = box_1->handle();
    PTC_iterator p_tet_2 = box_2->handle();
    const TSimplex& tet_1 = p_tet_1->simplex();
    const TSimplex& tet_2 = p_tet_2->simplex();

    // just to be safe, check that the tets aren't identical...
    CGAL_assertion(tet_1 != tet_2);

    // collect the shared vertices between both tets and get rid of any simple
    // case (one of tet's vertex is included -and not a vertex- of the other one)
    std::vector<std::size_t> shared_vertices_ids;
    for(std::size_t i=0; i<tet_1.size(); ++i)
    {
      Point_loc pl_1 = is_vertex_contained_in_tet(tet_1[i], tet_2);
      if(pl_1 == INSIDE)
      {
        // a vertex of tet_1 is inside tet_2 and is not a vertex of tet_2 => real intersection
        return mark_primal_tetrahedra_as_intersected(p_tet_1, p_tet_2);
      }
      else if(pl_1 == OUTSIDE)
        continue;
      else // pl_1 == FIRST/SECOND/THIRD or FOURTH_VERTEX
      {
        // here, the state pl_1 is the index in tet_2
        shared_vertices_ids.push_back(tet_2[pl_1]);
      }

      Point_loc pl_2 = is_vertex_contained_in_tet(tet_2[i], tet_1);
      if(pl_2 == INSIDE)
      {
        // a vertex of tet_2 is inside tet_1 and is not a vertex of tet_1 => real intersection
        return mark_primal_tetrahedra_as_intersected(p_tet_1, p_tet_2);
      }
      // no need to check the other states for pl_2, a shared vertex is
      // by definition shared, so if tet_2[i] were shared, we would have found
      // (or will find) it through pl_1.
      // (could be done anyway as a tiny optimization todo)
    }

#if (VERBOSITY > 25)
    std::cout << "operator() with tetrahedras: " << std::endl;
    std::cout << tet_1[0] << " " << tet_1[1] << " " << tet_1[2] << " " << tet_1[3] << std::endl;
    std::cout << tet_2[0] << " " << tet_2[1] << " " << tet_2[2] << " " << tet_2[3] << std::endl;
    std::cout << shared_vertices_ids.size() << " shared vertices" << std::endl;
    for(std::size_t i=0; i<shared_vertices_ids.size(); ++i)
      std::cout << shared_vertices_ids[i] << " ";
    std::cout << std::endl;
#endif

    for(std::size_t i=0; i<shared_vertices_ids.size(); ++i)
    {
      CGAL_assertion(std::find(tet_2.begin(), tet_2.end(),
                               shared_vertices_ids[i]) != tet_2.end());
    }

    // now, we know the shared vertices between the two tets, let's detail the cases :
    // - 3 shared vertices :
    //  In that case, if the last 2 vertices are on the same side of the plane
    //  determined by the 3 common vertices, we have an intersection

    // - 2 shared vertices :
    //  Inspect the opposite edges. If either intersects the other tetrahedron,
    //  we have an intersection : the intersection cannot be a vertex of the tetrahedron,
    //  otherwise we would have 3 shared vertices!

    // - 1 shared vertex :
    //  Inspect the opposite faces. If either intersects the other tetrahedron,
    //  we have an intersection (same reasons as above).

    // - 0 shared vertices
    //  Whether an intersection exists is simply given by do_intersect() since
    //  we have cleared off degenerate cases in the previous cases.

    if(shared_vertices_ids.size() == 3)
    {
      // trick to find the ids of the non shared vertices though the difference
      FT shared_ids_sum = shared_vertices_ids[0] + shared_vertices_ids[1] + shared_vertices_ids[2];
      std::size_t id_1 = tet_1[0] + tet_1[1] + tet_1[2] + tet_1[3] - shared_ids_sum;
      std::size_t id_2 = tet_2[0] + tet_2[1] + tet_2[2] + tet_2[3] - shared_ids_sum;

#if (VERBOSITY > 40)
      std::cout << "shared vertices 3" << std::endl;
      std::cout << "non shared ids : id1/2 :" << id_1 << " " << id_2 << std::endl;
#endif

      if(CGAL::orientation(canvas->seeds[shared_vertices_ids[0]],
                           canvas->seeds[shared_vertices_ids[1]],
                           canvas->seeds[shared_vertices_ids[2]],
                           canvas->seeds[id_1]) ==
         CGAL::orientation(canvas->seeds[shared_vertices_ids[0]],
                           canvas->seeds[shared_vertices_ids[1]],
                           canvas->seeds[shared_vertices_ids[2]],
                           canvas->seeds[id_2]))
      {
#if (VERBOSITY > 35)
        std::cout << "shared vertices 3 intersection" << std::endl;
#endif
        mark_primal_tetrahedra_as_intersected(p_tet_1, p_tet_2);
      }
      return;
    }
    else if(shared_vertices_ids.size() == 2)
    {
      const std::size_t shared_vertex_0 = shared_vertices_ids[0];
      const std::size_t shared_vertex_1 = shared_vertices_ids[1];

      const ESimplex opposite_edge_1 = get_opposite_edge(tet_1, shared_vertex_0, shared_vertex_1);
      const ESimplex opposite_edge_2 = get_opposite_edge(tet_2, shared_vertex_0, shared_vertex_1);

#if (VERBOSITY > 40)
      std::cout << "shared vertices 2" << std::endl;
      std::cout << "opposite edge 1 : " << opposite_edge_1[0] << " "
                << opposite_edge_1[1] << std::endl;
      std::cout << "opposite edge 2 : " << opposite_edge_2[0] << " "
                << opposite_edge_2[1] << std::endl;
#endif

      if(check_intersection_edge_tet(opposite_edge_1, tet_2) ||
         check_intersection_edge_tet(opposite_edge_2, tet_1))
      {
#if (VERBOSITY > 35)
        std::cout << "shared vertices 2 intersection" << std::endl;
#endif
        mark_primal_tetrahedra_as_intersected(p_tet_1, p_tet_2);
      }
      return;
    }
    else if(shared_vertices_ids.size() == 1)
    {
      const std::size_t shared_vertex = shared_vertices_ids[0];
      const TrSimplex opposite_face_1 = get_opposite_face(tet_1, shared_vertex);
      const TrSimplex opposite_face_2 = get_opposite_face(tet_2, shared_vertex);

#if (VERBOSITY > 40)
      std::cout << "shared vertices 2" << std::endl;
      std::cout << "opposite face 1: " << opposite_face_1[0] << " "
                << opposite_face_1[1] << " " << opposite_face_1[2] << std::endl;
      std::cout << "opposite face 2: " << opposite_face_2[0] << " "
                << opposite_face_2[1] << " " << opposite_face_2[2] << std::endl;
#endif

      if(check_intersection_triangle_tet(opposite_face_1, tet_2) ||
         check_intersection_triangle_tet(opposite_face_2, tet_1))
      {
#if (VERBOSITY > 35)
        std::cout << "shared vertices 1 intersection" << std::endl;
#endif

        mark_primal_tetrahedra_as_intersected(p_tet_1, p_tet_2);
      }
      return;
    }

    // at this point, we have cleared all the cases that are not really intersections
    // so if we actually have an intersection with CGAL::do_intersect(), we have a 'real'
    // intersection
    Tetrahedron_3 tetra_1(canvas->seeds[tet_1[0]], canvas->seeds[tet_1[1]],
                          canvas->seeds[tet_1[2]], canvas->seeds[tet_1[3]]);
    Tetrahedron_3 tetra_2(canvas->seeds[tet_2[0]], canvas->seeds[tet_2[1]],
                          canvas->seeds[tet_2[2]], canvas->seeds[tet_2[3]]);
    if(CGAL::do_intersect(tetra_1, tetra_2))
    {
#if (VERBOSITY > 35)
        std::cout << "tet-tet intersection" << std::endl;
#endif
      mark_primal_tetrahedra_as_intersected(p_tet_1, p_tet_2);
    }
  }

  void tests()
  {
    CGAL_assertion(is_vertex_contained_in_tet(Point_3(2,2,0), Point_3(2, 2, 0), Point_3(-2, 2, 0), Point_3(0, -2, 0), Point_3(0, 0, 3)) == FIRST_VERTEX);
    CGAL_assertion(is_vertex_contained_in_tet(Point_3(-2,2,0), Point_3(2, 2, 0), Point_3(-2, 2, 0), Point_3(0, -2, 0), Point_3(0, 0, 3)) == SECOND_VERTEX);
    CGAL_assertion(is_vertex_contained_in_tet(Point_3(0,-2,0), Point_3(2, 2, 0), Point_3(-2, 2, 0), Point_3(0, -2, 0), Point_3(0, 0, 3)) == THIRD_VERTEX);
    CGAL_assertion(is_vertex_contained_in_tet(Point_3(0,0,3), Point_3(2, 2, 0), Point_3(-2, 2, 0), Point_3(0, -2, 0), Point_3(0, 0, 3)) == FOURTH_VERTEX);
    CGAL_assertion(is_vertex_contained_in_tet(Point_3(0,0,4), Point_3(2, 2, 0), Point_3(-2, 2, 0), Point_3(0, -2, 0), Point_3(0, 0, 3)) == FOURTH_VERTEX);
    CGAL_assertion(is_vertex_contained_in_tet(Point_3(0,0,0), Point_3(2, 2, 0), Point_3(-2, 2, 0), Point_3(0, -2, 0), Point_3(0, 0, 3)) == INSIDE);
    CGAL_assertion(is_vertex_contained_in_tet(Point_3(0,2,0), Point_3(-2, 2, 0), Point_3(2, 2, 0), Point_3(0, -2, 0), Point_3(0, 0, 3)) == INSIDE);
    CGAL_assertion(is_vertex_contained_in_tet(Point_3(1,0,0), Point_3(-2, 2, 0), Point_3(2, 2, 0), Point_3(0, 0, 3), Point_3(0, -2, 0)) == INSIDE);
    CGAL_assertion(is_vertex_contained_in_tet(Point_3(1,0,-1e-7), Point_3(-2, 2, 0), Point_3(2, 2, 0), Point_3(0, 0, 3), Point_3(0, -2, 0)) == OUTSIDE);

    // todo: real tet-tet intersection tests by refactoring a bit the operator()
  }

  Intersect_tetrahedra(Canvas* canvas_) : canvas(canvas_) { }
};

}  // Anisotropic_mesh_3
}  // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_PRIMAL_SIMPLEX_H

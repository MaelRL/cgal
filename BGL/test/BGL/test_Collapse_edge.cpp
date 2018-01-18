#include "test_Prefix.h"
#include <boost/range/distance.hpp>
#include <CGAL/boost/graph/Euler_operations.h>

template < typename Mesh>
typename boost::graph_traits<Mesh>::
halfedge_descriptor find_halfedge(float x1, float y1,
                                  float x2, float y2,
                                  Mesh& m,
                                  bool is_border = false)
{
  typedef typename boost::property_map<Mesh, CGAL::vertex_point_t>::type VPMAP;
  typedef typename boost::property_traits<VPMAP>::value_type Point;

  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
  VPMAP vpmap = get(CGAL::vertex_point, m);
  BOOST_FOREACH(halfedge_descriptor h, halfedges(m))
  {
    if(get(vpmap, source(h, m)) == Point(x1,y1,0)
       && get(vpmap, target(h, m)) == Point(x2,y2,0))
    {
      if(is_border == CGAL::is_border(h, m))
        return h;
      else
        return opposite(h, m);
    }
  }
  return boost::graph_traits<Mesh>::null_halfedge();
}

template <typename Mesh>
void
collapse_edge_test()
{
  CGAL_GRAPH_TRAITS_MEMBERS(Mesh);
  typedef typename boost::graph_traits<Mesh>:: vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Mesh>:: halfedge_descriptor halfedge_descriptor;

  const std::string fname = "data/flat_hexahedron.off";
  Mesh m;
  if(!CGAL::read_off(fname, m)) {
    std::cout << "Error reading file: " << fname << std::endl;
  }
  assert(CGAL::is_valid(m));

  Mesh test_mesh;
  CGAL::copy_face_graph(m, test_mesh);
  assert(CGAL::is_valid(test_mesh));

  //case 1
  {
    halfedge_descriptor e = find_halfedge(-0.5,0,
                                          0.5,0,
                                          test_mesh);
    halfedge_descriptor en = next(e, test_mesh);
    halfedge_descriptor eno = opposite(en, test_mesh);
    halfedge_descriptor eno_prime = opposite(next(opposite(e, test_mesh), test_mesh), test_mesh);
    vertex_descriptor v1 = target(edge(e, test_mesh), test_mesh);
    assert(CGAL::Euler::collapse_edge(edge(e, test_mesh), test_mesh) == v1);
    char found = 0;
    BOOST_FOREACH(halfedge_descriptor it, CGAL::halfedges_around_target(v1,test_mesh))
    {
      if(it == eno
         || it == eno_prime){
        ++found;
      }
    }
    assert(found == 2);
    CGAL::clear(test_mesh);

  }
  //case 2
  {
    CGAL::copy_face_graph(m, test_mesh);
    assert(CGAL::is_valid(test_mesh));
    halfedge_descriptor e = find_halfedge(0,0.5,
                                          -0.75,0.5,
                                          test_mesh);
    CGAL::Euler::remove_face(e, test_mesh);

    e = find_halfedge(-0.5,0,
                      0.5,0,
                      test_mesh);
    halfedge_descriptor en = next(e, test_mesh);
    halfedge_descriptor eno = opposite(en, test_mesh);
    halfedge_descriptor eno_prime = opposite(next(opposite(e, test_mesh), test_mesh), test_mesh);
    vertex_descriptor v1 = target(edge(e, test_mesh), test_mesh);
    assert(CGAL::Euler::collapse_edge(edge(e, test_mesh), test_mesh) == v1);
    char found = 0;
    BOOST_FOREACH(halfedge_descriptor it, CGAL::halfedges_around_target(v1,test_mesh))
    {
      if(it == eno
         || it == eno_prime){
        ++found;
      }
    }
    assert(found == 2);
    CGAL::clear(test_mesh);
  }
  //case 3
  {
    CGAL::copy_face_graph(m, test_mesh);
    assert(CGAL::is_valid(test_mesh));
    halfedge_descriptor e = find_halfedge(1.5,0,
                                          0.75,0.5,
                                          test_mesh);
    CGAL::Euler::remove_face(e, test_mesh);

    e = find_halfedge(-0.5,0,
                      0.5,0,
                      test_mesh);
    halfedge_descriptor en = next(e, test_mesh);
    halfedge_descriptor eno = opposite(en, test_mesh);
    halfedge_descriptor eno_prime = opposite(next(opposite(e, test_mesh), test_mesh), test_mesh);
    vertex_descriptor v1 = target(edge(e, test_mesh), test_mesh);
    assert(CGAL::Euler::collapse_edge(edge(e, test_mesh), test_mesh) == v1);
    char found = 0;
    BOOST_FOREACH(halfedge_descriptor it, CGAL::halfedges_around_target(v1,test_mesh))
    {
      if(it == eno
         || it == eno_prime){
        ++found;
      }
    }
    assert(found == 2);
    CGAL::clear(test_mesh);
  }
  //case 4
  {
    CGAL::copy_face_graph(m, test_mesh);
    assert(CGAL::is_valid(test_mesh));
    halfedge_descriptor e = find_halfedge(-0.5, 0,
                                          0, -0.5,
                                          test_mesh);
    CGAL::Euler::remove_face(e, test_mesh);
    e = find_halfedge(0, -0.5,
                      -0.5, 0,
                      test_mesh);
    CGAL::Euler::remove_face(e, test_mesh);
    e = find_halfedge(0, -0.5,
                      0.75, -0.5,
                      test_mesh);
    CGAL::Euler::remove_face(e, test_mesh);


    e = find_halfedge(-0.5,0,
                      0.5,0,
                      test_mesh);
    halfedge_descriptor en = next(e, test_mesh);
    halfedge_descriptor eno = opposite(en, test_mesh);
    halfedge_descriptor ep_prime = prev(opposite(e, test_mesh), test_mesh);
    halfedge_descriptor eno_prime = opposite(next(opposite(e, test_mesh), test_mesh), test_mesh);
    vertex_descriptor v1 = target(edge(e, test_mesh), test_mesh);
    assert(CGAL::Euler::collapse_edge(edge(e, test_mesh), test_mesh) == v1);
    char found = 0;
    BOOST_FOREACH(halfedge_descriptor it, CGAL::halfedges_around_target(v1,test_mesh))
    {
      if(it == eno
         || it == eno_prime
         || it == ep_prime){
        ++found;
      }
    }
    assert(found == 3);
    CGAL::clear(test_mesh);
  }
  //case 5
  {
    CGAL::copy_face_graph(m, test_mesh);
    assert(CGAL::is_valid(test_mesh));
    halfedge_descriptor e = find_halfedge(0.75,0.5,
                                          1.5,0,
                                          test_mesh);
    CGAL::Euler::remove_face(e, test_mesh);
    e = find_halfedge(0.75,-0.5,
                      1.5,0,
                      test_mesh);
    CGAL::Euler::remove_face(e, test_mesh);
    e = find_halfedge(0,0.5,
                      0.5,0,
                      test_mesh);
    CGAL::Euler::remove_face(e, test_mesh);
    e = find_halfedge(0.5,0,
                      0,-0.5,
                      test_mesh);
    CGAL::Euler::remove_face(e, test_mesh);

    e = find_halfedge(-0.5,0,
                      0.5,0,
                      test_mesh);
    CGAL::Euler::remove_face(e, test_mesh);
    halfedge_descriptor ep = prev(e, test_mesh);
    halfedge_descriptor eno_prime = opposite(next(opposite(e, test_mesh), test_mesh), test_mesh);
    vertex_descriptor v1 = target(edge(e, test_mesh), test_mesh);
    assert(CGAL::Euler::collapse_edge(edge(e, test_mesh), test_mesh) == v1);
    char found = 0;
    BOOST_FOREACH(halfedge_descriptor it, CGAL::halfedges_around_target(v1,test_mesh))
    {
      if(it == ep)
        ++found;
      else if( it == eno_prime){
        ++found;
      }
    }
    assert(found == 2);
    CGAL::clear(test_mesh);
  }
}


int main()
{

  collapse_edge_test<Polyhedron>();
  collapse_edge_test<SM>();

#ifdef CGAL_USE_OPENMESH
  collapse_edge_test<OMesh>();
#endif

  std::cerr << "done\n";
  return 0;
}

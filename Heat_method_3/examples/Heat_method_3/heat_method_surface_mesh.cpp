#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Heat_method_3/Heat_method_3.h>
#include <CGAL/Timer.h>

#include <iostream>
#include <fstream>
#include <iostream>

#include <boost/foreach.hpp>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point_3;
typedef CGAL::Surface_mesh<Point_3>                          Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;
typedef Surface_mesh::Property_map<vertex_descriptor,double> Vertex_distance_map;
typedef CGAL::Heat_method_3::Heat_method_3<Surface_mesh> Heat_method;



int main(int argc, char* argv[])
{
  CGAL::Timer timer;
  //read in mesh
  Surface_mesh sm;
  const char* filename = (argc > 1) ? argv[1] : "./data/bunny.off";
  std::ifstream in(filename);
  in >> sm;
  //the heat intensity will hold the distance values from the source set
  Vertex_distance_map heat_intensity = sm.add_property_map<vertex_descriptor, double>("v:heat_intensity", 0).first;

  timer.start();
  Heat_method hm(sm);

  //add the first vertex as the source set
  vertex_descriptor source = *(vertices(sm).first);
  hm.add_source(source);
  hm.fill_distance_map(heat_intensity);
  timer.stop();
  std::cout << timer.time() << " sec"<< std::endl;
  timer.reset();

  Point_3 sp = sm.point(source);
  
  vertex_descriptor far;
  double sdistance = 0;
  
  BOOST_FOREACH(vertex_descriptor vd , vertices(sm)){
    // std::cout << vd << "  is at distance " << get(heat_intensity, vd) << " from " << source << std::endl;
    if(squared_distance(sp,sm.point(vd)) > sdistance){
      far = vd;
      sdistance = squared_distance(sp,sm.point(vd));
    }
  }
  timer.start();
  hm.add_source(far);
  hm.fill_distance_map(heat_intensity);
  timer.stop();
  std::cout << timer.time() << " sec"<< std::endl;

  /*
  BOOST_FOREACH(vertex_descriptor vd , vertices(sm)){
    std::cout << vd << "  is at distance " << get(heat_intensity, vd) << " " << std::endl;
  }
  */
  return 0;
}

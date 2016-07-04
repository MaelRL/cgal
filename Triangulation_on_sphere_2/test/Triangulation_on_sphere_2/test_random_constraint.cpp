#include <boost/iterator/transform_iterator.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_sphere_traits_2.h>
#include <CGAL/Projection_sphere_traits_3.h>
#include <CGAL/Delaunay_triangulation_sphere_2.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/point_generators_2.h>
#include <fstream>
#include <CGAL/Timer.h>
#include <CGAL/squared_distance_3.h>
#include <cmath>


//convex Hull
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/convex_hull_3_to_polyhedron_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/convex_hull_incremental_3.h>
#include <CGAL/Timer.h>
#include <vector>


#include <CGAL/Constrained_Delaunay_triangulation_sphere_2.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Delaunay_triangulation_sphere_traits_2<K>             Gt;
typedef CGAL::Projection_sphere_traits_3<K>							Gt2;
typedef CGAL::Constrained_Delaunay_triangulation_sphere_2<Gt>       CTOS;
typedef K::Point_3                                                 Point;
typedef  CTOS::Vertex_handle								Vertex_handle;




int main(){
	
	
	 int nu_of_pts;
	double radius;
	nu_of_pts =pow(2,15);
	radius=sqrt(110);
	CGAL::Timer time;

	CGAL::Random random(nu_of_pts);
	typedef CGAL::Creator_uniform_3<double, Point> Creator;
    CGAL::Random_points_on_sphere_3<Point, Creator> on_sphere(radius);
	
	
	std::vector<Point> points;
	std::vector<Point> points2;
		
	points2.push_back(Point(0,0,0));
	
	for (int count=0; count<nu_of_pts; count++) {
		Point p = *on_sphere;
		points.push_back(p);
		points2.push_back(p);
		on_sphere++;
	}
	//Delaunay_traits
	CTOS rtos;
	rtos.set_radius(radius);

	CTOS cdts;
	cdts.set_radius(10);
	
	
	std::cout<<" ***STARTING***"<<std::endl;
	time.start();
	Point p0 =Point(7, -5, 6);
    Point p3 =Point(-6, 5, -7);
	Point p2 =Point(5, 7, 6);
	Point p1 =Point(5, 6, 7);
	
	 rtos.insert_constraint(p0,p1);
	rtos.insert_constraint(p1,p2);
	rtos.insert_constraint(p2,p3);
	rtos.insert_constraint(p3,p0);

	rtos.insert(points.begin(),points.end());
	time.stop();
	rtos.is_valid();
	std::cout<<"triangulation sphere    "<< time.time()<<std::endl;
	std::cout<<"number of vertices   "<<rtos.number_of_vertices()<<std::endl;
	
	
	assert((int)rtos.number_of_vertices() == nu_of_pts+4);
	
	
	
}
	

	
	
	
	
	
	
	

//Claudia WERNER
#ifndef CGAL_DELAUNAY_MESH_SPHERE_TRAITS_2_H
#define CGAL_DELAUNAY_MESH_SPHERE_TRAITS_2_H

#include <CGAL/number_utils_classes.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/basic_constructions_2.h>
#include <CGAL/distance_predicates_2.h>

#include <CGAL/triangulation_assertions.h>
//Claudia WERNER


namespace CGAL { 
//projects the computed points (e.g. midpoint, circumcenter) on the sphere			
template <class K, class Predicate_>
class Project_on_sphere_adaptor{
public:
	
	typedef Predicate_ Predicate;
	typedef typename K::Point_3 Point;
	Project_on_sphere_adaptor(double radius,Point sphere):_radius(radius), _sphere(sphere){}
double _radius;
	Point _sphere;
	typedef typename Predicate::result_type result_type;
	
	
	result_type operator()(const Point& p0, const Point& p1)  {
		Point p = Predicate()(p0,p1);
		return project(p);
	}
		
	result_type operator()(const Point& p0, const Point& p1, const Point& p2)  {
		Point p = Predicate()(p0,p1,p2);
		return project(p);
	}
	
private:
	Point project(const Point& p)
		{return project(p.x(), p.y(), p.z());}
					   
	Point project(double x, double y, double z){
		double dist = x*x+y*y+z*z;
		if (dist ==0) return Point(0,0,0);
		double scale= _radius/sqrt(dist);
		return Point(x*scale, y*scale, z*scale);
	}
	

};
	
template < class Base >
class Delaunay_mesh_sphere_traits_2
  : public Base
{
public:
	typedef Delaunay_mesh_sphere_traits_2<Base> Self;
	typedef typename Base::Point_2 Point_2;
	typedef typename Base::Kernel Kernel;
	typedef typename Kernel::Vector_3 Vector_2;
	typedef typename Kernel::Compute_area_3 Compute_area_2;
	typedef typename Kernel::Angle_3 Angle_2;
	typedef typename Kernel::Compare_distance_3 Compare_distance_2;
	typedef typename Kernel::Construct_vector_3 Construct_vector_2;
	typedef typename Kernel::Construct_scaled_vector_3 Construct_scaled_vector_2;
	typedef typename Kernel::Construct_translated_point_3 Construct_translated_point_2;
	typedef typename Kernel::Construct_circumcenter_3 Circumcenter_2;
	
	typedef Project_on_sphere_adaptor<Kernel,typename Kernel::Construct_midpoint_3> Construct_midpoint_2;
	typedef Project_on_sphere_adaptor<Kernel,typename Kernel::Construct_circumcenter_3> Construct_circumcenter_2;
	typedef boost::true_type requires_test;
	double _radius;
	double _minDistSquared;
	void set_radius(double a)
	{_radius = a;
		_minDistSquared = (_radius*pow (2, -23))*(_radius*pow (2, -23));
	}
	double radius(){return _radius;}
	Delaunay_mesh_sphere_traits_2(double radius=1);
	
	Compare_distance_2
	compare_distance_2_object()const
	{return Compare_distance_2();}
	
	Compute_area_2
	compute_area_2_object()const
	{return Compute_area_2();}
		
	Angle_2
	angle_2_object()const
	{return Angle_2();}
	
	Construct_vector_2
	construct_vector_2_object() const
	{return Construct_vector_2();}
	
	Construct_scaled_vector_2
	construct_scaled_vector_2_object()const
	{return Construct_scaled_vector_2();}
	
	Construct_translated_point_2
	construct_translated_point_2_object()const
	{return Construct_translated_point_2();}
	
		
	Construct_midpoint_2
	construct_midpoint_2_object()const
	{return Construct_midpoint_2(_radius, Base::_sphere);}
	
	Construct_circumcenter_2
	construct_circumcenter_2_object()const
	{return Construct_circumcenter_2(_radius, Base::_sphere);}
			
protected:
	Point_2 _sphere;
};	

template <class R>
	Delaunay_mesh_sphere_traits_2<R>::
	Delaunay_mesh_sphere_traits_2(double radius)
	:_radius(radius){}
	
	
	
} //namespace CGAL

#endif // CGAL_Reg_TRIANGULATION_SPHERE_TRAITS_2_H














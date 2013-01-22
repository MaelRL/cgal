
#ifndef CGAL_DELAUNAY_MESH_SPHERE_TRAITS_2_H
#define CGAL_DELAUNAY_MESH_SPHERE_TRAITS_2_H

#include <CGAL/number_utils_classes.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Line_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/basic_constructions_2.h>
#include <CGAL/distance_predicates_2.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Segment_2_Segment_2_intersection.h>


namespace CGAL { 
	/*
template <class Base>
class Compute_area_2{
public:
	typedef typename Base::Point_2   Point_2;
	typedef typename Base::Kernel Kernel;
	typedef double result_type;
	result_type 
	operator()( const Point_2& p, const Point_2& q, const Point_2& r ) const
	{return compute_area_2_object()(p,q);}
	
	
//protected:
	//Point_2 _sphere;
};
*/
/*
template <class Base>
class Construct_midpoint_2{
public:
	typedef typename Base::Point_2   Point_2;
	
	typedef Point_2 result_type;
	result_type
	operator()(const Point_2& p, const Point_2&q)
	{
		
		return midpoint(p,q);}
	
	};*/
	
	
	
template < class Base >
class Delaunay_mesh_sphere_traits_2
  : public Base
{
public:
	typedef Delaunay_mesh_sphere_traits_2<Base> Self;
	//typedef typename CGAL::Compute_area_2<Self> Compute_area_2;
	//typedef typename CGAL::Construct_midpoint_2<Self> Construct_midpoint_2;
	typedef typename Base::Point_2 Point_2;
		typedef typename Base::Kernel Kernel;
	typedef typename Kernel::Vector_3 Vector_2;
	typedef typename Kernel::Compute_area_3 Compute_area_2;
	typedef typename Kernel::Construct_midpoint_3 Construct_midpoint_2;
	typedef typename Kernel::Angle_3 Angle_2;
	typedef typename Kernel::Construct_vector_3 Construct_vector_2;
	typedef typename Kernel::Construct_scaled_vector_3 Construct_scaled_vector_2;
	typedef typename Kernel::Construct_translated_point_3 Construct_translated_point_2;
	Compute_area_2
	compute_area_2_object()const
	{return Compute_area_2();}
	
	Construct_midpoint_2
	construct_midpoint_2_object()const
	{return Construct_midpoint_2();}
	
	
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
	
	/*Compute_area_2
	compute_area_2_object()const
	{return Compute_area_2();}
	
	Construct_midpoint_2
	construct_midpoint_2_object()const
	{return Construct_midpoint_2();}*/
protected:
	Point_2 _sphere;
};	
	
/*
template < class R >
Delaunay_mesh_sphere_traits_2<R> ::
Delaunay_mesh_sphere_traits_2(const Point_2& sphere)
: _sphere(sphere)
	{}
*/
	
} //namespace CGAL

#endif // CGAL_Reg_TRIANGULATION_SPHERE_TRAITS_2_H














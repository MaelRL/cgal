
#ifndef CGAL_REAL_REGULAR_TRIANGULATION_SPHERE_TRAITS_2_H
#define CGAL_REAL_REGULAR_TRIANGULATION_SPHERE_TRAITS_2_H

//#include <CGAL/Weighted_point.h>
#include <CGAL/number_utils_classes.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Delaunay_triangulation_sphere_traits_2.h>


namespace CGAL { 

template<class Kernel>
class Power_test_on_sphere{
  typedef typename Kernel::Weighted_point_3    Weighted_point_3;
  typedef typename Kernel::FT                  FT;
  
  FT radius_;
  
public:
  typedef Sign result_type;

  Power_test_on_sphere(const FT& radius=1):radius_(radius){}

  result_type
  operator()(const Weighted_point_3& p,
             const Weighted_point_3& q,
             const Weighted_point_3& r,
             const Weighted_point_3& s ) const 
  {
    FT cp=2*radius_-p.weight();
    FT cq=2*radius_-q.weight();
    FT cr=2*radius_-r.weight();
    FT cs=2*radius_-s.weight();
    
    return
      sign (determinant(
              q.x()*cp -p.x()*cq , r.x()*cp -p.x()*cr , s.x()*cp -p.x()*cs,
              q.y()*cp -p.y()*cq , r.y()*cp -p.y()*cr , s.y()*cp -p.y()*cs,
              q.z()*cp -p.z()*cq , r.z()*cp -p.z()*cr , s.z()*cp -p.z()*cs
            )
      );
    
  }

  result_type
  operator()(const Weighted_point_3& p,
             const Weighted_point_3& q,
             const Weighted_point_3& r
            ) const 
  {
    return ON_ORIENTED_BOUNDARY;
  }
  
  result_type
  operator()(const Weighted_point_3& p,
             const Weighted_point_3& q
            ) const 
  {
    return ON_ORIENTED_BOUNDARY;
  }  
};


template < class R, class Wpoint >
class Real_regular_triangulation_sphere_traits_2
  : public R
{
public:
  typedef Delaunay_triangulation_sphere_traits_2<R>                     Base;
  typedef Wpoint                                               Point_2; 
  typedef Wpoint                                               Weighted_point_2;
                      

  typedef Real_regular_triangulation_sphere_traits_2<R,Wpoint>   Self;
  typedef CGAL::Power_test_on_sphere<Self>    Power_test_2;
  typedef CGAL::Orientation_sphere_2<Self>    Orientation_2;
  typedef CGAL::Coradial_sphere_2<Self>       Coradial_sphere_2;
  typedef CGAL::Inside_cone_2<Self>           Inside_cone_2;
  typedef CGAL::Orientation_sphere_1<Self>    Orientation_1;

  double _radius;

  Real_regular_triangulation_sphere_traits_2(const Point_2& sphere=typename Point_2::Point(0,0,0));
  
	void set_radius(double radius){}
	typedef boost::true_type requires_test;
	
  Orientation_2
  orientation_2_object()const
  {return Orientation_2(_sphere);}

  Orientation_1
  orientation_1_object() const {
    return Orientation_1(_sphere);
  }

  Power_test_2 
  power_test_2_object() const
    {  return Power_test_2();}

  Coradial_sphere_2
  coradial_sphere_2_object() const
  {return Coradial_sphere_2(_sphere);}

  Inside_cone_2
  inside_cone_2_object() const {
    return Inside_cone_2(_sphere);
  }

protected :
  Point_2 _sphere;

};

template < class R,class W >
Real_regular_triangulation_sphere_traits_2<R,W> ::
Real_regular_triangulation_sphere_traits_2(const Point_2& sphere)
: _sphere(sphere), _radius(1)
{}
 


} //namespace CGAL

#endif // CGAL_REAL_REGULAR_TRIANGULATION_SPHERE_TRAITS_2_H

#ifndef CONSTRAIN_SURFACE_3_H
#define CONSTRAIN_SURFACE_3_H

#include <CGAL/enum.h>
#include <CGAL/Object.h>

template<typename Kernal, typename Point_container>
class Constrain_surface_3{
	typedef typename Kernal::Point_3	  Point_3;
	typedef typename Kernal::Object_3	  Object_3;
	typedef typename Kernal::Segment_3	  Segment_3;
	typedef typename Kernal::Ray_3		  Ray_3;
	typedef typename Kernal::FT			  FT;

	public:
	
	typedef typename CGAL::Oriented_side	  Oriented_side;  
	
	///\param obj can be a Segment_3 or a Ray_3
	///a function that returns the intersection point of the object
	/// and the surface. 
	/// Precondition: the intersection exist.
	/// TODO Change so that it allows return null intersection.
	virtual Object_3 intersection(const Object_3 &obj)const=0;

	virtual FT get_bounding_radius() const=0;

	///return if the point is inside the constraint (ON_POSITIVE_SIDE)
	/// on the bordary of constraint (ON_ORIENTED_BOUNDARY)
	/// or out of the bordary of constraint (ON_NEGATIVE_SIDE)
	virtual Oriented_side side_of_constraint(const Point_3 &p)const=0;

	virtual Point_container initial_points()const=0;

	virtual ~Constrain_surface_3(){};
};

#endif //CONSTRAIN_SURFACE_3_H

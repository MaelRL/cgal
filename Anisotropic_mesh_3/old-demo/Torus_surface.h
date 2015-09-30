#ifndef TORUS_SURFACE_H
#define TORUS_SURFACE_H

#include <CGAL/enum.h>
#include <cmath>
#include <vector>

#ifdef DEBUG_SURFACE
#include <iostream>
#endif //DEBUG_SURFACE

#include <CGAL/Constrain_surface_3.h>	//extended incl
#include <CGAL/IO/Debug_print.h>		//extented incl

#define PI (2.*acos(0.))

/** A 3-D torus within X-Z plain, with assigned
 *  main_radius and minor_radius
 */
template<typename K>
class Torus_surface : Constrain_surface_3<K,std::vector<typename K::Point_3> >{

	public:

		typedef typename K::Point_3				  Point_3;
		typedef typename K::Object_3			  Object_3;
		typedef typename K::Segment_3			  Segment_3;
		typedef typename K::Ray_3				  Ray_3;
		typedef typename K::Vector_3			  Vector_3;
		typedef typename K::FT					  FT;
		typedef typename std::vector<Point_3>		  Point_container;
		typedef typename CGAL::Oriented_side	  Oriented_side;
		/**construction:
		 * main_r and minor_r must be positive
		 * and main_r > minor_r
		 * */
		Torus_surface(FT main_r,FT minor_r):
			main_radius(main_r),minor_radius(minor_r){
			}

		/**construction*/
		Torus_surface(const Torus_surface &cs){
			main_radius=cs.get_main_radius();
			minor_radius=cs.get_minor_radius();
		}


		Object_3 intersection(const Object_3 &obj)const{
			Segment_3 seg;
			Ray_3 ray;

			if(CGAL::assign(seg,obj)){
				return intersection_of_segment(seg);
			}else if(CGAL::assign(ray,obj)){
				return intersection_of_ray(ray);
			}else{
				return Object_3();
			}
		}

		FT	get_bounding_radius()const {
			return main_radius + minor_radius;
		} 

	private:
		/// \return The intersection of the ray with the surface.
		///  return Object() if intersection does not exist.
		Object_3 intersection_of_ray(const Ray_3 &ray) const {

			Point_3 source_point = ray.source();
			Vector_3 direction = ray.to_vector();
			FT vec_length = sqrt(direction.squared_length());
			FT factor = 2 * (main_radius + minor_radius)/vec_length;
			Point_3 outer_point = source_point + (factor * direction);

#ifdef DEBUG_SURFACE
			std::cout<<"ray intersection ";
#endif
			return intersection_of_segment(Segment_3(source_point,outer_point));
		}


		///\return Intersection of the segment p,q and the surface
		/// Precondition: the oriented_side of p and q are different
		Object_3 intersection_of_segment(const Segment_3 &seg) const {
			Point_3 ip = seg.source();
			Point_3 op = seg.target();
#ifdef DEBUG_SURFACE
			std::cout<<"main_radius = "<<main_radius<<"  minor_radius = "<<minor_radius;
			std::cout<<" between points:";
			print_point(ip);
			std::cout<<"and";
			print_point(op);
			std::cout.flush();
#endif
			Oriented_side osi = side_of_constraint(ip);
			Oriented_side oso = side_of_constraint(op);
			if(oso==CGAL::ON_ORIENTED_BOUNDARY){
			  return make_object(op);
			}
			if(osi==CGAL::ON_ORIENTED_BOUNDARY){
			  return make_object(ip);
			}
			if(osi == oso){
				return Object_3();
			}
			if(osi == CGAL::ON_NEGATIVE_SIDE){
				Point_3 temp = ip;
				ip = op;
				op = temp;
			}
			while(distance_bound(op,ip) > 0){
#ifdef DEBUG_SURFACE
				std::cout<<"between \n";
				print_point(ip);
				print_point(op);
#endif
			  Point_3 m=mid_point(op,ip);
			  Oriented_side os = side_of_constraint(m);
			  if(os==CGAL::ON_ORIENTED_BOUNDARY){
				return make_object(m);
			  }else if(os==CGAL::ON_POSITIVE_SIDE){
				if(ip == m){
					return make_object(ip);
				}
				ip = m;
			  }else{
				if(op == m){
					return make_object(op);
				}
				op = m;
			  }
			}//while
			return make_object((ip));
		}

		Point_3 mid_point(const Point_3 &p,const Point_3 &q)const{
		  return Point_3((p.x()+q.x())/2,(p.y()+q.y())/2,(p.z()+q.z())/2);
		}

	public:
		/** points can be expressed with 2 variable:
		 * the main angle \f$ \phi \f$  and the minor angle \f$ \theta \f$ 
		 * 
		 * projection on the surface:
		 *
		 * \f{array}{(c)c(c)\}}{
		 * \x &   & cos(\phi)R+cos(\theta)cos(\phi)r	\\
		 * \y &	= & sin(\theta)r					    \\
		 * \z &   & sin(\phi)R+cos(\theta)sin(\phi)r
		 * \f}
		 *
		 */
		Point_3 projection(const Point_3 &p) const {	
			FT phi,theta;
			surface_coordinate(p,phi,theta);
			Point_3 projected;
			cardisian_coordinate(phi,theta,projected);
			return projected;
		}

		///return if the point is inside the constraint (ON_POSITIVE_SIDE)
		/// on the bordary of constraint (ON_ORIENTED_BOUNDARY)
		/// or out of the bordary of constraint (ON_NEGATIVE_SIDE)
		Oriented_side side_of_constraint(const Point_3 &p)const{
			FT phi,theta;
			surface_coordinate(p,phi,theta);
			Point_3 center_torus(cos(phi)*main_radius,0,sin(phi)*main_radius);

			FT diff_x=p.x()-center_torus.x();
			FT diff_y=p.y()-center_torus.y();
			FT diff_z=p.z()-center_torus.z();
			FT sq_dist=diff_x*diff_x + diff_y*diff_y + diff_z*diff_z;
			if(sq_dist < minor_radius*minor_radius){
				return CGAL::ON_POSITIVE_SIDE;
			}else if (sq_dist > minor_radius*minor_radius){
				return CGAL::ON_NEGATIVE_SIDE;
			}else{
				return CGAL::ON_ORIENTED_BOUNDARY;
			}
		}

		/**Return the Point_container of initial points
		 * enough to construct a homogenius surface
		 */
		Point_container initial_points()const{
			Point_container	l_points;
			Point_3 p;
			int i=0;
			for(FT phi=-PI; phi<PI; phi+=PI/6+1e-2){
				FT theta=-PI;
				if(i%2==0){
					theta += PI/8;
				}
				for(; theta < PI; theta += PI/4){
					cardisian_coordinate(phi,theta,p);
					l_points.push_back(p);
				}
				i++;
			}
			return l_points;
		}

		FT get_main_radius() const {
			return main_radius;
		}

		FT get_minor_radius() const {
			return minor_radius;
		}

	private:

		/** Given a point, calculate its main angle phi and minor angle theta
		 * \f[ \phi\in (-\pi,\pi], \theta\in (-\pi,\pi] \f]
		 */
		void surface_coordinate(Point_3 p,FT &phi,FT &theta)const{
			phi = atan2(p.z(),p.x());
			Point_3 center_torus(cos(phi)*main_radius,0,sin(phi)*main_radius);
			FT diff_x=p.x()-center_torus.x();
			FT diff_z=p.z()-center_torus.z();
			FT dist = sqrt(diff_x * diff_x + diff_z * diff_z);
			theta = atan(p.y()/dist);
			if(fabs(p.x())<fabs(center_torus.x())){
				theta = PI - theta;
				if(theta>PI){
					theta -= 2* PI;
				}
			}
		}

		/**Given a set of angles \f$ \theta \phi\f$ calculate the 
		 * corresponding point on the surface;
		 **/
		void cardisian_coordinate(FT phi,FT theta,Point_3 &p)const{
			Point_3 center_torus(cos(phi)*main_radius,0,sin(phi)*main_radius);
			Point_3 vec(cos(theta)*minor_radius*cos(phi),
					sin(theta)*minor_radius,
					cos(theta)*minor_radius*sin(phi));
			p=Point_3(center_torus.x()+vec.x(),
					center_torus.y()+vec.y(),
					center_torus.z()+vec.z());
		}
	
		/** The maximun difference between the coordinates of two points
		*	which can be used as the approximation of the distance between two points
		*   \return \f$ max(|x_p-x_q|,|y_p-y_q|,|z_p-z_q|)\f$
		*/
		FT	distance_bound(const Point_3 &p,const Point_3 &q)const {
			FT max= fabs(p.x()-q.x());
			FT dis= fabs(p.y()-q.y());
			if(dis>max){
			  max = dis;
			}
			dis = fabs(p.z()-q.z());
			if(dis>max){
			  return dis;
			}
			return max;
		}


		FT main_radius;
		FT minor_radius;
};//class torus

#undef PI

#endif //TORUS_SURFACE_H

// Copyright (c) 2003-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Laurent RINEAU, Claudia WERNER

#ifndef CGAL_DELAUNAY_MESH_SPHERE_SIZE_CRITERIA_2_H
#define CGAL_DELAUNAY_MESH_SPHERE_SIZE_CRITERIA_2_H

#include <CGAL/Mesh_2/Face_badness.h>
#include <CGAL/Delaunay_mesh_criteria_2.h>
#include <utility>
#include <ostream>

namespace CGAL {
	
	template <class CDT>
	class Delaunay_mesh_sphere_size_criteria_2 : 
    public virtual Delaunay_mesh_criteria_2<CDT>
	{
	protected:
		typedef typename CDT::Geom_traits Geom_traits;
		double sizebound;
		Geom_traits traits;
		
	public:
		typedef Delaunay_mesh_criteria_2<CDT> Base;
		
		Delaunay_mesh_sphere_size_criteria_2(const double aspect_bound = 0.125, 
											 const double size_bound = 0,
											 const Geom_traits& traits = Geom_traits())
		: Base(aspect_bound, traits), sizebound(size_bound),traits(traits) {}
		
		inline
		double size_bound() const { return sizebound; }
		
		inline
		void set_size_bound(const double sb) { sizebound = sb; }
		
		// first: squared_minimum_sine
		// second: size
		
		struct Quality : public std::pair<double, double>
		{
			typedef std::pair<double, double> Base;
			bool irreparable;
			
			Quality() : Base() {irreparable=false;};
			Quality(double _sine, double _size) : Base(_sine, _size) {}
			
			const double& size() const { return second; }
			const double& sine() const { return first; }
			
			// q1<q2 means q1 is prioritised over q2
			// ( q1 == *this, q2 == q )
			bool operator<(const Quality& q) const
			{
				if( size() > 1 )
					if( q.size() > 1 )
						return ( size() > q.size() );
					else
						return true; // *this is big but not q
					else
						if( q.size() >  1 )
							return false; // q is big but not *this
				return( sine() < q.sine() );
			}
			
			std::ostream& operator<<(std::ostream& out) const
			{
				return out << "(size=" << size()
				<< ", sine=" << sine() << ")";
			}
		};
		
		class Is_bad: public Base::Is_bad
		{
		protected:
			const double squared_size_bound; // squared size bound on edge length
			const Geom_traits& traits;
		public:
			typedef typename Base::Is_bad::Point_2 Point_2;
			typedef typename CDT::Vertex_handle Vertex_handle;
			typedef typename CDT::Face_handle Face_handle;
			
			Is_bad(const double aspect_bound,
				   const double size_bound,
				   const Geom_traits& traits)
			: Base::Is_bad(aspect_bound, traits),
			squared_size_bound(size_bound * size_bound),traits(traits) {init(traits._radius);}
			double _radius;
			double _minDist;
			void init(double r){
				_radius = r;
				_minDist = _radius*pow(2,-23)*_radius*pow(2,-23);
			}
			
			
						
			Vertex_handle
			nearest_vertex(const Point_2& p, Face_handle f) const
			{
				typename Geom_traits::Compare_distance_2 
				compare_distance =  traits.compare_distance_2_object();
				Vertex_handle nn =  f->vertex(0);
				if ( compare_distance(p,  f->vertex(1)->point(),nn->point()) == SMALLER) 
					nn=f->vertex(1);
				if ( compare_distance(p,f->vertex(2)->point(), nn->point()) == SMALLER) 
					nn=f->vertex(2);
				
				look_nearest_neighbor(p,f,0,nn);
				look_nearest_neighbor(p,f,1,nn);
				look_nearest_neighbor(p,f,2,nn);
				return nn;
			}
			
			
			
			
			void
			look_nearest_neighbor(const Point_2& p, Face_handle f, int i,
								  Vertex_handle& nn) const
			{
				Face_handle  ni=f->neighbor(i);
				if ( ON_NEGATIVE_SIDE == traits.orientation_2_object()(ni->vertex(0)->point(),ni->vertex(1)->point(),ni->vertex(2)->point(),p) ) return;
				
				typename Geom_traits::Compare_distance_2 
				compare_distance =  traits.compare_distance_2_object();
				i = ni->index(f);
				if ( compare_distance(p,  ni->vertex(i)->point(), nn->point())  == SMALLER)  nn=ni->vertex(i);
				
				// recursive exploration of triangles whose circumcircle contains p
				look_nearest_neighbor(p, ni, CDT::ccw(i), nn);
				look_nearest_neighbor(p, ni, CDT::cw(i), nn);
			} 
			
						
			
			
			Mesh_2::Face_badness operator()(const Quality q) const
			{
				if(q.irreparable){
					return Mesh_2::NOT_BAD;
				}
				if( q.size() > 1 )
					return Mesh_2::IMPERATIVELY_BAD;
				if( q.sine() < this->B )
					return Mesh_2::BAD;
				else
					return Mesh_2::NOT_BAD;
			}
			
			Mesh_2::Face_badness operator()(const typename CDT::Face_handle& fh,
											Quality& q) const
			{
				
				typedef typename CDT::Geom_traits Geom_traits;
				typedef typename Geom_traits::Compute_area_2 Compute_area_2;
				typedef typename Geom_traits::Compute_squared_distance_2
								 Compute_squared_distance_2;
				typedef typename Geom_traits::Construct_circumcenter_2
				                 Construct_circumcenter_2;
				typedef typename Geom_traits::Construct_midpoint_2
								 Construct_midpoint_2;
				
				Compute_squared_distance_2 squared_distance = 
				traits.compute_squared_distance_2_object();
				
				const Point_2& pa = fh->vertex(0)->point();
				const Point_2& pb = fh->vertex(1)->point();
				const Point_2& pc = fh->vertex(2)->point();
				
				double
				a = CGAL::to_double(squared_distance(pb, pc)),
				b = CGAL::to_double(squared_distance(pc, pa)),
				c = CGAL::to_double(squared_distance(pa, pb));
				
				
				double max_sq_length; // squared max edge length
				double second_max_sq_length;
				
				
								
				if(a<b)
				{
					if(b<c) {
						max_sq_length = c;
						second_max_sq_length = b;
					}
					else { // c<=b
						max_sq_length = b;
						second_max_sq_length = ( a < c ? c : a );
					}
				}
				else // b<=a
				{
					if(a<c) {
						max_sq_length = c;
						second_max_sq_length = a;
					}
					else { // c<=a
						max_sq_length = a;
						second_max_sq_length = ( b < c ? c : b );
					}
				}
					
				Construct_circumcenter_2 circumcenter_2 = traits.construct_circumcenter_2_object();
				Point_2 cc = circumcenter_2(pa,pb,pc);
				Vertex_handle nearest = nearest_vertex(cc, fh);
				
				//compute the distance of circumcenter of to the nearest neighbor. Since the triangulation
				//is constrained this can be a vertex not belonging to fh
				//if distance to smal ->NOT_BAD (no point can be inserted)
				double circumradius= CGAL::to_double(squared_distance(cc, nearest->point()));
				if (circumradius<=_minDist){
					q.irreparable = true;
					return Mesh_2::NOT_BAD;
				}
				
				
				
				//circumcenter  in the triangle?
				Orientation o1 = traits.orientation_2_object()(pa, pb,cc);
				Orientation o2 = traits.orientation_2_object()(pb, pc,cc);
				Orientation o3 = traits.orientation_2_object()(pc,pa, cc);
				
				//if circumcenter lies outside the triangle:
				if (o1!=POSITIVE || o2 != POSITIVE || o3!=POSITIVE){
					//longest edge constrained -> NOT_BAD, else circumcenter would lie outside 
					//the domain. Notice: this is possible since on a sphere it is not always possible
					// to make a constrained edge Gabriel conform
					if ((max_sq_length ==a && fh->is_constrained(0))||
						(max_sq_length ==b && fh->is_constrained(1))||
						(max_sq_length ==c && fh->is_constrained(2))){
						q.irreparable = true;
						return Mesh_2:: NOT_BAD;
					}
					
					//has the neighbor of fh (adjacent to the longest edge) a constrained edge? If yes ->NOT_BAD
					// else there could be an encroached edge which can not be split ->infinity loop
					if((max_sq_length == a&& fh->neighbor(0)->is_constrained(0) || fh->neighbor(0)->is_constrained(1) || fh->neighbor(0)->is_constrained(2))||
					   (max_sq_length == b&& fh->neighbor(1)->is_constrained(0) || fh->neighbor(1)->is_constrained(1) || fh->neighbor(1)->is_constrained(2))||
					   (max_sq_length == b&& fh->neighbor(2)->is_constrained(0) || fh->neighbor(2)->is_constrained(1) || fh->neighbor(2)->is_constrained(2))){
						q.irreparable = true;
						return Mesh_2:: NOT_BAD;
					   }

					}
				
				//We can be sure that there are not more problems with to closed vertices
				//go on with the same computation as in Mesh_2
				
				
				
				q.second = 0;
				if( squared_size_bound != 0 )
				{
					//	  std::cerr << squared_size_bound << std::endl;
					q.second = max_sq_length / squared_size_bound;
					// normalized by size bound to deal
					// with size field
					if( q.size() > 1 )
					{
						q.first = 1; // (do not compute sine)
						return Mesh_2::IMPERATIVELY_BAD;
					}
				}
				
				Compute_area_2 area_2 = traits.compute_area_2_object();
				
				double area = 2*CGAL::to_double(area_2(pa, pb, pc));
				
				q.first = (area * area) / (max_sq_length * second_max_sq_length); // (sine)*/
								
				if( q.sine() < this->B )
					return Mesh_2::BAD;
				
				else
					return Mesh_2::NOT_BAD;
					
			
			}
		};
		
		Is_bad is_bad_object() const
		{ return Is_bad(this->bound(), size_bound(), 
						this->traits /* from the bad class */); }
	};
	
} // end namespace CGAL

#endif

#ifndef STRETCHED_TRAITS_3
#define STRETCHED_TRAITS_3

/**Stretched_traits_3.h: 
  * This file defines the stretched geometrical traits of 3 dimension meshing.
*/

#include <CGAL/basic.h>
#include <CGAL/Robust_circumcenter_traits_3.h>

///The first parameter is a model of Kernel, and the second parameter is a 
/// model of a Metric_field.
template < typename K, typename Mfield>/**metric field*/
struct Stretched_traits_3 {

	public:


		/**field number type*/
		typedef typename K::FT                           FT;

		/**the objects required by TriangulationTriaits_3
		 * @ref cgal_manual/Triangulation_3_ref/Concept_TriangulationTraits_3.html
		 */
		typedef typename K::Point_3                      Point_3;
		typedef typename K::Line_3						 Line_3;
		typedef typename K::Ray_3	                     Ray_3;
		typedef typename K::Object_3                     Object_3;
		typedef typename K::Segment_3                    Segment_3;
		typedef typename K::Tetrahedron_3				 Tetrahedron_3;
		typedef typename K::Triangle_3                   Triangle_3;

		/**Predicates*/
		typedef typename K::Less_x_3                     Less_x_3;
		typedef typename K::Less_y_3                     Less_y_3;
		typedef typename K::Less_z_3					 Less_z_3;
		typedef typename K::Compare_x_3                  Compare_x_3;
		typedef typename K::Compare_y_3                  Compare_y_3;
		typedef typename K::Compare_z_3                  Compare_z_3;
		typedef typename K::Compare_xyz_3				 Compare_xyz_3;
		typedef typename K::Construct_object_3			 Construct_object_3;
		typedef typename K::Construct_segment_3			 Construct_segment_3;
		typedef typename K::Construct_ray_3				 Construct_ray_3;
		typedef typename K::Construct_equidistant_line_3 Construct_equidistant_line_3;

		/**decide for 4 vertex (p,q,r,s) the position of s correspond
		 * to the plane defined by (p,q,r) in a Mfield
		 * POSITIVE / NEGATIVE / COPLANAR
		 *
		 * Given a center fixed, the metric transformation does not change the orientation 
		 * of the points
		 * */

		/** Constructions : */
		
		Stretched_traits_3(unsigned _index, Point_3 _center, Mfield &_mf): 
			index(_index), center(_center), mf(_mf) {}
		
		Stretched_traits_3(const Stretched_traits_3  & t):
			index(t.get_index()), center(t.get_center()), mf(t.get_field()) {}
		
		Stretched_traits_3 &operator= (const Stretched_traits_3 &t) {
			index = t.get_index();
			center = t.get_center();
			mf = t.get_field();
			return *this;
		}



		/** see the side of oriented sphere:
		 *  positive/ negative/ on the bord
		 *
		 *  @return returns the relative position of point test 
		 *   to the oriented sphere defined by p, q, r and s. 
		 *   The order of the points p, q, r, and s is important, 
		 *   since it determines the orientation of the implicitly constructed sphere. 
		 *   If the points p, q, r and s  are positive oriented, positive 
		 *   side is the bounded interior of the sphere.
		 *  
		 * Precondition:  p, q, r, and s are coplanar and p, q, and r are not collinear.
		 */
		class Side_of_oriented_sphere_3 {
			typedef typename K::Point_3        Point_3;
			Mfield &mf;
			Point_3 center;
			public:
			typedef typename K::Oriented_side    result_type;

			Side_of_oriented_sphere_3(Point_3 _center, Mfield& _mf): 
				mf(_mf), center(_center){}

			result_type
				operator()( const Point_3& p, const Point_3& q,
						const Point_3& r, const Point_3& s,const Point_3& test) const { 
					typename K::Side_of_oriented_sphere_3 soos;
					return soos(mf.transfo(center, p),
							mf.transfo(center, q),
							mf.transfo(center, r),
							mf.transfo(center, s),
							mf.transfo(center, test));
				}
		};

		class Side_of_bounded_sphere_3 {
			typedef typename K::Point_3        Point_3;
			Mfield &mf;
			Point_3 center;
			public:
			typedef typename K::Bounded_side    result_type;

			Side_of_bounded_sphere_3(Point_3 _center, Mfield& _mf): 
				mf(_mf), center(_center){}
			
			result_type
				operator()( const Point_3& p, const Point_3& q,const Point_3& test) const { 
					typename K::Side_of_bounded_sphere_3 sobs;
					return sobs(mf.transfo(center, p),
							mf.transfo(center, q),
							mf.transfo(center, test));
				}
			
			result_type
				operator()( const Point_3& p, const Point_3& q,
						const Point_3& r,const Point_3& test) const { 
					typename K::Side_of_bounded_sphere_3 sobs;
					return sobs(mf.transfo(center, p),
							mf.transfo(center, q),
							mf.transfo(center, r),
							mf.transfo(center, test));
				}
			
			result_type
				operator()( const Point_3& p, const Point_3& q,
						const Point_3& r, const Point_3& s,const Point_3& test) const { 
					typename K::Side_of_bounded_sphere_3 sobs;
					return sobs(mf.transfo(center, p),
							mf.transfo(center, q),
							mf.transfo(center, r),
							mf.transfo(center, s),
							mf.transfo(center, test));
				}
		};


		class Transformation_metric_3 {
		  typedef typename K::Point_3		  Point_3;
		  Mfield &mf;
		  Point_3 center;
		  public:
		  typedef Point_3			  result_type;

		  Transformation_metric_3(Point_3 _center,Mfield& _mf):
			mf(_mf),center(_center){}

		  result_type
			operator() (const Point_3& p) const {
			  return mf.transfo(center,p);
			}
		};

			/** A predicate object that must provide the function operator
		 * Bounded_side operator()(Point p, Point q, Point r, Point s),
		 * which determines the bounded side of the circle defined by p, q, and r on which s lies.
		 * Precondition:  p, q, r, and s are coplanar and p, q, and r are not collinear.
		 */
		class Coplanar_orientation_3{
			typedef typename K::Point_3		  Point_3;
			Mfield &mf;
			Point_3 center;
			public:
			typedef typename K::Orientation	result_type;

			Coplanar_orientation_3(Point_3 _center,  Mfield& _mf):
				mf(_mf),center(_center){}

			result_type
				operator()( const Point_3& p, const Point_3&q,
						const Point_3& r) const {
					typename K::Coplanar_orientation_3 co;
					return co(mf.transfo(center,p),
							mf.transfo(center,q),
							mf.transfo(center,r));
				}//operator
		};

	
	
		/** see the side of oriented sphere:
		 *  positive/ negative/ on the bord
		 *
		 *  @return returns the relative position of point test 
		 *   to the oriented sphere defined by p, q, r and s. 
		 *   The order of the points p, q, r, and s is important, 
		 *   since it determines the orientation of the implicitly constructed sphere. 
		 *   If the points p, q, r and s  are positive oriented, positive 
		 *   side is the bounded interior of the sphere.
		 */  
		class Coplanar_side_of_bounded_circle_3{
			typedef typename K::Point_3		  Point_3;
			Mfield &mf;
			Point_3 center;
			public:
			typedef typename K::Bounded_side	result_type;

			Coplanar_side_of_bounded_circle_3(Point_3 _center,  Mfield& _mf):
				mf(_mf),center(_center){}

			result_type
				operator()( const Point_3& p, const Point_3&q,
						const Point_3& r, const Point_3& s) const {
					typename K::Coplanar_side_of_bounded_circle_3 csbc;
					return csbc(mf.transfo(center,p),
							mf.transfo(center,q),
							mf.transfo(center,r),
							mf.transfo(center,s));
				}//operator
		};
	
		/**	 A predicate object that must provide the function operator
		 * Orientation operator()(Point p, Point q, Point r, Point s),
		 * which returns POSITIVE, if s lies on the positive side of the oriented 
		 * plane h defined by p, q, and r, returns NEGATIVE if s lies on the 
		 * negative side of h, and returns COPLANAR if s lies on h.
		 */
		class Orientation_3{
			typedef typename K::Point_3		  Point_3;
			Mfield &mf;
			Point_3 center;
			public:
			typedef typename K::Orientation	result_type;

			Orientation_3(Point_3 _center, Mfield& _mf):
				mf(_mf),center(_center){}

			result_type
				operator()( const Point_3& p, const Point_3&q,
						const Point_3& r, const Point_3& s) const {
					typename K::Orientation_3 orient;
					return orient(mf.transfo(center,p),
						mf.transfo(center,q),
						mf.transfo(center,r)
						,mf.transfo(center,s));
				}//operator
		};

	
		/** A predicate object that must provide the function operator
		 * Comparison_result operator()(Point p, Point q, Point r),
		 * which compares the distance between p and q to the distance between p and r.
		 * @return enum{SMALL, EQUAL, LARGE}
		 * */
		class Compare_distance_3 {
			typedef typename K::Point_3   Point_3;
			Mfield &mf;
			Point_3 center; 
			public:
			typedef typename K::Comparison_result     result_type;

			Compare_distance_3(Point_3 _center, Mfield& _mf): 
				mf(_mf),center(_center){}

			result_type
				operator()(const Point_3& p, const Point_3& q, const Point_3& r) const { 
					typename K::Compare_distance_3 cd;
					return cd(mf.transfo(center, p),
							mf.transfo(center, q),
							mf.transfo(center, r));
				}
		};

		/** Computes the stretched distance between two points
		 *
		 * */
		class Compute_squared_distance_3 {
			typedef typename K::Point_3   Point_3;
			Mfield &mf;
			Point_3 center;
			public:
			typedef typename K::FT			    result_type;

			Compute_squared_distance_3(Point_3 _center, Mfield& _mf):
				mf(_mf), center(_center){}

			result_type
				operator()(const Point_3& p, const Point_3& q) const {
					typename K::Compute_squared_distance_3 csd;
					return csd(mf.transfo(center, p),
							mf.transfo(center, q));
				}
		};

		/** Decide whether in triple (p,q,r), r encroaches the edge (p,q)
		 *  with st provided distance, this function can return stretched 
		 *  judgement of encroachment
		 */
		class Encroachment_edge_3 {
			typedef typename K::Point_3   Point_3;
			const Stretched_traits_3<K,Mfield>* st;
			public:
			typedef bool                        result_type;

			Encroachment_edge_3(const Stretched_traits_3<K,Mfield>* _st):
				st(_st){}

			result_type
				operator()(const Point_3& p, const Point_3& q, const Point_3& r) const {
					Compute_squared_distance_3 csd = st->compute_squared_distance_3_object();
					return csd(p, q) > csd(p, r) + csd(q, r);
				}
		};

		/** Decide where test encroaches (p,q,r)
		 * i.e. if c is the circumcenter of (p,q,r), r the circumradius
		 * then return bool distance(c,test)<r
		 * */
		class Encroachment_face_3 {
			typedef typename K::Point_3	  Point_3;
			const Stretched_traits_3<K,Mfield>* st;
			public:
			typedef bool						result_type;

			Encroachment_face_3(const Stretched_traits_3<K,Mfield>* _st):
				st(_st){}

			result_type
				operator()(const Point_3 &p,const Point_3& q,
						const Point_3& r,const Point_3& test) const {
					Construct_circumcenter_3 cc=st->construct_circumcenter_3_object();
					Compare_distance_3 cd=st->compare_distance_3_object();
					Point_3 c=cc(p,q,r);
					return cd(c,p,test)==CGAL::LARGER;
				}
		};//class Encroachment_face_3


		/**construct the circumcenter according to 4 point of input
		*/
		class Construct_circumcenter_3 {
			typedef typename K::Point_3   Point_3;
			Mfield &mf;
			Point_3 center;
			public:
			typedef Point_3                     result_type;
			typename CGAL::Robust_construct_circumcenter_3<K> cc;

			Construct_circumcenter_3(Point_3 _center, Mfield& _mf):
				mf(_mf), center(_center){}

			result_type
				operator()(const Point_3& p, const Point_3& q, 
						const Point_3& r, const Point_3& s) const {
					return mf.transfo_inv(center, cc(mf.transfo(center, p),
								mf.transfo(center, q),
								mf.transfo(center, r),
								mf.transfo(center, s)));
				}

			result_type
				operator()(const Point_3& p, const Point_3& q, const Point_3& r) const {
					return mf.transfo_inv(center, cc(mf.transfo(center, p),
								mf.transfo(center, q),
								mf.transfo(center, r)));
				}
		};

		/** Given 4 points, Conpute the distortion of a tretahedron
		 * i.e. compute the maximum distortion of each pairs of the 4 vertice
		 * */
		class Compute_distortion_3 {
			typedef typename K::Point_3 Point_3;
			Mfield &mf;
			public:
			typedef FT					  result_type;

			Compute_distortion_3(Mfield& _mf):
				mf(_mf){}

			private:
			result_type
				compute_distortion_3(const Point_3 &p,const Point_3 &q,const Point_3&r)const{
					FT max;
					FT distortion;
				    max = distortion = mf.Compute_distortion(p,q);
					distortion = mf.Compute_distortion(p,r);
					if(distortion>max){max=distortion;}
					distortion = mf.Compute_distortion(q,r);
					if(distortion>max){max=distortion;}
					return max;
				}

			public:
			result_type
				operator()(const Point_3& p,const Point_3& q,const Point_3& r) const{
					return compute_distortion_3(p,q,r);
				}
	
			result_type
				operator()(const Point_3& p,const Point_3& q,const Point_3& r,const Point_3 &s) const{
					FT max;
					FT distortion;
				    max = distortion = mf.Compute_distortion(p,q);
					distortion = compute_distortion_3(p,r,s);
					if(distortion>max){max=distortion;}
					distortion = compute_distortion_3(q,r,s);
					if(distortion>max){max=distortion;}
					return max;
				}

		};//class Compute_distortion_3


		FT Absolute_distortion() const{
			return mf.compute_absolute_distortion(center);
		}

		Less_x_3
			less_x_3_object() const
			{ return Less_x_3();}

		Less_y_3
			less_y_3_object() const
			{ return Less_y_3();}

		Less_z_3
			less_z_3_object() const
			{ return Less_z_3();}

		Compare_x_3
			compare_x_3_object() const
			{ return Compare_x_3();}

		Compare_y_3
			compare_y_3_object() const
			{ return Compare_y_3();}

		Compare_z_3
			compare_z_3_object() const
			{ return Compare_z_3();}

		Compare_xyz_3
		    compare_xyz_3_object() const
			{ return Compare_xyz_3();}

		Orientation_3
			orientation_3_object() const
			{ return Orientation_3(center, mf);}
	
		Coplanar_side_of_bounded_circle_3
			coplanar_side_of_bounded_circle_3_object() const
			{ return Coplanar_side_of_bounded_circle_3(center,mf);}

		Coplanar_orientation_3
			coplanar_orientation_3_object() const
			{ return Coplanar_orientation_3(center,mf);}

		Transformation_metric_3
		    transformation_metric_3_object() const
			{ return Transformation_metric_3(center,mf);}

		Encroachment_edge_3
			encroachment_edge_3_object() const
			{ return Encroachment_edge_3(this);}

		Encroachment_face_3
			encroachment_face_3_object() const
			{ return Encroachment_face_3(this);}

		Side_of_oriented_sphere_3
			side_of_oriented_sphere_3_object() const
			{return Side_of_oriented_sphere_3(center, mf);}

		
		Side_of_bounded_sphere_3
			side_of_bounded_sphere_3_object() const
			{return Side_of_bounded_sphere_3(center,mf);}

		Construct_circumcenter_3
			construct_circumcenter_3_object() const
			{ return Construct_circumcenter_3(center,mf);}

		Compute_distortion_3
			compute_distortion_3_object() const
			{ return Compute_distortion_3(mf);}

		Compute_squared_distance_3
			compute_squared_distance_3_object() const
			{return Compute_squared_distance_3(center, mf);}

		Compare_distance_3
			compare_distance_3_object() const
			{return Compare_distance_3(center, mf);}

		Construct_object_3
			construct_object_3_object() const
			{return Construct_object_3();}

		Construct_segment_3
			construct_segment_3_object() const
			{return Construct_segment_3();}

		Construct_ray_3
			construct_ray_3_object() const
			{return Construct_ray_3();}

		Construct_equidistant_line_3
			construct_equidistant_line_3_object() const
			{return Construct_equidistant_line_3();}

		unsigned get_index() const {
			return index;
		}

		Point_3 get_center() const {
			return center;
		}

		Mfield & get_field() const {
			return mf;
		}

	private:
		unsigned index;		///< the index of corresponding index center

		Point_3 center;	///< center for metric

		Mfield  &mf;		///< metric field
};

#endif//STRETCHED_TRAITS_3

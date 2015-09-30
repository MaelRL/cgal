#ifndef SURFACE_CRITERIA
#define SURFACE_CRITERIA

/** Criteria of surface
 *
 * Surfaces (triangles) are judged by the following criteria:
 * 1. size
 * 2. distortion
 * 3. coherence
 * 4. distance to the given surface
 *
 * Triangles will be judged and scored according to their "quality" of each item.
 * 
 * All triangles have scores greater than 0.
 * 
 * The Umbrella_set is required to have the concept of type of 
 *  - Traits: geometrical traits
 *  - Constrain_surface:  
 *
 *
 * */


template<class Umbrella_set,class Constrain_surface>
class Surface_criteria_for_umbrella_set{

	public:

	typedef typename Umbrella_set::Traits				  Traits;
	typedef typename Traits::FT							  FT;
	typedef typename Traits::Point_3					  Point_3;
	typedef typename Traits::Line_3                       Line_3;
	typedef typename Traits::Ray_3                        Ray_3;
	typedef typename Traits::Object_3                     Object_3;
	typedef typename Traits::Segment_3                    Segment_3;
	typedef typename Traits::Tetrahedron_3                Tetrahedron_3;
	typedef typename Traits::Triangle_3                   Triangle_3;
	typedef typename Traits::Construct_circumcenter_3	  Construct_circumcenter_3;
	typedef typename Traits::Compute_squared_distance_3	  Compute_squared_distance_3;
	typedef typename Traits::Compute_distortion_3		  Compute_distortion_3;
	typedef typename Traits::Transformation_metric_3	  Transformation_metric_3;

	Surface_criteria_for_umbrella_set(const Constrain_surface *cs, FT max_ratio,
			FT _max_distortion, FT max_size, FT max_surface_distance):
		constrain_surface(cs),
		max_squared_ratio(max_ratio * max_ratio),
		max_distortion(_max_distortion),
		max_squared_size(max_size * max_size),
		max_squared_surface_distance(max_surface_distance * max_surface_distance)
	{}

	bool is_bad(const Point_3 &a,const Point_3 &b, const Point_3 &c,const Traits &traits)const {
	
		FT over_size = squared_size(a,b,c,traits) - max_squared_size;
		if(over_size > 0.){
#ifdef DEBUG_SURFACE_CRITERIA
			std::cout<<"oversized:"<<over_size<<std::endl;
#endif
			return true;
		}
		FT over_distortion = distortion(a,b,c,traits) - max_distortion;
		if(over_distortion > 0.){
#ifdef DEBUG_SURFACE_CRITERIA
			std::cout<<"over_distortioned:"<<over_distortion<<std::endl;
#endif
		}
		FT over_ratio = squared_ratio(a,b,c,traits) - max_squared_ratio;
		if(over_ratio > 0.){
#ifdef DEBUG_SURFACE_CRITERIA
			std::cout<<"over_ratio:"<<over_ratio<<std::endl;
#endif
		}
		FT over_surface_distance = squared_surface_distace(a,b,c,traits) - max_squared_surface_distance;
		if(over_surface_distance>0){
#ifdef DEBUG_SURFACE_CRITERIA
			std::cout<<"over_surface_distance:"<<over_surface_distance<<std::endl;
#endif
		}

		return false;
	}

	FT score(const Point_3&a, const Point_3 &b, const Point_3 &c,const Traits &traits) const{
	  
		FT size_rate = squared_size(a,b,c,traits) / max_squared_size;
		FT distortion_rate = distortion(a,b,c,traits) / max_distortion;
		FT ratio_rate = squared_ratio(a,b,c,traits) / max_squared_ratio;
		FT surface_distance_rate = squared_surface_distace(a,b,c,traits) / max_squared_surface_distance;

		FT total_score = size_rate + distortion_rate + ratio_rate + surface_distance_rate;

		return total_score;
	}

//==================================================================================================
//								  The calculation functions
//==================================================================================================

	/** The squared size in space of transformed metric:
	 *	  S = ( |ab|^2 * |ac|^2 - (ab.ac)^2 )/4;
	 */
	FT squared_size(const Point_3 &a,const Point_3 &b,const Point_3 &c,
			const Traits &traits)const{
		
		Compute_squared_distance_3 csd = traits.compute_squared_distance_3_object();
		Transformation_metric_3 transfo = traits.transformation_metric_3_object();
		Point_3 ta=transfo(a);
		Point_3 tb=transfo(b);
		Point_3 tc=transfo(c);
		FT x[3];
		FT y[3];
		x[0]=tb.x()-ta.x(); y[0]=tc.x()-ta.x();
		x[1]=tb.y()-ta.y(); y[1]=tc.y()-ta.y();
		x[2]=tb.z()-ta.z(); y[2]=tc.z()-ta.z();
		FT product = 0;
		for(int i=0;i<3;i++){
			product += x[i]*y[i];
		}

		return (csd(a,b)*csd(a,c)-product*product)/4;
	}


	FT distortion(const Point_3 &a,const Point_3&b, const Point_3 &c,
			const Traits &traits)const{
		Compute_distortion_3  cd = traits.compute_distortion_3_object();
		return cd(a,b,c);
	}
	
	FT squared_ratio(const Point_3 &p,const Point_3 &q,const Point_3 &r,
			const Traits &traits)const{

		Construct_circumcenter_3  cc = traits.construct_circumcenter_3_object();
		Compute_squared_distance_3 csd = traits.compute_squared_distance_3_object();
		Point_3 center = cc(p,q,r);
		FT radius_square = csd(p,center);
		FT max_squared_distance = csd(p,q);
		FT squared_distance  = csd(p,r);
		if(max_squared_distance < squared_distance){
			max_squared_distance = squared_distance;
		}
		squared_distance  = csd(q,r);
		if(max_squared_distance < squared_distance){
			max_squared_distance = squared_distance;
		}
		return radius_square/max_squared_distance;
	}

	//TODO:  Questions:
	// 1. the distance between surface and the barycenter or the circumcenter of the triangle ?
	// 2. the projected point or the intersection point
	FT squared_surface_distace(const Point_3 &p,const Point_3 &q,const Point_3 &r, 
			const Traits &traits) const{
		Point_3 barycenter((p.x()+q.x()+r.x())/3,(p.y()+q.y()+r.y())/3,(p.z()+q.z()+r.z())/3);
		Point_3 projected_center = constrain_surface->projection(barycenter);
		Compute_squared_distance_3 csd = traits.compute_squared_distance_3_object();
		return csd(barycenter,projected_center);
	}

	private:

	const Constrain_surface	*constrain_surface;
	FT	max_squared_ratio;
	FT	max_distortion;
	FT	max_squared_size;
	FT	max_squared_surface_distance;
};


#endif //SURFACE_CRITERIA

#ifndef VOLUME_CRITERIA
#define VOLUME_CRITERIA

/** Criteria of surface
 *
 * Surfaces (triangles) are judged by the following criteria:
 * 1. volume
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


template<class Umbrella_set>
class Volume_criteria_for_umbrella_set{

	public:

	typedef typename Umbrella_set::Traits				  Traits;
	typedef typename Umbrella_set::Cell_handle			  Cell_handle;
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

	Volume_criteria_for_umbrella_set(FT max_ratio,FT _max_distortion, FT _max_volume):
		max_squared_ratio(max_ratio * max_ratio),
		max_distortion(_max_distortion),
		max_volume(_max_volume)
	{}

	bool is_bad(const Cell_handle &cell,const Traits &traits)const {
	
		FT over_volume = volume(cell,traits) - max_volume;
		if(over_volume > 0.){
#ifdef DEBUG_VOLUME_CRITERIA
			std::cout<<"over_volumed:"<<over_volume<<std::endl;
#endif
			return true;
		}
		FT over_distortion = distortion(cell,traits) - max_distortion;
		if(over_distortion > 0.){
#ifdef DEBUG_VOLUME_CRITERIA
			std::cout<<"over_distortioned:"<<over_distortion<<std::endl;
#endif
		}
		FT over_ratio = squared_ratio(cell,traits) - max_squared_ratio;
		if(over_ratio > 0.){
#ifdef DEBUG_VOLUME_CRITERIA
			std::cout<<"over_ratio:"<<over_ratio<<std::endl;
#endif
		}

		return false;
	}

	FT score(const Cell_handle	&cell,const Traits &traits) const{
	  
		FT volume_rate = volume(cell,traits) / max_volume;
		FT distortion_rate = distortion(cell,traits) / max_distortion;
		FT ratio_rate = squared_ratio(cell,traits) / max_squared_ratio;

		FT total_score = volume_rate + distortion_rate + ratio_rate ;

		return total_score;
	}

//==================================================================================================
//								  The calculation functions
//==================================================================================================

	/** The squared volume in space of transformed metric:
	 *	  S = ( |ab|^2 * |ac|^2 - (ab.ac)^2 )/4;
	 */
	FT volume(const Cell_handle & cell, const Traits &traits)const{
		
		Transformation_metric_3 transfo = traits.transformation_metric_3_object();
		Point_3 tp[4];
		for(int i=0; i<4; i++){
			tp[i] = transfo(cell->vertex(i)->point());
		}
		FT x[3];
		FT y[3];
		FT z[3];
		for(int i=1;i<3;i++){
			x[i]=tp[i].x()-tp[0].x();
			y[i]=tp[i].y()-tp[0].y();
			z[i]=tp[i].z()-tp[0].z();
		}
		FT determinant = x[0]*y[1]*z[2] + x[1]*y[2]*z[0] + x[2]*y[0]*z[1] 
			- x[0]*y[2]*z[1] - x[2]*y[1]*z[0] - x[1]*y[0]*z[2];

		return fabs(determinant/6);
	}


	FT distortion(const Cell_handle & cell,const Traits &traits)const{
		Compute_distortion_3  cd = traits.compute_distortion_3_object();
		return cd(cell->vertex(0)->point(),
			cell->vertex(1)->point(),
			cell->vertex(2)->point(),
			cell->vertex(3)->point());
	}
	
	FT squared_ratio(const Cell_handle & cell,const Traits &traits)const{

		Compute_squared_distance_3 csd = traits.compute_squared_distance_3_object();
		Point_3 p[4];
		for(int i=0;i<4;i++){
			p[i]=cell->vertex(i)->point();
		}
		Point_3 center = cell->circumcenter(traits);
		FT radius_square = csd(p[0],center);
		FT max_squared_distance = csd(p[0],p[1]);
		FT squared_distance;
		for(int i=2;i<4;i++){
			for(int j=0;j<i;j++){
				squared_distance = csd(p[j],p[i]);
				if(max_squared_distance < squared_distance){
					max_squared_distance = squared_distance;
				}
			}
		}
		return radius_square/max_squared_distance;
	}

	private:

	FT	max_squared_ratio;
	FT	max_distortion;
	FT	max_volume;
};


#endif //VOLUME_CRITERIA

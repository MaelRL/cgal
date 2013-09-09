#ifndef CIRCULAR_METRIC_FIELD_3
#define CIRCULAR_METRIC_FIELD_3

#include <CGAL/Metric_field_3.h>	//extended incl


/** Defines the interface of anisotopic
* geometrical metric of 3 dimensional geometry
*/
#include<cmath>
#include<cassert>

#define PI (2*acos(0.))
//================================================================================
//      Circular metric field: (along x-z plane)
//================================================================================


template <typename Point_3,typename FT=double>
class Circular_metric_field_3 : public Metric_field_spherical <FT, Point_3> {

	public:

	Circular_metric_field_3():stretch(1){
#ifndef NDEBUG
		std::cerr<<"\n\n\n\n\na mf of stretch 1 is created here."<<std::endl;
#endif
	}
	Circular_metric_field_3(FT _stretch):stretch(_stretch){
	}
	Circular_metric_field_3(const Circular_metric_field_3 & c):stretch(c.get_stretch()){
	}

	Circular_metric_field_3 & operator=(const Circular_metric_field_3 &c){
		stretch = c.get_stretch();
		return *this;
	}

	FT get_theta1(const Point_3 &p){
		return atan2(p.z(),p.x());
	}//get_theta1

	FT get_phi1(const Point_3 &p){
		return PI;
	}//get_phi1

	~Circular_metric_field_3(){}

	FT get_theta2(const Point_3 &){return acos(0.);}
	FT get_val1(const Point_3 &){return stretch;}
	FT get_val2(const Point_3 &){return 1.;}
	FT get_val3(const Point_3 &){return 1.;}

	FT get_stretch()const {
		return stretch;
	}

	private:

	FT stretch;
};

#undef PI	//3.1415925357

#endif//CIRCULAR_METRIC_FIELD_3


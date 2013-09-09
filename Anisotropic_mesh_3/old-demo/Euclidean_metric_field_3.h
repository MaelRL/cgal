
#ifndef EUCLIDEAN_METRIC_FIELD_3
#define EUCLIDEAN_METRIC_FIELD_3

#include <CGAL/Metric_field_3.h>

//================================================================================
//      Eucliean metric field
//================================================================================

//Standard Euclidean metric field with character
// 1 0
// 0 1
// default metric field
template < typename Point_3 >
struct Euclidean_metric_field_3 :public Metric_field_spherical < double, Point_3 > {
  	double get_theta1(const Point_3 &) {return acos(0.);};
 	double get_phi1(const Point_3 &) {return 0.;};
 	double get_theta2(const Point_3 &) {return acos(0.);};
	double get_val1(const Point_3&){return 1.;}
	double get_val2(const Point_3&){return 1.;}
	double get_val3(const Point_3&){return 1.;}

};

#endif//EUCLIDEAN_METRIC_FIELD_3


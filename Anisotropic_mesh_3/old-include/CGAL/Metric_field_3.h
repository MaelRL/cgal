
#ifndef METRIC_FIELD_3
#define METRIC_FIELD_3

/** Defines the interface of anisotopic
* geometrical metric of 3 dimensional geometry
*/
#include<cmath>
#include<cassert>

#define PI (2*acos(0.))

//================================================================================
//      spherical metric field(base type)
//================================================================================


/**FT field type*/
template < typename FT, typename Point_3 >
class Metric_field_spherical {
	/**Metric defined by 3 eigen value and its corresponding vector
	*The vectors must be orthogonal, so:
	*--- the first unit vector will be initialized by phi and theta in 
	*spherical coordinates. the second vector will only be specified its
	*theta, the other parameters can be deducted.
	*
	*      ref (spherical_coordinates):
	*http://en.wikipedia.org/wiki/Spherical_coordinates 
	*/

	public:

	///TODO maybe should copy all the matrice?
	Metric_field_spherical(){}

	virtual FT get_theta1(const Point_3&) = 0;
	virtual FT get_phi1  (const Point_3&) =0;
	virtual FT get_theta2(const Point_3&) =0;

	/** eigenvalue*/
	virtual FT get_val1(const Point_3&) = 0;
	virtual FT get_val2(const Point_3&) = 0;
	virtual FT get_val3(const Point_3&) = 0;

	virtual ~Metric_field_spherical(){};

	private:

	/**given theta1 phi1 and theta2,
	*solve the orthogonal angle of phi2 in spherical
	*coordination
	*/
	FT solve_phi2(){
		assert((theta1>=0)&&(theta2>=0));

		FT c=cos(theta1)*cos(theta2);
		FT s=sin(theta1)*sin(theta2);

		if(fabs(c)+fabs(s)<1e-10){
			return phi1+acos(0.);
		}
		if((c>s) || (c<-s)){  
			// (-1 > c/s) or (1 < c/s)
#ifndef NDEBUG_METRIC_VALIDITY 
			std::cerr<<"warning: theta2 value invalid:"
				<<c<<":"<<s<<" \
				resetting theta2 to pi/2."<<std::endl;
#endif
			theta2 = acos(0.);
			return phi1+acos(0.);
		}
		else{
			FT v=-c/s;
			assert((v>=-1)&&(v<=1));
			return phi1+acos(v);
		}
	}

	/**Given the unit spherical coordination of 
	//the first two eigen vectros theta1(t1) phi1(p1) theta2(t2) phi2(p2)
	//find out the orthogonal matrix \f$O\f$
	//
	\f{array}{cccccccccccccc}
    O & = & [& o00& o01& o02& ]& = & [& sin(t_1)cos(p_1)& sin(t_1)sin(p_1)& cos(t_1)& ]\\
	    & & [& o10& o11& o12& ]&   & [& sin(t_2)cos(p_2)& sin(t_2)sin(p_2)& cos(t_2)& ]\\
	    & & [& o20& o21& o22& ]&   & [&    ?           &     ?            &    ?    & ]
	* \f
	*/
	void find_orthogonal(){
		o[0][0] = sin(theta1)*cos(phi1);
		o[0][1] = sin(theta1)*sin(phi1);
		o[0][2] = cos(theta1);

		o[1][0] = sin(theta2)*cos(phi2);
		o[1][1] = sin(theta2)*sin(phi2);
		o[1][2] = cos(theta2);

		o[2][0] = o[1][2] * o[0][1] - o[1][1] * o[0][2];
		o[2][1] = o[1][0] * o[0][2] - o[0][0] * o[1][2];
		o[2][2] = o[0][0] * o[1][1] - o[1][0] * o[0][1];

		FT f=sqrt(o[2][0]*o[2][0] + o[2][1]*o[2][1]+o[2][2]*o[2][2]);
		
		for(int i=0;i<3;i++){
			o[2][i] /= f;
		}//for

	}//find_orthogonal

	/** set center(origin) of metric
	*calculate phi2 accordingly
	*   and calculate the square root of eigenvalue
	*   and then the square root of matrix
	*/
	void set_center(const Point_3 &center){
		// 0 <= theta1 < PI
		// 0 <=  phi1  < PI*2
		/*if(center == current_center){
		  return;
		}else{
		  current_center = center;
		}*/
		theta1=	fmod(fmod(get_theta1(center),PI)+PI,PI);
		phi1  =	fmod(fmod(get_phi1(center),PI*2)+PI*2,PI*2);
		theta2=	fmod(fmod(get_theta2(center),PI)+PI,PI);
		phi2  = solve_phi2();

		// find the orthogonal matrix accordingly (o[0][0]-o[2][2])
		find_orthogonal();
		calculate_sqrt_matrix(center);
		calculate_inv_sqrt_matrix(center);
	
	}//set_center

	private:

	void print_matrix(FT m[3][3]){
		std::cout<<"******************"<<std::endl;
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				if(fabs(m[i][j])<1e-10){
					std::cout<<0<<"\t";
				}else{
					std::cout<<m[i][j]<<"\t";
				}//if ~ 0
			}//for
			std::cout<<std::endl;
		}//for
	}//print_matrix

	void calculate_inv_sqrt_matrix(const Point_3 &center){

		FT is_v[3];//inverse sqrt eigenvalue
#ifndef NDEBUG_METRIC_VALIDITY
		if(get_val1(center)<=1e-9 || get_val2(center)<1e-9 || 
				get_val3(center)<=1e-9){
			std::cerr<<"eigenvalue positive required!"<<std::endl;
		}
#endif
		is_v[0]=1/sqrt(fabs(get_val1(center)));
		is_v[1]=1/sqrt(fabs(get_val2(center)));
		is_v[2]=1/sqrt(fabs(get_val3(center)));

		orth_eigen_multi(is_v,m_i);

		//std::cout<<"M_I ";
		//print_matrix(m_i);

	}//calculate_inv_sqrt_matrix


	void calculate_sqrt_matrix(const Point_3 &center){

		FT sqrt_v[3];//sqrt eigenvalue
#ifndef NDEBUG_METRIC_VALIDITY
		if(get_val1(center)<0 || get_val2(center)<0 || get_val3(center)<0){
			std::cerr<<"eigenvalue positive required!"<<std::endl;
		}
#endif
		sqrt_v[0]=sqrt(fabs(get_val1(center)));
		sqrt_v[1]=sqrt(fabs(get_val2(center)));
		sqrt_v[2]=sqrt(fabs(get_val3(center)));
		
		orth_eigen_multi(sqrt_v,m);

		//std::cout<<"M ";
		//print_matrix(m);
	}//calculate_sqrt_matrix

	//Calcul mat = O' * dia(v) * O
	void orth_eigen_multi(FT v[3],FT mat[3][3]){

		FT temp[3][3];
		for(int i=0 ;i<3 ;i++){
			for(int j=0; j<3; j++){
				temp[i][j] = o[i][j] * v[i];
			}//for
		}//for

		multi_orth(temp,mat);
	
	}//orth_eigen_multi

	//calculate mat = O' * right
    void multi_orth(FT right[3][3],FT mat[3][3]){
		for(int i=0;i<3;i++){
			for(int j=0; j<3; j++){
				mat[i][j] =0;
				for(int k=0;k<3;k++){
					//Multiply transverse of O
					mat[i][j] += o[k][i] * right[k][j];
				}//for
			}//for
		}//for	
	}//multi_orth

	public:
	/** transfer p's axis according to the center point metric
	*/
	Point_3 transfo(Point_3 center, Point_3 p) {
		set_center(center);
		return Point_3(m[0][0]*p.x() + m[0][1]*p.y() + m[0][2]*p.z(),
				m[1][0]*p.x() + m[1][1]*p.y() + m[1][2]*p.z(),
				m[2][0]*p.x() + m[2][1]*p.y() + m[2][2]*p.z());
	}

	/** inverse transfer of p's axis according to center point metric
	*/
	Point_3 transfo_inv(Point_3 center, Point_3 p) {
		set_center(center);
		return Point_3(m_i[0][0]*p.x() + m_i[0][1]*p.y() + m_i[0][2]*p.z(),
				m_i[1][0]*p.x() + m_i[1][1]*p.y() + m_i[1][2]*p.z(),
				m_i[2][0]*p.x() + m_i[2][1]*p.y() + m_i[2][2]*p.z());
	}

	/**Distortion of the two metrics of two points are defined as
	 * \f$ max(\|M_xM^{-1}_y)\|_1,\|(M_yM^{-1}_x)\|_1)
	 * \f$ 
	 * */
	FT Compute_distortion(const Point_3 &px,const Point_3 &py) {
	
		//calculate and record the metric matrice of x
		set_center(px);

		FT mx[3][3];
		FT mx_i[3][3];
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				mx[i][j]=m[i][j];
				mx_i[i][j]=m_i[i][j];
			}
		}//for

		//calculate metric matrice of y
		set_center(py);

		FT d1,d2;		//recorder of elements of multiplied matrix
		FT max_sum = 0.;
		FT sum1,sum2;
		for(int i=0;i<3;i++){
			sum1=sum2=0;
			for(int j=0;j<3;j++){
				d1=d2=0;
				for(int k=0;k<3;k++){
					d1+=mx[i][k]*m_i[k][j];
					d2+=mx_i[i][k]*m[k][j];
				}//for
				sum1+=fabs(d1);
				sum2+=fabs(d2);
			}
			if(sum1>max_sum)
				max_sum=sum1;
			if(sum2>max_sum)
				max_sum=sum2;
		}//for
		return max_sum;
	}

	FT compute_absolute_distortion(const Point_3 &p){

		FT sum1,sum2;
		FT max_sum = 0;

		set_center(p);

		for(int i=0;i<3;i++){
			sum1 = sum2 =0;
			for(int j=0;j<3;j++){
				sum1 += fabs(m[i][j]);
				sum2 += fabs(m_i[i][j]);
			}
			if(sum1>max_sum)
				max_sum=sum1;
			if(sum2>max_sum)
				max_sum=sum2;
		}//for
		return max_sum;
	}

	private:

	Point_3	  current_center;

	//spherical coordination 
	FT theta1;
	FT phi1;
	FT theta2;
	FT phi2;

	//The 9 elements of the orthogonal matrix
	FT o[3][3];
	
	//The 9 elements of squre root matrix:
	FT m[3][3];

	//The 9 elements of the inverse of square root matrix:
	FT m_i[3][3];

};

//================================================================================
//      Eucliean metric field
//================================================================================

//Standard Euclidean metric field with character
// 1 0
// 0 1
// default metric field
template < typename Point_3 >
struct Euclidean_mf_3 :public Metric_field_spherical < double, Point_3 > {
  	double get_theta1(const Point_3 &) {return acos(0.);};
 	double get_phi1(const Point_3 &) {return 0.;};
 	double get_theta2(const Point_3 &) {return acos(0.);};
	double get_val1(const Point_3&){return 1.;}
	double get_val2(const Point_3&){return 1.;}
	double get_val3(const Point_3&){return 1.;}

};

//================================================================================
//      Circular metric field: (along x-z plane)
//================================================================================


template <typename Point_3>
class Circular_mf_3 : public Metric_field_spherical <double, Point_3> {

	public:

	Circular_mf_3():stretch(1){
#ifndef NDEBUG_METRIC_VALIDITY
		std::cerr<<"\n\n\n\n\na mf of stretch 1 is created here."<<std::endl;
#endif
	}
	Circular_mf_3(double _stretch):stretch(_stretch){
	}
	Circular_mf_3(const Circular_mf_3 & c):stretch(c.get_stretch()){
	}

	Circular_mf_3 & operator=(const Circular_mf_3 &c){
		stretch = c.get_stretch();
		return *this;
	}

	double get_theta1(const Point_3 &p){
		return atan2(p.z(),p.x());
	}//get_theta1

	double get_phi1(const Point_3 &p){
		return PI;
	}//get_phi1

	~Circular_mf_3(){}

	double get_theta2(const Point_3 &){return acos(0.);}
	double get_val1(const Point_3 &){return stretch;}
	double get_val2(const Point_3 &){return 1.;}
	double get_val3(const Point_3 &){return 1.;}

	double get_stretch()const {
		return stretch;
	}

	private:

	double stretch;
};

#undef PI	//3.1415925357

#endif//METRIC_FIELD_3


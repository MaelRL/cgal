#ifndef CGAL_ANISOTROPIC_MESH_3_ELLIPSOID_METRIC_FIELD
#define CGAL_ANISOTROPIC_MESH_3_ELLIPSOID_METRIC_FIELD

#include <CGAL/Metric_field.h>

using namespace CGAL::Anisotropic_mesh_3;


template<typename K>
class Ellipsoid_metric_field : public Metric_field<K> {
public:
  typedef Metric_field<K>        Base;
  typedef typename Base::FT      FT;
  typedef typename Base::Metric  Metric;
  typedef typename Base::Point_3 Point_3;
  typedef typename Base::Vector_3 Vector_3;

public:
  FT a, b, c;

public:
  virtual void report(typename std::ofstream &fx) const {
    fx << "type:   ellipsoid" << std::endl;
    fx << "a, b, c :  " << a << " " << b << " " << c << std::endl;
    fx << "epsilon: " << this->epsilon << std::endl;
  }

  virtual Metric compute_metric(const Point_3 &p) const
  {
    FT x = p.x(), y = p.y(), z = p.z();
    FT R = std::sqrt(x*x + y*y + z*z);
    FT sq_a = a*a, sq_b = b*b, sq_c = c*c;

    Vector_3 g(x/(a*a),y/(b*b),z/(c*c));
    FT n_norm = std::sqrt( x*x/(sq_a*sq_a) + y*y/(sq_b*sq_b) + z*z/(sq_c*sq_c) );
    FT h = 1./n_norm;
    Vector_3 n(g/n_norm);

    //eigenvalues
    FT sq_h = h*h;
    FT denom = 1.0/(2.0 * sq_a * sq_b * sq_c);
    FT delta_1 = sq_h * h * (sq_a + sq_b + sq_c - (R * R) );
    FT delta_2 = 4.0 * sq_h * sq_h * sq_a * sq_b * sq_c;
    FT delta = delta_1 * delta_1 - delta_2;
    delta = std::sqrt(delta);
    FT e1 = denom * (delta + delta_1);
    FT e2 = denom * (-delta + delta_1);

    FT en = a * c/(sq_b * b); //curvature at the umbilic == max_curv ? probably wrong

    //eigenvectors

    // ----------------------------------------------------------
    // WARNING : ONLY CORRECT FOR PROLATE SPHEROIDS (A > B = C)
        //normal at P in the ellipse formed by the inter of yOz and the ellipsoid
    FT v1_norm = std::sqrt( z*z/(sq_c*sq_c) + y*y/(sq_b*sq_b) );
    Vector_3 v1(0,-z/(sq_c*v1_norm), y/(sq_b*v1_norm));
    Vector_3 v2 = CGAL::cross_product(n, v1);
    // ----------------------------------------------------------

    return this->build_metric(n, v1, v2, en, e1, e2);
  }

  Ellipsoid_metric_field(FT a_, FT b_, FT c_, FT epsilon_ = 1e-3)
    :Metric_field<K>(epsilon_), a(a_), b(b_), c(c_){ }
};


#endif

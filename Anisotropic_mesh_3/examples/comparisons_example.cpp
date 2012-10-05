#include <CGAL/Timer.h>
#include <CGAL/Default_configuration.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Euclidean_metric_field.h>
#include <CGAL/Implicit_curvature_metric_field.h>
#include <Domain/Constrain_surface_3_torus.h>

#include <Metric_field/Cylinder_metric_field.h>
#include <Domain/Constrain_surface_3_cylinder.h>


#include <CGAL/Anisotropic_surface_mesher_3.h>


using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel	K;
typedef CGAL::Timer Timer;


template<typename PlaneType, typename PointType>
struct Is_between
{
public:
  bool operator()(const PointType& p) const
  {
    if(xcondition && p.x() < 0) 
      return false;
    PointType pp1 = plane1.projection(p);
    PointType pp2 = plane2.projection(p);
    double scal = (p.x() - pp1.x()) * (p.x() - pp2.x())
      + (p.y() - pp1.y()) * (p.y() - pp2.y())
      + (p.z() - pp1.z()) * (p.z() - pp2.z());
    return (scal <= 0.);
  }
public:
  Is_between(const PlaneType& p1, const PlaneType& p2, const bool xc)
    : plane1(p1), plane2(p2), xcondition(xc){}
  Is_between(const Is_between& ib)
    : plane1(ib.plane1), plane2(ib.plane2), xcondition(ib.xcondition){}
private:
  PlaneType plane1;
  PlaneType plane2;
  bool xcondition;
};


int main(int argc, char* argv[])
{
  std::ofstream fx("compare_iso_aniso.txt");
  CGAL::default_random = CGAL::Random(0);

  K::FT R = (argc > 1) ? atof(argv[1]) : 10.;
  K::FT r = 1.;
  K::FT r0 = (argc > 2) ? atof(argv[2]) : 1.0/*default*/;
  K::FT gamma0 = (argc > 3) ? atof(argv[3]) : 1.5;
  int nb = (argc > 4) ? atoi(argv[4]) : 10;
  K::FT approx = (argc > 5) ? atof(argv[5]) : 0.02;
  K::FT lambda = (argc > 6) ? atof(argv[6]) : (r / R);

  fx << "Torus :" << std::endl;
  fx << "\tR = " << R << std::endl;
  fx << "\tr = " << r << std::endl;
  fx << "\teps = " << lambda << std::endl;
  fx << "\tr0 = " << r0 << std::endl;
  fx << "\tgamma0 = " << gamma0 << std::endl;

  Constrain_surface_3_torus<K>* pdomain
    = new Constrain_surface_3_torus<K>(R, r);

  K::FT y = (argc > 7) ? atof(argv[7]) : (R + r + 1);
  int xcondition = (argc > 8) ? atoi(argv[8]) : -1;//default : no condition on x
  K::Plane_3 plane1(0., 1., 0., -y);
  K::Plane_3 plane2(0., 1., 0., y);

  typedef Is_between<K::Plane_3, K::Point_3> RCondition;
  RCondition condition(plane1, plane2, (xcondition == 1));
  
  Criteria_base<K> criteria_curvature(3.0, //radius_edge_ratio_
    0.2, //sliverity_
    r0, //circumradius_ 0.1
    gamma0, //distortion_ 1.3
    2.5, //beta_
    0.3, //delta_
    60, //max tries in picking region
    approx);//
  fx << std::endl << "Criteria for curvature metric field : " << std::endl;
  criteria_curvature.report(fx);
  fx << std::endl;

  Implicit_curvature_metric_field<K> curvature_mf(*pdomain, lambda);
  Surface_star_set_3<K, RCondition> curvature_mesh(criteria_curvature, curvature_mf, pdomain, nb, condition);
  curvature_mesh.refine_all();
  fx << "Aniso :\n";
  fx << "  nb stars = " << curvature_mesh.number_of_surface_stars() << std::endl;
  //approx = curvature_mesh.compute_approximation_error();
  fx << "  approx   = " << approx << std::endl;
  curvature_mesh.output("mesh_curvature.off");
  curvature_mesh.clear();


  Criteria_base<K> criteria_euclidean(3.0, //radius_edge_ratio_
    0.2, //sliverity_
    r0, //circumradius_ 0.1
    gamma0, //distortion_ 1.3
    2.5, //beta_
    0.3, //delta_
    60, //max tries in picking region
    approx);//
  fx << std::endl << "Criteria for euclidean metric field : " << std::endl;
  criteria_euclidean.report(fx);
  fx << std::endl;
  
  Euclidean_metric_field<K> euclid_mf;
  Surface_star_set_3<K, RCondition> euclidean_mesh(criteria_euclidean, euclid_mf, pdomain, nb, condition);
  euclidean_mesh.refine_all(); 
  fx << "Iso :\n";
  fx << "  nb stars = " << euclidean_mesh.number_of_surface_stars() << std::endl;
  fx << "  approx   = " << euclidean_mesh.compute_approximation_error() << std::endl;
  euclidean_mesh.output("mesh_euclidean.off");
  euclidean_mesh.clear();

  delete pdomain;
  return 0;
}

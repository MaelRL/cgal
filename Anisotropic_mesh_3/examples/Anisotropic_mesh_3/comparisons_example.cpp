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
  CGAL::default_random = CGAL::Random(0);
  int n = 1;

//geometry
  K::FT R = (argc > n) ? atof(argv[n++]) : 10.;
  K::FT r = (argc > n) ? atof(argv[n++]) : 1.;

//metric field
  K::FT epsilon = (argc > n) ? atof(argv[n++]) : 1e-6;

//facet criteria
  K::FT approx = (argc > n) ? atof(argv[n++]) : 0.1;
  K::FT gamma0 = (argc > n) ? atof(argv[n++]) : 1.5;
  K::FT f_rho0 = (argc > n) ? atof(argv[n++]) : 3.0;
  K::FT f_r0 = (argc > n) ? atof(argv[n++]) : 0.1;

//cell criteria
  K::FT sliverity = (argc > n) ? atof(argv[n++]) : 0.;
  K::FT c_rho0 = (argc > n) ? atof(argv[n++]) : 0.;
  K::FT c_r0 = (argc > n) ? atof(argv[n++]) : 0.;
  bool c_consistency = false;

//misc
  K::FT beta = (argc > n) ? atof(argv[n++]) : 2.5;
  K::FT delta = (argc > n) ? atof(argv[n++]) : 0.3;
  int max_times_to_try = (argc > n) ? atoi(argv[n++]) : 60;
  int nb = (argc > n) ? atoi(argv[n++]) : 20;

  Criteria_base<K>* criteria = new Criteria_base<K>(approx, gamma0, f_rho0, f_r0,
                                                    sliverity, c_rho0, c_r0,
                                                    c_consistency, beta, delta,
                                                    max_times_to_try);

  K::FT y = (argc > n) ? atof(argv[n++]) : (R + r + 1);
  int xcondition = (argc > n) ? atoi(argv[n++]) : -1;//default : no condition on x
  K::Plane_3 plane1(0., 1., 0., -y);
  K::Plane_3 plane2(0., 1., 0., y);

  typedef Is_between<K::Plane_3, K::Point_3> RCondition;
  RCondition condition(plane1, plane2, (xcondition == 1));

  Constrain_surface_3_torus<K>* pdomain =
      new Constrain_surface_3_torus<K>(R, r);
  Implicit_curvature_metric_field<K>* curvature_mf =
      new Implicit_curvature_metric_field<K>(*pdomain, epsilon);
  Euclidean_metric_field<K>* euclid_mf =
      new Euclidean_metric_field<K>();

  Surface_star_set_3<K, RCondition> curvature_mesh(criteria, curvature_mf, pdomain, nb, 0 /*pass_n*/, condition);
  curvature_mesh.refine_all();
  curvature_mesh.output("mesh_curvature.off");

  Surface_star_set_3<K, RCondition> euclidean_mesh(criteria, euclid_mf, pdomain, nb, 0 /*pass_n*/, condition);
  euclidean_mesh.refine_all(); 
  euclidean_mesh.output("mesh_euclidean.off");

  delete pdomain;
  return 0;
}

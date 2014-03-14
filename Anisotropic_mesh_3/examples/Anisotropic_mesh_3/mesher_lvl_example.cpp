#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Euclidean_metric_field.h>
#include <Domain/Constrain_surface_3_ellipse.h>

#include <CGAL/IO/Star_set_output.h>

#include <CGAL/Stretched_Delaunay_3.h>
#include <CGAL/Anisotropic_mesher_3.h>

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef typename K::FT                                       FT;
typedef CGAL::Timer                                          Timer;
typedef Stretched_Delaunay_3<K>                              Star;
typedef std::vector<Star*>                                   Star_vector;

std::string output_filename(const double& a, const double& b, const double &c)
{
  std::ostringstream oss;
  oss << "ellipsoid_" << a << "_" << b << "_" << c << ".off";
  return oss.str();
}

int main(int argc, char** argv)
{
  Timer timer;
  CGAL::default_random = CGAL::Random(0);

  int n = 1;

//geometry
  FT a = (argc > n) ? atof(argv[n++]) : 10.;
  FT b = (argc > n) ? atof(argv[n++]) : 1.;
  FT c = (argc > n) ? atof(argv[n++]) : 1.;

//metric field
  FT epsilon = (argc > n) ? atof(argv[n++]) : 1e-6;

//facet criteria
  FT approx = (argc > n) ? atof(argv[n++]) : 1.0;
  FT gamma0 = (argc > n) ? atof(argv[n++]) : 1.5;
  FT f_rho0 = (argc > n) ? atof(argv[n++]) : 3.0;
  FT f_r0 = (argc > n) ? atof(argv[n++]) : 1.0;

//cell criteria
  FT sliverity = (argc > n) ? atof(argv[n++]) : 0.2;
  FT c_rho0 = (argc > n) ? atof(argv[n++]) : 3.0;
  FT c_r0 = (argc > n) ? atof(argv[n++]) : 1.0;
  bool c_consistency = true;

//misc
  FT beta = (argc > n) ? atof(argv[n++]) : 2.5;
  FT delta = (argc > n) ? atof(argv[n++]) : 0.3;
  int max_times_to_try = (argc > n) ? atoi(argv[n++]) : 60;
  int nb = (argc > n) ? atoi(argv[n++]) : 20;

  Criteria_base<K>* criteria = new Criteria_base<K>(approx, gamma0, f_rho0, f_r0,
                                                    sliverity, c_rho0, c_r0,
                                                    c_consistency, beta, delta,
                                                    max_times_to_try);

  timer.start();

  Constrain_surface_3_ellipse<K>* pdomain =
      new Constrain_surface_3_ellipse<K>(a, b, c);
  Euclidean_metric_field<K>* metric_field = new Euclidean_metric_field<K>();
  Star_vector stars;

  Anisotropic_mesher_3<K> mesher(stars, pdomain, criteria, metric_field);
  double elapsed_time = mesher.refine_mesh();
  std::cout << "elapsed time: " << elapsed_time << std::endl;
  mesher.report();

  std::ofstream out("bambimboum.mesh");
  output_medit(stars, out);
  std::ofstream out_facet("bambimboum_surf.mesh");
  output_surface_medit(stars, out_facet);
  std::ofstream out_off("bambimboum.off");
  output_off(stars, out_off);
  std::ofstream out_surface_off("bambimboum_surf.off");
  output_surface_off(stars, out_surface_off);

  delete criteria;
  delete pdomain;
  delete metric_field;

  return 0;
}


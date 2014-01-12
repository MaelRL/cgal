//#define ANISO_USE_INSIDE_EXACT
//#define ANISO_DEBUG

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

//geometry
  FT a = 10.;
  FT b = 1.;
  FT c = 1.;

//metric field
  FT epsilon = 0.;

//facet criteria
  FT approx = 0.05;
  FT gamma0 = 2.;
  FT f_rho0 = 3.0;
  FT f_r0 = 1.0;

//cell criteria
  FT sliverity = 0.2;
  FT c_rho0 = 3.0;
  FT c_r0 = 1.0;
  bool c_consistency = true;

//misc
  int nb = 20;
  FT beta = 2.5;
  FT delta = 0.3;
  int max_times_to_try = 60;

  Criteria_base<K>* criteria = new Criteria_base<K>(approx, gamma0, f_rho0, f_r0, sliverity,
                                                    c_rho0, c_r0, c_consistency, beta, delta,
                                                    max_times_to_try);

  timer.start();

  Constrain_surface_3_ellipse<K>* pdomain = new Constrain_surface_3_ellipse<K>(a, b, c);
  Euclidean_metric_field<K>* metric_field = new Euclidean_metric_field<K>();
  Star_vector stars;

  Anisotropic_mesher_3<K> mesher(stars, pdomain, criteria, metric_field);
  double elapsed_time = mesher.refine_mesh();
  std::cout << elapsed_time << std::endl;

  std::ofstream out("mshlvl.mesh");
  output_medit(stars, out);
  std::ofstream out_facet("mshlvl_surf.mesh");
  output_surface_medit(stars, out_facet);
  std::ofstream out_off("mshlvl.off");
  output_off(stars, out_off);
  std::ofstream out_surface_off("mshlvl_surf.off");
  output_surface_off(stars, out_surface_off);

  delete criteria;
  delete pdomain;
  delete metric_field;

  return 0;
}


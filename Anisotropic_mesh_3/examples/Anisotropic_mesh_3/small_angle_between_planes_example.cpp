#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <Metric_field/Custom_metric_field.h>
#include <Metric_field/Hyperbolic_shock_metric_field.h>
#include <Metric_field/Euclidean_metric_field.h>
#include <CGAL/Constrain_surface_3_polyhedral.h>

#include <CGAL/IO/Star_set_IO.h>

#include <CGAL/Stretched_Delaunay_3.h>
#include <CGAL/Anisotropic_mesher_3.h>

#include <fstream>

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef typename K::FT                                       FT;
typedef Stretched_Delaunay_3<K>                              Star;
typedef std::vector<Star*>                                   Star_vector;

void build_geometry(const double angle, const double base_length)
{
  double height = std::tan(angle)*base_length;

  std::ofstream out("small_angle_between_planes_input.off");

  /*       _ F
         _//|
       _/ / | h
     _/ _/  |_________________ A
 E  /  /   _/\__ B          _/
   |  /  _/     \__       _/
  h| / _/          \__  _/  l
   |/_/_______________\/ D
  C         l

   */

  out << "OFF" << std::endl;
  out << "6 8 0" << std::endl << std::endl;

  //vertices
  out << "0 0 0" << std::endl;
  out << base_length << " 0 0" << std::endl;
  out << base_length << " " << base_length << " 0" << std::endl;
  out << "0 " << base_length << " 0" << std::endl;
  out << base_length << " " << base_length << " " << height << std::endl;
  out << base_length << " 0 " << height << std::endl;

  //facets
  out << "3 1 0 3" << std::endl; // ABD
  out << "3 0 5 4" << std::endl; // AFE
  out << "3 0 4 3" << std::endl; // AED
  out << "3 1 5 0" << std::endl; // AFB
  out << "3 3 2 1" << std::endl; // BCD
  out << "3 2 5 1" << std::endl; // BFC
  out << "3 4 5 2" << std::endl; // CFE
  out << "3 2 3 4" << std::endl; // CDE
}

int main(int argc, char** argv)
{
  CGAL::default_random = CGAL::Random(0);

  int n = 2;

//geometry
  FT angle = (argc > n) ? atof(argv[n++]) : 0.05;
  FT bl = 1./std::tan(angle); //length of the side of the "base" ABCD

//metric field
  FT epsilon = (argc > n) ? atof(argv[n++]) : 1e-6;
  FT en_factor = (argc > n) ? atof(argv[n++]) : 0.999;

//facet criteria
  FT approx = (argc > n) ? atof(argv[n++]) : 0.1;
  FT gamma0 = (argc > n) ? atof(argv[n++]) : 1.5;
  FT f_rho0 = (argc > n) ? atof(argv[n++]) : 3.0;
  FT f_r0 = (argc > n) ? atof(argv[n++]) : 1.0;

//cell criteria
  FT sliverity = (argc > n) ? atof(argv[n++]) : 0.2;
  FT c_rho0 = (argc > n) ? atof(argv[n++]) : 3.0;
  FT c_r0 = (argc > n) ? atof(argv[n++]) : 1.0;

//misc
  FT beta = (argc > n) ? atof(argv[n++]) : 2.5;
  FT delta = (argc > n) ? atof(argv[n++]) : 0.3;
  int max_times_to_try = (argc > n) ? atoi(argv[n++]) : 60;
  int nb = (argc > n) ? atoi(argv[n++]) : 20;

  build_geometry(angle, 2*bl);
  Criteria_base<K>* criteria = new Criteria_base<K>(approx, f_r0, f_rho0,
                                                    c_r0, c_rho0, sliverity,
                                                    gamma0, beta, delta, nb,
                                                    max_times_to_try);

  Constrain_surface_3_polyhedral<K>* pdomain = new Constrain_surface_3_polyhedral<K>(argv[1]);
  Custom_metric_field<K>* metric_field = new Custom_metric_field<K>();
//  Euclidean_metric_field<K>* metric_field = new Euclidean_metric_field<K>();
//  Hyperbolic_shock_metric_field<K>* metric_field = new Hyperbolic_shock_metric_field<K>(0.6, epsilon, en_factor);

  Starset<K> starset;

  Anisotropic_mesher_3<K> mesher(starset, pdomain, criteria, metric_field);
  double elapsed_time = mesher.refine_mesh();
  std::cout << elapsed_time << std::endl;

  std::ofstream out("wedge.mesh");
  output_medit(starset, out);
  std::ofstream out_facet("wedge_surf.mesh");
  output_surface_medit(starset, out_facet);

  delete criteria;
  delete pdomain;
  delete metric_field;

  return 0;
}


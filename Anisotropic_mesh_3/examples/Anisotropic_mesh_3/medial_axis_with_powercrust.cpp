
#include <CGAL/Default_configuration.h>
#define ANISO_OUTPUT_POINTS_FOR_MEDIAL_AXIS

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Robust_circumcenter_traits_3.h>

#include <CGAL/Constrain_surface_3_polyhedral.h>
#include <CGAL/Polyhedral_curvature_metric_field.h>
#include <CGAL/Anisotropic_surface_mesher_3.h>

#include <CGAL/Euclidean_metric_field.h>
#include <Metric_field/Arctan_metric_field.h>


#define CGAL_PROFILE

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel	Epick;
typedef CGAL::Robust_circumcenter_traits_3<Epick> K;

typedef Anisotropic_surface_mesher_3<K> Mesher;


std::string output_filename(const std::string& filename)
{
  std::size_t found = filename.find_last_of("/\\");
  std::string file = filename.substr(found+1, filename.length()-5);
  file.append("_output.off");
  return file;
}

int main(int argc, char *argv[])
{
  CGAL::default_random = CGAL::Random(0);
  std::string file(argv[1]);
  std::cout << "Let's mesh : " << file << std::endl;

  double r0 = (argc > 2) ? atof(argv[2]) : 1.;
  double gamma0 = (argc > 3) ? atof(argv[3]) : 1.8;
  double approx = (argc > 4) ? atof(argv[4]) : 0.;
  double epsilon = (argc > 5) ? atof(argv[5]) : (0.1);
  
  Criteria_base<K>* criteria = new Criteria_base<K>(0., //radius_edge_ratio_ //3.0
                                                    0., //sliverity_ (not used)
                                                    r0, //circumradius_
                                                    gamma0, //distortion_
                                                    2.5, //beta_
                                                    0.3, //delta_
                                                    60, //nb tries in pick_valid
                                                    approx);

  Constrain_surface_3_polyhedral<K>* pdomain
    = new Constrain_surface_3_polyhedral<K>(argv[1], epsilon);

  Polyhedral_curvature_metric_field<K>* metric_field =
    new Polyhedral_curvature_metric_field<K>(*pdomain, epsilon);
  
  Surface_star_set_3<K> starset(criteria, metric_field, pdomain);

//  starset.refine_all();

//  std::string outfile = output_filename(file);
//  starset.output(outfile);
  
  delete pdomain;
  return 0;
}


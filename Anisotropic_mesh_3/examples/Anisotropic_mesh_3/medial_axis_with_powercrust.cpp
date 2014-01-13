
#include <CGAL/Default_configuration.h>
#define ANISO_GIVE_POINTS_FOR_MEDIAL_AXIS

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Robust_circumcenter_traits_3.h>

#include <CGAL/Constrain_surface_3_polyhedral.h>
#include <CGAL/Polyhedral_curvature_metric_field.h>
#include <CGAL/Anisotropic_surface_mesher_3.h>

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel	Epick;
typedef CGAL::Robust_circumcenter_traits_3<Epick> K;

typedef Anisotropic_surface_mesher_3<K> Mesher;


std::string output_filename(const std::string& filename)
{
  std::size_t found = filename.find_last_of("/\\");
  std::string file = filename.substr(found+1, filename.length()-4);
  file.append("_output.off");
  return file;
}


int main(int argc, char *argv[])
{
  CGAL::default_random = CGAL::Random(0);
  std::string file(argv[1]);
  std::cout << "Let's mesh : " << file << std::endl;

  char* poles_filename;
  int n = 2;
  bool poles_given = false;
  if(argc > 3)
  {
    std::string poles(argv[2]);
    if(poles.find("-poles") != std::string::npos)
    {
      poles_filename = argv[3];
      std::cout << "poles filename : " << poles_filename << std::endl;
      n = n + 2;
      poles_given = true;
    }
    else std::cerr << "argc is " << argc << std::endl;
  }

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

  Constrain_surface_3_polyhedral<K>* pdomain
    = new Constrain_surface_3_polyhedral<K>(file.data(), epsilon);

  Polyhedral_curvature_metric_field<K>* metric_field =
    new Polyhedral_curvature_metric_field<K>(*pdomain, epsilon);
  
  Surface_star_set_3<K> starset(criteria, metric_field, pdomain, 10, 1, No_condition<K::Point_3>(), poles_given, poles_filename);

  starset.refine_all();

  std::string outfile = output_filename(file);
  starset.output(outfile.c_str());
  
  delete pdomain;
  return 0;
}


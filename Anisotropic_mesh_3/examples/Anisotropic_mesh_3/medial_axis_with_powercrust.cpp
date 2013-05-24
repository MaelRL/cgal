
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

void read_poles(const char* filename,/*.off file*/
                std::set<K::Point_3>& poles)
{
  std::string line;
  std::ifstream polesfile(filename);
  if(polesfile.is_open())
  {
    bool nbv_found = false;
    std::size_t nbv = 0;
    while(polesfile.good())
    {
      std::getline(polesfile,line);
      if(line.find("OFF") != std::string::npos)
        continue;
      else if(line.find("#") != std::string::npos)//comments
        continue;
      else if(!nbv_found)
      {
        polesfile >> nbv;
        std::cout << "Nb poles is : " << nbv << std::endl;
        nbv_found = true;
      }
      else if(poles.size() < nbv)
      {
        K::FT x,y,z;
        polesfile >> x >> y >> z;
        poles.insert(K::Point_3(x,y,z));
      }
      else break;
    }
    polesfile.close();
  }
  else 
    std::cout << "Unable to open the file containing poles : " << filename << "\n"; 
}

int main(int argc, char *argv[])
{
  CGAL::default_random = CGAL::Random(0);
  std::string file(argv[1]);
  std::cout << "Let's mesh : " << file << std::endl;

  char* poles_filename;
  if(argc > 3 && argv[2] == "-poles")
  {
    poles_filename = argv[3];
    std::cout << "poles filename : " << poles_filename << std::endl;
  }
  int i=4;
  double r0 = (argc > i) ? atof(argv[i++]) : 1.;
  double gamma0 = (argc > i) ? atof(argv[i++]) : 1.8;
  double approx = (argc > i) ? atof(argv[i++]) : 0.;
  double epsilon = (argc > i) ? atof(argv[i++]) : (0.1);
  
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
  
  std::set<K::Point_3> poles;
  read_poles(poles_filename, poles);

  Surface_star_set_3<K> starset(criteria, metric_field, pdomain,
    10, 2., No_condition<K::Point_3>(), poles);

//  starset.refine_all();

//  std::string outfile = output_filename(file);
//  starset.output(outfile);
  
  delete pdomain;
  return 0;
}


// #define NO_USE_AABB_TREE_OF_BBOXES
// #define LIFT_AABB_TREE
// #define ANISO_OUTPUT_WIP
// #define REJECT_FAILED_ELEMENTS
// #define ANISO_USE_CUSTOM_CONSISTENCY
// #define ANISO_USE_ENCROACH_SKIP
// #define ANISO_NO_CONSISTENCY

// should use at most one (if using QUEUE, using SKIP is pointless)
// #define ANISO_USE_DISTORTION_QUEUE
#define SKIP_PICK_VALID_DISTORTION

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Starset.h>
#include <CGAL/Stretched_Delaunay_2.h>
#include <CGAL/Anisotropic_mesher_2.h>

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <Domain/Rectangle_domain.h>

#include <CGAL/IO/Star_set_output.h>

#include <boost/timer/timer.hpp>
#include <cmath>

using namespace CGAL::Anisotropic_mesh_2;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef typename K::FT                                       FT;
typedef CGAL::Timer                                          Timer;
typedef Stretched_Delaunay_2<K>                              Star;
typedef typename Star::Point_2                               Point_2;

int main(int argc, char** argv)
{
//  std::freopen("log.txt", "w", stdout);
  boost::timer::auto_cpu_timer t;

  std::streambuf * old = std::cout.rdbuf();
//  std::cout.rdbuf(0);

  Timer timer;
  CGAL::default_random = CGAL::Random(0);

  int n = 1;

//geometry
  FT a = (argc > n) ? atof(argv[n++]) : 2.0; // half the side
  FT b = (argc > n) ? atof(argv[n++]) : 2.0;
  Point_2 center(0., 0.);

//metric field
  FT epsilon = (argc > n) ? atof(argv[n++]) : 1e-6;

//face criteria
  FT f_r0 = (argc > n) ? atof(argv[n++]) : 1.0;
  FT f_rho0 = (argc > n) ? atof(argv[n++]) : 3.0;

//misc
  FT gamma = (argc > n) ? atof(argv[n++]) : 2;
  FT beta = (argc > n) ? atof(argv[n++]) : 3;
  FT delta = (argc > n) ? atof(argv[n++]) : 0.3;
  int max_times_to_try = (argc > n) ? atoi(argv[n++]) : 60;
  int nb = (argc > n) ? atoi(argv[n++]) : 20;

  Criteria_base<K>* criteria = new Criteria_base<K>(f_r0, f_rho0,
                                                    gamma, beta, delta, nb,
                                                    max_times_to_try);

  timer.start();

  //----------- pick a domain! ----------
  Rectangle_domain<K>* pdomain = new Rectangle_domain<K>(a, b, center);

  //----------- pick a metric field! ----
//  Euclidean_metric_field<K>* metric_field = new Euclidean_metric_field<K>(1000., 10.);
  Custom_metric_field<K>* metric_field = new Custom_metric_field<K>(epsilon);

  Starset<K> starset;

  //----------- pick a mesher! -----------
  Anisotropic_mesher_2<K> mesher(starset, pdomain, criteria, metric_field);

#if 1
  double elapsed_time = mesher.refine_mesh();

  std::cout.rdbuf(old);
  std::cout << "elapsed time: " << elapsed_time << std::endl;
  mesher.report();

  std::ofstream out("bambimboum.mesh");
  output_medit(starset, out, false/*consistent_only*/);

//  starset.draw_metric_vector_field();
//  starset.draw_metric_vector_field_2(metric_field);
//  starset.output_radius_of_stars();

#else
// IF YOU WANT TO CONSTRAIN STARS (sweep)
//  starset.constrain();
//  face_edge_length_histogram(starset);

// IF YOU WANT TO RESUME FROM AN EXISTING MESH
  std::cout << "resuming!" << std::endl;
  starset.clear();
  mesher.resume_from_mesh_file("resume_in.mesh");
#endif

  delete metric_field;
  delete pdomain;
  delete criteria;

  return 0;
}


// Copyright (c) 2013 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

//#define ANISO_USE_INSIDE_EXACT

#include <CGAL/Timer.h>
#include <CGAL/Default_configuration.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Implicit_curvature_metric_field.h>

#include <Domain/Constrain_surface_3_double_ellipsoid.h>

#include <CGAL/Anisotropic_surface_mesher_3.h>

#include "refinement_condition_is_between_planes.h"

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel	K;
typedef CGAL::Timer Timer;

std::string output_filename(const double& a, const double& b, const double &c)
{
  std::ostringstream oss;
  oss << "ellipsoid_" << a << "_" << b << "_" << c << ".off";
  return oss.str();
}

int main(int argc, char* argv[])
{
  std::ofstream fx("ellipsoid_timings.txt");

  Timer timer;
  CGAL::default_random = CGAL::Random(0);

  K::FT a = (argc > 1) ? atof(argv[1]) : 0.15;
  K::FT b = (argc > 2) ? atof(argv[2]) : 0.25;
  K::FT c = (argc > 3) ? atof(argv[3]) : 0.35;
  K::FT d = (argc > 4) ? atof(argv[4]) : 0.1;
  K::FT e = (argc > 5) ? atof(argv[5]) : 0.2;
  K::FT f = (argc > 6) ? atof(argv[6]) : 0.3;
  K::Point_3 s_center(0.25, 0.1, 0.1);

  K::FT epsilon = (argc > 7) ? atof(argv[7]) : 1.0;

  K::FT r0 = (argc > 8) ? atof(argv[8]) : 1.0/*default*/;
  K::FT gamma0 = (argc > 9) ? atof(argv[9]) : 1.5;

  int nb = (argc > 10) ? atoi(argv[10]) : 10;

  K::FT beta = (argc > 11) ? atof(argv[11]) : 2.5;
  K::FT delta = (argc > 12) ? atof(argv[12]) : 0.3;

  fx << "Ellipsoid :" << std::endl;
  fx << "\ta = " << a << std::endl;
  fx << "\tb = " << b << std::endl;
  fx << "\tc = " << c << std::endl;
  fx << "\teps = " << epsilon << std::endl;
  fx << "\tr0 = " << r0 << std::endl;
  fx << "\tgamma0 = " << gamma0 << std::endl;
  fx << "\tbeta = " << beta << std::endl;
  fx << "\tdelta = " << delta << std::endl;
  
  Criteria_base<K> criteria(3.0, //radius_edge_ratio_
    0.2,    //sliverity_
    r0,     //circumradius_ 0.1
    gamma0, //distortion_ 1.3
    beta,   //beta_ 2.5
    delta); //delta_ 0.3

  fx << std::endl << "nbV" << "\t" << "time" << std::endl;
  timer.start();

  Constrain_surface_3_double_ellipsoid<K>* pdomain =
          new Constrain_surface_3_double_ellipsoid<K>(a, b, c, d, e, f, s_center);

  Implicit_curvature_metric_field<K> metric_field(*pdomain, epsilon);
  //Ellipsoid_metric_field<K> metric_field(a, b, c, epsilon);

  int xcondition = (argc > 13) ? atoi(argv[13]) : -1;//default : no condition on x
  K::Plane_3 plane1(1., 0., 0., -4);
  K::Plane_3 plane2(1., 0., 0., 4);

  typedef Is_between<K::Plane_3, K::Point_3> RCondition;

  RCondition condition(plane1, plane2, (xcondition == 1));
  //Surface_star_set_3<K, RCondition> starset(criteria, metric_field, pdomain, nb, condition);
  Surface_star_set_3<K> starset(criteria, metric_field, pdomain, nb);

  timer.stop();
  starset.output("base_double_ellipse.off", false);

  std::cout << "n : " << starset.number_of_stars() << std::endl;
  std::cout << "surf : " << starset.number_of_surface_stars() << std::endl;
  std::cout << "restricted facets : " << starset.count_restricted_facets() << std::endl;

  fx << starset.number_of_stars() << "\t" <<  timer.time() << std::endl;
  starset.refine_all();

  starset.output("double_ellipse.off");

  delete pdomain;
  return 0;

}


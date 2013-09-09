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
//#define ANISO_DEBUG
#undef ANISO_VERBOSE

#include <CGAL/Default_configuration.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Implicit_curvature_metric_field.h>
#include <Domain/Constrain_surface_3_ellipse.h>

#include <CGAL/Anisotropic_surface_mesher_3.h>

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel	K;

int main(int argc, char* argv[])
{
  std::ofstream fx("ellipsoid_parameters.txt");

  CGAL::default_random = CGAL::Random(0);

  K::FT a = (argc > 1) ? atof(argv[1]) : 10.;
  K::FT b = (argc > 2) ? atof(argv[2]) : 1.;
  K::FT c = (argc > 3) ? atof(argv[3]) : 1.;

  K::FT epsilon = (argc > 4) ? atof(argv[4]) : 0.;

  K::FT r0 = (argc > 5) ? atof(argv[5]) : 0./*default*/;
  K::FT gamma0 = (argc > 6) ? atof(argv[6]) : 0.;
  K::FT approx = (argc > 7) ? atof(argv[7]) : 0.;

  int nb = (argc > 8) ? atoi(argv[8]) : 10;

  K::FT beta = (argc > 9) ? atof(argv[9]) : 2.5;
  K::FT delta = (argc > 10) ? atof(argv[10]) : 0.3;

  K::FT sigma0 = 0.;

  fx << "Ellipsoid :" << std::endl;
  fx << "\ta = " << a << std::endl;
  fx << "\tb = " << b << std::endl;
  fx << "\tc = " << c << std::endl;
  fx << "\teps = " << epsilon << std::endl;
  fx << "\tr0 = " << r0 << std::endl;
  fx << "\tgamma0 = " << gamma0 << std::endl;
  fx << "\tsigma0 = " << sigma0 << std::endl;
  fx << "\tapprox = " << approx << std::endl;
  fx << "\tbeta = " << beta << std::endl;
  fx << "\tdelta = " << delta << std::endl;
  
  // RUN ANISO MESHER
  Criteria_base<K>* criteria = new Criteria_base<K>(0., //radius_edge_ratio_
                                                    sigma0,    //sliverity_
                                                    r0,     //circumradius_ 0.1
                                                    gamma0, //distortion_ 1.3
                                                    beta,   //beta_ 2.5
                                                    delta, //delta_ 0.3
                                                    60,
                                                    approx);//approx

  Constrain_surface_3_ellipse<K>* pdomain = new Constrain_surface_3_ellipse<K>(a, b, c);

  Implicit_curvature_metric_field<K>* metric_field =
    new Implicit_curvature_metric_field<K>(*pdomain, epsilon);

  Surface_star_set_3<K> starset(criteria, metric_field, pdomain, nb);
  starset.refine_all();
  std::cout << "Star set has " << starset.number_of_surface_stars() << " surface vertices.\n";

  // RUN ISO MESHER
  approx = starset.compute_approximation_error();
  pdomain->run_mesh_3(approx, true/*verbose*/);

  delete pdomain;
  return 0;
}


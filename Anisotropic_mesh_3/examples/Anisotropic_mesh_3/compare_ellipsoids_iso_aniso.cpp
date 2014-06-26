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

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

int main(int argc, char* argv[])
{
  CGAL::default_random = CGAL::Random(0);
  int n = 1;

//geometry
  K::FT a = (argc > n) ? atof(argv[n++]) : 10.;
  K::FT b = (argc > n) ? atof(argv[n++]) : 1.;
  K::FT c = (argc > n) ? atof(argv[n++]) : 1.;

//metric field
  K::FT epsilon = (argc > n) ? atof(argv[n++]) : 1e-6;

//facet criteria
  K::FT approx = (argc > n) ? atof(argv[n++]) : 0.1;
  K::FT f_r0 = (argc > n) ? atof(argv[n++]) : 0.1;
  K::FT f_rho0 = (argc > n) ? atof(argv[n++]) : 3.0;

//cell criteria
  K::FT c_r0 = (argc > n) ? atof(argv[n++]) : 0.;
  K::FT c_rho0 = (argc > n) ? atof(argv[n++]) : 0.;
  K::FT sliverity = (argc > n) ? atof(argv[n++]) : 0.;

//misc
  K::FT gamma0 = (argc > n) ? atof(argv[n++]) : 1.5;
  K::FT beta = (argc > n) ? atof(argv[n++]) : 2.5;
  K::FT delta = (argc > n) ? atof(argv[n++]) : 0.3;
  int max_times_to_try = (argc > n) ? atoi(argv[n++]) : 60;
  int nb = (argc > n) ? atoi(argv[n++]) : 20;

  Criteria_base<K>* criteria = new Criteria_base<K>(approx, f_r0, f_rho0,
                                                    c_r0, c_rho0, sliverity,
                                                    gamma0, beta, delta, nb,
                                                    max_times_to_try);

  Constrain_surface_3_ellipse<K>* pdomain = new Constrain_surface_3_ellipse<K>(a, b, c);

  Implicit_curvature_metric_field<K>* metric_field =
    new Implicit_curvature_metric_field<K>(*pdomain, epsilon);

  //aniso
  Starset<K> starset;
  Anisotropic_surface_mesher_3<K> mesher(starset, pdomain, criteria, metric_field);
  mesher.refine_mesh();
  std::cout << "Star set has " << starset.number_of_surface_stars() << " surface vertices.\n";

  //iso
//  approx = starset.compute_approximation_error(); fixme
  pdomain->run_mesh_3(approx, true/*verbose*/);

  delete metric_field;
  delete pdomain;
  delete criteria;

  return 0;
}


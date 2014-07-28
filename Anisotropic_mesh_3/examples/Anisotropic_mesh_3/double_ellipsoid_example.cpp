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

#include <CGAL/Default_configuration.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Anisotropic_surface_mesher_3.h>
#include <CGAL/Starset.h>

#include <CGAL/Implicit_curvature_metric_field.h>
#include <Domain/Constrain_surface_3_double_ellipsoid.h>

#include "refinement_condition_is_between_planes.h"

#include <string>

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel	K;

std::string output_filename(const double& a, const double& b, const double &c)
{
  std::ostringstream oss;
  oss << "ellipsoid_" << a << "_" << b << "_" << c << ".off";
  return oss.str();
}

int main(int argc, char* argv[])
{
  CGAL::default_random = CGAL::Random(0);
  int n = 1;

//geometry
  K::FT a = (argc > n) ? atof(argv[n++]) : 0.15;
  K::FT b = (argc > n) ? atof(argv[n++]) : 0.25;
  K::FT c = (argc > n) ? atof(argv[n++]) : 0.35;
  K::FT d = (argc > n) ? atof(argv[n++]) : 0.1;
  K::FT e = (argc > n) ? atof(argv[n++]) : 0.2;
  K::FT f = (argc > n) ? atof(argv[n++]) : 0.3;
  K::Point_3 s_center(0.25, 0.1, 0.1);

//metric field
  K::FT epsilon = (argc > n) ? atof(argv[n++]) : 1e-6;
  K::FT en_factor = (argc > n) ? atof(argv[n++]) : 0.999;

//facet criteria
  K::FT approx = (argc > n) ? atof(argv[n++]) : 1.0;
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

  Constrain_surface_3_double_ellipsoid<K>* pdomain =
          new Constrain_surface_3_double_ellipsoid<K>(a, b, c, d, e, f, s_center);

  Implicit_curvature_metric_field<K>* metric_field =
    new Implicit_curvature_metric_field<K>(*pdomain, epsilon, en_factor);

//  int xcondition = (argc > n) ? atoi(argv[n++]) : -1; //default : no condition on x
//  K::Plane_3 plane1(1., 0., 0., -4);
//  K::Plane_3 plane2(1., 0., 0., 4);
//  typedef Is_between<K::Plane_3, K::Point_3> RCondition;
//  RCondition condition(plane1, plane2, (xcondition == 1));

  Starset<K> starset;
  Anisotropic_surface_mesher_3<K> mesher(starset, pdomain, criteria, metric_field);

  output_surface_off(starset, "base_double_ellipse.off", false);

  mesher.refine_mesh();
  output_surface_off(starset, "double_ellipse.off");

  delete metric_field;
  delete pdomain;
  delete criteria;

  return 0;
}


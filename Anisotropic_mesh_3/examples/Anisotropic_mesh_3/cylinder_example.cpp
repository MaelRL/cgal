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
#include <Metric_field/Cylinder_metric_field.h>
#include <Domain/Constrain_surface_3_cylinder.h>

#include "refinement_condition_is_between_planes.h"

#include <string>

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel	K;
typedef CGAL::Timer Timer;

std::string output_filename(const double& r, const double& h)
{
  std::ostringstream oss;
  oss << "cylinder_" << r << "_" << h << ".off";
  return oss.str();
}

int main(int argc, char* argv[])
{
  CGAL::default_random = CGAL::Random(0);
  int n = 1;

//geometry
  K::FT r = (argc > n) ? atof(argv[n++]) : 1.;
  K::FT h = (argc > n) ? atof(argv[n++]) : 10.;

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

  Constrain_surface_3_cylinder<K>* pdomain = new Constrain_surface_3_cylinder<K>(r, h);

  //Implicit_curvature_metric_field<K>* metric_field =
  //  new Implicit_curvature_metric_field<K>(*pdomain, epsilon);
  Cylinder_metric_field<K>* metric_field = new Cylinder_metric_field<K>(r, h, epsilon, en_factor);

//  int xcondition = (argc > n) ? atoi(argv[n++]) : -1;//default : no condition on x
//  K::Plane_3 plane1(0., 0., 1., -0.01*h);
//  K::Plane_3 plane2(0., 0., 1., -0.99*h);
//  typedef Is_between<K::Plane_3, K::Point_3> RCondition;
//  RCondition condition(plane1, plane2, (xcondition == 1));

  Starset<K> starset;
  Anisotropic_surface_mesher_3<K> mesher(starset, pdomain, criteria, metric_field);

  output_surface_off(starset, "base_closed_cylinder.off", false/*consistency*/);

  mesher.refine_mesh();
  output_surface_off(starset, "cylinder.off");

  delete metric_field;
  delete pdomain;
  delete criteria;

  return 0;
}


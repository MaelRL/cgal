// Copyright (c) 2011  INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
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
//
// Author(s) : Kan-Le Shi

#include <CGAL/Timer.h>
#include <CGAL/Default_configuration.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Implicit_curvature_metric_field.h>
#include <CGAL/Euclidean_metric_field.h>
#include <Metric_field/Torus_metric_field.h>

#include <Domain/Constrain_surface_3_torus.h>
#include <CGAL/Criteria.h>
#include <CGAL/Surface_star_set_3.h>

#include "refinement_condition_is_between_planes.h"

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Timer Timer;

std::string output_filename(const double& R,
                            const double& r)
{
  std::ostringstream oss;
  oss << "torus_" << R << "_" << r << ".off";
  return oss.str();
}

int main(int argc, char* argv[])
{
  CGAL::default_random = CGAL::Random(0);
  int n = 1;

//geometry
  K::FT R = (argc > n) ? atof(argv[n++]) : 10.;
  K::FT r = (argc > n) ? atof(argv[n++]) : 1.;

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

  Constrain_surface_3_torus<K>* pdomain
    = new Constrain_surface_3_torus<K>(R, r);
  Implicit_curvature_metric_field<K>* metric_field =
    new Implicit_curvature_metric_field<K>(*pdomain, epsilon);
  //Torus_metric_field<K>* metric_field = new Torus_metric_field<K>(R, r, epsilon);
  //Euclidean_metric_field<K>* metric_field = new Euclidean_metric_field<K>();

  std::string file = output_filename(R, r);

  //K::FT y = (argc > n) ? atof(argv[n++]) : (0.5*(R + r));
  //int xcondition = (argc > n) ? atoi(argv[n++]) : -1;//default : no condition on x
  K::Plane_3 plane1(0., 1., 0., 0.);
  double angle = M_PI/16.;
  K::Plane_3 plane2(K::Point_3(0.,0.,0.), K::Point_3(0.,0.,1.), K::Point_3(std::cos(angle), std::sin(angle),0.));
  typedef Is_between<K::Plane_3, K::Point_3> RCondition;
  RCondition condition(plane1, plane2, 1);//(xcondition == 1));

  //Surface_star_set_3<K, RCondition> starset(criteria, metric_field, pdomain, nb, 10, condition);
  Surface_star_set_3<K> starset(criteria, metric_field, pdomain, nb);

  starset.refine_all();
  starset.output_surface_star_set("torus_example.mesh");//mesh
  starset.output(file.c_str());

  delete pdomain;
  return 0;
}

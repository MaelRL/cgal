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

#include <CGAL/Timer.h>
#include <CGAL/Default_configuration.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Euclidean_metric_field.h>

#include <Domain/Constrain_surface_3_ellipse.h>

#include <CGAL/Anisotropic_tet_mesher_3.h>

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Timer Timer;

std::string output_filename(const double& a, const double& b, const double &c)
{
  std::ostringstream oss;
  oss << "ellipsoid_" << a << "_" << b << "_" << c << ".off";
  return oss.str();
}

int main(int argc, char* argv[])
{
  Timer timer;
  CGAL::default_random = CGAL::Random(0);
  int n = 1;

//geometry
  K::FT a = (argc > n) ? atof(argv[n++]) : 10.;
  K::FT b = (argc > n) ? atof(argv[n++]) : 1.;
  K::FT c = (argc > n) ? atof(argv[n++]) : 1.;

//metric field
  K::FT epsilon = (argc > n) ? atof(argv[n++]) : 1e-6;

//facet criteria
  K::FT approx = (argc > n) ? atof(argv[n++]) : 0.;
  K::FT gamma0 = (argc > n) ? atof(argv[n++]) : 1.5;
  K::FT f_rho0 = (argc > n) ? atof(argv[n++]) : 0.;
  K::FT f_r0 = (argc > n) ? atof(argv[n++]) : 0.;

//cell criteria
  K::FT sliverity = (argc > n) ? atof(argv[n++]) : 0.2;
  K::FT c_rho0 = (argc > n) ? atof(argv[n++]) : 3.0;
  K::FT c_r0 = (argc > n) ? atof(argv[n++]) : 0.1;
  bool c_consistency = true;

//misc
  K::FT beta = (argc > n) ? atof(argv[n++]) : 2.5;
  K::FT delta = (argc > n) ? atof(argv[n++]) : 0.3;
  int max_times_to_try = (argc > n) ? atoi(argv[n++]) : 60;
  int nb = (argc > n) ? atoi(argv[n++]) : 20;

  Criteria_base<K>* criteria = new Criteria_base<K>(approx, gamma0, f_rho0, f_r0,
                                                    sliverity, c_rho0, c_r0,
                                                    c_consistency, beta, delta,
                                                    max_times_to_try);

  timer.start();
  Constrain_surface_3_ellipse<K>* pdomain = new Constrain_surface_3_ellipse<K>(a, b, c);
  Euclidean_metric_field<K>* metric_field = new Euclidean_metric_field<K>();

  Cell_star_set_3<K> starset(criteria, metric_field, pdomain, nb);
  starset.dump();
  timer.stop();

  starset.refine_all();
  std::ofstream fmedit("cell_mesher_elli.mesh");
  starset.output_medit(fmedit);

  delete pdomain;
  return 0;
}


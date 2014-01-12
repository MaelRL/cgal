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
  std::ofstream fx("ellipsoid_timings.txt");

  Timer timer;
  CGAL::default_random = CGAL::Random(0);

  K::FT a = (argc > 1) ? atof(argv[1]) : 10.;
  K::FT b = (argc > 2) ? atof(argv[2]) : 1.;
  K::FT c = (argc > 3) ? atof(argv[3]) : 1.;

  K::FT epsilon = (argc > 4) ? atof(argv[4]) : 0.;

  K::FT r0 = (argc > 5) ? atof(argv[5]) : 1.0;
  K::FT gamma0 = (argc > 6) ? atof(argv[6]) : 5.;
  K::FT rho0 = (argc > 7) ? atof(argv[7]) : 3.0;
  K::FT approx = (argc > 8) ? atof(argv[8]) : 10.;

  int nb = (argc > 9) ? atoi(argv[9]) : 10;

  K::FT beta = (argc > 10) ? atof(argv[10]) : 2.5;
  K::FT delta = (argc > 11) ? atof(argv[11]) : 0.3;

  fx << "Ellipsoid :" << std::endl;
  fx << "\ta = " << a << std::endl;
  fx << "\tb = " << b << std::endl;
  fx << "\tc = " << c << std::endl;
  fx << "\teps = " << epsilon << std::endl;
  fx << "\tr0 = " << r0 << std::endl;
  fx << "\tgamma0 = " << gamma0 << std::endl;
  fx << "\tapprox = " << approx << std::endl;
  fx << "\tbeta = " << beta << std::endl;
  fx << "\tdelta = " << delta << std::endl;
  
  Criteria_base<K>* criteria = new Criteria_base<K>(rho0, //radius_edge_ratio_
                                                    0.2,    //sliverity_
                                                    r0,     //circumradius_ 0.1
                                                    gamma0, //distortion_ 1.3
                                                    beta,   //beta_ 2.5
                                                    delta, //delta_ 0.3
                                                    60,
                                                    approx);//approx

  fx << std::endl << "nbV" << "\t" << "time" << std::endl;
  timer.start();

  Constrain_surface_3_ellipse<K>* pdomain = new Constrain_surface_3_ellipse<K>(a, b, c);
  Euclidean_metric_field<K>* metric_field = new Euclidean_metric_field<K>();

  Cell_star_set_3<K> starset(criteria, metric_field, pdomain, nb);
  starset.dump();
  timer.stop();

  fx << starset.number_of_stars() << "\t" <<  timer.time() << std::endl;
  starset.refine_all();
  std::ofstream fmedit("cell_mesher_elli.mesh");
  starset.output_medit(fmedit);

  delete pdomain;
  return 0;
}


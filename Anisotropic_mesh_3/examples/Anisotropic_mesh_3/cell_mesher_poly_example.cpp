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

#include <CGAL/Default_configuration.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Constrain_surface_3_polyhedral.h>
#include <CGAL/Euclidean_metric_field.h>
#include <Metric_field/Hyperbolic_shock_metric_field.h>

#include <CGAL/Anisotropic_tet_mesher_3.h>

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel	K;

int main(int argc, char* argv[])
{
  CGAL::default_random = CGAL::Random(0);

  K::FT epsilon = (argc > 2) ? atof(argv[2]) : 0.;

  K::FT r0 = (argc > 3) ? atof(argv[3]) : 0.01;
  K::FT gamma0 = (argc > 4) ? atof(argv[4]) : 2.;
  K::FT rho0 = (argc > 5) ? atof(argv[5]) : 3.0;
  K::FT approx = (argc > 6) ? atof(argv[6]) : 10.;

  int nb = (argc > 7) ? atoi(argv[7]) : 20;

  K::FT beta = (argc > 8) ? atof(argv[8]) : 2.5;
  K::FT delta = (argc > 9) ? atof(argv[9]) : 0.3;
  
  Criteria_base<K>* criteria = new Criteria_base<K>(rho0, //radius_edge_ratio_
                                                    0.2,    //sliverity_
                                                    r0,     //circumradius_ 0.1
                                                    gamma0, //distortion_ 1.3
                                                    beta,   //beta_ 2.5
                                                    delta, //delta_ 0.3
                                                    60,
                                                    approx);//approx


  Constrain_surface_3_polyhedral<K>* pdomain = new Constrain_surface_3_polyhedral<K>(argv[1], epsilon);
  Hyperbolic_shock_metric_field<K>* metric_field = new Hyperbolic_shock_metric_field<K>(0.6, 0.001);

  Cell_star_set_3<K> starset(criteria, metric_field, pdomain, nb);
  starset.dump();
  starset.refine_all();
  std::ofstream fx("cell_mesher_poly.mesh");
  starset.output_medit(fx);

  delete pdomain;
  return 0;
}


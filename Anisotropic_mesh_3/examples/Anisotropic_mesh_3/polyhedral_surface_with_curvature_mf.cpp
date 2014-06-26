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

//#define CGAL_PROFILE

#include <CGAL/Default_configuration.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Robust_circumcenter_traits_3.h>

#include <CGAL/Anisotropic_surface_mesher_3.h>
#include <CGAL/Starset.h>

#include <CGAL/Constrain_surface_3_polyhedral.h>
#include <CGAL/Euclidean_metric_field.h>
#include <CGAL/Polyhedral_curvature_metric_field.h>

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel	Epick;
typedef CGAL::Robust_circumcenter_traits_3<Epick> K;

typedef Anisotropic_surface_mesher_3<K> Mesher;

std::string output_filename(const std::string& filename)
{
  std::size_t found = filename.find_last_of("/\\");
  std::string file = filename.substr(found+1, filename.length()-5);
  file.append("_output.off");
  return file;
}

int main(int argc, char *argv[])
{
  CGAL::default_random = CGAL::Random(0);
  std::string file(argv[1]);
  std::cout << "Let's mesh : " << file << std::endl;

  int n = 2;

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
  const double en_factor = (argc > n) ? atof(argv[n++]) : 0.999;

  Criteria_base<K>* criteria = new Criteria_base<K>(approx, f_r0, f_rho0,
                                                    c_r0, c_rho0, sliverity,
                                                    gamma0, beta, delta, nb,
                                                    max_times_to_try);

  Constrain_surface_3_polyhedral<K>* pdomain = new Constrain_surface_3_polyhedral<K>(argv[1]);

  Polyhedral_curvature_metric_field<K>* metric_field =
    new Polyhedral_curvature_metric_field<K>(*pdomain, epsilon, en_factor);
  //Euclidean_metric_field<K>* metric_field = new Euclidean_metric_field<K>();

  Starset<K> starset;
  Anisotropic_surface_mesher_3<K> mesher(starset, pdomain, criteria, metric_field);

  mesher.refine_mesh();

  std::string outfile = output_filename(file);
  output_surface_off(starset, outfile.c_str());
  
  delete metric_field;
  delete pdomain;
  delete criteria;

  return 0;
}


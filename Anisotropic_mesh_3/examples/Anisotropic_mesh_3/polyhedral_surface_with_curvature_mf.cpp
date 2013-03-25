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

#include <CGAL/Default_configuration.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrain_surface_3_polyhedral.h>
#include <CGAL/Polyhedral_curvature_metric_field.h>
#include <CGAL/Anisotropic_surface_mesher_3.h>

#include <CGAL/Euclidean_metric_field.h>
#include <Metric_field/Arctan_metric_field.h>


#define CGAL_PROFILE

using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel	K;
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

  double r0 = (argc > 2) ? atof(argv[2]) : 1.;
  double gamma0 = (argc > 3) ? atof(argv[3]) : 1.8;
  double approx = (argc > 4) ? atof(argv[4]) : 0.;
  double epsilon = (argc > 5) ? atof(argv[5]) : (0.1);
  
  Criteria_base<K>* criteria = new Criteria_base<K>(0., //radius_edge_ratio_ //3.0
                                                    0., //sliverity_ (not used)
                                                    r0, //circumradius_
                                                    gamma0, //distortion_
                                                    2.5, //beta_
                                                    0.3, //delta_
                                                    60, //nb tries in pick_valid
                                                    approx);

  Constrain_surface_3_polyhedral<K>* pdomain
    = new Constrain_surface_3_polyhedral<K>(argv[1], epsilon);

  Polyhedral_curvature_metric_field<K>* metric_field =
    new Polyhedral_curvature_metric_field<K>(*pdomain, epsilon);
  
  //Criteria_base<K>* criteria = new Criteria_base<K>(3.0, 0.2, r0, gamma0, 2.5, 0.3, approx);
  //Euclidean_metric_field<K>* metric_field = new Euclidean_metric_field<K>();

  Surface_star_set_3<K> starset(criteria, metric_field, pdomain);

  std::cout << "Max curvature : " << pdomain->global_max_curvature() << std::endl;
  std::cout << "Min curvature : " << pdomain->global_min_curvature() << std::endl;

  Mesher mesher(starset);
  mesher.refine_all();

  std::string outfile = output_filename(file);
  mesher.output(outfile);
  
  delete pdomain;
  return 0;
}


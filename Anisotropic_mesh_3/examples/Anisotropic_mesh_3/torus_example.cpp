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
  std::ofstream fx("torus_info.txt");

  Timer timer;
  CGAL::default_random = CGAL::Random(0);

  if(argc == 2)
  {
    std::cout << "torus_example parameters :" << std::endl;
    std::cout << "R=10 r=1 epsilon=1 r0=1 gamma0=1.5 rho0=3 approx=0";
    std::cout << " nb=20 beta=2.5 delta=0.3 gamma1=2 y=5.5 xcondition=-1 " << std::endl;
    return 0;
  }
  
  K::FT R = (argc > 1) ? atof(argv[1]) : 10.;
  K::FT r = (argc > 2) ? atof(argv[2]) : 1.;

  K::FT epsilon = (argc > 3) ? atof(argv[3]) : 1.;

  K::FT r0 = (argc > 4) ? atof(argv[4]) : 1.0/*default*/;
  K::FT gamma0 = (argc > 5) ? atof(argv[5]) : 1.5;
  K::FT rho0 = (argc > 6) ? atof(argv[6]) : 3.0;
  K::FT approx = (argc > 7) ? atof(argv[7]) : 0.;

  int nb = (argc > 8) ? atoi(argv[8]) : 20;

  K::FT beta = (argc > 9) ? atof(argv[9]) : 2.5;
  K::FT delta = (argc > 10) ? atof(argv[10]) : 0.3;
  //m_distortion_bound_avoid_pick_valid
  K::FT distortion_bound = (argc > 11) ? atof(argv[11]) : 2.;

  Criteria_base<K>* criteria = new Criteria_base<K>(rho0, //radius_edge_ratio_
                                                    0.2,    //sliverity_
                                                    r0,     //circumradius_ 0.1
                                                    gamma0, //distortion_ 1.3
                                                    beta,   //beta_ 2.5
                                                    delta,  //delta_ 0.3
                                                    60,     //max_times_to_try_in_picking_region_
                                                    approx); //approximation

  fx << "Torus :" << std::endl;
  fx << "\tR = " << R << std::endl;
  fx << "\tr = " << r << std::endl;
  fx << "\teps = " << epsilon << std::endl;
  fx << "\tr0 = " << r0 << std::endl;
  fx << "\tgamma0 = " << gamma0 << std::endl;
  fx << "\trho0 = " << rho0 << std::endl;
  fx << "\tapprox = " << approx << std::endl;
  fx << "\tbeta = " << beta << std::endl;
  fx << "\tdelta = " << delta << std::endl;
  fx << "\tgamma1 = " << distortion_bound << std::endl;
  
  Constrain_surface_3_torus<K>* pdomain
    = new Constrain_surface_3_torus<K>(R, r);
  Implicit_curvature_metric_field<K>* metric_field =
    new Implicit_curvature_metric_field<K>(*pdomain, epsilon);
  //Torus_metric_field<K>* metric_field = new Torus_metric_field<K>(R, r, epsilon);
  //Euclidean_metric_field<K>* metric_field = new Euclidean_metric_field<K>();

  std::string file = output_filename(R, r);
  if(argc > 12)
  {
    K::FT y = (argc > 12) ? atof(argv[12]) : (0.5*(R + r));
    int xcondition = (argc > 13) ? atoi(argv[13]) : -1;//default : no condition on x
    K::Plane_3 plane1(0., 1., 0., 0.);
    K::Plane_3 plane2(1.0, 0.0001, 0., 0.);
    typedef Is_between<K::Plane_3, K::Point_3> RCondition;
    RCondition condition(plane1, plane2, (xcondition == 1));

    fx << "\ty_condition = " << y << std::endl;
    fx << "\tx_condition = " << (xcondition == 1) << std::endl;
    fx << std::endl << "nbV" << "\t" << "time" << std::endl;
    timer.start();
    Surface_star_set_3<K, RCondition> starset(criteria, metric_field, pdomain, nb,                                               
                                              distortion_bound, condition); 
    timer.stop();
    fx << starset.number_of_stars() << "\t" <<  timer.time() << std::endl;
  
    timer.start();
    starset.refine_all(/*max nb of points*/);
    timer.stop();
    fx << starset.number_of_stars() << "\t" <<  timer.time() << std::endl;
  
    starset.output_surface_star_set("torus_example.mesh");//mesh
    starset.output(file.c_str());//off
  }
  else
  {
    timer.start();
    Surface_star_set_3<K> starset(criteria, metric_field, pdomain, nb, distortion_bound); 
    timer.stop();
    fx << starset.number_of_stars() << "\t" <<  timer.time() << std::endl;
  
    timer.start();
    starset.refine_all(/*max nb of points*/);
    timer.stop();
    fx << starset.number_of_stars() << "\t" <<  timer.time() << std::endl;
    
    starset.output(file.c_str());
  }

  delete pdomain;
  return 0;
}

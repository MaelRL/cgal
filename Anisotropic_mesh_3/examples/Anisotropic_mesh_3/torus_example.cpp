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
//#include <CGAL/Euclidean_metric_field.h>
//#include <Metric_field/Torus_metric_field.h>
//#include <Metric_field/Cylinder_metric_field.h>

//#include <Domain/Constrain_surface_3_torus.h>
//#include <Domain/Constrain_surface_3_sphere.h>
//#include <Domain/Constrain_surface_3_ellipse.h>
//#include <Domain/Constrain_surface_3_tanglecube.h>
#include <Domain/Constrain_surface_3_cylinder.h>

#include <CGAL/Anisotropic_surface_mesher_3.h>


using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel	K;
typedef CGAL::Timer Timer;


std::string output_filename(const double& R,
                            const double& r)
{
  std::ostringstream oss;
  oss << "torus_" << R << "_" << r << ".off";
  return oss.str();
}


template<typename PlaneType, typename PointType>
struct Is_between
{
public:
  bool operator()(const PointType& p) const
  {
    if(xcondition && p.x() < 0) 
      return false;
    PointType pp1 = plane1.projection(p);
    PointType pp2 = plane2.projection(p);
    double scal = (p.x() - pp1.x()) * (p.x() - pp2.x())
      + (p.y() - pp1.y()) * (p.y() - pp2.y())
      + (p.z() - pp1.z()) * (p.z() - pp2.z());
    return (scal <= 0.);
  }
public:
  Is_between(const PlaneType& p1, const PlaneType& p2, const bool xc)
    : plane1(p1), plane2(p2), xcondition(xc){}
  Is_between(const Is_between& ib)
    : plane1(ib.plane1), plane2(ib.plane2), xcondition(ib.xcondition){}
private:
  PlaneType plane1;
  PlaneType plane2;
  bool xcondition;
};


int main(int argc, char* argv[])
{
  std::ofstream fx("torus_timings.txt");

#ifdef ANISO_USE_EIGEN
  std::cout << "Use Eigen" << std::endl;
#else
  std::cout << "Don't use Eigen" << std::endl;
#endif
  Timer timer;
  CGAL::default_random = CGAL::Random(0);

  K::FT R = (argc > 1) ? atof(argv[1]) : 10.;
  K::FT r = (argc > 2) ? atof(argv[2]) : 1.;

  K::FT lambda = (argc > 3) ? atof(argv[3]) : 1.;

  K::FT r0 = (argc > 4) ? atof(argv[4]) : 1.0/*default*/;
  K::FT gamma0 = (argc > 5) ? atof(argv[5]) : 1.5;

  int nb = (argc > 6) ? atoi(argv[6]) : 10;

  fx << "Torus :" << std::endl;
  fx << "\tR = " << R << std::endl;
  fx << "\tr = " << r << std::endl;
  fx << "\teps = " << lambda << std::endl;
  fx << "\tr0 = " << r0 << std::endl;
  fx << "\tgamma0 = " << gamma0 << std::endl;

  Criteria_base<K> criteria(3.0, //radius_edge_ratio_
    0.2, //sliverity_
    r0, //circumradius_ 0.1
    gamma0, //distortion_ 1.3
    2.5, //beta_
    0.3);//delta_

  fx << std::endl << "nbV" << "\t" << "time" << std::endl;
  timer.start();
  //Constrain_surface_3_torus<K>* pdomain
  //  = new Constrain_surface_3_torus<K>(R, r);
//Constrain_surface_3_sphere<K>* pdomain
//  = new Constrain_surface_3_sphere<K>();
//  Constrain_surface_3_ellipse<K>* pdomain
//    = new Constrain_surface_3_ellipse<K>(R/*a*/, r/*b*/, 1./*c*/);
//  Constrain_surface_3_tanglecube<K>* pdomain
//    = new Constrain_surface_3_tanglecube<K>();

  K::FT h = r;
  Constrain_surface_3_cylinder<K>* pdomain
     = new Constrain_surface_3_cylinder<K>(R/*r*/, h);

  Implicit_curvature_metric_field<K> metric_field(*pdomain, lambda);
//  Torus_metric_field<K> metric_field(R, r, lambda);
//  Euclidean_metric_field<K> metric_field;
//  Cylinder_metric_field<K> metric_field(R/*r*/, r/*h*/, lambda/*epsilon*/);

//  K::FT y = (argc > 7) ? atof(argv[7]) : (R + r + 1);
  int xcondition = (argc > 8) ? atoi(argv[8]) : -1;//default : no condition on x
  
  K::Plane_3 plane1(0., 0., 1., -0.1*h);
  K::Plane_3 plane2(0., 0., 1., -0.9*h);

  typedef Is_between<K::Plane_3, K::Point_3> RCondition;
  //typedef Anisotropic_surface_mesher_3<K, RCondition> Mesher;

  RCondition condition(plane1, plane2, (xcondition == 1));
  Surface_star_set_3<K, RCondition> starset(criteria, metric_field, pdomain, nb, condition);
  //Surface_star_set_3<K> starset(criteria, metric_field, pdomain, nb);
  timer.stop();
  
  fx << starset.number_of_stars() << "\t" <<  timer.time() << std::endl;
  starset.refine_all();

  //starset.refine_all(fx, starttime);

  //timer.stop();
  //std::cerr << timer.time() << " seconds" << std::endl;

  //std::string file = output_filename(R, r);
  //mesher.output(file);
  starset.output("closed_cylinder.off");

  delete pdomain;
  return 0;
}

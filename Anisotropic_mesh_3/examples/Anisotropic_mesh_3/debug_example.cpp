#include <CGAL/Default_configuration.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

//#include <Domain/Constrain_surface_3_torus.h>
//#include <CGAL/Implicit_curvature_metric_field.h>
//#include <Metric_field/Torus_metric_field.h>
#include <CGAL/Polyhedral_curvature_metric_field.h>

#include <iostream>
#include <fstream>


using namespace CGAL::Anisotropic_mesh_3;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;


Point_3 point_on_torus(double R, double r, //radii
                       double u, double v) //angles
{
  return Point_3((R + r * std::cos(v)) * std::cos(u),
                 (R + r * std::cos(v)) * std::sin(u),
                 r * std::sin(v));
}


int main()
{
  K::FT R = 1;
  K::FT r = 0.25;
  K::FT lambda = 1.0;
  //Constrain_surface_3_torus<K>* pdomain
  //  = new Constrain_surface_3_torus<K>(R, r);

  char* filename = "fertility.off";//"torus_original.off";

  /*torus polyhedron*/
  std::ifstream input(filename);
  if(!input)
  {
    std::cout << "\nError : file does not exist\n";
    return 0;
  }    
  Polyhedron poly;
  input >> poly;


  /*Estimation of curvature on implicit torus*/
  //Implicit_curvature_metric_field<K> mf_implicit(*pdomain, lambda);

  //std::ofstream file_implicit_min("vectors_implicit_min.txt");
  //std::ofstream file_implicit_max("vectors_implicit_max.txt");

  /*Explicit curvature on implicit torus*/
  //Torus_metric_field<K> mf_torus(R, r, lambda);

  //std::ofstream file_torus_min("vectors_torus_min.txt");
  //std::ofstream file_torus_max("vectors_torus_max.txt");

  /*Estimation of curvature on polyhedral torus*/
  Constrain_surface_3_polyhedral<K>* p_poly
    = new Constrain_surface_3_polyhedral<K>(filename);    
  Polyhedral_curvature_metric_field<K> mf_poly(*p_poly, lambda);

  std::ofstream file_poly_min("vectors_poly_min.txt");
  std::ofstream file_poly_max("vectors_poly_max.txt");
  std::ofstream file_poly_third("vectors_poly_third.txt");

  /*iterator on polyhedron vertices*/
  std::cout << "Polyhedron : ok\n";
  std::cerr << "Compute vectors...";
  for(Polyhedron::Vertex_iterator v = poly.vertices_begin(); 
      v != poly.vertices_end(); 
      ++v)
  {
    std::cout << ".";
    const Point_3& p = v->point();
    Vector_3 eigen;
    //---------- estimation of curvature (implicit) ------------------
    //Implicit_curvature_metric_field<K>::Metric mf1 = mf_implicit.compute_metric(p);
    //mf1.get_min_eigenvector(eigen);
    //file_implicit_min << "2 " << p << " " << (p + 0.1*eigen) << std::endl;

    //mf1.get_max_eigenvector(eigen);
    //file_implicit_max << "2 " << p << " " << (p + 0.1*eigen) << std::endl;

    //---------- explicit torus curvature ------------------
    //Torus_metric_field<K>::Metric mf2 = mf_torus.compute_metric(p);
    //mf2.get_min_eigenvector(eigen);
    //file_torus_min << "2 " << p << " " << (p + 0.1*eigen) << std::endl;

    //mf2.get_max_eigenvector(eigen);
    //file_torus_max << "2 " << p << " " << (p + 0.1*eigen) << std::endl;

    //---------- estimation of curvature (polyhedron) ------------------
    Polyhedral_curvature_metric_field<K>::Metric mf3 = mf_poly.compute_metric(p);
    mf3.get_min_eigenvector(eigen);
    file_poly_min << "2 " << p << " " << (p + 0.1*eigen) << std::endl;

    mf3.get_max_eigenvector(eigen);
    file_poly_max << "2 " << p << " " << (p + 0.1*eigen) << std::endl;

    mf3.get_third_eigenvector(eigen);
    file_poly_third << "2 " << p << " " << (p + 0.1*eigen) << std::endl;
  }
  std::cout << "done." << std::endl;

  //file_implicit_min.close();
  //file_implicit_max.close();
  //file_torus_min.close();
  //file_torus_max.close();
  file_poly_min.close();
  file_poly_max.close();
  file_poly_third.close();

//  std::cout << "\nEigenvalues :\n";
//  Point_3 pn = point_on_torus(R, r, 0.1, 0.5*CGAL_PI);
//  Point_3 ps = point_on_torus(R, r, 0.1, 1.5*CGAL_PI);
//  Point_3 pe = point_on_torus(R, r, 0.1, 0.);
//  Point_3 po = point_on_torus(R, r, 0.1, CGAL_PI);
//
//  Implicit_curvature_metric_field<K>::Metric mfa
//    = mf_implicit.compute_metric(pn);
//  std::cout << "North : \n";
//  std::cout << "   c_min = " <<  ((mfa.get_min_eigenvalue() - 1.) / lambda) << std::endl;
//  std::cout << "   c_max = " <<  ((mfa.get_max_eigenvalue() - 1.) / lambda) << std::endl;
//
//  Implicit_curvature_metric_field<K>::Metric mfb
//    = mf_implicit.compute_metric(ps);
//  std::cout << "South : \n";
//  std::cout << "   c_min = " <<  ((mfb.get_min_eigenvalue() - 1.) / lambda) << std::endl;
//  std::cout << "   c_max = " <<  ((mfb.get_max_eigenvalue() - 1.) / lambda) << std::endl;
//
//  Implicit_curvature_metric_field<K>::Metric mfc
//   = mf_implicit.compute_metric(pe);
//  std::cout << "East : \n";
//  std::cout << "   c_min = " <<  ((mfc.get_min_eigenvalue() - 1.) / lambda) << std::endl;
//  std::cout << "   c_max = " <<  ((mfc.get_max_eigenvalue() - 1.) / lambda) << std::endl;
//
//  Implicit_curvature_metric_field<K>::Metric mfd
//   = mf_implicit.compute_metric(po);
//  std::cout << "West : \n";
//  std::cout << "   c_min = " <<  ((mfd.get_min_eigenvalue() - 1.) / lambda) << std::endl;
//  std::cout << "   c_max = " <<  ((mfd.get_max_eigenvalue() - 1.) / lambda) << std::endl;

  return 0;
}

#include <CGAL/Default_configuration.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <Domain/Constrain_surface_3_torus.h>
#include <CGAL/Implicit_curvature_metric_field.h>
#include <Metric_field/Torus_metric_field.h>
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

Point_3 point_on_ellipsoid  (double a, double b, double c, //radii
                       double u, double v) //angles
{
  return Point_3(a*std::cos(u)*std::cos(v),
                 b*std::cos(u)*std::sin(v),
                 c*std::sin(u));
}

int main()
{
  K::FT R = 10.;
  K::FT r = 1.;
  K::FT epsilon = 0.;

  Constrain_surface_3_torus<K>* pdomain = new Constrain_surface_3_torus<K>(R, r);

  //char* filename = "fertility.off";//"torus_original.off";

  /*torus polyhedron*/
  //std::ifstream input(filename);
  //if(!input)
  //{
  //  std::cout << "\nError : file does not exist\n";
  //  return 0;
  //}
  //Polyhedron poly;
  //input >> poly;


  /*Estimation of curvature on implicit torus*/
  Implicit_curvature_metric_field<K> mf_implicit(*pdomain, epsilon);

  //std::ofstream file_implicit_min("vectors_implicit_min.txt");
  //std::ofstream file_implicit_max("vectors_implicit_max.txt");

  /*Explicit curvature on implicit torus*/
  Torus_metric_field<K> mf_torus(R, r, epsilon);

  //std::ofstream file_torus_min("vectors_torus_min.txt");
  //std::ofstream file_torus_max("vectors_torus_max.txt");

  /*Estimation of curvature on polyhedral torus*/

/*
  ------------------------ BAD STUFF ------------------------------

  Constrain_surface_3_polyhedral<K>* p_poly = new Constrain_surface_3_polyhedral<K>(filename);
  Polyhedral_curvature_metric_field<K> mf_poly(*p_poly, epsilon);

  std::ofstream file_poly_min("vectors_poly_min.txt");
  std::ofstream file_poly_max("vectors_poly_max.txt");
  std::ofstream file_poly_third("vectors_poly_third.txt");

  //------------ iterator on polyhedron vertices ----------------
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

  file_implicit_min.close();
  file_implicit_max.close();
  file_torus_min.close();
  file_torus_max.close();
  file_poly_min.close();
  file_poly_max.close();
  file_poly_third.close();
*/

  // --------------------------------- TORUS STUFF --------------------------------

  std::cout << "\nEigenvalues, Eigenvectors & stuff :\n";

  std::cout << "R : " << R << " r : " << r << " eps : " << epsilon << std::endl;
  std::cout << "1/R : " << 1./R << " | 1/r : " << 1./r << std::endl;


  int point_nb = 10;
  double shift = 2.*CGAL_PI/(point_nb + 1.);

  std::cout << "shift : " << shift << std::endl;

  std::vector<Point_3> points(point_nb);
  std::vector<Implicit_curvature_metric_field<K>::Metric> point_imp_metrics(point_nb);
  std::vector<Torus_metric_field<K>::Metric> point_exp_metrics(point_nb);

  Vector_3 v_min, v_max, v_en;

  for(int i=0; i<point_nb; ++i){

    points[i] = point_on_torus(R, r, 2.*i*CGAL_PI/point_nb, 0.*CGAL_PI);

    std::cout << "Point n° " << i << " : " << R << " " << r << " 0.1 " << 2*i*CGAL_PI/point_nb << std::endl;
    std::cout << " coordinates : " << points[i].x() << " " << points[i].y() << " " << points[i].z() << std::endl;

//    point_imp_metrics[i] = mf_implicit.compute_metric(points[i]);
    point_exp_metrics[i] = mf_torus.compute_metric(points[i]);

//    std::cout << "Implicit : " << std::endl;
//    std::cout << " e_min = " <<  (point_imp_metrics[i]).get_min_eigenvalue() << std::endl;
//    std::cout << " e_max = " <<  (point_imp_metrics[i]).get_max_eigenvalue() << std::endl;
//    std::cout << " e_max = " <<  (point_imp_metrics[i]).get_third_eigenvalue() << std::endl;
//    (point_imp_metrics[i]).get_min_eigenvector(v_min);
//    (point_imp_metrics[i]).get_max_eigenvector(v_max);
//    (point_imp_metrics[i]).get_third_eigenvector(v_en);
//    std::cout << " v_min = " << v_min << std::endl;
//    std::cout << " v_max = " << v_max << std::endl;
//    std::cout << " v_en = " << v_en << std::endl;

    std::cout << "Explicit : " << std::endl;
    std::cout << " e_min = " <<  (point_exp_metrics[i]).get_min_eigenvalue() << std::endl;
    std::cout << " e_max = " <<  (point_exp_metrics[i]).get_max_eigenvalue() << std::endl;
    std::cout << " e_max = " <<  (point_exp_metrics[i]).get_third_eigenvalue() << std::endl;
    (point_exp_metrics[i]).get_min_eigenvector(v_min);
    (point_exp_metrics[i]).get_max_eigenvector(v_max);
    (point_exp_metrics[i]).get_third_eigenvector(v_en);
    std::cout << " v_min = " << v_min << std::endl;
    std::cout << " v_max = " << v_max << std::endl;
    std::cout << " v_en = " << v_en << std::endl;
  }

//  std::cout << std::endl << "Distortion Computations" << std::endl;
//  std::cout << "Implicit : " << std::endl;
//  for(int i=0; i<point_nb; ++i){
//    int next_i = (i == (point_nb-1))?0:(i+1);
//    double d = (point_imp_metrics[i]).compute_distortion( (point_imp_metrics[next_i]) );
//    std::cout << i << " " << next_i << " | " << d << std::endl;
//  }

  std::cout << "Explicit : " << std::endl;
  for(int i=0; i<point_nb; ++i){
    int next_i = (i == (point_nb-1))?0:(i+1);
    double d = (point_exp_metrics[i]).compute_distortion( (point_exp_metrics[next_i]) );
    std::cout << i << " " << next_i << " | " << d << std::endl;
  }

  // --------------------------------- ELLI STUFF --------------------------------

/*
  K::FT a = 10.;
  K::FT b = 1.;
  K::FT c = 1.;
  K::FT epsilon = 0.;

  Constrain_surface_3_ellipse<K>* pdomain = new Constrain_surface_3_ellipse<K>(a, b, c);
  Implicit_curvature_metric_field<K> mf_implicit(*pdomain, epsilon);

  std::cout << "\nEigenvalues, Eigenvectors & stuff :\n";

  std::cout << "a : " << a << " b : " << b << " c : " << c << " eps : " << epsilon << std::endl;

  int point_nb = 10;
  std::vector<Point_3> points(point_nb);
  std::vector<Implicit_curvature_metric_field<K>::Metric> point_imp_metrics(point_nb);

  Vector_3 v_min, v_max, v_en;

  for(int i=0; i<point_nb; ++i){
    points[i] = point_on_ellipsoid(a, b, c,-0.5*CGAL_PI + i/point_nb*CGAL_PI, -CGAL_PI + 0.3);
    //u is in [-pi/2;pi/2] & v is in [-pi;pi]

    std::cout << "Point n° " << i << " : " << R << " " << r << " 0.1 " << 2*i*CGAL_PI/point_nb << std::endl;
    std::cout << " coordinates : " << points[i].x() << " " << points[i].y() << " " << points[i].z() << std::endl;

    point_imp_metrics[i] = mf_implicit.compute_metric(points[i]);

    std::cout << "Implicit : " << std::endl;
    std::cout << " e_min = " <<  (point_imp_metrics[i]).get_min_eigenvalue() << std::endl;
    std::cout << " e_max = " <<  (point_imp_metrics[i]).get_max_eigenvalue() << std::endl;
    std::cout << " e_max = " <<  (point_imp_metrics[i]).get_third_eigenvalue() << std::endl;
    (point_imp_metrics[i]).get_min_eigenvector(v_min);
    (point_imp_metrics[i]).get_max_eigenvector(v_max);
    (point_imp_metrics[i]).get_third_eigenvector(v_en);
    std::cout << " v_min = " << v_min << std::endl;
    std::cout << " v_max = " << v_max << std::endl;
    std::cout << " v_en = " << v_en << std::endl;
  }

  std::cout << std::endl << "Distortion Computations" << std::endl;
  std::cout << "Implicit : " << std::endl;
  for(int i=0; i<point_nb; ++i){
    int next_i = (i == (point_nb-1))?0:(i+1);
    double d = (point_imp_metrics[i]).compute_distortion( (point_imp_metrics[next_i]) );
    std::cout << i << " " << next_i << " | " << d << std::endl;
  }
*/

  return 0;
}

#ifndef CGAL_ANISOTROPIC_MESH_3_CANVAS_MESH_3_H
#define CGAL_ANISOTROPIC_MESH_3_CANVAS_MESH_3_H

#include <CGAL/Canvas/canvas_triangulation_io.h>

#include <CGAL/Metric.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/make_mesh_3.h>

#include <iostream>
#include <fstream>
#include <list>
#include <vector>

namespace CGAL
{
// functions to generate the canvas using Mesh_3 (Aniso_mesh_3 is a turtle).
// Need to define our own refinement criterion based on metric distortion

template<typename K>
typename K::FT cube_function(const typename K::Point_3& p)
{
  if( p.x() >= 0. && p.x() <= 2. &&
      p.y() >= 0. && p.y() <= 2. &&
      p.z() >= 0. && p.z() <= 2. )
    return -1.;
  return 1.;
}

template<typename K>
void create_polylines(std::list<std::vector<typename K::Point_3> >& polylines)
{
  typedef typename K::FT                      FT;
  typedef typename K::Point_3                 Point_3;
  typedef std::vector<Point_3>                Polyline_3;

  FT side_x = 2.; // fixme don't hardcode stuff like that...
  FT side_y = 2.;
  FT side_z = 2.;

  // bot 4 edges
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(0, 0, 0));
    polyline.push_back(Point_3(side_x, 0, 0));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(side_x, 0, 0));
    polyline.push_back(Point_3(side_x, side_y, 0));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(side_x, side_y, 0));
    polyline.push_back(Point_3(0, side_y, 0));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(0, side_y, 0));
    polyline.push_back(Point_3(0, 0, 0));
    polylines.push_back(polyline);
  }

  // vertical edges
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(0, 0, 0));
    polyline.push_back(Point_3(0, 0, side_z));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(side_x, 0, 0));
    polyline.push_back(Point_3(side_x, 0, side_z));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(0, side_y, 0));
    polyline.push_back(Point_3(0, side_y, side_z));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(side_x, side_y, 0));
    polyline.push_back(Point_3(side_x, side_y, side_z));
    polylines.push_back(polyline);
  }

  // top 4 edges
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(0, 0, side_z));
    polyline.push_back(Point_3(side_x, 0, side_z));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(side_x, 0, side_z));
    polyline.push_back(Point_3(side_x, side_y, side_z));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(side_x,side_y,side_z));
    polyline.push_back(Point_3(0,side_y,side_z));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(0,side_y,side_z));
    polyline.push_back(Point_3(0,0,side_z));
    polylines.push_back(polyline);
  }
}

template <typename Tr>
class Mesh_edge_criteria_w_distortion_3 :
    public Mesh_edge_criteria_3<Tr>

{
  //todo
};

template<typename Tr>
class Mesh_facet_criteria_w_distortion_3 :
    public Mesh_facet_criteria_3<Tr>
{
  //todo
};

template<typename Tr>
class Mesh_cell_criteria_w_distortion_3 :
    public Mesh_cell_criteria_3<Tr>
{
  //todo
};

// Sizing field
template<typename Tr, typename MF>
struct Geo_sizing_field
{
  typedef typename Tr::Geom_traits                               Gt;
  typedef typename Gt::FT                                        FT;
  typedef typename Gt::Point_3                                   Point_3;

  typedef FT (Function)(const Point_3&);
  typedef Implicit_mesh_domain_3<Function, Gt>                   Mesh_domain;
  typedef typename Mesh_domain::Index                            Index;

  typedef typename MF::Metric                                    Metric;

  const MF& mf;

  FT operator()(const Point_3& p, const int, const Index&) const
  {
    const FT base = 10.;

    const Metric& m = mf.compute_metric(p);
    const FT e_max = m.get_max_eigenvalue();
    const FT e_min = m.get_min_eigenvalue();
    const FT e_third = m.get_third_eigenvalue();
    CGAL_precondition(e_max >= e_min);
    const FT max = (std::max)(e_max, e_third);
    const FT width = 1. / max;

    const FT discretization = 6.;
    const FT metric_based_size = width / discretization;

//    std::cout << metric_based_size << std::endl;

    return (std::min)(base, metric_based_size);
  }

  Geo_sizing_field(const MF& mf_) : mf(mf_) { }
};

template<typename K, typename MF>
void generate_canvas(const MF& mf)
{
  typedef typename K::FT                                         FT;
  typedef typename K::Point_3                                    Point_3;

  typedef std::vector<Point_3>                                   Polyline_3;
  typedef std::list<Polyline_3>                                  Polylines;

  typedef FT (Function)(const Point_3&);
  typedef Implicit_mesh_domain_3<Function, K>                    Mesh_domain;
  typedef Mesh_domain_with_polyline_features_3<Mesh_domain>      Mesh_domain_with_features;

  typedef typename Mesh_triangulation_3<Mesh_domain>::type       Tr;
  typedef Mesh_complex_3_in_triangulation_3<Tr>                  C3t3;

  typedef Mesh_criteria_3<Tr>                                    Mesh_criteria;
  typedef typename Mesh_criteria::Edge_criteria                  Edge_criteria;
  typedef typename Mesh_criteria::Facet_criteria                 Facet_criteria;
  typedef typename Mesh_criteria::Cell_criteria                  Cell_criteria;

  // Domain
  Mesh_domain_with_features domain(cube_function<K>,
                                   typename K::Sphere_3(Point_3(1., 1., 1.), 5.*5.),
                                   1e-6);

  // Polylines
  Polylines polylines;
  create_polylines<K>(polylines);
  domain.add_features(polylines.begin(), polylines.end());

  // Set mesh criteria
  Geo_sizing_field<Tr, MF> size(mf);
  Edge_criteria edge_criteria(size);
  Facet_criteria facet_criteria(30, size, 0.05); // angle, size, approximation
  Cell_criteria cell_criteria(2., size); // radius-edge ratio, size
  Mesh_criteria criteria(edge_criteria, facet_criteria, cell_criteria);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

  std::cout << "Number of vertices: " << c3t3.triangulation().number_of_vertices() << std::endl;

  // Output
  std::ofstream medit_file("super_dense_base_mesh.mesh");
  output_to_medit_all_cells<C3t3, true/*rebind*/, true/*no patch*/>(medit_file, c3t3);
}

}  // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_CANVAS_MESH_3_H

#ifndef CGAL_ANISOTROPIC_MESH_3_C3T3_CANVAS_GENERATOR_H
#define CGAL_ANISOTROPIC_MESH_3_C3T3_CANVAS_GENERATOR_H

#include <CGAL/Canvas/canvas_config.h>
#include <CGAL/Canvas/canvas_triangulation_io.h>

#include <CGAL/Metric.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/make_mesh_3.h>

#include <Domain/Mesh_3/surface_3_cube.h>

#include <iostream>
#include <fstream>
#include <list>
#include <vector>

namespace CGAL
{
// functions to generate the canvas using Mesh_3
// todo: refinement criterion based on the metric field's distortion

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

  const MF* mf;

  FT operator()(const Point_3& p, const int, const Index&) const
  {
    const FT base = 1.;

    const Metric& m = mf->compute_metric(p);
    const FT e_max = m.get_max_eigenvalue();
    const FT e_min = m.get_min_eigenvalue();
    const FT e_third = m.get_third_eigenvalue();
    CGAL_precondition(e_max >= e_min);
    const FT max = (std::max)(e_max, e_third);
    const FT width = 1. / max;

    const FT discretization = 6.;
    const FT metric_based_size = width / discretization;

    return (std::min)(base, metric_based_size);
  }

  Geo_sizing_field(const MF* mf_) : mf(mf_) { }
};

template<typename K, typename MF>
void generate_canvas(const MF* mf)
{
  typedef typename K::FT                                         FT;
  typedef typename K::Point_3                                    Point;

  typedef FT (Function)(const Point&);
  typedef CGAL::Mesh_domain_with_polyline_features_3<
                CGAL::Implicit_mesh_domain_3<Function,K> >       Mesh_domain;

#ifdef CGAL_CONCURRENT_MESH_3
  typedef CGAL::Mesh_triangulation_3<Mesh_domain,
                                     CGAL::Kernel_traits<Mesh_domain>::Kernel, // Same as sequential
                                     CGAL::Parallel_tag                        // Tag to activate parallelism
                                    >::type                      Tr;
#else
  typedef typename CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
#endif

typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr,
                                                typename Mesh_domain::Corner_index,
                                                typename Mesh_domain::Curve_segment_index> C3t3;

  typedef Mesh_criteria_3<Tr>                                    Mesh_criteria;
  typedef typename Mesh_criteria::Edge_criteria                  Edge_criteria;
  typedef typename Mesh_criteria::Facet_criteria                 Facet_criteria;
  typedef typename Mesh_criteria::Cell_criteria                  Cell_criteria;

  // Domain
  Mesh_domain const * const domain = cube_domain_with_features<Mesh_domain>();

  // Set mesh criteria
  Geo_sizing_field<Tr, MF> size(mf);
  Edge_criteria edge_criteria(size);
  Facet_criteria facet_criteria(30, size, 0.05); // angle, size, approximation
  Cell_criteria cell_criteria(2., size); // radius-edge ratio, size
  Mesh_criteria criteria(edge_criteria, facet_criteria, cell_criteria);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(*domain, criteria);

  std::cout << "Number of vertices: "
            << c3t3.triangulation().number_of_vertices() << std::endl;

  // Output
  std::ofstream medit_file("input_c3t3.mesh");
  output_to_medit_all_cells<C3t3, true/*rebind*/, true/*no patch*/>(medit_file, c3t3);

  delete domain;
}

}  // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_C3T3_CANVAS_GENERATOR_H

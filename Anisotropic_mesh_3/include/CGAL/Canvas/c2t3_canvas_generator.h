#ifndef CGAL_ANISOTROPIC_MESH_3_C2T3_CANVAS_GENERATOR_H
#define CGAL_ANISOTROPIC_MESH_3_C2T3_CANVAS_GENERATOR_H

#include <CGAL/Canvas/canvas_config.h>
#include <CGAL/Canvas/canvas_triangulation_io.h>

#include <CGAL/Metric.h>
#include <Domain/Mesh_3/surface_3_cube.h>
#include <Domain/Mesh_3/surface_3_sphere.h>
#include <Domain/Mesh_3/surface_3_chair.h>

#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>

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
  // todo
};

template<typename Tr>
class Mesh_facet_criteria_w_distortion_3 :
    public Mesh_facet_criteria_3<Tr>
{
  // todo
};

template<typename Tr>
class Mesh_cell_criteria_w_distortion_3 :
    public Mesh_cell_criteria_3<Tr>
{
  // todo
};

// Sizing field
template<typename Tr, typename MF, typename Domain>
struct Geo_sizing_field
{
  typedef typename Tr::Geom_traits                               Gt;
  typedef typename Gt::FT                                        FT;
  typedef typename Gt::Point_3                                   Point_3;

  typedef typename Domain::Index                                 Index;

  typedef typename MF::Metric                                    Metric;

  const MF* mf;

  FT operator()(const Point_3& p, const int, const Index&) const
  {
    const FT base = 1.;

    const Metric& m = mf->compute_metric(p); // could be kept in memory todo
    const FT l_max = 1. / std::sqrt( m.get_max_eigenvalue() );
    const FT l_min = 1. / std::sqrt( m.get_min_eigenvalue() );
    const FT l_third = 1. / std::sqrt( m.get_third_eigenvalue() );
    const FT min_width = std::min(l_max, l_min);

    const FT discretization = 50.;
    const FT metric_based_size = min_width / discretization;

    // fixme needs to take into account the curvature of the domain here
    return (std::min)(base, metric_based_size);
  }

  Geo_sizing_field(const MF* mf_) : mf(mf_) { }
};

template<typename C3t3, typename MF, typename Domain>
void generate_canvas(C3t3& c3t3,
                     const MF* metric_field,
                     const Domain& domain)
{
  typedef typename C3t3::Triangulation                           Triangulation;
  typedef Mesh_criteria_3<Triangulation>                         Mesh_criteria;
  typedef typename Mesh_criteria::Edge_criteria                  Edge_criteria;
  typedef typename Mesh_criteria::Facet_criteria                 Facet_criteria;
  typedef typename Mesh_criteria::Cell_criteria                  Cell_criteria;

  // Set mesh criteria
  Geo_sizing_field<Triangulation, MF, Domain> size(metric_field);
  Edge_criteria edge_criteria(size);
  Facet_criteria facet_criteria(30, size, 0.1); // angle, size, approximation

  // large radius-edge ratio & radius so the cells aren't refined
  Cell_criteria cell_criteria(Anisotropic_mesh_3::FT_inf,
                              Anisotropic_mesh_3::FT_inf); // radius-edge ratio, size
  Mesh_criteria criteria(edge_criteria, facet_criteria, cell_criteria);

  // Mesh generation
  c3t3 = make_mesh_3<C3t3>(domain, criteria,
                           parameters::no_perturb(),
                           parameters::no_exude());

  std::cout << "Number of corners: "
            << c3t3.number_of_vertices_in_complex() << std::endl;

  // Output
  std::ofstream medit_file("input_c2t3.mesh");
  output_to_medit_all_cells<C3t3, true/*rebind*/, false/*no patch*/>(medit_file, c3t3);
}

template<typename C3t3, typename MF, typename MD>
void generate_canvas(C3t3& c3t3, const MF* metric_field)
{
#ifdef IMPLICIT_DOMAIN
  MD const * const domain = sphere_domain<MD>();
//  MD const * const domain = cube_domain_with_features<MD>();
//  MD const * const domain = chair_domain<MD>();
#else
  // Polyhedral
    // bit ugly, would be nicer to use the constrain_surface_3_poly_domain
    // but it uses enriched polyhedron items...
  typedef typename MD::Polyhedron_type Polyhedron;
  std::ifstream input("/home/mrouxell/Data/OFF/fandisk.off");
  if(!input)
    std::cout << "\nWarning : file does not exist" << std::endl;
  Polyhedron poly;
  input >> poly;
  MD * const domain = new MD(poly);
  domain->detect_features();
#endif

  generate_canvas<C3t3, MF, MD>(c3t3, metric_field, *domain);
  delete domain;
}

} // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_C2T3_CANVAS_GENERATOR_H

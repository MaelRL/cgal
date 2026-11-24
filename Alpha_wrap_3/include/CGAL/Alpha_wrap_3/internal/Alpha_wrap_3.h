// Copyright (c) 2019-2023 Google LLC (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Pierre Alliez
//                 Cedric Portaneri,
//                 Mael Rouxel-Labbé
//                 Andreas Fabri
//                 Michael Hemmer
//
#ifndef CGAL_ALPHA_WRAP_3_INTERNAL_ALPHA_WRAP_3_H
#define CGAL_ALPHA_WRAP_3_INTERNAL_ALPHA_WRAP_3_H

#ifdef CGAL_AW3_DEBUG_PP
 #ifndef CGAL_AW3_DEBUG
  #define CGAL_AW3_DEBUG
  #define CGAL_AW3_DEBUG_INITIALIZATION
  #define CGAL_AW3_DEBUG_STEINER_COMPUTATION
  #define CGAL_AW3_DEBUG_QUEUE
  #define CGAL_AW3_DEBUG_FACET_STATUS
  #define CGAL_AW3_DEBUG_MANIFOLDNESS
 #endif
#endif

#include <CGAL/license/Alpha_wrap_3.h>

#include <CGAL/Alpha_wrap_3/internal/Alpha_wrap_triangulation_cell_base_3.h>
#include <CGAL/Alpha_wrap_3/internal/Alpha_wrap_triangulation_vertex_base_3.h>
#include <CGAL/Alpha_wrap_3/internal/Alpha_wrap_AABB_geom_traits.h>
#include <CGAL/Alpha_wrap_3/internal/gate_priority_queue.h>
#include <CGAL/Alpha_wrap_3/internal/geometry_utils.h>
#include <CGAL/Alpha_wrap_3/internal/oracles.h>

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Regular_triangulation_vertex_base_3.h>
#include <CGAL/Regular_triangulation_cell_base_with_weighted_circumcenter_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Robust_weighted_circumcenter_filtered_traits_3.h>

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Default.h>
#include <CGAL/Named_function_parameters.h>
#ifdef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
 #include <CGAL/Modifiable_priority_queue.h>
#endif
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup_extension.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h> // only if non-manifoldness is not treated
#include <CGAL/property_map.h>
#include <CGAL/Real_timer.h>

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <queue>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace CGAL {
namespace Alpha_wraps_3 {
namespace internal {

namespace {

namespace AW3i = ::CGAL::Alpha_wraps_3::internal;

} // unnamed namespace

struct Wrapping_default_visitor
{
  Wrapping_default_visitor() { }

  template <typename AlphaWrapper>
  void on_alpha_wrapping_begin(const AlphaWrapper&) const { }

  template <typename AlphaWrapper>
  void on_flood_fill_begin(const AlphaWrapper&) const { }

  // Whether the flood filling process should continue
  template <typename Wrapper>
  constexpr bool go_further(const Wrapper&) const { return true; }

  template <typename AlphaWrapper, typename Facet>
  bool consider_facet(const AlphaWrapper&, const Facet&) const { return true; }

  template <typename AlphaWrapper, typename Gate>
  void before_facet_treatment(const AlphaWrapper&, const Gate&) const { }

  template <typename AlphaWrapper, typename Point>
  void before_Steiner_point_insertion(const AlphaWrapper&, const Point&) const { }

  template <typename AlphaWrapper, typename VertexHandle>
  void after_Steiner_point_insertion(const AlphaWrapper&, VertexHandle) const { }

  template <typename AlphaWrapper>
  void on_flood_fill_end(const AlphaWrapper&) const { }

  template <typename AlphaWrapper>
  void on_alpha_wrapping_end(const AlphaWrapper&) const { };
};

template <typename Oracle_,
          typename Triangulation_ = CGAL::Default>
class Alpha_wrapper_3
{
  using Oracle = Oracle_;

  // Triangulation
  using Base_GT = typename Oracle::Geom_traits;
  using Default_Gt = CGAL::Robust_circumcenter_filtered_traits_3<Base_GT>;

  using Default_Vbb = Regular_triangulation_vertex_base_3<Default_Gt>;
  using Default_Cbb = Regular_triangulation_cell_base_with_weighted_circumcenter_3<Default_Gt>;
  using Default_Vb = Alpha_wrap_triangulation_vertex_base_3<Default_Gt, Default_Vbb>;
  using Default_Cb = Alpha_wrap_triangulation_cell_base_3<Default_Gt, Default_Cbb>;
  using Default_Cbt = Cell_base_with_timestamp<Default_Cb>; // for determinism
  using Default_Tds = CGAL::Triangulation_data_structure_3<Default_Vb, Default_Cbt>;
  using Default_Triangulation = CGAL::Regular_triangulation_3<Default_Gt, Default_Tds>;

public:
  using Triangulation = typename Default::Get<Triangulation_, Default_Triangulation>::type;

  // Use the geom traits from the triangulation, and trust the (advanced) user that provided it
  using Geom_traits = typename Triangulation::Geom_traits;

private:
  using Cell_handle = typename Triangulation::Cell_handle;
  using Facet = typename Triangulation::Facet;
  using Vertex_handle = typename Triangulation::Vertex_handle;
  using Locate_type = typename Triangulation::Locate_type;

  using Gate = internal::Gate<Triangulation>;

  // A sorted queue is a priority queue sorted by circumradius, and is experimentally significantly
  // slower. However, intermediate results are usable: at each point of the algorithm, the wrap
  // has a somewhat uniform mesh element size.
  //
  // An unsorted queue is a LIFO queue, and is experimentally much faster (~35%),
  // but intermediate results are not useful: a LIFO carves deep before than wide,
  // and can thus for example leave scaffolding faces till almost the end of the refinement.
#ifdef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
  using Alpha_PQ = Modifiable_priority_queue<Gate, Less_gate, Gate_ID_PM<Triangulation>, CGAL_BOOST_PAIRING_HEAP>;
#else
  using Alpha_PQ = std::stack<Gate>;
#endif

  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Weighted_point_3 = typename Geom_traits::Weighted_point_3;
  using Segment_3 = typename Geom_traits::Segment_3;
  using Vector_3 = typename Geom_traits::Vector_3;
  using Sphere_3 = typename Geom_traits::Sphere_3;
  using Ball_3 = typename Geom_traits::Ball_3;
  using Triangle_3 = typename Geom_traits::Triangle_3;
  using Iso_cuboid_3 = typename Geom_traits::Iso_cuboid_3;

  using SC = Simple_cartesian<double>;
  using SC_Point_3 = SC::Point_3;
  using SC_Vector_3 = SC::Vector_3;
  using SC_Iso_cuboid_3 = SC::Iso_cuboid_3;
  using SC2GT = Cartesian_converter<SC, Geom_traits>;

  using Seeds = std::vector<Point_3>;

protected:
  mutable Oracle m_oracle;
  SC_Iso_cuboid_3 m_bbox;

  FT m_alpha = FT(-1), m_sq_alpha = FT(-1);
  FT m_offset = FT(-1), m_sq_offset = FT(-1);

  Seeds m_seeds;

  Triangulation m_tr;

  Alpha_PQ m_queue;

public:
  Alpha_wrapper_3()
#ifdef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
      // '4096' is an arbitrary, not-too-small value for the largest ID in queue initialization
    : m_queue(4096)
#endif
  {
    // Due to the Steiner point computation being a dichotomy, the algorithm is inherently inexact
    // and passing exact kernels is explicitly disabled to ensure no misunderstanding.
    static_assert(std::is_floating_point<FT>::value);
  }

  Alpha_wrapper_3(const Oracle& oracle)
    :
      m_oracle(oracle),
      m_tr(Geom_traits(oracle.geom_traits()))
#ifdef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
      , m_queue(4096)
#endif
  {
    static_assert(std::is_floating_point<FT>::value);
  }

  void clear()
  {
    m_oracle.clear();
    m_bbox = {};
    m_alpha = m_sq_alpha = FT(-1);
    m_offset = m_sq_offset = FT(-1);
    m_seeds.clear();
    m_tr.clear();
    m_queue.clear();
  }

public:
  const Geom_traits& geom_traits() const { return m_tr.geom_traits(); }
  Oracle& oracle() { return m_oracle; }
  const Oracle& oracle() const { return m_oracle; }
  Triangulation& triangulation() { return m_tr; }
  const Triangulation& triangulation() const { return m_tr; }
  Alpha_PQ& queue() { return m_queue; }
  const Alpha_PQ& queue() const { return m_queue; }
  const FT& alpha() const { return m_alpha; }
  const FT& offset() const { return m_offset; }

  double default_alpha() const
  {
    const Bbox_3 bbox = m_oracle.bbox();
    const double diag_length = std::sqrt(square(bbox.xmax() - bbox.xmin()) +
                                         square(bbox.ymax() - bbox.ymin()) +
                                         square(bbox.zmax() - bbox.zmin()));

    return diag_length / 20.;
  }

private:
  const Point_3& circumcenter(const Cell_handle c) const
  {
    return c->circumcenter(geom_traits());
  }

public:
  template <typename OutputMesh,
            typename InputNamedParameters = parameters::Default_named_parameters,
            typename OutputNamedParameters = parameters::Default_named_parameters>
  void operator()(const double alpha, // = default_alpha()
                  const double offset, // = alpha / 30.
                  OutputMesh& output_mesh,
                  const InputNamedParameters& in_np = parameters::default_values(),
                  const OutputNamedParameters& out_np = parameters::default_values())
  {
    namespace PMP = Polygon_mesh_processing;

    using parameters::get_parameter;
    using parameters::get_parameter_reference;
    using parameters::choose_parameter;

    //
    using OVPM = typename CGAL::GetVertexPointMap<OutputMesh, OutputNamedParameters>::type;
    OVPM ovpm = choose_parameter(get_parameter(out_np, internal_np::vertex_point),
                                 get_property_map(vertex_point, output_mesh));

    //
    using Visitor = typename internal_np::Lookup_named_param_def<
                      internal_np::visitor_t,
                      InputNamedParameters,
                      Wrapping_default_visitor // default
                    >::reference;

    Wrapping_default_visitor default_visitor;
    Visitor visitor = choose_parameter(get_parameter_reference(in_np, internal_np::visitor), default_visitor);

    // Points used to create initial cavities
    m_seeds = choose_parameter(get_parameter_reference(in_np, internal_np::seed_points), Seeds());

    // Whether or not some cells should be reflagged as "inside" after the refinement+carving loop
    // as ended, as to ensure that the result is not only combinatorially manifold, but also
    // geometrically manifold.
    //
    // -- Warning --
    // Manifoldness postprocessing will be performed even if the wrapping is interrupted (and
    // this option is enabled).
    const bool do_enforce_manifoldness = choose_parameter(get_parameter(in_np, internal_np::do_enforce_manifoldness), true);

    // Whether to keep pockets of "outside" cells that are not connected to the exterior (or to the
    // initial cavities, if used).
    //
    // -- Warning --
    // Pockets of "outside" cells will be purged even if the wrapping is interrupted (and
    // this option is enabled).
    const bool keep_inner_ccs = choose_parameter(get_parameter(in_np, internal_np::keep_inner_connected_components), false);

    // This parameter enables avoiding recomputing the triangulation from scratch when wrapping
    // the same input for multiple values of alpha (and typically the same offset values).
    //
    // -- Warning --
    // If this is enabled, the 3D triangulation will NOT be re-initialized at launch.
    // This means that the triangulation is NOT cleared, even if:
    // - you use an alpha value that is greater than what was used in a previous run; you will
    //   obtain the same result as the last run.
    // - you use a different offset value between runs, you might then get points that are not
    //   on the offset surface corresponding to that corresponding to the latter offset value.
    const bool refining = choose_parameter(get_parameter(in_np, internal_np::refine_triangulation), false);

#ifdef CGAL_AW3_TIMER
    CGAL::Real_timer t;
    t.start();
#endif

    visitor.on_alpha_wrapping_begin(*this);

    if(!initialize(alpha, offset, refining))
      return;

#ifdef CGAL_AW3_TIMER
    t.stop();
    std::cout << "Initialization took: " << t.time() << " s." << std::endl;
    t.reset();
    t.start();
#endif

#ifdef CGAL_AW3_DEBUG_DUMP_INTERMEDIATE_WRAPS
    dump_triangulation_faces("starting_wrap.off", true /*only_boundary_faces*/);
#endif

    alpha_flood_fill(visitor);

#ifdef CGAL_AW3_DEBUG_DUMP_INTERMEDIATE_WRAPS
    dump_triangulation_faces("flood_filled_wrap.off", true /*only_boundary_faces*/);
#endif

#ifdef CGAL_AW3_TIMER
    t.stop();
    std::cout << "Flood filling took: " << t.time() << " s." << std::endl;
    t.reset();
    t.start();
#endif

    if(do_enforce_manifoldness)
    {
      make_manifold();

#ifdef CGAL_AW3_DEBUG_DUMP_INTERMEDIATE_WRAPS
      dump_triangulation_faces("manifold_wrap.off", true /*only_boundary_faces*/);
#endif

#ifdef CGAL_AW3_TIMER
      t.stop();
      std::cout << "Manifoldness postprocessing took: " << t.time() << " s." << std::endl;
      t.reset();
      t.start();
#endif
    }

    if(!keep_inner_ccs)
    {
      // We could purge *before* manifold enforcement, but making the mesh manifold is
      // very cheap in most cases, so it is better to keep the code simpler.
      purge_inner_connected_components();

#ifdef CGAL_AW3_DEBUG_DUMP_INTERMEDIATE_WRAPS
      dump_triangulation_faces("purged_wrap.off", true /*only_boundary_faces*/);
#endif
    }

 #ifdef CGAL_AW3_DEBUG_DUMP_INTERMEDIATE_WRAPS
    IO::write_polygon_mesh("final_wrap.off", output_mesh, CGAL::parameters::stream_precision(17));
    dump_triangulation_faces("final_tr.off", false /*only_boundary_faces*/);
 #endif

    extract_surface(output_mesh, ovpm, !do_enforce_manifoldness);

#ifdef CGAL_AW3_TIMER
    t.stop();
    std::cout << "Surface extraction took: " << t.time() << " s." << std::endl;
#endif

#ifdef CGAL_AW3_DEBUG
    std::cout << "Alpha wrap vertices:  " << vertices(output_mesh).size() << std::endl;
    std::cout << "Alpha wrap faces:     " << faces(output_mesh).size() << std::endl;

    std::cout << "CREASE VERTICES" << std::endl;
    for(const Vertex_handle v : m_tr.finite_vertex_handles())
      if(v->type() == Vertex_type::CREASE_VERTEX)
        std::cout << " " << m_tr.point(v) << std::endl;

    std::cout << "CORNER VERTICES" << std::endl;
    for(const Vertex_handle v : m_tr.finite_vertex_handles())
      if(v->type() == Vertex_type::CORNER_VERTEX)
        std::cout << " " << m_tr.point(v) << std::endl;
#endif

    visitor.on_alpha_wrapping_end(*this);
  }

  // Convenience overloads
  template <typename OutputMesh>
  void operator()(const double alpha,
                  OutputMesh& output_mesh)
  {
    return operator()(alpha, alpha / 30. /*offset*/, output_mesh);
  }

  template <typename OutputMesh>
  void operator()(OutputMesh& output_mesh)
  {
    return operator()(default_alpha(), output_mesh);
  }

  // This function is public only because it is used in the tests
  SC_Iso_cuboid_3 construct_bbox(const double offset)
  {
    // Input axis-aligned bounding box
    SC_Iso_cuboid_3 bbox = m_oracle.bbox();
    const SC_Point_3 bbox_centroid = midpoint((bbox.min)(), (bbox.max)());

    // Scale a bit to create the initial points not too close to the input
    double scaling = 1.2;
    CGAL::Aff_transformation_3<SC> scale(SCALING, scaling);
    bbox = SC_Iso_cuboid_3(scale.transform((bbox.min)()), scale.transform((bbox.max)()));

    // Translate bbox back to initial centroid
    const SC_Point_3 bbox_transformed_centroid = midpoint((bbox.min)(), (bbox.max)());
    const SC_Vector_3 diff_centroid = bbox_centroid - bbox_transformed_centroid;
    CGAL::Aff_transformation_3<SC> centroid_translate(TRANSLATION, diff_centroid);
    bbox = SC_Iso_cuboid_3(centroid_translate.transform((bbox.min)()),
                           centroid_translate.transform((bbox.max)()));

    // Add the offset
    SC_Vector_3 offset_ext = std::sqrt(3.) * offset * SC_Vector_3(1, 1, 1);
    CGAL::Aff_transformation_3<SC> translate_m(TRANSLATION, - offset_ext);
    CGAL::Aff_transformation_3<SC> translate_M(TRANSLATION,   offset_ext);
    bbox = SC_Iso_cuboid_3(translate_m.transform((bbox.min)()), translate_M.transform((bbox.max)()));

    return bbox;
  }

private:
  // The distinction between inside boundary and outside boundary is the presence of cells
  // being flagged for manifoldness: inside boundary considers those outside, and outside
  // boundary considers them inside.
  bool is_on_inside_boundary(Cell_handle ch, Cell_handle nh) const
  {
    return (ch->is_inside() != nh->is_inside()); // one is "inside", the other is not
  }

  bool is_on_outside_boundary(Cell_handle ch, Cell_handle nh) const
  {
    return (ch->is_outside() != nh->is_outside()); // one is "outside", the other is not
  }

private:
  // Adjust the bbox & insert its corners to construct the starting triangulation
  void insert_bbox_corners()
  {
    m_bbox = construct_bbox(CGAL::to_double(m_offset));

#ifdef CGAL_AW3_DEBUG_INITIALIZATION
    std::cout << "Insert Bbox vertices" << std::endl;
#endif

    // insert in dt the eight corner vertices of the input loose bounding box
    for(int i=0; i<8; ++i)
    {
      const Point_3 bp = SC2GT()(m_bbox.vertex(i));
      Vertex_handle bv = m_tr.insert(bp);
#ifdef CGAL_AW3_DEBUG_INITIALIZATION
      std::cout << "\t" << bp << std::endl;
#endif
      bv->type() = AW3i::Vertex_type::BBOX_VERTEX;
    }
  }

  // Two criteria:
  // - Cells that are intersecting the input are inside
  // - Cells whose circumcenter is in the offset volume are inside: this is because
  // we need to have outside cell circumcenters outside of the volume to ensure
  // that the refinement point is separated from the existing point set.
  Cell_label cavity_cell_label(const Cell_handle ch)
  {
    CGAL_precondition(!m_tr.is_infinite(ch));

    const Tetrahedron_with_outside_info<Geom_traits> tet(ch, geom_traits());
    if(m_oracle.do_intersect(tet))
      return Cell_label::INSIDE;

    const Point_3& ch_cc = circumcenter(ch);
    typename Geom_traits::Construct_ball_3 ball = geom_traits().construct_ball_3_object();
    const Ball_3 ch_cc_offset_ball = ball(ch_cc, m_sq_offset);
    const bool is_cc_in_offset = m_oracle.do_intersect(ch_cc_offset_ball);

    return is_cc_in_offset ? Cell_label::INSIDE : Cell_label::OUTSIDE;
  }

  // Create a cavity using seeds rather than starting from the infinity.
  //
  // For each seed, insert the seeds and its neighboring vertices on a regular icosahedron.
  // The idea behind an icosahedron rather than e.g. a simple tetrahedron is twofold:
  // - Improve the odds of both starting the algorithm, but also of not immediately stopping,
  //   which could happen if a conflict zone included all the initial cavity when using a simple cavity.
  // - Base the cavity diameter on alpha, and allow its cells to intersect the input. If a single
  //   tetrahedron is used, one is forced to make it small-enough such that it does not intersect
  //   the input, and then it could be forced to be smaller than 'alpha' and the algorithm cannot start,
  //   for example if the seed is close to the offset. If we have many tetrahedra, this is far less
  //   likely to happen.
  //
  // For an edge length a, the radius of an icosahedron is
  //   r = 0.5 * a * sqrt(phi * sqrt(5)), phi being the golden ratio
  // which yields
  //   r = a * sin(2pi / 5) =~ 0.95105651629515353 * a
  // Faces of an icosahedron are equilateral triangles with size a and circumradius a / sqrt(3)
  // Since we want faces of the icosahedron to be traversable, we want `a` such that
  //   a / sqrt(3) > alpha
  // Hence r such that
  //   r / (sqrt(3) * sin(2pi/5)) > alpha
  //   r > alpha * sqrt(3) * sin(2pi / 5) ~= 1.6472782070926637 * alpha
  //
  // Furthermore, the triangles between edges of the icosahedron and the center of the icosahedron
  // are not equilateral triangles since a is slightly bigger than r. They are
  // slightly flattened isocele triangles with base `a` and the circumradius is smaller.
  // The circumradius is
  //   r_iso = r² / (2 * h) = r² / (2 * sqrt(r² - (a / 2)²))
  // Since r = a * sin(2pi / 5)
  //  r_iso = r² / (2 * sqrt(r² - (r / 2 * sin(2pi / 5))²))
  //        = r / (2 * sqrt(1 - 1 / (2 * sin(2pi / 5))²))
  // So we want r_iso > alpha, i.e.
  //   r > (2 * sqrt(1 - 1 / (2 * sin(2pi / 5))²)) * alpha ~= 1.7013016167040798 * alpha
  //
  // Another way is to simply make faces incident to the seed always traversable, and then
  // we only have to ensure faces opposite of the seed are traversable (i.e., radius ~= 1.65 * alpha)
  bool initialize_with_cavities()
  {
#ifdef CGAL_AW3_DEBUG_INITIALIZATION
    std::cout << "> Dig cavities" << std::endl;
    std::cout << m_seeds.size() << " seed(s)" << std::endl;
#endif

    CGAL_precondition(!m_seeds.empty());

    // Get a double value approximating the scaling factors
//    std::cout << sqrt(3) * sin(2pi / 5) << std::endl;
//    std::cout << (2. * std::sqrt(1. - 1. / square(2 * std::sin(2 * CGAL_PI / 5)))) << std::endl;

    Iso_cuboid_3 bbox = SC2GT()(m_bbox);

    std::vector<Vertex_handle> seed_vs;
    for(const Point_3& seed_p : m_seeds)
    {
#ifdef CGAL_AW3_DEBUG_INITIALIZATION
      std::cout << "Initialize from seed " << seed_p << std::endl;
#endif

      if(bbox.has_on_unbounded_side(seed_p))
      {
#ifdef CGAL_AW3_DEBUG_INITIALIZATION
        std::cerr << "Warning: seed " << seed_p << " is outside the bounding box" << std::endl;
#endif
        continue;
      }

      // get the closest point on the input
      const Point_3 closest_pt = m_oracle.closest_point(seed_p);
      const FT sq_d_to_closest = geom_traits().compute_squared_distance_3_object()(seed_p, closest_pt);

      if(sq_d_to_closest < m_sq_offset)
      {
#ifdef CGAL_AW3_DEBUG_INITIALIZATION
        std::cerr << "Warning: seed " << seed_p << " is in the offset" << std::endl;
#endif
        continue;
      }

      // Mark the seeds and icosahedron vertices as "scaffolding" vertices such that the facets
      // incident to these vertices are always traversable regardless of their circumcenter.
      // This is done because otherwise some cavities can appear on the mesh: non-traversable facets
      // with two vertices on the offset, and the third being a deeper inside seed / ico_seed.
      // Making them traversable will cause more refinement than "alpha", but they will eventually
      // not appear anymore in the inside/outside boundary and the surface will look smoother.
      //
      // This problem only appears when the seed and icosahedron vertices are close to the offset surface,
      // which usually happens for large alpha values.

      Vertex_handle seed_v = m_tr.insert(seed_p);
      seed_v->type() = AW3i::Vertex_type::SEED_VERTEX;
      seed_vs.push_back(seed_v);

      // Icosahedron vertices (see also BGL::make_icosahedron())
      const Point_3 center = seed_p;
      const FT radius = 1.65 * m_alpha;
      const FT phi = (FT(1) + approximate_sqrt(FT(5))) / FT(2);
      const FT t = radius / approximate_sqrt(1 + square(phi));
      const FT t_phi = t * phi;

      std::array<Point_3, 12> ico_ps =
      {
        Point_3(center.x(), center.y() + t, center.z() + t_phi),
        Point_3(center.x(), center.y() + t, center.z() - t_phi),
        Point_3(center.x(), center.y() - t, center.z() + t_phi),
        Point_3(center.x(), center.y() - t, center.z() - t_phi),

        Point_3(center.x() + t, center.y() + t_phi, center.z()),
        Point_3(center.x() + t, center.y() - t_phi, center.z()),
        Point_3(center.x() - t, center.y() + t_phi, center.z()),
        Point_3(center.x() - t, center.y() - t_phi, center.z()),

        Point_3(center.x() + t_phi, center.y(), center.z() + t),
        Point_3(center.x() + t_phi, center.y(), center.z() - t),
        Point_3(center.x() - t_phi, center.y(), center.z() + t),
        Point_3(center.x() - t_phi, center.y(), center.z() - t)
      };

      for(const Point_3& seed_neighbor_p : ico_ps)
      {
#ifdef CGAL_AW3_DEBUG_PP
        std::cout << seed_neighbor_p << std::endl;
#endif
        if(bbox.has_on_unbounded_side(seed_neighbor_p))
          continue;

        Vertex_handle ico_v = m_tr.insert(seed_neighbor_p, seed_v /*hint*/);
        ico_v->type() = AW3i::Vertex_type::SEED_VERTEX;
      }
    }

    if(seed_vs.empty())
    {
#ifdef CGAL_AW3_DEBUG_INITIALIZATION
      std::cerr << "Error: no acceptable seed was provided" << std::endl;
#endif
      return false;
    }

#ifdef CGAL_AW3_DEBUG_INITIALIZATION
    std::cout << m_tr.number_of_vertices() - 8 /*bbox*/ << " vertice(s) due to seeds" << std::endl;
#endif

    for(Vertex_handle seed_v : seed_vs)
    {
      std::vector<Cell_handle> inc_cells;
      inc_cells.reserve(64);
      m_tr.incident_cells(seed_v, std::back_inserter(inc_cells));

      for(Cell_handle ch : inc_cells)
        ch->set_label(cavity_cell_label(ch));
    }

    // Should be cheap enough to go through the full triangulation as only seeds have been inserted
    for(Cell_handle ch : m_tr.all_cell_handles())
    {
      if(ch->is_inside())
        continue;

      // When the algorithm starts from a manually dug hole, infinite cells are initialized
      // as "inside" such that they do not appear on the boundary
      CGAL_assertion(!m_tr.is_infinite(ch));

      for(int i=0; i<4; ++i)
        push_facet(std::make_pair(ch, i));
    }

    if(m_queue.empty())
    {
#ifdef CGAL_AW3_DEBUG_INITIALIZATION
      std::cerr << "Could not initialize the algorithm with these seeds, and alpha|offset values" << std::endl;
#endif
      return false;
    }

    return true;
  }

  // tag all infinite cells "outside" and all finite cells "inside"
  // init queue with all convex hull facets
  bool initialize_from_infinity()
  {
    for(Cell_handle ch : m_tr.all_cell_handles())
    {
      if(m_tr.is_infinite(ch))
      {
        ch->set_label(Cell_label::OUTSIDE);
        const int inf_index = ch->index(m_tr.infinite_vertex());
        push_facet(std::make_pair(ch, inf_index));
      }
      else
      {
        ch->set_label(Cell_label::INSIDE);
      }
    }

    return true;
  }

  void reset_manifold_labels()
  {
    // No erase counter increment, or it will mess up with a possibly non-empty queue.
    for(Cell_handle ch : m_tr.all_cell_handles())
      if(ch->label() == Cell_label::MANIFOLD)
        ch->set_label(Cell_label::OUTSIDE);
  }

  // This function is used in the case of resumption of a previous run: m_tr is not cleared,
  // and we fill the queue with the new parameters.
  bool initialize_from_existing_triangulation()
  {
#ifdef CGAL_AW3_DEBUG_INITIALIZATION
    std::cout << "Restart from a DT of " << m_tr.number_of_cells() << " cells" << std::endl;
#endif

    for(Cell_handle ch : m_tr.all_cell_handles())
    {
      if(ch->is_inside())
        continue;

      for(int i=0; i<4; ++i)
      {
        if(ch->neighbor(i)->is_inside())
          push_facet(std::make_pair(ch, i));

      }
    }

    return true;
  }

private:
  // Manifoldness is tolerated while debugging and extracting at intermediate states
  // Not the preferred way because it uses 3*nv storage
  template <typename OutputMesh, typename OVPM>
  void extract_possibly_non_manifold_surface(OutputMesh& output_mesh,
                                             OVPM ovpm) const
  {
    namespace PMP = Polygon_mesh_processing;

#ifdef CGAL_AW3_DEBUG
    std::cout << "> Extract possibly non-manifold wrap... ()" << std::endl;
#endif

    remove_all_elements(output_mesh);

    CGAL_assertion_code(for(auto cit=m_tr.finite_cells_begin(), cend=m_tr.finite_cells_end(); cit!=cend; ++cit))
    CGAL_assertion(cit->tds_data().is_clear());

    for(auto cit=m_tr.finite_cells_begin(), cend=m_tr.finite_cells_end(); cit!=cend; ++cit)
    {
      Cell_handle seed = cit;
      if(seed->is_outside() || seed->tds_data().processed())
        continue;

      std::queue<Cell_handle> to_visit;
      to_visit.push(seed);

      std::vector<Point_3> points;
      std::vector<std::vector<size_t> > faces;
      std::size_t idx = 0;

      while(!to_visit.empty())
      {
        const Cell_handle cell = to_visit.front();
        CGAL_assertion(!cell->is_outside() && !m_tr.is_infinite(cell));

        to_visit.pop();

        if(cell->tds_data().processed())
          continue;
        cell->tds_data().mark_processed();

        for(int fid=0; fid<4; ++fid)
        {
          const Cell_handle neighbor = cell->neighbor(fid);
          if(neighbor->is_outside())
          {
            // There shouldn't be any artificial vertex on the inside/outside boundary
            // (past initialization)
//            CGAL_assertion(cell->vertex((fid + 1)&3)->type() == AW3i::Vertex_type::DEFAULT);
//            CGAL_assertion(cell->vertex((fid + 2)&3)->type() == AW3i::Vertex_type::DEFAULT);
//            CGAL_assertion(cell->vertex((fid + 3)&3)->type() == AW3i::Vertex_type::DEFAULT);

            points.push_back(m_tr.point(cell, Triangulation::vertex_triple_index(fid, 0)));
            points.push_back(m_tr.point(cell, Triangulation::vertex_triple_index(fid, 1)));
            points.push_back(m_tr.point(cell, Triangulation::vertex_triple_index(fid, 2)));
            faces.push_back({idx, idx + 1, idx + 2});
            idx += 3;
          }
          else
          {
            to_visit.push(neighbor);
          }
        }
      }

      CGAL_assertion(PMP::is_polygon_soup_a_polygon_mesh(faces));
      PMP::polygon_soup_to_polygon_mesh(points, faces, output_mesh,
                                        CGAL::parameters::default_values(),
                                        CGAL::parameters::vertex_point_map(ovpm));

      PMP::stitch_borders(output_mesh, CGAL::parameters::vertex_point_map(ovpm));
      CGAL_assertion(is_closed(output_mesh));
    }

    for(auto cit=m_tr.finite_cells_begin(), cend=m_tr.finite_cells_end(); cit!=cend; ++cit)
      cit->tds_data().clear();

    CGAL_postcondition(!is_empty(output_mesh));
    CGAL_postcondition(is_valid_polygon_mesh(output_mesh));
    CGAL_postcondition(is_closed(output_mesh));

    PMP::orient_to_bound_a_volume(output_mesh, CGAL::parameters::vertex_point_map(ovpm));

    collect_garbage(output_mesh);
  }

  template <typename OutputMesh, typename OVPM>
  void extract_manifold_surface(OutputMesh& output_mesh,
                                OVPM ovpm) const
  {
    namespace PMP = Polygon_mesh_processing;

#ifdef CGAL_AW3_DEBUG
    std::cout << "> Extract manifold wrap... ()" << std::endl;
#endif

    CGAL_assertion_code(for(Vertex_handle v : m_tr.finite_vertex_handles()))
    CGAL_assertion(!is_non_manifold(v));

    remove_all_elements(output_mesh);

    std::vector<Point_3> points;
    std::vector<std::array<std::size_t, 3> > faces;

    std::unordered_map<Vertex_handle, std::size_t> vertex_to_id;

    for(auto fit=m_tr.all_facets_begin(), fend=m_tr.all_facets_end(); fit!=fend; ++fit)
    {
      Facet f = *fit;
      if(!f.first->is_outside())
        f = m_tr.mirror_facet(f);

      const Cell_handle ch = f.first;
      const int s = f.second;
      const Cell_handle nh = ch->neighbor(s);
      if(!is_on_outside_boundary(ch, nh))
        continue;

      std::array<std::size_t, 3> ids;
      for(int pos=0; pos<3; ++pos)
      {
        Vertex_handle vh = ch->vertex(Triangulation::vertex_triple_index(s, pos));
        auto insertion_res = vertex_to_id.emplace(vh, vertex_to_id.size());
        if(insertion_res.second) // successful insertion, never-seen-before vertex
          points.push_back(m_tr.point(vh));

        ids[pos] = insertion_res.first->second;
      }

      faces.emplace_back(std::array<std::size_t, 3>{ids[0], ids[1], ids[2]});
    }

#ifdef CGAL_AW3_DEBUG
    std::cout << "\t" << points.size() << " points" << std::endl;
    std::cout << "\t" << faces.size() << " faces" << std::endl;
#endif

    if(faces.empty())
    {
#ifdef CGAL_AW3_DEBUG
      std::cerr << "Warning: empty wrap?..." << std::endl;
#endif
      return;
    }

    if(!PMP::is_polygon_soup_a_polygon_mesh(faces))
    {
      CGAL_warning_msg(false, "Failed to extract a manifold boundary!");
      return;
    }

    PMP::polygon_soup_to_polygon_mesh(points, faces, output_mesh,
                                      CGAL::parameters::default_values(),
                                      CGAL::parameters::vertex_point_map(ovpm));

    CGAL_postcondition(!is_empty(output_mesh));
    CGAL_postcondition(is_valid_polygon_mesh(output_mesh));
    CGAL_postcondition(is_closed(output_mesh));
    CGAL_postcondition(PMP::does_bound_a_volume(output_mesh, CGAL::parameters::vertex_point_map(ovpm)));

    PMP::orient_to_bound_a_volume(output_mesh, CGAL::parameters::vertex_point_map(ovpm));
  }

public:
  template <typename OutputMesh, typename OVPM>
  void extract_surface(OutputMesh& output_mesh,
                       OVPM ovpm,
                       const bool tolerate_non_manifoldness = false) const
  {
    if(tolerate_non_manifoldness)
      extract_possibly_non_manifold_surface(output_mesh, ovpm);
    else
      extract_manifold_surface(output_mesh, ovpm);
  }

private:
  bool is_traversable(const Facet& f) const
  {
    return less_squared_radius_of_min_empty_sphere(m_sq_alpha, f, m_tr);
  }

  enum class Steiner_construction_rule
  {
    NONE = 0,
    R1 = 1,
    R2 = 2
  };

  enum class Steiner_location
  {
    // do not change values without changing the places where they are compared
    ON_CORNER = 0,
    ON_CREASE = 1,
    ON_FACE = 2,
    UNOPTIMIZED = 3,
    MINIMIZER = 4 // sample point that minimizes the QEM error of all sample points
  };

  Steiner_status compute_steiner_point(const Cell_handle ch,
                                       const Cell_handle neighbor,
                                       Point_3& steiner_point) const
  {
    CGAL_precondition(!m_tr.is_infinite(neighbor));

    typename Geom_traits::Construct_ball_3 ball = geom_traits().construct_ball_3_object();
    typename Geom_traits::Construct_vector_3 vector = geom_traits().construct_vector_3_object();
    typename Geom_traits::Construct_translated_point_3 translate = geom_traits().construct_translated_point_3_object();
    typename Geom_traits::Construct_scaled_vector_3 scale = geom_traits().construct_scaled_vector_3_object();

    const Point_3& neighbor_cc = circumcenter(neighbor);
    const Ball_3 neighbor_cc_offset_ball = ball(neighbor_cc, m_sq_offset);
    const bool is_neighbor_cc_in_offset = m_oracle.do_intersect(neighbor_cc_offset_ball);

#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
    std::cout << "Compute_steiner_point(" << &*ch << ", " << &*neighbor << ")" << std::endl;

    const Point_3& chc = circumcenter(ch);
    std::cout << "CH" << std::endl;
    std::cout << "\t" << ch->vertex(0)->point() << std::endl;
    std::cout << "\t" << ch->vertex(1)->point() << std::endl;
    std::cout << "\t" << ch->vertex(2)->point() << std::endl;
    std::cout << "\t" << ch->vertex(3)->point() << std::endl;
    std::cout << "cc: " << chc << std::endl;
    std::cout << "CC Distance to input: " << CGAL::squared_distance(chc, m_oracle.closest_point(chc)) << std::endl;

    std::cout << "NCH" << std::endl;
    std::cout << "\t" << neighbor->vertex(0)->point() << std::endl;
    std::cout << "\t" << neighbor->vertex(1)->point() << std::endl;
    std::cout << "\t" << neighbor->vertex(2)->point() << std::endl;
    std::cout << "\t" << neighbor->vertex(3)->point() << std::endl;
    std::cout << "ncc: " << neighbor_cc << std::endl;
    std::cout << "NC Distance to input: " << CGAL::squared_distance(neighbor_cc, m_oracle.closest_point(neighbor_cc)) << std::endl;
    std::cout << "is_neighbor_cc_in_offset = " << is_neighbor_cc_in_offset << std::endl;

    std::cout << "squared offset " << m_sq_offset << std::endl;
#endif

    // ch's circumcenter should not be within the offset volume
    CGAL_assertion_code(const Point_3& ch_cc = circumcenter(ch);)
    CGAL_assertion_code(const Ball_3 ch_cc_offset_ball = ball(ch_cc, m_sq_offset);)
    CGAL_assertion(!m_oracle.do_intersect(ch_cc_offset_ball));

    if(is_neighbor_cc_in_offset)
    {
      const Point_3& ch_cc = circumcenter(ch);

      // If the voronoi edge intersects the offset, the steiner point is the first intersection
      if(m_oracle.first_intersection(ch_cc, neighbor_cc, steiner_point, m_offset))
      {
#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
        std::cout << "Steiner found through first_intersection(): " << steiner_point << std::endl;
#endif
        return Steiner_status::RULE_1;
      }
    }

    Tetrahedron_with_outside_info<Geom_traits> tet(neighbor, geom_traits());
    if(m_oracle.do_intersect(tet))
    {
      // steiner point is the closest point on input from cell centroid with offset
      const Point_3 closest_pt = m_oracle.closest_point(neighbor_cc);
#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
      std::cout << "Steiner found through neighboring triangle intersecting the input" << std::endl;
      std::cout << "Closest point: " << closest_pt << std::endl;
#endif
      CGAL_assertion(closest_pt != neighbor_cc);

      Vector_3 unit = vector(closest_pt, neighbor_cc);

      // PMP::internal::normalize() requires sqrt()
      unit = scale(unit, m_offset / approximate_sqrt(geom_traits().compute_squared_length_3_object()(unit)));
      steiner_point = translate(closest_pt, unit);

#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
      std::cout << "Direction: " << unit << std::endl;
      std::cout << "Steiner: " << steiner_point << std::endl;
#endif

      return Steiner_status::RULE_2;
    }

#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
    std::cout << "No Steiner point" << std::endl;
#endif

    return Steiner_status::NO_STEINER_POINT;
  }

  Steiner_status compute_steiner_point(const Gate& gate,
                                       Point_3& steiner_point) const
  {
    const Facet& f = gate.facet();
    const Cell_handle ch = f.first;
    const int s = f.second;
    const Cell_handle nh = ch->neighbor(s);
    return compute_steiner_point(ch, nh, steiner_point);
  }

    // given qem, move it so that it is 'offset' away from the closest input point
  bool project_to_offset(const Point_3& p,
                         Point_3& offset_pt,
                         Vector_3& normal) const
  {
    std::cout << "Projecting " << p << std::endl;

    const Point_3 closest_pt = m_oracle.closest_point(p);

    // @fixme this is not correct and can give points arbitrarily close to the input (concave zones)

    // if(squared_distance(closest_pt, p) < m_sq_offset) // @todo give a bit of tolerance?
    // {
    //   std::cerr << "Warning: Projection of point within offset: " << p << std::endl;
    //   return false;
    // }

    normal = p - closest_pt;
    const FT sqnorm = normal.squared_length();
    if(sqnorm != 0.0)
      normal = normal / CGAL::approximate_sqrt(sqnorm);

    offset_pt = closest_pt + m_offset * normal;
    std::cout << " to " << offset_pt << std::endl;

    return true;
  }

  bool project_to_offset(const Point_3& p,
                         Point_3& offset_pt) const
  {
    Vector_3 unused_normal;
    return project_to_offset(p, offset_pt, unused_normal);
  }

  std::pair<Point_3, Steiner_location> qem_error_minimizer(const std::vector<Point_3>& p_proj,
                                                           const std::vector<Vector_3>& n_proj,
                                                           const std::vector<FT>& w,
                                                           const std::vector<Point_3>& candidates) const
  {
    Eigen::MatrixXd Q(4, 4);
    Q.setZero();

    for(int i=0; i<p_proj.size(); ++i)
    {
      const Vector_3& normal = n_proj[i];
      double area = w[i];

      double a = normal.x();
      double b = normal.y();
      double c = normal.z();
      double d = CGAL::scalar_product(normal, (p_proj[i] - CGAL::ORIGIN));

      Eigen::MatrixXd local_Q(4, 4);
      local_Q << area * a * a, area * a * b, area * a * c, area * a * d,
                 area * a * b, area * b * b, area * b * c, area * b * d,
                 area * a * c, area * b * c, area * c * c, area * c * d,
                 area * a * d, area * b * d, area * c * d, area * d * d;

      Q += local_Q;
    }

    // returns the sample point that minimizes the QEM error among all samples
    FT min_error = std::numeric_limits<FT>::max();
    Point_3 p_min = candidates[0];
    for(int i=0; i<candidates.size(); ++i)
    {
      Eigen::VectorXd p(4);
      p << candidates[i].x(), candidates[i].y(), candidates[i].z(), 1.;
      FT error = p.transpose() * Q * p;
      if(error < min_error)
      {
        min_error = error;
        p_min = candidates[i];
      }
    }

    return std::make_pair(p_min, Steiner_location::MINIMIZER);
  }

  std::pair<Point_3, Steiner_location> qem_weighted(const Point_3& p,
                                                    const std::vector<Point_3>& p_proj,
                                                    const std::vector<Vector_3>& n_proj,
                                                    const std::vector<FT>& w) const
  {
    Eigen::MatrixXd A(3, 3);
    A.setZero();
    Eigen::VectorXd B(3);
    B.setZero();

    for(int i=0; i<p_proj.size(); ++i)
    {
      const Vector_3& normal = n_proj[i];
      double area = w[i];

      double a = normal.x();
      double b = normal.y();
      double c = normal.z();
      double d = CGAL::scalar_product(normal, (p_proj[i] - CGAL::ORIGIN));

      Eigen::MatrixXd local_A(3, 3);
      local_A << area * a * a, area * a * b, area * a * c,
                 area * a * b, area * b * b, area * b * c,
                 area * a * c, area * b * c, area * c * c;

      Eigen::VectorXd local_B(3);
      local_B << area * a * d, area * b * d, area * c * d;

      A += local_A;
      B += local_B;
    }

    // analyze the eigenvalues of matrix A to tell us the nature of the QEM point
    const auto evs = A.eigenvalues();
    const auto evsr = evs.real();
    const double ev_min = std::min( { evsr[0], evsr[1], evsr[2] } );
    const double ev_max = std::max( { evsr[0], evsr[1], evsr[2] } );
    const double ev_mid = std::max( { std::min(evsr[0], evsr[1]),
                                      std::min(evsr[0], evsr[2]),
                                      std::min(evsr[1], evsr[2]) } );

    const double min_over_max = ev_min / ev_max;
    const double mid_over_max = ev_mid / ev_max;

    Steiner_location type;
    if(min_over_max > 0.1) // all eigenvalues roughly equal
      type = Steiner_location::ON_CORNER;
    else if(mid_over_max > 0.1) // max&mid are roughly equal
      type = Steiner_location::ON_CREASE;
    else
      type = Steiner_location::ON_FACE;

    std::cout << "EVS " << evsr[0] << " " << evsr[1] << " " << evsr[2] << std::endl;
    std::cout << "type ==> " << static_cast<int>(type) << std::endl;

    // check rank
    Eigen::FullPivLU<Eigen::MatrixXd> solver(A);
    solver.setThreshold(1e-5);

    Eigen::VectorXd x_hat(3);
    x_hat(0) = p.x();
    x_hat(1) = p.y();
    x_hat(2) = p.z();

    Eigen::JacobiSVD<Eigen::MatrixXd> svd_decomp(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    svd_decomp.setThreshold(1e-5);

    const Eigen::VectorXd optim = x_hat + svd_decomp.solve(B - A * x_hat);
    const Point_3& p_new = Point_3(optim(0), optim(1), optim(2));
    return std::make_pair(p_new, type);
  }

  std::pair<Point_3, Steiner_location> qem_weighted(const Point_3& p,
                                                    const std::vector<Triangle_3>& tris,
                                                    const std::vector<FT>& ws) const
  {
    // noet yet implemented
    CGAL_assertion(false);
    return std::make_pair(Point_3(), Steiner_location::UNOPTIMIZED);
  }

  bool compute_steiner_point_QEM(const Facet f,
                                 Point_3& steiner_point,
                                 Steiner_location& steiner_location) const
  {
#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
    std::cout << " ======== computing Steiner point w/ QEM" << std::endl;
#endif

    // get steiner point first
    Steiner_construction_rule steiner_rule;
    bool successful_computation = compute_steiner_point(f, steiner_point, steiner_rule);

    if(!successful_computation)
      return false;

#ifdef CGAL_AW3_ORIGINAL_STEINER_CONSTRUCTION
    return true;
#endif

#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
    std::cout << "Original Steiner point is: " << steiner_point << " found by R" << static_cast<int>(steiner_rule) << "\n";
#endif

    Point_3 optimized_steiner_point;
#ifdef CGAL_AW3_SHARPEN_WITH_SAMPLING
    Sampling_optimization_status res =
      optimize_with_sampling(f, steiner_point, steiner_rule,
                             optimized_steiner_point, steiner_location);
    bool successful_optimization = (res == Sampling_optimization_status::FOUND_OPTIMAL_POINT);
#elif defined(CGAL_AW3_SHARPEN_WITH_ITERATIVE_SAMPLING)
    bool successful_optimization =
      optimize_with_iterative_sampling(f, steiner_point, steiner_rule,
                                       optimized_steiner_point, steiner_location);
#elif defined(CGAL_AW3_SHARPEN_WITH_ITERATIVE_BALL_QUERIES)
    bool successful_optimization =
      optimize_with_iterative_ball_queries(f, steiner_point, 2 /*iter*/,
                                           m_offset + 1.5 * m_alpha /*radius*/, // @fixme silly for offset >>> alpha
                                           optimized_steiner_point, steiner_location);
#elif defined(CGAL_AW3_SHARPEN_WITH_DYNAMIC_BALL_QUERIES)
      optimize_with_dynamic_ball_queries(f, steiner_point, optimized_steiner_point, steiner_location);
#endif

#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
    std::cout << "\tsuccessful optimization = " << successful_optimization << "\n";
#endif

    if(successful_optimization)
      steiner_point = optimized_steiner_point;

    return true;
  }

  bool seek_closer_existing_vertex(const Point_3& steiner_point,
                                   const FT sq_radius) const
  {
    Vertex_handle vh = m_tr.nearest_vertex(steiner_point);
    if(CGAL::squared_distance(steiner_point, m_tr.point(vh)) < sq_radius)
    {
      std::cout << "Closest point is too close to: " << m_tr.point(vh) << std::endl;
      return true;
    }

    return false;
  }

  bool optimize_with_iterative_ball_queries(const Facet& f,
                                            const Point_3& steiner_point,
                                            const int num_iter,
                                            const FT radius,
                                            Point_3& optimized_steiner_point,
                                            Steiner_location& steiner_location) const
  {
    const Cell_handle ch = f.first;
    const int id = f.second;
    const Cell_handle neighbor = ch->neighbor(id);

    Point_3 intermediate_steiner_point = steiner_point;

    for(int i=0; i<num_iter; ++i)
    {
      // snap_point_by_ball snap_point_by_ball_all_combinations ctrl+C helper
      if(snap_point_by_ball_all_combinations(f, intermediate_steiner_point, radius, optimized_steiner_point, steiner_location))
      {
        std::cout << "snap by ball point, iteration #" << i << " --> new point: " << optimized_steiner_point << std::endl;
        intermediate_steiner_point = optimized_steiner_point;
      }
      else
      {
        std::cout << "snap by ball, iteration #" << i << " --> no point" << std::endl;
        break;
      }
    }

    return (steiner_location != Steiner_location::UNOPTIMIZED);
  }

  bool check_does_break_gate(const Facet& f,
                             const Point_3& query_point) const
  {
    const Cell_handle ch = f.first;
    const int id = f.second;
    const Cell_handle neighbor = ch->neighbor(id);

    typename Geom_traits::Construct_sphere_3 sphere = geom_traits().construct_sphere_3_object();
    bool breaks_gate = false;

    if((!m_tr.is_infinite(ch) &&
         m_tr.side_of_sphere(ch, query_point, false /*perturb*/) == ON_BOUNDED_SIDE) ||
       (!m_tr.is_infinite(neighbor) &&
         m_tr.side_of_sphere(neighbor, query_point, false /*perturb*/) == ON_BOUNDED_SIDE))
    {
      breaks_gate = true;
    }

    // logging result
    std::cout << "\tquery point " << query_point
              << " does" << (breaks_gate ? " " : " NOT ") << "break the gate" << std::endl;

    return breaks_gate;
  }

  std::vector<Point_3> get_neighbor_samples(const Point_3& og_point) const
  {
    // find the closest primitive of the og point
    const auto closest_pair = m_oracle.tree().closest_point_and_primitive(og_point);
    auto dpm = m_oracle.get_m_dpm();
    const Triangle_3 p_tri = get(dpm, closest_pair.second);
    Vector_3 axis = p_tri.supporting_plane().orthogonal_vector();

    // find four vectors orthogonal to the normal, 90 degrees apart from each other
    const FT mag = axis.squared_length();
    if(mag != 0.0)
      axis /= CGAL::sqrt(CGAL::to_double(mag));

    const FT u_x = axis.x();
    const FT u_y = axis.y();
    const FT u_z = axis.z();
    // vs SC = simple_Cartesion<double>? RT vs FT for args? will have type errors
    CGAL::Aff_transformation_3<Geom_traits> rotate(
      u_x*u_x      , u_x*u_y - u_z, u_x*u_z + u_y,
      u_y*u_x + u_z, u_y*u_y      , u_y*u_z - u_x,
      u_z*u_x - u_y, u_z*u_y + u_x, u_z*u_z
    );

    std::vector<Vector_3> samples;
    Vector_3 sample = Vector_3(0, 1, 0);
    if(u_y != 0.0)
      sample = Vector_3(1, CGAL::to_double(u_x)/CGAL::to_double(u_y), 0); // perpendicular to axis

    samples.push_back(sample);
    std::cout << sample << std::endl;
    for(int i=0; i<3; ++i)
    {
      sample = rotate.transform(sample);
      samples.push_back(sample);
      std::cout << sample << std::endl;
    }

    // resize vector magnitudes to be m_alpha
    const double circumradius = CGAL::sqrt(geom_traits().compute_squared_radius_3_object()(p_tri.vertex(0), p_tri.vertex(1), p_tri.vertex(2)));
    for(auto& sample : samples)
    {
      const FT sample_mag = sample.squared_length();
      if(sample_mag != 0.0)
      {
        sample /= CGAL::sqrt(CGAL::to_double(sample_mag));
        sample *= m_alpha;
        // sample *= circumradius/5;
      }
    }

    // add these vectors to original point to get new sample points
    std::vector<Point_3> points;
    for(auto& sample : samples)
    {
      const Point_3 new_point = og_point + sample;
      points.push_back(new_point);
      std::cout << new_point << std::endl;
    }

    return points;
  }

  // given a point on offset surface, probe its neighbors and refind another qem point
  // returns true if a new point was found, false otherwise
  bool snap_point_by_neighbors(const Point_3& initial_point,
                               Point_3& optimized_point,
                               Steiner_location& optimized_point_location) const
  {
    // you need to port a lot of improvements from the other function first (like seek closer,
    // QEM type 2 ignoring etc.)
    CGAL_assertion(false);

    // find the samples
    const std::vector<Point_3> samples = get_neighbor_samples(initial_point);

    std::vector<Point_3> p_proj;
    std::vector<Vector_3> n_proj;
    std::vector<FT> w;
    std::vector<std::pair<int, int> > primitives;

    auto dpm = m_oracle.get_m_dpm();

    // for each sample, find their closest primitive
    for(int i=0; i<samples.size(); ++i)
    {
      const Point_3 sample = samples[i];
      const auto closest_pair = m_oracle.tree().closest_point_and_primitive(sample);
      const Point_3 closest_pt = closest_pair.first;
      std::cout << "\tNew snap primitive id: " << closest_pair.second.first << ", " << closest_pair.second.second << std::endl;

#ifdef CGAL_AW3_USE_NORMALS_OF_CLOSEST_PRIMITIVE
      const Triangle_3& p_tri = get(dpm, closest_pair.second);
      Vector_3 normal = p_tri.supporting_plane().orthogonal_vector(); // should always point outward?
      const FT area = CGAL::sqrt(p_tri.squared_area());
#else
      Vector_3 normal = sample - closest_pt;
      const FT area = 1;
#endif

      FT mag_norm = normal.squared_length();
      if(mag_norm != 0.0)
        normal = normal / CGAL::approximate_sqrt(mag_norm);

      p_proj.push_back(closest_pt + normal*m_offset);
      n_proj.push_back(normal);
      w.push_back(area);
      primitives.push_back(closest_pair.second);
    }

    if(!p_proj.empty())
    {
      Point_3 qem_point;
      std::tie(qem_point, optimized_point_location) = qem_weighted(initial_point, p_proj, n_proj, w);

      bool successful_projection = project_to_offset(qem_point, optimized_point);
      if(!successful_projection)
        return false;

      return (initial_point != optimized_point);
    }

    return false;
  }

  // very similar to snap_point_by_neighbors, but acquiring the nearby primitives in a different way
  // for now, assumes ball will intersect input when center is on offset (by choice of radius)
  // initial_point could be the original steiner point, or the snapped steiner, or the qem
  bool snap_point_by_ball(const Facet& f,
                          const Point_3& initial_point,
                          const FT radius,
                          Point_3& optimized_point,
                          Steiner_location& optimized_point_location) const
  {
    std::cout << "  Snap point " << initial_point << std::endl;

    // you need to port a lot of improvements from the other function first (like seek closer,
    // QEM type 2 ignoring etc.)
    CGAL_assertion(false);

    // construct the ball
    typename Geom_traits::Construct_ball_3 ball = geom_traits().construct_ball_3_object();

    const Ball_3 neighborhood_ball = ball(initial_point, square(radius));

    std::cout << "query ball: " << initial_point << " radius " << radius << std::endl;

    // see where the ball intersects the input
    std::list<std::pair<int, int> > primitives;
    m_oracle.tree().all_intersected_primitives(neighborhood_ball, std::back_inserter(primitives));

    std::cout << primitives.size() << " intersected primitives" << std::endl;

    // given those primitives, compute some qem points (for now, using all of them)
    // determine which qem point to return, if any
    std::vector<Point_3> p_proj;
    std::vector<Vector_3> n_proj;
    std::vector<FT> w;
    std::vector<Plane_3> planes; // @todo make sure it isn't a single plane?

    const auto dpm = m_oracle.get_m_dpm();
    for(auto primitive : primitives)
    {
      const Triangle_3& p_tri = get(dpm, primitive);

      std::cout << "\t\tintersected primitive: " << primitive.first << ", " << primitive.second << std::endl;
      std::cout << "\t\tdatum: " << p_tri << std::endl;

      const Plane_3 plane = p_tri.supporting_plane();
      Vector_3 normal = plane.orthogonal_vector();
      const FT area = 1; // CGAL::approximate_sqrt(p_tri.squared_area()); // @tmp
      FT mag_norm = normal.squared_length();
      if(mag_norm != 0.)
        normal /= CGAL::approximate_sqrt(mag_norm);

      p_proj.push_back(p_tri.vertex(0) + normal * m_offset);
      n_proj.push_back(normal);
      w.push_back(area);
      planes.push_back(plane);
    }

    if(p_proj.empty())
      return false;

    Point_3 qem_point;
    Steiner_location qem_location;
    std::tie(qem_point, qem_location) = qem_weighted(initial_point, p_proj, n_proj, w);

    std::cout << "QEM point: " << qem_point << " type " << static_cast<int>(qem_location) << std::endl;

    Point_3 projected_qem_point;
    bool successful_projection = project_to_offset(qem_point, projected_qem_point);
    if(!successful_projection)
      return false;

    if(check_does_break_gate(f, projected_qem_point))
    {
      optimized_point = projected_qem_point;
      std::cout << optimized_point << " breaks the gate" << std::endl;
    }
    else
    {
      std::cout << "optimized steiner point did not break gate" << std::endl;
    }

    return (initial_point != optimized_point);
  }

  // very similar to snap_point_by_ball(), but at a given iteration, keep the best QEM point
  // from all QEM points obtained by combinaisons of 2 or 3 primitives.
  // initial_point could be the original steiner point, or the snapped steiner, or the qem point
  bool snap_point_by_ball_all_combinations(const Facet& f,
                                           const Point_3& initial_point,
                                           const FT radius,
                                           Point_3& optimized_point,
                                           Steiner_location& optimized_point_location)  const
  {
    std::cout << "  Snap point " << initial_point << std::endl;

    // construct the ball
    typename Geom_traits::Construct_ball_3 ball = geom_traits().construct_ball_3_object();

    const Ball_3 neighborhood_ball = ball(initial_point, square(radius));

    std::cout << "query ball: " << initial_point << " radius " << radius << std::endl;

    // see where the ball intersects the input
    using PID = std::pair<int, int>;
    std::vector<PID> primitives;
    m_oracle.tree().all_intersected_primitives(neighborhood_ball, std::back_inserter(primitives));

    std::cout << primitives.size() << " intersected primitives" << std::endl;

    // get all combinations of 2 or 3 primitives within the vector
    // @todo expensive, and not sufficient for complex corners
    std::vector<std::vector<PID> > primitive_combinations;
    for(int i=0; i<primitives.size(); ++i)
      for(int j=i+1; j<primitives.size(); ++j)
        primitive_combinations.push_back({primitives[i], primitives[j]});

    for(int i=0; i<primitives.size(); ++i)
      for(int j=i+1; j<primitives.size(); ++j)
        for(int k=j+1; k<primitives.size(); ++k)
          primitive_combinations.push_back({primitives[i], primitives[j], primitives[k]});

    std::cout << primitive_combinations.size() << " combinations" << std::endl;

    primitive_combinations.push_back(primitives);

    const auto dpm = m_oracle.get_m_dpm();

    // now for all combinations, compute some QEM points
    std::vector<std::pair<Point_3, Steiner_location> > candidate_QEM_points;
    for(const std::vector<PID>& combination : primitive_combinations)
    {
      // given those primitives, compute some qem points (for now, using all of them)
      std::vector<Point_3> p_proj;
      std::vector<Vector_3> n_proj;
      std::vector<FT> w;
      std::vector<Plane_3> planes; // @todo make sure it isn't a single plane?

      // @todo how to select meaningful combinations? E.g. discard a set of parallel planes
      std::cout << "-- new combination" << std::endl;

      for(const auto& primitive : combination)
      {
        const Triangle_3& p_tri = get(dpm, primitive);

        std::cout << "\t\tintersected primitive: " << primitive.first << ", " << primitive.second << std::endl;
        std::cout << "\t\tdatum: " << p_tri << std::endl;

        const Plane_3 plane = p_tri.supporting_plane();
        Vector_3 normal = plane.orthogonal_vector();
        const FT area = 1; // CGAL::approximate_sqrt(p_tri.squared_area()); // @tmp
        FT mag_norm = normal.squared_length();
        if(mag_norm != 0.)
          normal /= CGAL::approximate_sqrt(mag_norm);

        p_proj.push_back(p_tri.vertex(0) + normal * m_offset);
        n_proj.push_back(normal);
        w.push_back(area);
        planes.push_back(plane);
      }

      if(!p_proj.empty())
      {
        Point_3 qem_point;
        Steiner_location qem_location;
        std::tie(qem_point, qem_location) = qem_weighted(initial_point, p_proj, n_proj, w);

        std::cout << "QEM point: " << qem_point << " type " << static_cast<int>(qem_location) << std::endl;

        // Ignore QEM giving us an optimal point on a face: we only care to snap to something sharp
        if(qem_location == Steiner_location::ON_FACE)
          continue;

        // Reject if the optimal point is too far from the planes
        bool reject_candidate = false;
        for(int i=0; i<planes.size(); ++i)
        {
          const Plane_3& plane = planes[i];

          // Try to avoid points being on the wrong offset surface
          // if(plane.oriented_side(qem_point) != CGAL::ON_POSITIVE_SIDE) // @todo too harsh?
          // {
          //   reject_candidate = true;
          //   break;
          // }

          // The idea is that at a convex corner, the intersection of the offset planes is
          // at distance (sqrt(2) - 1) * offset ~= 0.5 * offset of the offset surface,
          // and it shouldn't be rejected.
          const FT sq_dist = squared_distance(qem_point, plane);
          if(sq_dist > 2.25 * m_sq_offset) // @tmp 2.25; 0.25 ?
          {
            reject_candidate = true;
            break;
          }
        }

        if(reject_candidate)
          continue;

        Point_3 projected_qem_point;
        bool successful_projection = project_to_offset(qem_point, projected_qem_point);
        if(!successful_projection)
          return false;

        const FT sq_near_radius = square(0.5 * m_alpha);
        if(check_does_break_gate(f, projected_qem_point) &&
           !seek_closer_existing_vertex(projected_qem_point, sq_near_radius))
        {
          candidate_QEM_points.emplace_back(projected_qem_point, qem_location);
        }
      }
    }

    std::cout << candidate_QEM_points.size() << " QEM candidates" << std::endl;
    for(const auto& cp : candidate_QEM_points)
      std::cout << "\t" << cp.first << " " << static_cast<int>(cp.second) << std::endl;

    if(candidate_QEM_points.empty())
      return false;

    for(const auto& candidate : candidate_QEM_points)
    {
      // This is just saying: "we should prioritize QEM types that have the lowest dimension:
      // corners ideally, or creases otherwise"
      if(candidate.second != optimized_point_location)
      {
        if(candidate.second < optimized_point_location)
        {
          optimized_point = candidate.first;
          optimized_point_location = candidate.second;
        }

        continue;
      }

      // Out of the QEM points that break the gate, get the best QEM point available:
      // for now arbitraly defined as the closest from the CC of the outside cell
      const Point_3& cc = circumcenter(f.first);
      if(squared_distance(cc, candidate.first) < squared_distance(cc, optimized_point))
      {
        optimized_point = candidate.first;
        optimized_point_location = candidate.second;
      }
    }

    std::cout << "Best QEM candidate: " << optimized_point << " type " << static_cast<int>(optimized_point_location) << std::endl;

    return (initial_point != optimized_point);
  }

  std::pair<bool, int> optimize_with_dynamic_ball_queries(const Facet& f,
                                                          const Point_3& initial_point,
                                                          const Steiner_construction_rule steiner_rule,
                                                          Point_3& optimized_point) const
  {
    std::cout << "Dynamic snapping, og opoint is: " << initial_point << std::endl;

    // you need to port a lot of improvements from the other function first (like seek closer,
    // QEM type 2 ignoring etc.)
    CGAL_assertion(false);

    const Cell_handle ch = f.first;
    const int id = f.second;
    const Cell_handle neighbor = ch->neighbor(id);

    const Point_3& neighbor_a = neighbor->vertex(0)->point();
    const Point_3& neighbor_b = neighbor->vertex(1)->point();
    const Point_3& neighbor_c = neighbor->vertex(2)->point();
    const Point_3& neighbor_d = neighbor->vertex(3)->point();

    const Sphere_3 delaunay_ball = Sphere_3(neighbor_a, neighbor_b, neighbor_c, neighbor_d);
    const double max_radius = CGAL::approximate_sqrt(CGAL::to_double(delaunay_ball.squared_radius()));

    FT radius = 0.25 * m_alpha;
    int degen = 3;

    const auto dpm = m_oracle.get_m_dpm();

    Point_3 intermed_point;
    std::set<int> probed_primitives;
    std::set<std::pair<int, int> > new_primitives;
    std::set<std::pair<int, int> > old_primitives;

    while(degen != 0 && radius < max_radius)
    {
      // same as snap_point_by_ball, but use the degen value returned from qem_weighted
      // construct the ball -- center point is projection of initial_point onto the input surface
      typename Geom_traits::Construct_ball_3 ball = geom_traits().construct_ball_3_object();
      const Ball_3 neighborhood_ball = ball(m_oracle.closest_point(initial_point), square(radius));

      // see where the ball intersects the input
      std::list<std::pair<int, int> > primitives;
      m_oracle.tree().all_intersected_primitives(neighborhood_ball, std::back_inserter(primitives));

      // given those primitives, compute some qem points
      // determine which qem point to return, if any
      std::vector<Point_3> p_proj;
      std::vector<Vector_3> n_proj;
      std::vector<FT> w;

      bool found_new_primitives = false;
      for(auto primitive : primitives)
      {
        if(probed_primitives.find(primitive.second) == probed_primitives.end())
        {
          found_new_primitives = true;
          probed_primitives.emplace(primitive.second);
          new_primitives.emplace(primitive);
          std::cout << "\t\tthe new ball intersected primitives are: " << primitive.second << std::endl;
        }

        const Triangle_3& p_tri = get(dpm, primitive);
        Vector_3 normal = p_tri.supporting_plane().orthogonal_vector();
        const FT area = CGAL::sqrt(p_tri.squared_area());
        FT mag_norm = normal.squared_length();
        if(mag_norm != 0.0)
          normal /= CGAL::approximate_sqrt(mag_norm);

        p_proj.push_back(p_tri.vertex(0) + normal*m_offset);
        n_proj.push_back(normal);
        w.push_back(1); // used to be area
      }

      if(!p_proj.empty() && found_new_primitives)
      {
        const auto intermed_point_pair = qem_weighted(initial_point, p_proj, n_proj, w);
        std::cout << "in loop of dynamic, after qem weighted, og point is: " << initial_point << std::endl;

        bool successful_projection = project_to_offset(intermed_point_pair.first, intermed_point);
        if(!successful_projection)
          return std::make_pair(false, 3);

        std::cout << "\tnewly projected point is: " << intermed_point << std::endl;
        std::cout << "\tsq distance from new point to input is: " << CGAL::approximate_sqrt(m_oracle.squared_distance(intermed_point)) << std::endl;
        std::cout << "\thas degen: " << intermed_point_pair.second << std::endl;
        std::cout << "\twith radius: " << radius << " vs radius: " << max_radius << std::endl;

        if(intermed_point_pair.second == 0)
        {
          const bool found_true_corner = find_best_intermed_corner_point(old_primitives, new_primitives, intermed_point, optimized_point, initial_point, f, steiner_rule);
          if(found_true_corner)
            degen = 0;

          std::cout << "the final point used is: " << optimized_point
                  << " while the og point is: " << initial_point
                  << " and are they equal? " << (optimized_point == initial_point) << std::endl;

          return std::make_pair(optimized_point != initial_point, degen);
        }
        else if(intermed_point_pair.second < degen &&
                check_does_break_gate(f, intermed_point))
        {
          std::cout << "decreased degen to: " << intermed_point_pair.second << " and updated point (which breaks gate) to: " << intermed_point << std::endl;
          std::cout << "meanwhile og point is: " << initial_point << std::endl;
          degen = intermed_point_pair.second;
          optimized_point = intermed_point;
        }
        // else if((intermed_point_pair.second == degen) && (optimized_point != intermed_point))
        // {
        //   std::cout << "changed point but did not decrease rank -- changed to bad point, end. snap point: " << optimized_point << " vs new intermed point: " << intermed_point << std::endl;
        //   return std::make_pair(optimized_point != initial_point, degen);
        // }
      }
      // updates for next while loop
      old_primitives.insert(new_primitives.begin(), new_primitives.end());
      new_primitives.clear();
      radius += 0.25 * m_alpha;
    }

    if(degen == 3)
    {
      std::cout << "no radius was successful rip, don't expect this. degen must be < 3" << std::endl;

      return std::make_pair(false, 3);
    }

    std::cout << "did not reach corner but made improvement (interior or edge, expect this), with degen: " << std::endl;

    return std::make_pair(optimized_point != initial_point, degen);
  }

  bool find_best_intermed_corner_point(const std::set<std::pair<int, int> >& old_primitives,
                                       const std::set<std::pair<int, int> >& new_primitives,
                                       const Point_3& all_point,
                                       const Point_3& initial_point,
                                       const Facet& f,
                                       const Steiner_construction_rule steiner_rule,
                                       Point_3& optimized_point) const
  {
    const Cell_handle ch = f.first;
    const int id = f.second;
    const Cell_handle neighbor = ch->neighbor(id);

    std::cout << "in find best intermed corner point, new_primitives has size: " << new_primitives.size() << std::endl;
    std::cout << "old primitives has size: " << old_primitives.size() << std::endl;
    const auto dpm = m_oracle.get_m_dpm();
    std::set<Point_3> corner_points;

    const std::vector<std::pair<int, int> > new_prims_vector(new_primitives.begin(), new_primitives.end());
    const auto all_new_subsets = get_all_subsets(new_prims_vector);
    for(auto new_prim_subset : all_new_subsets)
    {
      std::vector<Point_3> p_proj;
      std::vector<Vector_3> n_proj;
      std::vector<FT> w;

      std::set<std::pair<int, int> > primitives = old_primitives;
      std::copy(new_prim_subset.begin(), new_prim_subset.end(), std::inserter(primitives, primitives.end()));

      for(auto prim : primitives)
      {
        const Triangle_3 p_tri = get(dpm, prim);
        Vector_3 normal = p_tri.supporting_plane().orthogonal_vector();
        const FT area = CGAL::sqrt(p_tri.squared_area());
        FT mag_norm = normal.squared_length();
        if(mag_norm != 0.0)
          normal /= CGAL::approximate_sqrt(mag_norm);

        p_proj.push_back(p_tri.vertex(0) + normal*m_offset);
        n_proj.push_back(normal);
        w.push_back(1);
      }

      if(!p_proj.empty())
      {
        const auto point_pair = qem_weighted(all_point, p_proj, n_proj, w);

        Point_3 corner_candidate;
        bool successful_projection = project_to_offset(point_pair.first, corner_candidate);

        if(successful_projection &&
           point_pair.second == 0 &&
           is_close_enough_to_offset_surface(corner_candidate) &&
           check_does_break_gate(f, corner_candidate))
        {
          std::cout << corner_candidate << " broke the gate related to og steiner point: " << initial_point << std::endl;
          corner_points.emplace(corner_candidate);
        }
      }
    }

    std::cout << "\tfound all the candidates, corner points has size: " << corner_points.size() << std::endl;
    if(corner_points.size() == 0)
    {
      std::cout << "\tin find corner point, no intermed (including all_point)\n";
      return false;
    }

    std::cout << "\tin find corner point, there are: " << corner_points.size() << " potential corners\n";

    // given all potential points, return the "best" according to faithfulness
    std::map<double, Point_3> evaluate_map; // could have same double values though
    for(auto point : corner_points)
    {
      // from gate priority code
      int li, lj = 0;
      Locate_type lt;
      const Cell_handle conflict_cell = m_tr.locate(point, lt, li, lj, neighbor);
      // CGAL_assertion(lt != Triangulation::VERTEX);
      if(lt == Triangulation::VERTEX)
        continue;

      std::vector<Facet> boundary_facets;
      std::vector<Cell_handle> conflict_zone;
      boundary_facets.reserve(32);
      conflict_zone.reserve(32);
      m_tr.find_conflicts(point, conflict_cell,
                          std::back_inserter(boundary_facets),
                          std::back_inserter(conflict_zone));

      double max_faithful_facet = 0.;
      for(const Facet& facet : boundary_facets)
      {
        const Cell_handle ch_facet = facet.first;
        const int id_facet = facet.second;
        for(int i=0; i<3; ++i)
        {
          const Point_3& first = m_tr.point(ch_facet, (id_facet + 1 + i)&3);
          const Point_3& second = m_tr.point(ch_facet, (id_facet + 1 + (i+1)%3)&3);
          CGAL_assertion(first != m_tr.point(ch_facet, id_facet));
          CGAL_assertion(second != m_tr.point(ch_facet, id_facet));

          const Triangle_3 triangle_facet(point, first , second);
          std::queue<Triangle_3> triangles_queue;
          triangles_queue.push(triangle_facet);
          // cout_triangle(triangle_facet);
          // std::cout << "about to start faithful queue\n";
          if(is_faithful_triangles_queue(triangles_queue, m_parsimony_tolerance))
          {
            const double potential_gate_area = triangle_facet.squared_area();
            max_faithful_facet += potential_gate_area;
          }
        }
      }

      std::cout << "\t\tpoint: " << point << " gives faithfulness: " << max_faithful_facet << std::endl;
      evaluate_map.insert({max_faithful_facet, point});
    }

    if(evaluate_map.size() == 0)
    {
      if(is_close_enough_to_offset_surface(all_point))
      {
        optimized_point = all_point;
        std::cout << "\tin find corner point, no new intermed corners but og is close enough\n";
        return true;
      }

      std::cout << "\tin find corner point, no new intermed and og not close enough\n";
      return false;
    }

    const double max_value = evaluate_map.rbegin()->first;
    optimized_point = evaluate_map.find(max_value)->second;
    std::cout << "\tfinal 'true' corner point is: " << optimized_point << " with faithful value: " << max_value << std::endl;

    return true;
  }

  std::vector<std::vector<std::pair<int, int> > > get_all_subsets(const std::vector<std::pair<int, int> > elements) const
  {
    std::vector<std::vector<std::pair<int, int> > > all_subsets;
    std::vector<std::pair<int, int> > current_subset;
    recurse_all_subsets(elements, current_subset, 0, all_subsets);
    return all_subsets;
  }

  void recurse_all_subsets(const std::vector<std::pair<int, int> >& elements,
                           std::vector<std::pair<int, int> >& current_subset,
                           const int index,
                           std::vector<std::vector<std::pair<int, int> > >& all_subsets) const
  {
    all_subsets.push_back(current_subset);

    for(int i=index; i<elements.size(); ++i)
    {
      current_subset.push_back(elements[i]);
      recurse_all_subsets(elements, current_subset, i+1, all_subsets);
      current_subset.pop_back();
    }
  }

  bool is_close_enough_to_offset_surface(const Point_3& q_point) const
  {
    const FT sq_dist = m_oracle.squared_distance(q_point);
    return (sq_dist < 3 * m_sq_offset);
  }

  // Explore the neighborhood using random walks around the original Steiner point
  template <std::size_t n>
  void get_random_walk_samples(const Facet& f,
                               const Point_3& steiner_point,
                               const Steiner_construction_rule steiner_rule,
                               const FT seeking_radius,
                               std::array<Point_3,n>& samples) const
  {
    const Cell_handle ch = f.first;
    const int id = f.second;
    const Point_3& cc = circumcenter(ch); // @fixme handle inf cells
    const FT ball_sq_radius = squared_distance(cc, m_tr.point(ch, 0));
    const FT scaled_ball_sq_radius = FT(1.25) * ball_sq_radius;

    const Cell_handle nc = f.first->neighbor(id);
    const Point_3& ncc = circumcenter(nc);
    const FT neighbor_ball_sq_radius = squared_distance(ncc, m_tr.point(nc, 0));
    const FT scaled_neighbor_ball_sq_radius = FT(1.25) * neighbor_ball_sq_radius;

    const FT precision = 0.01 * m_offset;
    const FT sq_offset_minus_precision = CGAL::square(m_offset - precision);
    const FT sq_offset_plus_precision = CGAL::square(m_offset + precision);

    const FT seeking_sq_radius = CGAL::square(seeking_radius);

    Vector_3 dir; // points away from steiner_point

    if(steiner_rule == Steiner_construction_rule::R1)
    {
      // not dir = { steiner_point, cc } because ch might be infinite, in which case
      // cc is (awkwardly) set to the center of the diametral ball of the facet,
      // and if we use this, the dual might be inverted.
      const Point_3& a = m_tr.point(ch, Triangulation::vertex_triple_index(id, 0));
      const Point_3& b = m_tr.point(ch, Triangulation::vertex_triple_index(id, 1));
      const Point_3& c = m_tr.point(ch, Triangulation::vertex_triple_index(id, 2));
      const Vector_3 edge_1 {a, b};
      const Vector_3 edge_2 {a, c};
      dir = CGAL::cross_product(edge_1, edge_2);
    }
    else
    {
      dir = { steiner_point, ncc };
    }

    CGAL_assertion(dir != CGAL::NULL_VECTOR);
    dir = dir / CGAL::approximate_sqrt(dir.squared_length());

    const Point_3 source = steiner_point + 0.5 * seeking_radius * dir;

    std::cout << "c:\n"
              << "\t" << ch->vertex(0)->point() << "\n"
              << "\t" << ch->vertex(1)->point() << "\n"
              << "\t" << ch->vertex(2)->point() << "\n"
              << "\t" << ch->vertex(3)->point() << std::endl;
    std::cout << "cc: " << cc << std::endl;
    std::cout << "nc:\n"
              << "\t" << nc->vertex(0)->point() << "\n"
              << "\t" << nc->vertex(1)->point() << "\n"
              << "\t" << nc->vertex(2)->point() << "\n"
              << "\t" << nc->vertex(3)->point() << std::endl;
    std::cout << "ncc: " << ncc << std::endl;
    std::cout << "ball_sq_radius: " << ball_sq_radius << std::endl;
    std::cout << "neighbor_ball_sq_radius: " << neighbor_ball_sq_radius << std::endl;
    std::cout << "scaled_ball_sq_radius: " << scaled_ball_sq_radius << std::endl;
    std::cout << "scaled_neighbor_ball_sq_radius: " << scaled_neighbor_ball_sq_radius << std::endl;
    std::cout << "seeking_radius: " << seeking_radius << std::endl;
    std::cout << "steiner point: " << steiner_point << std::endl;
    std::cout << "source: " << source << std::endl;

    CGAL::Random rng(0);
    Random_points_on_sphere_3<Point_3> random_point_on_sphere(1., rng);

    for(int i=0; i<n; ++i)
    {
      std::cout << "Generate path of sample #" << i << std::endl;

      // Start from a random point on the initial sphere
      Point_3 current_pt = source;
      Point_3 closest_point = m_oracle.closest_point(current_pt);

      static int walk_id = 0;
      std::cout << "Walk ID " << m_tr.number_of_vertices() << " " << i << std::endl;
      std::ofstream walk_out("results/walk_" + std::to_string(m_tr.number_of_vertices()) + "_" + std::to_string(walk_id++) + ".cgal");
      walk_out.precision(17);

      for(;;)
      {
        // std::cout << "Current point: " << current_pt << std::endl;
        // std::cout << "Current closest point: " << closest_point << std::endl;
        FT sq_current_dist = squared_distance(current_pt, closest_point);

#ifdef CGAL_AW3_RANDOM_WALK_USE_WALK_ON_SPHERES
        // if the current point is close enough to the offset, we are done
        if((sq_current_dist > sq_offset_minus_precision) &&
           (sq_current_dist < sq_offset_plus_precision))
        {
          samples[i] = current_pt;
          break;
        }

        FT radius = CGAL::approximate_sqrt(sq_current_dist) - m_offset;
#elif defined(CGAL_AW3_RANDOM_WALK_USE_FIXED_STEP)
        // check if there is an intersection between the previous point and the current point
        #error // not yet implemented

        const FT min_step = 0.1 * m_alpha;
#else
# error
#endif

        for(;;) // until a good next point is found @todo safety step count?
        {
#ifdef CGAL_AW3_RANDOM_WALK_USE_WALK_ON_SPHERES
          // take a random point on the boundary of the sphere
          // of radius dist(current_pt, closest_pt) - offset
          Vector_3 rnd_dir { CGAL::ORIGIN, *random_point_on_sphere++ };
          Point_3 candidate_pt = current_pt + radius * rnd_dir;
#elif defined(CGAL_AW3_RANDOM_WALK_USE_FIXED_STEP)
          todo
#endif

//          std::cout << "Candidate: " << candidate_pt << std::endl;

          // allow to leave the Delaunay balls: we will make sure that the point breaks the facet later
          bool in_Delaunay_ball =
            (compare_squared_distance(cc, candidate_pt, scaled_ball_sq_radius) == CGAL::SMALLER) ||
            (compare_squared_distance(ncc, candidate_pt, scaled_neighbor_ball_sq_radius) == CGAL::SMALLER);
          bool in_local_ball = (compare_squared_distance(source, candidate_pt, seeking_sq_radius) == CGAL::SMALLER);
          if(!in_Delaunay_ball || !in_local_ball)
          {
//            std::cout << in_Delaunay_ball << " " << in_local_ball << std::endl;
            // std::cin.get();
            continue;
          }

          walk_out << "2 " << current_pt << " " << candidate_pt << "\n";
          current_pt = candidate_pt;
          closest_point = m_oracle.closest_point(current_pt);

          break;
        }
      }
    }
  }

  void spanwn_fountain(const Facet& f,
                       const Point_3& steiner_point,
                       const Point_3& source,
                       const Vector_3& dir,
                       const FT seeking_radius,
                       const std::size_t n,
                       std::vector<Point_3>& samples,
                       const bool can_spawn = true) const
  {
    static CGAL::Random rng(0);
    Random_points_on_sphere_3<Point_3> random_point_on_sphere(1., rng);

    // don't change this without changing it in the other functions
    const FT seeking_sq_radius = CGAL::square(seeking_radius);

    const FT precision = 0.01 * m_offset;
    const FT sq_offset_minus_precision = CGAL::square(m_offset - precision);
    const FT sq_offset_plus_precision = CGAL::square(m_offset + precision);

    const Plane_3 pl { steiner_point, dir };
    Vector_3 base1 = pl.base1();
    Vector_3 base2 = pl.base2();
    base1 = base1 / CGAL::approximate_sqrt(base1.squared_length());
    base2 = base2 / CGAL::approximate_sqrt(base2.squared_length());

    std::vector<Point_3> sub_sources;

    // to avoid too much structure
    FT lng_offset = rng.get_double(0., 2 * CGAL_PI);
    std::cout << "lng_offset = " << lng_offset << std::endl;

    for(unsigned int i=0; i<n; ++i)
    {
      const FT angle = lng_offset + FT(i) * 2 * CGAL_PI / FT(n);
      const FT cos_angle = std::cos(angle);
      const FT sin_angle = std::sin(angle);
      const Vector_3 lng = cos_angle * base1 + sin_angle * base2;

#ifndef CGAL_AW3_RANDOM_FOUNTAIN_STREAMS
      Vector_3 rot = CGAL::cross_product(dir, lng);
      rot = rot / CGAL::approximate_sqrt(rot.squared_length());
#endif

      // walk on spheres on the longitude given by the intersection of the seeking sphere and the plane defined by R and lng
      Point_3 current_pt = source;
      Point_3 closest_point = m_oracle.closest_point(current_pt);

      static int walk_id = 0;
      std::cout << "------------ Walk ID " << m_tr.number_of_vertices() << " " << i << std::endl;
      std::ofstream walk_out("results/walk_" + std::to_string(m_tr.number_of_vertices()) + "_" + std::to_string(walk_id++) + ".cgal");
      walk_out.precision(17);

      // until we have intersected the offset surface, or detected that there is no intersection
      FT distance_traveled = 0;
      bool has_spawned = false;

      for(;;)
      {
        std::cout << "current_pt = " << current_pt << " traveled " << distance_traveled << std::endl;

#ifndef CGAL_AW3_RANDOM_FOUNTAIN_STREAMS
        if(distance_traveled > 2 * CGAL_PI * seeking_radius) // could maybe avoid the 2 if n is even?
        {
          std::cout << "No intersection found along this great circle" << std::endl;
          break;
        }
#endif

#ifdef CGAL_AW3_SAMPLE_WITH_FOUNTAINS_RECURSIVE
        if(distance_traveled > 0.5 * seeking_radius) // sub fountain
        {
          if(!has_spawned && can_spawn)
          {
            sub_sources.push_back(current_pt);
            has_spawned = true;
          }
        }
#endif

        // if the current point is close enough to the offset, we are done
        FT sq_current_dist = squared_distance(current_pt, closest_point);
        if((sq_current_dist > sq_offset_minus_precision) &&
           (sq_current_dist < sq_offset_plus_precision))
        {
          samples.push_back(current_pt);
          break;
        }

        FT step = CGAL::approximate_sqrt(sq_current_dist) - m_offset;
        distance_traveled += step;

#ifdef CGAL_AW3_RANDOM_FOUNTAIN_STREAMS
        Vector_3 rot { CGAL::ORIGIN, *random_point_on_sphere++ };
#endif

        // create matrix of the rotation
        FT s = (std::sin)(2 * std::asin(step / (2. * seeking_radius)));
        FT c = CGAL::approximate_sqrt(1. - square(s));
        FT ux(rot.x()), uy(rot.y()), uz(rot.z());

        FT matrix[9] =
        {
          ux*ux*(1-c) + c,    ux*uy*(1-c) - uz*s, ux*uz*(1-c) + uy*s,
          ux*uy*(1-c) + uz*s, uy*uy*(1-c) + c,    uy*uz*(1-c) - ux*s,
          ux*uz*(1-c) - uy*s, uy*uz*(1-c) + ux*s, uz*uz*(1-c) + c,
        };

        // apply the rotation
        Point_3 tr_current_pt = current_pt - Vector_3(CGAL::ORIGIN, steiner_point);
        FT x = tr_current_pt.x(), y = tr_current_pt.y(), z = tr_current_pt.z();
        FT nx = matrix[0]*x + matrix[1]*y + matrix[2]*z;
        FT ny = matrix[3]*x + matrix[4]*y + matrix[5]*z;
        FT nz = matrix[6]*x + matrix[7]*y + matrix[8]*z;
        Point_3 tr_candidate_pt = Point_3(nx, ny, nz);
        Point_3 candidate_pt = tr_candidate_pt + Vector_3(CGAL::ORIGIN, steiner_point);

        // print everything
        std::cout << "step = " << step << std::endl;
        std::cout << "step / (2. * seeking_radius) = " << step / (2. * seeking_radius) << std::endl;
        std::cout << "s = " << s << std::endl;
        std::cout << "c = " << c << std::endl;
        std::cout << "rot = " << rot << std::endl;
        std::cout << "matrix = " << matrix[0] << " " << matrix[1] << " " << matrix[2] << std::endl;
        std::cout << "         " << matrix[3] << " " << matrix[4] << " " << matrix[5] << std::endl;
        std::cout << "         " << matrix[6] << " " << matrix[7] << " " << matrix[8] << std::endl;
        std::cout << "current_pt = " << current_pt << std::endl;
        std::cout << "tr_current_pt = " << tr_current_pt << std::endl;
        std::cout << "tr_candidate_pt = " << tr_candidate_pt << std::endl;
        std::cout << "candidate_pt = " << candidate_pt << std::endl;

        walk_out << "2 " << current_pt << " " << candidate_pt << "\n";
        current_pt = candidate_pt;
        closest_point = m_oracle.closest_point(current_pt);
      }
    }

#ifdef CGAL_AW3_SAMPLE_WITH_FOUNTAINS_RECURSIVE
    for(const Point_3& sub_source : sub_sources)
    {
      Vector_3 sub_dir = sub_source - steiner_point;
      spanwn_fountain(f, steiner_point, sub_source, sub_dir, seeking_radius,
                      n, samples, false /*no recursive spawning*/); // n -> 2 @tmp
    }
#endif
  }

  void get_fountain_samples(const Facet& f,
                            const Point_3& steiner_point,
                            const Steiner_construction_rule steiner_rule,
                            const FT seeking_radius,
                            const std::size_t n,
                            std::vector<Point_3>& samples) const
  {
    const Cell_handle ch = f.first;
    const int id = f.second;
    const Point_3& cc = circumcenter(ch); // @fixme inf cell...
    const FT ball_sq_radius = squared_distance(cc, m_tr.point(ch, 0));

    const Cell_handle nc = f.first->neighbor(id);
    const Point_3& ncc = circumcenter(nc);
    const FT neighbor_ball_sq_radius = squared_distance(ncc, m_tr.point(nc, 0));

    const FT seeking_sq_radius = CGAL::square(seeking_radius);

    Vector_3 dir; // points away from steiner_point

    if(steiner_rule == Steiner_construction_rule::R1)
    {
      // not dir = { steiner_point, cc } because ch might be infinite, in which case
      // cc is (awkwardly) set to the center of the diametral ball of the facet,
      // and if we use this, the dual might be inverted.
      const Point_3& a = m_tr.point(ch, Triangulation::vertex_triple_index(id, 0));
      const Point_3& b = m_tr.point(ch, Triangulation::vertex_triple_index(id, 1));
      const Point_3& c = m_tr.point(ch, Triangulation::vertex_triple_index(id, 2));
      const Vector_3 edge_1 {a, b};
      const Vector_3 edge_2 {a, c};
      dir = CGAL::cross_product(edge_1, edge_2);
    }
    else
    {
      dir = { steiner_point, ncc };
    }

    CGAL_assertion(dir != CGAL::NULL_VECTOR);
    dir = dir / CGAL::approximate_sqrt(dir.squared_length());

    const Point_3 source = steiner_point + seeking_radius * dir;

    std::cout << "c:\n"
              << "\t" << ch->vertex(0)->point() << "\n"
              << "\t" << ch->vertex(1)->point() << "\n"
              << "\t" << ch->vertex(2)->point() << "\n"
              << "\t" << ch->vertex(3)->point() << std::endl;
    std::cout << "cc: " << cc << std::endl;
    std::cout << "nc:\n"
              << "\t" << nc->vertex(0)->point() << "\n"
              << "\t" << nc->vertex(1)->point() << "\n"
              << "\t" << nc->vertex(2)->point() << "\n"
              << "\t" << nc->vertex(3)->point() << std::endl;
    std::cout << "ncc: " << ncc << std::endl;
    std::cout << "ball_radius: " << CGAL::approximate_sqrt(ball_sq_radius) << std::endl;
    std::cout << "ball_sq_radius: " << ball_sq_radius << std::endl;
    std::cout << "neighbor_ball_radius: " << CGAL::approximate_sqrt(neighbor_ball_sq_radius) << std::endl;
    std::cout << "neighbor_ball_sq_radius: " << neighbor_ball_sq_radius << std::endl;
    std::cout << "seeking_radius: " << seeking_radius << std::endl;
    std::cout << "steiner point: " << steiner_point << std::endl;
    std::cout << "seeking sphere center: " << steiner_point << std::endl;
    std::cout << "source: " << source << std::endl;

    spanwn_fountain(f, steiner_point, source, dir, seeking_radius, n, samples);
  }

  template <std::size_t n>
  void get_barycentric_samples(const Facet& f,
                               std::array<Point_3, n>& samples) const
  {
    std::array<Point_3, n> barycentric_coords = { Point_3( 1/3.,  1/3.,  1/3.),
                                                  Point_3( 4/5., 1/10., 1/10.),
                                                  Point_3(1/10.,  4/5., 1/10.),
                                                  Point_3(1/10., 1/10.,  4/5.) };

    const Cell_handle ch = f.first;
    const int id = f.second;

    const Point_3& a = m_tr.point(ch, Triangulation::vertex_triple_index(id, 0));
    const Point_3& b = m_tr.point(ch, Triangulation::vertex_triple_index(id, 1));
    const Point_3& c = m_tr.point(ch, Triangulation::vertex_triple_index(id, 2));
    const Vector_3 edge_1 {a, b};
    const Vector_3 edge_2 {a, c};
    const FT length_1 = CGAL::approximate_sqrt(edge_1.squared_length());
    const FT length_2 = CGAL::approximate_sqrt(edge_2.squared_length());

    for(int i=0; i<n; ++i)
    {
      const Point_3& b_coords = barycentric_coords[i];
      const FT mag_1 = length_1 * b_coords.x();
      const FT mag_2 = length_2 * b_coords.y();
      const Vector_3 vec_1 = mag_1 * edge_1/length_1;
      const Vector_3 vec_2 = mag_2 * edge_2/length_2;
      Point_3 sample = a;
      sample += vec_1;
      sample += vec_2;
      samples[i] = sample;
    }
  }

  // This removes samples that are redundant: if two samples (point + normal) represent
  // roughly the same plane, then only one of them is kept
  void filter_QEM_samples(std::vector<Point_3>& p_proj_subset,
                          std::vector<Vector_3>& n_proj_subset,
                          std::vector<FT>& w_subset)
  {

  }

  enum class Sampling_optimization_status
  {
    COULD_NOT_FIND_OPTIMAL_POINT = 0,
    OPTIMAL_POINT_IS_NOT_IN_CONFLICT,
    FOUND_OPTIMAL_POINT
  };

  Sampling_optimization_status optimize_with_sampling(const Facet& f,
                                                      const Point_3& steiner_point,
                                                      const Steiner_construction_rule steiner_rule,
                                                      const FT seeking_radius,
                                                      Point_3& optimized_point,
                                                      Steiner_location& optimized_point_location) const
  {
    namespace PMP = Polygon_mesh_processing;

    std::cout << "Sampled primitives are: \n";

#ifdef CGAL_AW3_SAMPLE_WITH_RANDOM_WALKS
    const int NUM_SAMPLES = 20;
    std::array<Point_3, NUM_SAMPLES> samples;
    get_random_walk_samples(f, steiner_point, steiner_rule, seeking_radius, samples);
#elif defined(CGAL_AW3_SAMPLE_WITH_FOUNTAINS)
    const int NUM_SAMPLES = 10;
    std::vector<Point_3> samples;
    get_fountain_samples(f, steiner_point, steiner_rule, seeking_radius, NUM_SAMPLES, samples);
#elif defined(CGAL_AW3_SAMPLE_WITH_BARYCENTRIC_SEEDS)
    // get the samples via barycentric points
    const int NUM_SAMPLES = 4;
    std::array<Point_3, NUM_SAMPLES> samples;
    get_barycentric_samples(f, samples);
#elif defined(CGAL_AW3_SAMPLE_WITH_RAY_INTERSECTIONS)
    const int NUM_SAMPLES = 4;
    std::vector<Point_3> samples;
    get_ray_shooting_samples(f, steiner_point, steiner_rule, seeking_radius, NUM_SAMPLES, samples);
#else
    #error "TBD"
#endif

    std::cout << "Samples:\n";
    for(const Point_3& smp : samples)
      std::cout << "\t" << smp << "\n";

    std::vector<Point_3> p_proj;
    std::vector<Vector_3> n_proj;
    std::vector<FT> w;

    // for each sample, find closest point on input and its normal, then lift to offset surface
    for(const Point_3& sample : samples)
    {
#if 1
     // The samples are already on the offset surface
      const auto closest_pair = m_oracle.tree().closest_point_and_primitive(sample);
      const Point_3 closest_pt = closest_pair.first;

      std::cout << "\tSample " << sample << " has primitive id: " << closest_pair.second.first << ", " << closest_pair.second.second << "\n";
      std::cout << "\tSample " << sample << " has closest point: " << closest_pt << "\n";

      // using the normal of the offset surface should give better stability
      Vector_3 normal = sample - closest_pt;
      const FT area = 1.;
#elif 0
      // closest point + normal of the primitive

      const auto closest_pair = m_oracle.tree().closest_point_and_primitive(sample);
      const Point_3 closest_pt = closest_pair.first;
      // given primitive id, get the primitive so you can get the datum
      std::cout << "\tSample " << sample << " has primitive id: " << closest_pair.second.first << ", " << closest_pair.second.second << "\n";

      auto dpm = m_oracle.get_m_dpm();
      const Triangle_3& p_tri = get(dpm, closest_pair.second);
      Vector_3 normal = p_tri.supporting_plane().orthogonal_vector();
      const FT area = CGAL::sqrt(p_tri.squared_area());
#elif 0
      // closest point + normal of the offset surface

      const auto closest_pair = m_oracle.tree().closest_point_and_primitive(sample);
      const Point_3 closest_pt = closest_pair.first;
      // given primitive id, get the primitive so you can get the datum
      std::cout << "\tSample " << sample << " has primitive id: " << closest_pair.second.first << ", " << closest_pair.second.second << "\n";

      // using the normal of the offset surface should give better stability
      Vector_3 normal = sample - closest_pt;
      const FT area = 1.;
#elif 0
      // Shoot a ray orthogonal to the facet f, and the normal of the offset surface

      // @todo doesn't depend on the sample, so it could be outside the loop
      const Cell_handle ch = f.first;
      const int id = f.second;
      const Point_3& a = m_tr.point(Triangulation::vertex_triple_index(id, 0));
      const Point_3& b = m_tr.point(Triangulation::vertex_triple_index(id, 1));
      const Point_3& c = m_tr.point(Triangulation::vertex_triple_index(id, 2));
      const Vector_3 edge_1 {a, b};
      const Vector_3 edge_2 {a, c};
      const Vector_3 shooting_dir = CGAL::cross_product(edge_1, edge_2);

      // @fixme "100x" isn't very robust...
      Point_3 ip;
      if(!m_oracle.first_intersection(sample, sample + 100 * shooting_dir, ip, m_offset))
        continue;

      std::cout << "\tSample " << sample << " has intersection point: " << ip << "\n";

      const auto closest_pair = m_oracle.tree().closest_point_and_primitive(ip);
      const Point_3 closest_pt = closest_pair.first;

      std::cout << "\tSample " << sample << " has primitive id: " << closest_pair.second.first << ", " << closest_pair.second.second << "\n";
      std::cout << "\tSample " << sample << " has closest point: " << closest_pt << "\n";

      // using the normal of the offset surface should give better stability
      Vector_3 normal = ip - closest_pt;
      const FT area = 1.;
#else
# error
#endif

      FT mag_norm = normal.squared_length();
      if(mag_norm != 0.)
        normal = normal / CGAL::approximate_sqrt(mag_norm);

      p_proj.push_back(closest_pt + m_offset*normal);
      n_proj.push_back(normal);
      w.push_back(area);
    }

    std::cout << p_proj.size() << " samples" << std::endl;
    for(int i=0; i<p_proj.size(); ++i)
      std::cout << "\t" << p_proj[i] << "\n";

    if(p_proj.empty())
    {
      std::cerr << "Warning: No Samples?" << std::endl;
      return Sampling_optimization_status::COULD_NOT_FIND_OPTIMAL_POINT;
    }

    std::cout << "Initial samples[" << m_tr.number_of_vertices() << "]:" << p_proj.size() << "\n";
    std::ofstream samples_out("results/samples_" + std::to_string(m_tr.number_of_vertices()) + ".xyz");
    samples_out.precision(17);
    for(int i=0; i<p_proj.size(); ++i)
    {
      std::cout << "\t" << p_proj[i] << "\n";
      samples_out << p_proj[i] << " " << n_proj[i] << "\n";
    }

#ifdef CGAL_AW3_ENHANCE_SAMPLES_WITH_QEM_POINTS
    std::vector<Point_3> additional_p_proj;
    std::vector<Vector_3> additional_n_proj;
    std::vector<FT> additional_w;
    std::vector<std::pair<Point_3, Steiner_location> > candidate_QEM_points;

    auto get_additional_QEM = [&](const std::vector<Point_3>& p_proj_subset,
                                  const std::vector<Vector_3>& n_proj_subset,
                                  const std::vector<FT>& w_subset) -> bool
    {
      std::cout << " ~~~" << std::endl;
      std::cout << "Additional QEM with " << p_proj_subset.size() << " samples" << std::endl;
      for(int i=0; i<p_proj_subset.size(); ++i)
        std::cout << "\t" << p_proj_subset[i] << "\n";

      Point_3 qem_point;
      Steiner_location qem_location;
      std::tie(qem_point, qem_location) = qem_weighted(steiner_point, p_proj_subset, n_proj_subset, w_subset);

      // Ignore QEM giving us an optimal point on a face: we only care to snap to something sharp
      if(qem_location == Steiner_location::ON_FACE)
        return false;

      for(int i=0; i<p_proj_subset.size(); ++i)
      {
        const Plane_3 pl = Plane_3(p_proj_subset[i], n_proj_subset[i]);

        // The idea is that at a convex corner, the intersection of the offset planes is
        // at distance (sqrt(2) - 1) * offset ~= 0.5 * offset of the offset surface,
        // and it shouldn't be rejected.
        const FT sq_dist = squared_distance(qem_point, pl);
        if(sq_dist > 0.25 * m_sq_offset)
          return false;
      }

      Point_3 proj_qem_point;
      Vector_3 proj_normal;
      bool successful_proj = project_to_offset(qem_point, proj_qem_point, proj_normal);
      if(!successful_proj)
        return false;

      if(!is_close_enough_to_offset_surface(proj_qem_point) ||
         !check_does_break_gate(f, proj_qem_point) ||
         seek_closer_existing_vertex(proj_qem_point, 0.25 * m_sq_alpha))
        return false;

      additional_p_proj.push_back(proj_qem_point);
      additional_n_proj.push_back(proj_normal);
      additional_w.push_back(1 /*area*/);

      candidate_QEM_points.emplace_back(proj_qem_point, qem_location);

      return true;
    };

    // for each pair of samples, compute the QEM and add it to the samples
    for(int i=0; i<p_proj.size(); ++i)
    {
      for(int j=i+1; j<p_proj.size(); ++j)
      {
        std::vector<Point_3> p_proj_subset = { p_proj[i], p_proj[j] };
        std::vector<Vector_3> n_proj_subset = { n_proj[i], n_proj[j] };
        std::vector<FT> w_subset = { w[i], w[j] };
        get_additional_QEM(p_proj_subset, n_proj_subset, w_subset);
      }
    }

    // same, for combinaisons of three samples
    for(int i=0; i<p_proj.size(); ++i)
    {
      for(int j=i+1; j<p_proj.size(); ++j)
      {
        for(int k=j+1; k<p_proj.size(); ++k)
        {
          std::vector<Point_3> p_proj_subset = { p_proj[i], p_proj[j], p_proj[k] };
          std::vector<Vector_3> n_proj_subset = { n_proj[i], n_proj[j], n_proj[k] };
          std::vector<FT> w_subset = { w[i], w[j], w[k] };
          get_additional_QEM(p_proj_subset, n_proj_subset, w_subset);
        }
      }
    }

    // QEM of all samples
    get_additional_QEM(p_proj, n_proj, w);

    //
    std::cout << "Added " << additional_p_proj.size() << " additional samples\n";
    p_proj.insert(p_proj.end(), additional_p_proj.begin(), additional_p_proj.end());
    n_proj.insert(n_proj.end(), additional_n_proj.begin(), additional_n_proj.end());
    w.insert(w.end(), additional_w.begin(), additional_w.end());
#endif

    std::cout << "Total samples[" << m_tr.number_of_vertices() << "]:" << p_proj.size() << "\n";
    std::ofstream all_samples_out("results/all_samples_" + std::to_string(m_tr.number_of_vertices()) + ".xyz");
    all_samples_out.precision(17);
    for(int i=0; i<p_proj.size(); ++i)
    {
      std::cout << "\t" << p_proj[i] << "\n";
      all_samples_out << p_proj[i] << " " << n_proj[i] << "\n";
    }

    // Samples and candidates have been computed; find the optimal point

#ifdef CGAL_AW3_STEINER_OPTIMIZATION_RETURNS_QEM_MINIMIZER
    std::tie(optimized_point, optimized_point_location) = qem_error_minimizer(p_proj, n_proj, w, p_proj);
#elif defined(CGAL_AW3_STEINER_OPTIMIZATION_RETURNS_CANDIDATE_CLOSEST_TO_CC)
    if(candidate_QEM_points.empty())
      return Sampling_optimization_status::COULD_NOT_FIND_OPTIMAL_POINT;

    std::cout << "Total candidates[" << m_tr.number_of_vertices() << "]:" << candidate_QEM_points.size() << "\n";
    std::ofstream all_candidates_out("results/all_candidates_" + std::to_string(m_tr.number_of_vertices()) + ".xyz");
    all_candidates_out.precision(17);
    for(std::size_t i=0; i<candidate_QEM_points.size(); ++i)
    {
      std::size_t gi = NUM_SAMPLES + i;
      std::cout << "\t" << p_proj[gi] << "\n";
      all_candidates_out << p_proj[gi] << " " << n_proj[gi] << "\n";
    }

    for(const auto& candidate : candidate_QEM_points)
    {
      // This is just saying: "we should prioritize QEM types that have the lowest dimension:
      // corners ideally, or creases otherwise"
      if(candidate.second != optimized_point_location)
      {
        if(candidate.second < optimized_point_location)
        {
          optimized_point = candidate.first;
          optimized_point_location = candidate.second;
        }

        continue;
      }

      // Out of the QEM points that break the gate, get the best QEM point available:
      // for now arbitraly defined as the closest from the CC of the outside cell
      const Point_3& cc = circumcenter(f.first);
      if(squared_distance(cc, candidate.first) < squared_distance(cc, optimized_point))
      {
        optimized_point = candidate.first;
        optimized_point_location = candidate.second;
      }
    }

    std::cout << "Best QEM candidate: " << optimized_point << " type " << static_cast<int>(optimized_point_location) << std::endl;

    if(steiner_point == optimized_point)
      return Sampling_optimization_status::COULD_NOT_FIND_OPTIMAL_POINT;
    return Sampling_optimization_status::FOUND_OPTIMAL_POINT;

#elif defined(CGAL_AW3_STEINER_OPTIMIZATION_RETURNS_QEM_OPTIMAL_POINT)
    Point_3 qem_point;
    std::tie(qem_point, optimized_point_location) = qem_weighted(steiner_point, p_proj, n_proj, w);
    std::cout << "\t original qem_point: " << qem_point << "\n";

    bool successful_projection = project_to_offset(qem_point, optimized_point);
    if(!successful_projection)
      return Sampling_optimization_status::COULD_NOT_FIND_OPTIMAL_POINT;
#endif

    std::cout << "\t offset qem point: " << optimized_point << "\n";
    std::cout << "\t is final equal to og steiner? " << (steiner_point == optimized_point) << "\n";

    const FT sq_near_radius = square(0.5 * m_alpha);
    if(!check_does_break_gate(f, optimized_point) ||
       seek_closer_existing_vertex(optimized_point, sq_near_radius))
    {
      return Sampling_optimization_status::OPTIMAL_POINT_IS_NOT_IN_CONFLICT;
    }

    if(steiner_point == optimized_point)
      return Sampling_optimization_status::COULD_NOT_FIND_OPTIMAL_POINT;
    return Sampling_optimization_status::FOUND_OPTIMAL_POINT;
  }

  Sampling_optimization_status optimize_with_sampling(const Facet& f,
                                                      const Point_3& steiner_point,
                                                      const Steiner_construction_rule steiner_rule,
                                                      Point_3& optimized_point,
                                                      Steiner_location& optimized_point_location) const
  {
    const FT seeking_radius = FT(1.5) * m_alpha;
    return optimize_with_sampling(f, steiner_point, steiner_rule, seeking_radius,
                                  optimized_point, optimized_point_location);
  }

  bool optimize_with_iterative_sampling(const Facet& f,
                                        const Point_3& steiner_point,
                                        const Steiner_construction_rule steiner_rule,
                                        Point_3& optimized_point,
                                        Steiner_location& optimized_point_location) const
  {
    // The idea is to first try with a larger seeking radius:
    // - If we find a good point, we can use it.
    // - If we find a good point but it's not in conflict, it tells us not to get too close
    // - If we don't find a good point, carry on with a smaller seeking radius
    FT seeking_radius = FT(1.5) * m_alpha;
    Sampling_optimization_status res = optimize_with_sampling(f, steiner_point, steiner_rule,
                                                              optimized_point, optimized_point_location);

    if(res == Sampling_optimization_status::FOUND_OPTIMAL_POINT)
      return true;

    seeking_radius = m_alpha;

    if(res == Sampling_optimization_status::OPTIMAL_POINT_IS_NOT_IN_CONFLICT)
    {
      Point_3 rejected_optimized_position = optimized_point;

      res = optimize_with_sampling(f, steiner_point, steiner_rule, seeking_radius,
                                   optimized_point, optimized_point_location);

      if(res != Sampling_optimization_status::FOUND_OPTIMAL_POINT)
        return false;

      if(squared_distance(optimized_point, rejected_optimized_position) < 0.25 * m_sq_alpha)
        return false;

      return true;
    }

    res = optimize_with_sampling(f, steiner_point, steiner_rule, seeking_radius,
                                 optimized_point, optimized_point_location);

    return (res == Sampling_optimization_status::FOUND_OPTIMAL_POINT);
  }

private:
  enum Facet_queue_status
  {
    IRRELEVANT = 0,
    ARTIFICIAL_FACET,
    TRAVERSABLE
  };

  // @speed some decent time may be spent constructing Facet (pairs) for no reason as it's always
  // just to grab the .first and .second as soon as it's constructed, and not due to API requirements
  // e.g. from DT3
  Facet_queue_status facet_status(const Facet& f) const
  {
    CGAL_precondition(!m_tr.is_infinite(f));

#ifdef CGAL_AW3_DEBUG_FACET_STATUS
    std::cout << "facet status: "
              << f.first->vertex((f.second + 1)&3)->point() << " "
              << f.first->vertex((f.second + 2)&3)->point() << " "
              << f.first->vertex((f.second + 3)&3)->point() << std::endl;
#endif

    // skip if neighbor is OUTSIDE or infinite
    const Cell_handle ch = f.first;
    const int id = f.second;
    const Cell_handle nh = ch->neighbor(id);

    if(nh->is_outside())
    {
#ifdef CGAL_AW3_DEBUG_FACET_STATUS
      std::cout << "Neighbor already outside" << std::endl;
#endif
      return IRRELEVANT;
    }

    if(m_tr.is_infinite(nh))
    {
#ifdef CGAL_AW3_DEBUG_FACET_STATUS
      std::cout << "Neighbor is infinite" << std::endl;
#endif
      return TRAVERSABLE;
    }

    // push if facet is connected to artificial vertices
    for(int i=0; i<3; ++i)
    {
      const Vertex_handle vh = ch->vertex(Triangulation::vertex_triple_index(id, i));
      if(vh->type() == AW3i::Vertex_type::BBOX_VERTEX ||
         vh->type() == AW3i::Vertex_type::SEED_VERTEX)
      {
#ifdef CGAL_AW3_DEBUG_FACET_STATUS
        std::cout << "artificial facet due to artificial vertex #" << i << std::endl;
#endif
        return ARTIFICIAL_FACET;
      }
    }

    if(m_parsimonious)
    {
      if(is_faithful(f, m_parsimony_tolerance))
      {
#ifdef CGAL_AW3_DEBUG_FACET_STATUS
        std::cout << "faithful facet" << std::endl;
#endif
        return IRRELEVANT;
      }
    }

    // skip if f min empty sphere radius is smaller than alpha
    if(is_traversable(f))
    {
#ifdef CGAL_AW3_DEBUG_FACET_STATUS
      std::cout << "traversable facet" << std::endl;
#endif
      return TRAVERSABLE;
    }

#if 0
    // try to reach stuff regardless of alpha
    if(ch->vertex(id)->type() == AW3i::Vertex_type::CREASE_VERTEX ||
       ch->vertex(id)->type() == AW3i::Vertex_type::CORNER_VERTEX)
    {
      std::cout << "Carvable for creases & corners" << std::endl;
      return TRAVERSABLE;
    }
#endif

#ifdef CGAL_AW3_DEBUG_FACET_STATUS
    std::cout << "not traversable facet" << std::endl;
#endif
    return IRRELEVANT;
  }

  double projected_distance_from_point(const Point_3& sample_point) const
  {
    const Point_3 closest_pt = m_oracle.closest_point(sample_point);
    Vector_3 vec = sample_point - closest_pt;

    const FT sqnorm = vec.squared_length();
    if(sqnorm != 0.0)
      vec /= CGAL::sqrt(CGAL::to_double(sqnorm));

    const Point_3 proj_point = closest_pt + vec * m_offset;
    const Vector_3 diff = proj_point - sample_point;
    return CGAL::approximate_sqrt(CGAL::to_double(diff.squared_length()));
  }

  bool valid_sample_point(Point_3& sample_point,
                          const double tolerance) const
  {
    const double dist = projected_distance_from_point(sample_point);
    return (dist <= tolerance);
  }

   // This function checks for the gate priority, before adding them to queue based
   // on the 4 classes we proposed.
  double gate_priority(const Facet& f) // const
  {
    CGAL_precondition(!m_tr.is_infinite(f));

    const Cell_handle ch = f.first;
    const int id = f.second;

    // std::cout << "in gate_priority() with facet " << &*f << " id: " << f.second << std::endl;

    // @tmp
    return DBL_MAX; // unsorted

    const Cell_handle neighbor = ch->neighbor(id);
    if(m_tr.is_infinite(neighbor))
      return DBL_MAX;

    // @tmp give absolute priority to gates involving an artificial vertex,
    for(int i=1; i<4; ++i)
    {
      if(ch->vertex((id+i)%4)->type() == AW3i::Vertex_type::BBOX_VERTEX ||
         ch->vertex((id+i)%4)->type() == AW3i::Vertex_type::SEED_VERTEX)
      {
        return DBL_MAX;
      }
    }

    Point_3 steiner_point;
    Steiner_location steiner_type;
    bool successful_computation = compute_steiner_point_QEM(f, steiner_point, steiner_type);
    if(!successful_computation)
      return DBL_MAX;

    const Point_3& p0 = m_tr.point(ch, (id+1)&3);
    const Point_3& p1 = m_tr.point(ch, (id+2)&3);
    const Point_3& p2 = m_tr.point(ch, (id+3)&3);
    if(Triangle_3(p0, p1, p2).is_degenerate() == 0)
      return 0;

    return geom_traits().compute_squared_radius_3_object()(p0, p1, p2);
  }

private:
  // A permissive gate is a gate that we can traverse without checking its circumradius
  enum class Facet_status
  {
    IRRELEVANT = 0,
    HAS_INFINITE_NEIGHBOR, // the cell incident to the mirrored facet is infinite (permissive)
    SCAFFOLDING, // incident to a SEED or BBOX vertex (permissive)
    TRAVERSABLE
  };

  inline const char* get_status_message(const Facet_status status)
  {
    constexpr std::size_t status_count = 4;

    // Messages corresponding to Error_code list above. Must be kept in sync!
    static const char* message[status_count] = {
      "Irrelevant facet",
      "Facet incident to infinite neighbor",
      "Facet with a bbox/seed vertex",
      "Traversable facet"
    };

    const std::size_t status_id = static_cast<std::size_t>(status);
    if(status_id > status_count || status_id < 0)
      return "Unknown status";
    else
      return message[status_id];
  }

public:
  // @speed some decent time may be spent constructing Facet (pairs) for no reason as it's always
  // just to grab the .first and .second as soon as it's constructed, and not due to API requirements
  Facet_status facet_status(const Facet& f) const
  {
    CGAL_precondition(!m_tr.is_infinite(f));

#ifdef CGAL_AW3_DEBUG_FACET_STATUS
    std::cout << "facet status: "
              << m_tr.point(f.first, Triangulation::vertex_triple_index(f.second, 0)) << " "
              << m_tr.point(f.first, Triangulation::vertex_triple_index(f.second, 1)) << " "
              << m_tr.point(f.first, Triangulation::vertex_triple_index(f.second, 2)) << std::endl;
#endif

    // skip if neighbor is "outside" or infinite
    const Cell_handle ch = f.first;
    const int id = f.second;
    CGAL_precondition(ch->label() == Cell_label::INSIDE || ch->label() == Cell_label::OUTSIDE);

    if(!ch->is_outside()) // "inside" or "manifold"
    {
#ifdef CGAL_AW3_DEBUG_FACET_STATUS
      std::cout << "Facet is inside" << std::endl;
#endif
      return Facet_status::IRRELEVANT;
    }

    const Cell_handle nh = ch->neighbor(id);
    CGAL_precondition(ch->label() == Cell_label::INSIDE || ch->label() == Cell_label::OUTSIDE);

    if(m_tr.is_infinite(nh))
    {
#ifdef CGAL_AW3_DEBUG_FACET_STATUS
      std::cout << "Infinite neighboring cell" << std::endl;
#endif
      return Facet_status::HAS_INFINITE_NEIGHBOR;
    }

    if(nh->is_outside())
    {
#ifdef CGAL_AW3_DEBUG_FACET_STATUS
      std::cout << "Neighbor already outside" << std::endl;
#endif
      return Facet_status::IRRELEVANT;
    }

    // push if facet is connected to scaffolding vertices
    for(int i=0; i<3; ++i)
    {
      const Vertex_handle vh = ch->vertex(Triangulation::vertex_triple_index(id, i));
      if(vh->type() == AW3i::Vertex_type::BBOX_VERTEX ||
         vh->type() == AW3i::Vertex_type::SEED_VERTEX)
      {
#ifdef CGAL_AW3_DEBUG_FACET_STATUS
        std::cout << "Scaffolding facet due to vertex #" << i << std::endl;
#endif
        return Facet_status::SCAFFOLDING;
      }
    }

    // skip if f min empty sphere radius is smaller than alpha
    if(is_traversable(f))
    {
#ifdef CGAL_AW3_DEBUG_FACET_STATUS
      std::cout << "traversable" << std::endl;
#endif
      return Facet_status::TRAVERSABLE;
    }

#ifdef CGAL_AW3_DEBUG_FACET_STATUS
    std::cout << "not traversable" << std::endl;
#endif
    return Facet_status::IRRELEVANT;
  }

private:
  bool push_facet(const Facet& f)
  {
    CGAL_precondition(f.first->is_outside());

#ifdef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
    // skip if f is already in queue
    if(m_queue.contains_with_bounds_check(Gate(f)))
      return false;
#endif

    // @todo could avoid useless facet_status() calls by doing it after the zombie check
    // for the unsorted priority queue, but AFAIR, it doesn't save noticeable time (and that
    // increases the queue size).
    const Facet_status status = facet_status(f);
    if(status == Facet_status::IRRELEVANT)
      return false;

#ifdef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
    const FT sqr = smallest_squared_radius_3(f, m_tr);
    const bool is_permissive = (status == Facet_status::HAS_INFINITE_NEIGHBOR ||
                                status == Facet_status::SCAFFOLDING);
    Gate new_gate(f, sqr, is_permissive);
#else
    Gate new_gate(f, m_tr);
#endif

#ifdef CGAL_AW3_COMPUTE_AND_STORE_STEINER_INFO_AT_GATE_CREATION
    if(status != Facet_status::HAS_INFINITE_NEIGHBOR)
    {
      Point_3 steiner_point;
      Steiner_status steiner_status = compute_steiner_point(new_gate, steiner_point);
      new_gate.m_steiner_status = steiner_status;
      if(steiner_status == Steiner_status::RULE_1 || steiner_status == Steiner_status::RULE_2)
        new_gate.m_steiner_point = steiner_point;
    }
#endif

#ifdef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
  m_queue.resize_and_push(new_gate);
#else
  m_queue.push(new_gate);
#endif

#ifdef CGAL_AW3_DEBUG_QUEUE
    const Cell_handle ch = f.first;
    const int s = f.second;
    const Point_3& p0 = m_tr.point(ch, Triangulation::vertex_triple_index(s, 0));
    const Point_3& p1 = m_tr.point(ch, Triangulation::vertex_triple_index(s, 1));
    const Point_3& p2 = m_tr.point(ch, Triangulation::vertex_triple_index(s, 2));

    static int gid = 0;
    std::cout << "Queue insertion #" << gid++ << "\n"
              << "  ch = " << &*ch << " (" << m_tr.is_infinite(ch) << ") " << "\n"
              << "\t" << p0 << "\n\t" << p1 << "\n\t" << p2 << std::endl;
    std::cout << "  Status: " << get_status_message(status) << std::endl;
 #ifdef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
    std::cout << "  SQR: " << sqr << std::endl;
    std::cout << "  Permissiveness: " << is_permissive << std::endl;

    CGAL_assertion(is_permissive || sqr >= m_sq_alpha);
 #endif // CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
#endif // CGAL_AW3_DEBUG_QUEUE

    return true;
  }

private:
  bool initialize(const double alpha,
                  const double offset,
                  const bool refining)
  {
#ifdef CGAL_AW3_DEBUG
    std::cout << "> Initialize..." << std::endl;
#endif

    const bool resuming = refining && (alpha == m_alpha) && (offset == m_offset);

#ifdef CGAL_AW3_DEBUG
    std::cout << "\tAlpha: " << alpha << std::endl;
    std::cout << "\tOffset: " << offset << std::endl;
    std::cout << "\tRefining? " << refining << std::endl;
    std::cout << "\tResuming? " << resuming << std::endl;
#endif

    if(!is_positive(alpha) || !is_positive(offset))
    {
#ifdef CGAL_AW3_DEBUG
      std::cerr << "Error: invalid input parameters: " << alpha << " and " << offset << std::endl;
#endif
      return false;
    }

#ifdef CGAL_AW3_DEBUG
    if(refining && alpha > m_alpha)
      std::cerr << "Warning: refining with an alpha greater than the last iteration's!" << std::endl;
    if(refining && offset != m_offset)
      std::cerr << "Warning: refining with a different offset value!" << std::endl;
#endif

    m_alpha = FT(alpha);
    m_sq_alpha = square(m_alpha);
    m_offset = FT(offset);
    m_sq_offset = square(m_offset);

    // Resuming means that we do not need to re-initialize the queue: either we have finished
    // and there is nothing to do, or the interruption was due to a user callback in the visitor,
    // and we can resume with the current queue
    if(resuming)
    {
#ifdef CGAL_AW3_DEBUG
      std::cout << "Resuming with a queue of size: " << m_queue.size() << std::endl;
#endif

      reset_manifold_labels();
      return true;
    }

#ifdef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
    m_queue.clear();
#else
    m_queue = { };
#endif

    if(refining)
    {
      // If we are reusing the triangulation, change the label of the extra elements
      // that we have added to ensure a manifold result back to external ("manifold" -> "outside")
      reset_manifold_labels();

      return initialize_from_existing_triangulation();
    }
    else
    {
      m_tr.clear();

      insert_bbox_corners();

      if(m_seeds.empty())
        return initialize_from_infinity();
      else
        return initialize_with_cavities();
    }
  }

  template <typename Visitor>
  bool alpha_flood_fill(Visitor& visitor)
  {
#ifdef CGAL_AW3_DEBUG
    std::cout << "> Flood fill..." << std::endl;
#endif

    visitor.on_flood_fill_begin(*this);

    // Explore all finite cells that are reachable from one of the initial outside cells.
    while(!m_queue.empty())
    {
      if(!visitor.go_further(*this))
        return false;

#ifdef CGAL_AW3_DEBUG_QUEUE_PP
      check_queue_sanity();
#endif

      // const& to something that will be popped, but safe as `ch` && `id` are extracted before the pop
      const Gate& gate = m_queue.top();

#ifndef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
      if(gate.is_zombie())
      {
        m_queue.pop();
        continue;
      }
#endif

      const Facet& f = gate.facet();
      const Cell_handle ch = f.first;
      const int s = f.second;
      const Cell_handle nh = ch->neighbor(s);

#ifdef CGAL_AW3_DEBUG_QUEUE
      static int fid = 0;
      std::cout << m_tr.number_of_vertices() << " triangulation vertices" << std::endl;
      std::cout << m_queue.size() << " facets in the queue" << std::endl;
      std::cout << "Face " << fid++ << " ["
                << m_tr.point(ch, Triangulation::vertex_triple_index(s, 0)) << " "
                << m_tr.point(ch, Triangulation::vertex_triple_index(s, 1)) << " "
                << m_tr.point(ch, Triangulation::vertex_triple_index(s, 2)) << "]" << std::endl;
      std::cout << "c = " << &*ch << " (inf: " << m_tr.is_infinite(ch) << ", label: " << ch->label() << "), n = "
                          << &*nh << " (inf: " << m_tr.is_infinite(nh) << ", label: " << nh->label() << ")" << "\n"

# ifdef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
      std::cout << "Priority: " << gate.priority() << " (sq alpha: " << m_sq_alpha << ")" << std::endl;
      std::cout << "Permissiveness: " << gate.is_permissive_facet() << std::endl;
# endif
#endif

#ifdef CGAL_AW3_DEBUG_DUMP_EVERY_STEP
      static int i = 0;
      std::string step_name = "results/steps/step_" + std::to_string(static_cast<int>(i)) + ".off";
      dump_triangulation_faces(step_name, true /*only_boundary_faces*/);

      std::string face_name = "results/steps/face_" + std::to_string(static_cast<int>(i++)) + ".xyz";
      std::ofstream face_out(face_name);
      face_out.precision(17);
      face_out << "3\n" << m_tr.point(ch, Triangulation::vertex_triple_index(s, 0)) << "\n"
                        << m_tr.point(ch, Triangulation::vertex_triple_index(s, 1)) << "\n"
                        << m_tr.point(ch, Triangulation::vertex_triple_index(s, 2)) << std::endl;
      face_out.close();
#endif

      CGAL_precondition(!m_tr.is_infinite(f));
      CGAL_precondition(ch->is_outside());
      CGAL_precondition(nh->label() == Cell_label::INSIDE || nh->label() == Cell_label::OUTSIDE);

      visitor.before_facet_treatment(*this, gate);

      m_queue.pop();

      if(m_tr.is_infinite(nh))
      {
        nh->set_label(Cell_label::OUTSIDE);
#ifndef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
        nh->increment_erase_counter();
#endif
        continue;
      }

      Point_3 steiner_point;
      Steiner_status ss = compute_steiner_point(ch, nh, steiner_point);
      if(ss != Steiner_status::NO_STEINER_POINT)
      {
//        std::cout << CGAL::abs(CGAL::approximate_sqrt(m_oracle.squared_distance(steiner_point)) - m_offset)
//                  << " vs " << 1e-2 * m_offset << std::endl;
        CGAL_assertion(CGAL::abs(CGAL::approximate_sqrt(m_oracle.squared_distance(steiner_point)) - m_offset) <= 1e-2 * m_offset);

        // locate cells that are going to be destroyed and remove their facet from the queue
        int li, lj = 0;
        Locate_type lt;
        const Cell_handle conflict_cell = m_tr.locate(steiner_point, lt, li, lj, nh);
        CGAL_assertion(lt != Triangulation::VERTEX);

        // Using small vectors like in Triangulation_3 does not bring any runtime improvement
        std::vector<Facet> boundary_facets;
        std::vector<Cell_handle> conflict_zone;
        boundary_facets.reserve(32);
        conflict_zone.reserve(32);

        m_tr.find_conflicts(steiner_point, conflict_cell,
                            std::back_inserter(boundary_facets),
                            std::back_inserter(conflict_zone));

#ifdef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
        // Purge the queue of facets that will be deleted/modified by the Steiner point insertion,
        // and which might have been gates
        for(const Cell_handle& cch : conflict_zone)
        {
          for(int i=0; i<4; ++i)
          {
            const Facet cf = std::make_pair(cch, i);
            if(m_queue.contains_with_bounds_check(Gate(cf)))
              m_queue.erase(Gate(cf));
          }
        }

        for(const Facet& f : boundary_facets)
        {
          const Facet mf = m_tr.mirror_facet(f); // boundary facets have incident cells in the CZ
          if(m_queue.contains_with_bounds_check(Gate(mf)))
            m_queue.erase(Gate(mf));
        }
#endif

        visitor.before_Steiner_point_insertion(*this, steiner_point);

        // Actual insertion of the Steiner point
        // We could use TDS functions to avoid recomputing the conflict zone, but in practice
        // it does not bring any runtime improvements
        Vertex_handle vh = m_tr.insert(steiner_point, lt, conflict_cell, li, lj);

        if(steiner_type == Steiner_location::ON_CORNER)
          vh->type() = AW3i::Vertex_type::CORNER_VERTEX;
        else if(steiner_type == Steiner_location::ON_CREASE)
          vh->type() = AW3i::Vertex_type::CREASE_VERTEX;
        else
          vh->type() = AW3i::Vertex_type::DEFAULT;

        visitor.after_Steiner_point_insertion(*this, vh);

        std::vector<Cell_handle> new_cells;
        new_cells.reserve(32);
        m_tr.incident_cells(vh, std::back_inserter(new_cells));
        for(const Cell_handle& new_ch : new_cells)
        {
          // std::cout << "new cell has time stamp " << new_ch->time_stamp() << std::endl;
          new_ch->set_label(m_tr.is_infinite(new_ch) ? Cell_label::OUTSIDE : Cell_label::INSIDE);
        }

        // Push all new boundary facets to the queue.
        // It is not performed by looking at the facets on the boundary of the conflict zones
        // because we need to handle internal facets, infinite facets, and also more subtle changes
        // such as a new cell being marked inside which now creates a boundary
        // with its incident "outside" flagged cell.
        for(Cell_handle new_ch : new_cells)
        {
          for(int i=0; i<4; ++i)
          {
            if(m_tr.is_infinite(new_ch, i))
              continue;

            const Cell_handle new_nh = new_ch->neighbor(i);
            if(new_nh->label() == new_ch->label()) // not on a boundary
              continue;

            const Facet boundary_f = std::make_pair(new_ch, i);
            if(new_ch->is_outside())
              push_facet(boundary_f);
            else
              push_facet(m_tr.mirror_facet(boundary_f));
          }
        }
      }
      else // no need for a Steiner point, carve through and continue
      {
        nh->set_label(Cell_label::OUTSIDE);
#ifndef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
        nh->increment_erase_counter();
#endif

        // for each finite facet of neighbor, push it to the queue
        const int mi = m_tr.mirror_index(ch, s);
        for(int i=1; i<4; ++i)
        {
          const Facet neighbor_f = std::make_pair(nh, (mi+i)&3);
          push_facet(neighbor_f);
        }
      }
    } // while(!queue.empty())

    visitor.on_flood_fill_end(*this);

    // Check that no useful facet has been ignored
    CGAL_postcondition_code(for(auto fit=m_tr.finite_facets_begin(), fend=m_tr.finite_facets_end(); fit!=fend; ++fit) {)
    CGAL_postcondition_code(  Cell_handle ch = fit->first; Cell_handle nh = fit->first->neighbor(fit->second); )
    CGAL_postcondition_code(  if(ch->label() == nh->label()) continue;)
    CGAL_postcondition_code(  Facet f = *fit;)
    CGAL_postcondition_code(  if(ch->is_inside()) f = m_tr.mirror_facet(f);)
    CGAL_postcondition(       facet_status(f) == Facet_status::IRRELEVANT);
    CGAL_postcondition_code(})

    return true;
  }

  // Any outside cell that isn't reachable from infinity is a cavity that can be discarded.
  std::size_t purge_inner_connected_components()
  {
#ifdef CGAL_AW3_DEBUG
    std::cout << "> Purge inner islands..." << std::endl;

    std::size_t pre_counter = 0;
    for(Cell_handle ch : m_tr.all_cell_handles())
      if(ch->is_outside())
        ++pre_counter;
    std::cout << pre_counter << " / " << m_tr.all_cell_handles().size() << " (pre purge)" << std::endl;
#endif

    std::size_t label_change_counter = 0;

    std::stack<Cell_handle> cells_to_visit;

    if(!m_seeds.empty())
    {
      for(const Point_3& seed : m_seeds)
      {
        Locate_type lt;
        int li, lj;
        Cell_handle ch = m_tr.locate(seed, lt, li, lj);

        if(!ch->is_outside())
        {
          std::cerr << "Warning: cell containing seed is not outside?!" << std::endl;
          continue;
        }

        cells_to_visit.push(ch);
      }
    }
    else // typical flooding from outside
    {
      cells_to_visit.push(m_tr.infinite_vertex()->cell());
    }

    while(!cells_to_visit.empty())
    {
      Cell_handle curr_c = cells_to_visit.top();
      cells_to_visit.pop();

      CGAL_assertion(curr_c->is_outside());

      if(curr_c->tds_data().processed())
        continue;
      curr_c->tds_data().mark_processed();

      for(int j=0; j<4; ++j)
      {
        Cell_handle neigh_c = curr_c->neighbor(j);
        if(neigh_c->tds_data().processed() || !neigh_c->is_outside())
          continue;

        cells_to_visit.push(neigh_c);
      }
    }

    for(Cell_handle ch : m_tr.all_cell_handles())
    {
      if(ch->tds_data().is_clear() && ch->is_outside())
      {
        ch->set_label(Cell_label::INSIDE);
#ifndef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
        ch->increment_erase_counter();
#endif
        ++label_change_counter;
      }
    }

    // reset the conflict flags
    for(Cell_handle ch : m_tr.all_cell_handles())
      ch->tds_data().clear();

#ifdef CGAL_AW3_DEBUG
    std::size_t post_counter = 0;
    for(Cell_handle ch : m_tr.all_cell_handles())
      if(ch->is_outside())
        ++post_counter;
    std::cout << post_counter << " / " << m_tr.all_cell_handles().size() << " (pre purge)" << std::endl;

    std::cout << label_change_counter << " label changes" << std::endl;
#endif

    return label_change_counter;
  }

private:
  bool is_non_manifold(Vertex_handle v) const
  {
    CGAL_precondition(!m_tr.is_infinite(v));

    bool is_non_manifold = false;

    std::vector<Cell_handle> inc_cells;
    inc_cells.reserve(64);
    m_tr.incident_cells(v, std::back_inserter(inc_cells));

    // Flood one inside and outside CC within the cell umbrella of the vertex.
    // Process both an inside and an outside CC to also detect edge pinching.
    // If there are still unprocessed afterwards, there is a non-manifoldness issue.
    //
    // Squat the conflict cell data to mark visits
    Cell_handle inside_start = Cell_handle();
    Cell_handle outside_start = Cell_handle();

    for(Cell_handle ic : inc_cells)
    {
      ic->tds_data().clear();
      if(ic->is_outside())
        outside_start = ic;
      else if(inside_start == Cell_handle())
        inside_start = ic; // can be "inside" or "manifold"
    }

    // fully inside / outside
    if(inside_start == Cell_handle() || outside_start == Cell_handle())
      return false;

    std::stack<Cell_handle> cells_to_visit;
    cells_to_visit.push(inside_start);
    cells_to_visit.push(outside_start);
    while(!cells_to_visit.empty())
    {
      Cell_handle curr_c = cells_to_visit.top();
      cells_to_visit.pop();

      if(curr_c->tds_data().processed())
        continue;
      curr_c->tds_data().mark_processed();

      int v_pos = -1;
      CGAL_assertion_code(bool res = ) curr_c->has_vertex(v, v_pos);
      CGAL_assertion(res);

      for(int j=0; j<4; ++j)
      {
        if(j == v_pos)
          continue;

        Cell_handle neigh_c = curr_c->neighbor(j);
        CGAL_assertion(neigh_c->has_vertex(v));

        if(neigh_c->tds_data().processed() ||
           is_on_outside_boundary(neigh_c, curr_c)) // do not cross the boundary
        {
          continue;
        }

        cells_to_visit.push(neigh_c);
      }
    }

    for(Cell_handle ic : inc_cells)
    {
      if(ic->tds_data().is_clear()) // <=> unvisited cell
      {
        is_non_manifold = true;
        break;
      }
    }

    // reset the conflict flags
    for(Cell_handle ic : inc_cells)
      ic->tds_data().clear();

    return is_non_manifold;
  }

  bool is_non_manifold(Cell_handle c) const
  {
    CGAL_precondition(!m_tr.is_infinite(c));

    for(int i=0; i<4; ++i)
    {
      Vertex_handle v = c->vertex(i);
      if(is_non_manifold(v))
        return true;
    }

    return false;
  }

  bool is_manifold() const
  {
    // Not the best complexity, but it's not important: this function is purely for information
    // Better complexity --> see PMP::non_manifold_vertices + throw
    for(const Vertex_handle v : m_tr.finite_vertex_handles())
      if(is_non_manifold(v))
        return true;

    return false;
  }

  // Remove bbox vertices, if they are not necessary (i.e., no "inside" incident cell)
  // This is to try and avoid having long tets with bbox vertices being tagged "inside" as part
  // of the manifold re-tagging
  bool remove_bbox_vertices()
  {
    bool do_remove = true;
    auto vit = m_tr.finite_vertices_begin();
    for(std::size_t i=0; i<8; ++i)
    {
      Vertex_handle v = vit++;

      std::vector<Cell_handle> inc_cells;
      inc_cells.reserve(64);
      m_tr.finite_incident_cells(v, std::back_inserter(inc_cells));

      for(Cell_handle c : inc_cells)
      {
        if(!c->is_outside())
        {
          do_remove = false;
          break;
        }
      }

      if(!do_remove)
        break;
    }

    std::cout << "Removing bbox vertices? " << do_remove << std::endl;
    if(!do_remove)
      return false;

    vit = m_tr.finite_vertices_begin();
    for(std::size_t i=0; i<8; ++i)
    {
      Vertex_handle v = vit++;
      m_tr.remove(v);
    }

    return true;
  }

public:
  // Not the best complexity, but it's very cheap compared to the rest of the algorithm.
  void make_manifold()
  {
#ifdef CGAL_AW3_DEBUG
    std::cout << "> Make manifold..." << std::endl;

    auto wrap_volume = [&]()
    {
      FT vol = 0;
      for(Cell_handle ch : m_tr.finite_cell_handles())
        if(!ch->is_outside())
          vol += volume(m_tr.point(ch, 0), m_tr.point(ch, 1), m_tr.point(ch, 2), m_tr.point(ch, 3));

      return vol;
    };

 #ifdef CGAL_AW3_DEBUG_DUMP_INTERMEDIATE_WRAPS
    dump_triangulation_faces("carved_tr.off", true /*only_boundary_faces*/);
 #endif

    FT base_vol = wrap_volume();
    if(!is_positive(base_vol))
      std::cerr << "Warning: wrap with non-positive volume?" << std::endl;
#endif

    // This ends up more harmful than useful after the priority queue has been introduced since
    // it usually results in a lot of flat cells into the triangulation, which then get added
    // to the mesh during manifoldness fixing.
//    remove_bbox_vertices();

    std::stack<Vertex_handle> non_manifold_vertices; // @todo sort somehow?
    for(Vertex_handle v : m_tr.finite_vertex_handles())
    {
      if(is_non_manifold(v))
      {
#ifdef CGAL_AW3_DEBUG_MANIFOLDNESS_PP
        std::cout << v->point() << " is non-manifold" << std::endl;
#endif
        non_manifold_vertices.push(v);
      }
    }

    // Some lambdas for the comparer
    auto has_scaffolding_vertex = [](Cell_handle c) -> bool
    {
      for(int i=0; i<4; ++i)
      {
        if(c->vertex(i)->type() == AW3i::Vertex_type::BBOX_VERTEX ||
           c->vertex(i)->type() == AW3i::Vertex_type::SEED_VERTEX)
        {
          return true;
        }
      }

      return false;
    };

    // This originally seemed like a good idea, but in the end it can have strong cascading issues,
    // whereas some cells with much smaller volume could have solved the non-manifoldness.
//    auto is_on_boundary = [](Cell_handle c, int i) -> bool
//    {
//      return is_on_outside_boundary(c, c->neighbor(i));
//    };
//
//    auto count_boundary_facets = [&](Cell_handle c, Vertex_handle v) -> int
//    {
//      const int v_index_in_c = c->index(v);
//      int boundary_facets = 0;
//      for(int i=0; i<3; ++i) // also take into account the opposite facet?
//      {
//        if(i == v_index_in_c)
//          continue;
//
//        boundary_facets += is_on_boundary(c, i);
//      }
//
//      return boundary_facets;
//    };

    // Experimentally, longest edge works better
//    auto sq_circumradius = [&](Cell_handle c) -> FT
//    {
//      const Point_3& cc = circumcenter(c);
//      return geom_traits().compute_squared_distance_3_object()(m_tr.point(c, 0), cc);
//    };

    // The reasoning behind using longest edge rather than volume is that we want to avoid
    // spikes (which would have a small volume), and can often happen since we do not spend
    // any care on the quality of tetrahedra.
    auto sq_longest_edge = [&](Cell_handle c) -> FT
    {
      return (std::max)({ squared_distance(m_tr.point(c, 0), m_tr.point(c, 1)),
                          squared_distance(m_tr.point(c, 0), m_tr.point(c, 2)),
                          squared_distance(m_tr.point(c, 0), m_tr.point(c, 3)),
                          squared_distance(m_tr.point(c, 1), m_tr.point(c, 2)),
                          squared_distance(m_tr.point(c, 1), m_tr.point(c, 3)),
                          squared_distance(m_tr.point(c, 2), m_tr.point(c, 3)) });
    };

#ifdef CGAL_AW3_DEBUG_MANIFOLDNESS_PP
    std::cout << non_manifold_vertices.size() << " initial NMV" << std::endl;
#endif

    while(!non_manifold_vertices.empty())
    {
#ifdef CGAL_AW3_DEBUG_MANIFOLDNESS_PP
      std::cout << non_manifold_vertices.size() << " NMV in queue" << std::endl;
#endif

      Vertex_handle v = non_manifold_vertices.top();
      non_manifold_vertices.pop();

      if(!is_non_manifold(v))
        continue;

      // Prioritize:
      // - cells without bbox vertices
      // - small cells when equal number of boundary facets
      //
      // Note that these are properties that do not depend on the cell labels, and so we only need
      // to sort once. However, if a criterion such as the number of incident inside cells were added,
      // we would need to sort after each modification of "inside"/"outside" labels.
      auto comparer = [&](Cell_handle l, Cell_handle r) -> bool
      {
        CGAL_precondition(!m_tr.is_infinite(l) && !m_tr.is_infinite(r));

        if(has_scaffolding_vertex(l) != has_scaffolding_vertex(r))
          return has_scaffolding_vertex(r);

        return sq_longest_edge(l) < sq_longest_edge(r);
      };

      std::vector<Cell_handle> inc_cells;
      inc_cells.reserve(64);
      m_tr.finite_incident_cells(v, std::back_inserter(inc_cells));

      std::vector<Cell_handle> finite_outside_inc_cells;
      finite_outside_inc_cells.reserve(64);
      std::copy_if(inc_cells.begin(), inc_cells.end(), std::back_inserter(finite_outside_inc_cells),
                   [&](Cell_handle c) -> bool { return !m_tr.is_infinite(c) && c->is_outside(); });

      // 'std::stable_sort' to have determinism without having to write something like:
      //     if(longest_edge(l) == longest_edge(r)) return ...
      // in the comparer. It's almost always a small range, so the extra cost does not matter.
      std::stable_sort(finite_outside_inc_cells.begin(), finite_outside_inc_cells.end(), comparer);

      for(Cell_handle ic : finite_outside_inc_cells)
      {
        CGAL_assertion(!m_tr.is_infinite(ic) && ic->is_outside());

        // This is where new material is added
        ic->set_label(Cell_label::MANIFOLD);

#ifdef CGAL_AW3_DEBUG_DUMP_EVERY_STEP
        static int i = 0;
        std::string step_name = "results/steps_manifold/step" + std::to_string(static_cast<int>(i++)) + ".off";
        dump_triangulation_faces(step_name, true /*only_boundary_faces*/);
#endif

        // @speed could update the manifold status while tagging
        if(!is_non_manifold(v))
          break;
      }

      CGAL_assertion(!is_non_manifold(v));

      // Check if the new material has not created a non-manifold configuration.
      // @speed this could be done on only the vertices of cells whose labels have changed.
      std::vector<Vertex_handle> adj_vertices;
      adj_vertices.reserve(64);
      m_tr.finite_adjacent_vertices(v, std::back_inserter(adj_vertices));

      for(Vertex_handle nv : adj_vertices)
        if(is_non_manifold(nv))
          non_manifold_vertices.push(nv);
    }

    CGAL_assertion_code(for(Vertex_handle v : m_tr.finite_vertex_handles()))
    CGAL_assertion(!is_non_manifold(v));

#ifdef CGAL_AW3_DEBUG
    std::size_t nm_cells_counter = 0;
    for(Cell_handle ch : m_tr.all_cell_handles())
      if(ch->label() == Cell_label::MANIFOLD)
        ++nm_cells_counter;
    std::cout << "Number of added cells: " << nm_cells_counter << std::endl;

    if(!is_zero(base_vol))
    {
      const FT manifold_vol = wrap_volume();
      const FT ratio = manifold_vol / base_vol;

      std::cout << "Volumes post-manifoldness fix:\n"
                << "before: " << base_vol << "\n"
                << "after:  " << manifold_vol << "\n"
                << "ratio:  " << ratio << std::endl;
      if(ratio > 1.1) // more than 10% extra volume
        std::cerr << "Warning: large increase of volume after manifoldness resolution" << std::endl;
    }
#endif
  }

private:
  void check_queue_sanity()
  {
    std::cout << "\t~~~ Check queue sanity ~~~" << std::endl;

    std::vector<Gate> queue_gates;
#ifdef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
    Gate previous_top_gate = m_queue.top();
#endif
    while(!m_queue.empty())
    {
      const Gate& current_gate = m_queue.top();
      queue_gates.push_back(current_gate);
      const Facet& current_f = current_gate.facet();
      const Cell_handle ch = current_f.first;
      const int id = current_f.second;
      const Point_3& p0 = m_tr.point(ch, Triangulation::vertex_triple_index(id, 0));
      const Point_3& p1 = m_tr.point(ch, Triangulation::vertex_triple_index(id, 1));
      const Point_3& p2 = m_tr.point(ch, Triangulation::vertex_triple_index(id, 2));

      std::cout << "Facet with VID " << get(Gate_ID_PM<Triangulation>(), current_gate) << "\n";
      std::cout << "\t" << p0 << "\n\t" << p1 << "\n\t" << p2 << "\n";

#ifdef CGAL_AW3_USE_SORTED_PRIORITY_QUEUE
      std::cout << "  Permissiveness: " << current_gate.is_permissive_facet() << "\n";
      std::cout << "  SQR: " << geom_traits().compute_squared_radius_3_object()(p0, p1, p2) << "\n";
      std::cout << "  Priority " << current_gate.priority() << std::endl;

      if(Less_gate()(current_gate, previous_top_gate))
        std::cerr << "Error: current gate has higher priority than the previous top" << std::endl;

      previous_top_gate = current_gate;
#else
      if(current_gate.is_zombie())
        std::cout << "Gate is a zombie!" << std::endl;
#endif

      m_queue.pop();
    }

    std::cout << "\t~~~ End queue sanity check ~~~" << std::endl;

    // Rebuild
    CGAL_assertion(m_queue.empty());
    for(const auto& g : queue_gates)
      m_queue.push(g); // no need for a resize here since the vector capacity is unchanged
  }

  void dump_triangulation_faces(const std::string filename,
                                bool only_boundary_faces = false)
  {
    std::stringstream vertices_ss;
    vertices_ss.precision(17);

    std::stringstream facets_ss;
    facets_ss.precision(17);

    std::unordered_map<Vertex_handle, std::size_t> vertex_to_id;
    std::size_t nv = 0;
    std::size_t nf = 0;

    for(auto fit=m_tr.finite_facets_begin(), fend=m_tr.finite_facets_end(); fit!=fend; ++fit)
    {
      Cell_handle ch = fit->first;
      int s = fit->second;

      Cell_handle nh = ch->neighbor(s);
      if(only_boundary_faces && ch->label() == nh->label())
        continue;

      std::array<std::size_t, 3> ids;
      for(std::size_t pos=0; pos<3; ++pos)
      {
        Vertex_handle v = ch->vertex((s+pos+1)&3);
        auto insertion_res = vertex_to_id.emplace(v, nv);
        if(insertion_res.second)
        {
          vertices_ss << m_tr.point(v) << "\n";
          ++nv;
        }

        ids[pos] = insertion_res.first->second;
      }

      facets_ss << "3 " << ids[0] << " " << ids[1] << " " << ids[2] << "\n";
      ++nf;
    }

    std::ofstream out(filename.c_str());
    out << "OFF\n" << nv << " " << nf << " 0\n";
    out << vertices_ss.str() << "\n" << facets_ss.str() << std::endl;
  }
};

} // namespace internal
} // namespace Alpha_wraps_3
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_3_INTERNAL_ALPHA_WRAP_3_H

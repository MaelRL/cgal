// Copyright (c) 2019-2022 Google LLC (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé
//
#ifndef CGAL_ALPHA_WRAP_3_INTERNAL_TRIANGLE_MESH_ORACLE_H
#define CGAL_ALPHA_WRAP_3_INTERNAL_TRIANGLE_MESH_ORACLE_H

#include <CGAL/license/Alpha_wrap_3.h>

#include <CGAL/Alpha_wrap_3/internal/Alpha_wrap_AABB_traits.h>
#include <CGAL/Alpha_wrap_3/internal/Oracle_base.h>
#include <CGAL/Alpha_wrap_3/internal/splitting_helper.h>
#include <CGAL/Alpha_wrap_3/internal/QEM_placer.h>

#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <algorithm>
#include <iostream>
#include <functional>
#include <memory>
#include <type_traits>
#include <vector>

namespace CGAL {
namespace Alpha_wraps_3 {
namespace internal {

// Just some typedefs for readability in the main oracle class
template <typename TriangleMesh, typename GT_>
struct TM_oracle_traits
{
  using Default_GT = typename GetGeomTraits<TriangleMesh>::type;

  using Base_GT = typename Default::Get<GT_, Default_GT>::type; // = Kernel, usually
  using Geom_traits = Alpha_wrap_AABB_traits<Base_GT>; // Wrap the kernel to add Ball_3 + custom Do_intersect_3
  using Point_3 = typename Geom_traits::Point_3;
  using AABB_traits = typename AABB_tree_splitter_traits<Point_3, Geom_traits>::AABB_traits;
  using AABB_tree = typename AABB_tree_splitter_traits<Point_3, Geom_traits>::AABB_tree;
};

// @speed could do a partial specialization 'subdivide = false' with simpler code for speed?
template <typename TriangleMesh,
          typename GT_ = CGAL::Default,
          typename BaseOracle = int,
          bool subdivide = true>
class Triangle_mesh_oracle
  : // this is the base that handles calls to the AABB tree
    public AABB_tree_oracle<typename TM_oracle_traits<TriangleMesh, GT_>::Geom_traits,
                            typename TM_oracle_traits<TriangleMesh, GT_>::AABB_tree,
                            typename std::conditional<
                                      /*condition*/subdivide,
                                      /*true*/Splitter_traversal_traits<typename TM_oracle_traits<TriangleMesh, GT_>::AABB_traits>,
                                      /*false*/Default_traversal_traits<typename TM_oracle_traits<TriangleMesh, GT_>::AABB_traits> >::type,
                            BaseOracle>,
    // this is the base that handles splitting input faces and inserting them into the AABB tree
    public AABB_tree_oracle_splitter<subdivide,
                                     typename TM_oracle_traits<TriangleMesh, GT_>::Point_3,
                                     typename TM_oracle_traits<TriangleMesh, GT_>::Geom_traits>
{
  using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;

  using TMOT = TM_oracle_traits<TriangleMesh, GT_>;
  using Base_GT = typename TMOT::Base_GT;

public:
  using Geom_traits = typename TMOT::Geom_traits;

private:
  using Point_3 = typename Geom_traits::Point_3;
  using Triangle_3 = typename Geom_traits::Triangle_3;

  using AABB_traits = typename TMOT::AABB_traits;
  using AABB_tree = typename TMOT::AABB_tree;
  using AABB_traversal_traits = typename std::conditional<
                                  /*condition*/subdivide,
                                  /*true*/Splitter_traversal_traits<AABB_traits>,
                                  /*false*/Default_traversal_traits<AABB_traits> >::type;

  using Oracle_base = AABB_tree_oracle<Geom_traits, AABB_tree, AABB_traversal_traits, BaseOracle>;
  using Splitter_base = AABB_tree_oracle_splitter<subdivide, Point_3, Geom_traits>;

public:
  // Constructors
  //
  // When using this constructor (and thus doing actual splitting), note that the oracle
  // will be adapted to this particular 'alpha', and so when calling again AW3(other_alpha)
  // the oracle might not have performed a split that is adapted to this other alpha value.
  Triangle_mesh_oracle(const double alpha,
                       const BaseOracle& base_oracle = BaseOracle(),
                       const Base_GT& gt = Base_GT())
    : Oracle_base(base_oracle, gt), Splitter_base(alpha)
  {
    Splitter_base::initialize_tree_property_maps(this->tree());
  }

  Triangle_mesh_oracle(const double alpha,
                       const Base_GT& gt,
                       const BaseOracle& base_oracle = BaseOracle())
    : Triangle_mesh_oracle(alpha, base_oracle, gt)
  { }

 Triangle_mesh_oracle(const BaseOracle& base_oracle,
                      const Base_GT& gt = Base_GT())
   : Triangle_mesh_oracle(0. /*alpha*/, base_oracle, gt)
 { }

 Triangle_mesh_oracle(const Base_GT& gt,
                      const BaseOracle& base_oracle = BaseOracle())
   : Triangle_mesh_oracle(0. /*alpha*/, base_oracle, gt)
 { }

 Triangle_mesh_oracle()
   : Triangle_mesh_oracle(0. /*alpha*/, BaseOracle(), Base_GT())
 { }

public:
  template <typename VPM, typename FT>
  auto get_projection(const Point_3& p,
                      const TriangleMesh& tmesh,
                      VPM vpm,
                      const FT alpha,
                      const FT offset) const
  {
    using Point_3 = typename boost::property_traits<VPM>::value_type;
    using K = typename CGAL::Kernel_traits<Point_3>::type;
    using Vector_3 = typename K::Vector_3;
    using Triangle_3 = typename K::Triangle_3;
    using Sphere_3 = typename K::Sphere_3;

    using Primitive_id = typename AABB_tree::Primitive_id;

    std::vector<Primitive_id> qem_primitives;
    this->tree().all_intersected_primitives(Sphere_3(p, CGAL::square(alpha)), // some padding @cache
                                            std::back_inserter(qem_primitives));

    std::cout << "QEM Faces intersection: " << qem_primitives.size() << std::endl;

    if(!qem_primitives.empty())
    {
      std::vector<Triangle_3> qem_triangles;
      qem_triangles.reserve(qem_primitives.size());
      for(const Primitive_id id : qem_primitives)
      {
        qem_triangles.push_back(get(this->m_dpm, id));
//        std::cout << qem_triangles.back() << std::endl;
      }

      auto optimal_point = solve_qem<Geom_traits>(p, qem_triangles.begin(), qem_triangles.end(),
                                                  this->tree(), alpha, offset);

      std::cout << "Sharpened SP: " << optimal_point.first << " type " << optimal_point.second << std::endl;
      return optimal_point;
    }

    return std::make_pair(p, 0 /*plane*/);
  }

  template <typename VPM, typename FT>
  auto get_splitting_position(typename boost::graph_traits<TriangleMesh>::edge_descriptor e,
                              const TriangleMesh& tmesh,
                              VPM vpm,
                              const FT alpha,
                              const FT offset) const
  {
    using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;
    using edge_descriptor = typename boost::graph_traits<TriangleMesh>::edge_descriptor;
    using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;

    using Point_3 = typename boost::property_traits<VPM>::value_type;
    using K = typename CGAL::Kernel_traits<Point_3>::type;
    using Vector_3 = typename K::Vector_3;
    using Ray_3 = typename K::Ray_3;
    using Sphere_3 = typename K::Sphere_3;

    using Primitive_id = typename AABB_tree::Primitive_id;

    std::cout << std::endl << "----" << std::endl;
    std::cout << "get_splitting_position of " << e << std::endl;

    auto get_cc = [&](const face_descriptor f) -> Point_3
    {
      halfedge_descriptor h = halfedge(f, tmesh);
      return CGAL::circumcenter(get(vpm, target(h, tmesh)),
                                get(vpm, target(next(h, tmesh), tmesh)),
                                get(vpm, source(h, tmesh)));
    };

//    auto get_steiner = [&](const face_descriptor f,
//                           Point_3& steiner_point) -> bool
//    {
//      Point_3 cc = get_cc(f);

//      Vector_3 v = Polygon_mesh_processing::compute_face_normal(f, tmesh);
//      v = -v; // to point inwards
//      Ray_3 r(cc, v);

//      std::cout << "CC of " << f << " is " << cc << std::endl;
//      std::cout << "Inverted normal: " << v << std::endl;

//      auto closest = this->tree().first_intersection(r);
//      if(!closest)
//      {
//        std::cerr << "Huh?" << std::endl;
//        return false;
//      }

//      const Point_3* closest_pt = boost::get<Point_3>( &closest->first );
//      if(!closest_pt)
//      {
//        std::cerr << "Hah?" << std::endl;
//        return false;
//      }

//      std::cout << "Intersection with input: " << *closest_pt << std::endl;

//      bool res = this->first_intersection(cc, *closest_pt, steiner_point, offset);
//      if(!res)
//      {
//        std::cerr << "Heh?" << std::endl;
//        std::exit(1);
//      }

//      std::cout << "Steiner: " << steiner_point << std::endl;

//      return true;
//    };

    auto get_qem_primitives = [&](const face_descriptor f,
                                  std::vector<Primitive_id>& intersections) -> void
    {
//      Point_3 steiner_point;
//      bool res = get_steiner(f, steiner_point);

//      if(!res)
//        return;

      // center the ball at the CC of the face
      Point_3 cc = get_cc(f);

      this->tree().all_intersected_primitives(Sphere_3(cc, CGAL::square(2 * alpha)), // some padding @cache
                                              std::back_inserter(intersections));

      std::cout << "Input P: " << cc << " Intersections: " << intersections.size() << std::endl;
    };

    // Get the new position
    face_descriptor f1 = face(halfedge(e, tmesh), tmesh);
    face_descriptor f2 = face(opposite(halfedge(e, tmesh), tmesh), tmesh);

    std::vector<Primitive_id> qem_primitives_1;
    get_qem_primitives(f1, qem_primitives_1);

    std::vector<Primitive_id> qem_primitives_2;
    get_qem_primitives(f2, qem_primitives_2);

    std::sort(qem_primitives_1.begin(), qem_primitives_1.end());
    std::sort(qem_primitives_2.begin(), qem_primitives_2.end());

    std::vector<Primitive_id> qem_primitives;
    std::set_intersection(qem_primitives_1.begin(), qem_primitives_1.end(),
                          qem_primitives_2.begin(), qem_primitives_2.end(),
                          back_inserter(qem_primitives));

    std::cout << "QEM Faces intersection: " << qem_primitives.size() << std::endl;

    Point_3 m = CGAL::midpoint(get(vpm, source(e, tmesh)),
                               get(vpm, target(e, tmesh)));

    if(!qem_primitives.empty())
    {
      std::vector<Triangle_3> qem_triangles;
      qem_triangles.reserve(qem_primitives.size());
      for(const Primitive_id id : qem_primitives)
      {
        qem_triangles.push_back(get(this->m_dpm, id));
//        std::cout << qem_triangles.back() << std::endl;
      }

      Point_3 opt_point = solve_qem<Geom_traits>(m, qem_triangles.begin(), qem_triangles.end(),
                                                 this->tree(), alpha, offset);

      std::cout << "Sharpened SP: " << opt_point << std::endl;
      return opt_point;
    }

    return m;
  }

public:
  template <typename NamedParameters = parameters::Default_named_parameters>
  void add_triangle_mesh(const TriangleMesh& tmesh,
                         const NamedParameters& np = CGAL::parameters::default_values())
  {
    using parameters::get_parameter;
    using parameters::choose_parameter;

    using VPM = typename GetVertexPointMap<TriangleMesh>::const_type;
    using Point_ref = typename boost::property_traits<VPM>::reference;

    CGAL_precondition(CGAL::is_triangle_mesh(tmesh));

    if(is_empty(tmesh))
    {
#ifdef CGAL_AW3_DEBUG
      std::cout << "Warning: Input is empty " << std::endl;
#endif
      return;
    }

#ifdef CGAL_AW3_DEBUG
    std::cout << "Insert into AABB tree (faces)..." << std::endl;
#endif

    VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_const_property_map(vertex_point, tmesh));
    CGAL_static_assertion((std::is_same<typename boost::property_traits<VPM>::value_type, Point_3>::value));

    Splitter_base::reserve(num_faces(tmesh));

    for(face_descriptor f : faces(tmesh))
    {
      if(Polygon_mesh_processing::is_degenerate_triangle_face(f, tmesh, np))
        continue;

      const Point_ref p0 = get(vpm, source(halfedge(f, tmesh), tmesh));
      const Point_ref p1 = get(vpm, target(halfedge(f, tmesh), tmesh));
      const Point_ref p2 = get(vpm, target(next(halfedge(f, tmesh), tmesh), tmesh));

      const Triangle_3 tr = this->geom_traits().construct_triangle_3_object()(p0, p1, p2);

      Splitter_base::split_and_insert_datum(tr, this->tree(), this->geom_traits());
    }

#ifdef CGAL_AW3_DEBUG
    std::cout << "Tree: " << this->tree().size() << " primitives (" << num_faces(tmesh) << " faces in input)" << std::endl;
#endif
  }
};

} // namespace internal
} // namespace Alpha_wraps_3
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_3_INTERNAL_TRIANGLE_MESH_ORACLE_H

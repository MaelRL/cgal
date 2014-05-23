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

#ifndef CGAL_ANISOTROPIC_MESH_3_CRITERIA_H
#define CGAL_ANISOTROPIC_MESH_3_CRITERIA_H

#include <iostream>
#include <fstream>
#include <utility>

#include <CGAL/Delaunay_traits_3.h>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

#define FOUR_POINTS(TYPE)  const TYPE &p, const TYPE &q, const TYPE &r, const TYPE &s
#define THREE_POINTS(TYPE)  const TYPE &p, const TYPE &q, const TYPE &r

template<typename K>
class Criteria_base
{
public:
  typedef typename K::FT  FT;

public:
//if 0. : criterion is not used
  const FT approximation;
  const FT squared_approximation;

  // facet
  const FT facet_circumradius;
  const FT squared_facet_circumradius;
  const FT facet_radius_edge_ratio;
  const FT squared_facet_radius_edge_radio;

  // cell
  const FT cell_circumradius;
  const FT squared_cell_circumradius;
  const FT cell_radius_edge_ratio;
  const FT squared_cell_radius_edge_radio;
  const FT sliverity;
  const bool c_consistency;

//can't be 0
  const FT distortion;
  const FT beta;
  const FT delta;

//misc
  const int nb_initial_points;
  const std::size_t max_times_to_try_in_picking_region;

  const int dimension;

public:
  void report() const
  {
    std::cout << "general criteria : " << std::endl;
    std::cout << "approximation error :               " << approximation << std::endl;
    std::cout << "maximum distortion (gamma_0):       " << distortion << std::endl;
    std::cout << "checking threshold (beta):          " << beta << std::endl;
    std::cout << "picking region radius (delta):      " << delta << std::endl;
    std::cout << "max. time to try in picking region: " << max_times_to_try_in_picking_region << std::endl;
    std::cout << "nb of initial points: " << nb_initial_points << std::endl;

    std::cout << "facet criteria : " << std::endl;
    std::cout << "maximum circumradius (r_0):         " << facet_circumradius << std::endl;
    std::cout << "radius-edge-ratio (rho_0):          " << facet_radius_edge_ratio << std::endl;

    std::cout << "cell criteria : " << std::endl;
    std::cout << "maximum circumradius (r_0):         " << cell_circumradius << std::endl;
    std::cout << "radius-edge-ratio (rho_0):          " << cell_radius_edge_ratio << std::endl;
    std::cout << "sliverity ratio (sigma_0):          " << sliverity << std::endl;
  }

  void report(std::ofstream& fx) const
  {
    fx << "general criteria : " << std::endl;
    fx << "approximation error :               " << approximation << std::endl;
    fx << "maximum distortion (gamma_0):       " << distortion << std::endl;
    fx << "checking threshold (beta):          " << beta << std::endl;
    fx << "picking region radius (delta):      " << delta << std::endl;
    fx << "max. time to try in picking region: " << max_times_to_try_in_picking_region << std::endl;
    fx << "nb of initial points: " << nb_initial_points << std::endl;

    fx << "facet criteria : " << std::endl;
    fx << "maximum circumradius (r_0):         " << facet_circumradius << std::endl;
    fx << "radius-edge-ratio (rho_0):          " << facet_radius_edge_ratio << std::endl;

    fx << "cell criteria : " << std::endl;
    fx << "maximum circumradius (r_0):         " << cell_circumradius << std::endl;
    fx << "radius-edge-ratio (rho_0):          " << cell_radius_edge_ratio << std::endl;
    fx << "sliverity ratio (sigma_0):          " << sliverity << std::endl;
  }

  bool no_cell_criteria()
  {
    return !sliverity && !cell_circumradius && !cell_radius_edge_ratio && !c_consistency;
  }

  Criteria_base(const FT approximation_ = 0.,
                const FT f_circumradius_ = 1.0,
                const FT f_radius_edge_ratio_ = 3.0,
                const FT c_circumradius_ = 1.0,
                const FT c_radius_edge_ratio_ = 3.0,
                const FT sliverity_ = 0.2,
                const FT distortion_ = 1.2,
                const FT beta_ = 2.5,
                const FT delta_ = 0.3,
                const int nb_initial_points_ = 4,
                const std::size_t max_times_to_try_in_picking_region_ = 60,
                const int dimension_ = 6
                )
    :
      approximation(approximation_),
      squared_approximation(approximation_*approximation_),
      facet_circumradius(f_circumradius_),
      squared_facet_circumradius(f_circumradius_ * f_circumradius_),
      facet_radius_edge_ratio(f_radius_edge_ratio_),
      squared_facet_radius_edge_radio(f_radius_edge_ratio_ * f_radius_edge_ratio_),
      cell_circumradius(c_circumradius_),
      squared_cell_circumradius(c_circumradius_ * c_circumradius_),
      cell_radius_edge_ratio(c_radius_edge_ratio_),
      squared_cell_radius_edge_radio(c_radius_edge_ratio_ * c_radius_edge_ratio_),
      sliverity(sliverity_),
      c_consistency(true), //this thing is useless todo
      distortion(distortion_),
      beta(beta_),
      delta(delta_),
      nb_initial_points(nb_initial_points_),
      max_times_to_try_in_picking_region(max_times_to_try_in_picking_region_),
      dimension(dimension_)
  { }
};

template<typename K, typename KExact = K>
class Stretched_criteria
{
public:
  //typedef Stretched_traits_3<K>  Traits;
  typedef Delaunay_traits_3<K, KExact>  Traits;
  typedef typename K::Point_3   Point_3;
  typedef typename K::FT        FT;
  typedef Criteria_base<K>      Criteria;

public:
  const Traits& traits;
  const Criteria* criteria;

public:
  FT compute_squared_shortest_edge(FOUR_POINTS(Point_3)) const
  {
    typename Traits::Compute_squared_distance_3 o = traits.compute_squared_distance_3_object();
    return min(min(min(min(min(o(p, q), o(p, r)), o(p, s)), o(q, r)), o(q, s)), o(r, s));
  }

  FT compute_squared_shortest_edge(THREE_POINTS(Point_3)) const
  {
    typename   Traits::Compute_squared_distance_3 o = traits.compute_squared_distance_3_object();
    return min(min(o(p, q), o(p, r)), o(q, r));
  }

  FT compute_squared_circumradius(FOUR_POINTS(Point_3)) const
  {
    typename   Traits::Compute_squared_distance_3 o = traits.compute_squared_distance_3_object();
    typename   Traits::Construct_circumcenter_3 c = traits.construct_circumcenter_3_object();
    return o(c(p, q, r, s), p);
  }

  FT compute_squared_circumradius(THREE_POINTS(Point_3)) const
  {
    typename Traits::Compute_squared_distance_3 sqd = traits.compute_squared_distance_3_object();
    typename Traits::Construct_circumcenter_3 cc = traits.construct_circumcenter_3_object();
    return sqd(cc(p, q, r), p);
  }

  FT compute_squared_radius_edge_ratio(FOUR_POINTS(Point_3)) const
  {
    return compute_squared_circumradius(p, q, r, s) / compute_squared_shortest_edge(p, q, r, s);
  }

  FT compute_squared_radius_edge_ratio(THREE_POINTS(Point_3)) const
  {
    return compute_squared_circumradius(p, q, r) / compute_squared_shortest_edge(p, q, r);
  }

  FT radius_edge_ratio_overflow(FOUR_POINTS(Point_3)) const
  {
    FT value = compute_squared_radius_edge_ratio(p, q, r, s) - criteria->squared_cell_radius_edge_radio;
    return (value < 0.0) ? 0.0 : value;
  }

  FT radius_edge_ratio_overflow(THREE_POINTS(Point_3)) const
  {
    FT value = compute_squared_radius_edge_ratio(p, q, r) - criteria->squared_facet_radius_edge_radio;
    return (value < 0.0) ? 0.0 : value;
  }

  FT compute_volume(FOUR_POINTS(Point_3)) const
  {
    typename K::Compute_volume_3 o;
    return std::abs(o(p, q, r, s));
  }

  FT compute_volume(THREE_POINTS(Point_3)) const
  {
    typename K::Compute_area_3 o;
    return std::abs(o(p, q, r));
  }

  FT compute_sliverity(FOUR_POINTS(Point_3)) const
  {
    return pow(compute_volume(p, q, r, s), 1.0 / 3.0) /
        sqrt(compute_squared_shortest_edge(p, q, r, s));
  }

  FT sliverity_overflow(FOUR_POINTS(Point_3)) const
  {
    FT value = criteria->sliverity - compute_sliverity(p, q, r, s);
    return (value < 0.0) ? 0.0 : value;
  }

  FT circumradius_overflow(FOUR_POINTS(Point_3)) const
  {
    FT value = compute_squared_circumradius(p, q, r, s) - criteria->squared_cell_circumradius;
    return (value < 0.0) ? 0.0 : value;
  }

  FT circumradius_overflow(THREE_POINTS(Point_3)) const
  {
    FT value = compute_squared_circumradius(p, q, r) - criteria->squared_facet_circumradius;
    return (value < 0.0) ? 0.0 : value;
  }

  FT element_quality(FOUR_POINTS(Point_3)) const
  {
    FT alpha = 216.*std::sqrt(3.);
    FT V = compute_volume(p, q, r, s);
    FT a = compute_volume(p, q, r);
    FT b = compute_volume(p, q, s);
    FT c = compute_volume(p, r, s);
    FT d = compute_volume(q, r, s);
    FT denom = (a+b+c+d);
    FT quality = alpha*V*V/(denom*denom*denom);

    return quality;
  }

  FT element_quality(THREE_POINTS(Point_3)) const
  {
    typename Traits::Compute_squared_distance_3 sqd =
        traits.compute_squared_distance_3_object();
    FT alpha = 4.*std::sqrt(3.);
    FT A = compute_volume(p, q, r);
    FT a = sqd(p, q);
    FT b = sqd(p, r);
    FT c = sqd(q, r);
    FT quality = alpha*A/(a+b+c);

    return quality;
  }

public:
  Stretched_criteria(const Traits &traits_, const Criteria* criteria_)
   : traits(traits_), criteria(criteria_) { }
};

#undef  FOUR_POINTS
#undef  TRHEE_POINTS

} // Anisotropic_mesh_3
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_CRITERIA_H

#ifndef CGAL_ANISOTROPIC_MESH_2_CRITERIA_H
#define CGAL_ANISOTROPIC_MESH_2_CRITERIA_H

#include <CGAL/Delaunay_traits_2.h>

#include <iostream>
#include <fstream>
#include <utility>

namespace CGAL
{
namespace Anisotropic_mesh_2
{
#define THREE_POINTS(TYPE)  const TYPE &p, const TYPE &q, const TYPE &r

template<typename K>
class Criteria_base
{
public:
  typedef typename K::FT  FT;

public:
//if 0. : criterion is not used
  // face
  const FT face_circumradius;
  const FT squared_face_circumradius;
  const FT face_radius_edge_ratio;
  const FT squared_face_radius_edge_radio;

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
    std::cout << "maximum distortion (gamma_0):       " << distortion << std::endl;
    std::cout << "checking threshold (beta):          " << beta << std::endl;
    std::cout << "picking region radius (delta):      " << delta << std::endl;
    std::cout << "max. time to try in picking region: " << max_times_to_try_in_picking_region << std::endl;
    std::cout << "nb of initial points: " << nb_initial_points << std::endl;

    std::cout << "face criteria : " << std::endl;
    std::cout << "maximum circumradius (r_0):         " << face_circumradius << std::endl;
    std::cout << "radius-edge-ratio (rho_0):          " << face_radius_edge_ratio << std::endl;
  }

  void report(std::ofstream& fx) const
  {
    fx << "general criteria : " << std::endl;
    fx << "maximum distortion (gamma_0):       " << distortion << std::endl;
    fx << "checking threshold (beta):          " << beta << std::endl;
    fx << "picking region radius (delta):      " << delta << std::endl;
    fx << "max. time to try in picking region: " << max_times_to_try_in_picking_region << std::endl;
    fx << "nb of initial points: " << nb_initial_points << std::endl;

    fx << "face criteria : " << std::endl;
    fx << "maximum circumradius (r_0):         " << face_circumradius << std::endl;
    fx << "radius-edge-ratio (rho_0):          " << face_radius_edge_ratio << std::endl;
  }

  Criteria_base(const FT f_circumradius_ = 1.0,
                const FT f_radius_edge_ratio_ = 3.0,
                const FT distortion_ = 1.2,
                const FT beta_ = 2.5,
                const FT delta_ = 0.3,
                const int nb_initial_points_ = 4,
                const std::size_t max_times_to_try_in_picking_region_ = 60,
                const int dimension_ = 6)
    :
      face_circumradius(f_circumradius_),
      squared_face_circumradius(f_circumradius_ * f_circumradius_),
      face_radius_edge_ratio(f_radius_edge_ratio_),
      squared_face_radius_edge_radio(f_radius_edge_ratio_ * f_radius_edge_ratio_),
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
  //typedef Stretched_traits_2<K>  Traits;
  typedef Delaunay_traits_2<K, KExact>  Traits;
  typedef typename K::Point_2           Point_2;
  typedef typename K::FT                FT;
  typedef Criteria_base<K>              Criteria;

public:
  const Traits& traits;
  const Criteria* criteria;

public:
  FT compute_squared_shortest_edge(THREE_POINTS(Point_2)) const
  {
    typename   Traits::Compute_squared_distance_2 o = traits.compute_squared_distance_2_object();
    return min(min(o(p, q), o(p, r)), o(q, r));
  }

  FT compute_squared_circumradius(THREE_POINTS(Point_2)) const
  {
    typename Traits::Compute_squared_distance_2 sqd = traits.compute_squared_distance_2_object();
    typename Traits::Construct_circumcenter_2 cc = traits.construct_circumcenter_2_object();
    return sqd(cc(p, q, r), p);
  }

  FT compute_squared_radius_edge_ratio(THREE_POINTS(Point_2)) const
  {
    return compute_squared_circumradius(p, q, r) / compute_squared_shortest_edge(p, q, r);
  }

  FT radius_edge_ratio_overflow(THREE_POINTS(Point_2)) const
  {
    FT value = compute_squared_radius_edge_ratio(p, q, r) - criteria->squared_face_radius_edge_radio;
    return (value < 0.0) ? 0.0 : value;
  }

  FT compute_volume(THREE_POINTS(Point_2)) const
  {
    typename K::Compute_area_2 o;
    return std::abs(o(p, q, r));
  }

  FT circumradius_overflow(THREE_POINTS(Point_2)) const
  {
    FT value = compute_squared_circumradius(p, q, r) - criteria->squared_face_circumradius;
    return (value < 0.0) ? 0.0 : value;
  }

  FT element_quality(THREE_POINTS(Point_2)) const
  {
    typename Traits::Compute_squared_distance_2 sqd =
        traits.compute_squared_distance_2_object();
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

#undef  TRHEE_POINTS

} // Anisotropic_mesh_2
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_2_CRITERIA_H

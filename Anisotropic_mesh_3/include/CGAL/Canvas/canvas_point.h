#ifndef CGAL_ANISOTROPIC_MESH_3_CANVAS_POINT_H
#define CGAL_ANISOTROPIC_MESH_3_CANVAS_POINT_H

#include <CGAL/Canvas/canvas_enum.h>

#include <CGAL/Metric.h>

#include <CGAL/assertions.h>

#include <Eigen/Dense>

#include <cstddef>
#include <iostream>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename Cp>
struct Canvas_point_comparer
{
  bool operator()(Cp const * const cp1, Cp const * const cp2)
  {
    return cp1->distance_to_closest_seed() > cp2->distance_to_closest_seed();
  }
};

template<typename K>
class Canvas_point
{
public:
  typedef typename K::FT                                    FT;
  typedef typename K::Point_3                               Point_3;
  typedef Metric_base<K>                                    Metric;
  typedef typename Eigen::Matrix<FT, 3, 1>                  Vector3d;

  Point_3 m_point;
  std::size_t m_index;
  FT m_distance_to_closest_seed;
  std::size_t m_closest_seed_id;
  FMM_state m_state;
  Metric m_metric;
  const Canvas_point* m_ancestor;

  Point_3& point() { return m_point; }
  const Point_3& point() const { return m_point; }
  std::size_t& index() { return m_index; }
  std::size_t index() const { return m_index; }
  FT& distance_to_closest_seed() { return m_distance_to_closest_seed; }
  const FT distance_to_closest_seed() const { return m_distance_to_closest_seed; }
  std::size_t& closest_seed_id() { return m_closest_seed_id; }
  std::size_t closest_seed_id() const { return m_closest_seed_id; }
  FMM_state& state() { return m_state; }
  FMM_state state() const { return m_state; }
  Metric& metric() { return m_metric; }
  const Metric& metric() const { return m_metric; }
  const Canvas_point* ancestor() { return m_ancestor; }
  const Canvas_point* ancestor() const { return m_ancestor; }

  void change_state(FMM_state new_state, std::size_t& known_count,
                    std::size_t& trial_count, std::size_t& far_count)
  {
    if(new_state == state())
      std::cerr << "WARNING: useless state change..." << std::endl;

    if(state() == FAR)
      --far_count;
    else if(state() == TRIAL)
      --trial_count;
    else if(state() == KNOWN)
      --known_count;

    if(new_state == KNOWN)
      ++known_count;
    else if(new_state == TRIAL)
      ++trial_count;
    else if(new_state == FAR)
      ++far_count;

    state() = new_state;
  }

  void initialize_from_point(const FT d,
                             const std::size_t seed_id)
  {
#if (verbosity > 15)
    std::cout << "initialize " << index() << " (" << point() << ")";
    std::cout << " at distance " << d << " from " << seed_id << std::endl;
#endif

    closest_seed_id() = seed_id;
    distance_to_closest_seed() = d;
    state() = TRIAL;
    m_ancestor = NULL;
  }

  void print_ancestor_tree() const
  {
    std::cout << "ancestors: ";
    const Canvas_point* anc = ancestor();
    while(anc)
    {
      std::cout << anc->index() << " ";
      anc = anc->ancestor();
    }
    std::cout << std::endl;
  }

  FT distortion_to_seed() const
  {
    FT gamma = 1.;
//    std::cout << "init gamma: " << gamma << " " << index << std::endl;
    const Canvas_point* curr = this;
    const Canvas_point* anc = ancestor();

    while(anc)
    {
      const Metric& m1 = anc->metric();
      const Metric& m2 = curr->metric();
      FT loc_gamma = m1.compute_distortion(m2);
//      std::cout << "loc: " << loc_gamma << " " << anc->index << std::endl;
#if 1
      gamma *= loc_gamma;
#else
      gamma = (std::max)(loc_gamma, gamma);
#endif
//      std::cout << "gamma: " << gamma << std::endl;
      anc = anc->ancestor();
    }

    gamma = (std::min)(25., gamma);

//    std::cout << "final gamma:" << gamma << std::endl;
//    exit(0);
    return gamma;
  }

  std::size_t count_ancestors() const
  {
    std::size_t i = 1;
    const Canvas_point* anc = ancestor();
    while(anc)
    {
      ++i;
      anc = anc->ancestor();
    }
    return i;
  }

  void reset_paint()
  {
    m_state = FAR;
    m_distance_to_closest_seed = FT_inf;
    m_closest_seed_id = -1;
    m_ancestor = NULL;
  }

  Canvas_point() { }

  template<typename MF>
  Canvas_point(const Point_3& point_, const std::size_t index_, const MF mf)
    :
      m_point(point_),
      m_index(index_),
      m_distance_to_closest_seed(FT_inf),
      m_closest_seed_id(-1),
      m_state(FAR),
      m_metric(mf->compute_metric(m_point)),
      m_ancestor(NULL)
  { }
};

}  // namespace Anisotropic_mesh_3
}  // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_CANVAS_POINT_H

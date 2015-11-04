#ifndef CGAL_ANISOTROPIC_MESH_3_CANVAS_POINT_H
#define CGAL_ANISOTROPIC_MESH_3_CANVAS_POINT_H

#include <CGAL/Canvas/canvas.h>
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

template<typename K, typename Canvas>
class Canvas_point
{
  typedef Canvas_point<K, Canvas>                           Self;
protected:
  typedef typename K::FT                                    FT;
  typedef typename K::Point_3                               Point_3;
  typedef Metric_base<K>                                    Metric;
  typedef typename Eigen::Matrix<FT, 3, 1>                  Vector3d;

  typedef boost::unordered_set<std::size_t>                 Point_set;

  Point_3 m_point;
  std::size_t m_index;
  FT m_distance_to_closest_seed;
  std::size_t m_closest_seed_id;
  FMM_state m_state;
  Metric m_metric;
  std::size_t m_ancestor;
  Point_set m_children;

  Canvas* m_canvas;

public:
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
  std::size_t& ancestor() { return m_ancestor; }
  std::size_t ancestor() const { return m_ancestor; }
  Canvas*& canvas() { return m_canvas; }
  const Canvas* canvas() const { return m_canvas; }
  Point_set& children() { return m_children; }
  const Point_set& children() const { return m_children; }

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

  void remove_from_children(const std::size_t c)
  {
    // 99% of the time, the child is found, but if during the initialization of
    // a new seed, you reset (at least) two points that have an ancestry relationship
    // then the child could have been reset already...
    typename Point_set::iterator it = m_children.find(c);
    if(it != m_children.end())
      m_children.quick_erase(it);
//    else
//    {
//      std::cout << "WARNING: call to remove_from_children didn't find the child (";
//      std::cout << c << " from " << index << ")" << std::endl;
//    }
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
    m_ancestor = -1;
  }

  void print_ancestor_tree() const
  {
    std::cout << index() << " has ancestors: ";
    std::size_t anc = ancestor();
    while(anc != static_cast<std::size_t>(-1))
    {
      std::cout << anc << " ";
      anc = m_canvas->canvas_points[anc].ancestor();
    }
    std::cout << std::endl;
  }

  FT distortion_to_seed() const
  {
    FT gamma = 1.;
//    std::cout << "init gamma: " << gamma << " " << index << std::endl;
    std::size_t curr = index();
    std::size_t anc = ancestor();

    while(anc != static_cast<std::size_t>(-1))
    {
      const Metric& m1 = m_canvas->canvas_points[anc].metric();
      const Metric& m2 = m_canvas->canvas_points[curr].metric();
      FT loc_gamma = m1.compute_distortion(m2);
//      std::cout << "loc: " << loc_gamma << " " << anc->index << std::endl;
#if 1
      gamma *= loc_gamma;
#else
      gamma = (std::max)(loc_gamma, gamma);
#endif
//      std::cout << "gamma: " << gamma << std::endl;
      anc = m_canvas->canvas_points[anc].ancestor();
    }

    gamma = (std::min)(25., gamma);

//    std::cout << "final gamma:" << gamma << std::endl;
//    exit(0);
    return gamma;
  }

  std::size_t count_ancestors() const
  {
    std::size_t i = 1;
    std::size_t anc = ancestor();
    while(anc != static_cast<std::size_t>(-1))
    {
      ++i;
      anc = m_canvas->get_point(anc).ancestor();
    }
    return i;
  }

  void reset_paint()
  {
    if(ancestor != static_cast<std::size_t>(-1))
      canvas()->get_point(m_ancestor).remove_from_children(index);
    typename Point_set::iterator cit = children.begin();
    typename Point_set::iterator end = children.end();
    for(; cit!=end; ++cit)
      canvas()->get_point(*cit).ancestor() = -1;

    m_state = FAR;
    m_distance_to_closest_seed = FT_inf;
    m_closest_seed_id = -1;
    m_ancestor = -1;
    children.clear();
  }

  Canvas_point() { }

  Canvas_point(const Point_3& point_, const std::size_t index_, Canvas* canvas)
    :
      m_point(point_),
      m_index(index_),
      m_distance_to_closest_seed(FT_inf),
      m_closest_seed_id(-1),
      m_state(FAR),
      m_metric(canvas->mf.compute_metric(m_point)),
      m_ancestor(-1),
      m_children(),
      m_canvas(canvas)
  { }
};

}  // namespace Anisotropic_mesh_3
}  // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_CANVAS_POINT_H

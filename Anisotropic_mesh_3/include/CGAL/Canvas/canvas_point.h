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
public:
  typedef Canvas_point<K, Canvas>                           Self;
  typedef typename K::FT                                    FT;
  typedef typename K::Point_3                               Point_3;
  typedef Metric_base<K>                                    Metric;
  typedef Eigen::Matrix<FT, 3, 1>                           Vector3d;

  typedef boost::unordered_set<std::size_t>                 Point_set;

protected:
  Canvas* m_canvas;

  Point_3 m_point;
  std::size_t m_index;
  Metric m_metric;
  bool m_is_on_domain_border;

  // stuff that depends on the seeds :
  FMM_state m_state;
  FT m_distance_to_closest_seed;
  std::size_t m_depth;
  int m_is_seed_holder; // ugly and awkward (should use the fact that it has no ancestor instead)
  std::size_t m_closest_seed_id;
  std::size_t m_ancestor;

  // 'children' needs to be 'mutable' because compute_closest_seed takes a const ref
  // to an ancestor (because some ancestors are temporaries whose lifestime I extend
  // through const refs) and we need to modify children in compute_closest_seed
  mutable Point_set m_children;

public:
  Point_3& point() { return m_point; }
  const Point_3& point() const { return m_point; }
  std::size_t& index() { return m_index; }
  std::size_t index() const { return m_index; }
  FT& distance_to_closest_seed() { return m_distance_to_closest_seed; }
  const FT distance_to_closest_seed() const { return m_distance_to_closest_seed; }
  int& is_seed_holder() { return m_is_seed_holder; }
  int is_seed_holder() const { return m_is_seed_holder; }
  std::size_t& depth() { return m_depth; }
  std::size_t depth() const { return m_depth; }
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
  bool border_info() const { return m_is_on_domain_border; }
  bool& border_info() { return m_is_on_domain_border; }

  void change_state(FMM_state new_state, std::size_t& known_count,
                    std::size_t& trial_count, std::size_t& far_count)
  {
    if(new_state == m_state)
      std::cout << "WARNING: useless state change..." << std::endl;

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

  void reset_descendants()
  {
    //  std::cout << "reset descendant at: " << index << std::endl;
    while(!m_children.empty())
    {
      Canvas_point& cp = m_canvas->get_point(*(m_children.begin()));
      cp.reset_descendants();
    }
    reset_paint();
  }

  void initialize_from_point(const FT d,
                             const std::size_t seed_id)
  {
#if (VERBOSITY > 15)
    std::cout << "initialize " << m_index << " (" << m_point << ")";
    std::cout << " at distance " << d << " from " << seed_id << std::endl;
#endif
    reset_paint();

    m_is_seed_holder = seed_id;
    m_closest_seed_id = seed_id;
    m_distance_to_closest_seed = d;
    m_state = TRIAL;
  }

  void print_ancestor_tree() const
  {
    std::cout << m_index << " has ancestors: ";
    std::size_t anc = m_ancestor;
    while(anc != static_cast<std::size_t>(-1))
    {
      std::cout << anc << " ";
      anc = m_canvas->get_point(anc).ancestor();
    }
    std::cout << std::endl;
  }

  FT distortion_to_seed() const
  {
    FT gamma = 1.;
//    std::cout << "init gamma: " << gamma << " " << index << std::endl;
    std::size_t curr = m_index;
    std::size_t anc = m_ancestor;

    while(anc != static_cast<std::size_t>(-1))
    {
      const Metric& m1 = m_canvas->get_point(anc).metric();
      const Metric& m2 = m_canvas->get_point(curr).metric();
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
    std::size_t anc = m_ancestor;
    while(anc != static_cast<std::size_t>(-1))
    {
      ++i;
      anc = m_canvas->get_point(anc).ancestor();
    }
    return i;
  }

  std::size_t ancestor_path_length() const
  {
    std::size_t i = 1;
    std::size_t n_anc = m_ancestor;
    while(n_anc != static_cast<std::size_t>(-1))
    {
      ++i;
      const Canvas_point& anc = m_canvas->get_point(n_anc);
      n_anc = anc.ancestor();

      CGAL_assertion(i < m_canvas->canvas_points.size());
    }
    return i;
  }

  void reset_paint()
  {
    if(m_ancestor != static_cast<std::size_t>(-1))
      m_canvas->get_point(m_ancestor).remove_from_children(m_index);

    typename Point_set::iterator cit = m_children.begin();
    typename Point_set::iterator end = m_children.end();
    for(; cit!=end; ++cit)
      m_canvas->get_point(*cit).ancestor() = -1;

    m_state = FAR;
    m_distance_to_closest_seed = FT_inf;
    m_depth = 0;
    m_is_seed_holder = -1;
    m_closest_seed_id = -1;
    m_ancestor = -1;
    m_children.clear();
  }

  Canvas_point() { }
  Canvas_point(const Point_3& point_, const std::size_t index_,
               Canvas* canvas, bool border_info_ = false)
    :
      m_canvas(canvas),
      m_point(point_),
      m_index(index_),
      m_metric(canvas->mf->compute_metric(m_point)),
      m_is_on_domain_border(border_info_),
      m_state(FAR),
      m_distance_to_closest_seed(FT_inf),
      m_depth(0),
      m_is_seed_holder(-1),
      m_closest_seed_id(-1),
      m_ancestor(-1),
      m_children()
  { }
};

}  // namespace Anisotropic_mesh_3
}  // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_CANVAS_POINT_H

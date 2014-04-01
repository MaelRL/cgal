#ifndef CGAL_ANISOTROPIC_MESH_3_MESHER_LEVEL_H
#define CGAL_ANISOTROPIC_MESH_3_MESHER_LEVEL_H

#include <CGAL/IO/Star_set_output.h>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

struct Null_anisotropic_mesher_level
{
  template<typename V>
  bool refine(V) { return true; }

  bool is_algorithm_done() { return true; }

  template<typename V>
  bool one_step(V) { return false; }

  template <typename P>
  bool test_point_conflict_from_superior(P) { return false; }

  template <typename I>
  void fill_ref_queues_from_superior(typename std::set<I>, I) { }
};

template < typename Star, /* Stretched DT */
           typename Derived, /* class that implements methods. */
           typename Previous /* = Null_anisotropic_mesher_level */>
class Anisotropic_mesher_level
{
public:
  typedef Anisotropic_mesher_level<Star, Derived, Previous>      Self;
  typedef typename Star::Point_3                                 Point;

private:
  Previous& previous_level;
  bool m_is_active;

private:
  Derived& derived() { return static_cast<Derived&>(*this); }
  const Derived& derived() const { return static_cast<const Derived&>(*this); }

public:
  const Previous& previous() const { return previous_level; }
  const bool is_active() const { return m_is_active; }
  bool& is_active() { return m_is_active; }

  void initialize()
  {
    derived().initialize_();
  }

  bool get_refinement_point_for_next_element(Point& p)
  {
    return derived().get_refinement_point_for_next_element_(p);
  }

  bool test_point_conflict_from_superior(const Point& p)
  {
    return derived().test_point_conflict_from_superior_(p);
  }

  bool is_point_in_conflict(const Point& p)
  {
    return previous_level.test_point_conflict_from_superior(p);
  }

  void fill_ref_queues_from_superior(const std::set<typename Star::Index>& modified_stars,
                                    typename Star::Index pid)
  {
    derived().fill_refinement_queue(modified_stars, pid);
    previous_level.fill_ref_queues_from_superior(modified_stars, pid);
  }

  void fill_previous_ref_queues(const std::set<typename Star::Index>& modified_stars,
                                typename Star::Index pid)
  {
    previous_level.fill_ref_queues_from_superior(modified_stars, pid);
  }

  template<typename Visitor>
  void fill_refinement_queues(const std::set<typename Star::Index>& modified_stars,
                              typename Star::Index pid,
                              const Visitor& visitor)
  {
    fill_previous_ref_queues(modified_stars, pid);
    derived().fill_refinement_queue(modified_stars, pid);
    visitor.fill_refinement_queue(modified_stars, pid);
  }

  template<typename Visitor>
  bool insert(const Point& p,
              const Visitor& visitor)
  {
    return derived().insert_(p, visitor);
  }

  template<typename Visitor>
  bool process_one_element(const Visitor& visitor)
  {
    Point p;
    if(!get_refinement_point_for_next_element(p))
      return false; // no next element

    if(!is_point_in_conflict(p)) //tmp
      return insert(p, visitor); // true for correct insertion, false if problems

    return true; // conflict and no point was inserted but the algorithm continues
  }

  bool is_algorithm_done()
  {
    return ( previous_level.is_algorithm_done() && derived().is_algorithm_done_() );
  }

  template<typename Visitor>
  bool refine(const Visitor& visitor) //boolean return type so that [stop at a prev level] => [immediate stop at all levels]
  {
    while(!is_algorithm_done() /*&& derived().continue_ smthg like that*/ )
    {
      if(!previous_level.is_algorithm_done())
        if(!previous_level.refine(visitor.previous()))
          return false;
      if(!process_one_element(visitor))
        return is_algorithm_done();
    }
    return true;
  }

  template<typename Visitor>
  bool one_step(const Visitor& visitor)
  {
    if(!previous_level.is_algorithm_done())
      return previous_level.one_step(visitor.previous());
    else
      return process_one_element(visitor);
  }

  Anisotropic_mesher_level(Previous& previous) : previous_level(previous), m_is_active(false) { }
}; // Anisotropic_mesher_level

} // Anisotropic_mesh_3
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_MESHER_LEVEL_H

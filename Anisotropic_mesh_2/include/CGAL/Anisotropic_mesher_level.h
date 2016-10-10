#ifndef CGAL_ANISOTROPIC_MESH_2_MESHER_LEVEL_H
#define CGAL_ANISOTROPIC_MESH_2_MESHER_LEVEL_H

#include <CGAL/IO/Star_set_output.h>

namespace CGAL
{
namespace Anisotropic_mesh_2
{

enum Refinement_point_status
{
  EMPTY_QUEUE = 0,
  POINT_IN_CONFLICT,
  PICK_VALID_FAILED,
  SUITABLE_POINT
};

struct Null_anisotropic_mesher_level
{
  template<typename V>
  bool refine(V) { return true; }

  bool is_algorithm_done() { return true; }

  template<typename V>
  bool one_step(V) { return false; }

  template <typename P>
  bool test_point_conflict_from_superior(P, bool, bool) { return false; }

  template <typename I>
  void fill_ref_queues_from_superior(I) { }

  Null_anisotropic_mesher_level() { }
};

template < typename Star, /* Stretched DT */
           typename Derived, /* class that implements methods. */
           typename Previous /* = Null_anisotropic_mesher_level */>
class Anisotropic_mesher_level
{
public:
  typedef Anisotropic_mesher_level<Star, Derived, Previous>      Self;
  typedef typename Star::Point_2                                 Point;

private:
  Previous& previous_level;
  bool m_is_active;

private:
  Derived& derived() { return static_cast<Derived&>(*this); }
  const Derived& derived() const { return static_cast<const Derived&>(*this); }

public:
  const Previous& previous() const { return previous_level; }
  const bool& is_active() const { return m_is_active; }
  bool& is_active() { return m_is_active; }

  void initialize()
  {
    derived().initialize_();
  }

  Refinement_point_status get_refinement_point_for_next_element(Point& p)
  {
    return derived().get_refinement_point_for_next_element_(p);
  }

  bool test_point_conflict_from_superior(const Point& p,
                                         const bool is_queue_updated = true,
                                         const bool need_picking_valid = false)
  {
    //can't test conflicts at all lower levels if we're interlacing surface
    //and mesher levels (since a facet refining point obviously encroaches a
    //facet) todo
    return derived().test_point_conflict_from_superior_(p, is_queue_updated,
                                                         need_picking_valid);
  }

  //this function potentially fills partially the conflict_zones:
  //it gives which stars are in conflict (without computing the conflict hole)
  bool is_point_in_conflict(const Point& p,
                            const bool is_queue_updated = true,
                            const bool need_picking_valid = false) const
  {
    return previous_level.test_point_conflict_from_superior(p, is_queue_updated,
                                                            need_picking_valid);
  }

  void fill_ref_queues_from_superior(typename Star::Index pid)
  {
    fill_previous_ref_queues(pid);
    derived().fill_refinement_queue(pid);
  }

  void fill_previous_ref_queues(typename Star::Index pid)
  {
    previous_level.fill_ref_queues_from_superior(pid);
  }

  template<typename Visitor>
  bool after_insertion(typename Star::Index pid,
                       Visitor& visitor)
  {
    fill_previous_ref_queues(pid);
    derived().fill_refinement_queue(pid);
    visitor.fill_refinement_queue(pid);

    derived().clear_conflict_zones();

    if(pid % 75 == 0)
      derived().clean_stars();

#if 0//def ANISO_DEBUG_QUEUE
    std::cout << "Enter fill ref queues debug. Filling with all stars" << std::endl;
    derived().fill_refinement_queue();
    derived().refine_queue().print();
    std::cout << "End fill ref queues debug" << std::endl;
#endif

    return true;
  }

  bool insert(const Point& p)
  {
    return derived().insert_(p);
  }

  template<typename Visitor>
  bool process_one_element(Visitor& visitor)
  {
    Point p;
    Refinement_point_status status = get_refinement_point_for_next_element(p);

    if(status == EMPTY_QUEUE)
      return false;
    if(status == POINT_IN_CONFLICT)
      return true; // no point was inserted but the algorithm continues (at the lower level)

    return (insert(p) && // true for correct insertion, false if problems
            after_insertion(derived().number_of_stars()-1, visitor)); //id of the last inserted star is size()-1
  }

  bool is_algorithm_done()
  {
    if(!m_is_active)
      initialize();

    return previous_level.is_algorithm_done() && derived().is_algorithm_done_();
  }

  //boolean return type so that [stop at a prev level] => [immediate stop at all levels]
  template<typename Visitor>
  bool refine(Visitor& visitor)
  {
    while(!is_algorithm_done())
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
  bool one_step(Visitor& visitor)
  {
    if(!previous_level.is_algorithm_done())
      return previous_level.one_step(visitor.previous());
    else
    {
      if(!m_is_active)
        initialize();

      return process_one_element(visitor);
    }
  }

  void resume_from_mesh_file(const char* filename)
  {
    derived().resume_from_mesh_file_(filename);
  }

  Anisotropic_mesher_level(Previous& previous)
    :
      previous_level(previous),
      m_is_active(false)
  { }
}; // Anisotropic_mesher_level

} // Anisotropic_mesh_2
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_2_MESHER_LEVEL_H

#ifndef CGAL_ANISOTROPIC_MESH_3_MESHER_VISITOR_H
#define CGAL_ANISOTROPIC_MESH_3_MESHER_VISITOR_H

namespace CGAL
{
namespace Anisotropic_mesh_3
{

struct Null_anisotropic_mesher_visitor
{
  typedef Null_anisotropic_mesher_visitor Self;

  const Self& previous() const { return *this; }

  void fill_refinement_queue(int) const { }

  Null_anisotropic_mesher_visitor() { }
};


template<typename Previous>
class Null_anisotropic_mesher_visitor_level
{
  Previous& previous_level;

public:
  const Previous& previous() const { return previous_level; }

  void fill_refinement_queue(int) const { }

  Null_anisotropic_mesher_visitor_level(Previous& previous_level_)
    : previous_level(previous_level_)
  { }
};

template <typename Mesher_level,
          typename Previous>
class Anisotropic_mesher_visitor
{
  Mesher_level& mesher_level;
  Previous& previous_level;
  bool m_is_active;

public:
  const Previous& previous() const { return previous_level; }
  const bool is_active() const { return m_is_active; }
  bool& is_active() { return m_is_active; }

  void fill_refinement_queue(int relative_point) const //should be Index (same above) todo
  {
    if(m_is_active)
      mesher_level.fill_refinement_queue(relative_point);
  }

  Anisotropic_mesher_visitor(Mesher_level& mesher_level_,
                             Previous& previous_level_)
    : mesher_level(mesher_level_),
      previous_level(previous_level_),
      m_is_active(false)
  { }

}; // Anisotropic_mesher_visitor

} // Anisotropic_mesh_3
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_MESHER_VISITOR_H

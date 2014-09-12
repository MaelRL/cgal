#ifndef CGAL_ANISOTROPIC_MESH_2_MESHER_VISITOR_H
#define CGAL_ANISOTROPIC_MESH_2_MESHER_VISITOR_H

#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_2
{

struct Null_anisotropic_mesher_visitor
{
  typedef Null_anisotropic_mesher_visitor Self;

  Self& previous() { return *this; }
  const Self& previous() const { return *this; }

  template<typename I>
  void fill_refinement_queue(I) const { }

  Null_anisotropic_mesher_visitor() { }
};


template<typename Previous>
class Null_anisotropic_mesher_visitor_level
{
  Previous& previous_level;

public:
  Previous& previous() { return previous_level; }
  const Previous& previous() const { return previous_level; }

  void fill_refinement_queue(int) const { }

  Null_anisotropic_mesher_visitor_level(Previous& previous_level_)
    : previous_level(previous_level_)
  { }
};

template <typename Trunk,
          typename Previous>
class Anisotropic_mesher_visitor
{
  Previous& previous_level;
  const std::vector<Trunk*> mesher_levels;

public:
  Previous& previous() { return previous_level; }
  const Previous& previous() const { return previous_level; }

  template<typename Index>
  void fill_refinement_queue(Index relative_point)
  {
    for(std::size_t i=0; i<mesher_levels.size(); ++i)
      if(mesher_levels[i]->is_active())
        mesher_levels[i]->fill_refinement_queue(relative_point);
  }

  Anisotropic_mesher_visitor(std::vector<Trunk*> mesher_level_,
                             Previous& previous_level_)
    :
      previous_level(previous_level_),
      mesher_levels(mesher_level_)
  { }

}; // Anisotropic_mesher_visitor


} // Anisotropic_mesh_2
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_2_MESHER_VISITOR_H

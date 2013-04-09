
#ifndef CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_SURFACE_MESHER_3_H
#define CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_SURFACE_MESHER_3_H

#include <CGAL/Surface_star_set_3.h>

namespace CGAL 
{
  namespace Anisotropic_mesh_3 
  {

    template<typename K, typename RefinementCondition = No_condition<typename K::Point_3> >
    class Anisotropic_surface_mesher_3 
    {
    public:
      typedef Surface_star_set_3<K, RefinementCondition>  Star_set;

    public:
      Star_set      &star_set;

    public:
      void refine_all(const int max_count = INT_MAX) 
      {
        star_set.refine_all(max_count);
      }

      const Star_set& starset() const 
      {
        return star_set;
      }

      void output(const std::string& filename) const
      {
        star_set.output(filename.data());
      }

#ifdef CGAL_MESH_3_MESHER_STATUS_ACTIVATED
      struct Mesher_status
      { 
        std::size_t vertices, facet_queue, cells_queue;

        Mesher_status(std::size_t v, std::size_t f, std::size_t c)
          : vertices(v), facet_queue(f), cells_queue(c) {}
      };

      Mesher_status status() const;
#endif


    public:
      Anisotropic_surface_mesher_3(Star_set &starset_)
        : star_set(starset_) 
      { }

      ~Anisotropic_surface_mesher_3() {}
    };
  }
}

#endif // CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_SURFACE_MESHER_3_H

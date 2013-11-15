
#ifndef CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_TET_MESHER_3_H
#define CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_TET_MESHER_3_H

#include <CGAL/Cell_star_set_3.h>

namespace CGAL 
{
  namespace Anisotropic_mesh_3 
  {

    template<typename K>
    class Anisotropic_tet_mesher_3 
    {
    public:
      typedef Cell_star_set_3<K>   Star_set;

    public:
      Star_set      &star_set;

    public:
      void refine_all()
      {
        star_set.refine_all();
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
      Anisotropic_tet_mesher_3(Star_set &star_set_)
        : star_set(star_set_)
        { }
      
      ~Anisotropic_tet_mesher_3() { }
    };
  }
}

#endif // CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_TET_MESHER_3_H


#ifndef CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_TET_MESHER_3_H
#define CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_TET_MESHER_3_H

#include <CGAL/Constrain_surface_3.h>
#include <CGAL/Metric_field.h>
#include <CGAL/Cell_star_set_3.h>
#include <CGAL/Criteria.h>

namespace CGAL 
{
  namespace Anisotropic_mesh_3 
  {

    template<typename K>
    class Anisotropic_tet_mesher_3 
    {
    public:
      typedef Constrain_surface_3<K>  Domain;
      typedef Metric_field<K>         Metric_field;
      typedef Criteria_base<K>        Criteria;
      typedef Cell_star_set_3<K>   Star_set;

    public:
      Star_set      star_set;
      Domain        &domain;
      Criteria      &criteria;
      Metric_field  &metric_field;

    public:
      void refine() 
      {
        star_set.refine_all();
      }

      const Star_set& starset() const 
      {
        return star_set;
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
      Anisotropic_tet_mesher_3(Domain &domain_, 
                               Metric_field &metric_field_, 
                               Criteria &criteria_)
        : domain(domain_), 
          metric_field(metric_field_), 
          criteria(criteria_),
          star_set(criteria_, metric_field_, &domain_, true) 
        { }
      
      ~Anisotropic_tet_mesher_3() { }
    };
  }
}

#endif // CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_TET_MESHER_3_H

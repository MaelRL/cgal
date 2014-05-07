#ifndef CGAL_ANISOTROPIC_MESH_3_MESHER_H
#define CGAL_ANISOTROPIC_MESH_3_MESHER_H

#include <CGAL/Starset.h>
#include <CGAL/Stretched_Delaunay_3.h>
#include <CGAL/Constrain_surface_3.h>
#include <CGAL/Metric_field.h>
#include <CGAL/Criteria.h>

#include <CGAL/Anisotropic_refine_facets_3.h>
#include <CGAL/Anisotropic_refine_cells_3.h>
#include <CGAL/Anisotropic_mesher_visitor.h>

#include <CGAL/IO/Star_set_output.h>

#include <CGAL/aabb_tree/aabb_tree_bbox.h>
#include <CGAL/aabb_tree/aabb_tree_bbox_primitive.h>

#include <CGAL/kd_tree/Kd_tree_for_star_set.h>

#include <CGAL/Timer.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/iterator.h>
#include <CGAL/helpers/timer_helper.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdio.h>
#include <utility>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename PointType>
struct No_condition
{
  bool operator()(const PointType& p) const
  {
    return true;
  }
};

template<typename K, typename RefinementCondition = No_condition<typename K::Point_3> >
class Anisotropic_mesher_3
{
public:
  //typedef typename CGAL::Exact_predicates_exact_constructions_kernel KExact;
  typedef K                                                 KExact;

  typedef Anisotropic_mesher_3<K, RefinementCondition>      Self;

  typedef Stretched_Delaunay_3<K, KExact>                   Star;
  typedef Star*                                             Star_handle;
  typedef typename Star::Base                               DT; // DT_3 with vertex_base_with_info
  typedef std::vector<Star_handle>                          Star_vector;
  typedef std::set<Star_handle>                             Star_set;
  typedef typename Star_vector::iterator                    Star_iterator;
  typedef typename Star::Index                              Index;
  typedef std::set<Index>                                   Index_set;

  typedef Constrain_surface_3<K>                            Constrain_surface;
  typedef CGAL::Anisotropic_mesh_3::Metric_field<K>         Metric_field;
  typedef typename Metric_field::Metric                     Metric;
  typedef Criteria_base<K>                                  Criteria;

  typedef CGAL::Anisotropic_mesh_3::Starset<K>              Starset;

//Filters
  typedef CGAL::AABB_tree_bbox<K, Star>                     AABB_tree;
  typedef CGAL::AABB_bbox_primitive<Star>                   AABB_primitive;

  typedef CGAL::Kd_tree_for_star_set<K, Star_handle>        Kd_tree;
  typedef typename Kd_tree::Traits                          Kd_traits;
  typedef typename Kd_tree::Box_query                       Kd_Box_query;
  typedef typename Kd_tree::key_type                        Kd_point_info;

  typedef CGAL::Anisotropic_mesh_3::Conflict_zone<K>        Conflict_zone;
  typedef CGAL::Anisotropic_mesh_3::Stars_conflict_zones<K> Stars_conflict_zones;

//Mesher levels
  typedef Null_anisotropic_mesher_level                     Null_mesher_level;
  typedef Anisotropic_refine_facets_3<K, Null_mesher_level> Facets_level;
  typedef Anisotropic_refine_cells_3<K, Facets_level>       Cells_level;

//Visitors
  typedef Null_anisotropic_mesher_visitor                   Null_mesher_visitor;
  typedef Anisotropic_mesher_visitor<Cells_level, Null_mesher_visitor>
                                                            Facets_visitor;
  typedef Null_anisotropic_mesher_visitor_level<Facets_visitor>
                                                            Cells_visitor;

private:
  // Star set
  Starset& m_starset;

  // Filters
  AABB_tree m_aabb_tree; //bboxes of stars
  DT m_ch_triangulation;
  Kd_tree m_kd_tree;     //stars* centers for box queries

  mutable Stars_conflict_zones m_star_czones; //conflict zones for the stars in conflict

  // Meshers
  Null_mesher_level m_null_mesher;
  Facets_level m_facet_mesher;
  Cells_level m_cell_mesher;

  // Visitors
  Null_mesher_visitor m_null_visitor;
  Facets_visitor m_facet_visitor;
  Cells_visitor m_cell_visitor;

  RefinementCondition m_refinement_condition; //todo

public:
  double refine_mesh()
  {
    CGAL::Timer timer;
    timer.start();
    double elapsed_time = 0.;

#if 1//ndef ANISO_VERBOSE
    // Scan surface and refine it
    m_facet_mesher.initialize();
    m_facet_mesher.refine(m_facet_visitor);

    m_facet_visitor.is_active() = true;
    m_facet_mesher.pick_valid_uses_3D_checks() = true;

    // Then scan volume and refine it
    m_cell_mesher.initialize();
    m_cell_mesher.refine(m_cell_visitor);

#else
    m_facet_mesher.initialize();

    elapsed_time += timer.time();
    timer.stop(); timer.reset(); timer.start();

    while (!m_facet_mesher.is_algorithm_done())
      m_facet_mesher.one_step(m_facet_visitor);

    std::cout << "Total refining surface time: " << timer.time() << "s" << std::endl;
    elapsed_time += timer.time();
    timer.stop(); timer.reset(); timer.start();

    m_facet_visitor.is_active() = true;
    m_facet_mesher.pick_valid_uses_3D_checks() = true;

    std::cout << "Start volume scan...";
    m_cell_mesher.initialize();
    std::cout << std::endl << "end scan" << std::endl;

    elapsed_time += timer.time();
    timer.stop(); timer.reset(); timer.start();

    while (!m_cell_mesher.is_algorithm_done())
      m_cell_mesher.one_step(m_cell_visitor);

    std::cout << "Total refining volume time: " << timer.time() << "s" << std::endl;
    std::cout << "Total refining time: " << timer.time()+elapsed_time << "s" << std::endl;
#endif

    timer.stop();
    elapsed_time += timer.time();
    return elapsed_time;
  }

  // Step-by-step methods
  void initialize()
  {
    m_facet_mesher.initialize();
  }

  void one_step()
  {
    if(!m_facet_mesher.is_active())
      m_facet_mesher.initialize();

    if(!m_facet_mesher.is_algorithm_done())
      m_facet_mesher.one_step(m_facet_visitor);
    else
    {
      if(!m_cell_mesher.is_active())
      {
        m_facet_visitor.is_active() = true;
        m_facet_mesher.pick_valid_uses_3D_checks() = true;

        m_cell_mesher.initialize();
      }
      m_cell_mesher.one_step(m_cell_visitor);
    }
  }

  bool is_algorithm_done()
  {
    return m_cell_mesher.is_algorithm_done();
  }

  void report()
  {
    m_cell_mesher.report();
    m_facet_mesher.report();
    std::cout << "consistency of EVERYTHING: ";
    std::cout << m_starset.is_consistent(true /*verbose*/) << std::endl;
  }

public:
  Anisotropic_mesher_3(Starset& starset_,
                       const Constrain_surface* pconstrain_,
                       const Criteria* criteria_,
                       const Metric_field* metric_field_)
    :
      m_starset(starset_),
      m_aabb_tree(100/*insertion buffer size*/),
      m_ch_triangulation(),
      m_kd_tree(m_starset.star_vector()),
      m_star_czones(m_starset),
      m_null_mesher(),
      m_facet_mesher(m_null_mesher, m_starset, pconstrain_, criteria_, metric_field_,
                      m_ch_triangulation, m_aabb_tree, m_kd_tree, m_star_czones),
      m_cell_mesher(m_facet_mesher, m_starset, pconstrain_, criteria_, metric_field_,
                     m_ch_triangulation, m_aabb_tree, m_kd_tree, m_star_czones),
      m_null_visitor(),
      m_facet_visitor(m_cell_mesher, m_null_visitor),
      m_cell_visitor(m_facet_visitor)
  { }

private:
  Anisotropic_mesher_3(const Self& src);
  Self& operator=(const Self& src);
}; // Anisotropic_mesher_3

}  // Anisotropic_mesh_3
}  // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_MESHER_H

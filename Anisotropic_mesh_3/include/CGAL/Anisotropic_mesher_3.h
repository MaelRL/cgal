#ifndef CGAL_ANISOTROPIC_MESH_3_MESHER_H
#define CGAL_ANISOTROPIC_MESH_3_MESHER_H

#include <CGAL/Stretched_Delaunay_3.h>
#include <CGAL/Constrain_surface_3.h>
#include <CGAL/Metric_field.h>
#include <CGAL/Criteria.h>

#include <CGAL/Anisotropic_refine_facets_3.h>
#include <CGAL/Anisotropic_refine_cells_3.h>

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

// Self
  typedef Anisotropic_mesher_3<K, RefinementCondition>      Self;

// Stars
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

// Filters
  typedef CGAL::AABB_tree_bbox<K, Star>                     AABB_tree;
  typedef CGAL::AABB_bbox_primitive<Star>                   AABB_primitive;

  typedef CGAL::Kd_tree_for_star_set<K, Star_handle>        Kd_tree;
  typedef typename Kd_tree::Traits                          Kd_traits;
  typedef typename Kd_tree::Box_query                       Kd_Box_query;
  typedef typename Kd_tree::key_type                        Kd_point_info;

  // Facets mesher level
  typedef CGAL::Anisotropic_mesh_3::Anisotropic_refine_facets_3<
      K,
      Null_mesher_level>                                    Facets_level;

  // Cells mesher level
  typedef CGAL::Anisotropic_mesh_3::Anisotropic_refine_cells_3<
      K,
      Facets_level>                                         Cells_level;

private:
  // Star set
  Star_vector& m_stars;

  // Filters
  AABB_tree m_aabb_tree; //bboxes of stars
  DT m_ch_triangulation;
  Kd_tree m_kd_tree;     //stars* centers for box queries

  // Meshers
  Null_mesher_level m_null_mesher;
  Facets_level m_facet_mesher;
  Cells_level m_cell_mesher;

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
    m_facet_mesher.refine();

    // Then scan volume and refine it
    m_cell_mesher.initialize();
    m_cell_mesher.refine();
#else
    std::cout << "Start surface scan...";
    m_facet_mesher.initialize();
    std::cout << std::endl << "end scan" << std::endl;

    std::ofstream out("initial.mesh");
    output_medit(m_stars, out);

    elapsed_time += timer.time();
    timer.stop(); timer.reset(); timer.start();

    while (!m_facet_mesher.is_algorithm_done())
      m_facet_mesher.one_step();

    std::cout << "Total refining surface time: " << timer.time() << "s" << std::endl;

    elapsed_time += timer.time();
    timer.stop(); timer.reset(); timer.start();

    std::cout << "Start volume scan...";
    m_cell_mesher.initialize();
    std::cout << std::endl << "end scan" << std::endl;

    elapsed_time += timer.time();
    timer.stop(); timer.reset(); timer.start();

    while (!m_cell_mesher.is_algorithm_done())
      m_cell_mesher.one_step();

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
      m_facet_mesher.one_step();
    else
    {
      if(!m_cell_mesher.is_active())
        m_cell_mesher.initialize();
      m_cell_mesher.one_step();
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
    std::cout << is_consistent(m_stars, true /*verbose*/) << std::endl;
  }

public:
  Anisotropic_mesher_3(Star_vector& stars,
                       const Constrain_surface* pconstrain_,
                       const Criteria* criteria_,
                       const Metric_field* metric_field_)
    :
      m_stars(stars),
      m_aabb_tree(100/*insertion buffer size*/),
      m_ch_triangulation(),
      m_kd_tree(m_stars),
      m_null_mesher(),
      m_facet_mesher(m_null_mesher, m_stars, pconstrain_, criteria_, metric_field_,
                      m_ch_triangulation, m_aabb_tree, m_kd_tree),
      m_cell_mesher(m_facet_mesher, m_stars, pconstrain_, criteria_, metric_field_,
                     m_ch_triangulation, m_aabb_tree, m_kd_tree)
  { }

  ~Anisotropic_mesher_3() {}

private:
  Anisotropic_mesher_3(const Self& src);
  Self& operator=(const Self& src);
}; // Anisotropic_mesher_3

}  // Anisotropic_mesh_3
}  // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_MESHER_H

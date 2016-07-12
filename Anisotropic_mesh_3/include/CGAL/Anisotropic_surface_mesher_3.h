#ifndef CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_SURFACE_MESHER_3_H
#define CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_SURFACE_MESHER_3_H

#include <CGAL/Starset.h>
#include <CGAL/Stretched_Delaunay_3.h>
#include <CGAL/Constrain_surface_3.h>
#include <CGAL/Metric_field.h>
#include <CGAL/Criteria.h>

#include <CGAL/Anisotropic_mesher_3_base.h>
#include <CGAL/Anisotropic_refine_facets_3.h>
#include <CGAL/Anisotropic_mesher_visitor.h>

#include <CGAL/IO/Star_set_IO.h>

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

template<typename K, typename RefinementCondition = No_condition<typename K::Point_3> >
class Anisotropic_surface_mesher_3 : public Anisotropic_mesher_3_base
{
  typedef Anisotropic_surface_mesher_3<K, RefinementCondition>  Self;
  typedef Anisotropic_mesher_3_base                             Base;

public:
  //typedef CGAL::Exact_predicates_exact_constructions_kernel KExact;
  typedef K                                                 KExact;

  typedef Stretched_Delaunay_3<K, KExact>                   Star;
  typedef Star*                                             Star_handle;

  typedef Constrain_surface_3<K>                            Constrain_surface;
  typedef CGAL::Anisotropic_mesh_3::Metric_field<K>         Metric_field;
  typedef typename Metric_field::Metric                     Metric;
  typedef Criteria_base<K>                                  Criteria;

  typedef CGAL::Anisotropic_mesh_3::Starset<K>              Starset;

  //Queues
  typedef CGAL::Anisotropic_mesh_3::Facet_refine_queue<K>   Facet_refine_queue;

  //Filters
  typedef CGAL::AABB_tree_bbox<K, Star>                     AABB_tree;
  typedef CGAL::Kd_tree_for_star_set<K, Star_handle>        Kd_tree;
  typedef CGAL::Anisotropic_mesh_3::Stars_conflict_zones<K> Stars_conflict_zones;

  //Mesher levels
  typedef Null_anisotropic_mesher_level                     Null_mesher_level;
  typedef Anisotropic_refine_facets_3<K, Null_mesher_level> Facets_level;

  //Refinement trunk
  typedef Anisotropic_refine_trunk<K>                       Trunk;

  //Visitors
  typedef Null_anisotropic_mesher_visitor                   Null_mesher_visitor;
  typedef Null_anisotropic_mesher_visitor_level<Null_mesher_visitor>
                                                            Facets_visitor;

private:
  // Star set
  Starset& m_starset;

  // Queues
  Facet_refine_queue m_facet_refine_queue;

  // Filters
  AABB_tree m_aabb_tree; //bboxes of stars
  Kd_tree m_kd_tree;     //stars* centers for box queries

  mutable Stars_conflict_zones m_star_czones; //conflict zones for the stars in conflict

  // Meshers
  Null_mesher_level m_null_mesher;
  Facets_level m_facet_mesher; // deals with all the rules of the facet level

  // Visitors
  Null_mesher_visitor m_null_visitor;
  Facets_visitor m_facet_visitor;

  // Refinement restrictions todo
  RefinementCondition m_refinement_condition;

public:
  double refine_mesh()
  {
    CGAL::Timer timer;
    timer.start();
    double elapsed_time = 0.;

#ifndef ANISO_VERBOSE
    // Scan surface and refine it
    m_facet_mesher.initialize();
    m_facet_mesher.refine(m_facet_visitor);
#else
    m_facet_mesher.initialize();

    elapsed_time += timer.time();
    timer.stop(); timer.reset(); timer.start();

    std::ofstream time_out("time.txt");
    while (!m_facet_mesher.is_algorithm_done())
    {
      m_facet_mesher.one_step(m_facet_visitor);
      if(m_starset.size()%10 == 0)
        time_out << m_starset.size() << " " << elapsed_time+timer.time() << std::endl;
    }

    std::cout << "Total refining time: " << timer.time() << "s" << std::endl;
    elapsed_time += timer.time();
    timer.stop(); timer.reset(); timer.start();
#endif

    timer.stop();
    elapsed_time += timer.time();
    return elapsed_time;
  }

  void resume_from_mesh_file(const char* filename)
  {
    m_facet_mesher.resume_from_mesh_file(filename);
  }

  void resume_from_dump_file(const char* filename)
  {
    m_facet_mesher.resume_from_dump_file(filename);
  }


  // Step-by-step methods
  void initialize()
  {
    m_facet_mesher.initialize();
  }

  void one_step()
  {
    m_facet_mesher.one_step(m_facet_visitor);
  }

  bool is_algorithm_done()
  {
    return m_facet_mesher.is_algorithm_done();
  }

  void report()
  {
    std::cout << m_starset.size() << " vertices" << std::endl;
    m_facet_mesher.report();
    std::cout << "consistency of EVERYTHING: ";
    std::cout << m_starset.is_consistent(true /*verbose*/) << std::endl;
  }

public:
  Anisotropic_surface_mesher_3(Starset& starset_,
                               const Constrain_surface* pconstrain_,
                               const Criteria* criteria_,
                               const Metric_field* metric_field_)
    :
      Base(),
      m_starset(starset_),
      m_facet_refine_queue(),
      m_aabb_tree(100/*insertion buffer size*/),
      m_kd_tree(m_starset.star_vector()),
      m_star_czones(m_starset),
      m_null_mesher(),
      m_facet_mesher(m_null_mesher, m_starset, pconstrain_, criteria_, metric_field_,
                     m_aabb_tree, m_kd_tree, m_star_czones,
                     m_facet_refine_queue, 0, 6, true/*use poles*/),
      m_null_visitor(),
      m_facet_visitor(m_null_visitor)
  { }

private:
  Anisotropic_surface_mesher_3(const Self& src);
  Self& operator=(const Self& src);
}; // Anisotropic_surface_mesher_3

}  // Anisotropic_mesh_3
}  // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_SURFACE_MESHER_3_H

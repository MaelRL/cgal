#ifndef CGAL_ANISOTROPIC_MESH_2_MESHER_H
#define CGAL_ANISOTROPIC_MESH_2_MESHER_H

#include <CGAL/Starset.h>
#include <CGAL/Stretched_Delaunay_2.h>
#include <CGAL/Domain_2.h>
#include <CGAL/Metric_field.h>
#include <CGAL/Criteria.h>

#include <CGAL/Anisotropic_mesher_2_base.h>
#include <CGAL/Anisotropic_refine_faces_2.h>
#include <CGAL/Anisotropic_mesher_visitor.h>
#include <CGAL/bbox.h>

#include <CGAL/Kd_tree_for_star_set.h>
#include <CGAL/aabb_tree/aabb_tree_bbox.h>
#include <CGAL/aabb_tree/aabb_tree_bbox_primitive.h>

#include <CGAL/IO/Star_set_output.h>

#include <CGAL/Timer.h>
#include <CGAL/iterator.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>

#include <boost/assign/list_of.hpp> // for 'list_of()'

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
namespace Anisotropic_mesh_2
{

template<typename K, typename RefinementCondition = No_condition<typename K::Point_2> >
class Anisotropic_mesher_2 : public Anisotropic_mesher_2_base
{
  typedef Anisotropic_mesher_2<K, RefinementCondition>      Self;
  typedef Anisotropic_mesher_2_base                         Base;

public:
  //typedef CGAL::Exact_predicates_exact_constructions_kernel KExact;
  typedef K                                                 KExact;

  typedef Stretched_Delaunay_2<K, KExact>                   Star;
  typedef Star*                                             Star_handle;

  typedef Domain_2<K>                                       Domain;
  typedef CGAL::Anisotropic_mesh_2::Metric_field<K>         Metric_field;
  typedef typename Metric_field::Metric                     Metric;
  typedef Criteria_base<K>                                  Criteria;

  typedef CGAL::Anisotropic_mesh_2::Starset<K>              Starset;

  //Queues
  typedef CGAL::Anisotropic_mesh_2::Face_refine_queue<K>    Face_refine_queue;

  //Filters
  typedef CGAL::AABB_tree_bbox<K, Star>                     AABB_tree;
  typedef CGAL::Kd_tree_for_star_set<K, Star_handle>        Kd_tree;
  typedef CGAL::Anisotropic_mesh_2::Stars_conflict_zones<K> Stars_conflict_zones;

  //Mesher levels
  typedef Null_anisotropic_mesher_level                     Null_mesher_level;
  typedef Anisotropic_refine_faces_2<K, Null_mesher_level>  Faces_level;
  typedef Anisotropic_refine_faces_2<K, Faces_level>        Faces_consistency_level;

  //Refinement trunk
  typedef Anisotropic_refine_trunk<K>                       Trunk;

  //Visitors
  typedef Null_anisotropic_mesher_visitor                         Null_mesher_visitor;
  typedef Anisotropic_mesher_visitor<Trunk, Null_mesher_visitor>  Faces_visitor;
  typedef Null_anisotropic_mesher_visitor_level<Faces_visitor>    Faces_consistency_visitor;

private:
  // Star set
  Starset& m_starset;

  // Queues
  Face_refine_queue m_face_refine_queue;

  // Filters
  AABB_tree m_aabb_tree; // bboxes of stars
  Kd_tree m_kd_tree; // stars* centers for box queries

  mutable Stars_conflict_zones m_star_czones; //conflict zones for the stars in conflict

  // Meshers
  Null_mesher_level m_null_mesher;
  Faces_level m_face_mesher; // deals with rules 0-3 of the face level
  Faces_consistency_level m_face_consistency_mesher; // deals with rules 4-5

  // Visitors
  Null_mesher_visitor m_null_visitor;
  Faces_visitor m_face_visitor;
  Faces_consistency_visitor m_face_consistency_visitor;

  // Refinement restrictions todo
  RefinementCondition m_refinement_condition;

public:
  double refine_mesh()
  {
    CGAL::Timer timer;
    timer.start();
    double elapsed_time = 0.;

#if 0//ndef ANISO_VERBOSE
    // Scan surface and refine it
    m_face_mesher.initialize();
    m_face_mesher.refine(m_face_visitor);

    // Then scan and solve face inconsistencies
#ifndef ANISO_NO_CONSISTENCY
    m_face_consistency_mesher.initialize();
    m_face_consistency_mesher.refine(m_face_consistency_visitor);
#endif
#else
    m_face_mesher.initialize();

    elapsed_time += timer.time();
    timer.stop(); timer.reset(); timer.start();

    std::ofstream time_out("time.txt");
    while (!m_face_mesher.is_algorithm_done())
    {
      m_face_mesher.one_step(m_face_visitor);
      if(m_starset.size()%10 == 0)
        time_out << m_starset.size() << " " << elapsed_time+timer.time() << '\n';
    }

    std::cout << "Total refining surface time: " << timer.time() << "s" << std::endl;
    elapsed_time += timer.time();
    timer.stop(); timer.reset(); timer.start();

    // ------------------------------------------------
    std::cout << "Start consistency surface scan...";
    m_face_consistency_mesher.initialize();
    std::cout << std::endl << "end scan" << std::endl;

    elapsed_time += timer.time();
    timer.stop(); timer.reset(); timer.start();

    while (!m_face_consistency_mesher.is_algorithm_done())
    {
      m_face_consistency_mesher.one_step(m_face_consistency_visitor);
      if(m_starset.size()%10 == 0)
        time_out << m_starset.size() << " " << elapsed_time+timer.time() << '\n';
    }

    std::cout << "Total refining consistency surface time: " << timer.time() << "s" << std::endl;
    elapsed_time += timer.time();
    timer.stop(); timer.reset(); timer.start();

    // ------------------------------------------------
    std::cout << "Total refining volume time: " << timer.time() << "s" << std::endl;
    std::cout << "Total refining time: " << timer.time()+elapsed_time << "s" << std::endl;
#endif

    report();
    timer.stop();
    elapsed_time += timer.time();
    time_out << std::endl;
    return elapsed_time;
  }

  void resume_from_mesh_file(const char* filename)
  {
#ifndef ANISO_NO_CONSISTENCY
    m_face_consistency_mesher.resume_from_mesh_file(filename);
#else
    m_face_mesher.resume_from_mesh_file(filename);
#endif
  }

  // Step-by-step methods
  void initialize()
  {
    m_face_mesher.initialize();
  }

  void one_step()
  {
#ifndef ANISO_NO_CONSISTENCY
    m_face_consistency_mesher.one_step(m_face_consistency_visitor);
#else
    m_face_mesher.one_step(m_face_visitor);
#endif
  }

  bool is_algorithm_done()
  {
#ifndef ANISO_NO_CONSISTENCY
    return m_face_consistency_mesher.is_algorithm_done();
#else
    return m_face_mesher.is_algorithm_done();
#endif
  }

  void report()
  {
    std::cout << m_starset.size() << " vertices" << std::endl;
    m_face_mesher.report();
#ifndef ANISO_NO_CONSISTENCY
    m_face_consistency_mesher.report();
#endif
    std::cout << "consistency of EVERYTHING: ";
    std::cout << m_starset.is_consistent(true /*verbose*/) << std::endl;
  }

public:
  Anisotropic_mesher_2(Starset& starset_,
                       const Domain* pdomain_,
                       const Criteria* criteria_,
                       const Metric_field* metric_field_)
    :
      Base(),
      m_starset(starset_),
      m_face_refine_queue(),
      m_aabb_tree(100/*insertion buffer size*/),
      m_kd_tree(m_starset.star_vector()),
      m_star_czones(m_starset),
      m_null_mesher(),
      m_face_mesher(m_null_mesher, m_starset, pdomain_, criteria_, metric_field_,
                    m_aabb_tree, m_kd_tree, m_star_czones, m_face_refine_queue, 0, 3),
      m_face_consistency_mesher(m_face_mesher, m_starset, pdomain_, criteria_,
                                metric_field_, m_aabb_tree, m_kd_tree, m_star_czones,
                                m_face_refine_queue, 4, 5),
      m_null_visitor(),
      m_face_visitor(boost::assign::list_of((Trunk*) &m_face_consistency_mesher), m_null_visitor),
      m_face_consistency_visitor(m_face_visitor)
  { }

private:
  Anisotropic_mesher_2(const Self& src);
  Self& operator=(const Self& src);
}; // Anisotropic_mesher_2

}  // Anisotropic_mesh_2
}  // CGAL

#endif // CGAL_ANISOTROPIC_MESH_2_MESHER_H

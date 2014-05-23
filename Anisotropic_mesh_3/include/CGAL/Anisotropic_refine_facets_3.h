#ifndef CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_REFINE_FACETS_H
#define CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_REFINE_FACETS_H

#include <set>
#include <vector>

#include <CGAL/Timer.h>

#include <CGAL/Facet_refine_queue.h>
#include <CGAL/Stretched_Delaunay_3.h>
#include <CGAL/Anisotropic_mesher_level.h>
#include <CGAL/Anisotropic_refine_trunk.h>

#include <CGAL/helpers/histogram_helper.h>
#include <CGAL/helpers/combinatorics_helper.h>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K, typename Previous_lvl>
class Anisotropic_refine_facets_3 :
    public Anisotropic_mesher_level<Stretched_Delaunay_3<K>,
                                    Anisotropic_refine_facets_3<K, Previous_lvl>,
                                    Previous_lvl>,
    public Anisotropic_refine_trunk<K>
{
private:
  typedef Anisotropic_refine_facets_3<K, Previous_lvl>             Self;
public:
  typedef Anisotropic_mesher_level<Stretched_Delaunay_3<K>,
                                   Anisotropic_refine_facets_3<K, Previous_lvl>,
                                   Previous_lvl>                   Mesher_lvl;
  typedef Anisotropic_refine_trunk<K>                              Trunk;

  typedef Stretched_Delaunay_3<K>                                  Star;
  typedef typename Star::FT                                        FT;
  typedef typename Star::Base                                      DT; // DT_3 with vertex_base_with_info
  typedef Star*                                                    Star_handle;
  typedef std::vector<Star_handle>                                 Star_vector;
  typedef typename Star_vector::iterator                           Star_iterator;
  typedef std::set<Star_handle>                                    Star_set;
  typedef typename Star::Index                                     Index;
  typedef std::set<Index>                                          Index_set;
  typedef typename Star::Point_3                                   Point_3;
  typedef typename Star::TPoint_3                                  TPoint_3;
  typedef std::set<Point_3>                                        Point_set;
  typedef typename Star::Vertex_handle                             Vertex_handle;
  typedef typename Star::Cell_handle                               Cell_handle;
  typedef typename Star::Facet_vector                              Facet_vector;
  typedef typename Star::Facet_handle                              Facet_handle;
  typedef typename Star::Facet_set_iterator                        Facet_set_iterator;
  typedef typename Star::Cell_handle_handle                        Cell_handle_handle;
  typedef typename Star::Facet                                     Facet;
  typedef typename Star::Edge                                      Edge;
  typedef typename Star::Vector_3                                  Vector_3;
  typedef typename Star::Constrain_surface                         Constrain_surface;
  typedef typename Star::Criteria                                  Criteria;

  typedef CGAL::Anisotropic_mesh_3::Starset<K>                     Starset;

  typedef typename CGAL::Anisotropic_mesh_3::Metric_field<K>       Metric_field;
  typedef typename Metric_field::Metric                            Metric;

  typedef CGAL::AABB_tree_bbox<K, Star>                            AABB_tree;
  typedef CGAL::AABB_bbox_primitive<Star>                          AABB_primitive;

  typedef CGAL::Kd_tree_for_star_set<K, Star_handle>               Kd_tree;
  typedef typename Kd_tree::Traits                                 Kd_traits;
  typedef typename Kd_tree::Box_query                              Kd_Box_query;
  typedef typename Kd_tree::key_type                               Kd_point_info;

  typedef CGAL::Anisotropic_mesh_3::Conflict_zone<K>               Conflict_zone;
  typedef CGAL::Anisotropic_mesh_3::Stars_conflict_zones<K>        Stars_conflict_zones;

  typedef CGAL::Anisotropic_mesh_3::Facet_refine_queue<K>          Refine_queue;
  typedef typename Refine_queue::Rfacet                            Refine_facet;
  typedef typename Refine_queue::Rfacet_set_iterator               Rfacet_set_iterator;

private:
  Refine_queue& m_refine_queue;

  // these two ints determine which queues (in the facet refinement queue) are
  // used by this level
  const int m_queue_ids_start;
  const int m_queue_ids_end;

  //this boolean is exactly m_facet_visitor.is_active(), but carrying the visitor
  //all the way down to is_valid_point() is even less pretty
  bool m_pick_valid_uses_3D_checks;

  //debug & info
  int m_pick_valid_succeeded;
  int m_pick_valid_failed;
  int m_pick_valid_skipped;
  int m_pick_valid_rejected;
  int m_pick_valid_max_failures;
  mutable int vertex_with_picking_count;
  mutable int vertex_without_picking_count;
  mutable int m_leak_counter;

public:
  mutable CGAL::Timer timer_pv;
  mutable CGAL::Timer timer_npv;

  bool& pick_valid_uses_3D_checks() { return m_pick_valid_uses_3D_checks; }
  const bool& pick_valid_uses_3D_checks() const { return m_pick_valid_uses_3D_checks; }

  bool& is_active() { return Mesher_lvl::is_active(); }
  const bool& is_active() const { return Mesher_lvl::is_active(); }

  bool is_criterion_tested(int queue_id) const
  {
    return (m_queue_ids_start <= queue_id && queue_id <= m_queue_ids_end);
  }

public:
//functions used in the Mesher_lvl class
  void initialize_()
  {
    std::cout << "Initializing facet level with ids: " << m_queue_ids_start;
    std::cout << " " << m_queue_ids_end << std::endl;

    if(this->m_starset.empty())
    {
      std::cout << "Starting with criteria: " << std::endl;
      this->m_criteria->report();

      Trunk::initialize_stars();
      this->m_ch_triangulation.infinite_vertex()->info() = -10;
      Trunk::build_aabb_tree();
    }

    fill_refinement_queue();
    Mesher_lvl::is_active() = true;

    std::ofstream out("initial.mesh");
    output_medit(this->m_starset, out);
    std::ofstream out_off("initial.off");
    output_off(this->m_starset, out_off);
  }

  bool is_algorithm_done_()
  {
    return m_refine_queue.empty(m_queue_ids_start, m_queue_ids_end);
  }

  Refinement_point_status get_refinement_point_for_next_element_(Point_3& steiner_point)
  {
    Rfacet_set_iterator bad_facet;
    bool need_picking_valid;
    Facet f; // to be refined

    //Top of the prio queue is not popped yet since the ref point could encroach.
    //It'll be popped at insertion if there is no encroachment.
    if(!next_refine_facet(bad_facet, f, need_picking_valid))
      return EMPTY_QUEUE;

#ifdef ANISO_DEBUG_REFINEMENT_PP
    Vertex_handle v1 = f.first->vertex((f.second+1)%4);
    Vertex_handle v2 = f.first->vertex((f.second+2)%4);
    Vertex_handle v3 = f.first->vertex((f.second+3)%4);
    Index ind_v1 = v1->info();
    Index ind_v2 = v2->info();
    Index ind_v3 = v3->info();
    std::cout << "Trying to refine : " << ind_v1 << " " << ind_v2 << " " << ind_v3 << std::endl;
    std::cout << "\tp"<< v1->info() <<" : " << bad_facet->star->metric().inverse_transform(v1->point()) << std::endl;
    std::cout << "\tp"<< v2->info() <<" : " << bad_facet->star->metric().inverse_transform(v2->point()) << std::endl;
    std::cout << "\tp"<< v3->info() <<" : " << bad_facet->star->metric().inverse_transform(v3->point()) << std::endl;
    std::cout << "\tcheck coordinates : from stars, " << std::endl;
    std::cout << "\t" << m_stars[ind_v1]->center_point() << " " << &(m_stars[ind_v1]) << std::endl;
    std::cout << "\t" << m_stars[ind_v2]->center_point() << " " << &(m_stars[ind_v2]) << std::endl;
    std::cout << "\t" << m_stars[ind_v3]->center_point() << " " << &(m_stars[ind_v3]) << std::endl;

    double deg_value = 1e-4;
    bool degenerate = ( (std::abs(v1->point().x()-v2->point().x()) < deg_value &&
                         std::abs(v1->point().y()-v2->point().y()) < deg_value &&
                         std::abs(v1->point().z()-v2->point().z()) < deg_value ) ||
                        (std::abs(v2->point().x()-v3->point().x()) < deg_value &&
                         std::abs(v2->point().y()-v3->point().y()) < deg_value &&
                         std::abs(v2->point().z()-v3->point().z()) < deg_value ) ||
                        (std::abs(v1->point().x()-v3->point().x()) < deg_value &&
                         std::abs(v1->point().y()-v3->point().y()) < deg_value &&
                         std::abs(v1->point().z()-v3->point().z()) < deg_value) );

    if(degenerate)
      std::cout << "trying to refine a degenerate (in the metric) facet" << std::endl;
#endif

    if(!this->m_criteria->max_times_to_try_in_picking_region || //skip pick_valid if the number of tries is set to 0
       (need_picking_valid &&
        this->m_starset.compute_distortion(f) > this->m_criteria->distortion)) //pick_valid trick #1: skip pick_valid if the distortion is too high
    {
      m_pick_valid_skipped++;
      need_picking_valid = false;
    }

    // note: failure in pick_valid AND a conflicting circumcenter gives POINT_IN_CONFLICT status (trick#2 is not applied)
    Refinement_point_status rp_status = compute_steiner_point(bad_facet->star, f,
                                                              need_picking_valid, steiner_point);

    if(need_picking_valid)
      pick_valid_output(rp_status); //counts point_in_conflict as fail (not sure if should)

    if(rp_status == POINT_IN_CONFLICT)
      return rp_status;

    // pick_valid trick #2: If an element fails a pick_valid test, put it at the end
    // of the (same) queue in hope that the (succesful) refinement of another element
    // will also solve the problem for the rejected element.
    if(rp_status == PICK_VALID_FAILED &&
       bad_facet->value != m_refine_queue.queue_min_value(bad_facet->queue_type) && //nothing to do if already last
       !bad_facet->prev_rejection) // only allow one rejection
    {
      Trunk::clear_conflict_zones();
      timer_npv.start();
      m_pick_valid_rejected++;
      m_refine_queue.reject_rfacet(bad_facet->queue_type);
      timer_npv.stop();
      return get_refinement_point_for_next_element_(steiner_point);
    }

    //We already know the conflict zones if it's a suitable point from pick_valid, but we need
    //to compute them for the other cases (won't cost anything if it's already known).
    Trunk::compute_conflict_zones(steiner_point);
    //same for elements needing checks
    this->m_stars_czones.compute_elements_needing_check();

    return SUITABLE_POINT;
  }

  bool test_point_conflict_from_superior_(const Point_3& p,
                                          const bool is_queue_updated = true,
                                          const bool need_picking_valid = false)
  {
    //return false; //if we want to ignore the surface level

      //this creates entries in the czones map
    this->m_stars_czones.conflicting_point_id() =  Trunk::simulate_insert_to_stars(p);
    typename Stars_conflict_zones::iterator czit = this->m_stars_czones.begin();
    typename Stars_conflict_zones::iterator czend = this->m_stars_czones.end();
    for(; czit!=czend; ++czit)
    {
      Index i = czit->first;
      Star_handle star = Trunk::get_star(i);

      Facet f;
      if(star->is_encroached(p, f))
      {
/*
        std::cout << "conflict" << std::endl;
        std::cout << p << " encroaches " << f.first->vertex((f.second+1)%4)->info();
        std::cout << " " << f.first->vertex((f.second+2)%4)->info() << " ";
        std::cout << f.first->vertex((f.second+3)%4)->info() << " up: ";
        std::cout << is_queue_updated << std::endl;
*/
        if(is_queue_updated)
        {
          //there can only be one encroached facet in the queue at a time
          //so we use the value to pass the boolean
          m_refine_queue.push(star, f, need_picking_valid, 0);
        }
        return true;
      }
    }
    return false;
  }

  bool insert_(const Point_3& steiner_point)
  {
    Refine_facet bad_facet;
    bool need_picking_valid;
    Facet f; // to be refined

    //At that point, recovering the bad_facet and everything else is just useful for debugging
    //Only queue's pop() is truly needed
    if(!next_refine_facet_pop(bad_facet, f, need_picking_valid))
      return false; //should never happen

    Vertex_handle v1 = f.first->vertex((f.second+1)%4);
    Vertex_handle v2 = f.first->vertex((f.second+2)%4);
    Vertex_handle v3 = f.first->vertex((f.second+3)%4);
    Index ind_v1 = v1->info();
    Index ind_v2 = v2->info();
    Index ind_v3 = v3->info();

    //todo re-add those conditions, (need to template the mesher lvl with it)
    //if(!m_refinement_condition(steiner_point))
    //  return true; //false would stop refinement

    Index pid = Trunk::insert(steiner_point, true/*conditional*/, true/*surface point*/);

    if(pid != static_cast<Index>(this->m_starset.size()-1))
      std::cout << "warning in insert_" << std::endl;

    //Debug: check if f has been destroyed -------------------------------------
    Cell_handle c;
    int i,j,k;
    if(bad_facet.star->is_facet(v1, v2, v3, c, i, j, k))
    {
      int index = 6 - i - j - k;
      Facet ff = bad_facet.star->make_canonical(Facet(c,index));

      typename K::Plane_3 fplane = bad_facet.star->triangle(ff).supporting_plane();

      Point_3 steiner2;
      if(bad_facet.star->is_restricted(ff, steiner2)
        && fplane.oriented_side(bad_facet.star->metric().transform(steiner_point))
           == fplane.oriented_side(bad_facet.star->metric().transform(steiner2)))
      {
        //if steiner_point and steiner2 are not on the same side from ff, it's legal

        Trunk::clear_conflict_zones();
        Trunk::pop_back_star();

        std::cout << "found & restricted" << std::endl;
        if(bad_facet.star->is_facet(v1, v2, v3, c, i, j, k))
        {
          index = 6 - i - j - k;
          ff = bad_facet.star->make_canonical(Facet(c,index));

          Trunk::compute_conflict_zones(steiner2);
          compute_exact_steiner_point(bad_facet.star, ff, need_picking_valid, steiner2);
          pid = Trunk::insert(steiner2, true/*conditional*/, true/*surface point*/);
        }

#if 1//def ANISO_DEBUG_REFINEMENT
        if(bad_facet.star->is_facet(v1, v2, v3, c, i, j, k))
        {
          index = 6 - i - j - k;
          ff = bad_facet.star->make_canonical(Facet(c,index));
          if(bad_facet.star->is_restricted(ff))
          {
            std::cout.precision(15);
            std::cout << "Bad facet still there. Bad_facet.star : " << bad_facet.star->index_in_star_set() << ". stars : " << this->m_starset.size() << std::endl;
            std::cout << "star's metric" << std::endl << bad_facet.star->metric().get_transformation() << std::endl;
            std::cout << "evs : " << bad_facet.star->metric().get_max_eigenvalue() << " ";
            std::cout << bad_facet.star->metric().get_min_eigenvalue() << " ";
            std::cout << bad_facet.star->metric().get_third_eigenvalue() << std::endl;
            std::cout << "\tp"<< v1->info() <<" : " << bad_facet.star->metric().inverse_transform(v1->point()) << std::endl;
            std::cout << "check : " << this->m_starset[v1->info()]->center_point() << std::endl;
            std::cout << "\tp"<< v2->info() <<" : " << bad_facet.star->metric().inverse_transform(v2->point()) << std::endl;
            std::cout << "check : " << this->m_starset[v2->info()]->center_point() << std::endl;
            std::cout << "\tp"<< v3->info() <<" : " << bad_facet.star->metric().inverse_transform(v3->point()) << std::endl;
            std::cout << "check : " << this->m_starset[v3->info()]->center_point() << std::endl;

            std::cout << "\tdim = " << bad_facet.star->dimension();
            std::cout << ", nbv = " << bad_facet.star->number_of_vertices();
            std::cout << ", pid = " << pid;
            std::cout << ", f(p) = " << this->m_pConstrain->side_of_constraint(steiner_point);
            std::cout << ", picking : " << need_picking_valid << std::endl;
            std::cout << "\tp = " << steiner_point  << " " << this->m_pConstrain->side_of_constraint(steiner_point) << std::endl;

            Cell_handle useless;
            std::cout << "checking conflicts with all three stars : ";
            std::cout << this->m_starset[ind_v1]->is_conflicted(this->transform_to_star_point(steiner_point, this->m_starset[ind_v1]), useless) << " ";
            std::cout << this->m_starset[ind_v2]->is_conflicted(this->transform_to_star_point(steiner_point, this->m_starset[ind_v2]), useless) << " ";
            std::cout << this->m_starset[ind_v3]->is_conflicted(this->transform_to_star_point(steiner_point, this->m_starset[ind_v3]), useless) << std::endl;

            std::cout << "checking needs_aabb_update in the tree of size " << this->m_aabb_tree.size() << " : ";
            std::cout << this->m_starset[ind_v1]->bbox_needs_aabb_update() << " ";
            std::cout << this->m_starset[ind_v2]->bbox_needs_aabb_update() << " ";
            std::cout << this->m_starset[ind_v3]->bbox_needs_aabb_update() << std::endl;

            bad_facet.star->debug_steiner_point(steiner_point, ff, true);

            //pop_back_star();
            return false;
          }
          else
            std::cout << "bad facet still in but not restricted anymore" << std::endl;
        }
        else
          std::cout << "Switching to exact succeeded" << std::endl;
#endif
      }
    }
    // --------------------------------------------------------------------------------------------------

    if(this->m_starset.size()%100 == 0) // TMP
    {
      std::ofstream out_med("bambimboum_wip.mesh");
      output_surface_medit(this->m_starset, out_med);
    }
    return true;
  }
//End of CRTP functions

  void report()
  {
    std::cout << "facet consistency : ";
    std::cout << this->m_starset.is_consistent(true /*verbose*/, FACETS_ONLY) << std::endl;
    all_facet_histograms(this->m_starset, this->m_pConstrain, this->m_criteria);
    std::cout << "FACET pick_valid stats: " << std::endl;
    std::cout << "tried: " << m_pick_valid_points_tried << " || ";
    std::cout << "skipped: " << m_pick_valid_skipped << " || ";
    std::cout << "succeeded: " << m_pick_valid_succeeded << " || ";
    std::cout << "rejected: " << m_pick_valid_rejected << " || ";
    std::cout << "failed: " << m_pick_valid_failed << std::endl;
    std::cout << "avg distortion: " << m_avg_dist/((double) (m_pick_valid_skipped+m_pick_valid_succeeded+m_pick_valid_failed)) << std::endl;
    std::cout << "facet pv: " << timer_pv.time() << " " << timer_pv.intervals() << std::endl;
    std::cout << "facet npv: " << timer_npv.time() << " " << timer_npv.intervals() << std::endl;
    std::cout << "facet leaking: " << m_leak_counter << std::endl;
  }

public:
  //Facet refinement function
  bool next_refine_facet(Rfacet_set_iterator& rfacet_it,
                         Facet& facet,
                         bool& need_picking_valid)
  {
    while(true)
    {
      if(!m_refine_queue.top(rfacet_it, m_queue_ids_start, m_queue_ids_end))
      {
        return false;
#if 0
        std::cout << "it says it's empty" << std::endl;
        fill_refinement_queue();
        if(!m_refine_queue.top(rfacet_it, m_queue_ids_start, m_queue_ids_end))
        {
          std::cout << "it ain't lying" << std::endl;
          return false;
        }
        else
        {
          std::cout << "it LIED" << std::endl;
          continue;
        }
#endif
      }

      if(rfacet_it->star->has_facet(rfacet_it->facet, facet))
      {
        need_picking_valid = m_refine_queue.need_picking_valid(rfacet_it->queue_type);

        //for encroachments, value == need_picking_valid
        if(rfacet_it->queue_type == m_refine_queue.encroachment_queue &&
           rfacet_it->value)
          need_picking_valid = true;

        return true;
      }
      else //top of the queue does not exist anymore
        m_refine_queue.pop(m_queue_ids_start, m_queue_ids_end);
    }
  }

  //Most of that function body/parameters is debug; it could simply be f(i){queues[i]->pop();}
  bool next_refine_facet_pop(Refine_facet& refine_facet,
                             Facet& facet,
                             bool& need_picking_valid)
  {
    Rfacet_set_iterator rfsit;
    if(!m_refine_queue.top(rfsit, m_queue_ids_start, m_queue_ids_end))
    {
      std::cout << "empty queue at facet pop time...?" << std::endl; //shouldn't happen
      return false;
    }

    refine_facet = *rfsit; //have to copy since the pop invalidates the iterator

    if(!rfsit->star->has_facet(rfsit->facet, facet))
    {
      std::cout << "problems at facet pop time" << std::endl;
      return false;
    }

    need_picking_valid = m_refine_queue.need_picking_valid(rfsit->queue_type);
    m_refine_queue.pop(m_queue_ids_start, m_queue_ids_end);
    return true;
  }

private:
  template<typename Facet_it>
  void test_facet(Star_handle star, Facet_it fi, // temporary, idea is to eventually break the criteria down like in /Mesh_3 (is_bad)
                  bool force_push = false, bool check_if_in = false)
  {
    //TODO (?) A possibility might be to check whether the facet already exists
    //in the queue (regardless of the queue/value/star). If this is the case, then we move on
    //to the next facet.

    //What is done atm is to update the queue IF HIGHER PRIO (new queue_type < old queue_type or
    // new_value > old_value).

    //Benefits would be the cost of checking the criteria + update cost but the queues would
    //not have the exact correct priority order they should have.

    // note : if a facet is encroached, the queue will be filled from test_conflict_from_superior()
    //        the code below can be used to brute force check that no facet is encroached
    // encroachment : 0
/*
    if(Trunk::is_encroached(star, *fi))
    {
      m_refine_queue.push_(star, *fi, star->compute_volume(*fi), 0, force_push);
      continue;
    }
*/

    // note : distortion is now used only to speed-up pick_valid (see pick_valid trick#1)
    // over distortion : 1
/*
    if(is_criterion_tested(m_refine_queue.over_distortion_queue) &&
       this->m_criteria->distortion > 0.)
    {
      FT over_distortion = this->m_starset.compute_distortion(*fi) - this->m_criteria->distortion;
      if(over_distortion > 0.)
      {
        if(!check_if_in || !m_refine_queue.is_facet_in(star, *fi, over_distortion, 1))
        {
          m_refine_queue.push(star, *fi, over_distortion, 1, force_push);
          if(check_if_in)
          {
            m_leak_counter++;
            std::cout << "not already in 1" << std::endl;
          }
        }
        return;
      }
    }
*/

    // size : 2
    if(is_criterion_tested(m_refine_queue.over_circumradius_queue) &&
       this->m_criteria->facet_circumradius > 0.)
    {
      FT over_circumradius = star->compute_circumradius_overflow(*fi);
      if (over_circumradius > 0)
      {
        if(!check_if_in || !m_refine_queue.is_facet_in(star, *fi, over_circumradius, 2))
        {
          m_refine_queue.push(star, *fi, over_circumradius, 2, force_push);
          if(check_if_in)
          {
            m_leak_counter++;
          }
        }
        return;
      }
    }

    // approx : 3
    if(is_criterion_tested(m_refine_queue.over_approximation_queue) &&
       this->m_criteria->approximation > 0.)
    {
      FT over_approx = Trunk::sq_distance_to_surface(*fi, star) - this->m_criteria->squared_approximation;
      if(over_approx > 0.)
      {
        if(!check_if_in || !m_refine_queue.is_facet_in(star, *fi, over_approx, 3))
        {
          m_refine_queue.push(star, *fi, over_approx, 3, force_push);
          if(check_if_in)
          {
            m_leak_counter++;
          }
        }
        return;
      }
    }

    // shape : 4
    if(is_criterion_tested(m_refine_queue.bad_shape_queue) &&
       this->m_criteria->facet_radius_edge_ratio > 0.)
    {
      FT over_radius_edge_ratio = star->compute_radius_edge_ratio_overflow(*fi);
      if (over_radius_edge_ratio > 0)
      {
        if(!check_if_in || !m_refine_queue.is_facet_in(star, *fi, over_radius_edge_ratio, 4))
        {
          m_refine_queue.push(star, *fi, over_radius_edge_ratio, 4, force_push);
          if(check_if_in)
          {
            std::cout << "not already in 4" << std::endl;
            std::cout << "star : " << star->index_in_star_set();
            std::cout << " || ids: ";
            std::cout << fi->first->vertex((fi->second+1)%4)->info() << " ";
            std::cout << fi->first->vertex((fi->second+2)%4)->info() << " ";
            std::cout << fi->first->vertex((fi->second+3)%4)->info() << std::endl;

            std::cout << "star vertices:" << std::endl;
            star->print_vertices();
          }
        }
        return;
      }
    }

    // consistency : 5
    if(is_criterion_tested(m_refine_queue.inconsistent_queue) &&
       !this->m_starset.is_consistent(*fi))
    {
      FT vol = star->compute_volume(*fi);
      if(!check_if_in || !m_refine_queue.is_facet_in(star, *fi, vol, 5))
      {
        m_refine_queue.push(star, *fi, vol, 5, force_push);
        if(check_if_in)
        {
          std::cout << "not already in 5" << std::endl;
          std::cout << "star : " << star->index_in_star_set();
          std::cout << " || ids: ";
          std::cout << fi->first->vertex((fi->second+1)%4)->info() << " ";
          std::cout << fi->first->vertex((fi->second+2)%4)->info() << " ";
          std::cout << fi->first->vertex((fi->second+3)%4)->info() << std::endl;

          std::cout << "star vertices:" << std::endl;
          star->print_vertices();
        }
      }
    }
  }

public:
  //facets here are necessarily inconsistent, but need to test other criteria too
  void fill_from_unmodified_stars()
  {
    typename std::map<Index, Facet_vector>::iterator mit = this->m_stars_czones.facets_to_check().begin();
    typename std::map<Index, Facet_vector>::iterator mend = this->m_stars_czones.facets_to_check().end();
    for(; mit!=mend; ++mit)
    {
      Index i = mit->first;
      Star_handle si = Trunk::get_star(i);

      Facet_handle fit = mit->second.begin();
      Facet_handle fend = mit->second.end();
      for(; fit!=fend; ++fit)
      {
        test_facet(si, fit, true/*force push*/);
      }
    }
  }

  void fill_refinement_queue(Index pid)
  {
    //fill from the stars that were directly in conflict with the inserted point
    fill_refinement_queue(this->m_stars_czones, pid);

    //fill from unmodified stars, where some facets could become inconsistent (or the
    //star ownership of the facet in the queue should change if the facet was already in another
    //queue).
    fill_from_unmodified_stars();

#if 0
    std::cout << "Enter fill f_ref_queue debug. Filling with all stars" << std::endl;
    fill_refinement_queue();
    m_refine_queue.print();
    std::cout << "End fill f_ref_queue debug" << std::endl;
#endif
    std::cout << this->m_starset.size() << " ";
    m_refine_queue.print();
  }

  void fill_refinement_queue()
  {
    std::cout << "fill from all facets" << std::endl;
    fill_refinement_queue(this->m_starset, -1);
  }

private:
  template<typename Stars>
  void fill_refinement_queue(const Stars& stars,
                             Index relative_point,
                             bool check_if_in = false)
  {
    /*
    std::cout << "fill facet call with relative " << relative_point << std::endl;
    typename Stars::const_iterator cit = stars.begin();
    typename Stars::const_iterator citend = stars.end();
    for (; cit != citend; cit++)
      std::cout << " " << this->get_star(cit)->index_in_star_set();
    std::cout << std::endl;
    */

#ifdef USE_ANISO_TIMERS
    std::clock_t start_time = clock();
#endif
    typename Stars::const_iterator si = stars.begin();
    typename Stars::const_iterator siend = stars.end();
    for (; si != siend; si++)
    {
      Star_handle star = Trunk::get_star(si);

      Facet_set_iterator fi = star->restricted_facets_begin();
      Facet_set_iterator fiend = star->restricted_facets_end();
      for (; fi != fiend; fi++)
      {
/*TODO
        Point_3 cc;
        star->compute_dual_intersection(*fi, cc);
        if(!m_refinement_condition(cc))
          continue;
*/
        Cell_handle cell = fi->first;
        int offset = fi->second;

#ifdef ANISO_DEBUG_REFINEMENT_PP
        //check for absurdities
        Vertex_handle v1 = cell->vertex((offset+1)%4);
        Vertex_handle v2 = cell->vertex((offset+2)%4);
        Vertex_handle v3 = cell->vertex((offset+3)%4);

        typename Star::Traits::Compute_squared_distance_3 csd =
            m_traits->compute_squared_distance_3_object();

        if(csd(star->metric().inverse_transform(v1->point()), m_stars[v1->info()]->center_point()) > 1e-10 ||
           csd(star->metric().inverse_transform(v2->point()), m_stars[v2->info()]->center_point()) > 1e-10 ||
           csd(star->metric().inverse_transform(v3->point()), m_stars[v3->info()]->center_point()) > 1e-10)
        {
          std::cout.precision(20);
          std::cout << "points differ in fill_ref_queue: " << v1->info() << " " << v2->info() << " " << v3->info() << std::endl;
          std::cout << "index in star set : " <<  star->index_in_star_set() << std::endl;
          std::cout << star->metric().inverse_transform(v1->point()) << std::endl;
          std::cout << star->metric().inverse_transform(v2->point()) << std::endl;
          std::cout << star->metric().inverse_transform(v3->point()) << std::endl;
          std::cout << m_stars[v1->info()]->center_point() << std::endl;
          std::cout << m_stars[v2->info()]->center_point() << std::endl;
          std::cout << m_stars[v3->info()]->center_point() << std::endl;
        }
#endif

        if(relative_point >= 0) // we do not consider not-relative facets
        {
          bool relative = false;
          for(int i = 1; i <= 3; i++)
          {
            if(relative_point == cell->vertex((offset + i) % 4)->info())
            {
              relative = true;
              break;
            }
          }
          if(!relative)
            continue;
        }

        test_facet(star, fi, false/*no force push*/, check_if_in);
      }
    }
  }

private:
  //Point computation
  Point_3 compute_random_steiner_point(const Star_handle star,
                                       const Facet& facet,
                                       const TPoint_3& tccf,
                                       const FT& circumradius) const
  {
    typename Star::Traits::Compute_random_point_3 random =
        star->traits()->compute_random_point_3_object();

    bool found_acceptable_steiner_point = false;
    int failures_count = 0;
    Point_3 steiner_point;
    TPoint_3 random_point_within_sphere_1;
    TPoint_3 random_point_within_sphere_2;

    while(!found_acceptable_steiner_point)
    {
      if(++failures_count > 100)
        std::cout << "failures_count is getting big : " << failures_count << std::endl;

      random_point_within_sphere_1 = random(tccf, circumradius);
      random_point_within_sphere_2 = random(tccf, circumradius);

      found_acceptable_steiner_point =
          star->compute_steiner_dual_intersection(steiner_point,
                                                  random_point_within_sphere_1,
                                                  random_point_within_sphere_2,
                                                  tccf,
                                                  facet,
                                                  circumradius,
                                                  failures_count);
    }

    return steiner_point;
  }

  void debug_is_cycle(const std::set<Edge_ij>& edges) const
  {
    std::map<int, int> m;
    typename std::set<Edge_ij>::const_iterator it;
    for(it = edges.begin(); it != edges.end(); it++)
    {
      int i = (*it).vertex(0);
      if(m.count(i) > 0) m.erase(i);
      else m[i] = 1;

      i = (*it).vertex(1);
      if(m.count(i) > 0) m.erase(i);
      else m[i] = 1;
    }
    if(!m.empty())
      std::cerr << "Error : edges_around_c is not a cycle!\n";
  }

  template<typename Edge_ijSet, typename IndexOutputIterator>
  void get_vertices(const Edge_ijSet& edges, IndexOutputIterator oit) const
  {
    typename Edge_ijSet::const_iterator eit;
    for(eit = edges.begin(); eit != edges.end(); eit++)
    {
      *oit++ = (*eit).vertex(0);
      *oit++ = (*eit).vertex(1);
    }
  }

  template<typename FacetVector, typename EdgeSet>
  void get_cycle(const FacetVector& bfacets,
                 EdgeSet& edges_around_c,
                 const int& center_id) const
  {
    typename FacetVector::const_iterator fit = bfacets.begin();
    for( ; fit != bfacets.end(); fit++)
    {
      Facet_ijk f(*fit);
      boost::array<int,2> others;
      if(f.has(center_id, others))
        edges_around_c.insert(Edge_ij(others[0], others[1]));
    }
    debug_is_cycle(edges_around_c);
  }

  void facets_created(Star_handle star,
                      std::map<Facet_ijk, int>& facets) const
  {
    typename Star::Facet_set_iterator it = star->restricted_facets_begin();
    typename Star::Facet_set_iterator iend = star->restricted_facets_end();
    for(; it != iend; it++)
      add_to_map(Facet_ijk(*it), facets);
  }

  // restricted facets(i,j,k) which would be created by the insertion
  // of 'p' in 'star' are added to 'facets'
  void facets_created(const Point_3& p, // new point, not yet in the star set
                      const Index p_id,   // index of p in star set
                      Star_handle star,
                      std::map<Facet_ijk, int>& facets) const
  {
    Index center_id = star->index_in_star_set();
    if(center_id == p_id)
      return facets_created(star, facets);

    // boundary facets of the conflict zone should already have been computed
    if(this->m_stars_czones.status() != Stars_conflict_zones::CONFLICT_ZONES_ARE_KNOWN)
      std::cout << "Zones should be known before entering facets_created..." << std::endl;
    if(this->m_stars_czones.conflict_zones().find(center_id) == this->m_stars_czones.end())
      std::cout << "Trying to compute facets_created for a star not in conflict" << std::endl;

    const Conflict_zone& star_cz = this->m_stars_czones.conflict_zone(center_id);
    const std::vector<Facet>& local_bfacets = star_cz.boundary_facets();

    int dim = star->dimension();
    // dimension 1
    if(dim == 1)
    {
      Edge_ij e(*(star->finite_edges_begin()));
      add_to_map(Facet_ijk(p_id, e.vertex(0), e.vertex(1)), facets);
      return;
    }

    std::set<Edge_ij> edges_around_c;
    if(dim == 2)
    {
      typename Star::Finite_edges_iterator eit = star->finite_edges_begin();
      typename Star::Finite_edges_iterator endit = star->finite_edges_begin();
      for(; eit != endit; eit++)
      {
        Edge_ij e(*eit);
        if(e.has(center_id))
          edges_around_c.insert(e);
      }
    }
    // dimension 3
    else if(dim == 3)
    { // get the circular set of edges/vertices around the center in 'boundary facets'
      get_cycle(local_bfacets, edges_around_c, center_id);
    }
    else std::cerr << "Warning : in 'facets_created', no conflict zone!\n";

    // get the set (circular if dim is 3) of vertices around the center
    std::set<int> vertices_around_c;
    get_vertices(edges_around_c, std::inserter(vertices_around_c, vertices_around_c.end()));

#ifdef ANISO_DEBUG_CONFLICT
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "enter facets_created with center_id!=pid" << std::endl;
    std::cout << "centerid/pid: " << center_id << " " << p_id << std::endl;
    std::cout << "dimension: " << star->dimension() << std::endl;
    std::cout << "get_cycle returns the following BOUNDARY facets:" << std::endl;
    for(std::size_t i=0; i<local_bfacets.size(); ++i)
    {
      Facet f = local_bfacets[i];
      std::cout << f.first->vertex((f.second+1)%4)->info() << " ";
      std::cout << f.first->vertex((f.second+2)%4)->info() << " ";
      std::cout << f.first->vertex((f.second+3)%4)->info() << std::endl;
    }
    std::cout << "edges around c:" << std::endl;
    for(typename std::set<Edge_ij>::iterator it=edges_around_c.begin(); it!=edges_around_c.end(); ++it)
    {
      Edge_ij e = *it;
      std::cout << e[0] << " " << e[1] << std::endl;
    }
    std::cout << "thus vertices around the center are: " << std::endl;
    for(std::set<int>::iterator it=vertices_around_c.begin(); it!=vertices_around_c.end(); ++it)
    {
      std::cout << *it << " ";
    }
    std::cout << std::endl;
#endif

    // facets to be added are : (p, center_i, vertices of 'vertices_around_c')
    // add them only if they are restricted to the constrain surface
    typename std::set<int>::iterator iit;
    for(iit = vertices_around_c.begin(); iit != vertices_around_c.end(); iit++)
    {
      int v_id = *iit;
      if(Trunk::is_infinite_vertex(v_id))
        continue; //it is an infinite vertex

      Edge_ij e1, e2; // edges incident to vertex with 'v_id' in edges_around_c
      bool one_found = false;
      typename std::set<Edge_ij>::iterator eit;
      for(eit = edges_around_c.begin(); eit != edges_around_c.end(); eit++)
      {
        Edge_ij tmp = *eit;
        if(tmp.has(v_id)) //v_id is shared by e1 and e2
        {
          if(!one_found)
          {
            e1 = tmp;
            one_found = true;
          }
          else
          {
            e2 = tmp;
            break;
          }
        }
      }

      //if(!e1.is_valid())
      //  std::cout << "Warning : e1 was not found!" << std::endl;
      if(!e2.is_valid())
      {
        std::cout << "Warning : e2 was not found!" << std::endl;
        continue;//TODO : check this is correct
      }
      bool is_finite_e1 = !e1.is_infinite();
      bool is_finite_e2 = !e2.is_infinite();
      bool is_inside_c1 = false;
      bool is_inside_c2 = false;
      if(is_finite_e1)
      {
        Point_3 cc = Trunk::compute_circumcenter(p, star->center_point(),
                                          this->m_starset[e1.vertex(0)]->center_point(),
                                          this->m_starset[e1.vertex(1)]->center_point(),
                                          star);
        is_inside_c1 = Trunk::is_inside_domain(cc);
      }
      if(is_finite_e2)
      {
        Point_3 cc = Trunk::compute_circumcenter(p, star->center_point(),
                                          this->m_starset[e2.vertex(0)]->center_point(),
                                          this->m_starset[e2.vertex(1)]->center_point(),
            star);
        is_inside_c2 = Trunk::is_inside_domain(cc);
      }

      if((is_finite_e1 && is_finite_e2 && is_inside_c1 != is_inside_c2) //finite case
          || (!is_finite_e1 && is_finite_e2 && is_inside_c2)  // finite+infinite case
          || (!is_finite_e2 && is_finite_e1 && is_inside_c1)) // finite+infinite case
        add_to_map(Facet_ijk(center_id, p_id, v_id), facets);
    }
  }

  bool check_consistency(Star_handle to_be_refined, //facet to be refined belongs to this star
                         Star_handle new_star,     //the newly created star
                         const double& sq_radius_bound) const
  {
    //list all facets that would be created by p's insertion
    Point_3 p = new_star->center_point();
    Index p_index = new_star->index_in_star_set();

    std::map<Facet_ijk, int> facets;
    facets_created(new_star, facets);

    typename Stars_conflict_zones::iterator czit = this->m_stars_czones.begin();
    typename Stars_conflict_zones::iterator czend = this->m_stars_czones.end();
    for(; czit!=czend; ++czit)
    {
      Index i = czit->first;
      Star_handle si = Trunk::get_star(i);
      facets_created(p, p_index, si, facets);
    }

    typename std::map<Facet_ijk, int>::iterator itf;
    for(itf = facets.begin(); itf != facets.end(); itf++)
    {
      std::size_t nmax = this->m_starset.size();
      if( (*itf).second != 3) // the face is not there 3 times
      {
        if((*itf).first.is_infinite()) // should not happen
          return false;

        TPoint_3 tp0, tp1, tp2;
        TPoint_3 tp = Trunk::transform_to_star_point(p, to_be_refined);

        tp0 = ((*itf).first.vertex(0) == nmax) ? tp
                                               : Trunk::transform_to_star_point(this->m_starset[(*itf).first.vertex(0)]->center_point(),to_be_refined);
        tp1 = ((*itf).first.vertex(1) == nmax) ? tp
                                               : Trunk::transform_to_star_point(this->m_starset[(*itf).first.vertex(1)]->center_point(),to_be_refined);
        tp2 = ((*itf).first.vertex(2) == nmax) ? tp
                                               : Trunk::transform_to_star_point(this->m_starset[(*itf).first.vertex(2)]->center_point(),to_be_refined);

        double sqr = to_be_refined->compute_squared_circumradius(tp0, tp1, tp2);
        if(sqr < sq_radius_bound)
        {
          // small inconsistency (radius) is forbidden. A big one is fine.
          return false;
        }
      }
    }

    return true;
  }

  bool is_valid_point_2D(const Point_3 &p,
                         const FT& sq_radius_bound, // in M_{star to_be_refined}
                         Star_handle to_be_refined,
                         Star_handle& new_star) const
  {
    Index id = Trunk::compute_conflict_zones(p);
    this->m_stars_czones.compute_elements_needing_check();

    if(!this->m_stars_czones.are_check_maps_empty())
    {
      //std::cout << "proposed FPV point creates inconsistencies in unmodified stars" << std::endl;
      return false;
    }

    if(id < 0) // no conflict
      return false;
    else if(id < (int) this->m_starset.size()) //already in star set
      return false;

    Trunk::create_star(p, id, new_star);

    bool is = check_consistency(to_be_refined, new_star, sq_radius_bound);

    if(!is)
      new_star->invalidate();

    return is;
  }

  bool is_valid_point(const Point_3 &p,
                      const FT& sq_radius_bound, // in M_{star to_be_refined}
                      Star_handle to_be_refined,
                      Star_handle& new_star) const
  {
    if(m_pick_valid_uses_3D_checks)
      return Trunk::is_valid_point_3D(p, sq_radius_bound, to_be_refined, new_star);
    else
      return is_valid_point_2D(p, sq_radius_bound, to_be_refined, new_star);
  }

  void pick_valid_output(const Refinement_point_status rp_status)
  {
    bool success = (rp_status == SUITABLE_POINT);
    if(success)
      m_pick_valid_succeeded++;
    else //POINT_IN_CONFLICT or PICK_VALID_FAILED (maybe distinguish between them?)
      m_pick_valid_failed++;

#ifdef ANISO_VERBOSE
    if((!success && m_pick_valid_failed % 100 == 0 && m_pick_valid_failed > 0) ||
       (success && m_pick_valid_succeeded % 100 == 0 && m_pick_valid_succeeded > 0))
    {
      std::cout << "Fpick_valid : ";
      std::cout << m_pick_valid_succeeded << " success and ";
      std::cout << m_pick_valid_failed << " failures" << std::endl;
    }
#endif
  }

  Refinement_point_status pick_valid(const Star_handle star, //to be refined
                                     const Facet& facet, //belongs to star and should be refined
                                     Point_3& p) const
  {
    timer_pv.start();
    // compute radius bound, in M_star
    Point_3 center;
#ifdef ANISO_USE_EXACT
    star->compute_exact_dual_intersection(facet, center);
#else
    star->compute_dual_intersection(facet, center);
#endif
    TPoint_3 tccf = Trunk::transform_to_star_point(center, star);

    // pick_valid region radius
    TPoint_3 tp2 = facet.first->vertex((facet.second + 1) % 4)->point();
    FT sq_circumradius = star->traits()->compute_squared_distance_3_object()(tccf, tp2);
    FT sq_radiusbound = this->m_criteria->beta * this->m_criteria->beta * sq_circumradius;
    FT circumradius = this->m_criteria->delta * std::sqrt(sq_circumradius);

    std::size_t tried_times = 0;
    Star_handle newstar = new Star(this->m_criteria, this->m_pConstrain, true/*surface*/);

    //possible trick#4, if(is_in_conflict(center)) then directly go to lower levels

    p = center; // first test is the circumcenter
    while(true)
    {
#ifdef ANISO_DEBUG_REFINEMENT
      star->debug_steiner_point(p, facet, false);
#endif
      if(this->m_stars_czones.status() != Stars_conflict_zones::CONFLICTS_UNKNOWN)
      {
        std::cout << "Conflicts zones should have been empty there" << std::endl;
        std::cout << this->m_stars_czones << std::endl;
      }
      this->m_stars_czones.conflicting_point() = p;

      // Pick_valid trick#3: check conflict (encroachment...) at lower levels
      //before testing the validity of the point.
      if(Mesher_lvl::is_point_in_conflict(p, false/*no insertion in lower level queue*/) ||
         !is_valid_point(p, sq_radiusbound, star, newstar))
      {
        Trunk::clear_conflict_zones();
      }
      else
      {
        timer_pv.stop();
        delete newstar;
        //here, the newstar was a fully created star that could have been pushed in the star_set.
        //Not done at the moment since it introduces non symmetrical operations and force us
        //to check whether we have n or n+1 stars at various point, etc.etc. TODO
        return SUITABLE_POINT;
      }

      if(++tried_times >= this->m_criteria->max_times_to_try_in_picking_region)
      {
        p = center;
        delete newstar;
        if(Mesher_lvl::is_point_in_conflict(p, true, true))
        {
          timer_pv.stop();
          Trunk::clear_conflict_zones();
          return POINT_IN_CONFLICT;
        }
        timer_pv.stop();
        return PICK_VALID_FAILED;
      }

      p = compute_random_steiner_point(star, facet, tccf, circumradius);
    }
  }

  Point_3 compute_insert_or_snap_point(const Star_handle star, const Facet &facet) const
  {
    Point_3 p;
#ifdef ANISO_USE_EXACT
    star->compute_exact_dual_intersection(facet, p);
#else
    star->compute_dual_intersection(facet, p);
#endif
    return p;
  }

  Refinement_point_status compute_exact_steiner_point(Star_handle to_be_refined,
                                                      const Facet& f, //facet to be refined
                                                      const bool need_picking_valid,
                                                      Point_3& steiner) const
  {
    if (need_picking_valid)
    {
      vertex_with_picking_count++;
      return pick_valid(to_be_refined, f, steiner); //not exact... TODO (boolean in pv parameters maybe)
    }
    else
    {
      vertex_without_picking_count++;
      to_be_refined->compute_exact_dual_intersection(f, steiner);

      if(Mesher_lvl::is_point_in_conflict(steiner, true, false))
        return POINT_IN_CONFLICT;

      return SUITABLE_POINT;
    }
  }

  Refinement_point_status compute_steiner_point(Star_handle to_be_refined,
                                                const Facet& f, //facet to be refined
                                                const bool need_picking_valid,
                                                Point_3& steiner) const
  {
    if(need_picking_valid)
    {
      vertex_with_picking_count++;
      return pick_valid(to_be_refined, f, steiner);
    }
    else
    {
      vertex_without_picking_count++;
      steiner = compute_insert_or_snap_point(to_be_refined, f);

      if(Mesher_lvl::is_point_in_conflict(steiner, true, false))
      {
        Trunk::clear_conflict_zones();
        return POINT_IN_CONFLICT;
      }
      return SUITABLE_POINT;
    }
  }

public:
  Anisotropic_refine_facets_3(Previous_lvl& previous,
                              Starset& starset_,
                              const Constrain_surface* pconstrain_,
                              const Criteria* criteria_,
                              const Metric_field* metric_field_,
                              DT& ch_triangulation_,
                              AABB_tree& aabb_tree_,
                              Kd_tree& kd_tree_,
                              Stars_conflict_zones& m_stars_czones_,
                              Refine_queue& refine_queue_,
                              int queue_ids_start_,
                              int queue_ids_end_)
    :
      Mesher_lvl(previous),
      Trunk(starset_, pconstrain_, criteria_, metric_field_,
           ch_triangulation_, aabb_tree_, kd_tree_, m_stars_czones_),
      m_refine_queue(refine_queue_),
      m_queue_ids_start(queue_ids_start_),
      m_queue_ids_end(queue_ids_end_),
      m_pick_valid_uses_3D_checks(false),
      m_pick_valid_succeeded(0),
      m_pick_valid_failed(0),
      m_pick_valid_skipped(0),
      m_pick_valid_rejected(0),
      m_pick_valid_max_failures(0),
      vertex_with_picking_count(0),
      vertex_without_picking_count(0),
      m_leak_counter(0),
      timer_pv(),
      timer_npv()
  {}

private:
  Anisotropic_refine_facets_3(const Self& src);
  Self& operator=(const Self& src);
}; // Anisotropic_refine_facets_3

} // Anisotropic_mesh_3
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_REFINE_FACETS_H

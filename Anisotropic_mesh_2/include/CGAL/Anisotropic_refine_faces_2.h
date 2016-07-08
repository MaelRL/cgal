#ifndef CGAL_ANISOTROPIC_MESH_2_ANISOTROPIC_REFINE_FACES_H
#define CGAL_ANISOTROPIC_MESH_2_ANISOTROPIC_REFINE_FACES_H

#include <set>
#include <vector>

#include <CGAL/Timer.h>

#include <CGAL/Face_refine_queue.h>
#include <CGAL/Stretched_Delaunay_2.h>
#include <CGAL/Anisotropic_mesher_level.h>
#include <CGAL/Anisotropic_refine_trunk.h>

#include <CGAL/helpers/combinatorics_helper.h>

#include <CGAL/IO/Star_set_output.h>

namespace CGAL
{
namespace Anisotropic_mesh_2
{

template<typename K, typename Previous_lvl>
class Anisotropic_refine_faces_2 :
    public Anisotropic_mesher_level<Stretched_Delaunay_2<K>,
                                    Anisotropic_refine_faces_2<K, Previous_lvl>,
                                    Previous_lvl>,
    public Anisotropic_refine_trunk<K>
{
private:
  typedef Anisotropic_refine_faces_2<K, Previous_lvl>             Self;
public:
  typedef Anisotropic_mesher_level<Stretched_Delaunay_2<K>,
                                   Anisotropic_refine_faces_2<K, Previous_lvl>,
                                   Previous_lvl>                   Mesher_lvl;
  typedef Anisotropic_refine_trunk<K>                              Trunk;

  typedef Stretched_Delaunay_2<K>                                  Star;
  typedef typename Star::FT                                        FT;
  typedef typename Star::Base                                      DT; // DT_2 with vertex_base_with_info
  typedef Star*                                                    Star_handle;
  typedef std::vector<Star_handle>                                 Star_vector;
  typedef typename Star_vector::iterator                           Star_iterator;
  typedef std::set<Star_handle>                                    Star_set;
  typedef typename Star::Index                                     Index;
  typedef std::set<Index>                                          Index_set;
  typedef typename Star::Point_2                                   Point_2;
  typedef typename Star::TPoint_2                                  TPoint_2;
  typedef std::set<Point_2>                                        Point_set;
  typedef typename Star::Vertex_handle                             Vertex_handle;
  typedef typename Star::Edge                                      Edge;
  typedef typename Star::Vector_2                                  Vector_2;
  typedef typename Star::Domain                                    Domain;
  typedef typename Star::Criteria                                  Criteria;
  typedef typename Star::Face                                      Face;
  typedef typename Star::Face_handle                               Face_handle;
  typedef typename Star::Face_handle_vector                        Face_handle_vector;
  typedef typename Star::Face_handle_handle                        Face_handle_handle;

  typedef CGAL::Anisotropic_mesh_2::Starset<K>                     Starset;

  typedef CGAL::Anisotropic_mesh_2::Metric_field<K>                Metric_field;
  typedef typename Metric_field::Metric                            Metric;

  typedef CGAL::Kd_tree_for_star_set<K, Star_handle>               Kd_tree;

  typedef CGAL::Anisotropic_mesh_2::Conflict_zone<K>               Conflict_zone;
  typedef CGAL::Anisotropic_mesh_2::Stars_conflict_zones<K>        Stars_conflict_zones;

  typedef CGAL::Anisotropic_mesh_2::Face_refine_queue<K>           Refine_queue;
  typedef typename Refine_queue::Rface                             Refine_face;
  typedef typename Refine_queue::Rface_set_iterator                Rface_set_iterator;

private:
  Refine_queue& m_refine_queue;

  // these two ints determine which queues (in the face refinement queue) are
  // used by this level
  const int m_queue_ids_start;
  const int m_queue_ids_end;

  //debug & info
  mutable int m_pick_valid_points_tried;
  int m_pick_valid_succeeded;
  int m_pick_valid_failed;
  int m_pick_valid_skipped;
  int m_pick_valid_rejected;
  mutable int vertex_with_picking_count;
  mutable int vertex_without_picking_count;
  mutable int m_leak_counter;

public:
  mutable CGAL::Timer timer_pv;
  mutable CGAL::Timer timer_npv;

  bool& is_active() { return Mesher_lvl::is_active(); }
  const bool& is_active() const { return Mesher_lvl::is_active(); }

  bool is_criterion_tested(int queue_id) const
  {
    return (m_queue_ids_start <= queue_id && queue_id <= m_queue_ids_end);
  }

  Refine_queue& refine_queue() { return m_refine_queue; }

public:
//functions used in the Mesher_lvl class
  void initialize_()
  {
    std::cout << "Initializing face level with ids: " << m_queue_ids_start;
    std::cout << " " << m_queue_ids_end << std::endl;

    if(this->m_starset.empty())
    {
      std::cout << "Starting with criteria: " << std::endl;
      this->m_criteria->report();

      Trunk::initialize_stars();
    }

    fill_refinement_queue();
    Mesher_lvl::is_active() = true;

#ifdef ANISO_DEBUG_REFINEMENT_PP

    Star_iterator sit = this->m_starset.begin();
    Star_iterator sitend = this->m_starset.end();
    for(; sit!=sitend; ++sit)
    {
      Star_handle star = Trunk::get_star(sit);
      star->print_vertices();
      star->print_faces();
    }
#endif

    std::ofstream out("initial.mesh");
    output_medit(this->m_starset, out);
  }

  bool is_algorithm_done_()
  {
#ifdef ANISO_EXPENSIVE_REFINE_QUEUE_CHECKS
    if(m_refine_queue.empty(m_queue_ids_start, m_queue_ids_end))
    {
      std::cout << "it says it's empty" << std::endl;
      fill_refinement_queue();
      if(m_refine_queue.empty(m_queue_ids_start, m_queue_ids_end))
      {
        std::cout << "it ain't lying" << std::endl;
        return true;
      }
      else
      {
        std::cout << "it LIED" << std::endl;
        return false;
      }
    }
#else
    return m_refine_queue.empty(m_queue_ids_start, m_queue_ids_end);
#endif
  }

  Refinement_point_status get_refinement_point_for_next_element_(Point_2& steiner_point)
  {
    Rface_set_iterator bad_face;
    bool need_picking_valid;
    Face_handle fh; // to be refined

    //Top of the prio queue is not popped yet since the ref point could encroach.
    //It'll be popped at insertion if there is no encroachment.
    if(!next_refine_face(bad_face, fh, need_picking_valid))
      return EMPTY_QUEUE;

#ifdef ANISO_DEBUG_REFINEMENT_PP
    Vertex_handle v1 = fh->vertex(0);
    Vertex_handle v2 = fh->vertex(1);
    Vertex_handle v3 = fh->vertex(2);
    Index ind_v1 = v1->info();
    Index ind_v2 = v2->info();
    Index ind_v3 = v3->info();
    std::cout << "Trying to refine : " << ind_v1 << " " << ind_v2 << " " << ind_v3 << std::endl;
    std::cout << "Bad_face belongs to: " << bad_face->star->index_in_star_set() << " npv: " << need_picking_valid << std::endl;
    std::cout << "\tp"<< v1->info() <<" : " << bad_face->star->metric().inverse_transform(v1->point()) << std::endl;
    std::cout << "\tp"<< v2->info() <<" : " << bad_face->star->metric().inverse_transform(v2->point()) << std::endl;
    std::cout << "\tp"<< v3->info() <<" : " << bad_face->star->metric().inverse_transform(v3->point()) << std::endl;
    std::cout << "\ttp"<< v1->info() <<" : " << v1->point() << std::endl;
    std::cout << "\ttp"<< v2->info() <<" : " << v2->point() << std::endl;
    std::cout << "\ttp"<< v3->info() <<" : " << v3->point() << std::endl;

    double deg_value = 1e-4;
    bool degenerate = ( (std::abs(v1->point().x()-v2->point().x()) < deg_value &&
                         std::abs(v1->point().y()-v2->point().y()) < deg_value) ||
                        (std::abs(v2->point().x()-v3->point().x()) < deg_value &&
                         std::abs(v2->point().y()-v3->point().y()) < deg_value ) ||
                        (std::abs(v1->point().x()-v3->point().x()) < deg_value &&
                         std::abs(v1->point().y()-v3->point().y()) < deg_value ) );

    if(degenerate)
      std::cout << "trying to refine a degenerate (in the metric) face" << std::endl;
#endif

    if(!this->m_criteria->max_times_to_try_in_picking_region || //skip pick_valid if the number of tries is set to 0
       (need_picking_valid &&
        this->m_starset.compute_distortion(fh) > this->m_criteria->distortion)) //pick_valid trick #1: skip pick_valid if the distortion is too high
    {
      m_pick_valid_skipped++;
      need_picking_valid = false;
    }

    // note: failure in pick_valid AND a conflicting circumcenter gives POINT_IN_CONFLICT status (trick#2 is not applied)
    Refinement_point_status rp_status = compute_steiner_point(bad_face->star, fh,
                                                              need_picking_valid, steiner_point);

    if(need_picking_valid)
      pick_valid_output(rp_status); //counts point_in_conflict as fail (not sure if should)

    if(rp_status == POINT_IN_CONFLICT)
      return rp_status;

    // pick_valid trick #2: If an element fails a pick_valid test, put it at the end
    // of the (same) queue in hope that the (successful) refinement of another element
    // will also solve the problem for the rejected element.
    if(0 && rp_status == PICK_VALID_FAILED &&
       bad_face->value != m_refine_queue.queue_min_value(bad_face->queue_type) && //nothing to do if already last
       !bad_face->prev_rejection) // only allow one rejection
    {
      Trunk::clear_conflict_zones();
      timer_npv.start();
      m_pick_valid_rejected++;
      m_refine_queue.reject_rface(bad_face->queue_type);
      timer_npv.stop();
      return get_refinement_point_for_next_element_(steiner_point);
    }

    // We already know the conflict zones if it's a suitable point from pick_valid,
    // but we need to compute them for the other cases (won't cost anything
    // if it's already known).
    Trunk::compute_conflict_zones(steiner_point);

    // same for elements needing checks
    this->m_stars_czones.compute_elements_needing_check();

    return SUITABLE_POINT;
  }

  bool test_point_conflict_from_superior_(const Point_2&, const bool, const bool)
  {
    return false;
  }

  bool insert_(const Point_2& steiner_point)
  {
    Refine_face bad_face;
    bool need_picking_valid;
    Face_handle fh; // to be refined

    //At that point, recovering the bad_face and everything else is just useful for debugging
    //Only queue's pop() is truly needed
    if(!next_refine_face_pop(bad_face, fh, need_picking_valid))
      return false; //should never happen

    Vertex_handle v1 = fh->vertex(0);
    Vertex_handle v2 = fh->vertex(1);
    Vertex_handle v3 = fh->vertex(2);

    //todo re-add those conditions, (need to template the mesher lvl with it)
    //if(!m_refinement_condition(steiner_point))
    //  return true; //false would stop refinement

    Index pid = Trunk::insert(steiner_point, true/*conditional*/);

    if(pid != static_cast<Index>(this->m_starset.size()-1))
      std::cout << "warning in insert_" << std::endl;

//Debug: check if f has been destroyed -------------------------------------
    //The face is not necessarily destroyed, its dual can simply be shortened

    Face_handle fh2;
    if(bad_face.star->is_face(v1, v2, v3, fh2))
    {
      std::cout << "Bad face still there. Bad_face.star : ";
      std::cout << bad_face.star->index_in_star_set() << ". ";
      std::cout << "starset size : " << this->m_starset.size() << std::endl;
    }
// --------------------------------------------------------------------------------------------------

#ifdef ANISO_OUTPUT_WIP
    if(this->m_starset.size()%1000 == 0)
    {
      std::ofstream out("bambimboum_wip.mesh");
      output_medit(this->m_starset, out);
    }
#endif

    return true;
  }
// End of CRTP functions

  void report()
  {
    std::cout << "face consistency : ";
    std::cout << this->m_starset.is_consistent(true /*verbose*/) << std::endl;
    //all_face_histograms(this->m_starset, this->m_pdomain, this->m_criteria);
    std::cout << "face pick_valid stats: " << std::endl;
    std::cout << "tried: " << m_pick_valid_points_tried << " || ";
    std::cout << "skipped: " << m_pick_valid_skipped << " || ";
    std::cout << "succeeded: " << m_pick_valid_succeeded << " || ";
    std::cout << "rejected: " << m_pick_valid_rejected << " || ";
    std::cout << "failed: " << m_pick_valid_failed << std::endl;
    std::cout << "face pv: " << timer_pv.time() << " " << timer_pv.intervals() << std::endl;
    std::cout << "face npv: " << timer_npv.time() << " " << timer_npv.intervals() << std::endl;
    std::cout << "face leaking: " << m_leak_counter << std::endl;
  }

public:
  //face refinement function
  bool next_refine_face(Rface_set_iterator& rface_it,
                         Face_handle& fh,
                         bool& need_picking_valid)
  {
    while(true)
    {
      if(!m_refine_queue.top(rface_it, m_queue_ids_start, m_queue_ids_end))
        return false;

      if(rface_it->star->has_face(rface_it->face, fh))
      {
        need_picking_valid = m_refine_queue.need_picking_valid(rface_it->queue_type);
        return true;
      }
      else // top of the queue does not exist anymore
        m_refine_queue.pop(m_queue_ids_start, m_queue_ids_end);
    }
  }

  // Most of that function body/parameters is debug... It could simply be :
  // f(i){queues[i]->pop();}
  bool next_refine_face_pop(Refine_face& refine_face,
                             Face_handle& fh,
                             bool& need_picking_valid)
  {
    Rface_set_iterator rfsit;
    if(!m_refine_queue.top(rfsit, m_queue_ids_start, m_queue_ids_end))
    {
      std::cout << "empty queue at face pop time...?" << std::endl; //shouldn't happen
      return false;
    }

    refine_face = *rfsit; // must copy since pop() invalidates the iterator

    if(!rfsit->star->has_face(rfsit->face, fh))
    {
      std::cout << "problems at face pop time" << std::endl;
      return false;
    }

    need_picking_valid = m_refine_queue.need_picking_valid(rfsit->queue_type);
    m_refine_queue.pop(m_queue_ids_start, m_queue_ids_end);
    return true;
  }

private:
  template<typename Fh>
  void test_face(Star_handle star, Fh fit,
                  bool force_push = false, bool check_if_in = false)
  {
    //TODO (?) A possibility might be to check whether the face already exists
    // in the queue (regardless of the queue/value/star). If this is the case,
    // then we move on to the next face.

    // What is done atm is to update the queue IF HIGHER PRIO
    // (new queue_type < old queue_type or new_value > old_value).

    // Benefits would be the cost of checking the criteria + update cost
    // but the queues would not have the exact correct priority order they should have.


#ifdef ANISO_DEBUG_REFINEMENT_PP
    std::cout << "testing face: ";
    std::cout << fit->vertex(0)->info() << " ";
    std::cout << fit->vertex(1)->info() << " ";
    std::cout << fit->vertex(2)->info() << " ";
    std::cout << " @ star: " << star->index_in_star_set() << std::endl;
#endif

    // note : distortion is now used only to speed-up pick_valid (see pick_valid trick#1)
    // over distortion : 1

    if(is_criterion_tested(m_refine_queue.over_distortion_queue) &&
       this->m_criteria->distortion > 1.)
    {
      FT over_distortion = this->m_starset.compute_distortion(fit) - this->m_criteria->distortion;
      if(over_distortion > 0.)
      {
        if(!check_if_in || !m_refine_queue.is_face_in(star, fit, over_distortion, 1))
        {
          m_refine_queue.push(star, fit, over_distortion, 1, force_push);
          if(check_if_in)
          {
            m_leak_counter++;
            std::cout << "not already in 1" << std::endl;
          }
        }
        return;
      }
    }

    // size : 2
    if(is_criterion_tested(m_refine_queue.over_circumradius_queue) &&
       this->m_criteria->face_circumradius > 0.)
    {
      FT over_circumradius = star->compute_circumradius_overflow(fit);
      if (over_circumradius > 0)
      {
        if(!check_if_in || !m_refine_queue.is_face_in(star, fit, over_circumradius, 2))
        {
          m_refine_queue.push(star, fit, over_circumradius, 2, force_push);
          if(check_if_in)
          {
            m_leak_counter++;
          }
        }
        return;
      }
    }

    // shape : 3
    if(is_criterion_tested(m_refine_queue.bad_shape_queue) &&
       this->m_criteria->face_radius_edge_ratio > 0.)
    {
      FT over_radius_edge_ratio = star->compute_radius_edge_ratio_overflow(fit);
      if (over_radius_edge_ratio > 0)
      {
        if(!check_if_in || !m_refine_queue.is_face_in(star, fit, over_radius_edge_ratio, 3))
        {
          m_refine_queue.push(star, fit, over_radius_edge_ratio, 3, force_push);
          if(check_if_in)
          {
            std::cout << "not already in 3" << std::endl;
            std::cout << "star : " << star->index_in_star_set();
            std::cout << " || ids: ";
            std::cout << fit->vertex(0)->info() << " ";
            std::cout << fit->vertex(1)->info() << " ";
            std::cout << fit->vertex(2)->info() << std::endl;

            std::cout << "star vertices:" << std::endl;
            star->print_vertices();
          }
        }
        return;
      }
    }

    // custom consistency : 4
#ifdef ANISO_USE_CUSTOM_CONSISTENCY
    if(is_criterion_tested(m_refine_queue.inconsistent_queue) &&
//     !this->m_starset.is_vertex_consistent(fit))
       !this->m_starset.is_flip_consistent(fit))
    {
      FT vol = star->compute_volume(fit);
      if(!check_if_in || !m_refine_queue.is_face_in(star, fit, vol, 4))
      {
        m_refine_queue.push(star, fit, vol, 4, force_push);
        if(check_if_in)
        {
          std::cout << "not already in 4" << std::endl;
          std::cout << "star : " << star->index_in_star_set();
          std::cout << " || ids: ";
          std::cout << fit->vertex(0)->info() << " ";
          std::cout << fit->vertex(1)->info() << " ";
          std::cout << fit->vertex(2)->info() << std::endl;

          std::cout << "star vertices:" << std::endl;
          star->print_vertices();
        }
      }
    }
#endif

    // consistency : 5
    if(is_criterion_tested(m_refine_queue.inconsistent_queue) &&
       !this->m_starset.is_consistent(fit))
    {
      FT vol = star->compute_volume(fit);
      if(!check_if_in || !m_refine_queue.is_face_in(star, fit, vol, 5))
      {
        m_refine_queue.push(star, fit, vol, 5, force_push);
        if(check_if_in)
        {
          std::cout << "not already in 5" << std::endl;
          std::cout << "star : " << star->index_in_star_set();
          std::cout << " || ids: ";
          std::cout << fit->vertex(0)->info() << " ";
          std::cout << fit->vertex(1)->info() << " ";
          std::cout << fit->vertex(2)->info() << std::endl;

          std::cout << "star vertices:" << std::endl;
          star->print_vertices();
        }
      }
    }
  }

public:
  //faces here are necessarily inconsistent, but need to test other criteria too
  void fill_from_unmodified_stars()
  {
    typename std::map<Index, Face_handle_vector>::iterator mit =
                                  this->m_stars_czones.faces_to_check().begin();
    typename std::map<Index, Face_handle_vector>::iterator mend =
                                  this->m_stars_czones.faces_to_check().end();
    for(; mit!=mend; ++mit)
    {
      Index i = mit->first;
      Star_handle si = Trunk::get_star(i);

      Face_handle_handle fit = mit->second.begin();
      Face_handle_handle fend = mit->second.end();
      for(; fit!=fend; ++fit)
        test_face(si, *fit, true/*force push*/);
    }
  }

  void fill_refinement_queue(Index pid)
  {
    std::cout << "fill_refinement_queue pid : " << pid << std::endl;
    std::cout << this->m_stars_czones.size() << " czones" << std::endl;

    if(this->m_stars_czones.size() == 0) // we inserted regardless of conflicts
      return fill_refinement_queue();
    else // fill from the stars that were directly in conflict with the inserted point
      fill_refinement_queue(this->m_stars_czones, pid);

    // fill from unmodified stars, where some faces could become inconsistent (or the
    // star ownership of the face in the queue should change if the face was already in another
    // queue).
    fill_from_unmodified_stars();

#if ANISO_DEBUG_REFINEMENT_PP
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
    fill_refinement_queue(this->m_starset, -1);
    std::cout << this->number_of_stars() << " ";
    m_refine_queue.print();
  }

private:
  template<typename Stars>
  void fill_refinement_queue(const Stars& stars,
                             Index relative_point,
                             bool check_if_in = false)
  {
#ifdef ANISO_DEBUG_REFINEMENT_PP
    std::cout << "fill face call with relative " << relative_point << std::endl;
    typename Stars::const_iterator cit = stars.begin();
    typename Stars::const_iterator citend = stars.end();
    for (; cit != citend; cit++)
      std::cout << " " << this->get_star(cit)->index_in_star_set();
    std::cout << std::endl;
#endif

#ifdef USE_ANISO_TIMERS
    std::clock_t start_time = clock();
#endif
    typename Stars::const_iterator si = stars.begin();
    typename Stars::const_iterator siend = stars.end();
    for (; si != siend; si++)
    {
      Star_handle star = Trunk::get_star(si);

      Face_handle_handle fit = star->finite_incident_faces_begin();
      Face_handle_handle fend = star->finite_incident_faces_end();
      for (; fit!=fend; fit++)
      {
        Face_handle fh = *fit;
/*TODO re-add refinement conditions
        Point_2 c = Trunk::compute_circumcenter(fh, star);
        if(!m_refinement_condition(c))
          continue;
*/

#ifdef ANISO_DEBUG_REFINEMENT_PP
        //check for absurdities
        Vertex_handle v1 = fh->vertex(0);
        Vertex_handle v2 = fh->vertex(1);
        Vertex_handle v3 = fh->vertex(2);
        std::cout << "check" << std::endl;
        std::cout << v1->info() << " " << v1->point() << std::endl;
        std::cout << v2->info() << " " << v2->point() << std::endl;
        std::cout << v3->info() << " " << v3->point() << std::endl;

        typename Star::Traits::Compute_squared_distance_2 csd =
            K().compute_squared_distance_2_object();

        if(csd(star->metric().inverse_transform(v1->point()), Trunk::get_star(v1->info())->center_point()) > 1e-10 ||
           csd(star->metric().inverse_transform(v2->point()), Trunk::get_star(v2->info())->center_point()) > 1e-10 ||
           csd(star->metric().inverse_transform(v3->point()), Trunk::get_star(v3->info())->center_point()) > 1e-10)
        {
          std::cout.precision(20);
          std::cout << "points differ in fill_ref_queue: " << v1->info() << " " << v2->info() << " " << v3->info() << std::endl;
          std::cout << "index in star set : " <<  star->index_in_star_set() << std::endl;
          std::cout << star->metric().inverse_transform(v1->point()) << std::endl;
          std::cout << star->metric().inverse_transform(v2->point()) << std::endl;
          std::cout << star->metric().inverse_transform(v3->point()) << std::endl;
          std::cout << Trunk::get_star(v1->info())->center_point() << std::endl;
          std::cout << Trunk::get_star(v2->info())->center_point() << std::endl;
          std::cout << Trunk::get_star(v3->info())->center_point() << std::endl;
        }
#endif
        if(!star->is_inside(fh)) //checks infinity
          continue;

        if(relative_point >= 0)
        {
          bool relative = false;
          for(int i=0; i<3; i++)
          {
            if(relative_point == fh->vertex(i)->info())
            {
              relative = true;
              break;
            }
          }

          // we do not consider non-relative faces
          if(!relative)
            continue;
        }

        test_face(star, fh, false/*no force push*/, check_if_in);
      }
    }
  }

private:
  void faces_created(Star_handle star,
                     std::map<Facet_ijk, int>& faces) const
  {
    Face_handle_handle fit = star->finite_incident_faces_begin();
    Face_handle_handle fend = star->finite_incident_faces_end();
    for(; fit!=fend; fit++)
      add_to_map(Facet_ijk(*fit), faces);
  }

  // restricted faces(i,j,k) which would be created by the insertion
  // of 'p' in 'star' are added to 'faces'
  void faces_created(const Point_2& p, // new point, not yet in the star set
                      const Index p_id,   // index of p in star set
                      Star_handle star,
                      std::map<Facet_ijk, int>& faces) const
  {
    Index center_id = star->index_in_star_set();
    if(center_id == p_id)
      return faces_created(star, faces);

    // boundary faces of the conflict zone should already have been computed
    if(this->m_stars_czones.status() != Stars_conflict_zones::CONFLICT_ZONES_ARE_KNOWN)
      std::cout << "Zones should be known before entering faces_created..." << std::endl;
    if(this->m_stars_czones.conflict_zones().find(center_id) == this->m_stars_czones.end())
      std::cout << "Trying to compute faces_created for a star not in conflict" << std::endl;

    const Conflict_zone& star_cz = this->m_stars_czones.conflict_zone(center_id);
    const std::vector<Edge>& bedges = star_cz.boundary_edges();

    typename std::vector<Edge>::const_iterator eit = bedges.begin();
    for( ; eit!=bedges.end(); eit++)
    {
      int id1 = eit->first->vertex((eit->second +1)%3)->info();
      int id2 = eit->first->vertex((eit->second +2)%3)->info();
      int sid = star->index_in_star_set();

      if(id1 == star->infinite_vertex_index() ||
         id2 == star->infinite_vertex_index())
        continue;

      if(id1 != sid && id2 != sid)
        continue;

      Point_2 c = Trunk::compute_circumcenter(this->m_starset[id1]->center_point(),
                                              this->m_starset[id2]->center_point(),
                                              p, star);
      if(Trunk::is_inside_domain(c))
        add_to_map(Facet_ijk(id1, id2, p_id), faces);
    }
  }

  bool check_consistency(Star_handle to_be_refined, //face to be refined belongs to this star
                         Star_handle new_star,     //the newly created star
                         const double& sq_radius_bound) const
  {
    // list all faces that would be created by p's insertion
    Point_2 p = new_star->center_point();
    Index p_index = new_star->index_in_star_set();

    std::map<Facet_ijk, int> faces;
    faces_created(new_star, faces);

    typename Stars_conflict_zones::iterator czit = this->m_stars_czones.begin();
    typename Stars_conflict_zones::iterator czend = this->m_stars_czones.end();
    for(; czit!=czend; ++czit)
    {
      Index i = czit->first;
      Star_handle si = Trunk::get_star(i);
      faces_created(p, p_index, si, faces);
    }

    typename std::map<Facet_ijk, int>::iterator itf;
    for(itf = faces.begin(); itf != faces.end(); itf++)
    {
      std::size_t nmax = this->m_starset.size();
      if( (*itf).second != 3) // the face is not there 3 times
      {
        if((*itf).first.is_infinite()) // should not happen
          return false;

        TPoint_2 tp0, tp1, tp2;
        TPoint_2 tp = Trunk::transform_to_star_point(p, to_be_refined);

        tp0 = (itf->first.vertex(0) == nmax) ? tp
                                             : Trunk::transform_to_star_point(this->m_starset[itf->first.vertex(0)]->center_point(),to_be_refined);
        tp1 = (itf->first.vertex(1) == nmax) ? tp
                                             : Trunk::transform_to_star_point(this->m_starset[itf->first.vertex(1)]->center_point(),to_be_refined);
        tp2 = (itf->first.vertex(2) == nmax) ? tp
                                             : Trunk::transform_to_star_point(this->m_starset[itf->first.vertex(2)]->center_point(),to_be_refined);

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

  bool is_valid_point_1D(const Point_2 &p,
                         const FT& sq_radius_bound, // in M_{star to_be_refined}
                         Star_handle to_be_refined,
                         Star_handle& new_star) const
  {
    // todo
    CGAL_assertion(false);
    return false;
  }

  bool is_valid_point(const Point_2 &p,
                      const FT& sq_radius_bound, // in M_{star to_be_refined}
                      Star_handle to_be_refined,
                      Star_handle& new_star) const
  {
    if(1) // todo switch? corresponds to 'is_3D_level' for aniso_3d
      return Trunk::is_valid_point_2D(p, sq_radius_bound, to_be_refined, new_star);
    else
      return is_valid_point_1D(p, sq_radius_bound, to_be_refined, new_star);
  }

  void pick_valid_output(const Refinement_point_status rp_status)
  {
    bool success = (rp_status == SUITABLE_POINT);
    if(success)
      m_pick_valid_succeeded++;
    else //POINT_IN_CONFLICT or PICK_VALID_FAILED (maybe distinguish between them?) todo
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
                                     const Face_handle& fh, //belongs to star and should be refined
                                     Point_2& p) const
  {
    timer_pv.start();
    // compute radius bound, in M_star
    TPoint_2 circumcenter = fh->circumcenter(*(star->traits()));
    Point_2 center = Trunk::transform_from_star_point(circumcenter, star);
    TPoint_2 tp2 = fh->vertex(0)->point();

    // pick_valid region radius
    FT sq_circumradius = star->traits()->compute_squared_distance_2_object()(circumcenter, tp2);
    FT sq_radiusbound = this->m_criteria->beta * this->m_criteria->beta * sq_circumradius;
    FT circumradius = this->m_criteria->delta * std::sqrt(sq_circumradius);

    typename Star::Traits::Compute_random_point_2 random =
                                star->traits()->compute_random_point_2_object();

    std::size_t tried_times = 0;
    Star_handle newstar = new Star(this->m_criteria, this->m_pdomain);

    p = center; // first test is the circumcenter
    while(true)
    {
#ifdef ANISO_DEBUG_REFINEMENT
      star->debug_steiner_point(p, face, false);
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
         (++m_pick_valid_points_tried &&
         !is_valid_point(p, sq_radiusbound, star, newstar)))
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

      TPoint_2 tp = random(circumcenter, circumradius);
      p = Trunk::transform_from_star_point(tp, star);
    }
  }

  Point_2 compute_insert_or_snap_point(const Star_handle star,
                                       const Face_handle& fh) const
  {
    Point_2 p = Trunk::compute_circumcenter(fh, star);
    return p;
  }

  Refinement_point_status compute_exact_steiner_point(Star_handle to_be_refined,
                                                      const Face_handle& fh, //face to be refined
                                                      const bool need_picking_valid,
                                                      Point_2& steiner) const
  {
    if (need_picking_valid)
    {
      vertex_with_picking_count++;
      return pick_valid(to_be_refined, fh, steiner); //not exact... TODO (boolean in pv parameters maybe)
    }
    else
    {
      vertex_without_picking_count++;
      to_be_refined->compute_exact_dual_intersection(fh, steiner);

      if(Mesher_lvl::is_point_in_conflict(steiner, true, false))
        return POINT_IN_CONFLICT;

      return SUITABLE_POINT;
    }
  }

  Refinement_point_status compute_steiner_point(Star_handle to_be_refined,
                                                const Face_handle& fh, //face to be refined
                                                const bool need_picking_valid,
                                                Point_2& steiner) const
  {
    if(need_picking_valid)
    {
      vertex_with_picking_count++;
      return pick_valid(to_be_refined, fh, steiner);
    }
    else
    {
      vertex_without_picking_count++;
      steiner = compute_insert_or_snap_point(to_be_refined, fh);

      if(Mesher_lvl::is_point_in_conflict(steiner, true, false))
      {
        Trunk::clear_conflict_zones();
        return POINT_IN_CONFLICT;
      }
      return SUITABLE_POINT;
    }
  }

public:
  Anisotropic_refine_faces_2(Previous_lvl& previous,
                              Starset& starset_,
                              const Domain* pdomain_,
                              const Criteria* criteria_,
                              const Metric_field* metric_field_,
                              Kd_tree& kd_tree_,
                              Stars_conflict_zones& m_stars_czones_,
                              Refine_queue& refine_queue_,
                              int queue_ids_start_,
                              int queue_ids_end_)
    :
      Mesher_lvl(previous),
      Trunk(starset_, pdomain_, criteria_, metric_field_, kd_tree_, m_stars_czones_),
      m_refine_queue(refine_queue_),
      m_queue_ids_start(queue_ids_start_),
      m_queue_ids_end(queue_ids_end_),
      m_pick_valid_points_tried(0),
      m_pick_valid_succeeded(0),
      m_pick_valid_failed(0),
      m_pick_valid_skipped(0),
      m_pick_valid_rejected(0),
      vertex_with_picking_count(0),
      vertex_without_picking_count(0),
      m_leak_counter(0),
      timer_pv(),
      timer_npv()
  {}

private:
  Anisotropic_refine_faces_2(const Self& src);
  Self& operator=(const Self& src);
}; // Anisotropic_refine_faces_2

} // Anisotropic_mesh_2
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_2_ANISOTROPIC_REFINE_faceS_H

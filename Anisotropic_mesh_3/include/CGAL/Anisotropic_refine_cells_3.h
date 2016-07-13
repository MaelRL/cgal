#ifndef CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_REFINE_CELLS_H
#define CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_REFINE_CELLS_H

#include <set>
#include <vector>

#include <CGAL/Timer.h>

#ifdef NAIVE_CELL_QUEUES
  #include <CGAL/Cell_refine_queue_naive.h>
#else
  #include <CGAL/Cell_refine_queue.h>
#endif

#include <CGAL/Stretched_Delaunay_3.h>
#include <CGAL/Anisotropic_mesher_level.h>
#include <CGAL/Anisotropic_refine_trunk.h>

#include <CGAL/helpers/combinatorics_helper.h>
#include <CGAL/helpers/histogram_helper.h>

#include <boost/chrono.hpp>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K, typename Previous_lvl>
class Anisotropic_refine_cells_3 :
    public Anisotropic_mesher_level<Stretched_Delaunay_3<K>,
                                    Anisotropic_refine_cells_3<K, Previous_lvl>,
                                    Previous_lvl>,
    public Anisotropic_refine_trunk<K>
{
private:
  typedef Anisotropic_refine_cells_3<K, Previous_lvl>              Self;
public:
  typedef Anisotropic_mesher_level<Stretched_Delaunay_3<K>,
                                   Anisotropic_refine_cells_3<K, Previous_lvl>,
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
  typedef typename Star::Cell_handle_vector                        Cell_handle_vector;
  typedef typename Star::Cell_handle_handle                        Cell_handle_handle;
  typedef typename Star::Facet_handle                              Facet_handle;
  typedef typename Star::Facet_set_iterator                        Facet_set_iterator;
  typedef typename Star::Facet                                     Facet;
  typedef typename Star::Edge                                      Edge;
  typedef typename Star::Vector_3                                  Vector_3;
  typedef typename Star::Constrain_surface                         Constrain_surface;
  typedef typename Star::Criteria                                  Criteria;

  typedef CGAL::Anisotropic_mesh_3::Starset<K>                     Starset;

  typedef CGAL::Anisotropic_mesh_3::Metric_field<K>                Metric_field;
  typedef typename Metric_field::Metric                            Metric;

  typedef CGAL::AABB_tree_bbox<K, Star>                            AABB_tree;
  typedef CGAL::AABB_bbox_primitive<Star>                          AABB_primitive;

  typedef CGAL::Kd_tree_for_star_set<K, Star_handle>               Kd_tree;

  typedef CGAL::Anisotropic_mesh_3::Conflict_zone<K>               Conflict_zone;
  typedef CGAL::Anisotropic_mesh_3::Stars_conflict_zones<K>        Stars_conflict_zones;

  typedef CGAL::Anisotropic_mesh_3::Cell_refine_queue<K>          Refine_queue;
  typedef typename Refine_queue::Rcell                            Refine_cell;
#ifndef NAIVE_CELL_QUEUES
  typedef typename Refine_queue::Rcell_set_iterator               Rcell_set_iterator;
#endif

private:
  Refine_queue& m_refine_queue;
  const int m_queue_ids_start;
  const int m_queue_ids_end;

  mutable int m_pick_valid_points_tried;
  int m_pick_valid_succeeded;
  int m_pick_valid_failed;
  int m_pick_valid_skipped;
  mutable int m_pick_valid_skipped_due_to_conflict;
  int m_pick_valid_rejected;
  mutable int vertex_with_picking_count;
  mutable int vertex_without_picking_count;
  mutable int m_leak_counter;

public:
  mutable CGAL::Timer timer_pv;
  mutable CGAL::Timer timer_npv;

#define OUTPUT_BOOST_TIMERS
#ifdef OUTPUT_BOOST_TIMERS
  mutable std::ofstream btimers_out; // n of stars, time of fill ref queue, time of pick valid, time of insert
#endif

public:
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
    std::cout << "Initializing cell level with ids: " << m_queue_ids_start;
    std::cout << " " << m_queue_ids_end << std::endl;

    if(this->m_starset.empty())
    {
      std::cout << "Starting with criteria: " << std::endl;
      this->m_criteria->report();

      Trunk::initialize_stars();
      Trunk::build_aabb_tree();
    }

    Trunk::switch_to_volume_bboxes();

    fill_refinement_queue();
    m_refine_queue.print();

    Mesher_lvl::is_active() = true;
    this->m_stars_czones.cells_need_checks() = true;

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
    boost::chrono::thread_clock::time_point start = boost::chrono::thread_clock::now();

#ifdef NAIVE_CELL_QUEUES
    Refine_cell const * bad_cell;
#else
    Rcell_set_iterator bad_cell;
#endif
    bool need_picking_valid;
    Cell_handle c;

    //Top of the prio queue is not popped yet since the ref point could encroach.
    //It'll be popped at insertion if there is no encroachment.
    if(!next_refine_cell(bad_cell, c, need_picking_valid))
      return EMPTY_QUEUE;

#ifdef ANISO_DEBUG_REFINEMENT_PP
    Vertex_handle v0 = c->vertex(0);
    Vertex_handle v1 = c->vertex(1);
    Vertex_handle v2 = c->vertex(2);
    Vertex_handle v3 = c->vertex(3);
    Index ind_v0 = v0->info();
    Index ind_v1 = v1->info();
    Index ind_v2 = v2->info();
    Index ind_v3 = v3->info();
    std::cout << "Trying to refine : " << ind_v0 << " " << ind_v1 << " " << ind_v2 << " " << ind_v3 << std::endl;
    std::cout << "Bad_cell belongs to: " << bad_cell->star->index_in_star_set() << " npv: " << need_picking_valid << std::endl;
    std::cout << "Queue type: " << bad_cell->queue_type << " and value " << bad_cell->value << std::endl;
    std::cout << "\tp"<< v0->info() <<" : " << bad_cell->star->metric().inverse_transform(v0->point()) << std::endl;
    std::cout << "\tp"<< v1->info() <<" : " << bad_cell->star->metric().inverse_transform(v1->point()) << std::endl;
    std::cout << "\tp"<< v2->info() <<" : " << bad_cell->star->metric().inverse_transform(v2->point()) << std::endl;
    std::cout << "\tp"<< v3->info() <<" : " << bad_cell->star->metric().inverse_transform(v3->point()) << std::endl;
#endif

    //std::cout << "cell to be refined has distortion: " << this->m_starset.compute_distortion(c) << std::endl;
    //is_consistent(this->m_stars, c, false/*verbose*/, bad_cell->star->index_in_star_set());

    if(!this->m_criteria->max_times_to_try_in_picking_region || //dodge Pick_valid if the number of tries is set to 0
       (need_picking_valid &&
        this->m_starset.compute_distortion(c) > this->m_criteria->distortion)) //Pick_valid trick #1
    {
      m_pick_valid_skipped++;
      need_picking_valid = false;
    }

#ifdef OUTPUT_BOOST_TIMERS
    btimers_out << this->number_of_stars() << " ";
#endif

    // note: failure in pick_valid AND a conflicting circumcenter gives POINT_IN_CONFLICT status (trick#2 is not applied)
    Refinement_point_status rp_status = compute_steiner_point(bad_cell->star, c,
                                                              need_picking_valid, steiner_point);

    if(need_picking_valid)
      pick_valid_output(rp_status);

    if(rp_status == POINT_IN_CONFLICT)
      return rp_status;

    // pick_valid trick #2: If an element fails a pick_valid test, put it at the end
    // of the (same) queue in hope that the (successful) refinement of another element
    // will also solve the problem for the rejected element.
#ifdef ANISO_REJECT_FAILED_PV
    if(rp_status == PICK_VALID_FAILED &&
       bad_cell->value != m_refine_queue.queue_min_value(bad_cell->queue_type) && //nothing to push if already last
       !bad_cell->prev_rejection) // if cell has not already been rejected
    {
      Trunk::clear_conflict_zones();
      timer_npv.start();
      m_pick_valid_rejected++;
      m_refine_queue.reject_rcell(bad_cell->queue_type);
      timer_npv.stop();
      return get_refinement_point_for_next_element_(steiner_point);
    }
#endif

    //We already know the conflict zones if it's a suitable point from pick_valid, but we need
    //to compute them for the other cases (won't cost anything if it's already known).
    Trunk::compute_conflict_zones(steiner_point);
    //same for elements needing checks
    this->m_stars_czones.compute_elements_needing_check();

#ifdef ANISO_DEBUG_REFINEMENT //need public m_stars_czones
    if(this->m_stars_czones.conflict_zones().find(bad_cell->star->index_in_star_set()) ==
       this->m_stars_czones.conflict_zones().end())
    {
      std::cout << "bad cell star not in conflict with refine point" << std::endl;
      Cell_handle useless;
      std::cout << "check the last statment: ";
      std::cout << bad_cell->star->is_conflicted(Trunk::transform_to_star_point(steiner_point,
                                                                                bad_cell->star),
                                                 useless) << std::endl;
      std::cout << "are bboxes doing silly stuff? ";
      std::cout << bad_cell->star->is_in_a_volume_delaunay_ball(Trunk::transform_to_star_point(steiner_point,
                                                                                               bad_cell->star),
                                                                useless) << std::endl;
      std::cout << "even after update? ";
      bad_cell->star->update_bbox();
      std::cout << bad_cell->star->is_conflicted(Trunk::transform_to_star_point(steiner_point,
                                                                                bad_cell->star),
                                                 useless) << std::endl;
    }
#endif

    std::cout << "duration of get_next_ref: "
              << boost::chrono::duration_cast<boost::chrono::microseconds>(boost::chrono::thread_clock::now() - start).count()
              << " micro s\n";

    return SUITABLE_POINT;
  }

  bool test_point_conflict_from_superior_(const Point_3& /*p*/,
                                          const bool /*is_queue_updated*/ = true,
                                          const bool /*need_picking_valid*/ = false) const
  {
    return false; //not used atm, but it would be [return (dist(tp, tc) < dist(tc, tv0))]
  }

  bool insert_(const Point_3& p)
  {
    boost::chrono::thread_clock::time_point start = boost::chrono::thread_clock::now();

    Refine_cell bad_cell;
    bool need_picking_valid;
    Cell_handle c;

    if(!next_refine_cell_pop(bad_cell, c, need_picking_valid))
      return false; // should never happen

    Vertex_handle v1 = c->vertex(0);
    Vertex_handle v2 = c->vertex(1);
    Vertex_handle v3 = c->vertex(2);
    Vertex_handle v4 = c->vertex(3);

    //todo re-add those conditions, (need to template the mesher lvl with it)
    //if(!m_refinement_condition(steiner_point))
    //  return true; //false would stop refinement

    Index pid;

    // The cell refpoint could be on the surface and should technically give a surface star...
    // commented code doesn't really work (not precise enough) todo
    if(0/*this->m_pConstrain->side_of_constraint(p) == CGAL::ON_ORIENTED_BOUNDARY*/) // tmp
      pid = Trunk::insert(p, true/*conditional*/, true/*surface point*/);
    else
      pid = Trunk::insert(p, true/*conditional*/);

    if(pid != static_cast<Index>(this->m_starset.size()-1))
      std::cout << "warning in insert_" << std::endl;

    //Debug: check if the cell c has been destroyed ----------------------------
    Cell_handle ctest;
    int i,j,k,l;
    if(bad_cell.star->is_cell(v1, v2, v3, v4, ctest, i, j, k, l))
    {
      std::cout << "bad cell still in..." << std::endl;
      std::cout << "star: " << bad_cell.star->index_in_star_set() << std::endl;
      std::cout << v1->info() << " " << v1->point() << std::endl;
      std::cout << v2->info() << " " << v2->point() << std::endl;
      std::cout << v3->info() << " " << v3->point() << std::endl;
      std::cout << v4->info() << " " << v4->point() << std::endl;

      this->m_starset[pid]->print_vertices();
      this->m_starset[v1->info()]->print_vertices();
      this->m_starset[v2->info()]->print_vertices();
      this->m_starset[v3->info()]->print_vertices();
      this->m_starset[v4->info()]->print_vertices();

      std::cout << "conflict check" << std::endl;
      Cell_handle useless;
      std::cout << this->m_starset[v1->info()]->bbox().xmin() << " " << this->m_starset[v1->info()]->bbox().xmax() << " ";
      std::cout << this->m_starset[v1->info()]->bbox().ymin() << " " << this->m_starset[v1->info()]->bbox().ymax() << " ";
      std::cout << this->m_starset[v1->info()]->bbox().zmin() << " " << this->m_starset[v1->info()]->bbox().zmax() << " ";
      std::cout << this->m_starset[v1->info()]->is_conflicted(this->transform_to_star_point(p, this->m_starset[v1->info()]), useless) << std::endl;
      std::cout << this->m_starset[v2->info()]->bbox().xmin() << " " << this->m_starset[v2->info()]->bbox().xmax() << " ";
      std::cout << this->m_starset[v2->info()]->bbox().ymin() << " " << this->m_starset[v2->info()]->bbox().ymax() << " ";
      std::cout << this->m_starset[v2->info()]->bbox().zmin() << " " << this->m_starset[v2->info()]->bbox().zmax() << " ";
      std::cout << this->m_starset[v2->info()]->is_conflicted(this->transform_to_star_point(p, this->m_starset[v2->info()]), useless) << std::endl;
      std::cout << this->m_starset[v3->info()]->bbox().xmin() << " " << this->m_starset[v3->info()]->bbox().xmax() << " ";
      std::cout << this->m_starset[v3->info()]->bbox().ymin() << " " << this->m_starset[v3->info()]->bbox().ymax() << " ";
      std::cout << this->m_starset[v3->info()]->bbox().zmin() << " " << this->m_starset[v3->info()]->bbox().zmax() << " ";
      std::cout << this->m_starset[v3->info()]->is_conflicted(this->transform_to_star_point(p, this->m_starset[v3->info()]), useless) << std::endl;
      std::cout << this->m_starset[v4->info()]->bbox().xmin() << " " << this->m_starset[v4->info()]->bbox().xmax() << " ";
      std::cout << this->m_starset[v4->info()]->bbox().ymin() << " " << this->m_starset[v4->info()]->bbox().ymax() << " ";
      std::cout << this->m_starset[v4->info()]->bbox().zmin() << " " << this->m_starset[v4->info()]->bbox().zmax() << " ";
      std::cout << this->m_starset[v4->info()]->is_conflicted(this->transform_to_star_point(p, this->m_starset[v4->info()]), useless) << std::endl;

      std::cout << this->m_starset[v1->info()]->metric().get_transformation() << std::endl;
      std::cout << this->m_starset[v2->info()]->metric().get_transformation() << std::endl;
      std::cout << this->m_starset[v3->info()]->metric().get_transformation() << std::endl;
      std::cout << this->m_starset[v4->info()]->metric().get_transformation() << std::endl;

      //pop back the star & switch to exact? TODO
      return false;
    }

    /*
    std::cout << "info on the newly created star: " << std::endl;
    std::map<Cell_ijkl, int> new_cells;
    Star_handle star = this->m_starset[pid];
    typename Star::Cell_handle_handle ci = star->finite_star_cells_begin();
    typename Star::Cell_handle_handle ciend = star->finite_star_cells_end();
    for (; ci != ciend; ++ci)
      if(!is_consistent(this->m_stars, *ci))
      {
        std::cout << "new incoherent cell with distortion: " << this->m_starset.compute_distortion(*ci) << std::endl;
        new_cells[Cell_ijkl(*ci)] = is_consistent(this->m_stars, *ci);
      }

    typename std::map<Cell_ijkl, int>::iterator itc;
    for(itc = new_cells.begin(); itc != new_cells.end(); ++itc)
    {
      std::cout << (*itc).first.vertex(0) << " ";
      std::cout << (*itc).first.vertex(1) << " ";
      std::cout << (*itc).first.vertex(2) << " ";
      std::cout << (*itc).first.vertex(3) << " count ";
      std::cout << (*itc).second << std::endl;
    }
    */

#ifdef ANISO_OUTPUT_WIP
    if(this->m_starset.size()%1000 == 0) //should be somewhere else todo
    {
      std::ofstream out_med("bambimboum_wip.mesh");
      output_medit(this->m_starset, out_med, false);
      std::ofstream out_dump("dump_wip.txt");
      dump(this->m_starset, out_dump);
    }
#endif

#ifdef OUTPUT_BOOST_TIMERS
    btimers_out << boost::chrono::duration_cast<boost::chrono::microseconds>(boost::chrono::thread_clock::now() - start).count() << " ";
#else
    std::cout << "duration of insert: "
              << boost::chrono::duration_cast<boost::chrono::microseconds>(boost::chrono::thread_clock::now() - start).count()
              << " micro s\n";
#endif

    return true;
  }
// end of CRTP functions

  void report()
  {
    std::cout << "cell consistency : ";
    std::cout << this->m_starset.is_consistent(true /*verbose*/, CELLS_ONLY) << std::endl;
    all_cell_histograms(this->m_starset, this->m_criteria);
    std::cout << "CELL pick_valid stats: " << std::endl;
    std::cout << "tried: " << m_pick_valid_points_tried << " || ";
    std::cout << "skipped: " << m_pick_valid_skipped << " || ";
    std::cout << "skipped (due to conflicts): " << m_pick_valid_skipped_due_to_conflict << " || ";
    std::cout << "succeeded: " << m_pick_valid_succeeded << " || ";
    std::cout << "rejected: " << m_pick_valid_rejected << " || ";
    std::cout << "failed: " << m_pick_valid_failed << std::endl;
    std::cout << "cell pv: " << timer_pv.time() << " " << timer_pv.intervals() << std::endl;
    std::cout << "cell npv: " << timer_npv.time() << " " << timer_npv.intervals() << std::endl;
    std::cout << "cell leaking: " << m_leak_counter << std::endl;
  }

private:
//Cell refinement queue functions
  //Most of that function body/parameters is debug; it could simply be f(i){queues[i]->pop();}
  bool next_refine_cell_pop(Refine_cell& refine_cell,
                            Cell_handle& cell,
                            bool& need_picking_valid)
  {
#ifdef NAIVE_CELL_QUEUES
    Refine_cell const * rcsit;
#else
    Rcell_set_iterator rcsit;
#endif
    if(!m_refine_queue.top(rcsit, m_queue_ids_start, m_queue_ids_end))
    {
      std::cout << "empty queue at cell pop time...?" << std::endl; //shouldn't happen
      return false;
    }

    refine_cell = *rcsit; //have to copy since the pop invalidates the iterator

    if(!rcsit->star->has_cell(rcsit->cell, cell))
    {
      std::cout << "problems at cell pop time" << std::endl;
      return false;
    }

    need_picking_valid = m_refine_queue.need_picking_valid(rcsit->queue_type);
    m_refine_queue.pop(m_queue_ids_start, m_queue_ids_end);
    return true;
  }

#ifdef NAIVE_CELL_QUEUES
  bool next_refine_cell(Refine_cell const *& refine_cell,
#else
  bool next_refine_cell(Rcell_set_iterator& refine_cell,
#endif
                        Cell_handle& cell,
                        bool& need_picking_valid)
  
  {
    while(true)
    {
      if(!m_refine_queue.top(refine_cell, m_queue_ids_start, m_queue_ids_end))
      {
#if 1
        std::cout << "it says it's empty" << std::endl;
        fill_refinement_queue(this->m_starset, -1);
        if(!m_refine_queue.top(refine_cell, m_queue_ids_start, m_queue_ids_end))
        {
          std::cout << "it ain't lying" << std::endl;
          return false;
        }
        else
        {
          std::cout << "it LIED" << std::endl;
          continue;
        }
#else
        return false;
#endif
      }

      if(refine_cell->star->has_cell(refine_cell->cell, cell))
      {
        if(!refine_cell->star->is_inside(cell))
        {
          m_refine_queue.pop(m_queue_ids_start, m_queue_ids_end);
          continue;
        }
        need_picking_valid = m_refine_queue.need_picking_valid(refine_cell->queue_type);
        return true;
      }
      else
        m_refine_queue.pop(m_queue_ids_start, m_queue_ids_end);
    }
  }

  void test_cell(Star_handle star, Cell_handle c,
                 bool force_push = false, bool check_if_in = false)
  {
    //note : distortion is now used only to speed-up pick_valid
    // over distortion 1

#ifndef ANISO_NO_DISTORTION_REFINEMENT
    if(is_criterion_tested(m_refine_queue.over_distortion_queue) &&
       this->m_criteria->distortion > 0.)
    {
      FT over_distortion = this->m_starset.compute_distortion(c) - this->m_criteria->distortion;
      if(over_distortion > 0.)
      {
#ifndef NAIVE_CELL_QUEUES
        if(!check_if_in || !m_refine_queue.is_cell_in(star, c, over_distortion, 1))
#endif
        {
          m_refine_queue.push(star, c, over_distortion, 1, force_push);
#ifndef NAIVE_CELL_QUEUES
          if(check_if_in)
          {
            m_leak_counter++;
            std::cout << "not already in 1" << std::endl;
          }
#endif
        }
        return;
      }
    }
#endif

    // too big 2
    if(is_criterion_tested(m_refine_queue.over_circumradius_queue) &&
       this->m_criteria->cell_circumradius > 0.)
    {
      FT over_circumradius = star->compute_circumradius_overflow(c);
      if(over_circumradius > 0)
      {
        //std::cout << "over circum : " <<  over_circumradius << std::endl;
#ifndef NAIVE_CELL_QUEUES
        if(!check_if_in || !m_refine_queue.is_cell_in(star, c, over_circumradius, 2))
#endif
        {
          m_refine_queue.push(star, c, over_circumradius, 2, force_push);
#ifndef NAIVE_CELL_QUEUES
          if(check_if_in)
          {
            m_leak_counter++;
            std::cout << "not already in 2" << std::endl;
            std::cout << "star : " << star->index_in_star_set();
            std::cout << " cell : " << c->vertex(0)->info() << " " << c->vertex(1)->info();
            std::cout << " " << c->vertex(2)->info() << " " << c->vertex(3)->info();
            std::cout << " " << over_circumradius << std::endl;
          }
#endif
        }
        return;
      }
    }

    // bad shape 3
    if(is_criterion_tested(m_refine_queue.bad_shape_queue) &&
       this->m_criteria->cell_radius_edge_ratio > 0.)
    {
      FT over_radius_edge_ratio = star->compute_radius_edge_ratio_overflow(c);
      if(over_radius_edge_ratio > 0)
      {
#ifndef NAIVE_CELL_QUEUES
        if(!check_if_in || !m_refine_queue.is_cell_in(star, c, over_radius_edge_ratio, 3))
#endif
        {
          m_refine_queue.push(star, c, over_radius_edge_ratio, 3, force_push);
#ifdef NAIVE_CELL_QUEUES
          if(check_if_in)
          {
            m_leak_counter++;
            std::cout << "not already in 3" << std::endl;
          }
#endif
        }
        return;
      }
    }

    // sliverity 4
    if(is_criterion_tested(m_refine_queue.sliver_queue) &&
       this->m_criteria->sliverity > 0.)
    {
      FT over_sliverity = star->compute_sliverity_overflow(c);
      if(over_sliverity > 0)
      {
#ifndef NAIVE_CELL_QUEUES
        if(!check_if_in || !m_refine_queue.is_cell_in(star, c, over_sliverity, 4))
#endif
        {
          m_refine_queue.push(star, c, over_sliverity, 4, force_push);
#ifndef NAIVE_CELL_QUEUES
          if(check_if_in)
          {
            m_leak_counter++;
            std::cout << "not already in 4" << std::endl;
          }
#endif
        }
        return;
      }
    }

#ifndef ANISO_NO_CONSISTENCY
    // consistency 5
    if(is_criterion_tested(m_refine_queue.inconsistent_queue) &&
       !this->m_starset.is_consistent(c))
    {
      FT vol = star->compute_volume(c);
#ifndef NAIVE_CELL_QUEUES
      if(!check_if_in || !m_refine_queue.is_cell_in(star, c, vol, 5))
#endif
      {
        m_refine_queue.push(star, c, vol, 5, force_push);
#ifdef NAIVE_CELL_QUEUES
        if(check_if_in)
        {
          m_leak_counter++;
          std::cout << "not already in 5" << std::endl;
          std::cout << "star : " << star->index_in_star_set();
          std::cout << " || ids: ";
          std::cout << c->vertex(0)->info() << " ";
          std::cout << c->vertex(1)->info() << " ";
          std::cout << c->vertex(2)->info() << " ";
          std::cout << c->vertex(3)->info();
          std::cout << " " << star->compute_volume(c) << std::endl;
        }
#endif
      }
    }
#endif
  }

public:
  //cells here are necessarily inconsistent, but need to test other criteria too
  void fill_from_unmodified_stars()
  {
    typename std::map<Index, Cell_handle_vector>::iterator mit = this->m_stars_czones.cells_to_check().begin();
    typename std::map<Index, Cell_handle_vector>::iterator mend = this->m_stars_czones.cells_to_check().end();
    for(; mit!=mend; ++mit)
    {
      Index i = mit->first;
      Star_handle si = Trunk::get_star(i);

      Cell_handle_handle chit = mit->second.begin();
      Cell_handle_handle chend = mit->second.end();
      for(; chit!=chend; ++chit)
      {
        test_cell(si, *chit, true/*force push*/);
      }
    }
  }

  template<typename Stars>
  void fill_refinement_queue(const Stars& stars,
                             Index relative_point,
                             bool check_if_in = false)
  {
/*
    std::cout << "fill cell call with relative " << relative_point << std::endl;
    typename Stars::const_iterator cit = stars.begin();
    typename Stars::const_iterator citend = stars.end();
    for (; cit != citend; cit++)
      std::cout << " " << Trunk::get_star(cit)->index_in_star_set();
    std::cout << std::endl;
*/

    typename Stars::const_iterator si = stars.begin();
    typename Stars::const_iterator siend = stars.end();
    for(; si != siend; si++)
    {
      Star_handle star = Trunk::get_star(si);

      /*
      typename Star::Facet_set_iterator fit = star->begin_restricted_facets();
      typename Star::Facet_set_iterator fend = star->end_restricted_facets();
      for(; fit != fend; fit++)
      {
        if((relative_point >= 0) && (fit->first->vertex(fit->second)->info() != relative_point))
          continue;
       //CHECK_FACE_ENCROACHMENT within the restricted facets?... Nothing should appear...
      }
      */

      Cell_handle_handle ci = star->finite_star_cells_begin();
      Cell_handle_handle ciend = star->finite_star_cells_end();
      for(; ci != ciend; ci++)
      {
        Cell_handle c = *ci;
        if(!star->is_inside(c)) //checks infinity
          continue;

        if(relative_point >= 0)
        {
          bool relative = false;
          for(int i = 0; i < 4; i++)
          {
            if(relative_point == c->vertex(i)->info())
            {
              relative = true;
              break;
            }
          }
          if(!relative)
            continue;
        }

        test_cell(star, c, false/*no force push*/, check_if_in);
      }
    }
  }

  void fill_refinement_queue(Index pid)
  {
    boost::chrono::thread_clock::time_point start = boost::chrono::thread_clock::now();

    fill_refinement_queue(this->m_stars_czones, pid);
    fill_from_unmodified_stars();

#if 0
    std::cout << "Enter fill c_ref_queue debug. Filling with all stars" << std::endl;
    fill_refinement_queue();
    m_refine_queue.print();
    std::cout << "End fill c_ref_queue debug" << std::endl;
#endif
    std::cout << this->m_starset.size() << " ";
    m_refine_queue.print();

#ifdef OUTPUT_BOOST_TIMERS
    btimers_out << boost::chrono::duration_cast<boost::chrono::microseconds>(boost::chrono::thread_clock::now() - start).count() << '\n';
#else
    std::cout << "duration of fill ref queue: "
              << boost::chrono::duration_cast<boost::chrono::microseconds>(boost::chrono::thread_clock::now() - start).count()
              << " micro s\n";
#endif
  }

  void fill_refinement_queue()
  {
    std::cout << "fill from all cell" << std::endl;
    return fill_refinement_queue(this->m_starset, -1);
  }

private:
//Point computation
  bool is_valid_point(const Point_3 &p,
                      const FT sq_radiusbound,
                      Star_handle to_be_refined,
                      Star_handle& new_star) const
  {
    return Trunk::is_valid_point_3D(p, sq_radiusbound, to_be_refined, new_star);
  }

  void pick_valid_output(const Refinement_point_status rp_status)
  {
    bool success = (rp_status == SUITABLE_POINT);
    if(success)
      m_pick_valid_succeeded++;
    else //POINT_IN_CONFLICT or PICK_VALID_FAILED
      m_pick_valid_failed++;

#ifdef ANISO_VERBOSE
    if((!success && m_pick_valid_failed % 1 == 0 && m_pick_valid_failed > 0) ||
       (success && m_pick_valid_succeeded % 1 == 0 && m_pick_valid_succeeded > 0))
    {
      std::cout << "Cpick_valid : ";
      std::cout << m_pick_valid_succeeded << " success and ";
      std::cout << m_pick_valid_failed << " failures" << std::endl;
    }
#endif
  }

  Refinement_point_status pick_valid(const Star_handle star, //to be refined
                                     const Cell_handle& cell, //belongs to star and should be refined
                                     Point_3& p) const
  {
    boost::chrono::thread_clock::time_point start = boost::chrono::thread_clock::now();

    timer_pv.start();
    TPoint_3 circumcenter = cell->circumcenter(*(star->traits()));
    TPoint_3 tp2 = cell->vertex(0)->point();

    Point_3 center = Trunk::transform_from_star_point(circumcenter, star);

    // pick_valid region radius
    FT sq_circumradius = star->traits()->compute_squared_distance_3_object() (circumcenter, tp2);
    FT sq_radiusbound = sq_circumradius *  this->m_criteria->beta * this->m_criteria->beta;
    //std::cout << "radius & beta: " << sq_circumradius << " " << this->m_criteria->beta << std::endl;
    // std::cout << "sqrb: " << sq_radiusbound << std::endl;
    FT circumradius = sqrt(sq_circumradius) * this->m_criteria->delta;
    typename Star::Traits::Compute_random_point_3 random =
                                star->traits()->compute_random_point_3_object();

    std::size_t tried_times = 0;
    Star_handle new_star = new Star(this->m_criteria, this->m_pConstrain, false/*surface star*/, Trunk::m_is_3D_level);
      //again here, the new_star is potentially a surface star TODO

    //possible trick#4, if(is_in_conflict(center)) then directly go to lower levels
//   bool is_center_in_conflict = Mesher_lvl::is_point_in_conflict(Trunk::transform_from_star_point(circumcenter, star), false, true);
//   Trunk::clear_conflict_zones();

    p = center;
    while(true)
    {
      if(this->m_stars_czones.status() != Stars_conflict_zones::CONFLICTS_UNKNOWN)
      {
        std::cout << "Conflicts zones should have been empty there" << std::endl;
        std::cout << this->m_stars_czones << std::endl;
      }
      this->m_stars_czones.conflicting_point() = p;

      // Pick_valid trick#3: check conflict (encroachment...) at lower levels
      //before testing the validity of the point.
      if((Mesher_lvl::is_point_in_conflict(p, false/*no insertion in lower level queue*/) &&
         ++m_pick_valid_skipped_due_to_conflict) ||
         (++m_pick_valid_points_tried &&
          !is_valid_point(p, sq_radiusbound, star, new_star)))
      {
        Trunk::clear_conflict_zones();
      }
      else
      {
        timer_pv.stop();
        delete new_star;

#ifdef OUTPUT_BOOST_TIMERS
        btimers_out << boost::chrono::duration_cast<boost::chrono::microseconds>(boost::chrono::thread_clock::now() - start).count() << " ";
#else
        std::cout << "duration of pv (success): "
                  << boost::chrono::duration_cast<boost::chrono::microseconds>(boost::chrono::thread_clock::now() - start).count()
                  << " micro s\n";
#endif
        return SUITABLE_POINT;
      }

      if(++tried_times >= this->m_criteria->max_times_to_try_in_picking_region)
      {
        p = center;
        delete new_star;
        if(Mesher_lvl::is_point_in_conflict(p, true/*modifies the lower level queue*/,
                                            true/*need_picking_valid*/))
        {
          timer_pv.stop();
          this->clear_conflict_zones();
          return POINT_IN_CONFLICT;
        }

        timer_pv.stop();

#ifdef OUTPUT_BOOST_TIMERS
        btimers_out << boost::chrono::duration_cast<boost::chrono::microseconds>(boost::chrono::thread_clock::now() - start).count() << " ";
#else
        std::cout << "duration of pv (failed): "
                  << boost::chrono::duration_cast<boost::chrono::microseconds>(boost::chrono::thread_clock::now() - start).count()
                  << " micro s\n";
#endif

        return PICK_VALID_FAILED;
      }

      TPoint_3 tp = random(circumcenter, circumradius);
      p = Trunk::transform_from_star_point(tp, star);
    }
  }

  Point_3 compute_insert_or_snap_point(const Star_handle star, const Cell_handle& cell) const
  {
    Point_3 p = Trunk::transform_from_star_point(star->compute_circumcenter(cell), star);
    /*
    std::cout << "in star " << star->index_in_star_set() << " cell : " << std::endl;
    std::cout << cell->vertex(0)->info() << " " << cell->vertex(0)->point() << std::endl;
    std::cout << cell->vertex(1)->info() << " " << cell->vertex(1)->point() << std::endl;
    std::cout << cell->vertex(2)->info() << " " << cell->vertex(2)->point() << std::endl;
    std::cout << cell->vertex(3)->info() << " " << cell->vertex(3)->point() << std::endl;
    std::cout << "insert or snap : " << p << std::endl;

    Cell_handle useless;
    std::cout << "check conflicts : " << std::endl;
    std::cout << m_stars[cell->vertex(0)->info()]->is_conflicted(p, useless) << " "; //no transf coz iso met
    std::cout << m_stars[cell->vertex(1)->info()]->is_conflicted(p, useless) << " "; //no transf coz iso met
    std::cout << m_stars[cell->vertex(2)->info()]->is_conflicted(p, useless) << " "; //no transf coz iso met
    std::cout << m_stars[cell->vertex(3)->info()]->is_conflicted(p, useless) << std::endl; //no transf coz iso met
    */
    return p;
  }

  Refinement_point_status compute_steiner_point(Star_handle to_be_refined,
                             const Cell_handle& c, //cell to be refined
                             const bool need_picking_valid,
                             Point_3& steiner) const
  {
    if(need_picking_valid)
    {
      vertex_with_picking_count++;
      return pick_valid(to_be_refined, c, steiner);
    }
    else
    {
      boost::chrono::thread_clock::time_point start = boost::chrono::thread_clock::now();

      vertex_without_picking_count++;
      steiner = compute_insert_or_snap_point(to_be_refined, c);

#ifdef OUTPUT_BOOST_TIMERS
      btimers_out << boost::chrono::duration_cast<boost::chrono::microseconds>(boost::chrono::thread_clock::now() - start).count() << " ";
#else
      std::cout << "duration of compute_insert_or_snap_point: "
                << boost::chrono::duration_cast<boost::chrono::microseconds>(boost::chrono::thread_clock::now() - start).count()
                << " micro s\n";
#endif

      if(Mesher_lvl::is_point_in_conflict(steiner, true /*lower level queue insertion*/,
                                          false/*need_picking_valid*/))
      {
        Trunk::clear_conflict_zones();
        return POINT_IN_CONFLICT;
      }
      return SUITABLE_POINT;
    }
  }

public:
  void attempt_consistency()
  {
    return Trunk::attempt_consistency(CELLS_ONLY);
  }

public:
  Anisotropic_refine_cells_3(Previous_lvl& previous,
                             Starset& starset_,
                             const Constrain_surface* pconstrain_,
                             const Criteria* criteria_,
                             const Metric_field* metric_field_,
                             AABB_tree& aabb_tree_,
                             Kd_tree& kd_tree_,
                             Stars_conflict_zones& m_stars_czones_,
                             Refine_queue& refine_queue_,
                             const int queue_ids_start_,
                             const int queue_ids_end_)
  :
    Mesher_lvl(previous),
    Trunk(starset_, pconstrain_, criteria_, metric_field_,
          aabb_tree_, kd_tree_, m_stars_czones_,
          true/*3D*/),
    m_refine_queue(refine_queue_),
    m_queue_ids_start(queue_ids_start_),
    m_queue_ids_end(queue_ids_end_),
    m_pick_valid_points_tried(0),
    m_pick_valid_succeeded(0),
    m_pick_valid_failed(0),
    m_pick_valid_skipped(0),
    m_pick_valid_skipped_due_to_conflict(0),
    m_pick_valid_rejected(0),
    vertex_with_picking_count(0),
    vertex_without_picking_count(0),
    m_leak_counter(0),
    timer_pv(),
    timer_npv(),
    btimers_out("boost_timers.txt")
  { 
    btimers_out.precision(20);
  }

private:
  Anisotropic_refine_cells_3(const Self& src);
  Self& operator=(const Self& src);
}; // Anisotropic_refine_cells_3

} // Anisotropic_mesh_3
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_REFINE_CELLS_H

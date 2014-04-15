#ifndef CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_REFINE_CELLS_H
#define CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_REFINE_CELLS_H

#include <set>
#include <vector>

#include <CGAL/Timer.h>

#include <CGAL/Cell_refine_queue.h>
#include <CGAL/Stretched_Delaunay_3.h>
#include <CGAL/Anisotropic_mesher_level.h>
#include <CGAL/Anisotropic_refine_trunk.h>
#include <CGAL/Star_consistency.h>

#include <CGAL/helpers/combinatorics_helper.h>

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

  typedef CGAL::Anisotropic_mesh_3::Cell_refine_queue<K>          Refine_queue;
  typedef typename Refine_queue::Rcell                            Refine_cell;
  typedef typename Refine_queue::Rcell_set_iterator               Rcell_set_iterator;

public: //tmp
  Refine_queue m_refine_queue;

private:
  int m_pick_valid_succeeded;
  int m_pick_valid_failed;
  int m_pick_valid_skipped;
  mutable int m_pick_valid_skipped_due_to_conflict;
  int m_pick_valid_rejected;
  int m_pick_valid_max_failures;
  mutable int vertex_with_picking_count;
  mutable int vertex_without_picking_count;

public:
  mutable CGAL::Timer timer_pv;
  mutable CGAL::Timer timer_npv;

public:
//functions used in the Mesher_lvl class
  void initialize_()
  {
    //remark #1
    //star set shouldn't be empty since we initialized with a few surface stars first,
    //even if we didn't refine for the facet level. Should probably check anyway. TODO

    std::cout << "\nInitializing Cells...\n";
    fill_refinement_queue(this->m_stars, -1);

    m_refine_queue.print();
    Mesher_lvl::is_active() = true;
  }

  bool is_algorithm_done_()
  {
    return m_refine_queue.empty();
  }

  Refinement_point_status get_refinement_point_for_next_element_(Point_3& steiner_point)
  {
    std::cout << "###############################################################" << std::endl;
    std::cout << "get ref point for next element cells" << std::endl;
    std::cout << "###############################################################" << std::endl;

    Rcell_set_iterator bad_cell;
    bool need_picking_valid;
    Cell_handle c;

    //Top of the prio queue is not popped yet since the ref point could encroach.
    //It'll be popped at insertion if there is no encroachment.
    if(!next_refine_cell(bad_cell, c, need_picking_valid))
      return EMPTY_QUEUE;

    //std::cout << "cell to be refined has distortion: " << Trunk::compute_distortion(c) << std::endl;
    //is_consistent(this->m_stars, c, false/*verbose*/, bad_cell->star->index_in_star_set());

    if(!this->m_criteria->max_times_to_try_in_picking_region || //dodge Pick_valid if the number of tries is set to 0
       (need_picking_valid &&
        Trunk::compute_distortion(c) > this->m_criteria->distortion)) //Pick_valid trick #1
    {
      m_pick_valid_skipped++;
      need_picking_valid = false;
    }

    // note: failure in pick_valid AND a conflicting circumcenter gives POINT_IN_CONFLICT status (trick#2 is not applied)
    Refinement_point_status rp_status = compute_steiner_point(bad_cell->star, c,
                                                              need_picking_valid, steiner_point);

    if(need_picking_valid)
      pick_valid_output(rp_status);

    if(rp_status == POINT_IN_CONFLICT)
      return rp_status;

    // pick_valid trick #2: If an element fails a pick_valid test, put it at the end
    // of the (same) queue in hope that the (succesful) refinement of another element
    // will also solve the problem for the rejected element.
    if(0 && rp_status == PICK_VALID_FAILED &&
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

    Trunk::compute_conflict_zones(steiner_point);

    return SUITABLE_POINT;
  }

  bool test_point_conflict_from_superior_(const Point_3& p,
                                          bool is_queue_updated = false) const
  {
    return false; //not used atm, but it would be [return (dist(tp, tc) < dist(tc, tv0))]
  }

  bool insert_(const Point_3& p)
  {

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

    // The cell refpoint could be on the surface and should technically give a surface star... TODO
    // commented code doesn't really work (not precise enough)
    if(0/*this->m_pConstrain->side_of_constraint(p) == CGAL::ON_ORIENTED_BOUNDARY*/)
      pid = Trunk::insert(p, true/*conditional*/, true/*surface point*/);
    else
      pid = Trunk::insert(p, true/*conditional*/);

    std::cout << "cell insertion : " << p << std::endl;

    if(pid != static_cast<Index>(this->m_stars.size()-1))
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

      this->m_stars[pid]->print_vertices();
      this->m_stars[v1->info()]->print_vertices();
      this->m_stars[v2->info()]->print_vertices();
      this->m_stars[v3->info()]->print_vertices();
      this->m_stars[v4->info()]->print_vertices();

      std::cout << "conflict check" << std::endl;
      Cell_handle useless;
      std::cout << this->m_stars[v1->info()]->bbox().xmin() << " " << this->m_stars[v1->info()]->bbox().xmax() << " ";
      std::cout << this->m_stars[v1->info()]->bbox().ymin() << " " << this->m_stars[v1->info()]->bbox().ymax() << " ";
      std::cout << this->m_stars[v1->info()]->bbox().zmin() << " " << this->m_stars[v1->info()]->bbox().zmax() << " ";
      std::cout << this->m_stars[v1->info()]->is_conflicted(this->transform_to_star_point(p, this->m_stars[v1->info()]), useless) << std::endl;
      std::cout << this->m_stars[v2->info()]->bbox().xmin() << " " << this->m_stars[v2->info()]->bbox().xmax() << " ";
      std::cout << this->m_stars[v2->info()]->bbox().ymin() << " " << this->m_stars[v2->info()]->bbox().ymax() << " ";
      std::cout << this->m_stars[v2->info()]->bbox().zmin() << " " << this->m_stars[v2->info()]->bbox().zmax() << " ";
      std::cout << this->m_stars[v2->info()]->is_conflicted(this->transform_to_star_point(p, this->m_stars[v2->info()]), useless) << std::endl;
      std::cout << this->m_stars[v3->info()]->bbox().xmin() << " " << this->m_stars[v3->info()]->bbox().xmax() << " ";
      std::cout << this->m_stars[v3->info()]->bbox().ymin() << " " << this->m_stars[v3->info()]->bbox().ymax() << " ";
      std::cout << this->m_stars[v3->info()]->bbox().zmin() << " " << this->m_stars[v3->info()]->bbox().zmax() << " ";
      std::cout << this->m_stars[v3->info()]->is_conflicted(this->transform_to_star_point(p, this->m_stars[v3->info()]), useless) << std::endl;
      std::cout << this->m_stars[v4->info()]->bbox().xmin() << " " << this->m_stars[v4->info()]->bbox().xmax() << " ";
      std::cout << this->m_stars[v4->info()]->bbox().ymin() << " " << this->m_stars[v4->info()]->bbox().ymax() << " ";
      std::cout << this->m_stars[v4->info()]->bbox().zmin() << " " << this->m_stars[v4->info()]->bbox().zmax() << " ";
      std::cout << this->m_stars[v4->info()]->is_conflicted(this->transform_to_star_point(p, this->m_stars[v4->info()]), useless) << std::endl;

      std::cout << this->m_stars[v1->info()]->metric().get_transformation() << std::endl;
      std::cout << this->m_stars[v2->info()]->metric().get_transformation() << std::endl;
      std::cout << this->m_stars[v3->info()]->metric().get_transformation() << std::endl;
      std::cout << this->m_stars[v4->info()]->metric().get_transformation() << std::endl;

      //pop back the star & switch to exact? todo
      return false;
    }

    /*
    std::cout << "info on the newly created star: " << std::endl;
    std::map<Cell_ijkl, int> new_cells;
    Star_handle star = this->m_stars[pid];
    typename Star::Cell_handle_handle ci = star->begin_finite_star_cells();
    typename Star::Cell_handle_handle ciend = star->end_finite_star_cells();
    for (; ci != ciend; ++ci)
      if(!is_consistent(this->m_stars, *ci))
      {
        std::cout << "new incoherent cell with distortion: " << Trunk::compute_distortion(*ci) << std::endl;
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

    if(this->m_stars.size()%100 == 0) //should be somewhere else todo
    {
      std::ofstream out_med("bambimboum_wip.mesh");
      output_medit(this->m_stars, out_med, false);
    }
    return true;
  }
// end of CRTP functions

  void report()
  {
    std::cout << "cell consistency : ";
    std::cout << is_consistent(this->m_stars, true /*verbose*/, CELLS_ONLY) << std::endl;
    all_cell_histograms(this->m_stars, this->m_criteria);
    std::cout << "CELL pick_valid stats: " << std::endl;
    std::cout << "skipped: " << m_pick_valid_skipped << " || ";
    std::cout << "skipped (due to conflicts): " << m_pick_valid_skipped_due_to_conflict << " || ";
    std::cout << "succeeded: " << m_pick_valid_succeeded << " || ";
    std::cout << "rejected: " << m_pick_valid_rejected << " || ";
    std::cout << "failed: " << m_pick_valid_failed << std::endl;
    std::cout << "cell pv: " << timer_pv.time() << " " << timer_pv.intervals() << std::endl;
    std::cout << "cell npv: " << timer_npv.time() << " " << timer_npv.intervals() << std::endl;
  }

private:
//Cell refinement queue functions
  //Most of that function body/parameters is debug; it could simply be f(i){queues[i]->pop();}
  bool next_refine_cell_pop(Refine_cell& refine_cell,
                             Cell_handle& cell,
                             bool& need_picking_valid)
  {
    Rcell_set_iterator rcsit;
    if(!m_refine_queue.top(rcsit))
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
    m_refine_queue.pop();
    return true;
  }

  bool next_refine_cell(Rcell_set_iterator& refine_cell,
                        Cell_handle& cell,
                        bool& need_picking_valid)
  {
    while(true)
    {
      //m_refine_queue.print();
      if(!m_refine_queue.top(refine_cell))
      {
        std::cout << "it says it's empty" << std::endl;
        fill_refinement_queue(this->m_stars, -1);
        if(!m_refine_queue.top(refine_cell))
        {
          std::cout << "it ain't lying" << std::endl;
          return false;
        }
        else
        {
          std::cout << "it LIED" << std::endl;
          continue;
        }
      }

      if(refine_cell->star->has_cell(refine_cell->cell, cell))
      {
        if(!refine_cell->star->is_inside(cell))
        {
          m_refine_queue.pop();
          continue;
        }
        need_picking_valid = m_refine_queue.need_picking_valid(refine_cell->queue_type);
        return true;
      }
      else
        m_refine_queue.pop();
    }
  }

  void test_cell(Star_handle star, Cell_handle c,
                 bool force_push = false, bool check_if_in = false)
  {
    //note : distortion is now used only to speed-up pick_valid
    // over distortion 1
    if(0 && this->m_criteria->distortion > 0.)
    {
      FT over_distortion = Trunk::compute_distortion(c) - this->m_criteria->distortion;
      if(over_distortion > 0.)
      {
        if(!check_if_in || !m_refine_queue.is_cell_in(star, c, over_distortion, 0))
        {
          m_refine_queue.push(star, c, over_distortion, 1, force_push);
          if(check_if_in)
            std::cout << "not already in 1" << std::endl;
        }
        return;
      }
    }

    // too big 2
    if(this->m_criteria->cell_circumradius > 0.)
    {
      FT over_circumradius = star->compute_circumradius_overflow(c);
      if(over_circumradius > 0)
      {
        //std::cout << "over circum : " <<  over_circumradius << std::endl;
        if(!check_if_in || !m_refine_queue.is_cell_in(star, c, over_circumradius, 1))
        {
          m_refine_queue.push(star, c, over_circumradius, 2, force_push);
          if(check_if_in)
          {
            std::cout << "not already in 2" << std::endl;
            std::cout << "star : " << star->index_in_star_set();
            std::cout << " cell : " << c->vertex(0)->info() << " " << c->vertex(1)->info();
            std::cout << " " << c->vertex(2)->info() << " " << c->vertex(3)->info();
            std::cout << " " << over_circumradius << std::endl;
          }
        }
        return;
      }
    }

    // bad shape 3
    if(this->m_criteria->cell_radius_edge_ratio > 0.)
    {
      FT over_radius_edge_ratio = star->compute_radius_edge_ratio_overflow(c);
      if(over_radius_edge_ratio > 0)
      {
        if(!check_if_in || !m_refine_queue.is_cell_in(star, c, over_radius_edge_ratio, 3))
        {
          m_refine_queue.push(star, c, over_radius_edge_ratio, 3, force_push);
          if(check_if_in)
            std::cout << "not already in 3" << std::endl;
        }
        return;
      }
    }

    // sliverity 4
    if(this->m_criteria->sliverity > 0.)
    {
      FT over_sliverity = star->compute_sliverity_overflow(c);
      if(over_sliverity > 0)
      {
        if(!check_if_in || !m_refine_queue.is_cell_in(star, c, over_sliverity, 4))
        {
          m_refine_queue.push(star, c, over_sliverity, 4, force_push);
          if(check_if_in)
            std::cout << "not already in 4" << std::endl;
        }
        return;
      }
    }

    // consistency 5
    if(!is_consistent(this->m_stars, c))
    {
      FT vol = star->compute_volume(c);
      if(!check_if_in || !m_refine_queue.is_cell_in(star, c, vol, 5))
      {
        m_refine_queue.push(star, c, vol, 5, force_push);
        if(check_if_in)
        {
          std::cout << "not already in 5" << std::endl;
          std::cout << "star : " << star->index_in_star_set();
          std::cout << " || ids: ";
          std::cout << c->vertex(0)->info() << " ";
          std::cout << c->vertex(1)->info() << " ";
          std::cout << c->vertex(2)->info() << " ";
          std::cout << c->vertex(3)->info();
          std::cout << " " << star->compute_volume(c) << std::endl;
        }
      }
    }
  }

public:
  //cells here are necessarily inconsistent, but need to test other criteria too
  void fill_from_unmodified_stars()
  {
    std::cout << "fill from unmodified @ cell level. Is empty: ";
    std::cout << this->m_stars_czones.cells_to_check().empty() << std::endl;

    typename std::map<Index, Cell_handle_vector>::const_iterator mit = this->m_stars_czones.cells_to_check().begin();
    typename std::map<Index, Cell_handle_vector>::const_iterator mend = this->m_stars_czones.cells_to_check().end();
    for(; mit!=mend; ++mit)
    {
      Index i = mit->first;
      Cell_handle_vector cell_handles = mit->second;
      Star_handle si = Trunk::get_star(i);

      Cell_handle_handle chit = cell_handles.begin();
      Cell_handle_handle chend = cell_handles.end();
      for(; chit!=chend; ++chit)
      {
        test_cell(si, *chit, true/*force push*/);
      }
    }
    m_refine_queue.print();
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

      Cell_handle_handle ci = star->begin_finite_star_cells();
      Cell_handle_handle ciend = star->end_finite_star_cells();
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

    std::cout << this->m_stars.size() << " ";
    m_refine_queue.print();
  }

  void fill_refinement_queue(Index pid)
  {
    fill_refinement_queue(this->m_stars_czones, pid);
    fill_from_unmodified_stars();
/*
    std::cout << "Enter fill c_ref_queue debug. Filling with all stars" << std::endl;
    fill_refinement_queue();
    m_refine_queue.print();
    std::cout << "End fill c_ref_queue debug" << std::endl;
*/
  }

  void fill_refinement_queue()
  {
    std::cout << "fill from all cell" << std::endl;
    return fill_refinement_queue(this->m_stars, -1, true);
  }

private:
//Point computation
  template<typename Element, typename OneMap>
  void add_to_map(const Element& e, OneMap& m) const
  {
    std::pair<typename OneMap::iterator, bool> is_insert_successful;
    is_insert_successful = m.insert(std::make_pair(e,1));
    if(!is_insert_successful.second)
      (is_insert_successful.first)->second += 1; // m[e] += 1
  }

  void cells_created(Star_handle star,
                      std::map<Cell_ijkl, int>& cells) const
  {
    Cell_handle_handle ci = star->begin_finite_star_cells();
    Cell_handle_handle ciend = star->end_finite_star_cells();
    for(; ci != ciend; ci++)
      add_to_map(Cell_ijkl(*ci), cells);
  }

  // finite cells(i,j,k,l) which would be created by the insertion
  // of 'p' in 'star' are added to 'cells'
  void cells_created(const Point_3& p, // new point, not yet in the star set
                      const int p_id,   // index of p in star set
                      Star_handle star,
                      std::map<Cell_ijkl, int>& cells) const
  {
    Index center_id = star->index_in_star_set();
    if(center_id == p_id)
      return cells_created(star, cells);

    // boundary facets of the conflict zone should already have been computed
    if(this->m_stars_czones.status() != Stars_conflict_zones::CONFLICT_ZONES_ARE_KNOWN)
      std::cout << "Zones should be known before entering facets_created..." << std::endl;
    if(this->m_stars_czones.conflict_zones().find(center_id) == this->m_stars_czones.end())
      std::cout << "Trying to compute cells_created for a star not in conflict" << std::endl;

    const Conflict_zone& star_cz = this->m_stars_czones.conflict_zone(center_id);
    const std::vector<Facet>& bfacets = star_cz.boundary_facets();

    typename std::vector<Facet>::const_iterator fit = bfacets.begin();
    for( ; fit != bfacets.end(); fit++)
    {
      int id1 = fit->first->vertex((fit->second+1)%4)->info();
      int id2 = fit->first->vertex((fit->second+2)%4)->info();
      int id3 = fit->first->vertex((fit->second+3)%4)->info();
      int sid = star->index_in_star_set();

      if(id1 == star->infinite_vertex_index() ||
         id2 == star->infinite_vertex_index() ||
         id3 == star->infinite_vertex_index())
        continue;

      if(id1 != sid && id2 != sid && id3 != sid)
        continue;

      Point_3 c = Trunk::compute_circumcenter(this->m_stars[id1]->center_point(),
                                             this->m_stars[id2]->center_point(),
                                             this->m_stars[id3]->center_point(),
                                             p, star);
      if(Trunk::is_inside_domain(c))
        add_to_map(Cell_ijkl(id1, id2, id3, p_id), cells);
    }
  }

  bool check_consistency_and_sliverity(Star_handle to_be_refined,
                                       const Star_handle& new_star,
                                       FT sq_radius_bound) const
  {
    Point_3 p = new_star->center_point();
    int p_index = new_star->index_in_star_set();

    std::map<Cell_ijkl, int> cells;
    cells_created(new_star, cells);

    typename Stars_conflict_zones::iterator czit = this->m_stars_czones.begin();
    typename Stars_conflict_zones::iterator czend = this->m_stars_czones.end();
    for(; czit!=czend; ++czit)
    {
      Index i = czit->first;
      Star_handle si = Trunk::get_star(i);
      cells_created(p, p_index, si, cells);
    }

    /*
    std::cout << "cells sum up predicts : " << std::endl;
    typename std::map<Cell_ijkl, int>::iterator itc;
    for(itc = cells.begin(); itc != cells.end(); ++itc)
    {
      std::cout << (*itc).first.vertex(0) << " ";
      std::cout << (*itc).first.vertex(1) << " ";
      std::cout << (*itc).first.vertex(2) << " ";
      std::cout << (*itc).first.vertex(3) << " count ";
      std::cout << (*itc).second << std::endl;
    }
    */

    typename std::map<Cell_ijkl, int>::iterator itc;
    for(itc = cells.begin(); itc != cells.end(); ++itc)
    {
      std::size_t nmax = Trunk::number_of_stars();
      int c0 = (*itc).first.vertex(0);
      int c1 = (*itc).first.vertex(1);
      int c2 = (*itc).first.vertex(2);
      int c3 = (*itc).first.vertex(3);
      if((*itc).second != 4 )
      {
        if((*itc).first.is_infinite())
          return false;

        TPoint_3 tp0, tp1, tp2, tp3;
        TPoint_3 tp = this->transform_to_star_point(p, to_be_refined);
        tp0 = (c0 == nmax) ? tp : this->transform_to_star_point(this->m_stars[c0]->center_point(),
                                                                to_be_refined);
        tp1 = (c1 == nmax) ? tp : this->transform_to_star_point(this->m_stars[c1]->center_point(),
                                                                to_be_refined);
        tp2 = (c2 == nmax) ? tp : this->transform_to_star_point(this->m_stars[c2]->center_point(),
                                                                to_be_refined);
        tp3 = (c3 == nmax) ? tp : this->transform_to_star_point(this->m_stars[c3]->center_point(),
                                                                to_be_refined);

        double sqr = to_be_refined->compute_squared_circumradius(tp0, tp1, tp2, tp3);
        if(sqr < sq_radius_bound)
        {
          //std::cout << "sqr sqrb: " << sqr << " " << sq_radius_bound << std::endl;
          //std::cout << "I" << std::endl;
          return false; // small inconsistency (radius) is forbidden.
        }
      }

      //should we finish the loop of consistency before looping the sliverity?
      //or the other way around? or test both at the same time? TODO

      //check sliverity
      for(int i=0; i<4; ++i)
      {
        Star_handle star;
        if((*itc).first.vertex(i) == nmax)
          star = new_star;
        else
          star = this->m_stars[(*itc).first.vertex(i)];

        Point_3 p0, p1, p2, p3;
        p0 = (c0 == nmax) ? p : this->m_stars[c0]->center_point();
        p1 = (c1 == nmax) ? p : this->m_stars[c1]->center_point();
        p2 = (c2 == nmax) ? p : this->m_stars[c2]->center_point();
        p3 = (c3 == nmax) ? p : this->m_stars[c3]->center_point();

        FT sliver_overflow = star->compute_sliverity_overflow(p0, p1, p2, p3);
        if (sliver_overflow > 0.0)
        {
          //std::cout << "sliver fail with overflow: " << sliver_overflow << std::endl;
          return false;
        }
        //add squared_radius_bound & volume? todo
      }
    }

    return true;
  }

  bool is_valid_point(const Point_3 &p,
                      const FT sq_radiusbound,
                      Star_handle to_be_refined,
                      Star_handle& new_star) const
  {
    Index id = Trunk::compute_conflict_zones(p);

    //TODO some kind of check that they already have been computed. Below isn't good enough
    if(0 && this->m_stars_czones.cells_to_check().empty())
      this->m_stars_czones.compute_elements_needing_check();

    this->m_stars_czones.compute_elements_needing_check();
    if(!this->m_stars_czones.cells_to_check().empty())
    {
      std::cout << "proposed pick_valid point creates cell inconsistencies in unmodified stars" << std::endl;
      return false;
    }

    if(id < 0) // No conflict
      return false;
    else if(id < (int) this->m_stars.size()) //already in star set
      return false;

    Trunk::create_star(p, id, new_star, new_star->is_surface_star());

    bool is = check_consistency_and_sliverity(to_be_refined, new_star, sq_radiusbound);

    if(!is)
      new_star->invalidate();

    return is;
  }

  void pick_valid_output(const Refinement_point_status rp_status)
  {
    bool success = (rp_status == SUITABLE_POINT);
    if(success)
      m_pick_valid_succeeded++;
    else //POINT_IN_CONFLICT or PICK_VALID_FAILED
      m_pick_valid_failed++;

#ifdef ANISO_VERBOSE
    if((!success && m_pick_valid_failed % 100 == 0 && m_pick_valid_failed > 0) ||
       (success && m_pick_valid_succeeded % 100 == 0 && m_pick_valid_succeeded > 0))
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
    std::cout << "pv cell" << std::endl;
    timer_pv.start();
    TPoint_3 circumcenter = cell->circumcenter(*(star->traits()));
    TPoint_3 tp2 = cell->vertex(0)->point();

    // pick_valid region radius
    FT sq_circumradius = star->traits()->compute_squared_distance_3_object() (circumcenter, tp2);
    FT sq_radiusbound = sq_circumradius *  this->m_criteria->beta * this->m_criteria->beta;
    //std::cout << "radius & beta: " << sq_circumradius << " " << this->m_criteria->beta << std::endl;
    // std::cout << "sqrb: " << sq_radiusbound << std::endl;
    FT circumradius = sqrt(sq_circumradius) * this->m_criteria->delta;
    typename Star::Traits::Compute_random_point_3 random =
                                star->traits()->compute_random_point_3_object();

    std::size_t tried_times = 0;
    Star_handle new_star = new Star(this->m_criteria, this->m_pConstrain, false/*surface star*/);
      //again here, the new_star is potentially a surface star TODO

    //possible trick#4, if(is_in_conflict(center)) then directly go to lower levels
//   bool is_center_in_conflict = Mesher_lvl::is_point_in_conflict(Trunk::transform_from_star_point(circumcenter, star), false);
//   Trunk::clear_conflict_zones();

    while(true)
    {
      TPoint_3 tp = random(circumcenter, circumradius);
      p = Trunk::transform_from_star_point(tp, star);

      if(this->m_stars_czones.status() != Stars_conflict_zones::CONFLICTS_UNKNOWN)
      {
        std::cout << "Conflicts zones should have been empty there" << std::endl;
        std::cout << this->m_stars_czones << std::endl;
      }
      this->m_stars_czones.conflicting_point() = p;

      // Pick_valid trick#3: check conflict (encroachment...) at lower levels
      //before testing the validity of the point.
      if((0 && Mesher_lvl::is_point_in_conflict(p, false/*no insertion in lower level queue*/) &&
          ++m_pick_valid_skipped_due_to_conflict) ||
         !is_valid_point(p, sq_radiusbound, star, new_star))
      {
        Trunk::clear_conflict_zones();
      }
      else
      {
//        if(is_center_in_conflict)
//          std::cout << "CENTER IN CONFLICT BUT THE MIGHTY PV PREVAILS" << std::endl;
        timer_pv.stop();
        delete new_star;
        return SUITABLE_POINT;
      }

      if(tried_times++ > this->m_criteria->max_times_to_try_in_picking_region)
      {
        p = Trunk::transform_from_star_point(circumcenter, star);
        delete new_star;
        if(Mesher_lvl::is_point_in_conflict(p)) //modifies the lower level queue
        {
          timer_pv.stop();
          this->clear_conflict_zones();
          return POINT_IN_CONFLICT;
        }
        timer_pv.stop();
        return PICK_VALID_FAILED;
      }
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
      vertex_without_picking_count++;
      steiner = compute_insert_or_snap_point(to_be_refined, c);

      if(Mesher_lvl::is_point_in_conflict(steiner))
      {
        Trunk::clear_conflict_zones();
        return POINT_IN_CONFLICT;
      }
      return SUITABLE_POINT;
    }
  }

public:
  Anisotropic_refine_cells_3(Previous_lvl& previous,
                             Star_vector& stars_,
                             const Constrain_surface* pconstrain_,
                             const Criteria* criteria_,
                             const Metric_field* metric_field_,
                             DT& ch_triangulation_,
                             AABB_tree& aabb_tree_,
                             Kd_tree& kd_tree_,
                             Stars_conflict_zones& m_stars_czones_)
  :
    Mesher_lvl(previous),
    Trunk(stars_, pconstrain_, criteria_, metric_field_,
          ch_triangulation_, aabb_tree_, kd_tree_, m_stars_czones_),
    m_refine_queue(),
    m_pick_valid_succeeded(0),
    m_pick_valid_failed(0),
    m_pick_valid_skipped(0),
    m_pick_valid_skipped_due_to_conflict(0),
    m_pick_valid_rejected(0),
    m_pick_valid_max_failures(0),
    vertex_with_picking_count(0),
    vertex_without_picking_count(0),
    timer_pv(),
    timer_npv()
  { }

private:
  Anisotropic_refine_cells_3(const Self& src);
  Self& operator=(const Self& src);
}; // Anisotropic_refine_cells_3

} // Anisotropic_mesh_3
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_REFINE_CELLS_H

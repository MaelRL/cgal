#ifndef CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_REFINE_CELLS_H
#define CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_REFINE_CELLS_H

#include <set>
#include <vector>

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
  typedef typename Star::Facet_handle                              Facet_handle;
  typedef typename Star::Facet_set_iterator                        Facet_set_iterator;
  typedef typename Star::Cell_handle_handle                        Cell_handle_handle;
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

  typedef typename CGAL::Anisotropic_mesh_3::Cell_refine_queue<K>  Refine_queue;
  typedef typename CGAL::Anisotropic_mesh_3::Refine_cell<K>        Refine_cell;

private:
  Refine_queue m_refine_queue;
  int m_pick_valid_succeeded;
  int m_pick_valid_failed;
  int m_pick_valid_skipped;
  bool m_pick_valid_causes_stop;
  int m_pick_valid_max_failures;
  mutable int vertex_with_picking_count;
  mutable int vertex_without_picking_count;

public:
//functions used in the Mesher_lvl class
  void initialize_()
  {
    //remark #1
    //star set shouldn't be empty since we initialized with a few surface stars first,
    //but maybe check anyway

    //remark #2
    //poles can be annoying, they should probably be removed at that point
    //is it as simple as remove(p1...pn)? what is the cost?

    std::cout << "\nInitializing Cells...\n";
    fill_refinement_queue(this->m_stars, -1);
    m_refine_queue.print();
    this->is_active() = true;
  }

  bool is_algorithm_done_()
  {
    //maybe if(empty), check that it really is empty with a fill
    return m_refine_queue.empty();
  }

  bool get_refinement_point_for_next_element_(Point_3& steiner_point)
  {
    int queue_type = 0;
    Refine_cell bad_cell;
    bool need_picking_valid;
    Cell_handle c;

    //Top of the prio queue is not popped yet since the ref point could encroach.
    //It'll be popped at insertion if there is no encroachment.
    if(!next_refine_cell(bad_cell, c, need_picking_valid, queue_type, false))
      return false;

    if(queue_type == 0)
    {
      std::cerr << "Error : encroachment queue for cells is not implemented.\n";
      return false;
    }
    else
    {
      if(!this->m_criteria->max_times_to_try_in_picking_region ||
         (need_picking_valid &&
          this->compute_distortion(c) > this->m_criteria->distortion))
      {
        m_pick_valid_skipped++;
        need_picking_valid = false;
      }

      bool success = compute_steiner_point(bad_cell.star, c, need_picking_valid, steiner_point);
      if(!pick_valid_output(need_picking_valid, success))
        return false;
    }

    std::cout << std::endl << this->number_of_stars() << " refine: " << steiner_point << " [" << c->vertex(0)->info() << "-";
    std::cout << c->vertex(1)->info() << "-" << c->vertex(2)->info() << "-";
    std::cout << c->vertex(3)->info() << "] " << "t:" << queue_type << " i:";
    std::cout << bad_cell.vertices[0] << " s:" << bad_cell.star->center()->info() << std::endl;
    return true;
  }

  bool test_point_conflict_from_superior_(const Point_3& p)
  {
    return false; //not used atm, but it would be [return (dist(tp, tc) < dist(tc, tv0))]
  }

  template<typename Visitor>
  bool insert_(const Point_3& p,
               const Visitor& visitor)
  {
    int queue_type = 0;
    Refine_cell bad_cell;
    bool need_picking_valid;
    Cell_handle c;

    if(!next_refine_cell(bad_cell, c, need_picking_valid, queue_type, true /*popping top*/))
      return false; // should never happen

    Vertex_handle v1 = c->vertex(0);
    Vertex_handle v2 = c->vertex(1);
    Vertex_handle v3 = c->vertex(2);
    Vertex_handle v4 = c->vertex(3);

    //todo readd those conditions, (need to template the mesher lvl with it)
    //if(!m_refinement_condition(steiner_point))
    //  return true; //false would stop refinement

    Index_set modified_stars;
    Index pid;

    //the point could accidentally be on the surface... todo...
    if(0)
      //this->m_pConstrain->side_of_constraint(p) == CGAL::ON_ORIENTED_BOUNDARY)
      // doesn't work for polyehdral
      pid = Trunk::insert(p, modified_stars, true);
    else
      pid = this->insert_in_domain(p, modified_stars, true);

    std::cout << "cell insertion : " << p << std::endl;
    if(pid != static_cast<Index>(this->m_stars.size()-1))
      std::cout << "warning in insert_" << std::endl;

    //check if the cell has been destroyed
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

    if(!modified_stars.empty())
    {
      this->fill_refinement_queues(modified_stars, pid, visitor);
    }
    if(this->m_stars.size()%100 == 0)
    {
      this->clean_stars();
      std::cout << this->m_stars.size() << " stars" << std::endl;
      std::ofstream out_off("bambimboum_wip.off");
      output_off(this->m_stars, out_off);
    }

    //debug ------------------------------------------------------------------------------
    std::cout << "fill from modified, badcellstar: ";
    std::cout << bad_cell.star->index_in_star_set() << " with relative " << pid << std::endl;
    typename Index_set::const_iterator cit = modified_stars.begin();
    typename Index_set::const_iterator citend = modified_stars.end();
    for (; cit != citend; cit++)
      std::cout << " " << *cit;
    std::cout << std::endl;

    m_refine_queue.print();
    int count1 = m_refine_queue.count();

    m_refine_queue.print_queue(2);

    std::cout << "fill check" << std::endl;
    fill_refinement_queue(m_stars, -1, true);

    std::cout << "post check queue : " << std::endl;
    m_refine_queue.print();

    int count2 = m_refine_queue.count();
    m_refine_queue.print_queue(2);

    std::cout << "double check" << std::endl;
    fill_refinement_queue(m_stars, -1, true);
    std::cout << "post double check queue : " << std::endl;
    m_refine_queue.print();
    int count3 = m_refine_queue.count();

    if(count1 != count2)
      m_refine_queue.clear();
    // --------------------------------------------------------------------------------
    */

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
    std::cout << "succeeded: " << m_pick_valid_succeeded << " || ";
    std::cout << "failed: " << m_pick_valid_failed << std::endl;
  }

private:
//Cell refinement queue functions
  bool next_refine_cell(Refine_cell &refine_cell,
                        Cell_handle &cell,
                        bool &need_picking_valid,
                        int &queue_type,
                        bool is_top_popped)
  {
    while(true)
    {
      //m_refine_queue.print();
      if(!m_refine_queue.top(refine_cell, queue_type))
      {
        std::cout << "it says it's empty" << std::endl;
        fill_refinement_queue(this->m_stars, -1);
        if(!m_refine_queue.top(refine_cell, queue_type))
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
      if(refine_cell.star->has_cell(refine_cell.vertices, cell))
      {
        if(!refine_cell.star->is_inside(cell))
        {
          m_refine_queue.pop();
          continue;
        }
        need_picking_valid = m_refine_queue.need_picking_valid(queue_type);
        if(is_top_popped)
          m_refine_queue.pop();
        return true;
      }
      else
        m_refine_queue.pop();
    }
  }

  template<typename Stars>
  void fill_refinement_queue(const Stars &modified_stars,
                             int relative_point = -1,
                             bool check_if_in = false)
  {
    typename Stars::const_iterator si = modified_stars.begin();
    typename Stars::const_iterator siend = modified_stars.end();
    for(; si != siend; si++)
    {
      Star_handle star = this->get_star(si);

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
        if(0 && !check_if_in)
        {
          std::cout << " cell : " << c->vertex(0)->info() << " " << c->vertex(1)->info();
          std::cout << " " << c->vertex(2)->info() << " " << c->vertex(3)->info();
        }
        if(!star->is_inside(c))
          continue;
        if(0 && !check_if_in)
          std::cout << " is go" << std::endl;
        assert(!star->is_infinite(c));

        if(relative_point >= 0)
        {
          bool relative = false;
          for(int i = 0; i < 4; i++)
            if(relative_point == c->vertex(i)->info())
            {
              relative = true;
              break;
            }
            if(!relative)
              continue;
        }

        /*
        std::cout << "fill ";
        std::cout << c->vertex(0)->info() << " ";
        std::cout << c->vertex(1)->info() << " ";
        std::cout << c->vertex(2)->info() << " ";
        std::cout << c->vertex(3)->info() << std::endl;
        */

//note : distortion is now used only to speed-up pick_valid
        // over distortion 1
//        if(m_criteria->distortion > 0.)
//        {
//          FT over_distortion = compute_distortion(c) - m_criteria->distortion;
//          if(over_distortion > 0.)
//          {
//            if(!check_if_in || !m_refine_queue.is_cell_in(star, c, over_distortion, 0, 1))
//            {
//              m_refine_queue.push_over_distortion(star, c, over_distortion, 0);
//              if(check_if_in)
//                std::cout << "not already in 1" << std::endl;
//            }
//            continue;
//          }
//        }

        // too big 2
        if(this->m_criteria->cell_circumradius > 0.)
        {
          FT over_circumradius = star->compute_circumradius_overflow(c);
          if(over_circumradius > 0)
          {
            //std::cout << "over circum : " <<  over_circumradius << std::endl;
            if(!check_if_in || !m_refine_queue.is_cell_in(star, c, over_circumradius, 0, 2))
            {
              m_refine_queue.push_over_circumradius(star, c, over_circumradius, 0);
              if(check_if_in)
              {
                std::cout << "not already in 2" << std::endl;
                std::cout << "star : " << star->index_in_star_set();
                std::cout << " cell : " << c->vertex(0)->info() << " " << c->vertex(1)->info();
                std::cout << " " << c->vertex(2)->info() << " " << c->vertex(3)->info();
                std::cout << " " << over_circumradius << std::endl;
              }
            }
            continue;
          }
        }

        // bad shape 3
        if(this->m_criteria->cell_radius_edge_ratio > 0.)
        {
          FT over_radius_edge_ratio = star->compute_radius_edge_ratio_overflow(c);
          if(over_radius_edge_ratio > 0)
          {
            if(!check_if_in || !m_refine_queue.is_cell_in(star, c, over_radius_edge_ratio, 0, 3))
            {
              m_refine_queue.push_bad_shape(star, c, over_radius_edge_ratio, 0);
              if(check_if_in)
                std::cout << "not already in 3" << std::endl;
            }
            continue;
          }
        }

        // sliverity 4
        if(this->m_criteria->sliverity > 0.)
        {
          FT over_sliverity = star->compute_sliverity_overflow(c);
          if(over_sliverity > 0)
          {
            if(!check_if_in || !m_refine_queue.is_cell_in(star, c, over_sliverity, 0, 4))
            {
              m_refine_queue.push_sliver(star, c, over_sliverity, 0);
              if(check_if_in)
                std::cout << "not already in 4" << std::endl;
            }
            continue;
          }
        }

        // consistency 5
        if(!is_consistent(this->m_stars, c))
        {
          if(!check_if_in || !m_refine_queue.is_cell_in(star, c, star->compute_volume(c), 0, 5))
          {
            m_refine_queue.push_inconsistent(star, c, star->compute_volume(c), 0);
            if(check_if_in)
            {
              std::cout << "not already in 5" << std::endl;
              std::cout << "star : " << star->index_in_star_set();
              std::cout << " " << star->compute_volume(c) << std::endl;
            }
          }
          continue;
        }
      } // for cells
    } // for stars

    std::cout << this->m_stars.size() << " ";
    m_refine_queue.print();
  }

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
    std::vector<Facet> bfacets;
    star->find_conflicts(p, std::back_inserter(bfacets));

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

      Point_3 c = this->compute_circumcenter(this->m_stars[id1]->center_point(),
                                             this->m_stars[id2]->center_point(),
                                             this->m_stars[id3]->center_point(),
                                             p, star);
      if(this->is_inside_domain(c))
        add_to_map(Cell_ijkl(id1, id2, id3, p_id), cells);
    }
  }

  bool check_consistency_and_sliverity(Star_handle to_be_refined,
                                       const Star_handle& new_star,
                                       const Index_set& modified_stars,
                                       FT sq_radius_bound) const
  {
    Point_3 p = new_star->center_point();
    int p_index = new_star->index_in_star_set();

    std::map<Cell_ijkl, int> cells;
    cells_created(new_star, cells);

    typename Index_set::iterator it;
    for(it = modified_stars.begin(); it != modified_stars.end(); ++it)
      cells_created(p, p_index, this->m_stars[(*it)], cells);

    typename std::map<Cell_ijkl, int>::iterator itc;
    for(itc = cells.begin(); itc != cells.end(); ++itc)
    {
      std::size_t nmax = this->number_of_stars();
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
          return false; // small inconsistency (radius) is forbidden.
      }

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
          return false;

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
    Index_set modified_stars;
    int id = this->simulate_insert_to_stars(p, modified_stars);
    if(id < 0)
      return false;
    else if(id < (int) this->m_stars.size())
    {
      new_star = this->m_stars[id]; //already in star set
      return false;
    }
    else
      this->create_star(p, id, modified_stars, new_star, new_star->is_surface_star());

    bool is = check_consistency_and_sliverity(to_be_refined, new_star, modified_stars, sq_radiusbound);
    if(!is)
    {
      new_star->invalidate();
      std::cout << ("C");
    }

    return is;
  }

  bool pick_valid_output(const bool need_picking_valid,
                         const bool success)
  {
    if(!success)
      m_pick_valid_failed++;
    else if(need_picking_valid)
      m_pick_valid_succeeded++;

#ifdef ANISO_VERBOSE
    if((!success && m_pick_valid_failed % 100 == 0 && m_pick_valid_failed > 0) ||
       (success && need_picking_valid && m_pick_valid_succeeded % 100 == 0 && m_pick_valid_succeeded > 0))
    {
      std::cout << "pick_valid : ";
      std::cout << m_pick_valid_succeeded << " success and ";
      std::cout << m_pick_valid_failed << " failures" << std::endl;
    }
#endif
    if(m_pick_valid_causes_stop && m_pick_valid_failed >= m_pick_valid_max_failures)
      return false;
    return true;
  }

  bool pick_valid(const Star_handle star, //to be refined
                  const Cell_handle &cell, //belongs to star and should be refined
                  Point_3& p) const
  {
    TPoint_3 circumcenter = cell->circumcenter(*(star->traits()));
    TPoint_3 tp2 = cell->vertex(0)->point();

    // pick_valid region radius
    FT sq_circumradius = star->traits()->compute_squared_distance_3_object() (circumcenter, tp2);
    FT sq_radiusbound = sq_circumradius *  this->m_criteria->beta * this->m_criteria->beta;
    FT circumradius = sqrt(sq_circumradius) * this->m_criteria->delta;
    typename Star::Traits::Compute_random_point_3 random =
             star->traits()->compute_random_point_3_object();

    Star_handle new_star = new Star(this->m_criteria, this->m_pConstrain, false/*surface star*/);
    std::size_t tried_times = 0;
    bool success = false;

    while(true)
    {
      TPoint_3 tp = random(circumcenter, circumradius);
      p = this->transform_from_star_point(tp, star);

      if(is_valid_point(p, sq_radiusbound, star, new_star))
      {
        success = true;
        break;
      }
      if(tried_times++ > this->m_criteria->max_times_to_try_in_picking_region)
      {
        std::cout << ("X");
        p = this->transform_from_star_point(circumcenter, star);
        break;
      }
    }
    delete new_star; //this delete could be avoided if success is true
    return success;
  }

  Point_3 compute_insert_or_snap_point(const Star_handle star, const Cell_handle& cell) const
  {
    Point_3 p = this->transform_from_star_point(star->compute_circumcenter(cell), star);
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

  bool compute_steiner_point(Star_handle to_be_refined,
                             const Cell_handle& c, //facet to be refined
                             const bool need_picking_valid,
                             Point_3& steiner) const
  {
    bool success = true;
    if(need_picking_valid)
    {
      vertex_with_picking_count++;
      success = pick_valid(to_be_refined, c, steiner);
    }
    else
    {
      vertex_without_picking_count++;
      steiner = compute_insert_or_snap_point(to_be_refined, c);
    }
    return success;
  }

public:
  Anisotropic_refine_cells_3(Previous_lvl& previous,
                             Star_vector& stars_,
                             const Constrain_surface* pconstrain_,
                             const Criteria* criteria_,
                             const Metric_field* metric_field_,
                             DT& ch_triangulation_,
                             AABB_tree& aabb_tree_,
                             Kd_tree& kd_tree_)
  :
    Mesher_lvl(previous),
    Trunk(stars_, pconstrain_, criteria_, metric_field_,
          ch_triangulation_, aabb_tree_, kd_tree_),
    m_refine_queue(),
    m_pick_valid_succeeded(0),
    m_pick_valid_failed(0),
    m_pick_valid_skipped(0),
    m_pick_valid_causes_stop(false),
    m_pick_valid_max_failures(0),
    vertex_with_picking_count(0),
    vertex_without_picking_count(0)
  { }

  ~Anisotropic_refine_cells_3() { }

private:
  Anisotropic_refine_cells_3(const Self& src);
  Self& operator=(const Self& src);
}; // Anisotropic_refine_cells_3

} // Anisotropic_mesh_3
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_REFINE_CELLS_H

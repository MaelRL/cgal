#ifndef CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_REFINE_FACETS_H
#define CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_REFINE_FACETS_H

#include <set>
#include <vector>

#include <CGAL/Facet_refine_queue.h>
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

  typedef CGAL::Anisotropic_mesh_3::Facet_refine_queue<K>          Refine_queue;
  typedef CGAL::Anisotropic_mesh_3::Refine_facet<K>                Refine_facet;

private:
  Refine_queue m_refine_queue;

  //debug & info
  int m_pick_valid_succeeded;
  int m_pick_valid_failed;
  int m_pick_valid_skipped;
  int m_pick_valid_rejected;
  bool m_pick_valid_causes_stop;
  int m_pick_valid_max_failures;
  mutable int vertex_with_picking_count;
  mutable int vertex_without_picking_count;

public:
  Refine_queue& facet_ref_queue() const { return m_refine_queue; }

public:
//functions used in the Mesher_lvl class

  void initialize_()
  {
    initialize_stars();
    this->m_ch_triangulation.infinite_vertex()->info() = -10;
    this->build_aabb_tree();

    std::cout << this->m_stars.size() << " stars and ";
    std::cout << this->count_restricted_facets() << " restricted facets" << std::endl;
    fill_refinement_queue(this->m_stars, -1);
    m_refine_queue.print();
    this->is_active() = true;

    std::ofstream out("initial.mesh");
    output_medit(this->m_stars, out);
    std::ofstream out_off("initial.off");
    output_off(this->m_stars, out_off);
  }

  bool is_algorithm_done_()
  {
    return m_refine_queue.empty();
  }

  bool get_refinement_point_for_next_element_(Point_3& steiner_point)
  {
    int queue_type = 0;
    Refine_facet bad_facet;
    bool need_picking_valid;
    Facet f; // to be refined

    //Top of the prio queue is not popped yet since the ref point could encroach.
    //It'll be popped at insertion if there is no encroachment.
    if (!next_refine_cell(bad_facet, f, need_picking_valid, queue_type, false))
      return false;

#ifdef ANISO_DEBUG_REFINEMENT_PP
    Vertex_handle v1 = f.first->vertex((f.second+1)%4);
    Vertex_handle v2 = f.first->vertex((f.second+2)%4);
    Vertex_handle v3 = f.first->vertex((f.second+3)%4);
    Index ind_v1 = v1->info();
    Index ind_v2 = v2->info();
    Index ind_v3 = v3->info();
    std::cout << "Trying to refine : " << ind_v1 << " " << ind_v2 << " " << ind_v3 << std::endl;
    std::cout << "\tp"<< v1->info() <<" : " << bad_facet.star->metric().inverse_transform(v1->point()) << std::endl;
    std::cout << "\tp"<< v2->info() <<" : " << bad_facet.star->metric().inverse_transform(v2->point()) << std::endl;
    std::cout << "\tp"<< v3->info() <<" : " << bad_facet.star->metric().inverse_transform(v3->point()) << std::endl;
    std::cout << "\tcheck coordinates : from stars, " << std::endl;
    std::cout << "\t" << m_stars[ind_v1]->center_point() << " " << &(m_stars[ind_v1]) << std::endl;
    std::cout << "\t" << m_stars[ind_v2]->center_point() << " " << &(m_stars[ind_v2]) << std::endl;
    std::cout << "\t" << m_stars[ind_v3]->center_point() << " " << &(m_stars[ind_v3]) << std::endl;
#endif

#ifdef ANISO_DEBUG
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

    if(!this->m_criteria->max_times_to_try_in_picking_region ||
       (need_picking_valid &&
        this->compute_distortion(f) > this->m_criteria->distortion))
    {
      need_picking_valid = false;
      m_pick_valid_skipped++;
    }

    bool success = compute_steiner_point(bad_facet.star, f, need_picking_valid, steiner_point);
    if(!pick_valid_output(need_picking_valid, success))
      return false;
    return true;
  }

  bool test_point_conflict_from_superior_(const Point_3& p)
  {
    Index_set modified_stars;
    this->simulate_insert_to_stars(p, modified_stars);
    typename Index_set::const_iterator si = modified_stars.begin();
    typename Index_set::const_iterator siend = modified_stars.end();
    for (; si != siend; si++)
    {
      Star_handle star = this->get_star(si);
      Facet f;
      if(star->is_encroached(p, f))
      {
        m_refine_queue.push_encroachment(star, f, star->compute_volume(f));
        return true;
      }
    }
    return false;
  }

  bool insert_(const Point_3& steiner_point)
  {
    int queue_type = 0;
    Refine_facet bad_facet;
    bool need_picking_valid;
    Facet f; // to be refined

    if (!next_refine_cell(bad_facet, f, need_picking_valid, queue_type, true/*popping top*/))
    {
      std::cout << "yeah buddy" << std::endl;
      return false; //should never happen
    }

    Vertex_handle v1 = f.first->vertex((f.second+1)%4);
    Vertex_handle v2 = f.first->vertex((f.second+2)%4);
    Vertex_handle v3 = f.first->vertex((f.second+3)%4);
    Index ind_v1 = v1->info();
    Index ind_v2 = v2->info();
    Index ind_v3 = v3->info();

    //todo readd those conditions, (need to template the mesher lvl with it)
    //if(!m_refinement_condition(steiner_point))
    //  return true; //false would stop refinement

    Index_set modified_stars;// stars that would be modified by p's insertion
    int pid = Trunk::insert(steiner_point, modified_stars, true/*conditional*/);

    std::cout << "facet insertion" << std::endl;

    if(pid != static_cast<Index>(this->m_stars.size()-1))
      std::cout << "warning in insert_" << std::endl;

#ifdef ANISO_DEBUG_REFINEMENT_PP
    std::cout << "p = " << steiner_point << std::endl;

    bool v1_in_modified = ( !m_stars[ind_v1]->has_facet(f) || std::find(modified_stars.begin(), modified_stars.end(), ind_v1) != modified_stars.end());
    bool v2_in_modified = ( !m_stars[ind_v2]->has_facet(f) || std::find(modified_stars.begin(), modified_stars.end(), ind_v2) != modified_stars.end());
    bool v3_in_modified = ( !m_stars[ind_v3]->has_facet(f) || std::find(modified_stars.begin(), modified_stars.end(), ind_v3) != modified_stars.end());

    if(!v1_in_modified || !v2_in_modified || !v3_in_modified)
      std::cout << "Some stars that contain the facet are not modified by the insertion of p" << std::endl;
#endif

    //check if f has been destroyed
    // debug -------------------------------------------------------------------------------------
    Cell_handle c;
    int i,j,k;
    if(bad_facet.star->is_facet(v1, v2, v3, c, i, j, k))
    {
      int index = 6 - i - j - k;
      Facet ff = bad_facet.star->make_canonical(Facet(c,index));

      typename K::Plane_3 fplane = bad_facet.star->triangle(ff).supporting_plane();

      Point_3 steiner2;
      if(bad_facet.star->is_restricted(ff,steiner2)
        && fplane.oriented_side(bad_facet.star->metric().transform(steiner_point))
           == fplane.oriented_side(bad_facet.star->metric().transform(steiner2)))
      {
        //if steiner_point and steiner2 are not on the same side from ff, it's legal

        modified_stars.clear();
        this->pop_back_star();
#ifdef ANISO_DEBUG_REFINEMENT
        std::cout << "found & restricted" << std::endl;
#endif
        if(bad_facet.star->is_facet(v1, v2, v3, c, i, j, k))
        {
          index = 6 - i - j - k;
          ff = bad_facet.star->make_canonical(Facet(c,index));

          compute_exact_steiner_point(bad_facet.star, ff, need_picking_valid, steiner2);
          pid = Trunk::insert(steiner2, modified_stars, true/*conditional*/);
        }

#if 1//def ANISO_DEBUG_REFINEMENT
        if(bad_facet.star->is_facet(v1, v2, v3, c, i, j, k))
        {
          index = 6 - i - j - k;
          ff = bad_facet.star->make_canonical(Facet(c,index));
          if(bad_facet.star->is_restricted(ff))
          {
            std::cout.precision(15);
            std::cout << "Bad facet still there. Bad_facet.star : " << bad_facet.star->index_in_star_set() << ". stars : " << this->m_stars.size() << std::endl;
            std::cout << "star's metric" << std::endl << bad_facet.star->metric().get_transformation() << std::endl;
            std::cout << "evs : " << bad_facet.star->metric().get_max_eigenvalue() << " ";
            std::cout << bad_facet.star->metric().get_min_eigenvalue() << " ";
            std::cout << bad_facet.star->metric().get_third_eigenvalue() << std::endl;
            std::cout << "\tp"<< v1->info() <<" : " << bad_facet.star->metric().inverse_transform(v1->point()) << std::endl;
            std::cout << "check : " << this->m_stars[v1->info()]->center_point() << std::endl;
            std::cout << "\tp"<< v2->info() <<" : " << bad_facet.star->metric().inverse_transform(v2->point()) << std::endl;
            std::cout << "check : " << this->m_stars[v2->info()]->center_point() << std::endl;
            std::cout << "\tp"<< v3->info() <<" : " << bad_facet.star->metric().inverse_transform(v3->point()) << std::endl;
            std::cout << "check : " << this->m_stars[v3->info()]->center_point() << std::endl;

#ifdef ANISO_DEBUG_REFINEMENT_COLORING
            refinement_debug_coloring = true;
            debug_coloring_p1 = bad_facet.star->metric().inverse_transform(v1->point());
            debug_coloring_p2 = bad_facet.star->metric().inverse_transform(v2->point());
            debug_coloring_p3 = bad_facet.star->metric().inverse_transform(v3->point());
            debug_coloring_p = steiner_point;
            debug_coloring_fc = bad_facet.star->metric().inverse_transform(bad_facet.star->compute_circumcenter(ff));
#endif

            std::cout << "\tdim = " << bad_facet.star->dimension();
            std::cout << ", nbv = " << bad_facet.star->number_of_vertices();
            std::cout << ", pid = " << pid;
            std::cout << ", f(p) = " << this->m_pConstrain->side_of_constraint(steiner_point);
            std::cout << ", picking : " << need_picking_valid << std::endl;
            std::cout << "\tp = " << steiner_point  << " " << this->m_pConstrain->side_of_constraint(steiner_point) << std::endl;

            Cell_handle useless;
            std::cout << "checking conflicts with all three stars : ";
            std::cout << this->m_stars[ind_v1]->is_conflicted(this->transform_to_star_point(steiner_point, this->m_stars[ind_v1]), useless) << " ";
            std::cout << this->m_stars[ind_v2]->is_conflicted(this->transform_to_star_point(steiner_point, this->m_stars[ind_v2]), useless) << " ";
            std::cout << this->m_stars[ind_v3]->is_conflicted(this->transform_to_star_point(steiner_point, this->m_stars[ind_v3]), useless) << std::endl;

            std::cout << "checking needs_aabb_update in the tree of size " << this->m_aabb_tree.size() << " : ";
            std::cout << this->m_stars[ind_v1]->bbox_needs_aabb_update() << " ";
            std::cout << this->m_stars[ind_v2]->bbox_needs_aabb_update() << " ";
            std::cout << this->m_stars[ind_v3]->bbox_needs_aabb_update() << std::endl;

            bad_facet.star->debug_steiner_point(steiner_point, ff, true);

            //pop_back_star();
            //return false;
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
    return true;
  }

//End of CRTP functions

  void report()
  {
    std::cout << "facet consistency : ";
    std::cout << is_consistent(this->m_stars, true /*verbose*/, FACETS_ONLY) << std::endl;
    std::cout << "FACET pick_valid stats: " << std::endl;
    std::cout << "skipped: " << m_pick_valid_skipped << " || ";
    std::cout << "succeeded: " << m_pick_valid_succeeded << " || ";
    std::cout << "rejected: " << m_pick_valid_rejected << " || ";
    std::cout << "failed: " << m_pick_valid_failed << std::endl;
  }

private:
//Initialization
  void initialize_medial_axis(Point_set& poles)
  {
#ifdef ANISO_VERBOSE
    std::cout << "Initialize medial axis..." << std::endl;
    std::cout << "Get poles..." << std::endl;
#endif
    this->m_pConstrain->compute_poles(poles);

    //TEMP ---------------------------------------------------------------------
    //ugly stuff to reduce the number of poles to the desired fixed amount
    std::vector<Point_3> poles_v;
    std::copy(poles.begin(), poles.end(), std::back_inserter(poles_v));
    std::random_shuffle(poles_v.begin(), poles_v.end());
    poles_v.resize(std::min((std::size_t) 500, poles_v.size()));
    //remark: random shuffle breaks the fact that only taking condition of 1/x
    //will work since the geometry will be covered properly...
    // -------------------------------------------------------------------------

#ifdef ANISO_VERBOSE
    //insert them in all stars
    std::cout << "Insert poles..." << std::endl;
#endif

    unsigned int i = 1;
    unsigned int done = 0;
    typename std::vector<Point_3>::const_iterator it;
    for(it = poles_v.begin(); it != poles_v.end(); ++it, ++i)
    {
      if(i % 100 == 0)
        this->clean_stars();

      bool conditional = (i % 50 != 0); //1/50 with no condition
      if(true) //m_refinement_condition(*it)) todo
      {
        Index_set modified_stars;
        //this->insert_in_domain(*it, modified_stars, conditional);
        ++done;
      }
    }
    this->clean_stars();

#ifdef ANISO_VERBOSE
    std::cout << done << " points." << std::endl;
#endif
  }

  void initialize_stars(const int nb = 100)
  {
#ifdef ANISO_VERBOSE
    std::cout << "Initialize "<< nb << " stars..." << std::endl;
#endif
    double approx = this->m_criteria->approximation/this->m_pConstrain->get_bounding_radius();
    approx = 1e-4;

    typename Constrain_surface::Pointset initial_points = this->m_pConstrain->get_surface_points(2*nb);

    //this should -probably- be skipped if we don't have any criterion requiring facet refinement
    Point_set poles;
    initialize_medial_axis(poles); // poles

#ifdef ANISO_VERBOSE
    std::cout << "Picked " << initial_points.size() << " initial points" << std::endl;
#endif
    typename Constrain_surface::Pointset::iterator pi = initial_points.begin();
    typename Constrain_surface::Pointset::iterator pend = initial_points.end();
    int nbdone = 0;
    for (; pi != pend && nbdone < nb; pi++)
    {
      if(nbdone % 100 == 0)
        this->clean_stars();

      std::size_t this_id = this->m_stars.size();
      int id = -1;

      //if(m_refinement_condition(*pi)) todo
        id = Trunk::insert(*pi, false/*under no condition*/);

      if(this_id == id)
        nbdone++;
    }
    this->clean_stars();

#ifdef ANISO_VERBOSE
    std::cout << "done (" << this->m_stars.size() << " stars)." << std::endl;
#endif
  }

  //Facet refinement function
  bool next_refine_cell(Refine_facet &refine_facet,
                        Facet &facet,
                        bool &need_picking_valid,
                        int &queue_type,
                        bool is_top_popped)
  {
    while(true)
    {
      if (!m_refine_queue.top(refine_facet, queue_type))
      {
        std::cout << "it says it's empty" << std::endl;
        fill_refinement_queue(this->m_stars, -1);
        if (!m_refine_queue.top(refine_facet, queue_type))
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
      if (refine_facet.star->has_facet(refine_facet.vertices, facet))
      {
        need_picking_valid = m_refine_queue.need_picking_valid(queue_type);
        if(is_top_popped)
          m_refine_queue.pop(); //double pop can't happen since
        return true;            //if the facet has been selected for pt insertion, has_facet is true
      }
      else
        m_refine_queue.pop();
    }
  }

public:
  template<typename Stars>
  void fill_refinement_queue(const Stars& modified_stars,
                             int relative_point = -1)
  {
#ifdef USE_ANISO_TIMERS
    std::clock_t start_time = clock();
#endif
    typename Stars::const_iterator si = modified_stars.begin();
    typename Stars::const_iterator siend = modified_stars.end();
    for (; si != siend; si++)
    {
      Star_handle star = this->get_star(si);
      Facet_set_iterator fi = star->begin_restricted_facets();
      Facet_set_iterator fiend = star->end_restricted_facets();
      for (; fi != fiend; fi++)
      {
        Point_3 cc;
        star->compute_dual_intersection(*fi, cc);
//        if(!m_refinement_condition(cc)) todo
//          continue;

        Cell_handle cell = fi->first;
        int offset = fi->second;

#ifdef ANISO_DEBUG_REFINEMENT_PP
        //check for absurdities
        Vertex_handle v1 = cell->vertex((offset+1)%4);
        Vertex_handle v2 = cell->vertex((offset+2)%4);
        Vertex_handle v3 = cell->vertex((offset+3)%4);

        if(CGAL::squared_distance(star->metric().inverse_transform(v1->point()), m_stars[v1->info()]->center_point()) > 1e-10 ||
           CGAL::squared_distance(star->metric().inverse_transform(v2->point()), m_stars[v2->info()]->center_point()) > 1e-10 ||
           CGAL::squared_distance(star->metric().inverse_transform(v3->point()), m_stars[v3->info()]->center_point()) > 1e-10)
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

        if (relative_point >= 0) // we don't consider not-relative facets
        {
          bool relative = false;
          for (int i = 1; i <= 3; i++)
            if (relative_point == cell->vertex((offset + i) % 4)->info())
            {
              relative = true;
              break;
            }
          if (!relative)
            continue;
        }

        //note : if a facet is encroached, the queue will be filled from conflict_from_superior
        // encroachment : 0
        //if(this->is_encroached(star, *fi))
        //{
        //  m_refine_queue.push_encroachment(star, *fi, star->compute_volume(*fi));
        //  continue;
        //}

        //note : distortion is now used only to speed-up pick_valid
        // // over distortion : 1
        //if(m_criteria->distortion > 0.)
        //{
        //  FT over_distortion = compute_distortion(*fi) - m_criteria->distortion;
        //  if(over_distortion > 0.)
        //  {
        //    m_refine_queue.push_over_distortion(star, *fi, over_distortion);
        //    continue;
        //  }
        //}

        // too big : 2
        if(this->m_criteria->facet_circumradius > 0.)
        {
          FT over_circumradius = star->compute_circumradius_overflow(*fi);
          if (over_circumradius > 0)
          {
            m_refine_queue.push_over_circumradius(star, *fi, over_circumradius);
            continue;
          }
        }
        // bad approx : 3
        if(this->m_criteria->approximation > 0.)
        {
          FT over_approx = this->sq_distance_to_surface(*fi, star) - this->m_criteria->squared_approximation;
          if(over_approx > 0.)
          {
            m_refine_queue.push_bad_approximation(star, *fi, over_approx);
            continue;
          }
        }
        // bad shape : 4
        if(this->m_criteria->facet_radius_edge_ratio > 0.)
        {
          FT over_radius_edge_ratio = star->compute_radius_edge_ratio_overflow(*fi);
          if (over_radius_edge_ratio > 0)
          {
            m_refine_queue.push_bad_shape(star, *fi, over_radius_edge_ratio);
            continue;
          }
        }
        // inconsistency : 5
        if(!is_consistent(this->m_stars, *fi))
          m_refine_queue.push_inconsistent(star, *fi, star->compute_volume(*fi));

      } // facet
    } // star
#ifdef USE_ANISO_TIMERS
    m_fill_queue_timer += duration(start_time);
#endif
    m_refine_queue.print();
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

  template<typename Element, typename OneMap>
  void add_to_map(const Element& e, OneMap& m) const
  {
    std::pair<typename OneMap::iterator, bool> is_insert_successful;
    is_insert_successful = m.insert(std::make_pair(e,1));
    if(!is_insert_successful.second)
      (is_insert_successful.first)->second += 1; // m[e] += 1
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
    typename Star::Facet_set_iterator it = star->begin_restricted_facets();
    typename Star::Facet_set_iterator iend = star->end_restricted_facets();
    for(; it != iend; it++)
      add_to_map(Facet_ijk(*it), facets);
  }

  // restricted facets(i,j,k) which would be created by the insertion
  // of 'p' in 'star' are added to 'facets'
  void facets_created(const Point_3& p, // new point, not yet in the star set
                      const int p_id,   // index of p in star set
                      Star_handle star,
                      std::map<Facet_ijk, int>& facets) const
  {
    int center_id = star->index_in_star_set();
    if(center_id == p_id)
      return facets_created(star, facets);

    // find boundary facets of the conflict zone
    std::vector<Facet> local_bfacets;
    int dim = star->dimension();
    // dimension 1
    if(dim == 1)
    {
      Edge_ij e(*(star->finite_edges_begin()));
      add_to_map(Facet_ijk(p_id, e.vertex(0), e.vertex(1)), facets);
      return;
    }

    std::set<Edge_ij> edges_around_c;
    // dimension 2
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
    else if(dim == 3 && star->find_conflicts(p, std::back_inserter(local_bfacets)))
    { // get the circular set of edges/vertices around the center in 'boundary facets'
      get_cycle(local_bfacets, edges_around_c, center_id);
    }
    else std::cerr << "Warning : in 'facets_created', no conflict zone!\n";

    // get the set (circular if dim is 3) of vertices around the center
    std::set<int> vertices_around_c;
    get_vertices(edges_around_c, std::inserter(vertices_around_c, vertices_around_c.end()));

    // facets to be added are : (p, center_i, vertices of 'vertices_around_c')
    // add them only if they are restricted to the constrain surface
    typename std::set<int>::iterator iit;
    for(iit = vertices_around_c.begin(); iit != vertices_around_c.end(); iit++)
    {
      int v_id = *iit;
      if(this->is_infinite_vertex(v_id))
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
        continue;//todo : check this is correct
      }
      bool is_finite_e1 = !e1.is_infinite();
      bool is_finite_e2 = !e2.is_infinite();
      bool is_inside_c1 = false;
      bool is_inside_c2 = false;
      if(is_finite_e1)
      {
        Point_3 cc = this->compute_circumcenter(p, star->center_point(),
                                          this->m_stars[e1.vertex(0)]->center_point(),
                                          this->m_stars[e1.vertex(1)]->center_point(),
                                          star);
        is_inside_c1 = this->is_inside_domain(cc);
      }
      if(is_finite_e2)
      {
        Point_3 cc = this->compute_circumcenter(p, star->center_point(),
                                          this->m_stars[e2.vertex(0)]->center_point(),
                                          this->m_stars[e2.vertex(1)]->center_point(),
            star);
        is_inside_c2 = this->is_inside_domain(cc);
      }

      if((is_finite_e1 && is_finite_e2 && is_inside_c1 != is_inside_c2) //finite case
          || (!is_finite_e1 && is_finite_e2 && is_inside_c2)  // finite+infinite case
          || (!is_finite_e2 && is_finite_e1 && is_inside_c1)) // finite+infinite case
        add_to_map(Facet_ijk(center_id, p_id, v_id), facets);
    }
  }

  bool check_consistency(Star_handle to_be_refined, //facet to be refined belongs to this star
                         Star_handle new_star,     //the newly created star
                         const Index_set& modified_stars,
                         const double& sq_radius_bound) const
  {
    //list all facets that would be created by p's insertion
    Point_3 p = new_star->center_point();
    int p_index = new_star->index_in_star_set();

    std::map<Facet_ijk, int> facets;
    facets_created(new_star, facets);

    typename Index_set::const_iterator it;
    for(it = modified_stars.begin(); it != modified_stars.end(); it++)
      facets_created(p, p_index, this->m_stars[(*it)], facets);

    typename std::map<Facet_ijk, int>::iterator itf;
    for(itf = facets.begin(); itf != facets.end(); itf++)
    {
      std::size_t nmax = this->m_stars.size();
      if( (*itf).second != 3) // the face is not there 3 times
      {
        if((*itf).first.is_infinite()) // should not happen
          return false;

        TPoint_3 tp0, tp1, tp2;
        TPoint_3 tp = this->transform_to_star_point(p, to_be_refined);

        tp0 = ((*itf).first.vertex(0) == nmax) ? tp
                                               : this->transform_to_star_point(this->m_stars[(*itf).first.vertex(0)]->center_point(),to_be_refined);
        tp1 = ((*itf).first.vertex(1) == nmax) ? tp
                                               : this->transform_to_star_point(this->m_stars[(*itf).first.vertex(1)]->center_point(),to_be_refined);
        tp2 = ((*itf).first.vertex(2) == nmax) ? tp
                                               : this->transform_to_star_point(this->m_stars[(*itf).first.vertex(2)]->center_point(),to_be_refined);

        double sqr = to_be_refined->compute_squared_circumradius(tp0, tp1, tp2);
        if(sqr < sq_radius_bound)
        {
          // small inconsistency (radius) is forbidden. A big one is fine.
          CGAL_PROFILER("[is_valid failure : small inconsistency]");
          return false;
        }
      }
    }

    CGAL_PROFILER("[is_valid success]");
    return true;
  }

  bool is_valid_point(const Point_3 &p,
                      const FT& sq_radius_bound, // in M_{star to_be_refined}
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
      this->create_star(p, id, modified_stars, new_star);

    bool is = check_consistency(to_be_refined, new_star, modified_stars, sq_radius_bound);

    if(!is)
      new_star->invalidate();

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
                  const Facet &facet, //belongs to star and should be refined
                  Point_3& p) const
  {
    // compute radius bound, in M_star
    Point_3 center;
#ifdef ANISO_USE_EXACT
    star->compute_exact_dual_intersection(facet, center);
#else
    star->compute_dual_intersection(facet, center);
#endif
    TPoint_3 tccf = this->transform_to_star_point(center, star);

    // pick_valid region radius
    TPoint_3 tp2 = facet.first->vertex((facet.second + 1) % 4)->point();
    FT sq_circumradius = star->traits()->compute_squared_distance_3_object()(tccf, tp2);
    FT sq_radiusbound = this->m_criteria->beta * this->m_criteria->beta * sq_circumradius;
    FT circumradius = this->m_criteria->delta * std::sqrt(sq_circumradius);

    std::size_t tried_times = 0;
    bool success = false;

    Star_handle newstar = new Star(this->m_criteria, this->m_pConstrain, true/*surface*/);

    while(true)
    {
      p = compute_random_steiner_point(star, facet, tccf, circumradius);

#ifdef ANISO_DEBUG_REFINEMENT
      //std::cout << "<? check pick_valid\n";
      star->debug_steiner_point(p, facet, false);
      //std::cout << "?>\n";
#endif
      if(is_valid_point(p, sq_radiusbound, star, newstar))
      {
        success = true;
        CGAL_HISTOGRAM_PROFILER("[iterations for pick_valid]", tried_times);
        break;
      }
      if(tried_times++ > this->m_criteria->max_times_to_try_in_picking_region)
      {
        p = center;
        CGAL_HISTOGRAM_PROFILER("[iterations for pick_valid]", tried_times);
        break;
      }
    }
    delete newstar; //this delete could be avoided if success is true
#ifdef USE_ANISO_TIMERS
    m_pick_valid_timer += duration(start_time);
#endif
    return success;
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

  bool compute_exact_steiner_point(Star_handle to_be_refined,
                                   const Facet& f, //facet to be refined
                                   const bool need_picking_valid,
                                   Point_3& steiner) const
  {
    bool success = true;
    if (need_picking_valid)
    {
      vertex_with_picking_count++;
      success = pick_valid(to_be_refined, f, steiner); //not exactly (!) exact...
    }
    else
    {
      vertex_without_picking_count++;
      to_be_refined->compute_exact_dual_intersection(f, steiner);
    }
    return success;
  }

  bool compute_steiner_point(Star_handle to_be_refined,
                             const Facet& f, //facet to be refined
                             const bool need_picking_valid,
                             Point_3& steiner) const
  {
    bool success = true;
    if (need_picking_valid)
    {
      vertex_with_picking_count++;
      success = pick_valid(to_be_refined, f, steiner);
    }
    else
    {
      vertex_without_picking_count++;
      steiner = compute_insert_or_snap_point(to_be_refined, f);
    }
    return success;
  }

public:
  Anisotropic_refine_facets_3(Previous_lvl& previous,
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
      m_pick_valid_rejected(0),
      m_pick_valid_causes_stop(false),
      m_pick_valid_max_failures(0),
      vertex_with_picking_count(0),
      vertex_without_picking_count(0)
  {}

  ~Anisotropic_refine_facets_3() { }

private:
  Anisotropic_refine_facets_3(const Self& src);
  Self& operator=(const Self& src);
}; // Anisotropic_refine_facets_3

} // Anisotropic_mesh_3
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_ANISOTROPIC_REFINE_FACETS_H

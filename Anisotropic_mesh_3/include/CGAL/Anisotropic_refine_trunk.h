#ifndef CGAL_ANISOTROPIC_MESH_3_REFINE_TRUNK_H
#define CGAL_ANISOTROPIC_MESH_3_REFINE_TRUNK_H

#include <CGAL/Stretched_Delaunay_3.h>
#include <CGAL/Constrain_surface_3.h>
#include <CGAL/Metric_field.h>
#include <CGAL/Criteria.h>

#include <CGAL/aabb_tree/aabb_tree_bbox.h>
#include <CGAL/aabb_tree/aabb_tree_bbox_primitive.h>

#include <CGAL/kd_tree/Kd_tree_for_star_set.h>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K>
class Anisotropic_refine_trunk
{
private:
  typedef Anisotropic_refine_trunk<K>                              Self;
public:
  //typedef typename CGAL::Exact_predicates_exact_constructions_kernel KExact;
  typedef K                                                        KExact;

  typedef Stretched_Delaunay_3<K>                                  Star;
  typedef typename Star::FT                                        FT;
  typedef Star*                                                    Star_handle;
  typedef typename Star::Base                                      DT; // DT_3 with vertex_base_with_info
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

private:
  typedef typename KExact::Point_3                                 Exact_Point_3;
  typedef typename KExact::Point_3                                 Exact_TPoint_3;
  typedef CGAL::Cartesian_converter<K, KExact>                     To_exact;
  typedef CGAL::Cartesian_converter<KExact, K>                     Back_from_exact;

  To_exact to_exact;
  Back_from_exact back_from_exact;

protected:
  Star_vector& m_stars;

  const Constrain_surface* m_pConstrain;
  const Criteria* m_criteria;
  const Metric_field* m_metric_field;

  DT& m_ch_triangulation;
  AABB_tree& m_aabb_tree;
  Kd_tree& m_kd_tree;

public:
  double duration(const time_t& start) const
  {
    return ((clock() - start + 0.) / ((double)CLOCKS_PER_SEC));
  }

public:
  Star_set& stars() const { return m_stars; }
  const Constrain_surface* constrain_surface() const { return m_pConstrain; }
  const Criteria* criteria() const { return m_criteria; }
  const Metric_field* metric_field() const { return m_metric_field; }

  bool empty() const { return m_stars.empty(); }
  bool is_infinite_vertex(int i) const { return (i == Star::infinite_vertex_index()); }
  std::size_t number_of_stars() const { return m_stars.size(); }

  unsigned int number_of_surface_stars() const
  {
    unsigned int nb = 0;
    for(unsigned int i = 0; i < m_stars.size(); ++i)
      if(m_stars[i]->is_surface_star())
        nb++;
    return nb;
  }

  std::size_t total_number_of_vertices() const
  {
    std::size_t nbv = 0;
    for(unsigned int i = 0; i < m_stars.size(); i++)
      nbv += m_stars[i]->number_of_vertices();
    return nbv;
  }

  Star_handle get_star(Star_handle s) const { return s; }
  Star_handle get_star(std::size_t i) const { return m_stars[i]; }
  Star_handle get_star(int i) const         { return m_stars[i]; }
  Star_handle get_star(typename Index_set::const_iterator it) const   { return m_stars[*it]; }
  Star_handle get_star(typename Star_set::const_iterator it) const    { return *it; }
  Star_handle get_star(typename Star_vector::const_iterator it) const { return *it; }

  template <typename StarIterator>
  Star_set get_stars(StarIterator begin, StarIterator end) const
  {
    Star_set stars;
    for(StarIterator it = begin; it != end; ++it)
      stars.insert(get_star(it));
    return stars;
  }

  template<typename OutputIterator>
  void all_stars(OutputIterator oit) const
  {
    for(std::size_t i = 0; i < m_stars.size(); i++)
      *oit++ = m_stars[i];
  }

public:
  TPoint_3 transform_to_star_point(const Point_3& p, Star_handle star) const
  {
    return star->metric().transform(p);
  }
  Point_3 transform_from_star_point(const TPoint_3& p, Star_handle star) const
  {
    return star->metric().inverse_transform(p);
  }

  bool is_inside_domain(const Point_3& p) const
  {
    return (m_pConstrain->side_of_constraint(p) == ON_POSITIVE_SIDE);
  }

  Point_3 compute_circumcenter(const Facet& f, Star_handle here) const
  {
    return transform_from_star_point(here->compute_circumcenter(f), here);
  }
  Point_3 compute_circumcenter(Cell_handle cell, Star_handle here) const
  {
    return transform_from_star_point(here->compute_circumcenter(cell), here);
  }

  Point_3 compute_circumcenter(const Point_3& p0, const Point_3& p1,
                               const Point_3& p2, Star_handle here) const
  {
#ifdef ANISO_USE_CC_EXACT
    return back_from_exact(
       transform_from_star_point(here->compute_circumcenter(
              transform_to_star_point(to_exact(p0), here),
              transform_to_star_point(to_exact(p1), here),
              transform_to_star_point(to_exact(p2), here)), here));
#else
    return transform_from_star_point(here->compute_circumcenter(
              transform_to_star_point(p0, here),
              transform_to_star_point(p1, here),
              transform_to_star_point(p2, here)), here);
#endif
  }

  Point_3 compute_circumcenter(const Point_3& p0, const Point_3& p1,
                               const Point_3& p2, const Point_3& p3,
                               Star_handle here) const
  {
#ifdef ANISO_USE_CC_EXACT
    return back_from_exact(
       transform_from_star_point(here->compute_circumcenter(
              transform_to_star_point(to_exact(p0), here),
              transform_to_star_point(to_exact(p1), here),
              transform_to_star_point(to_exact(p2), here),
              transform_to_star_point(to_exact(p3), here)), here));
#else
    return transform_from_star_point(here->compute_circumcenter(
              transform_to_star_point(p0, here),
              transform_to_star_point(p1, here),
              transform_to_star_point(p2, here),
              transform_to_star_point(p3, here)), here);
#endif
  }

  bool is_encroached(Star_handle star, const Facet &facet)
  {
    Index p1 = facet.first->vertex((facet.second + 1)%4)->info();
    Index p2 = facet.first->vertex((facet.second + 2)%4)->info();
    Index p3 = facet.first->vertex((facet.second + 3)%4)->info();

    Star_iterator sit = m_stars.begin();
    Star_iterator sitend = m_stars.end();
    for(; sit!=sitend; ++sit)
    {
      Star_handle si = get_star(sit);

      Index nsi = si->index_in_star_set();
      if(nsi == p1 || nsi == p2 || nsi == p3)
        continue;

      if(star->is_facet_encroached(transform_to_star_point(si->center_point(), star), facet))
      {
#ifdef BLABLA //todo
        std::cout << "Facet ";
        std::cout << facet.first->vertex((facet.second + 1)%4)->info() << " ";
        std::cout << facet.first->vertex((facet.second + 2)%4)->info() << " ";
        std::cout << facet.first->vertex((facet.second + 3)%4)->info() << " ";
        std::cout << "is encroached by the star: " << si->index_in_star_set() << std::endl;
#endif
        return true;
      }
    }
    return false;
  }

  Point_3 barycenter(const Facet& f) const
  {
    Point_3 p1 = m_stars[f.first->vertex((f.second + 1) % 4)->info()]->center_point();
    Point_3 p2 = m_stars[f.first->vertex((f.second + 2) % 4)->info()]->center_point();
    Point_3 p3 = m_stars[f.first->vertex((f.second + 3) % 4)->info()]->center_point();
    FT third = 1./3.;
    return Point_3(third * (p1.x() + p2.x() + p3.x()),
                   third * (p1.y() + p2.y() + p3.y()),
                   third * (p1.z() + p2.z() + p3.z()));
  }

  FT sq_distance_to_surface(const Facet& f, const Star_handle s) const
  {
    //Point_3 steiner;
    //s->compute_dual_intersection(f, steiner);
    //Point_3 cc = transform_from_star_point(s->compute_circumcenter(f), s);
    //return CGAL::squared_distance(steiner, cc);

    return m_pConstrain->compute_sq_approximation(barycenter(f));
  }

  bool is_surface_facet(const Facet& f) const
  {
    return m_stars[f.first->vertex((f.second + 1) % 4)->info()]->is_surface_star()
        && m_stars[f.first->vertex((f.second + 2) % 4)->info()]->is_surface_star()
        && m_stars[f.first->vertex((f.second + 3) % 4)->info()]->is_surface_star();
  }

  std::size_t count_restricted_facets() const
  {
    typename std::set<Facet_ijk> facets;
    for(unsigned int i = 0; i < m_stars.size(); i++)
    {
      typename Star::Facet_set_iterator fit = m_stars[i]->begin_restricted_facets();
      typename Star::Facet_set_iterator fitend = m_stars[i]->end_restricted_facets();
      for(; fit != fitend; fit++)
        facets.insert(Facet_ijk(*fit));
    }
    return facets.size();
  }

//Basic distortion functions
  FT compute_distortion(const Facet& f) const
  {
    FT distortion = 1.;
    int index = f.second;
    Cell_handle c = f.first;
    for (int i = 0; i < 3; i++)
    {
      int i1 = (index + i + 1) % 4;
      int i2 = (index + (i + 1) % 3 + 1) % 4;
      distortion = (std::max)(distortion,
                              m_stars[c->vertex(i1)->info()]->metric().compute_distortion(
                     m_stars[c->vertex(i2)->info()]->metric()));
    }
    return distortion;
  }

  FT compute_distortion(const Cell_handle& c) const
  {
    FT distortion = 1.0;
    for(int i = 0; i < 4; i++)
    {
      for(int j = i + 1; j < 4; j++)
      {
        distortion = (std::max)(distortion,
                                m_stars[c->vertex(i)->info()]->metric().compute_distortion(
                       m_stars[c->vertex(j)->info()]->metric()));
      }
    }
    return distortion;
  }

  void debug_show_distortions() const
  {
    for(std::size_t i = 0; i < m_stars.size(); ++i)
    {
      Star_handle s = m_stars[i];
      if(!s->is_surface_star())
        continue;
      std::cout << "  " << i << " : ";
      typename Star::Facet_set_iterator fit = s->begin_restricted_facets();
      typename Star::Facet_set_iterator fend = s->end_restricted_facets();
      for(; fit != fend; ++fit)
      {
        std::cout << "(";
        Facet f = *fit;
        if(is_surface_facet(f))
        {
          for (int i = 0; i < 3; i++)
          {
            int index_1 = (f.second + i + 1) % 4;
            int index_2 = (f.second + (i + 1) % 3 + 1) % 4;
            FT distortion =
              m_stars[f.first->vertex(index_1)->info()]->metric().compute_distortion(
              m_stars[f.first->vertex(index_2)->info()]->metric());

            FT costheta = std::abs(
              m_stars[f.first->vertex(index_1)->info()]->metric().get_vmin()
              * m_stars[f.first->vertex(index_2)->info()]->metric().get_vmin());

            std::cout << (distortion - 1./costheta) << "\t";
          }
        }
        std::cout << ")";
      }
      std::cout << std::endl;
    }
  }

  double facet_distortion(Facet& f) const
  {
    FT max_distortion = 0.;
    for (int i = 0; i < 3; i++)
    {
      int index_1 = (f.second + i + 1) % 4;
      int index_2 = (f.second + (i + 1) % 3 + 1) % 4;
      FT distortion = m_stars[f.first->vertex(index_1)->info()]->metric().compute_distortion(
                        m_stars[f.first->vertex(index_2)->info()]->metric());
      max_distortion = (std::max)(distortion, max_distortion);
    }

    return max_distortion;
  }

  double average_facet_distortion(bool verbose = true) const
  {
    double avg_coh_dis = 0., min_coh_dis = 1e30, max_coh_dis = -1e30;
    double avg_incoh_dis = 0., min_incoh_dis = 1e30, max_incoh_dis = -1e30;
    int nb_coh = 0, nb_incoh = 0;

    std::set<Facet_ijk> done;
    for(std::size_t i = 0; i <number_of_stars(); i++)
    {
      FT max_distortion = 0.;
      Star_handle star = m_stars[i];
      if(!star->is_surface_star())
        continue;

      typename Star::Facet_set_iterator fit = star->begin_restricted_facets();
      typename Star::Facet_set_iterator fitend = star->end_restricted_facets();
      for(; fit != fitend; fit++)
      {
        Facet f = *fit;
        std::pair<typename std::set<Facet_ijk>::iterator, bool> is_insert_successful;
        is_insert_successful = done.insert(Facet_ijk(f));
        if(!is_insert_successful.second)
          continue;

        max_distortion = (std::max)(facet_distortion(f), max_distortion);

        if(is_consistent(f))
        {
          nb_coh++;
          avg_coh_dis += max_distortion;
          if(max_distortion > max_coh_dis)
            max_coh_dis = max_distortion;
          if(max_distortion < min_coh_dis)
            min_coh_dis = max_distortion;
        }
        else
        {
          nb_incoh++;
          avg_incoh_dis += max_distortion;
          if(max_distortion > max_incoh_dis)
            max_incoh_dis = max_distortion;
          if(max_distortion < min_incoh_dis)
            min_incoh_dis = max_distortion;
        }
      }
    }

    if(verbose)
    {
      std::cout << "SUM UP OF THE DISTORTION (FACET):" << std::endl;
      std::cout << nb_coh << " coherent facets with avg: " << avg_coh_dis/(double) nb_coh << std::endl;
      std::cout << "min: " << min_coh_dis << " max: " << max_coh_dis << std::endl;
      std::cout << nb_incoh << " incoherent facets with avg " << avg_incoh_dis/(double) nb_incoh << std::endl;
      std::cout << "min: " << min_incoh_dis << " max: " << max_incoh_dis << std::endl;
      std::cout << "---------------------------------------------" << std::endl;
    }

    return avg_coh_dis/(double) nb_coh;
  }

  double star_distortion(Star_handle star) const
  {
    FT max_distortion = 0.;

    /*
  typename Star::Facet_set_iterator fit = star->begin_restricted_facets();
  typename Star::Facet_set_iterator fitend = star->end_restricted_facets();
  for(; fit != fitend; fit++)
  {
    Facet f = *fit;
    max_distortion = (std::max)(max_distortion, facet_distortion(f));
  }
  return max_distortion;
  */

    typename std::vector<Vertex_handle>::iterator vit = star->begin_neighboring_vertices();
    typename std::vector<Vertex_handle>::iterator vend = star->end_neighboring_vertices();
    for(; vit!=vend; ++vit)
    {
      if(is_infinite_vertex((*vit)->info()))
        continue;
      if(!(m_stars[(*vit)->info()]->is_surface_star())) //todo this is bad
        continue;

      FT distortion = m_stars[(*vit)->info()]->metric().compute_distortion(star->metric());
      max_distortion = (std::max)(distortion, max_distortion);
    }

    return max_distortion;
  }

  double average_star_distortion(bool verbose = true) const
  {
    double avg_coh_dis = 0., min_coh_dis = 1e30, max_coh_dis = -1e30;
    double avg_incoh_dis = 0., min_incoh_dis = 1e30, max_incoh_dis = -1e30;
    int nb_coh = 0, nb_incoh = 0;

    for(std::size_t i=0; i<number_of_stars(); ++i)
    {
      Star_handle si = get_star(i);
      if(!(si->is_surface_star()))
        continue;

      FT max_distortion = star_distortion(si);

      if(is_consistent(si))
      {
        nb_coh++;
        avg_coh_dis += max_distortion;
        if(max_distortion > max_coh_dis)
          max_coh_dis = max_distortion;
        if(max_distortion < min_coh_dis)
          min_coh_dis = max_distortion;
      }
      else
      {
        nb_incoh++;
        avg_incoh_dis += max_distortion;
        if(max_distortion > max_incoh_dis)
          max_incoh_dis = max_distortion;
        if(max_distortion < min_incoh_dis)
          min_incoh_dis = max_distortion;
      }
    }

    if(verbose)
    {
      std::cout << "SUM UP OF THE DISTORTION (STAR) :" << std::endl;
      std::cout << nb_coh << " coherent (surface) stars with avg: " << avg_coh_dis/(double) nb_coh << std::endl;
      std::cout << "min: " << min_coh_dis << " max: " << max_coh_dis << std::endl;
      std::cout << nb_incoh << " incoherent (surface) stars with avg " << avg_incoh_dis/(double) nb_incoh << std::endl;
      std::cout << "min: " << min_incoh_dis << " max: " << max_incoh_dis << std::endl;
      std::cout << "---------------------------------------------" << std::endl;
    }

    return avg_coh_dis/(double) nb_coh;
  }

//aabb functions
public:
  void update_aabb_tree(Star_handle star) const
  {
    m_aabb_tree.update_primitive(AABB_primitive(star));
  }

  void build_aabb_tree()
  {
    m_aabb_tree.rebuild(m_stars.begin(), m_stars.end());
  }

  void update_bboxes() const
  {
#ifdef DEBUG_UPDATE_AABB_TREE
    std::cout << "updating bboxes. tree : " << m_aabb_tree.size() << " " << m_aabb_tree.m_insertion_buffer_size() << std::endl;
#endif
    std::size_t i;
    std::size_t N = m_stars.size();
    for(i = 0; i < N; i++)
    {
      m_stars[i]->update_bbox();
#ifndef NO_USE_AABB_TREE_OF_BBOXES
      if(m_aabb_tree.m_insertion_buffer_size() == 1) // tree has just been rebuilt
        m_stars[i]->bbox_needs_aabb_update() = false;
      if(m_stars[i]->bbox_needs_aabb_update() && i<m_aabb_tree.size())
      {
        m_stars[i]->bbox_needs_aabb_update() = false;
        update_aabb_tree(m_stars[i]);
      }
#endif
    }

#ifndef NO_USE_AABB_TREE_OF_BBOXES
    for(i = 0; i < N; i++)
      if(m_stars[i]->bbox_needs_aabb_update() && i<m_aabb_tree.size())
        std::cout << "forgot some stars in the update" << std::endl;
#endif
  }

  template<typename OutputIterator>
  void finite_stars_in_conflict(const Point_3& p,
                                OutputIterator oit) const
  {
    update_bboxes();
#ifndef NO_USE_AABB_TREE_OF_BBOXES
    //bigger set of stars
    m_aabb_tree.all_intersected_primitives(p, oit);
#else
    //exact set
    for(unsigned int i = 0; i < m_stars.size(); i++)
    {
      Star_handle s = m_stars[i];
      Cell_handle ch;
      if(s->is_conflicted(transform_to_star_point(p, s), ch))
        *oit++ = s;
    }
#endif
  }

  template<typename OutputIterator>
  bool infinite_stars_in_conflict(const Point_3& p,
                                  OutputIterator oit,
                                  const bool collect = true) const
  {
    if(m_ch_triangulation.dimension() < 3)
    {
      all_stars(oit);
      return true;
    }

    int li, lj;
    typename DT::Locate_type lt;
    typename DT::Cell_handle c = m_ch_triangulation.locate(p, lt, li, lj);

    if(lt != DT::OUTSIDE_CONVEX_HULL && lt != DT::OUTSIDE_AFFINE_HULL)
      return false;

    std::vector<typename DT::Facet> bfacets;
    m_ch_triangulation.find_conflicts(p, c, std::back_inserter(bfacets), Emptyset_iterator());

    if(!collect)
       return !bfacets.empty();

    typename std::vector<typename DT::Facet>::iterator fit = bfacets.begin();
    typename std::vector<typename DT::Facet>::iterator fitend = bfacets.end();
    for(; fit != fitend; fit++)
    {
      typename DT::Facet f = *fit;
      if(m_ch_triangulation.is_infinite(f))
        continue;

      *oit++ = m_stars[f.first->vertex((f.second + 1) % 4)->info()];
      *oit++ = m_stars[f.first->vertex((f.second + 2) % 4)->info()];
      *oit++ = m_stars[f.first->vertex((f.second + 3) % 4)->info()];
    }
    return true;
  }

  template<typename StarConstIterator>
  void remove_from_stars(const Index& id,
                         StarConstIterator begin,
                         StarConstIterator end)
  {
    StarConstIterator it;
    for(it = begin; it != end; it++)
    {
      Star_handle s = get_star(it);
      if(id != s->index_in_star_set())
        s->remove(id);
    }
  }

  void pop_back_star()
  {
    Star_handle last = m_stars.back();
    delete last;

    m_stars.pop_back();
    Index id = static_cast<Index>(m_stars.size());
    remove_from_stars(id, m_stars.begin(), m_stars.end());
    m_kd_tree.remove_last();
#ifndef NO_USE_AABB_TREE_OF_BBOXES
    m_aabb_tree.remove_last();
#endif
  }

  void clean_stars()
  {
    std::cout << "clean call" << std::endl;
    for(std::size_t i = 0; i < m_stars.size(); i++)
      m_stars[i]->clean();
  }

//Point insertion functions
  Index simulate_insert_to_stars(const Point_3& p,
                                 Index_set& modified_stars) const
  {
    Index this_id = static_cast<Index>(m_stars.size());

    // find conflicted stars
    Star_set stars;
    finite_stars_in_conflict(p, std::inserter(stars, stars.end())); //aabb tree
    infinite_stars_in_conflict(p, std::inserter(stars, stars.end()));//convex hull

    typename Star_set::iterator it = stars.begin();
    typename Star_set::iterator itend = stars.end();
    for(; it != itend; it++)
    {
      Star_handle si = get_star(it);
      int id = si->simulate_insert_to_star(p, this_id);

      if(id == -1) // no conflict
        continue;
      else if(id < (int)this_id) // already in star set
        return id;
      else // to be inserted, standard configuration
        modified_stars.insert(si->index_in_star_set());
    }
    return this_id;
  }

  void create_star(const Point_3 &p,
                   int pid,
                   const Index_set& modified_stars,//by p's insertion
                   Star_handle& star,
                   const bool surface_star = true,
                   const bool metric_reset = true) const
  {
#ifdef USE_ANISO_TIMERS
    std::clock_t start_time = clock();
#endif
    if(1) //metric_reset
    { //todododododdodododdotodododododdodododdotodododododdodododdotodododododdodododdo
      if(1) //surface_star
      {
        Metric m_p;

        /*
        if(m_use_c3t3_colors || m_compute_metric_from_triangle)
        {
          Vector_3 normal;
          std::vector<Vector_3> v(3);
          std::vector<FT> e(3);
          Eigen::Matrix3d M;
          if(m_use_c3t3_colors)
            m_poly_painter.get_color_from_poly(p, M);
          else
            M = compute_metric_from_index_set(p, modified_stars);

          get_eigen_vecs_and_vals<K>(M, v[0], v[1], v[2], e[0], e[1], e[2]);
          int n_ind = -1;
          m_p.get_third_eigenvector(normal);
          get_metric_normal_index<K>(normal, v, n_ind);
          int id_1 = (n_ind + 1) % 3;
          int id_2 = (n_ind + 2) % 3;
          m_p = metric_field()->build_metric(v[n_ind], v[id_1], v[id_2],
                                             std::sqrt(e[n_ind]),
                                             std::sqrt(e[id_1]),
                                             std::sqrt(e[id_2]));
        }
        else //standard surface case
        */
          m_p = metric_field()->compute_metric(p);

        star->reset(p, pid, m_p, surface_star);
      }
      else
        star->reset(p, pid, metric_field()->uniform_metric(p), surface_star);
    }

    typename Index_set::const_iterator si = modified_stars.begin();
    typename Index_set::const_iterator siend = modified_stars.end();
    for (; si != siend; si++)
    {
      Star_handle star_i = get_star(si);
      star->insert_to_star(star_i->center_point(), star_i->index_in_star_set(), false);
        //no condition because they should be there for consistency
    }

    typename Star::Bbox bbox = star->bbox(); // volume bbox + surface bbox when star is not a topo_disk
    Point_3 pmin(bbox.xmin(), bbox.ymin(), bbox.zmin());
    Point_3 pmax(bbox.xmax(), bbox.ymax(), bbox.zmax());
    Kd_Box_query query(pmin, pmax, /*3=dim,*/ 0./*epsilon*/, typename Kd_tree::Star_pmap(m_stars));
    std::set<Kd_point_info> indices;
    m_kd_tree.search(std::inserter(indices, indices.end()), query);

    if(indices.size() != modified_stars.size())// because 'indices' contains 'modified_stars'
    {
      Index_set diff;
      std::set_difference(indices.begin(), indices.end(),
                          modified_stars.begin(), modified_stars.end(), std::inserter(diff, diff.end()));
      typename Index_set::iterator it = diff.begin();
      while(it != diff.end())
      {
        Star_handle si = get_star(it++);
        star->insert_to_star(si->center_point(), si->index_in_star_set(), true/*conditional*/);
      }
    }
  }

  Star_handle create_star(const Point_3 &p,
                          int pid,
                          const Index_set& modified_stars)
  {
    Star_handle star = new Star(m_criteria, m_pConstrain, true /*surface star*/);
    create_star(p, pid, modified_stars, star, true /*surface star*/);
    return star;
  }

  Star_handle create_inside_star(const Point_3 &p,
                                 int pid,
                                 const Index_set& modified_stars) const //by p's insertion
  {
    Star_handle star = new Star(m_criteria, m_pConstrain, false/*surface*/);
    create_star(p, pid, modified_stars, star, false/*surface*/);
    return star;
  }

  // inserts p in target_stars ('conditional' to be in conflict or not)
  // modified_stars get filled with the stars in which p is actually inserted
  // triangulation of CH is updated accordingly
  Index perform_insertions(const Point_3& p,
                           const Index& this_id,
                           const Star_set& target_stars,
                           Index_set& modified_stars,
                           const bool conditional,
                           const bool infinite_stars = false)
  {
    Vertex_handle v_in_ch = Vertex_handle();
    typename Star_set::const_iterator it = target_stars.begin();
    typename Star_set::const_iterator itend = target_stars.end();
    for(; it != itend; it++)
    {
      Star_handle si = *it;
      Vertex_handle vi = si->insert_to_star(p, this_id, conditional);

      if(vi == Vertex_handle())  //no conflict
        continue;
      else if(vi->info() < this_id) // already in star set (should not happen)
      {
        std::cout << "Warning! Insertion of p"<< this_id
                  << " (" << p << ") in S" << si->index_in_star_set()
                  << " failed. vi->info() :"<< vi->info() << std::endl;
        if(v_in_ch != Vertex_handle())
          m_ch_triangulation.remove(v_in_ch);
        remove_from_stars(this_id, target_stars.begin(), ++it);

        si->print_vertices(true);
        si->print_faces();
        std::cout << "Metric : \n" << si->metric().get_transformation() << std::endl;
        return vi->info();
      }
      else // inserted, standard configuration
      {
        modified_stars.insert(si->index_in_star_set());
        // update triangulation of convex hull
        if(infinite_stars)//(si->is_infinite() && v_in_ch != Vertex_handle())
        {
          v_in_ch = m_ch_triangulation.insert(p);
          v_in_ch->info() = this_id;
        }
      }
    }
    return this_id;
  }

  Index insert_to_stars(const Point_3& p,
                        Index_set& modified_stars,
                        const bool conditional)
  {
    Index this_id = static_cast<Index>(m_stars.size());

    // find conflicted stars & insert p to these stars
    Index id;
    Star_set target_stars;
    if(conditional)
      finite_stars_in_conflict(p, std::inserter(target_stars, target_stars.end())); // aabb tree/exhaustive
    else
      all_stars(std::inserter(target_stars, target_stars.end()));

    //std::cout << target_stars.size() << " stars found by the conflict find" << std::endl;
    id = perform_insertions(p, this_id, target_stars, modified_stars, conditional, false);

    target_stars.clear();
    infinite_stars_in_conflict(p, std::inserter(target_stars, target_stars.end()));//convex hull
    id = perform_insertions(p, this_id, target_stars, modified_stars, conditional, true);

    return id;
  }

  Index insert(const Point_3& p, const bool conditional)
  {
    Index_set modified_stars;
    return insert(p, modified_stars, conditional);
  }

  Index insert(const Point_3 &p,
               Index_set& modified_stars,
               const bool conditional)
  {
    Index id = insert_to_stars(p, modified_stars, conditional);
    if(id < 0 || id < (int)m_stars.size())
      return id;

    Star_handle star = create_star(p, id, modified_stars); // implies surface star

    if(star->index_in_star_set() != m_stars.size())
      std::cout << "WARNING in insert..." << std::endl;

    m_stars.push_back(star);
    modified_stars.insert(star->index_in_star_set());
    m_kd_tree.insert(star->index_in_star_set());
#ifndef NO_USE_AABB_TREE_OF_BBOXES
    m_aabb_tree.insert(AABB_primitive(star));
#endif
    return id;
  }

  Index insert_in_domain(const Point_3& p,
                         Index_set& modified_stars,
                         const bool conditional)
  {
    Index id = insert_to_stars(p, modified_stars, conditional);

    if(id < 0 || id < (int)m_stars.size())
      return id;

    Star_handle star = create_inside_star(p, id, modified_stars);
    m_stars.push_back(star);
    modified_stars.insert(star->index_in_star_set());
    m_kd_tree.insert(star->index_in_star_set());
#ifndef NO_USE_AABB_TREE_OF_BBOXES
    m_aabb_tree.insert(AABB_primitive(star));
#endif
    return id;
  }


//Constructors
  Anisotropic_refine_trunk(Star_vector& m_stars_,
                           const Constrain_surface* m_pConstrain_,
                           const Criteria* m_criteria_,
                           const Metric_field* m_metric_field_,
                           DT& m_ch_triangulation_,
                           AABB_tree& m_aabb_tree_,
                           Kd_tree& m_kd_tree_)
  :
    m_stars(m_stars_),
    m_pConstrain(m_pConstrain_),
    m_criteria(m_criteria_),
    m_metric_field(m_metric_field_),
    m_ch_triangulation(m_ch_triangulation_),
    m_aabb_tree(m_aabb_tree_),
    m_kd_tree(m_kd_tree_)
  { }

};  // Anisotropic_refine_trunk

}  // Anisotropic_mesh_3
}  // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_REFINE_trunk_H

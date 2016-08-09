#ifndef CGAL_ANISOTROPIC_MESH_2_REFINE_TRUNK_H
#define CGAL_ANISOTROPIC_MESH_2_REFINE_TRUNK_H

//#define ANISO_DEBUG_CZONES
#define ANISO_USE_BOUNDING_BOX_VERTICES_AS_POLES

#include <CGAL/Criteria.h>
#include <CGAL/Domain_2.h>
#include <CGAL/Kd_tree_for_star_set.h>
#include <CGAL/Metric_field.h>
#include <CGAL/Stretched_Delaunay_2.h>
#include <CGAL/Starset.h>

#include <CGAL/Kd_tree_for_star_set.h>
#include <CGAL/aabb_tree/aabb_tree_bbox.h>
#include <CGAL/aabb_tree/aabb_tree_bbox_primitive.h>

#include <CGAL/IO/Star_set_output.h>

#include <omp.h>

#include <deque>
#include <set>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_2
{

//these classes are used to keep in memory the various conflict zones for each modified star.
//This is done for multiple reasons :
// * Avoid computing it multiple times (encroachment tests, pick_valid tests, point insertion...)
// * Keep in memory the faces/cells that will be destroyed to be able to track elements that might
//   need to be inserted in the pqueues (especially those not in modified stars)
template<typename K>
class Conflict_zone
{
  typedef Conflict_zone<K>                       Self;

public:
  typedef Stretched_Delaunay_2<K>                Star;
  typedef typename Star::Face                    Face;
  typedef typename Star::Face_handle             Face_handle;
  typedef typename Star::Face_handle_vector      Face_handle_vector;
  typedef typename Star::Edge                    Edge;
  typedef std::vector<Edge>                      Edge_vector;

private:
  //Conflict zones defining vectors. Valid only until the insertion of the point!
  Edge_vector m_boundary_edges;
  Edge_vector m_internal_edges;
  Face_handle_vector m_faces;

public:
  Edge_vector& boundary_edges() { return m_boundary_edges; }
  const Edge_vector& boundary_edges() const { return m_boundary_edges; }
  Edge_vector& internal_edges() { return m_internal_edges; }
  const Edge_vector& internal_edges() const { return m_internal_edges; }
  Face_handle_vector& faces() { return m_faces; }
  const Face_handle_vector& faces() const { return m_faces; }

  bool empty()
  {
    return m_boundary_edges.empty() && m_internal_edges.empty() && m_faces.empty();
  }

  //Constructors
  Conflict_zone() : m_boundary_edges(), m_internal_edges(), m_faces() { }
};

template<typename K, typename KExact = K>
class Stars_conflict_zones
{
  typedef Stars_conflict_zones<K, KExact>                  Self;

public:
  typedef Stretched_Delaunay_2<K, KExact>                  Star;
  typedef CGAL::Anisotropic_mesh_2::Starset<K, KExact>     Starset;
  typedef Star*                                            Star_handle;
  typedef std::vector<Star_handle>                         Star_vector;
  typedef CGAL::Anisotropic_mesh_2::Conflict_zone<K>       Czone;
  typedef typename Star::Index                             Index;
  typedef typename Star::Point_2                           Point;
  typedef typename KExact::Point_2                         Exact_Point;
  typedef typename Star::Vertex                            Vertex;
  typedef typename Star::Face                              Face;
  typedef typename Star::Face_handle                       Face_handle;
  typedef typename std::map<Index, Czone>::iterator        iterator;
  typedef typename std::map<Index, Czone>::const_iterator  const_iterator;

public:
  enum Conflict_zones_status
  {
    CONFLICTS_UNKNOWN = 0,
    CONFLICTING_POINT_IS_KNOWN,
    CONFLICTING_STARS_ARE_KNOWN, // conflicting_point is known too
    CONFLICT_ZONES_ARE_KNOWN // everything is known
  };

private:
  const Starset& m_starset;
  std::map<Index, Czone> m_conflict_zones;
  Point m_conflict_p; // point causing the conflict (only needed to debug / assert)
  Index m_conflict_p_id; //id of the point

  bool m_faces_need_check;
  bool m_are_checks_computed; //already computed or not

  //we need to keep track of the "soon to be" destroyed faces/cells SEEN FROM UNMODIFIED STARS
  //not lose consistency in the unmodified stars and add them to the queue if needed
  std::map<Index, std::vector<Face_handle> > m_faces_to_check;

public:
  const std::map<Index, Czone>& conflict_zones() const { return m_conflict_zones; }
  Czone& conflict_zone(Index i) { return m_conflict_zones[i]; }
  const Czone& conflict_zone(Index i) const { return m_conflict_zones[i]; }
  Point& conflicting_point() { return m_conflict_p; }
  const Point& conflicting_point() const { return m_conflict_p; }
  Index& conflicting_point_id() { return m_conflict_p_id; }
  const Index& conflicting_point_id() const { return m_conflict_p_id; }
  std::map<Index, std::vector<Face_handle> >& faces_to_check() { return m_faces_to_check; }
  const bool& are_checks_computed() const { return m_are_checks_computed; }

  const_iterator begin() const { return m_conflict_zones.begin(); }
  const_iterator end() const { return m_conflict_zones.end(); }
  iterator begin() { return m_conflict_zones.begin(); }
  iterator end() { return m_conflict_zones.end(); }
  std::size_t size() const { return m_conflict_zones.size(); }

  bool empty() const { return m_conflict_zones.empty(); }
  void clear()
  {
    m_conflict_zones.clear();
    m_conflict_p = Point(1e17, 1e17, 1e17); // todo change that to something better
    m_conflict_p_id = -1;
    m_faces_to_check.clear();
    m_are_checks_computed = false;
  }

  Conflict_zones_status status()
  {
    if(m_conflict_p == Point(1e17, 1e17, 1e17))
      return CONFLICTS_UNKNOWN;
    if(m_conflict_zones.empty())
      return CONFLICTING_POINT_IS_KNOWN;
    if((m_conflict_zones.begin()->second).empty())
      return CONFLICTING_STARS_ARE_KNOWN;
    return CONFLICT_ZONES_ARE_KNOWN;
  }

  bool are_check_maps_empty()
  {
    return m_faces_to_check.empty();
  }

  void check_from_face(Face_handle fit,
                        const std::vector<Index>& indices)
  { }

  void check_from_internal_edges(iterator mit,
                   const std::map<Edge_ij, std::size_t>& internal_edges_counter)
  { }

  void check_from_boundary_edges(iterator mit)
  { }

  void compute_elements_needing_check()
  {
    // code for this function is not finished so if you're calling insert_in_stars
    // with conditional = true, you're going to miss out on something inconsistencies
//    std::cerr << "compute_elements_needing_check() call, see comment" << std::endl;
    return;

//todo -----------------------
    //WHEN REFINEMENT_CONDITIONS ARE ADDED AGAIN THERE NEEDS TO BE A CHECK HERE
    //NOT TO CONSIDER USELESS FACES/CELLS
//-----------------------

    //What is being done here:
    //we need to make sure that no inconsistency appears in stars that are not
    //in conflict with the point. This might appear when a simplex is in conflict
    //in one star, but not in the others.

    if(m_are_checks_computed)
      return;

    //Count the elements. In the general case, an internal edge(face) is two(three)
    //times in conflict, and there is thus nothing to check.
    std::map<Edge_ij, std::size_t> internal_edges_counter;
    std::map<Facet_ijk, int> faces_counter;

    iterator mit = m_conflict_zones.begin();
    iterator mend = m_conflict_zones.end();
    for(; mit!=mend; ++mit)
    {
      Star_handle si = m_starset[mit->first];
      Czone& czi = mit->second;
      if(czi.empty())
        continue;

      // internal edges
      typename std::vector<typename Star::Edge>::iterator eit =
                                           czi.internal_edges().begin();
      typename std::vector<typename Star::Edge>::iterator eend =
                                           czi.internal_edges().end();
      for(; eit!=eend; ++eit)
      {
        const typename Star::Edge& e = *eit;
        if(e.first->vertex((e.second + 1)%3)->info() == mit->first || // incident edges only
           e.first->vertex((e.second + 2)%3)->info() == mit->first)
          add_to_map(Edge_ij(e), internal_edges_counter);
      }

      if(!m_faces_need_check)
        continue;

      // faces
      typename Star::Face_handle_handle fit = czi.faces().begin();
      typename Star::Face_handle_handle fend = czi.faces().end();
      for(; fit!=fend; ++fit)
      {
        Face_handle fh = *fit;
        if(si->is_inside(fh) &&
           (fh->vertex(0)->info() == mit->first || // incident cells only!
            fh->vertex(1)->info() == mit->first ||
            fh->vertex(2)->info() == mit->first))
          add_to_map(Facet_ijk(*fit), faces_counter);
      }
    }

    mit = m_conflict_zones.begin();
    for(; mit!=mend; ++mit)
    {
      check_from_internal_edges(mit, internal_edges_counter);
    }

    //TODO : std unique on both if not empty since we can have duplicates
    //       make sure we actually have duplicates and they are actually
    //       adjacent first..

#ifdef ANISO_DEBUG_UNMODIFIED_STARS
    if(!m_faces_to_check.empty())
    {
      std::cout << "faces_to_check is not empty. Size (n° of stars): ': " << m_faces_to_check.size() << std::endl;
      std::cout << "Detail: " << std::endl;

      typename std::map<Index, std::vector<face> >::iterator mvfit = m_faces_to_check.begin();
      typename std::map<Index, std::vector<face> >::iterator mvfend = m_faces_to_check.end();
      for(; mvfit!= mvfend; ++mvfit)
      {
        Index sid = mvfit->first;
        std::cout << "-- Star " << sid << std::endl;
        std::vector<face> vf = mvfit->second;

        face_handle fit = vf.begin();
        face_handle fend = vf.end();
        for(; fit!=fend; ++fit)
        {
          std::cout << fit->first->vertex((fit->second+1)%4)->info() << " ";
          std::cout << fit->first->vertex((fit->second+2)%4)->info() << " ";
          std::cout << fit->first->vertex((fit->second+3)%4)->info() << " || ";
          std::cout << fit->first->vertex(fit->second)->info() << std::endl;
        }
      }
    }
#endif
    m_are_checks_computed = true;
  }

//Constructors
  Stars_conflict_zones(const Starset& m_starset_)
    :
      m_starset(m_starset_),
      m_conflict_zones(),
      m_conflict_p(Point(1e17, 1e17, 1e17)),
      m_conflict_p_id(-1),
      m_faces_need_check(false),
      m_are_checks_computed(false),
      m_faces_to_check()
  { }
};

template<typename K>
inline std::ostream& operator<<(std::ostream& os, Stars_conflict_zones<K>& src)
{
  typedef typename Stars_conflict_zones<K>::Star     Star;
  typedef typename Star::Face_handle                 Face_handle;
  typedef typename Star::Face_handle_handle          Face_handle_handle;

  os << "printing current czones:" << std::endl;
  os << "point : " << src.conflicting_point() << std::endl;
  os << "id: " << src.conflicting_point_id() << std::endl;
  os << "size: " << src.size() << std::endl;
  os << "zones: " << std::endl;

  typename Stars_conflict_zones<K>::iterator it = src.begin();
  typename Stars_conflict_zones<K>::iterator iend = src.end();
  for(; it!=iend; ++it)
  {
    typename Stars_conflict_zones<K>::Czone& czi = it->second;
    os << "S " << it->first << " ";
//    os << czi.boundary_edges().size() << " ";
    os << czi.faces().size() << " ";
    os << std::endl;

    //    src.m_starset[it->first]->print_faces();

    os << "-*-*-*-*-" << std::endl;
    os << "-*-*-*-*-" << std::endl;
    if(1)
      continue;

/*
    typename Stars_conflict_zones<K>::Edge_handle eit = czi.boundary_edge().begin();
    typename Stars_conflict_zones<K>::Edge_handle eend = czi.boundary_edges().end();
    for(; eit!=eend; ++eit)
    {
      os << " " << eit->first->vertex((eit->second+1)%3)->info() << " ";
      os << eit->first->vertex((eit->second+2)%3)->info() << std::endl;
    }
    os << "-----" << std::endl;
*/

    Face_handle_handle fit = czi.faces().begin();
    Face_handle_handle fend = czi.faces().end();
    for(; fit!=fend; ++fit)
    {
      Face_handle fh = *fit;
      os << fh->vertex(0)->info() << " ";
      os << fh->vertex(1)->info() << " ";
      os << fh->vertex(2)->info() << std::endl;
    }
    os << "-*-*-*-*-" << std::endl;
    os << "-*-*-*-*-" << std::endl;
  }

  return os;
}

template<typename K>
class Anisotropic_refine_trunk
{
private:
  typedef Anisotropic_refine_trunk<K>                              Self;
public:
  //typedef CGAL::Exact_predicates_exact_constructions_kernel KExact;
  typedef K                                                        KExact;

  typedef Stretched_Delaunay_2<K>                                  Star;
  typedef typename Star::FT                                        FT;
  typedef Star*                                                    Star_handle;
  typedef typename Star::Base                                      DT; // DT_2 with vertex_base_with_info
  typedef std::vector<Star_handle>                                 Star_vector;
  typedef typename Star_vector::iterator                           Star_iterator;
  typedef std::set<Star_handle>                                    Star_set;
  typedef std::deque<Star_handle>                                  Star_deque;
  typedef typename Star::Index                                     Index;
  typedef std::set<Index>                                          Index_set;
  typedef typename Star::Point_2                                   Point_2;
  typedef typename Star::Point_3                                   Point_3;
  typedef typename Star::TPoint_2                                  TPoint_2;
  typedef std::set<Point_2>                                        Point_set;
  typedef typename Star::Edge                                      Edge;
  typedef typename Star::Vertex_handle                             Vertex_handle;
  typedef typename Star::Face_handle_handle                        Face_handle_handle;
  typedef typename Star::Face                                      Face;
  typedef typename Star::Face_handle                               Face_handle;
  typedef typename Star::Vector_2                                  Vector_2;
  typedef typename Star::Domain                                    Domain;
  typedef typename Star::Criteria                                  Criteria;

  typedef CGAL::Anisotropic_mesh_2::Starset<K>                     Starset;

  typedef CGAL::Anisotropic_mesh_2::Metric_field<K>                Metric_field;
  typedef typename Metric_field::Metric                            Metric;

  typedef CGAL::AABB_tree_bbox<K, Star>                            AABB_tree;
  typedef CGAL::AABB_bbox_primitive<Star>                          AABB_primitive;

  typedef CGAL::Kd_tree_for_star_set<K, Star_handle>               Kd_tree;
  typedef typename Kd_tree::Traits                                 Kd_traits;
  typedef typename Kd_tree::Box_query                              Kd_Box_query;
  typedef typename Kd_tree::key_type                               Kd_point_info;

  typedef CGAL::Anisotropic_mesh_2::Conflict_zone<K>               Conflict_zone;
  typedef CGAL::Anisotropic_mesh_2::Stars_conflict_zones<K>        Stars_conflict_zones;

private:
  typedef typename KExact::Point_2                                 Exact_Point_2;
  typedef typename KExact::Point_2                                 Exact_TPoint_2;
  typedef CGAL::Cartesian_converter<K, KExact>                     To_exact;
  typedef CGAL::Cartesian_converter<KExact, K>                     Back_from_exact;

  To_exact to_exact;
  Back_from_exact back_from_exact;

protected:
  Starset& m_starset;

  const Domain* m_pdomain;
  const Criteria* m_criteria;
  const Metric_field* m_metric_field;

  AABB_tree& m_aabb_tree;
  Kd_tree& m_kd_tree;

  Stars_conflict_zones& m_stars_czones;

public:
  //virtual functions used in the visitors
  //We could avoid the virtual by declaring all the different types of visitors manually...
  virtual const bool& is_active() const = 0;
  virtual void fill_refinement_queue(Index) = 0;

public:
  Star_vector& stars() const { return m_starset.star_vector(); }
  const Domain* domain() const { return m_pdomain; }
  const Criteria* criteria() const { return m_criteria; }
  const Metric_field* metric_field() const { return m_metric_field; }

  bool empty() const { return m_starset.empty(); }
  bool is_infinite_vertex(int i) const { return (i == Star::infinite_vertex_index()); }
  std::size_t number_of_stars() const { return m_starset.size(); }

  Star_handle get_star(Star_handle s) const { return s; }
  Star_handle get_star(std::size_t i) const { return m_starset[i]; }
  Star_handle get_star(int i) const         { return m_starset[i]; }
  Star_handle get_star(typename Index_set::const_iterator it) const   { return m_starset[*it]; }
  Star_handle get_star(typename Star_set::const_iterator it) const    { return *it; }
  Star_handle get_star(typename Star_vector::const_iterator it) const { return *it; }
  Star_handle get_star(typename Stars_conflict_zones::const_iterator it) const { return m_starset[it->first]; }

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
    for(std::size_t i = 0; i < number_of_stars(); i++)
      *oit++ = m_starset[i];
  }

  double duration(const time_t& start) const
  {
    return ((clock() - start + 0.) / ((double)CLOCKS_PER_SEC));
  }

public:
  TPoint_2 transform_to_star_point(const Point_2& p, Star_handle star) const
  {
    return star->metric().transform(p);
  }
  Point_2 transform_from_star_point(const TPoint_2& p, Star_handle star) const
  {
    return star->metric().inverse_transform(p);
  }

  bool is_inside_domain(const Point_2& p) const
  {
    return (m_pdomain->side_of_constraint(p) == ON_POSITIVE_SIDE);
  }

  Point_2 compute_circumcenter(const Face_handle& fh, Star_handle star) const
  {
    return transform_from_star_point(star->compute_circumcenter(fh), star);
  }

  Point_2 compute_circumcenter(const Point_2& p0, const Point_2& p1,
                               const Point_2& p2, Star_handle here) const
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

  Point_2 compute_circumcenter(const Point_2& p0, const Point_2& p1,
                               const Point_2& p2, const Point_2& p3,
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

  Point_2 barycenter(const Face_handle& fh) const
  {
    Point_2 p1 = get_star(fh->vertex(0)->info())->center_point();
    Point_2 p2 = get_star(fh->vertex(1)->info())->center_point();
    Point_2 p3 = get_star(fh->vertex(2)->info())->center_point();
    FT third = 1./3.;
    return Point_2(third * (p1.x() + p2.x() + p3.x()),
                   third * (p1.y() + p2.y() + p3.y()),
                   third * (p1.z() + p2.z() + p3.z()));
  }

public:
  //aabb functions
  void update_aabb_tree(Star_handle star) const
  {
    m_aabb_tree.update_primitive(AABB_primitive(star));
  }

  void build_aabb_tree()
  {
    std::cout << "build star set with " << m_starset.size() << " stars" << std::endl;
    m_aabb_tree.rebuild(m_starset.begin(), m_starset.end());
  }

  void update_bboxes() const
  {
    std::size_t i;
    std::size_t N = number_of_stars();
    for(i = 0; i < N; i++)
    {
      Star_handle si = get_star(i);
      si->update_bbox();
#ifndef NO_USE_AABB_TREE_OF_BBOXES
      if(m_aabb_tree.m_insertion_buffer_size() == 1) // tree has just been rebuilt
        si->bbox_needs_aabb_update() = false;
      if(si->bbox_needs_aabb_update() && i<m_aabb_tree.size())
      {
        si->bbox_needs_aabb_update() = false;
        update_aabb_tree(si);
      }
#endif
    }

#ifndef NO_USE_AABB_TREE_OF_BBOXES
    for(i = 0; i < N; i++)
      if(get_star(i)->bbox_needs_aabb_update() && i<m_aabb_tree.size())
        std::cout << "forgot some stars in the update" << std::endl;
#endif
  }

  template<typename OutputIterator>
  void finite_stars_in_conflict(const Point_2& p,
                                OutputIterator oit) const
  {
    // computes the exact set of stars in conflict

    update_bboxes();

#ifndef NO_USE_AABB_TREE_OF_BBOXES
  FT x = p.x(), y = p.y();
  #ifdef LIFT_AABB_TREE
    const Point_3 query(x, y, x);
  #else
    const Point_3 query(x, y, 0);
  #endif
    m_aabb_tree.all_intersected_primitives(query, oit);
#else
    const std::size_t ss = number_of_stars();
    for(std::size_t i = 0; i<ss; i++)
    {
      Star_handle s = get_star(i);
      Face_handle fh;
      if(s->is_conflicted(transform_to_star_point(p, s), fh))
        *oit++ = s;
    }
#endif
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

  void remove_from_stars(const Index& id)
  {
    Star_iterator sit = m_starset.begin();
    Star_iterator sitend = m_starset.end();
    remove_from_stars(id, sit, sitend);
  }

  void pop_back_star()
  {
    Star_handle last = m_starset.back();
    delete last;

    m_starset.pop_back();
    Index id = static_cast<Index>(number_of_stars());
    remove_from_stars(id, m_starset.begin(), m_starset.end());
    m_kd_tree.remove_last();
#ifndef NO_USE_AABB_TREE_OF_BBOXES
    m_aabb_tree.remove_last();
#endif
  }

  void clean_stars()
  {
    std::cout << "clean call, size: " << number_of_stars() << std::endl;
    for(std::size_t i = 0; i < number_of_stars(); i++)
      get_star(i)->clean();
  }

// Pick valid point validity checks for 3D (used by cell level
// and facet level when the cell level is active)
  void facets_created(Star_handle star,
                     std::map<Facet_ijk, int>& facets) const
  {
    Face_handle_handle fit = star->finite_incident_faces_begin();
    Face_handle_handle fiend = star->finite_incident_faces_end();
    for(; fit != fiend; fit++)
      add_to_map(Facet_ijk(*fit), facets);
  }

  // finite facets(i,j,k) which would be created by the insertion
  // of 'p' in 'star' are added to 'facets'
  void facets_created(const Point_2& p, // new point, not yet in the star set
                      const int p_id,   // index of p in star set
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
    const std::vector<Edge>& bedges = star_cz.boundary_edges();

    typename std::vector<Edge>::const_iterator bit = bedges.begin();
    typename std::vector<Edge>::const_iterator bend = bedges.end();
    for( ; bit!=bend; bit++)
    {
      int id1 = bit->first->vertex((bit->second + 1)%3)->info();
      int id2 = bit->first->vertex((bit->second + 2)%3)->info();
      int sid = star->index_in_star_set();

      if(id1 == star->infinite_vertex_index() ||
         id2 == star->infinite_vertex_index())
        continue;

      if(id1 != sid && id2 != sid)
        continue;

      Point_2 c = compute_circumcenter(this->m_starset[id1]->center_point(),
                                       this->m_starset[id2]->center_point(),
                                       p, star);
      if(is_inside_domain(c))
        add_to_map(Facet_ijk(id1, id2, p_id), facets);
    }
  }

  bool check_consistency_and_sliverity(Star_handle to_be_refined,
                                       const Star_handle& new_star,
                                       FT sq_radius_bound) const
  {

    Point_2 p = new_star->center_point();
    int p_index = new_star->index_in_star_set();

    std::map<Facet_ijk, int> facets;
    facets_created(new_star, facets);

    typename Stars_conflict_zones::iterator czit = this->m_stars_czones.begin();
    typename Stars_conflict_zones::iterator czend = this->m_stars_czones.end();
    for(; czit!=czend; ++czit)
    {
      Index i = czit->first;
      Star_handle si = get_star(i);
      facets_created(p, p_index, si, facets);
    }

    typename std::map<Facet_ijk, int>::iterator itc;
    for(itc = facets.begin(); itc != facets.end(); ++itc)
    {
      int nmax = number_of_stars();
      int c0 = (*itc).first.vertex(0);
      int c1 = (*itc).first.vertex(1);
      int c2 = (*itc).first.vertex(2);
      if((*itc).second != 3 ) // inconsistency
      {
        if((*itc).first.is_infinite())
          return false;

        TPoint_2 tp0, tp1, tp2;
        TPoint_2 tp = this->transform_to_star_point(p, to_be_refined);
        tp0 = (c0 == nmax) ? tp : this->transform_to_star_point(this->m_starset[c0]->center_point(),
                                                                to_be_refined);
        tp1 = (c1 == nmax) ? tp : this->transform_to_star_point(this->m_starset[c1]->center_point(),
                                                                to_be_refined);
        tp2 = (c2 == nmax) ? tp : this->transform_to_star_point(this->m_starset[c2]->center_point(),
                                                                to_be_refined);

        double sqr = to_be_refined->compute_squared_circumradius(tp0, tp1, tp2);
        if(sqr < sq_radius_bound)
        {
          return false; // small inconsistency (radius) is forbidden.
        }
      }
    }

    return true;
  }

  bool is_valid_point_2D(const Point_2 &p,
                         const FT sq_radiusbound,
                         Star_handle to_be_refined,
                         Star_handle& new_star) const
  {
    Index id = compute_conflict_zones(p);
    this->m_stars_czones.compute_elements_needing_check();

    if(!this->m_stars_czones.are_check_maps_empty())
    {
      //std::cout << "proposed CPV point creates inconsistencies in unmodified stars" << std::endl;
      return false;
    }

    if(id < 0) // No conflict
      return false;
    else if(id < (int) this->m_starset.size()) //already in star set
      return false;

    create_star(p, id, new_star);

    bool is = check_consistency_and_sliverity(to_be_refined, new_star, sq_radiusbound);
    if(!is)
      new_star->invalidate();

    return is;
  }

// Conflict zones
  Index simulate_insert_to_stars(const Point_2& p) const
  {
    const Index this_id = static_cast<Index>(number_of_stars());

    // find conflicted stars
    Star_deque stars;
    finite_stars_in_conflict(p, std::inserter(stars, stars.end()));

//    bool abort = false;
    std::size_t ds = stars.size();
//#pragma omp parallel for shared (abort) // omp does not accelerate (on the contrary)
    for(std::size_t i=0; i<ds; i++)
    {
//#pragma omp flush (abort) // make sure 'abort' is up to date
//      if(!abort)
//      {
        Star_handle si = get_star(stars[i]);
        const int id = si->simulate_insert_to_star(p, this_id);

        if(id == -1) // no conflict
          continue;
        else if(id < static_cast<int>(this_id)) // already in star set
        {
//#pragma omp critical
{
          std::cout << "Warning in simulate_insert_to_stars" << std::endl;
//          abort = true;
          return id;
}
        }
        else // to be inserted, standard configuration
        {
//#pragma omp critical
{
            //add them to the map of stars in conflict (conflict zones are not computed yet)
            m_stars_czones.conflict_zone(si->index_in_star_set());
}
        }
//      }
    }

    return this_id;
  }

  //those functions maybe could (should) be in the Czones class. But then it means
  //that all the filters (aabb, kd, etc.) need also a ref in the czone class(es) TODO (?)
  Index compute_conflict_zones(const Point_2& p) const
  {
    if(m_stars_czones.status() == Stars_conflict_zones::CONFLICT_ZONES_ARE_KNOWN
       && p == m_stars_czones.conflicting_point())
    {
      return m_stars_czones.conflicting_point_id();
    }

    if(m_stars_czones.status() == Stars_conflict_zones::CONFLICTS_UNKNOWN)
      m_stars_czones.conflicting_point() = p;

    if(p != m_stars_czones.conflicting_point())
    {
      std::cout << "Points should have been equal in compute_conflict_zones." << std::endl;
      std::cout << p << " " << m_stars_czones.conflicting_point() << std::endl;
    }

    if(m_stars_czones.status() == Stars_conflict_zones::CONFLICTING_POINT_IS_KNOWN)
      m_stars_czones.conflicting_point_id() = simulate_insert_to_stars(p);

    typename Stars_conflict_zones::iterator czit = m_stars_czones.begin();
    typename Stars_conflict_zones::iterator czend = m_stars_czones.end();
    for(; czit!=czend; ++czit)
    {
      Index i = czit->first;
      Star_handle si = get_star(i);
      Conflict_zone& ci = czit->second;

      si->find_conflicts(p, std::back_inserter(ci.faces()),
                            std::back_inserter(ci.boundary_edges()));
      CGAL_assertion(!ci.faces().empty());
    }

    for(std::size_t i=0; i<number_of_stars(); ++i)
    {
      Star_handle si = get_star(i);
      TPoint_2 tp = si->metric().transform(p);
      Face_handle fh;
      if(si->is_conflicted(tp, fh))
        CGAL_assertion(m_stars_czones.conflict_zones().find(i) !=
                       m_stars_czones.conflict_zones().end());
      else
        CGAL_assertion(m_stars_czones.conflict_zones().find(i) ==
            m_stars_czones.conflict_zones().end());

    }

    return m_stars_czones.conflicting_point_id();
  }

  void clear_conflict_zones() const
  {
    m_stars_czones.clear();
  }

  void insert_from_kd_tree(Star_handle star) const
  {
    int counter = 0;
    bool grew = true;
    while(grew)
    {
      counter++;

      const typename Star::Bbox& bbox = star->bbox();
      Point_2 pmin(bbox.xmin(), bbox.ymin());
      Point_2 pmax(bbox.xmax(), bbox.ymax());

      Kd_Box_query query(pmin, pmax, /*3=dim,*/ 0./*epsilon*/,
                         typename Kd_tree::Star_pmap(m_starset.star_vector()));
      std::set<Kd_point_info> indices;
      m_kd_tree.search(std::inserter(indices, indices.end()), query);

#ifdef ANISO_DEBUG_KDTREE
      std::cout << "bbox given to kd tree: " << std::endl;
      std::cout << pmin << std::endl << pmax << std::endl;
      std::cout << "kd tree finds : " << indices.size() << " pts out of "
                << this->number_of_stars() << std::endl;
      std::cout << "Brute force kd check...";
      Star_iterator sit = m_starset.begin();
      Star_iterator siend = m_starset.end();
      for(;sit!=siend;++sit)
        if(query.contains((*sit)->index_in_star_set()) &&
           indices.find((*sit)->index_in_star_set()) == indices.end())
        {
          std::cout << "the kd tree missed the star: "
                    << (*sit)->index_in_star_set() << std::endl;
        }
      std::cout << "end of brute force check" << std::endl;
#endif

      int old_nv = star->number_of_vertices();

      typename Index_set::iterator it = indices.begin();
      typename Index_set::iterator iend = indices.end();
      for(; it!=iend; ++it)
      {
        Star_handle si = get_star(it);
        star->insert_to_star(si->center_point(), si->index_in_star_set(),
                             true/*conditional*/);
      }

      int new_nv = star->number_of_vertices();
      grew = (old_nv != new_nv);
    }
  }

  void create_star_from_all_stars(Star_handle& star) const
  {
    // Increasingly slow but will produce the correct star
    typename Starset::iterator sit = m_starset.begin();
    typename Starset::iterator send = m_starset.end();
    for(; sit!=send; ++sit)
    {
      Star_handle star_i = get_star(sit);
      star->insert_to_star(star_i->center_point(), star_i->index_in_star_set(),
                           false/*conditional*/);
      // "false": no condition because they all need to be considered to get the
      // proper star. Not being in conflict with the star when your number comes up
      // does _NOT_ mean you're not part of the first ring of the new star in the end...
    }
  }

  void create_star(const Point_2 &p,
                   int pid,
                   Star_handle& star,
                   const bool reset_star = true) const
  {
#ifdef USE_ANISO_TIMERS
    std::clock_t start_time = clock();
#endif
    if(reset_star)
    {
      Metric m_p = metric_field()->compute_metric(p);
      star->reset(p, pid, m_p);
#ifdef ANISO_DEBUG
      std::cout << "compute metric" << std::endl;
      std::cout << "m_p:"  << std::endl << m_p.get_transformation() << std::endl;
      std::cout << "check:"  << std::endl << star->metric().get_transformation() << std::endl;
#endif
    }

#ifdef ANISO_USE_BOUNDING_BOX_VERTICES_AS_POLES
    //adding the 4 bounding box vertices (they are the first 4 stars)
    if(this->m_starset.size() > 4)
    {
      for(int i=0; i<4; ++i)
      {
        Star_handle star_i = get_star(i);
        star->insert_to_star(star_i->center_point(),
                             star_i->index_in_star_set(),
                             false /*no condition*/);
      }
    }
#endif

    CGAL_precondition((pid == 0 || !m_stars_czones.empty()) &&
                      "empty conflict map at the creation of the new star");

    typename Stars_conflict_zones::iterator czit = m_stars_czones.begin();
    typename Stars_conflict_zones::iterator czend = m_stars_czones.end();
    for(; czit!=czend; ++czit)
    {
      Index i = czit->first;
      Star_handle star_i = get_star(i);
      star->insert_to_star(star_i->center_point(), star_i->index_in_star_set(), false);
        //"false": no condition because they should be there for consistency
    }

    insert_from_kd_tree(star);
    star->clean();
  }

  Star_handle create_star(const Point_2 &p,
                          int pid)
  {
    Star_handle star = new Star(m_criteria, m_pdomain);
    create_star(p, pid, star);
    return star;
  }

  //this version of "perform_insertions()" uses the conflict zones previous computed
  Index perform_insertions(const Point_2& p,
                           const Index& this_id)
  {
    typename Stars_conflict_zones::iterator czit = m_stars_czones.begin();
    typename Stars_conflict_zones::iterator czend = m_stars_czones.end();
    for(; czit!=czend; ++czit)
    {
      Index i = czit->first;
      Star_handle si = get_star(i);
      Vertex_handle vi;

      Conflict_zone& si_czone = m_stars_czones.conflict_zone(i);

//      vi = si->insert_to_star(p, this_id, true/*conditional*/);
      vi = si->insert_to_star_hole(p, this_id,
                                   si_czone.boundary_edges(),
                                   si_czone.faces());

      if(vi == Vertex_handle()) // no conflict (should not happen)
      {
        std::cout << "No conflict in the insertion of a point in a conflict star" << std::endl;
        continue;
      }
      else if(vi->info() < this_id) // already in star set (should not happen)
      {
        std::cout << "Warning! Insertion of p"<< this_id
                  << " (" << p << ") in S" << si->index_in_star_set()
                  << " failed. vi->info() :"<< vi->info() << std::endl;

        si->print_vertices(true);
        si->print_faces();
        return vi->info();
      }
      else
      {
        // the insertion was done correctly... But we need to check for OLDER points
        // that could become in conflict with the star after the insertion of
        // the new point...
        insert_from_kd_tree(si);
      }
    }
    return this_id;
  }

  // this version is called when conditional is false : the point is inserted
  // in target_stars regardless of conflicts
  template<typename Stars>
  Index perform_insertions(const Point_2& p,
                           const Index& this_id,
                           const Stars& target_stars)
  {
    typename Stars::const_iterator it = target_stars.begin();
    typename Stars::const_iterator itend = target_stars.end();
    for(; it!=itend; it++)
    {
      Star_handle si = *it;
      Index i = si->index_in_star_set();
      Vertex_handle vi = si->insert_to_star(p, this_id, false/*conditional*/);
        // equivalent to directly calling si->base::insert(tp)

      if(vi == Vertex_handle()) // no conflict
        continue;
      else if(vi->info() < this_id) // already in star set (should not happen)
      {
        std::cout << "Warning! Insertion of p"<< this_id
                  << " (" << p << ") in S" << si->index_in_star_set()
                  << " failed. vi->info() :"<< vi->info() << std::endl;
        remove_from_stars(this_id, target_stars.begin(), ++it);

        si->print_vertices(true);
        si->print_faces();
        return vi->info();
      }
      else // inserted, standard configuration
      {
        // Conflict zones are not computed for the target_stars, so we need
        // to create entries in the conflict zones map (for fill_ref_queue)
        m_stars_czones.conflict_zone(i);

        // the insertion was done correctly... But we need to check for OLDER points
        // that could become in conflict with the star after the insertion of
        // the new point...
        insert_from_kd_tree(si);
      }
    }
    return this_id;
  }

  Index insert_to_stars(const Point_2& p,
                        const bool conditional)
  {
    Index this_id = static_cast<Index>(number_of_stars());

    Index id;
    if(conditional)
      id = perform_insertions(p, this_id); // stars in conflict
    else
      id = perform_insertions(p, this_id, m_starset); // insert in all stars

    return id;
  }

  Index insert(const Point_2& p,
               const bool conditional)
  {
    if(conditional &&
       m_stars_czones.status() != Stars_conflict_zones::CONFLICT_ZONES_ARE_KNOWN)
      std::cout << "Conflict zones unknown at insertion time...insert()" << std::endl;

    std::cout << "insert with p : "<< p << std::endl;

    Index id = insert_to_stars(p, conditional);
    if(id < 0 || id < (int) number_of_stars())
      return id;

    Star_handle star;
    star = create_star(p, id);

    if(star->index_in_star_set() != static_cast<Index>(number_of_stars()))
      std::cout << "WARNING in insert..." << std::endl;

    if(conditional)
    {
      // need to have the id of the new star in the list of ids to check in fill_ref_queue
      m_stars_czones.conflict_zone(number_of_stars());
    }
    else
    {
      std::cout << "clearing czones" << std::endl;
      // if(!conditional), the queue is filled from all stars anyway so the czones
      // have no more utility and should be cleared ASAP.
      m_stars_czones.clear();
    }

    m_starset.push_back(star);
    m_kd_tree.insert(star->index_in_star_set());
#ifndef NO_USE_AABB_TREE_OF_BBOXES
    m_aabb_tree.insert(AABB_primitive(star));
#endif
    return id;
  }

private:
  //Adding the 4 points of the domain's bbox to avoid infinite cells
  void initialize_bounding_box_vertices()
  {
    std::cout << "ini bounding box vertices" << std::endl;

    std::vector<Point_2> bbox_vertices;
    Bbox_2 bbox = this->m_pdomain->get_bbox();

    FT xmin = bbox.xmin(), xmax = bbox.xmax();
    FT ymin = bbox.ymin(), ymax = bbox.ymax();

    bbox_vertices.push_back(Point_2(xmin, ymin));
    bbox_vertices.push_back(Point_2(xmax, ymin));
    bbox_vertices.push_back(Point_2(xmax, ymax));
    bbox_vertices.push_back(Point_2(xmin, ymax));

    typename std::vector<Point_2>::iterator it = bbox_vertices.begin();
    for(; it!=bbox_vertices.end(); ++it)
      insert(*it, false /*no condition*/);
    std::cout << "done ini bounding vertices: " << bbox_vertices.size() << std::endl;
  }

protected:
  void initialize_stars(const int nb = 50)
  {
#ifdef ANISO_VERBOSE
    std::cout << "Initializse "<< nb << " stars..." << std::endl;
#endif

#ifdef ANISO_USE_BOUNDING_BOX_VERTICES_AS_POLES
    initialize_bounding_box_vertices();
#endif

    //The initial points need to be picked more cleverly as they completely ignore
    //the input metric field right now TODO
    typename Domain::Pointset initial_points = this->m_pdomain->get_boundary_points(2*nb);

#ifdef ANISO_VERBOSE
    std::cout << "Picked " << initial_points.size() << " initial points" << std::endl;
    for(typename Domain::Pointset::iterator it=initial_points.begin(); it!=initial_points.end(); ++it)
      std::cout << it->x() << " " << it->y() << " " << std::endl;
#endif

    typename Domain::Pointset::iterator pi = initial_points.begin();
    typename Domain::Pointset::iterator pend = initial_points.end();
    int nbdone = 0;
    for (; pi != pend && nbdone < nb; pi++)
    {
      if(nbdone > 0 && nbdone % 100 == 0)
        clean_stars();

      std::size_t this_id = this->m_starset.size();
      int id = -1;

      //if(m_refinement_condition(*pi)) TODO
      id = insert(*pi, false/*under no condition*/);

      if(static_cast<Index>(this_id) == id)
        nbdone++;
    }
    clean_stars();

#ifdef ANISO_VERBOSE
    std::cout << "done (" << this->m_starset.size() << " stars)." << std::endl;
#endif
  }

public:
  //Resuming
  void resume_from_mesh_file_(const char* filename)
  {
    if(!m_starset.empty())
    {
      std::cout << "resuming with a non empty star set... ?" << std::endl;
      return;
    }

    std::ifstream in(filename);
    std::string word;
    int useless, nv;
    FT x,y;

    in >> word >> useless; //MeshVersionFormatted i
    in >> word >> useless; //Dimension d
    in >> word >> nv;

    std::cout << "filename: " << filename << std::endl;
    std::cout << nv << " vertices in file" << std::endl;

    for(int i=0; i<nv; ++i)
    {
      in >> x >> y >> useless;
      Point_2 p(x,y);
      if(m_starset.size() < 10)
      {
        insert(p, false /*no condition*/);
      }
      else
      {
        if(m_starset.size() == 10)
          build_aabb_tree();

        compute_conflict_zones(p);
        m_stars_czones.compute_elements_needing_check();
        insert(p, true/*conditional*/);
      }
    }
    clean_stars();

    std::cout << "consistency of the triangulation: " << m_starset.is_consistent(true);
    std::cout << " " << m_starset.is_consistent(true) << std::endl;

    std::ofstream out("resumed.mesh");
    output_medit(m_starset, out);
    std::ofstream outoff("resumed.off");
    output_off(m_starset, outoff);
  }

protected:
  void switch_to_volume_bboxes()
  {
    Star_iterator sit = m_starset.begin();
    Star_iterator sitend = m_starset.end();
    for(; sit!=sitend; ++sit)
    {
      Star_handle star_i = get_star(sit);
      if(!star_i->is_in_2D_mesh())
      {
        star_i->is_in_2D_mesh() = true;
        star_i->invalidate_bbox_cache();
        star_i->update_bbox();
      }
    }
  }

//Constructors
public:
  Anisotropic_refine_trunk(Starset& starset_,
                           const Domain* pConstrain_,
                           const Criteria* criteria_,
                           const Metric_field* metric_field_,
                           AABB_tree& aabb_tree_,
                           Kd_tree& kd_tree_,
                           Stars_conflict_zones& stars_czones_)
  :
    m_starset(starset_),
    m_pdomain(pConstrain_),
    m_criteria(criteria_),
    m_metric_field(metric_field_),
    m_aabb_tree(aabb_tree_),
    m_kd_tree(kd_tree_),
    m_stars_czones(stars_czones_)
  { }

};  // Anisotropic_refine_trunk

}  // Anisotropic_mesh_2
}  // CGAL

#endif // CGAL_ANISOTROPIC_MESH_2_REFINE_trunk_H

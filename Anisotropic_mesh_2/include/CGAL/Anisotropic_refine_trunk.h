#ifndef CGAL_ANISOTROPIC_MESH_2_REFINE_TRUNK_H
#define CGAL_ANISOTROPIC_MESH_2_REFINE_TRUNK_H

//#define ANISO_DEBUG_CZONES
#define ANISO_USE_BOUNDING_BOX_VERTICES_AS_POLES

#include <CGAL/Domain_2.h>
#include <CGAL/Criteria.h>
#include <CGAL/Metric_field.h>
#include <CGAL/Stretched_Delaunay_2.h>
#include <CGAL/Starset.h>

#include <CGAL/IO/Star_set_output.h>

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
  typedef typename std::vector<Edge>             Edge_vector;

private:
  //Conflict zones defining vectors. Valid only until the insertion of the point!
  Edge_vector m_boundary_edges;
  Face_handle_vector m_internal_faces;

public:
  Edge_vector& boundary_edges() { return m_boundary_edges; }
  const Edge_vector& boundary_edges() const { return m_boundary_edges; }
  Face_handle_vector& internal_faces() { return m_internal_faces; }
  const Face_handle_vector& internal_faces() const { return m_internal_faces; }

  bool empty()
  {
    return m_internal_faces.empty();
  }

  //Constructors
  Conflict_zone() : m_boundary_edges(), m_internal_faces() { }
};

template<typename K, typename KExact = K>
class Stars_conflict_zones
{
  typedef Stars_conflict_zones<K, KExact>                  Self;

public:
  typedef Stretched_Delaunay_2<K, KExact>                  Star;
  typedef CGAL::Anisotropic_mesh_2::Starset<K, KExact>     Starset;
  typedef Star*                                            Star_handle;
  typedef typename std::vector<Star_handle>                Star_vector;
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

  bool is_empty() const { return m_conflict_zones.empty(); }
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
  {

  }

  void check_from_internal_faces(iterator mit,
                                  const std::map<Facet_ijk, int>& internal_faces_counter)
  {
  }

  void check_from_boundary_edges(iterator mit)
  {
  }

  void compute_elements_needing_check()
  {
//todo -----------------------
    //WHEN REFINEMENT_CONDITIONS ARE ADDED AGAIN THERE NEEDS TO BE A CHECK HERE
    //NOT TO CONSIDER USELESS faceS/CELLS
//-----------------------

    //What is being done here:
    //we need to make sure that no inconsistency appears in stars that are not
    //in conflict with the point. This might appear when a simplex is in conflict
    //in one star, but not in the others.

    if(m_are_checks_computed)
      return;

    //Count the elements. In the general case, an internal face(cell) is three(four)
    //times in conflict, and there is thus nothing to check.
    std::map<Facet_ijk, int> internal_faces_counter;

    iterator mit = m_conflict_zones.begin();
    iterator mend = m_conflict_zones.end();
    for(; mit!=mend; ++mit)
    {
      Czone& czi = mit->second;
      if(czi.empty())
        continue;

      //internal faces
      typename Star::Face_handle_handle fit = czi.internal_faces().begin();
      typename Star::Face_handle_handle fend = czi.internal_faces().end();
      for(; fit!=fend; ++fit)
      {
        Face_handle fh = *fit;
        if(fh->vertex(0)->info() == mit->first || //incident faces only
           fh->vertex(1)->info() == mit->first ||
           fh->vertex(2)->info() == mit->first)
          add_to_map(Facet_ijk(fh), internal_faces_counter);
      }
    }

    mit = m_conflict_zones.begin();
    for(; mit!=mend; ++mit)
    {
      check_from_internal_faces(mit, internal_faces_counter);
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
    os << czi.internal_faces().size() << " ";
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

    Face_handle_handle fit = czi.internal_faces().begin();
    Face_handle_handle fend = czi.internal_faces().end();
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
  //typedef typename CGAL::Exact_predicates_exact_constructions_kernel KExact;
  typedef K                                                        KExact;

  typedef Stretched_Delaunay_2<K>                                  Star;
  typedef typename Star::FT                                        FT;
  typedef Star*                                                    Star_handle;
  typedef typename Star::Base                                      DT; // DT_2 with vertex_base_with_info
  typedef std::vector<Star_handle>                                 Star_vector;
  typedef typename Star_vector::iterator                           Star_iterator;
  typedef std::set<Star_handle>                                    Star_set;
  typedef typename Star::Index                                     Index;
  typedef std::set<Index>                                          Index_set;
  typedef typename Star::Point_2                                   Point_2;
  typedef typename Star::TPoint_2                                  TPoint_2;
  typedef std::set<Point_2>                                        Point_set;
  typedef typename Star::Vertex_handle                             Vertex_handle;
  typedef typename Star::Face                                      Face;
  typedef typename Star::Face_handle                               Face_handle;
  typedef typename Star::Vector_2                                  Vector_2;
  typedef typename Star::Domain                         Domain;
  typedef typename Star::Criteria                                  Criteria;

  typedef CGAL::Anisotropic_mesh_2::Starset<K>                     Starset;

  typedef typename CGAL::Anisotropic_mesh_2::Metric_field<K>       Metric_field;
  typedef typename Metric_field::Metric                            Metric;

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
  void update_bboxes() const
  {
    std::size_t i;
    std::size_t N = number_of_stars();
    for(i = 0; i < N; i++)
    {
      Star_handle si = get_star(i);
      si->update_bbox();
    }
  }

  template<typename OutputIterator>
  void finite_stars_in_conflict(const Point_2& p,
                                OutputIterator oit) const
  {
    update_bboxes();

    //exact set
    for(std::size_t i = 0; i < number_of_stars(); i++)
    {
      Star_handle s = get_star(i);
      Face_handle fh;
      if(s->is_conflicted(transform_to_star_point(p, s), fh))
        *oit++ = s;
    }
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
  }

  void clean_stars()
  {
    std::cout << "clean call, size: " << number_of_stars() << std::endl;
    for(std::size_t i = 0; i < number_of_stars(); i++)
      get_star(i)->clean();
  }

  bool check_consistency_and_sliverity(Star_handle to_be_refined,
                                       const Star_handle& new_star,
                                       FT sq_radius_bound) const
  {

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

//Conflict zones
  Index simulate_insert_to_stars(const Point_2& p) const
  {
    CGAL::Timer t;
    t.start();

    Index this_id = static_cast<Index>(number_of_stars());

    // find conflicted stars
    Star_set stars;
    finite_stars_in_conflict(p, std::inserter(stars, stars.end()));

    typename Star_set::iterator it = stars.begin();
    typename Star_set::iterator iend = stars.end();
    for(; it != iend; it++)
    {
      Star_handle si = get_star(it);
      int id = si->simulate_insert_to_star(p, this_id);

      if(id == -1) // no conflict
        continue;
      else if(id < (int)this_id) // already in star set
      {
        std::cout << "Warning in simulate_insert_to_stars" << std::endl;
        return id;
      }
      else // to be inserted, standard configuration
      {
        //add them to the map of stars in conflict (conflict zones are not computed yet)
        m_stars_czones.conflict_zone(si->index_in_star_set());
      }
    }

    //std::cerr << s_in_conflict << std::endl;
    t.stop();
    //std::cout << "simulate: " << t.time() << std::endl;

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

      si->find_conflicts(p, std::back_inserter(ci.internal_faces()),
                            std::back_inserter(ci.boundary_edges()));
    }

    return m_stars_czones.conflicting_point_id();
  }

  void clear_conflict_zones() const
  {
    m_stars_czones.clear();
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
    }

#ifdef ANISO_SCHWARZENEGGER_CREATE_STAR
    Star_iterator csit = m_starset.begin();
    Star_iterator csend = m_starset.end();
    for(; csit!=csend; ++csit)
    {
      Star_handle si = get_star(csit);
      star->insert_to_star(si->center_point(), si->index_in_star_set(), false);
      si->insert_to_star(star->center_point(), star->index_in_star_set(), false);
    }
    return;
#endif

#ifdef ANISO_DEBUG_CREATE_STAR
    std::cout << "CONAN INSERT POINT YARR" << std::endl;

    Star_handle star_check = new Star(m_criteria, m_pdomain,
                                      surface_star, m_is_2D_level);

    Metric m_p = metric_field()->compute_metric(p);
    star_check->reset(p, pid, m_p, surface_star);

    //full rambo
    Star_iterator csit = m_starset.begin();
    Star_iterator csend = m_starset.end();
    for(; csit!=csend; ++csit)
    {
      Star_handle si = get_star(csit);
      star_check->insert_to_star(si->center_point(), si->index_in_star_set(), false);
    }

    int conan_rface_size = star_check->count_restricted_faces();
    star_check->clean(true);
    star_check->print_vertices();
    star_check->print_restricted_faces();
    star_check->clear();
    delete star_check;
#endif

    if(m_stars_czones.is_empty())
      std::cout << "Warning: empty conflict map at the creation of the new star" << std::endl;

#ifdef ANISO_USE_BOUNDING_BOX_VERTICES_AS_POLES
    //adding the 4 bounding box vertices (they are the first 8 stars)
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

    typename Stars_conflict_zones::iterator czit = m_stars_czones.begin();
    typename Stars_conflict_zones::iterator czend = m_stars_czones.end();
    for(; czit!=czend; ++czit)
    {
      Index i = czit->first;
      Star_handle star_i = get_star(i);
      star->insert_to_star(star_i->center_point(), star_i->index_in_star_set(), false);
        //"false": no condition because they should be there for consistency
    }
    star->clean();

#ifdef ANISO_DEBUG_CREATE_STAR
    star->print_vertices();
    star->print_restricted_faces();
    std::cout << "valid check @ create_star @ " << star->index_in_star_set() << " val: " << star->is_valid() << std::endl;

    if(conan_rface_size != normal_rface_size)
    {
      std::cout << "conan wins: " << conan_rface_size << " " << normal_rface_size << std::endl;
      assert(1==2);
    }
#endif
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
    std::cout << "no conflict zone perform insertion" << std::endl;
    return -1;
  }

  //this version is called when conditional is false (insert in target_stars regardless of conflicts)
  template<typename Stars>
  Index perform_insertions(const Point_2& p,
                           const Index& this_id,
                           const Stars& target_stars)
  {
    typename Stars::const_iterator it = target_stars.begin();
    typename Stars::const_iterator itend = target_stars.end();
    for(; it != itend; it++)
    {
      Star_handle si = *it;
      Index i = si->index_in_star_set();
      Vertex_handle vi = si->insert_to_star(p, this_id, false);
        // equivalent to directly calling si->base::insert(tp)

      if(vi == Vertex_handle())  //no conflict
        continue;
      else if(vi->info() < this_id) // already in star set (should not happen)
      {
        std::cout << "Warning! Insertion of p"<< this_id
                  << " (" << p << ") in S" << si->index_in_star_set()
                  << " failed. vi->info() :"<< vi->info() << std::endl;
        remove_from_stars(this_id, target_stars.begin(), ++it);

        si->print_vertices(true);
        si->print_faces();
        std::cout << "Metric : \n" << si->metric().get_transformation() << std::endl;
        return vi->info();
      }
      else // inserted, standard configuration
      {
        //Conflict zones are not computed for the target_stars, so we need to create
        //entries in the conflict zones map (for fill_ref_queue)
        m_stars_czones.conflict_zone(i);
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
      id = perform_insertions(p, this_id); //stars in conflict
    else
      id = perform_insertions(p, this_id, m_starset); // insert in all stars

    return id;
  }

  Index insert(const Point_2& p,
               const bool conditional)
  {
    if(conditional && m_stars_czones.status() != Stars_conflict_zones::CONFLICT_ZONES_ARE_KNOWN)
      std::cout << "Conflict zones unknown at insertion time...insert()" << std::endl;

    std::cout << "insert with p : "<< p << std::endl;

    Index id = insert_to_stars(p, conditional);
    if(id < 0 || id < (int) number_of_stars())
      return id;

    Star_handle star;
    star = create_star(p, id);

    if(star->index_in_star_set() != number_of_stars())
      std::cout << "WARNING in insert..." << std::endl;

    if(conditional)
    {
      //need to have the id of the new star in the list of ids to check in fill_ref_queue
      m_stars_czones.conflict_zone(number_of_stars());
    }
    else
    {
      //if(!conditional), the queue is filled from all stars anyway so the czones
      //have no more utility and should be cleared ASAP.
      m_stars_czones.clear();
    }

    m_starset.push_back(star);
    return id;
  }

private:
  //Adding the 8 points of the domain's bbox to avoid infinite cells
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
    std::cout << "done inibv" << std::endl;
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

      if(this_id == id)
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
      insert(p, false /*no condition*/);
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
                           Stars_conflict_zones& stars_czones_)
  :
    m_starset(starset_),
    m_pdomain(pConstrain_),
    m_criteria(criteria_),
    m_metric_field(metric_field_),
    m_stars_czones(stars_czones_)
  { }

};  // Anisotropic_refine_trunk

}  // Anisotropic_mesh_2
}  // CGAL

#endif // CGAL_ANISOTROPIC_MESH_2_REFINE_trunk_H

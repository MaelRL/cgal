#ifndef CGAL_ANISOTROPIC_MESH_2_STRETCHED_DELAUNAY_2
#define CGAL_ANISOTROPIC_MESH_2_STRETCHED_DELAUNAY_2

#include <CGAL/Profile_counter.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_domain_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Random.h>

#include <CGAL/Kd_tree_for_star_set.h>
#include <CGAL/aabb_tree/aabb_tree_bbox.h>
#include <CGAL/aabb_tree/aabb_tree_bbox_primitive.h>
#include <CGAL/aabb_tree/aniso_bbox_3.h>

#include <CGAL/bbox.h>
#include <CGAL/Metric.h>
#include <CGAL/Criteria.h>
#include <CGAL/Delaunay_traits_2.h>
#include <CGAL/Domain_2.h>
#include <CGAL/Timer.h>

#include <CGAL/helpers/combinatorics_helper.h>

#include <iostream>
#include <fstream>
#include <utility>
#include <time.h>

namespace CGAL
{
namespace Anisotropic_mesh_2
{
struct Star_index
{
  int info_;

  operator int () const { return info_;}
  operator int& () { return info_;}
  Star_index& operator=(const int i)
  {
    info_ = i;
    return *this;
  }

  Star_index(){ info_ = -1; }
};

template<typename K, typename KExact = K>
class Stretched_Delaunay_2
#ifdef CGAL_USE_BASIC_ANISO
    : public CGAL::Delaunay_triangulation_2<
                    Delaunay_traits_2<K, KExact>,
                    CGAL::Triangulation_data_structure_2<
                      CGAL::Triangulation_vertex_base_with_info_2<
                        Star_index,
                        Delaunay_traits_2<K, KExact> >,
                      CGAL::Triangulation_face_base_with_domain_info_2<
                        Delaunay_traits_2<K, KExact>,
                        Domain_2<K> > > >
#else
    : public CGAL::Constrained_Delaunay_triangulation_2<
                    Delaunay_traits_2<K, KExact>,
                    CGAL::Triangulation_data_structure_2<
                      CGAL::Triangulation_vertex_base_with_info_2<
                        Star_index,
                        Delaunay_traits_2<K, KExact> >,
                      CGAL::Triangulation_face_base_with_domain_info_2<
                        Delaunay_traits_2<K, KExact>,
                        Domain_2<K> > > >
#endif
{
  typedef Stretched_Delaunay_2<K, KExact>                       Self;
public:
  typedef K                                                     Kernel;
  typedef Delaunay_traits_2<K, KExact>                          Traits;
  typedef CGAL::Triangulation_data_structure_2<
               CGAL::Triangulation_vertex_base_with_info_2<Star_index,
                                                           Traits >,
               CGAL::Triangulation_face_base_with_domain_info_2<Traits,
                                                                Domain_2<K> >
                                                >                DS;

#ifdef CGAL_USE_BASIC_ANISO
  typedef CGAL::Delaunay_triangulation_2<Traits, DS>             Base;
#else
  typedef CGAL::Constrained_Delaunay_triangulation_2<Traits, DS> Base;
#endif

  typedef int Index;

  typedef Domain_2<K>                               Domain;
  typedef Metric_base<K, KExact>                    Metric;
  typedef typename K::FT                            FT;
  typedef typename K::Point_2                       Point_2;
  typedef typename K::Point_3                       Point_3;
  typedef typename K::Point_2                       TPoint_2; //Transformed_point_2
  typedef typename K::Vector_2                      Vector_2;
  typedef typename K::Line_2                        Line_2;
  typedef typename K::Ray_2                         Ray_2;
  typedef typename K::Circle_2                      Circle_2;
  typedef typename Base::Segment                    Segment;
  typedef typename Base::Triangle                   Triangle;
  typedef typename Base::Vertex                     Vertex;
  typedef typename Base::Edge                       Edge;
  typedef typename Base::Face                       Face;
  typedef typename DS::Vertex_handle                Vertex_handle;
  typedef typename DS::Face_handle                  Face_handle;
  typedef ::CGAL::Anisotropic_mesh_2::Stretched_criteria<K, KExact> Stretched_criteria;
  typedef Criteria_base<K>                          Criteria;
  typedef CGAL::Bbox<K>                             Bbox;
  typedef CGAL::Aniso_bbox_3<K>                     A_Bbox_3;

  typedef std::vector<Point_2>                      Point_vector;
  typedef std::vector<Vertex_handle>                Vertex_handle_vector;
  typedef std::vector<Face_handle>                  Face_handle_vector;
  typedef typename Face_handle_vector::iterator     Face_handle_handle;
  typedef std::set<Face_handle>                     Face_handle_set;
  typedef typename Vertex_handle_vector::iterator   Vertex_handle_handle;

  typedef typename KExact::Point_2                  Exact_Point_2;
  typedef typename KExact::Point_2                  Exact_TPoint_2;
  typedef typename KExact::Vector_2                 Exact_Vector_2;
  typedef CGAL::Cartesian_converter<K, KExact>      To_exact;
  typedef CGAL::Cartesian_converter<KExact, K>      Back_from_exact;

private:
  static const Index index_of_infinite_vertex = -10;

  Traits* m_traits;

  Point_2 m_center_point; // before transformation by m_metric.transform
  Vertex_handle m_center; // after transformation by m_metric.transform

  To_exact to_exact;
  Back_from_exact back_from_exact;
  Metric m_metric;
  const Domain* m_pdomain;
  Stretched_criteria* m_criteria;

  mutable bool m_is_topological_disk;
  mutable bool m_is_valid_topo_disk;
  mutable Bbox m_bbox;
  mutable A_Bbox_3 m_a_bbox_3;
  mutable bool m_is_valid_bbox;
  mutable bool m_bbox_needs_aabb_update;
  mutable bool m_metric_needs_update;

  mutable bool is_cache_dirty;
  mutable Vertex_handle_vector finite_adjacent_vertices_cache;
  mutable Face_handle_vector finite_incident_faces_cache;

public:
  mutable bool m_active; // used for removing points at the end of the process

public:
  static Index infinite_vertex_index() { return index_of_infinite_vertex; }

  bool is_infinite() const
  {
    return (infinite_vertex_index() == index_in_star_set());
  }

  inline bool is_infinite_vertex(const Vertex_handle& vh) const
  {
    return Base::is_infinite(vh);
  }

  inline bool is_infinite_face(const Face_handle& fh) const
  {
    return Base::is_infinite(fh);
  }

public:
  Point_2& center_point()             { return m_center_point; }
  const Point_2& center_point() const { return m_center_point; }

  Index& index_in_star_set()             { return m_center->info(); }
  const Index& index_in_star_set() const { return m_center->info(); }

public:
  Vertex_handle center() const      { return m_center; }
  Metric metric() const             { return m_metric; }
  Metric& metric()                  { return m_metric; }
  Traits* traits() const            { return m_traits; }
  const Domain* domain() const { return m_pdomain; }
  const Stretched_criteria* criteria() const { return m_criteria; }
  void set_criteria(const Criteria* criteria_)
  {
    delete m_criteria;
    m_criteria = new Stretched_criteria(*m_traits, criteria_);
  }

public:
  const Bbox& bbox(const bool verbose = false) const
  {
    update_bbox(verbose);
    return m_bbox;
  }

// to abuse the AABB structure with a 2D star set
  const A_Bbox_3& AABB_bbox(const bool verbose = false) const
  {
    update_bbox(verbose);
    return m_a_bbox_3;
  }

  Point_3 AABB_center_point() const
  {
    FT x = center_point().x();
    FT y = center_point().y();
#ifdef LIFT_AABB_TREE
    return Point_3(x, y, x);
#else
    return Point_3(x,y,0.);
#endif
  }

  bool& bbox_needs_aabb_update(){ return m_bbox_needs_aabb_update; }
  const bool& bbox_needs_aabb_update() const { return m_bbox_needs_aabb_update; }
  bool& metric_needs_update(){ return m_metric_needs_update; }
  const bool& metric_needs_update() const { return m_metric_needs_update; }

public:
  bool has_vertex(const int i) const
  {
    Vertex_handle_handle vit = finite_adjacent_vertices_begin();
    Vertex_handle_handle vend = finite_adjacent_vertices_end();
    for(; vit != vend; vit++)
    {
      Vertex_handle vh = *vit;
      if(vh->info() == i)
        return true;
    }
    return false;
  }

  bool has_vertex(const int i, Vertex_handle& v) const
  {
    Vertex_handle_handle vit = finite_adjacent_vertices_begin();
    Vertex_handle_handle vend = finite_adjacent_vertices_end();
    for(; vit != vend; vit++)
    {
      Vertex_handle vh = *vit;
      if(vh->info() == i)
      {
        v = vh;
        return true;
      }
    }
    return false;
  }

  //this function checks if f is incident to the center
  bool is_in_star(const Face_handle& f) const
  {
    for(int i = 0; i < 3; i++)
    {
      Vertex_handle v = f->vertex(i);
      if(v->info() == m_center->info())
        return true;
    }
    return false;
  }

  bool has_face(const Face_handle& fh) const
  {
    Face_handle_handle fit = finite_incident_faces_begin();
    Face_handle_handle fend = finite_incident_faces_end();
    for(; fit!=fend; fit++)
      if(Facet_ijk(*fit) == Facet_ijk(fh))
        return true;
    return false;
  }

  bool has_face(const Facet_ijk& f, Face_handle& fh) const
  {
    boost::array<int, 3> dids;
    Face_handle_handle fit = finite_incident_faces_begin();
    Face_handle_handle fend = finite_incident_faces_end();
    for (; fit!=fend; fit++)
    {
      Face_handle fhi = *fit;
      for (int i=0; i<3; i++)
        dids[i] = fhi->vertex(i)->info();
      std::sort(dids.begin(), dids.end());
      if(std::equal(dids.begin(), dids.end(), f.vertices().begin()))
      {
        fh = fhi;
        return true;
      }
    }
    return false;
  }

  bool has_face(int *cids, Face_handle& fh) const
  {
    int dids[3];
    Face_handle_handle fit = finite_incident_faces_begin();
    Face_handle_handle fend = finite_incident_faces_end();
    for (; fit != fend; fit++)
    {
      Face_handle fhi = *fit;
      for (int i=0; i<3; i++)
        dids[i] = fhi->vertex(i)->info();
      if (is_same_ids<3>(cids, dids))
      {
        fh = fhi;
        return true;
      }
    }
    return false;
  }

  // CACHE RELATED FUNCTIONS
public:
  inline void update_star_caches() const
  {
    if(!is_cache_dirty)
      return;

    CGAL_PROFILER("[update_star_caches]");
    finite_adjacent_vertices_cache.clear();
    finite_incident_faces_cache.clear();

    // update neighboring vertices
    typename Base::Vertex_circulator vc = Base::incident_vertices(m_center), done_v(vc);
    if(vc != 0)
    {
      do
      {
        if(!is_infinite_vertex(vc))
          finite_adjacent_vertices_cache.push_back(vc);
      }
      while(++vc!=done_v);
    }

    typename Base::Face_circulator fc = Base::incident_faces(m_center), done_f(fc);
    if(fc != 0)
    {
      do
      {
        if(!is_infinite_face(fc))
          finite_incident_faces_cache.push_back(fc);
      }
      while(++fc!=done_f);
    }

    is_cache_dirty = false;
  }

  void invalidate_bbox_cache() const
  {
    m_is_valid_bbox = false;
  }

  void invalidate_cache() const
  {
    is_cache_dirty = true;
    m_is_valid_bbox = false;
    m_is_valid_topo_disk = false;
  }

  inline Vertex_handle_handle finite_adjacent_vertices_begin() const
  {
    update_star_caches();
    return finite_adjacent_vertices_cache.begin();
  }
  inline Vertex_handle_handle finite_adjacent_vertices_end() const
  {
    return finite_adjacent_vertices_cache.end();
  }

  inline Face_handle_handle finite_incident_faces_begin() const
  {
    update_star_caches();
    return finite_incident_faces_cache.begin();
  }
  inline Face_handle_handle finite_incident_faces_end() const
  {
    return finite_incident_faces_cache.end();
  }

public:
  bool is_topological_disk() const
  {
    update_is_topological_disk();
    return m_is_topological_disk;
  }

  //need to be outside of the function below since it is used as a template...
  struct Vertex_neighborhood
  {
  public:
    Vertex_handle m_fn, m_sn; //first and second (potential) neighbors
    bool visited; //used during cycle check

    Vertex_neighborhood(Vertex_handle fn_) :
      m_fn(fn_), m_sn(Vertex_handle()), visited(false)
    { }
  };

  void update_is_topological_disk() const
  {
    if(m_is_valid_topo_disk)
      return;

    CGAL_PROFILER("[Update is_topological_disk]");
    Face_handle_handle fit = finite_incident_faces_begin();
    Face_handle_handle fend = finite_incident_faces_end();
    if(fit == fend) // no faces
    {
      m_is_topological_disk = false;
      m_is_valid_topo_disk = true;
      return;
    }

    //for each vertex V of the star S, the map contains its n neighbors (except
    //the center of the star).
    typedef std::map<Vertex_handle, Vertex_neighborhood> Cycle_map;
    Cycle_map cycle_map;

    struct Cycle_map_insertor
    {
      Cycle_map& m_cmap;

      bool operator()(Vertex_handle v1, Vertex_handle v2)
      {
        //insert v2 as neighbor of v1
        std::pair<typename Cycle_map::iterator, bool> is_insert_successful;
        is_insert_successful = m_cmap.insert(
                                 std::make_pair(v1, Vertex_neighborhood(v2)));
        if(!is_insert_successful.second) // v1 already is in the cycle map
        {
          Vertex_neighborhood& vhn = is_insert_successful.first->second;
          if(vhn.m_sn != Vertex_handle())
          {
            // edge [center,v1] has more than 2 incident facets
            return false;
          }
          vhn.m_sn = v2; // setting v1's second neighbor
        }
        return true;
      }

      Cycle_map_insertor(Cycle_map& cmap_) : m_cmap(cmap_) { }
    };

    Cycle_map_insertor cm_insertor(cycle_map);

    //fill the map from the facets
    for(; fit != fend; fit++)
    {
      std::vector<Vertex_handle> vns; // the 2 vertices of f that are not the center of the star
      for(int i = 0; i < 3; i++)
      {
        Vertex_handle v = (*fit)->vertex(i);
        if(v != m_center)
          vns.push_back(v);
      }

      //insert front() and back() into each other's neighborhoods
      if( !cm_insertor(vns.front(), vns.back()) ||
          !cm_insertor(vns.back(), vns.front()) )
      {
        //front() or back() already had 2 neighbors
        m_is_topological_disk = false;
        m_is_valid_topo_disk = true;
        return;
      }
    }

#ifdef ANISO_DEBUG_TOPO_DISK
    std::cout << "printing cycle map. Size : " << cycle_map.size() << std::endl;
    for(typename Cycle_map::iterator it=cycle_map.begin(); it!=cycle_map.end(); ++it)
    {
      std::cout << it->first->info() << " || ";
      std::cout << it->second.m_fn->info() << " , ";
      if(it->second.m_sn != Vertex_handle())
        std::cout << it->second.m_sn->info();
      std::cout << std::endl;
    }
#endif

    //check that all the vertices have 2 neighbors
    typename Cycle_map::iterator cmit = cycle_map.begin();
    for(; cmit!=cycle_map.end(); ++cmit)
    {
      //if in the map, m_fn is set properly so it is sufficient to check m_sn
      if(cmit->second.m_sn == Vertex_handle())
      {
        m_is_topological_disk = false;
        m_is_valid_topo_disk = true;
        return;
      }
    }

    //We then check that the vertices form exactly one cycle (we could have n=2
    //for all the vertices, but more than one cycle).
    typename Cycle_map::iterator it = cycle_map.begin();
    Vertex_handle starting_vh = it->first;
    Vertex_handle current_vh = starting_vh;
    Vertex_handle next_vh = it->second.m_fn;
    it->second.visited = true;

    while(true)
    {
      current_vh = next_vh;
      Vertex_neighborhood& vhn = cycle_map.find(current_vh)->second;

      if(current_vh == starting_vh)
        break;

      if(vhn.visited)
      {
        // next vertex was already visited but is not the starting vertex
        m_is_topological_disk = false;
        m_is_valid_topo_disk = true;
        return;
      }
      vhn.visited = true;

      if((vhn.m_fn == starting_vh || vhn.m_sn == starting_vh) &&
         cycle_map.find(vhn.m_fn)->second.visited &&
         cycle_map.find(vhn.m_sn)->second.visited )
        break; // next_vh will be starting_vh, which means we have a cycle

      if(!(cycle_map.find(vhn.m_fn)->second.visited))
        next_vh = vhn.m_fn;
      else
        next_vh = vhn.m_sn; // might already have been visited, we'll check that later
    }

    //we're back to the starting vertex, let's check that everything was visited
    it = cycle_map.begin();
    for(; it!=cycle_map.end(); ++it)
    {
      if(!it->second.visited)
      {
        //didn't visit everything
        m_is_topological_disk = false;
        m_is_valid_topo_disk = true;
        return;
      }
    }

    //All good, it's a topo disk!
    m_is_valid_topo_disk = true;
    m_is_topological_disk = true;
  }

  void build_a_bbox_from_bbox() const
  {
    // the bbox used in the AABB tree.
    // can try to lift the bboxes to see if it's better for the AABB tree...

    FT xmin = m_bbox.xmin(), xmax = m_bbox.xmax(),
       ymin = m_bbox.ymin(), ymax = m_bbox.ymax();

#ifdef LIFT_AABB_TREE
    m_a_bbox_3 = A_Bbox_3(xmin, ymin, xmin,
                          xmax, ymax, xmax); // bboxes do not live in the same plane
#else
    m_a_bbox_3 = A_Bbox_3(xmin, ymin, -1., xmax, ymax, 1.); // all the bboxes have the same z coordinates
#endif

    return;
  }

  void update_bbox(const bool verbose = false) const
  {
    if(m_is_valid_bbox)
      return;

    Bbox_2 old_bbox = m_bbox;

    if(verbose)
      std::cout << "Update Bbox...";
    CGAL_PROFILER("[update_bbox]");

    m_bbox = this->compute_bbox();

    m_is_valid_bbox = true;

#ifndef NO_USE_AABB_TREE_OF_BBOXES
   //checking if the bbox has grown
    if(m_bbox.xmin() < old_bbox.xmin() || m_bbox.xmax() > old_bbox.xmax() ||
       m_bbox.ymin() < old_bbox.ymin() || m_bbox.ymax() > old_bbox.ymax())
    {
#ifdef DEBUG_UPDATE_AABB_TREE
      std::cout << "growing bbox in star : " << index_in_star_set() << std::endl;
      std::cout << m_bbox.xmin() << " " << m_bbox.xmax() << " " << m_bbox.ymin();
      std::cout << " " << m_bbox.ymax() << std::endl;
#endif
      m_bbox_needs_aabb_update = true;
   }
#endif

    build_a_bbox_from_bbox();

    if(verbose)
      std::cout << "done." << std::endl;
  }

  Bbox compute_bbox() const // compute bbox of incident surface Delaunay balls
  {
    typename Traits::Compute_squared_distance_2 csd
        = m_traits->compute_squared_distance_2_object();

    Bbox bb = m_center->point().bbox();
    Face_handle_handle fit = finite_incident_faces_begin();
    Face_handle_handle fend = finite_incident_faces_end();
    for(; fit != fend; fit++)
    {
      Face_handle fh = *fit;
      TPoint_2 tp = compute_circumcenter(fh);
      FT squared_radius = csd(tp, m_center->point());
      Circle_2 c(tp, squared_radius);
      bb = bb + c.bbox();
    }
    return m_metric.inverse_transform(bb);
  }

  bool is_inside_bbox(const Point_2& p) const
  {
    update_bbox();
    return p.x() >= m_bbox.xmin() && p.x() <= m_bbox.xmax()
        && p.y() >= m_bbox.ymin() && p.y() <= m_bbox.ymax();
  }

public:
  bool is_face_visited(const Face_handle& fh) const
  {
    return fh->is_face_visited();
  }

public:
  inline FT compute_volume(const Face_handle& fh) const
  {
    return m_criteria->compute_volume(fh->vertex(0)->point(),
                                      fh->vertex(1)->point(),
                                      fh->vertex(2)->point());
  }

  inline FT compute_radius_edge_ratio_overflow(const Face_handle& fh) const
  {
    return m_criteria->radius_edge_ratio_overflow(fh->vertex(0)->point(),
                                                  fh->vertex(1)->point(),
                                                  fh->vertex(2)->point());
  }

  inline FT compute_squared_radius_edge_ratio(const Face_handle& fh) const
  {
    return m_criteria->compute_squared_radius_edge_ratio(fh->vertex(0)->point(),
                                                         fh->vertex(1)->point(),
                                                         fh->vertex(2)->point());
  }

  inline FT compute_squared_radius_edge_ratio(const TPoint_2& tp1,
                                              const TPoint_2& tp2,
                                              const TPoint_2& tp3) const
  {
    return m_criteria->compute_squared_radius_edge_ratio(tp1, tp2, tp3);
  }

  inline FT compute_circumradius_overflow(const Face_handle& fh) const
  {
    TPoint_2 tp = compute_circumcenter(fh);
    TPoint_2 tp2 = fh->vertex(0)->point();
    FT sqr = m_traits->compute_squared_distance_2_object()(tp, tp2);
    return sqr - m_criteria->criteria->squared_face_circumradius;
  }

  inline FT compute_squared_circumradius(const Face_handle& fh) const
  {
    return m_criteria->compute_squared_circumradius(fh->vertex(0)->point(),
                                                    fh->vertex(1)->point(),
                                                    fh->vertex(2)->point());
  }

  inline FT compute_squared_circumradius(const TPoint_2& p1,
                                         const TPoint_2& p2,
                                         const TPoint_2& p3) const
  {
    return m_criteria->compute_squared_circumradius(p1, p2, p3);
  }

  inline FT compute_element_quality(const TPoint_2& p1,
                                    const TPoint_2& p2,
                                    const TPoint_2& p3) const
  {
    return m_criteria->element_quality(p1, p2, p3);
  }

  inline bool is_inside(const Face_handle& fh) const
  {
    if(Base::is_infinite(fh))
      return false;
    return fh->template is_inside<Self>(*this, *domain(), *traits());
  }

#ifdef ANISO_USE_INSIDE_EXACT
  inline bool is_inside_exact(const Face_handle& fh) const
  {
    if(this->is_infinite(fh))
      return false;
    return fh->template is_inside_exact<Self>(*this, *domain(), *traits());
  }
#endif

public:
  bool is_in_delaunay_ball(const TPoint_2& tp,
                           const Face_handle& fh) const
  {
    CGAL_PROFILER("[conflict test]");
    return (Base::side_of_oriented_circle(fh, tp) == CGAL::ON_POSITIVE_SIDE);
  }

  inline bool is_in_a_delaunay_ball(const TPoint_2& tp,
                                    Face_handle& in_which_face) const
  {
    CGAL_PROFILER("[is_in_a_volume_delaunay_ball]");
    Face_handle_handle fit = finite_incident_faces_begin();
    Face_handle_handle fend = finite_incident_faces_end();
    for(; fit!=fend; fit++)
    {
      Face_handle fh = *fit;
      if(is_in_delaunay_ball(tp, fh))
      {
        in_which_face = fh;
        return true;
      }
    }
    return false;
  }

  inline bool is_conflicted(const TPoint_2& tp,
                            Face_handle& in_which_face) const
  {
    if(!is_inside_bbox(m_metric.inverse_transform(tp)))
      return false;
    else
      return is_in_a_delaunay_ball(tp, in_which_face);
  }

  template<typename FacesOutputIterator,
           typename BoundaryEdgesOutputIterator>
  bool find_conflicts(const Point_2& p,
                      FacesOutputIterator foit,
                      BoundaryEdgesOutputIterator beoit)
  {
    int dim = Base::dimension();
    if(dim <= 1)
      return false;

    Face_handle fh;
    TPoint_2 tp = m_metric.transform(p);
    if(!is_conflicted(tp, fh))
    {
      std::cout << "no conflict in find_conflict.........(bad)" << std::endl;
      return false;
    }
    Base::get_conflicts_and_boundary(tp, foit, beoit, fh);

    return true;
  }

  template<typename FacesOutputIterator>
  bool find_conflicts(const Point_2& p, FacesOutputIterator foit) const
  {
    return find_conflicts(p, foit, Emptyset_iterator());
  }

  inline void get_inverse_transformed_points(Point_2 *ps,
                                             const Face_handle& fh) const
  {
    int k = 0;
    for (int i = 0; i < 3; i++)
      ps[k++] = m_metric.inverse_transform(fh->vertex(i)->point());
  }

  inline void get_points(TPoint_2 *ps, const Face_handle &fh) const
  {
    int k = 0;
    for (int i = 0; i < 3; i++)
      ps[k++] = fh->vertex(i)->point();
  }

  TPoint_2 compute_circumcenter(const TPoint_2& p0,
                                const TPoint_2& p1,
                                const TPoint_2& p2) const
  {
#ifdef ANISO_USE_CC_EXACT
    return back_from_exact(m_traits->construct_circumcenter_2_object()(
                             to_exact(p0), to_exact(p1), to_exact(p2)));
#else
    return m_traits->construct_circumcenter_2_object()(p0, p1, p2);
#endif
  }

  TPoint_2 compute_circumcenter(const Face_handle& fh) const
  {
#ifdef ANISO_USE_CC_EXACT
    return back_from_exact(compute_exact_circumcenter(fh));
#else
    return compute_exact_circumcenter(fh);
#endif
  }

  Exact_TPoint_2 compute_exact_circumcenter(const Face_handle& fh) const
  {
    TPoint_2 ps[3];
    get_points(ps, fh);
#ifdef ANISO_USE_CC_EXACT
    return m_traits->construct_circumcenter_2_object()(
          to_exact(ps[0]), to_exact(ps[1]), to_exact(ps[2]));
#else
    return m_traits->construct_circumcenter_2_object()(ps[0], ps[1], ps[2]);
#endif
  }

public:
  std::size_t clean(bool verbose = false) //remove non-adjacent vertices
  {
    typedef std::pair<TPoint_2, int>          PPoint;
    std::vector<PPoint> backup;
    std::size_t nbv = this->number_of_vertices();

    Vertex_handle_handle vit = finite_adjacent_vertices_begin();
    Vertex_handle_handle vend = finite_adjacent_vertices_end();
    if(vend-vit == nbv-1) // -1 due to center
      return 0; // all the vertices are already adjacent

    /*
    if(vend-vit > 0.8*(nbv-1))
      return 0;
*/

    // backup star vertices
    backup.push_back(PPoint(m_center->point(), index_in_star_set()));
    for(; vit != vend; vit++)
      backup.push_back(PPoint((*vit)->point(), (*vit)->info()));

    this->clear();
    this->infinite_vertex()->info() = index_of_infinite_vertex;

    // re-insert
    m_center = insert(backup[0].first);
    m_center->info() = backup[0].second;
    //other vertices
    for(std::size_t i = 1; i < backup.size(); i++)
      Base::insert(backup[i].first)->info() = backup[i].second;

    invalidate_cache();

    if(verbose)
    {
      std::cout << "clean (old/new) @ : " << index_in_star_set() << " :: ";
      std::cout << nbv << " " << this->number_of_vertices() << std::endl;
    }
    return (nbv - backup.size());
  }

  void remove(Vertex_handle v)
  {
    Base::remove(v);
    invalidate_cache();
  }

  void remove(const Index& i)
  {
    typename Base::Finite_vertices_iterator vit = this->finite_vertices_begin();
    typename Base::Finite_vertices_iterator vend = this->finite_vertices_end();
    for(; vit != vend; vit++)
    {
      if(vit->info() == i)
      {
        this->remove(vit);
        break; //warning : vit is broken
      }
    }
  }

  // warning : point should be transformed before this insert function is called
  Vertex_handle insert(const TPoint_2& tp, const Face_handle& fh)
  {
    Vertex_handle retval = Base::insert(tp, fh);
    invalidate_cache();
    return retval;
  }

  // warning : point should be transformed before this insert function is called
  Vertex_handle insert(const TPoint_2& tp)
  {
    Vertex_handle retval = Base::insert(tp);
    invalidate_cache();
    return retval;
  }

  // warning : point should be transformed before this insert function is called
  Vertex_handle insert_to_star_(const TPoint_2& tp,
                                const bool conditional)
  {
    if(Base::dimension() < 2 || !conditional)
      return this->insert(tp);

    Face_handle fh;
    if(is_conflicted(tp, fh))
      return this->insert(tp, fh);

    else // point is not in conflict with cells incident to center
      return Vertex_handle();
  }

public:
  Vertex_handle insert_to_star(const Point_2& p,
                               const int id,
                               const bool conditional)
  {
    TPoint_2 tp = m_metric.transform(p);
#ifdef ANISO_DEUBG
    std::cout << "insert to star: " << std::endl;
    std::cout << "p: " << p << std::endl;
    std::cout << "m" << std::endl << m_metric.get_transformation() << std::endl;
    std::cout << "tp: " << tp << std::endl;
#endif

    Vertex_handle v = insert_to_star_(tp, conditional);
    if(v != center() && v != Vertex_handle())
      v->info() = id;
    return v;
  }

  Vertex_handle insert_to_star_hole(const Point_2& point,
                                    const Index id,
                                    std::vector<Edge>& boundary_edges,
                                    std::vector<Face_handle>& faces)
  {
    Vertex_handle v = Base::star_hole(m_metric.transform(point),
                                      boundary_edges.begin(),
                                      boundary_edges.end(),
                                      faces.begin(),
                                      faces.end());
    invalidate_cache();
    if(v != center() && v != Vertex_handle())
      v->info() = id;
    else
      std::cout << "problem in insert_to_star_hole..." << std::endl;
    return v;
  }

  int simulate_insert_to_star(const Point_2& p, const int id)
  {
    CGAL_HISTOGRAM_PROFILER("V", this->number_of_vertices());
    TPoint_2 tp = m_metric.transform(p);
    bool found_vertex = false;
    Vertex_handle vh;

#ifndef ANISO_BRUTE_FORCE_SIMULATE_INSERT_TO_STAR
    int li;
    typename Base::Locate_type lt;
    Face_handle fh = Base::locate(tp, lt, li, m_center->face());
    if(lt ==  Base::VERTEX)
    {
      found_vertex = true;
      vh = fh->vertex(li);
    }
#else
    for(typename Base::Finite_vertices_iterator it =  this->finite_vertices_begin();
        (! found_vertex) && (it != this->finite_vertices_end());
        ++it)
    {
      if(p == it->point())
      {
        found_vertex = true;
        vh = it;
      }
    }
#endif
    if(found_vertex) // already in star
    {
      std::cout << "Warning : simulate_insert_to_star re-inserts same point" << std::endl;
      std::cout << this->index_in_star_set() << " p already there: " << vh->info() << std::endl;
      return -1;
    }
    else // in conflict
    {
      Face_handle fh;
      if(is_conflicted(tp, fh))
        return id;
    }
    return -1; // no conflict
  }

  // INFO ABOUT THE VERTICES / FACETS / CELLS
public:
  void print_vertices(const bool print_points = false) const
  {
    std::cout << "\t Star_" << index_in_star_set();
    std::cout << " vertices ("<<(this->number_of_vertices())<<" v.)\t (";
    typename std::set<Vertex_handle> vertices;

    Face_handle_handle fit = finite_incident_faces_begin();
    Face_handle_handle fend = finite_incident_faces_end();
    for(; fit!=fend; fit++)
    {
      Face_handle fh = *fit;
      vertices.insert(fh->vertex(0));
      vertices.insert(fh->vertex(1));
      vertices.insert(fh->vertex(2));
    }
    for(typename std::set<Vertex_handle>::iterator vit = vertices.begin();
        vit != vertices.end(); vit++)
    {
      if(print_points) std::cout << "\n\t";
      std::cout << (*vit)->info() << " ";
      if(print_points) std::cout << " : " << (*vit)->point();
    }

    std::cout << ")\t(";
    typename Base::Finite_vertices_iterator vit = this->finite_vertices_begin();
    typename Base::Finite_vertices_iterator vend = this->finite_vertices_end();
    for(; vit != vend; vit++)
      std::cout << vit->info() << " ";
    std::cout << ")\n";
  }

  void print_faces() const
  {
    std::cout << "\t Star_" << index_in_star_set() << " faces :\n";
    Face_handle_handle fit = finite_incident_faces_begin();
    Face_handle_handle fend = finite_incident_faces_end();
    for(; fit!=fend; ++fit)
    {
      Face_handle fh = *fit;
      std::cout << "\t\t f(" << fh->vertex(0)->info()
                << " " << fh->vertex(1)->info()
                << " " << fh->vertex(2)->info() << ") : "
                << is_inside(fh);
      std::cout << ",\n";
    }
  }

  void face_indices(const Face_handle& fh) const
  {
    std::cout << "Facet[";
    std::vector<int> indices;
    indices.push_back(fh->vertex(0)->info());
    indices.push_back(fh->vertex(1)->info());
    indices.push_back(fh->vertex(2)->info());
    std::sort(indices.begin(), indices.end());
    std::cout << indices[0] << " " << indices[1] << " " << indices[2];
    std::cout << "]";
  }

public:
  void invalidate()
  {
    this->clear();
    m_center_point = CGAL::ORIGIN;
    m_center = Vertex_handle();
    invalidate_cache();
  }

  bool is_valid(const bool verbose = false) const
  {
    return (m_center != Vertex_handle()) && Base::is_valid(verbose);
  }

  void reset(const Point_2& centerpoint,
             const int index,
             const Metric& metric_)
  {
    this->clear();
    m_center_point = centerpoint;
    m_metric = metric_;
    m_center = Base::insert(m_metric.transform(centerpoint));
    m_center->info() = index;
    this->infinite_vertex()->info() = index_of_infinite_vertex;

    invalidate_cache();
  }

  void reset(const Point_2& new_centerpoint)
  {
    Index index = m_center->info();

    this->clear();
    m_center_point = new_centerpoint;
    m_center = Base::insert(m_metric.transform(new_centerpoint));
    m_center->info() = index;
    this->infinite_vertex()->info() = index_of_infinite_vertex;

    invalidate_cache();
  }

  void reset() { reset(center_point()); }

  void update_metric(const Metric& metric_)
  {
    Index id = m_center->info();
    this->clear();
    m_metric = metric_;
    m_center = Base::insert(m_metric.transform(m_center_point));
    m_center->info() = id;
    invalidate_cache();
  }

public:
  Stretched_Delaunay_2(const Criteria* criteria_,
                       const Domain* pdomain)
    :
      Base(*(m_traits = new Traits())),
      m_center_point(CGAL::ORIGIN),
      m_metric(),
      m_pdomain(pdomain),
      m_criteria(new Stretched_criteria(*m_traits, criteria_)),
      m_is_topological_disk(false),
      m_is_valid_topo_disk(false),
      m_is_valid_bbox(false),
      m_bbox_needs_aabb_update(false),
      m_metric_needs_update(false),
      is_cache_dirty(true),
      finite_adjacent_vertices_cache(),
      finite_incident_faces_cache(),
      m_active(true)
  {
    m_center = Vertex_handle();
    m_bbox = m_pdomain->get_bbox(); // in M_euclidean
    build_a_bbox_from_bbox();
    this->infinite_vertex()->info() = index_of_infinite_vertex;
  }

  Stretched_Delaunay_2(const Criteria* criteria_,
                       const Point_2 &centerpoint,
                       const int& index,
                       const Metric& metric_,
                       const Domain* pdomain)
    :
      Base(*(m_traits = new Traits())),
      m_center_point(centerpoint),
      m_metric(metric_),
      m_pdomain(pdomain),
      m_criteria(new Stretched_criteria(*m_traits, criteria_)),
      m_is_topological_disk(false),
      m_is_valid_topo_disk(false),
      m_is_valid_bbox(false),
      m_bbox_needs_aabb_update(false),
      m_metric_needs_update(false),
      is_cache_dirty(true),
      finite_adjacent_vertices_cache(),
      finite_incident_faces_cache(),
      m_active(true)
  {
    m_center = Base::insert(m_metric.transform(centerpoint));
    m_center->info() = index;
    m_bbox = m_pdomain->get_bbox(); // in M_euclidean
    build_a_bbox_from_bbox();
    this->infinite_vertex()->info() = index_of_infinite_vertex;
  }

  ~Stretched_Delaunay_2()
  {
    delete m_traits;
    delete m_criteria;
  }
};

} // Anisotropic_mesh_2
} // CGAL

#endif //CGAL_ANISOTROPIC_MESH_2_STRETCHED_DELAUNAY_2

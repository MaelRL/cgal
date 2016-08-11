#ifndef CGAL_ANISOTROPIC_MESH_3_STRETCHED_DELAUNAY_3
#define CGAL_ANISOTROPIC_MESH_3_STRETCHED_DELAUNAY_3

#include <CGAL/Profile_counter.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Random.h>

#include <CGAL/bbox.h>
#include <CGAL/Metric.h>
#include <CGAL/Criteria.h>
#include <CGAL/Delaunay_traits_3.h>
#include <CGAL/Constrain_surface_3.h>
#include <CGAL/Triangulation_cell_base_with_domain_info_3.h>

#include <CGAL/helpers/combinatorics_helper.h>
#include <CGAL/gl_draw/drawing_helper.h>

#include <CGAL/Timer.h>

#include <iostream>
#include <fstream>
#include <utility>
#include <time.h>

namespace CGAL
{
namespace Anisotropic_mesh_3
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
class Stretched_Delaunay_3
    : public CGAL::Delaunay_triangulation_3<
               Delaunay_traits_3<K, KExact>,
               CGAL::Triangulation_data_structure_3<
                  CGAL::Triangulation_vertex_base_with_info_3<
                    Star_index,
                    Delaunay_traits_3<K, KExact> >,
                  CGAL::Triangulation_cell_base_with_domain_info_3<
                    Delaunay_traits_3<K, KExact>,
                    Constrain_surface_3<K> > > >
{
  typedef Stretched_Delaunay_3<K, KExact>                       Self;
public:
  typedef K                                                     Kernel;
  typedef Delaunay_traits_3<K, KExact>                          Traits;
  typedef CGAL::Triangulation_data_structure_3<
            CGAL::Triangulation_vertex_base_with_info_3<
              Star_index,
              Traits >,
            CGAL::Triangulation_cell_base_with_domain_info_3<
              Traits, Constrain_surface_3<K> >
          >                                                     DS;
  typedef CGAL::Delaunay_triangulation_3<Traits, DS>            Base;

  typedef int Index;

  typedef Constrain_surface_3<K>                    Constrain_surface;
  typedef Metric_base<K, KExact>                    Metric;
  typedef typename K::FT                            FT;
  typedef typename K::Point_3                       Point_3;
  typedef typename K::Point_3                       TPoint_3; //Transformed_point_3
  typedef typename K::Sphere_3                      Sphere;
  typedef typename K::Vector_3                      Vector_3;
  typedef typename K::Line_3                        Line_3;
  typedef typename K::Ray_3                         Ray_3;
  typedef typename K::Plane_3                       Plane_3;
  typedef typename Base::Segment                    Segment;
  typedef typename Base::Triangle                   Triangle;
  typedef typename Base::Tetrahedron                Tetrahedron;
  typedef typename Base::Vertex                     Vertex;
  typedef typename Base::Cell                       Cell;
  typedef typename Base::Edge                       Edge;
  typedef typename Base::Facet                      Facet;
  typedef typename DS::Vertex_handle                Vertex_handle;
  typedef typename DS::Cell_handle                  Cell_handle;
  typedef ::CGAL::Anisotropic_mesh_3::Stretched_criteria<K, KExact>
                                                    Stretched_criteria;
  typedef Criteria_base<K>                          Criteria;
  typedef CGAL::Bbox<K>                             Bbox;

  typedef std::vector<Vertex_handle>                Vertex_handle_vector;
  typedef std::vector<Facet>                        Facet_vector;
  typedef std::vector<Cell_handle>                  Cell_handle_vector;
  typedef std::set<Facet>                           Facet_set;
  typedef std::vector<Point_3>                      Point_vector;
  typedef typename Facet_set::iterator              Facet_set_iterator;
  typedef typename Facet_vector::iterator           Facet_handle;
  typedef typename Cell_handle_vector::iterator     Cell_handle_handle;
  typedef typename Vertex_handle_vector::iterator   Vertex_handle_handle;

  typedef typename KExact::Point_3                  Exact_Point_3;
  typedef typename KExact::Point_3                  Exact_TPoint_3;
  typedef typename KExact::Vector_3                 Exact_Vector_3;
  typedef CGAL::Cartesian_converter<K, KExact>      To_exact;
  typedef CGAL::Cartesian_converter<KExact, K>      Back_from_exact;

private:
  static const Index index_of_infinite_vertex = -10;

  Traits* m_traits;

  Point_3 m_center_point; // before transformation by m_metric.transform
  Vertex_handle m_center; // after transformation by m_metric.transform
  bool m_is_surface_star;

  To_exact to_exact;
  Back_from_exact back_from_exact;
  Metric m_metric;
  const Constrain_surface* m_pConstrain;
  Stretched_criteria* m_criteria;

  mutable bool m_is_in_3D_mesh; // if we care about cells, we need 3D pick_valid and volume bboxes
  mutable bool m_is_topological_disk;
  mutable bool m_is_valid_topo_disk;
  mutable Bbox m_bbox;
  mutable bool m_is_valid_bbox;
  mutable bool m_bbox_needs_aabb_update;
  mutable bool m_metric_needs_update;

  mutable bool is_cache_dirty;
  mutable Facet_set restricted_facets_cache; // restricted and incident to m_center
  mutable Cell_handle_vector incident_cells_cache;
  mutable Cell_handle_vector finite_incident_cells_cache;
  mutable Vertex_handle_vector finite_adjacent_vertices_cache;
  mutable Vertex_handle_vector finite_adjacent_restricted_vertices_cache;

public:
  mutable bool m_active; // used for removing points at the end of the process

public:
  static Index infinite_vertex_index() { return index_of_infinite_vertex; }

  bool is_infinite() const
  {
    return (infinite_vertex_index() == index_in_star_set());
  }

  inline bool is_infinite_vertex(const Vertex_handle &v) const
  {
    return Base::is_infinite(v);
  }

  inline bool is_infinite(const Cell_handle& cell) const
  {
    return Base::is_infinite(cell);
  }

public:
  Point_3& center_point()             { return m_center_point; }
  const Point_3& center_point() const { return m_center_point; }

  Index& index_in_star_set()             { return m_center->info(); }
  const Index& index_in_star_set() const { return m_center->info(); }

  bool& is_surface_star()             { return m_is_surface_star; }
  const bool& is_surface_star() const { return m_is_surface_star; }

  bool& is_in_3D_mesh() { return m_is_in_3D_mesh; }
  const bool& is_in_3D_mesh() const { return m_is_in_3D_mesh; }

public:
  Vertex_handle center() const      { return m_center; }
  Metric metric() const             { return m_metric; }
  Metric& metric()                  { return m_metric; }
  Traits* traits() const            { return m_traits; }
  const Constrain_surface* constrain() const { return m_pConstrain; }
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

  bool& bbox_needs_aabb_update(){ return m_bbox_needs_aabb_update; }
  const bool& bbox_needs_aabb_update() const { return m_bbox_needs_aabb_update; }
  bool& metric_needs_update(){ return m_metric_needs_update; }
  const bool& metric_needs_update() const { return m_metric_needs_update; }

public:
  template<typename I>
  bool has_vertex(const I i) const
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

  template<typename I>
  bool has_vertex(const I i, Vertex_handle& v) const
  {
    Vertex_handle_handle vit = finite_adjacent_vertices_begin();
    Vertex_handle_handle vend = finite_adjacent_vertices_end();
    for(; vit != vend; vit++)
    {
      Vertex_handle vh = *vit;
      if(vit->info() == i)
      {
        v = vh;
        return true;
      }
    }
    return false;
  }

  //this function checks if f is incident to the center
  bool is_in_star(const Facet& f) const
  {
    for(int i = 1; i < 4; i++)
    {
      Vertex_handle v = f.first->vertex((f.second+i)%4);
      if(v->info() == m_center->info())
        return true;
    }
    return false;
  }

  bool has_facet(const Facet &facet) const
  {
    Facet_set_iterator fi = restricted_facets_begin();
    Facet_set_iterator fend = restricted_facets_end();
    for(; fi != fend; fi++)
      if(Facet_ijk(*fi) == Facet_ijk(facet))
        return true;
    return false;
  }

  bool has_facet(const Facet_ijk& f, Facet& facet) const
  {
    boost::array<int, 3> dids;
    Facet_set_iterator fi = restricted_facets_begin();
    Facet_set_iterator fiend = restricted_facets_end();
    for (; fi != fiend; fi++)
    {
      for (int i=0; i<3; i++)
        dids[i] = fi->first->vertex((fi->second+i+1) % 4)->info();
      std::sort(dids.begin(), dids.end());
      if(std::equal(dids.begin(), dids.end(), f.vertices().begin()))
      {
        facet = *fi;
        return true;
      }
    }
    return false;
  }

  bool has_facet(int *cids, Facet &facet) const
  {
    int dids[3];
    Facet_set_iterator fi = restricted_facets_begin();
    Facet_set_iterator fiend = restricted_facets_end();
    for (; fi != fiend; fi++)
    {
      for (int i=0; i<3; i++)
        dids[i] = fi->first->vertex((fi->second+i+1) % 4)->info();
      if (is_same_ids<3>(cids, dids))
      {
        facet = *fi;
        return true;
      }
    }
    return false;
  }

  bool has_cell(const Cell_handle &cell) const
  {
    int cids[4], dids[4];
    for (int i = 0; i < 4; i++)
      cids[i] = cell->vertex(i)->info();

    Cell_handle_handle ci = star_cells_begin();
    Cell_handle_handle cend = star_cells_end();
    for (; ci != cend; ci++)
    {
      for (int i = 0; i < 4; i++)
        dids[i] = (*ci)->vertex(i)->info();
      if (is_same_ids<4>(cids, dids))
        return true;
    }
    return false;
  }

  bool has_cell(int *vertices, Cell_handle &cell) const
  {
    int dids[4];
    Cell_handle_handle ci = star_cells_begin();
    Cell_handle_handle cend = star_cells_end();
    for (; ci != cend; ci++) {
      for (int i = 0; i < 4; i++)
        dids[i] = (*ci)->vertex(i)->info();
      if (is_same_ids<4>(vertices, dids))
      {
        cell = *ci;
        return true;
      }
    }
    return false;
  }

  bool has_cell(const Cell_ijkl& c, Cell_handle& cell) const
  {
    boost::array<int, 4> dids;
    Cell_handle_handle ci = star_cells_begin();
    Cell_handle_handle cend = star_cells_end();
    for (; ci != cend; ci++)
    {
      for (int i=0; i<4; i++)
        dids[i] = (*ci)->vertex(i)->info();
      std::sort(dids.begin(), dids.end());
      if(std::equal(dids.begin(), dids.end(), c.vertices().begin()))
      {
        cell = *ci;
        return true;
      }
    }
    return false;
  }

  // CACHE RELATED FUNCTIONS
public:
  void set_facet_cache(const Facet& facet, const Point_3& p) const
  {
    Facet mf = this->mirror_facet(facet);
    facet.first->set_facet_surface_center(facet.second, p);
    mf.first->set_facet_surface_center(mf.second, p);
    facet.first->set_facet_visited(facet.second);
    mf.first->set_facet_visited(mf.second);
  }

  inline void update_star_caches() const
  {
    if(!is_cache_dirty)
      return;

    CGAL_PROFILER("[update_star_caches]");
    incident_cells_cache.clear();
    finite_incident_cells_cache.clear();
    restricted_facets_cache.clear();
    finite_adjacent_vertices_cache.clear();

    // update neighboring vertices
    std::back_insert_iterator<Vertex_handle_vector>
        finite_adjacent_vertices_insertor(finite_adjacent_vertices_cache);
    Base::finite_adjacent_vertices(m_center, finite_adjacent_vertices_insertor);

    if(Base::dimension() == 2)
    {
      typename Base::Finite_facets_iterator fit = this->finite_facets_begin();
      typename Base::Finite_facets_iterator fend = this->finite_facets_end();
      for(; fit != fend; ++fit)
        if(is_in_star(*fit) && is_restricted(*fit))
          restricted_facets_cache.insert(this->make_canonical(*fit));
    }
    else if(Base::dimension() > 2)
    {
      // update neighboring cells
      std::back_insert_iterator<Cell_handle_vector> incident_cells_insertor(incident_cells_cache);
      Base::incident_cells(m_center, incident_cells_insertor);

      std::back_insert_iterator<Cell_handle_vector>
        finite_incident_cells_insertor(finite_incident_cells_cache);
      Base::finite_incident_cells(m_center, finite_incident_cells_insertor);

      Cell_handle_handle ci = finite_incident_cells_cache.begin();
      Cell_handle_handle cend = finite_incident_cells_cache.end();
      for (; ci != cend; ci++)
      {
        int center_index = (*ci)->index(m_center);
        for(int i = 0; i < 4; i++)
        {
          if(i == center_index)
            continue;
          Facet f = this->make_canonical(Facet(*ci, i));
          Point_3 p;
          if(is_restricted(f, p, true)) //updates surface delaunay center, if needed
            restricted_facets_cache.insert(f);
        }
      }
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

  template<typename OutputIndexIterator>
  void finite_star_vertices(OutputIndexIterator oit) const
  {
    std::vector<Vertex_handle> vertices;
    finite_adjacent_vertices(m_center, std::back_inserter(vertices));
    Vertex_handle_handle it = vertices.begin();
    Vertex_handle_handle itend = vertices.end();
    for(; it != itend; it++)
      *oit++ = (*it)->info();
  }

  inline Facet_set_iterator restricted_facets_begin() const
  {
    update_star_caches();
    return restricted_facets_cache.begin();
  }
  inline Facet_set_iterator restricted_facets_end() const
  {
    return restricted_facets_cache.end();
  }
  inline Cell_handle_handle star_cells_begin() const
  {
    update_star_caches();
    return incident_cells_cache.begin();
  }
  inline Cell_handle_handle star_cells_end() const
  {
    return incident_cells_cache.end();
  }
  inline Cell_handle_handle finite_star_cells_begin() const
  {
    update_star_caches();
    return finite_incident_cells_cache.begin();
  }
  inline Cell_handle_handle finite_star_cells_end() const
  {
    return finite_incident_cells_cache.end();
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
  inline int get_finite_star_cell_count() const
  {
    update_star_caches();
    return (int)finite_incident_cells_cache.size();
  }
  inline bool is_boundary_star() const
  {
    update_star_caches();
    return (restricted_facets_cache.size() > 0);
  }

  void compute_adjacent_restricted_vertices()
  {
    if(!is_boundary_star())
      return;

    finite_adjacent_restricted_vertices_cache.clear();

    std::vector<Vertex_handle> tmp_vertices;

    Facet_set_iterator it = restricted_facets_begin();
    Facet_set_iterator end = restricted_facets_end();
    for(; it!=end; ++it)
    {
      for(std::size_t i=1; i<4; ++i)
      {
        Vertex_handle vh = it->first->vertex((it->second + i)%4);
        if(!vh->visited_for_vertex_extractor)
        {
          vh->visited_for_vertex_extractor = true;
          finite_adjacent_restricted_vertices_cache.push_back(vh);
          tmp_vertices.push_back(vh);
        }
      }
    }

    for(std::size_t i=0; i<tmp_vertices.size(); ++i)
      tmp_vertices[i]->visited_for_vertex_extractor = false;

  }

  inline Vertex_handle_handle finite_adjacent_restricted_vertices_begin() const
  {
    // should put a boolean for safety...
    return finite_adjacent_restricted_vertices_cache.begin();
  }

  inline Vertex_handle_handle finite_adjacent_restricted_vertices_end() const
  {
    return finite_adjacent_restricted_vertices_cache.end();
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
    if(!is_surface_star() || this->dimension() < 3)
    {
      m_is_topological_disk = false;
      m_is_valid_topo_disk = true;
      return;
    }

    Facet_set_iterator fit = this->restricted_facets_begin();
    Facet_set_iterator fend = this->restricted_facets_end();
    if(fit == fend) // no restricted facets
    {
      m_is_topological_disk = false;
      m_is_valid_topo_disk = true;
      return;
    }

    //for each vertex V of the star S, the map contains its n neighbors (except
    //the center of the star) that belong to restricted facets. n can have any value
    //but if we have n!=2, the star is not a topological disk (each edge must be
    //incident to exactly 2 restricted facets), thus we simply keep the first 2
    //and exit if n>2 or n<2.
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
            // edge [center,v1] has more than 2 incident restricted facets
            return false;
          }
          vhn.m_sn = v2; // setting v1's second neighbor
        }
        return true;
      }

      Cycle_map_insertor(Cycle_map& cmap_) : m_cmap(cmap_) { }
    };

    Cycle_map_insertor cm_insertor(cycle_map);

    //fill the map from the restricted facets
    for(; fit != fend; fit++)
    {
      Facet f = *fit;
      std::vector<Vertex_handle> vns; // the 2 vertices of f that are not the center of the star
      for(int i = 0; i < 4; i++)
      {
        Vertex_handle v = f.first->vertex(i);
        if(i != f.second && v != m_center)
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

    //check that all the vertices have 2 neighbors that live on restricted facets
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

  void update_bbox(const bool verbose = false) const
  {
    if(m_is_valid_bbox)
      return;

    Bbox_3 old_bbox = m_bbox;

    if(verbose)
      std::cout << "Update Bbox...";
    CGAL_PROFILER("[update_bbox]");

    if(this->dimension() < 3)
    {
      m_bbox = m_pConstrain->get_bbox(); // must be found by every request to aabb_tree
    }
    else if(!m_is_in_3D_mesh)
    {
      m_bbox = this->surface_bbox();
      if(!is_topological_disk())
        m_bbox = m_bbox + this->volume_bbox();
    }
    else
      m_bbox = this->volume_bbox();

    m_is_valid_bbox = true;

#ifndef NO_USE_AABB_TREE_OF_BBOXES
   //checking if the bbox has grown
    if(m_bbox.xmin() < old_bbox.xmin() || m_bbox.xmax() > old_bbox.xmax() ||
       m_bbox.ymin() < old_bbox.ymin() || m_bbox.ymax() > old_bbox.ymax() ||
       m_bbox.zmin() < old_bbox.zmin() || m_bbox.zmax() > old_bbox.zmax())
    {
#ifdef DEBUG_UPDATE_AABB_TREE
      std::cout << "growing bbox in star : " << index_in_star_set() << std::endl;
      std::cout << m_bbox.xmin() << " " << m_bbox.xmax() << " " << m_bbox.ymin();
      std::cout << " " << m_bbox.ymax() << " " << m_bbox.zmin() << " " << m_bbox.zmax() << std::endl;
#endif
      m_bbox_needs_aabb_update = true;
   }
#endif

    if(verbose)
      std::cout << "done." << std::endl;
  }

  Bbox surface_bbox() const // compute bbox of incident surface Delaunay balls
  {
    typename Traits::Compute_squared_distance_3 csd
      = m_traits->compute_squared_distance_3_object();

    Bbox bb = m_center->point().bbox();
    Facet_set_iterator fi = restricted_facets_begin();
    Facet_set_iterator fend = restricted_facets_end();
    for(; fi != fend; fi++)
    {
      Point_3 p;
#ifdef ANISO_USE_EXACT
      compute_exact_dual_intersection(*fi, p);
#else
      compute_dual_intersection(*fi, p);
#endif
      TPoint_3 tp = m_metric.transform(p);
      FT squared_radius = csd(tp, m_center->point());
      Sphere s(tp, squared_radius);
      bb = bb + s.bbox();
    }
    return m_metric.inverse_transform(bb);
  }

  Bbox volume_bbox() const // compute bbox of incident cells circumspheres
  {
    typename Traits::Compute_squared_distance_3 csd
      = m_traits->compute_squared_distance_3_object();

    Bbox bb = m_center->point().bbox();
    Cell_handle_handle ci = finite_star_cells_begin();
    Cell_handle_handle cend = finite_star_cells_end();
    for (; ci != cend; ci++)
    {
      Cell_handle ch = *ci;
      TPoint_3 cc = ch->circumcenter(*m_traits);
      FT squared_radius = csd(cc, m_center->point());

#ifdef ANISO_DEBUG_BBOX
      std::cout << "computing circumcenter with points: " << std::endl;
      std::cout << ch->vertex(0)->info() << " " << ch->vertex(0)->point() << std::endl;
      std::cout << ch->vertex(1)->info() << " " << ch->vertex(1)->point() << std::endl;
      std::cout << ch->vertex(2)->info() << " " << ch->vertex(2)->point() << std::endl;
      std::cout << ch->vertex(3)->info() << " " << ch->vertex(3)->point() << std::endl;
      std::cout << "cc: " << cc << " radius is: " << squared_radius << std::endl;
#endif

#ifdef ANISO_APPROXIMATE_SPHERE_BBOX
      double radius = sqrt(squared_radius);
      Bbox bbs(cc.x()-radius, cc.y()-radius, cc.z()-radius,
              cc.x()+radius, cc.y()+radius, cc.z()+radius);
      bb = bb + bbs;
#else
     Sphere s(cc, squared_radius);
     bb = bb + s.bbox();
#endif
    }
    return m_metric.inverse_transform(bb);
  }

  bool is_inside_bbox(const Point_3& p) const
  {
    update_bbox();
    return p.x() >= m_bbox.xmin() && p.x() <= m_bbox.xmax()
        && p.y() >= m_bbox.ymin() && p.y() <= m_bbox.ymax()
        && p.z() >= m_bbox.zmin() && p.z() <= m_bbox.zmax();
  }

public:
  Facet make_canonical(const Facet& f) const
  {
    //returns always the same (c,i) for f and its mirror facet
    if(this->dimension() < 3)
      return f;

    Facet f2 = this->mirror_facet(f);
    Index i_f = f.first->vertex(f.second)->info();
    Index i_f2 = f2.first->vertex(f2.second)->info();
    return (i_f > i_f2) ? f : f2;
    // note : infinite vertex has index -10
    // we choose the one with the greatest index to get the facet on the finite side
  }

  bool is_facet_visited(const Facet& f) const
  {
    if(this->dimension() < 3)
      return f.first->is_facet_visited(f.second);

    Facet f2 = this->mirror_facet(f);
    return f.first->is_facet_visited(f.second) && f2.first->is_facet_visited(f2.second);
  }

public:
  inline FT compute_volume(const Cell_handle &cell) const
  {
    return m_criteria->compute_volume(
          cell->vertex(0)->point(), cell->vertex(1)->point(),
          cell->vertex(2)->point(), cell->vertex(3)->point());
  }

  inline FT compute_volume(const Facet &facet) const
  {
    return m_criteria->compute_volume(
      facet.first->vertex((facet.second + 1) % 4)->point(),
      facet.first->vertex((facet.second + 2) % 4)->point(),
      facet.first->vertex((facet.second + 3) % 4)->point());
  }

  inline FT compute_radius_edge_ratio_overflow(const Cell_handle &cell) const
  {
    return m_criteria->radius_edge_ratio_overflow(
      cell->vertex(0)->point(), cell->vertex(1)->point(),
      cell->vertex(2)->point(), cell->vertex(3)->point());
  }

  inline FT compute_radius_edge_ratio_overflow(const Facet &facet) const
  {
    return m_criteria->radius_edge_ratio_overflow(
      facet.first->vertex((facet.second + 1) % 4)->point(),
      facet.first->vertex((facet.second + 2) % 4)->point(),
      facet.first->vertex((facet.second + 3) % 4)->point());
  }

  inline FT compute_squared_radius_edge_ratio(const Facet& facet) const
  {
    return m_criteria->compute_squared_radius_edge_ratio(
      facet.first->vertex((facet.second + 1) % 4)->point(),
      facet.first->vertex((facet.second + 2) % 4)->point(),
      facet.first->vertex((facet.second + 3) % 4)->point());
  }

  inline FT compute_squared_radius_edge_ratio(const TPoint_3& tp1,
                                              const TPoint_3& tp2,
                                              const TPoint_3& tp3) const
  {
    return m_criteria->compute_squared_radius_edge_ratio(tp1, tp2, tp3);
  }

  inline FT compute_squared_radius_edge_ratio(const Cell_handle ch) const
  {
    return m_criteria->compute_squared_radius_edge_ratio(
          ch->vertex(0)->point(), ch->vertex(1)->point(),
          ch->vertex(2)->point(), ch->vertex(3)->point());
  }

  inline FT compute_squared_radius_edge_ratio(const TPoint_3& tp1,
                                              const TPoint_3& tp2,
                                              const TPoint_3& tp3,
                                              const TPoint_3& tp4) const
  {
    return m_criteria->compute_squared_radius_edge_ratio(tp1, tp2, tp3, tp4);
  }

  inline FT compute_circumradius_overflow(const Cell_handle &cell) const
  {
    return m_criteria->circumradius_overflow(
      cell->vertex(0)->point(), cell->vertex(1)->point(),
      cell->vertex(2)->point(), cell->vertex(3)->point());
  }

  inline FT compute_circumradius_overflow(const Facet &facet) const
  {
    // surface Delaunay ball radius
    Point_3 p;
    compute_dual_intersection(facet,p);
    TPoint_3 tp = m_metric.transform(p);
    TPoint_3 tp2 = facet.first->vertex((facet.second + 1) % 4)->point();
    FT sqr = m_traits->compute_squared_distance_3_object()(tp, tp2);
    return sqr - m_criteria->criteria->squared_facet_circumradius;
  }

  inline FT compute_squared_circumradius(const Facet &facet) const
  {
    return m_criteria->compute_squared_circumradius(
      facet.first->vertex((facet.second + 1) % 4)->point(),
      facet.first->vertex((facet.second + 2) % 4)->point(),
      facet.first->vertex((facet.second + 3) % 4)->point());
  }

  inline FT compute_squared_circumradius(const TPoint_3& p1,
                                         const TPoint_3& p2,
                                         const TPoint_3& p3) const
  {
    return m_criteria->compute_squared_circumradius(p1, p2, p3);
  }

  inline FT compute_sliverity_overflow(const Cell_handle &cell) const
  {
    return m_criteria->sliverity_overflow(
      cell->vertex(0)->point(), cell->vertex(1)->point(),
      cell->vertex(2)->point(), cell->vertex(3)->point());
  }

  inline FT compute_squared_circumradius(const Cell_handle ch) const
  {
    return m_criteria->compute_squared_circumradius(
          ch->vertex(0)->point(), ch->vertex(1)->point(),
          ch->vertex(2)->point(), ch->vertex(3)->point());
  }

  inline FT compute_squared_circumradius(const TPoint_3& p1,
                                         const TPoint_3& p2,
                                         const TPoint_3& p3,
                                         const TPoint_3& p4) const
  {
    return m_criteria->compute_squared_circumradius(p1, p2, p3, p4);
  }

  inline FT compute_sliverity_overflow(const Point_3& p1,
                                       const Point_3& p2,
                                       const Point_3& p3,
                                       const Point_3& p4) const
  {
    return m_criteria->sliverity_overflow(
      m_metric.transform(p1), m_metric.transform(p2),
      m_metric.transform(p3), m_metric.transform(p4));
  }

  inline FT compute_element_quality(const Cell_handle & cell) const
  {
    return m_criteria->element_quality(
      cell->vertex(0)->point(), cell->vertex(1)->point(),
      cell->vertex(2)->point(), cell->vertex(3)->point());
  }

  inline FT compute_element_quality(const TPoint_3& p1,
                                    const TPoint_3& p2,
                                    const TPoint_3& p3,
                                    const TPoint_3& p4) const
  {
    return m_criteria->element_quality(p1, p2, p3, p4);
  }


  inline FT compute_element_quality(const TPoint_3& p1,
                                    const TPoint_3& p2,
                                    const TPoint_3& p3) const
  {
    return m_criteria->element_quality(p1, p2, p3);
  }

  inline bool is_sliver(const TPoint_3& p1,
                        const TPoint_3& p2,
                        const TPoint_3& p3,
                        const TPoint_3& p4) const
  {
    return (compute_sliverity_overflow(p1, p2, p3, p4) > 0.);
  }

public:
  inline bool is_inside(const Cell_handle &cell) const
  {
    if(this->is_infinite(cell))
      return false;
    return cell->template is_inside<Self>(*this, *constrain(), *traits());
  }

#ifdef ANISO_USE_INSIDE_EXACT
  inline bool is_inside_exact(const Cell_handle &cell) const
  {
    if(this->is_infinite(cell))
      return false;
    return cell->template is_inside_exact<Self>(*this, *constrain(), *traits());
  }
#endif

public:
  unsigned int nb_restricted_facets() const
  {
    unsigned int count = 0;
    Facet_set_iterator fit = restricted_facets_begin();
    Facet_set_iterator fend = restricted_facets_end();
    for(; fit != fend; ++fit)
      count++;
    return count;
  }

  bool is_restricted(const Facet& f,
                     Point_3& surface_delaunay_ball_center,
                     const bool compute_intersection) const
  {
    if(Base::is_infinite(f))
      return false;
    else if(Base::dimension() == 2)
      return is_restricted_2_in_3(f, surface_delaunay_ball_center, true/*=exact*/);
    else if(Base::dimension() == 3)
    {
      Cell_handle c1 = f.first;
      Cell_handle c2 = f.first->neighbor(f.second);
#ifdef ANISO_USE_INSIDE_EXACT
      bool inside1 = is_inside_exact(c1);
      bool inside2 = is_inside_exact(c2);
#else
      bool inside1 = is_inside(c1);
      bool inside2 = is_inside(c2);
#endif
      if((inside1 && !inside2) || (!inside1 && inside2))
      {
        if(!compute_intersection)
          return true;
#ifdef ANISO_USE_EXACT
        if(!compute_exact_dual_intersection(f, surface_delaunay_ball_center))
          std::cout << "Warning : is_restricted exact can't find an intersection point.\n";
#else
        if(!compute_dual_intersection(f, surface_delaunay_ball_center))
          std::cout << "Warning : is_restricted can't find an intersection point.\n";
#endif
        return true;
      }
      else return false;
    }
    return false;
  }

  bool is_restricted(const Facet& f,
                     Point_3& surface_delaunay_ball_center) const
  {
    return is_restricted(f, surface_delaunay_ball_center, true);
  }

  bool is_restricted(const Facet& f) const
  {
    Point_3 c;
    return is_restricted(f,c,false/*do not compute intersection*/);
  }

public:
  bool is_in_surface_delaunay_ball(const TPoint_3& tp,
                                   const Facet& f,
                                   Cell_handle& in_which_cell) const
  {
    CGAL_PROFILER("[conflict surfaceDB test]");
    if(Base::dimension() < 3)
      return true;

    typename Traits::Compute_squared_distance_3 csd
      = m_traits->compute_squared_distance_3_object();
    Point_3 center;
#ifdef ANISO_USE_EXACT
    compute_exact_dual_intersection(f, center);
#else
    compute_dual_intersection(f, center);
#endif
    TPoint_3 tcenter = m_metric.transform(center);

    FT sq_radius = csd(tcenter, m_center->point()); // M_star
    if(csd(tcenter, tp) >= sq_radius)
      return false;

    Cell_handle c1 = f.first;
    Cell_handle c2 = f.first->neighbor(f.second);

#ifdef ANISO_DEBUG_DELAUNAY_BALLS
    std::cout << "is in surf dball @ " << index_in_star_set() << " fsecond: " << f.second << std::endl;
    std::cout << c1->vertex(0)->point() << " || " << c1->vertex(0)->info() << std::endl;
    std::cout << c1->vertex(1)->point() << " || " << c1->vertex(1)->info() << std::endl;
    std::cout << c1->vertex(2)->point() << " || " << c1->vertex(2)->info() << std::endl;
    std::cout << c1->vertex(3)->point() << " || " << c1->vertex(3)->info() << std::endl;
    std::cout << Base::is_valid(c1, true) << std::endl;

    std::cout << c2->vertex(0)->point() << " || " << c2->vertex(0)->info() << std::endl;
    std::cout << c2->vertex(1)->point() << " || " << c2->vertex(1)->info() << std::endl;
    std::cout << c2->vertex(2)->point() << " || " << c2->vertex(2)->info() << std::endl;
    std::cout << c2->vertex(3)->point() << " || " << c2->vertex(3)->info() << std::endl;
    std::cout << Base::is_valid(c2, true) << std::endl;
    std::cout << "end" << std::endl;
#endif

    if(Base::side_of_sphere(c1, tp) == CGAL::ON_BOUNDED_SIDE)
    {
      in_which_cell = c1;
      return true;
    }
    else if(Base::side_of_sphere(c2, tp) == CGAL::ON_BOUNDED_SIDE)
    {
      in_which_cell = c2;
      return true;
    }
    return false;
  }

  inline bool is_in_a_surface_delaunay_ball(const TPoint_3& tp,
                                           Cell_handle& in_which_cell) const
  {
    CGAL_PROFILER("[is_in_a_surface_delaunay_ball]");
    Facet_set_iterator fi = restricted_facets_begin();
    Facet_set_iterator fend = restricted_facets_end();
    for(; fi != fend; fi++)
      if(is_in_surface_delaunay_ball(tp, *fi, in_which_cell))
        return true;
    return false;
  }

  bool is_in_delaunay_ball(const TPoint_3& tp,
                           const Cell_handle& c) const
  {
    CGAL_PROFILER("[conflict volumeDB test]");
    return (Base::side_of_sphere(c, tp) == CGAL::ON_BOUNDED_SIDE);
  }

  inline bool is_in_a_volume_delaunay_ball(const TPoint_3& tp,
                                           Cell_handle& in_which_cell) const
  {
    CGAL_PROFILER("[is_in_a_volume_delaunay_ball]");
    if(this->dimension() < 3)
      return true;
    Cell_handle_handle ci = star_cells_begin();
    Cell_handle_handle cend = star_cells_end();
    for(; ci!=cend; ci++)
    {
      if(is_in_delaunay_ball(tp, *ci))
      {
        in_which_cell = *ci;
        return true;
      }
    }
    return false;
  }

  inline bool is_conflicted(const TPoint_3 &tp,
                            Cell_handle& in_which_cell) const
  {
    if(Base::dimension() < 3)
      return true;
    else if(!is_inside_bbox(m_metric.inverse_transform(tp)))
      return false;
    else if(!m_is_in_3D_mesh)
    {
      if(is_topological_disk())
        return is_in_a_surface_delaunay_ball(tp, in_which_cell);
      else
        return (is_in_a_volume_delaunay_ball(tp, in_which_cell) ||
                is_in_a_surface_delaunay_ball(tp, in_which_cell) );
    }
    else
      return is_in_a_volume_delaunay_ball(tp, in_which_cell);
  }

  template<typename BFacetsOutputIterator,
           typename CellsOutputIterator,
           typename IFacetsOutputIterator>
  bool find_conflicts(const Point_3& p, BFacetsOutputIterator bfoit /*boundary*/,
                                        CellsOutputIterator coit /*cells*/,
                                        IFacetsOutputIterator ifoit /*internal*/) const
  {
    int dim = Base::dimension();
    if(dim <= 1)
      return false;
    else if(dim == 2)
    {
      Facet_set_iterator fit = restricted_facets_begin();
      Facet_set_iterator fend = restricted_facets_end();
      for(; fit != fend; fit++)
        if(is_restricted(*fit))
          *bfoit++ = *fit;
      return true;
    }
    else
    {
      TPoint_3 tp = m_metric.transform(p);
      Cell_handle ch; // this cell has already been computed previously and could be kept in stars_czones memory todo
      if(!is_conflicted(tp, ch))
      {
        std::cout << "no conflict in find_conflict.........(bad)" << std::endl;
        return false;
      }

      Base::find_conflicts(tp, ch, bfoit, coit, ifoit);
      return true;
    }
    return false;
  }

  template<typename BoundaryFacetsOutputIterator>
  bool find_conflicts(const Point_3& p, BoundaryFacetsOutputIterator oit) const
  {
    return find_conflicts(p, oit, Emptyset_iterator(), Emptyset_iterator());
  }

  bool can_form_sliver(const TPoint_3 &p, const FT &squared_radius_bound) const
  {
    Cell_handle_handle ci = star_cells_begin();
    Cell_handle_handle cend = star_cells_end();
    for (; ci != cend; ci++)
    {
      if (Base::side_of_sphere((*ci), p) == CGAL::ON_UNBOUNDED_SIDE)
        continue;
      Cell_handle c = *ci;
      TPoint_3 points[3];
      TPoint_3 centerpoint = m_center->point();
      int k = 0;
      for (int i = 0; i < 4; i++)
      {
        if (c->vertex(i) != m_center)
          points[k++] = c->vertex(i)->point();
      }
      assert(k == 3);
      for (int i = 0; i < 3; i++)
      {
        FT sliver_overflow = m_criteria->sliverity_overflow(
          centerpoint, points[i], points[(i + 1) % 3], p);
        if (sliver_overflow <= 0.0)
          continue;
        if (m_criteria->compute_volume(centerpoint, points[i], points[(i + 1) % 3], p) < 1e-7)
          continue; //should be return true? todo
        if (m_criteria->compute_squared_circumradius(centerpoint, points[i],
          points[(i + 1) % 3], p) >= squared_radius_bound)
          continue;
        return true;
      }
    }
    return false;
  }

  inline void get_inverse_transformed_points(Point_3 *ps, const Facet &facet) const
  {
    int k = 0;
    for (int i = 0; i < 4; i++)
      if (i != facet.second)
        ps[k++] = m_metric.inverse_transform(facet.first->vertex(i)->point());
  }

  inline void get_points(TPoint_3 *ps, const Facet &facet) const
  {
    int k = 0;
    for (int i = 0; i < 4; i++)
      if (i != facet.second)
        ps[k++] = facet.first->vertex(i)->point();
  }

  inline void get_points(TPoint_3 *ps, Cell_handle cell) const
  {
    int k = 0;
    for (int i = 0; i < 4; i++)
      ps[k++] = cell->vertex(i)->point();
  }

  TPoint_3 compute_circumcenter(const TPoint_3& p0,
                                const TPoint_3& p1,
                                const TPoint_3& p2) const
  {
#ifdef ANISO_USE_CC_EXACT
    return back_from_exact(m_traits->construct_circumcenter_3_object()(
                to_exact(p0), to_exact(p1), to_exact(p2)));
#else
    return m_traits->construct_circumcenter_3_object()(p0, p1, p2);
#endif
  }

  TPoint_3 compute_circumcenter(const Facet &facet) const
  {
#ifdef ANISO_USE_CC_EXACT
    return back_from_exact(compute_exact_circumcenter(facet));
#else
    return compute_exact_circumcenter(facet);
#endif
  }

  Exact_TPoint_3 compute_exact_circumcenter(const Facet& facet) const
  {
    TPoint_3 ps[3];
    get_points(ps, facet);
#ifdef ANISO_USE_CC_EXACT
    return m_traits->construct_circumcenter_3_object()(
               to_exact(ps[0]), to_exact(ps[1]), to_exact(ps[2]));
#else
    return m_traits->construct_circumcenter_3_object()(ps[0], ps[1], ps[2]);
#endif
  }

  Exact_TPoint_3 compute_exact_circumcenter(Cell_handle cell) const
  {
#ifdef ANISO_USE_CC_EXACT
    TPoint_3 ps[4];
    get_points(ps, cell);
    return m_traits->construct_circumcenter_3_object()(
             to_exact(ps[0]), to_exact(ps[1]), to_exact(ps[2]), to_exact(ps[3]));
#else
    return cell->circumcenter(*m_traits);
#endif
  }

  TPoint_3 compute_circumcenter(Cell_handle cell) const
  {
    return compute_exact_circumcenter(cell);
  }

  template<typename TPoint> //exact or not
  TPoint compute_circumcenter(const TPoint &p0, const TPoint &p1,
                              const TPoint &p2, const TPoint &p3) const
  {
    return m_traits->construct_circumcenter_3_object()(p0, p1, p2, p3);
  }

  bool constrain_ray_intersection(const Point_3 &p1,
                                  const Point_3 &p2,
                                  Point_3 &res) const
  {
    return constrain_ray_intersection(p1, p2, res, p1); //taking default ref point as p1
  }

  bool constrain_ray_intersection(const Point_3 &p1, //source
                                  const Point_3 &p2, //target
                                  Point_3 &res, //intersection (if any)
                                  const Point_3 &ref) const //reference point
  {
    if(p1 == p2)
    {
#ifdef ANISO_VERBOSE
      std::cout << "CONSTRAIN_RAY_INTERSECTION : source cannot be == target ";
      std::cout.precision(20);
      std::cout << " ("<< p1 << " | " << p2 << ")" << std::endl;
#endif
      return false;
    }
    Vector_3 v(p1, p2);
    return (m_pConstrain->intersection(p1, p1+v*(m_pConstrain->get_bounding_radius()*2./std::sqrt(v*v)), ref).assign(res));
  }

  bool constrain_segment_intersection(const Point_3 &p1,
                                      const Point_3 &p2,
                                      Point_3 &res) const
  {
    return constrain_segment_intersection(p1, p2, res, p1); //taking default ref point as p1
  }

  bool constrain_segment_intersection(const Point_3 &p1, //point1
                                      const Point_3 &p2, //point2
                                      Point_3 &res, //intersection (if any)
                                      const Point_3 &ref) const //reference point
  {
    if(p1 == p2)
    {
#ifdef ANISO_VERBOSE
      std::cout << "CONSTRAIN_SEGMENT_INTERSECTION : source cannot be == target ";
      std::cout << " ("<< p1 << ")" << std::endl;
#endif
      return false;
    }

    return (m_pConstrain->intersection(p1, p2, ref).assign(res));
  }

  bool is_restricted_2_in_3(const Facet& f,
                            Point_3& surface_delaunay_ball_center,
                            const bool exact = true) const
  {
    if(Base::dimension() != 2)
      return false;
    CGAL_PROFILER("[is_restricted_2_in_3]");

    Triangle tr = m_metric.inverse_transform(Base::triangle(f));
    Vector_3 n = tr.supporting_plane().orthogonal_vector();
    if(exact)
    {
#ifdef ANISO_USE_CC_EXACT
     Point_3 fc = back_from_exact(m_metric.inverse_transform(compute_exact_circumcenter(f)));
#else
     Point_3 fc = m_metric.inverse_transform(compute_circumcenter(f));
#endif
      return (constrain_ray_intersection(fc, fc-n, surface_delaunay_ball_center)
           || constrain_ray_intersection(fc, fc+n, surface_delaunay_ball_center));
    }
    else
    {
      Point_3 fc = m_metric.inverse_transform(compute_circumcenter(f));
      return (constrain_ray_intersection(fc, fc-n, surface_delaunay_ball_center)
           || constrain_ray_intersection(fc, fc+n, surface_delaunay_ball_center));
    }
  }

  typename K::Object_3 dual(const Facet& facet) const
  {
    if(Base::dimension() == 2)
    {
      Point_3 fc = m_metric.inverse_transform(compute_circumcenter(facet));
      Triangle tr = m_metric.inverse_transform(Base::triangle(facet));
      return make_object(Line_3(fc, tr.supporting_plane().orthogonal_vector()));
    }

    Point_3 ps[3];
    get_inverse_transformed_points(ps, facet);

    Cell_handle c1 = facet.first;
    Cell_handle c2 = c1->neighbor(facet.second);
    bool f1 = !is_infinite(c1);
    bool f2 = !is_infinite(c2);

    Index offset = facet.second;
    if(!f1 && f2)
    {
      offset = c2->index(c1);
      Cell_handle tmp = c2;
      c2 = c1;    f2 = false;
      c1 = tmp;   f1 = true;
    }
    TPoint_3 offset_point = c1->vertex(offset)->point();

    Point_3 fc = m_metric.inverse_transform(compute_circumcenter(facet));
    Vector_3 t_n = this->dual_support(c1, offset).to_vector();
    Vector_3 n = m_metric.inverse_transform(t_n);
    n = std::sqrt(1./(n*n)) * n;

    if(f1)
    {
      if (f2)
      {
        Point_3 cp1 = m_metric.inverse_transform(c1->circumcenter(*(m_traits)));
        Point_3 cp2 = m_metric.inverse_transform(c2->circumcenter(*(m_traits)));
        return make_object(Segment(cp1, cp2));
      }
      else // !f2
      {
        Point_3 cp = m_metric.inverse_transform(c1->circumcenter(*(m_traits)));
        Point_3 ps3 = m_metric.inverse_transform(offset_point);

        CGAL::Orientation o1 = CGAL::orientation(ps[0], ps[1], ps[2], ps3);
        CGAL::Orientation o2 = CGAL::orientation(ps[0], ps[1], ps[2], fc + n);
        CGAL::Orientation o3 = CGAL::orientation(ps[0], ps[1], ps[2], cp);

        if(o3 == COPLANAR) //fc = cp
        {
          if(fc != cp)
          {
            std::cout << "fc != cp with cp on the facet..." << std::endl;
            std::cout << "fc : " << fc << std::endl;
            std::cout << "cp : " << cp << std::endl;
          }
          if(o1 == o2) //n points towards ps3
            n = -n;
          return make_object(Ray_3(cp, cp + n));
        }
        else if(o1 == o3)
          return make_object(Ray_3(cp, fc));
        else if(o1 != o3)
          return make_object(Ray_3(cp, Point_3(cp + Vector_3(fc, cp))));
      }
    }
    return typename K::Object_3();
  }

  bool compute_exact_dual_intersection(const Facet &facet,
                                       Point_3& p,
                                       const bool use_cache = true,
                                       const bool verbose = true,
                                       const bool super_verbose = false) const
  {
#ifdef USE_ANISO_TIMERS
    std::clock_t start = clock();
#endif
    if(!is_restricted(facet))// point is not computed here
      return false;

    if(super_verbose)
      facet_indices(facet);

    if(use_cache && is_facet_visited(facet))
    {
      p = facet.first->get_facet_surface_center(facet.second);
      return true;
    }
    bool ret_val = false;
    CGAL_PROFILER("[compute_exact_dual_intersection]");
    if(Base::dimension() == 2)
      ret_val = compute_exact_dual_intersection_2(facet, p);
    else
    {
      Point_3 ps[3];
      get_inverse_transformed_points(ps, facet);
      FT third = 1./3.;
      Point_3 facet_bar = CGAL::barycenter(ps[0], third, ps[1], third, ps[2], third);

      Cell_handle c1 = facet.first;
      Cell_handle c2 = c1->neighbor(facet.second);
      bool f1 = !is_infinite(c1);
      bool f2 = !is_infinite(c2);

      if(super_verbose)
        std::cout << "(case " << f1 << " " << f2 << ")";

      Index offset = facet.second;
      if(!f1 && f2)
      {
        offset = c2->index(c1);
        Cell_handle tmp = c2;
        c2 = c1;    f2 = false;
        c1 = tmp;   f1 = true;
      }
      TPoint_3 offset_point = c1->vertex(offset)->point();

      Exact_Point_3 fc = m_metric.inverse_transform(compute_exact_circumcenter(facet));
      Exact_Vector_3 t_n = to_exact(this->dual_support(c1, offset).to_vector());
      Exact_Vector_3 n = m_metric.inverse_transform(t_n);
      n = std::sqrt(1./(n*n)) * n;

      if(f1)
      {
        if (f2)
        {
          Exact_Point_3 cp1 = m_metric.inverse_transform(compute_exact_circumcenter(c1));
          Exact_Point_3 cp2 = m_metric.inverse_transform(compute_exact_circumcenter(c2));
          typename Traits::Compute_squared_distance_3 o = m_traits->compute_squared_distance_3_object();
          if(o(facet_bar, cp1) > o(facet_bar, cp2))
          {
            Exact_Point_3 temp = cp1;
            cp1 = cp2;
            cp2 = temp;
          }

          if(constrain_segment_intersection(back_from_exact(cp1), back_from_exact(cp2), p, back_from_exact(fc)))
            ret_val = true;
        }
        else // !f2
        {
          Exact_Point_3 cp = m_metric.inverse_transform(compute_exact_circumcenter(c1));
          Point_3 ps3 = m_metric.inverse_transform(offset_point);

          CGAL::Orientation o1 = CGAL::orientation(to_exact(ps[0]), to_exact(ps[1]),
                                                   to_exact(ps[2]), to_exact(ps3));
          CGAL::Orientation o2 = CGAL::orientation(to_exact(ps[0]), to_exact(ps[1]),
                                                   to_exact(ps[2]), fc + n);
          CGAL::Orientation o3 = CGAL::orientation(to_exact(ps[0]), to_exact(ps[1]),
                                                   to_exact(ps[2]), cp);

          if(o3 == COPLANAR) //fc = cp
          {
            if(fc != cp)
            {
              std::cout << "fc != cp with cp on the facet...[exact "<<f1<< " "<<f2 <<"]" << std::endl;
              std::cout << "fc : " << fc << std::endl;
              std::cout << "cp : " << cp << std::endl;
            }
            if(o1 == o2) //n points towards ps3
              n = -n;
            Exact_Point_3 ep = cp + n;
            if(constrain_ray_intersection(back_from_exact(cp), back_from_exact(ep), p, back_from_exact(fc)))
              ret_val = true;
          }
          else if(o1 == o3 && constrain_ray_intersection(back_from_exact(cp), back_from_exact(fc), p, back_from_exact(fc)))
            ret_val = true;
          else if(o1 != o3)
          {
            Exact_Point_3 target(cp + Exact_Vector_3(fc, cp));
            if(constrain_ray_intersection(back_from_exact(cp), back_from_exact(target), p, back_from_exact(fc)))
              ret_val = true;
          }
          if(super_verbose)
          {
            std::cout << "\t(o " << o1 << " " << o2 << " " << o3 << ")\n";
            std::cout << "\t(cp " << cp << ")\n";
            std::cout << "\t(fc " << fc << ")\n";
            std::cout << "\t(returns " << ret_val << ")\n";
          }
        }
      }

      if(ret_val)//do this in dimension 3 only (mirror_facet could crash, otherwise)
        set_facet_cache(facet, p);
      else if(verbose)
      {
        std::cout.precision(15);
        std::cout << "Oops! no intersection [exact]" << std::endl;
        std::cout << "Case:  " << f1 << " " << f2 << std::endl;
        std::cout << "Star:  " << m_center->info() << std::endl;
        std::cout << "Facet 1: " << std::endl;
        std::cout << c1->vertex((offset + 1) % 4)->info() << " ";
        std::cout << to_exact(m_metric.inverse_transform(c1->vertex((offset + 1) % 4)->point())) << std::endl;
        std::cout << c1->vertex((offset + 2) % 4)->info() << " ";
        std::cout << to_exact(m_metric.inverse_transform(c1->vertex((offset + 2) % 4)->point())) << std::endl;
        std::cout << c1->vertex((offset + 3) % 4)->info() << " ";
        std::cout << to_exact(m_metric.inverse_transform(c1->vertex((offset + 3) % 4)->point())) << std::endl;
        std::cout << c1->vertex(offset)->info() << " ";
        std::cout << to_exact(m_metric.inverse_transform(c1->vertex(offset)->point()));
        std::cout << " (.second)" << std::endl;

        std::cout << "Cell 2: " << std::endl;
        std::cout << c2->vertex(0)->info() << " ";
        std::cout << to_exact(m_metric.inverse_transform(c2->vertex(0)->point())) << std::endl;
        std::cout << c2->vertex(1)->info() << " ";
        std::cout << to_exact(m_metric.inverse_transform(c2->vertex(1)->point())) << std::endl;
        std::cout << c2->vertex(2)->info() << " ";
        std::cout << to_exact(m_metric.inverse_transform(c2->vertex(2)->point())) << std::endl;
        std::cout << c2->vertex(3)->info() << " ";
        std::cout << to_exact(m_metric.inverse_transform(c2->vertex(3)->point())) << std::endl;
      }
      if(super_verbose)
        std::cout << "Normal ("<< n <<")" << std::endl;
    }
    return ret_val;
  }

  //compute_exact_dual_intersection if triangulation is of dimension 2
  bool compute_exact_dual_intersection_2(const Facet &facet,
                                         Point_3& p) const
  {
    bool ret_val = false;
    Point_3 p1, p2;
    TPoint_3 fc = back_from_exact(m_metric.inverse_transform(compute_exact_circumcenter(facet)));
    Triangle tr = m_metric.inverse_transform(Base::triangle(facet));
    Vector_3 n = tr.supporting_plane().orthogonal_vector();

    bool b1 = false; bool b2 = false;
    if(constrain_ray_intersection(fc, fc-n, p1)) b1 = true;
    if(constrain_ray_intersection(fc, fc+n, p2)) b2 = true;
    if(b1 && b2)
    {
      typename Traits::Compute_squared_distance_3 o = m_traits->compute_squared_distance_3_object();
      if (o(fc, p1) > o(fc, p2))
        p = p2;
      else p = p1;
      ret_val = true;
    }
    else if(b1 && !b2) { p = p1 ; ret_val = true; }
    else if(!b1 && b2) { p = p2 ; ret_val = true; }
    return ret_val;
  }

  bool compute_dual_intersection(const Facet &facet,
                                 Point_3& p,
                                 const bool use_cache = true,
                                 const bool verbose = true,
                                 const bool super_verbose = false) const
  {
#ifdef USE_ANISO_TIMERS
    std::clock_t start = clock();
#endif
    if(!is_restricted(facet))// point is not computed here
      return false;

    if(use_cache && is_facet_visited(facet))
    {
      p = facet.first->get_facet_surface_center(facet.second);
      return true;
    }
    bool ret_val = false;
    typename Traits::Compute_squared_distance_3 csd = m_traits->compute_squared_distance_3_object();
    CGAL_PROFILER("[compute_dual_intersection]");

    if(Base::dimension() == 2)
    {
      Point_3 p1, p2;
      Point_3 fc = m_metric.inverse_transform(compute_circumcenter(facet));
      Triangle tr = m_metric.inverse_transform(Base::triangle(facet));
      Vector_3 n = tr.supporting_plane().orthogonal_vector();

      bool b1 = false;
      bool b2 = false;
      if(constrain_ray_intersection(fc, fc - n, p1)) b1 = true;
      if(constrain_ray_intersection(fc, fc + n, p2)) b2 = true;
      if(b1 && b2)
      {
        if (csd(fc, p1) > csd(fc, p2))
          p = p2;
        else
          p = p1;
        ret_val = true;
      }
      else if(b1 && !b2) { p = p1 ; ret_val = true; }
      else if(!b1 && b2) { p = p2 ; ret_val = true; }
    }
    else
    {
      Point_3 ps[3];
      get_inverse_transformed_points(ps, facet);
      FT third = 1./3.;
      Point_3 facet_bar = CGAL::barycenter(ps[0], third, ps[1], third, ps[2], third);

      Cell_handle c1 = facet.first;
      Cell_handle c2 = c1->neighbor(facet.second);
      bool f1 = !is_infinite(c1);
      bool f2 = !is_infinite(c2);

      if(super_verbose)
        std::cout << "(case " << f1 << " " << f2 << ")" << std::endl;

      Index offset = facet.second;
      if(!f1 && f2)
      {
        offset = c2->index(c1);
        Cell_handle tmp = c2;
        c2 = c1;    f2 = false;
        c1 = tmp;   f1 = true;
      }
      TPoint_3 offset_point = c1->vertex(offset)->point();

      Point_3 fc = m_metric.inverse_transform(compute_circumcenter(facet));
      Vector_3 t_n = this->dual_support(c1, offset).to_vector();
      Vector_3 n = m_metric.inverse_transform(t_n);
      n = std::sqrt(1./(n*n)) * n;

      if(f1)
      {
        if (f2)
        {
          Point_3 cp1 = m_metric.inverse_transform(c1->circumcenter(*(m_traits)));
          Point_3 cp2 = m_metric.inverse_transform(c2->circumcenter(*(m_traits)));
          if(csd(facet_bar, cp1) > csd(facet_bar, cp2))
          {
            Point_3 temp = cp1;
            cp1 = cp2;
            cp2 = temp;
          }

          if(cp1 == cp2)
          {
            p = cp1;
            ret_val = true;
          }
          else if(constrain_segment_intersection(cp1, cp2, p, fc))
            ret_val = true;
        }
        else // !f2
        {
          Point_3 cp = m_metric.inverse_transform(c1->circumcenter(*(m_traits)));
          Point_3 ps3 = m_metric.inverse_transform(offset_point);

          CGAL::Orientation o1 = CGAL::orientation(ps[0], ps[1], ps[2], ps3);
          CGAL::Orientation o2 = CGAL::orientation(ps[0], ps[1], ps[2], fc + n);
          CGAL::Orientation o3 = CGAL::orientation(ps[0], ps[1], ps[2], cp);

          if(o3 == COPLANAR) //fc = cp
          {
            if(fc != cp)
            {
              std::cout << "fc != cp with cp on the facet..." << std::endl;
              std::cout << "fc : " << fc << std::endl;
              std::cout << "cp : " << cp << std::endl;
            }
            if(o1 == o2) //n points towards ps3
              n = -n;
            if(constrain_ray_intersection(cp, cp + n, p, facet_bar))
              ret_val = true;
          }
          else if(o1 == o3 && constrain_ray_intersection(cp, fc, p, fc))
            ret_val = true;
          else if(o1 != o3 && constrain_ray_intersection(cp, Point_3(cp + Vector_3(fc, cp)), p, fc))
            ret_val = true;

          if(super_verbose)
          {
            std::cout << "(o " << o1 << " " << o2 << " " << o3 << ")" << std::endl;
            std::cout << "(cp " << cp << ")" << std::endl;
          }
        }
      }

      if(ret_val)//do this in dimension 3 only (mirror_facet could crash, otherwise)
        set_facet_cache(facet, p);
      else if(verbose)
      {
        std::cout.precision(15);
        std::cout << "Oops! no intersection!" << std::endl;
        std::cout << "Case:  " << f1 << " " << f2 << std::endl;
        std::cout << "Star:  " << m_center->info() << std::endl;
        std::cout << "Facet 1: " << std::endl;
        std::cout << c1->vertex((offset + 1) % 4)->info() << " ";
        std::cout << m_metric.inverse_transform(c1->vertex((offset + 1) % 4)->point()) << std::endl;
        std::cout << c1->vertex((offset + 2) % 4)->info() << " ";
        std::cout << m_metric.inverse_transform(c1->vertex((offset + 2) % 4)->point()) << std::endl;
        std::cout << c1->vertex((offset + 3) % 4)->info() << " ";
        std::cout << m_metric.inverse_transform(c1->vertex((offset + 3) % 4)->point()) << std::endl;
        std::cout << c1->vertex(offset)->info() << " ";
        std::cout << m_metric.inverse_transform(c1->vertex(offset)->point()) << " (.second)" << std::endl;
        std::cout << "Cell 2: " << std::endl;
        std::cout << c2->vertex(0)->info() << " ";
        std::cout << m_metric.inverse_transform(c2->vertex(0)->point()) << std::endl;
        std::cout << c2->vertex(1)->info() << " ";
        std::cout << m_metric.inverse_transform(c2->vertex(1)->point()) << std::endl;
        std::cout << c2->vertex(2)->info() << " ";
        std::cout << m_metric.inverse_transform(c2->vertex(2)->point()) << std::endl;
        std::cout << c2->vertex(3)->info() << " ";
        std::cout << m_metric.inverse_transform(c2->vertex(3)->point()) << std::endl;
        std::cout << "cp1/cp2" << std::endl;
        if(f1 && f2)
        {
          Point_3 cp1 = m_metric.inverse_transform(c1->circumcenter(*(m_traits)));
          Point_3 cp2 = m_metric.inverse_transform(c2->circumcenter(*(m_traits)));
          std::cout << cp1 << " " << cp2 << std::endl;
          std::cout << "is inside : " << is_inside(c1) << " " << is_inside(c2) << std::endl;
          std::cout << "check : " << m_pConstrain->side_of_constraint(cp1) << " " << m_pConstrain->side_of_constraint(cp2) << std::endl;
        }
        else
        {
          std::cout << m_metric.inverse_transform(c1->circumcenter(*(m_traits))) << std::endl;
          std::cout << m_metric.inverse_transform(offset_point) << std::endl;
        }
      }
      if(super_verbose)
        std::cout << " Normal ("<< n <<")" << std::endl;
    }
    return ret_val;
  }

  bool compute_steiner_dual_intersection(Point_3& p,
                                         const TPoint_3& tcandidate_1, // first point within the ball
                                         const TPoint_3& tcandidate_2, // second point within the ball
                                         const TPoint_3& tccf,
                                         const Facet& facet,
                                         const FT& circumradius,
                                         const int failures_count = 0) const
  {
    CGAL_PROFILER("[compute_steiner_dual_intersection]");

    typename Traits::Compute_squared_distance_3 csd =
        m_traits->compute_squared_distance_3_object();

    Point_3 candidate_1 = m_metric.inverse_transform(tcandidate_1);
    Point_3 candidate_2 = m_metric.inverse_transform(tcandidate_2);

    if(Base::dimension() == 2) //need to check if correct todo
    {
      TPoint_3 tfacetp = compute_circumcenter(facet);
      Point_3 facetp = m_metric.inverse_transform(tfacetp);
      Triangle tr = m_metric.inverse_transform(Base::triangle(facet));
      Vector_3 n = tr.supporting_plane().orthogonal_vector();
      return (constrain_ray_intersection(facetp, facetp - n, p)
             || constrain_ray_intersection(facetp, facetp + n, p));
    }
    else
    {
      if(candidate_1 == candidate_2)
        return false;

      if(constrain_segment_intersection(candidate_1, candidate_2, p))
      {
        //check if not too far away from the intersection of the dual & surface (should never happen)
        TPoint_3 tp = m_metric.transform(p);
        if(std::sqrt(csd(tccf, tp)) > circumradius)
        {
          std::cout << "intersection point outside of picking ball in compute_steiner_dual_intersection" << std::endl;
#if 1//def ANISO_DEBUG_STEINER_DUAL
          std::cout << "that should not happen !" << std::endl;
          std::cout << "candidates : " << std::endl;
          std::cout << tcandidate_1 << std::endl << tcandidate_2 << std::endl;
          std::cout << "tccf : " << tccf << std::endl;
          std::cout << "tp : " << tp << std::endl;
          std::cout << "checking dist tccf-tcandi1 : " << std::sqrt(csd(tccf, tcandidate_1)) << std::endl;
          std::cout << "checking dist tccf-tcandi2 : " << std::sqrt(csd(tccf, tcandidate_2)) << std::endl;
          std::cout << "dist tccf-tp : " << std::sqrt(csd(tccf, tp)) << " and circumradius : " << circumradius << std::endl;
#endif
          return false;
        }
        else
          return true;
      }
      else
      {
        if(failures_count > 100)
        {
          Point_3 ccf = m_metric.inverse_transform(tccf);
#if 1//def ANISO_DEBUG_STEINER_DUAL
          std::cout << "failures : " << failures_count << std::endl;
          double eval1 = m_pConstrain->side_of_constraint(candidate_1);
          double eval2 = m_pConstrain->side_of_constraint(candidate_2);
          double evalc = m_pConstrain->side_of_constraint(ccf);
          std::cout << "no intersection..." << std::endl;
          std::cout << "tcandidate_1 : " << tcandidate_1 << std::endl;
          std::cout << "tcandidate_2 : " << tcandidate_2 << std::endl;
          std::cout << "candidate_1 : " << candidate_1 << " side of constraint : " << eval1 << std::endl;
          std::cout << "candidate_2 : " << candidate_2 << " side of constraint : " << eval2 << std::endl;
          std::cout << "checking dist 1 : " << std::sqrt(csd(tccf, tcandidate_1)) << std::endl;
          std::cout << "checking dist 2 : " << std::sqrt(csd(tccf, tcandidate_2)) << std::endl;
          std::cout << "tccf : " << tccf << " side of constraint for ccf : " << evalc << std::endl;
          std::cout << "circum : " << circumradius << std::endl;
#endif
          std::cout << "had to use p = ccf in compute_steiner_dual (not good)." << std::endl;
          p = ccf;
          return true;
        }
        return false;
      }
    }
  }

// ENCROACHMENT
  bool is_facet_encroached(const TPoint_3 &tp, const Facet &facet) const
  {
    Point_3 c;
#ifdef ANISO_USE_EXACT
    compute_exact_dual_intersection(facet, c);
#else
    compute_dual_intersection(facet, c);
#endif
    TPoint_3 tc = m_metric.transform(c);
    TPoint_3 tq = facet.first->vertex((facet.second + 1) % 4)->point();
    typename Traits::Compute_squared_distance_3 o = m_traits->compute_squared_distance_3_object();
    return o(tc, tp) < o(tc, tq);
  }

  bool is_encroached(const Point_3 &p, Facet &facet) const
  {
    TPoint_3 tp = m_metric.transform(p);

    Facet_set_iterator fi = restricted_facets_begin();
    Facet_set_iterator fend = restricted_facets_end();
    for(; fi != fend; fi++)
    {
      if(is_facet_encroached(tp, *fi))
      {
        facet = *fi;
        return true;
      }
    }
    return false;
  }

public:
  std::size_t clean(bool verbose = false) //remove non-adjacent vertices
  {
    typedef std::pair<TPoint_3, int>          PPoint;
    std::vector<PPoint> backup;
    std::size_t nbv = this->number_of_vertices();

    Vertex_handle_handle vit = finite_adjacent_vertices_begin();
    Vertex_handle_handle vend = finite_adjacent_vertices_end();
    if(static_cast<std::size_t>(vend-vit) == nbv-1) // -1 due to center
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
      insert(backup[i].first)->info() = backup[i].second;

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
  Vertex_handle insert(const TPoint_3 &tp, const Cell_handle &cell)
  {
    Vertex_handle retval = Base::insert(tp, cell);
    invalidate_cache();
    return retval;
  }

  // warning : point should be transformed before this insert function is called
  Vertex_handle insert(const TPoint_3 &tp)
  {
    Vertex_handle retval = Base::insert(tp);
    invalidate_cache();
    return retval;
  }

  // warning : point should be transformed before this insert function is called
  Vertex_handle insert_to_star_(const TPoint_3 &tp,
                                const bool conditional)
  {
    if(Base::dimension() < 3 || !conditional)
      return this->insert(tp);

    Cell_handle cell;
    if(is_conflicted(tp, cell))
      return this->insert(tp, cell);

    else // point is not in conflict with cells incident to center
      return Vertex_handle();
  }

public:
  Vertex_handle insert_to_star(const Point_3& p,
                               const int id,
                               const bool conditional)
  {
    Vertex_handle v = insert_to_star_(m_metric.transform(p), conditional);
    if(v != center() && v != Vertex_handle())
      v->info() = id;
    return v;
  }

  Vertex_handle insert_to_star_hole(const Point_3& point,
                                    const Index id,
                                    std::vector<Cell_handle>& cells,
                                    Facet facet)
  {
    Vertex_handle v = Base::insert_in_hole(m_metric.transform(point),
                                           cells.begin(),
                                           cells.end(),
                                           facet.first,
                                           facet.second);
    invalidate_cache();
    if(v != center() && v != Vertex_handle())
      v->info() = id;
    else
      std::cout << "problem in insert_to_star_hole..." << std::endl;
    return v;
  }

  //could return "cell" and avoid a call of "is_conflicted()" in find_conflict() by giving the
  //correct cell hint directly TODO
  int simulate_insert_to_star(const Point_3& p, const int id)
  {
    CGAL_HISTOGRAM_PROFILER("V", this->number_of_vertices());
    TPoint_3 tp = m_metric.transform(p);
    bool found_vertex = false;
    Vertex_handle vh;

#ifndef ANISO_BRUTE_FORCE_SIMULATE_INSERT_TO_STAR
    int li, lj;
    typename Base::Locate_type lt;
    Cell_handle c = Base::locate(tp, lt, li, lj, m_center->cell());
    if(lt ==  Base::VERTEX)
    {
      found_vertex = true;
      vh = c->vertex(li);
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
      Cell_handle cell;
      if(is_conflicted(tp, cell))
        return id;
    }
    return -1; // no conflict
  }

//Shaking
  void shake_center(FT radius)
  {
    std::cout << "shaking " << index_in_star_set() << "'s center" << std::endl;
    typename Traits::Compute_random_point_3 random =
        traits()->compute_random_point_3_object();

    Point_3 new_center;

    // if surface star, get 2 random pts in a sphere and find the intersection with the surface
    if(is_surface_star())
    {
      bool found_new_point = false;
      while(!found_new_point)
      {
        TPoint_3 tcandidate_1 = random(center()->point(), radius);
        TPoint_3 tcandidate_2 = random(center()->point(), radius);
        Point_3 candidate_1 = m_metric.inverse_transform(tcandidate_1);
        Point_3 candidate_2 = m_metric.inverse_transform(tcandidate_2);

        found_new_point = constrain_segment_intersection(candidate_1, candidate_2, new_center);
      }
    }
    else
      new_center = m_metric.inverse_transform(random(center()->point(), radius));

    //can't reset the metric yet, even if we wanted to
    reset(new_center, this->index_in_star_set(), m_metric, is_surface_star());
  }

  FT shake_center()
  {
    // find closest point in the first ring
    FT min_sq_dist = 1e30;

    Vertex_handle_handle favit = finite_adjacent_vertices_begin();
    Vertex_handle_handle favend = finite_adjacent_vertices_end();
    for(; favit!=favend; ++favit)
    {
      Vertex_handle vit = *favit;
      FT sq_dist = m_traits->compute_squared_distance_3_object()(center()->point(), vit->point());
      if(sq_dist < min_sq_dist)
        min_sq_dist = sq_dist;
    }

    FT radius = CGAL::sqrt(min_sq_dist)/3.;
    shake_center(radius);
    return radius;
  }

// INFO ABOUT THE VERTICES / FACETS / CELLS
public:
  void print_vertices(const bool print_points = false) const
  {
    std::cout << "\t Star_" << index_in_star_set();
    std::cout << " vertices ("<<(this->number_of_vertices())<<" v.)\t (";
    typename std::set<Vertex_handle> vertices;

    if(restricted_facets_begin() == restricted_facets_end())
      std::cout << "empty";
    else
    {
      for(Facet_set_iterator fit = restricted_facets_begin();
          fit != restricted_facets_end();
          fit++)
      {
        Facet f = *fit;
        vertices.insert(f.first->vertex((f.second + 1) % 4));
        vertices.insert(f.first->vertex((f.second + 2) % 4));
        vertices.insert(f.first->vertex((f.second + 3) % 4));
      }
      for(typename std::set<Vertex_handle>::iterator vit = vertices.begin();
          vit != vertices.end(); vit++)
      {
        if(print_points) std::cout << "\n\t";
        std::cout << (*vit)->info() << " ";
        if(print_points) std::cout << " : " << (*vit)->point();
      }
    }

    std::cout << ")\t(";
    typename Base::Finite_vertices_iterator vit = this->finite_vertices_begin();
    typename Base::Finite_vertices_iterator vend = this->finite_vertices_end();
    for(; vit != vend; vit++)
      std::cout << vit->info() << " ";
    std::cout << ")\n";
  }

  std::size_t count_restricted_facets() const
  {
    update_star_caches();
    return restricted_facets_cache.size();
  }

  void print_restricted_facets() const
  {
    std::cout << "Restricted facets of " << m_center->info() << std::endl;
    Facet_set_iterator fi = restricted_facets_begin();
    Facet_set_iterator fend = restricted_facets_end();
    for(; fi != fend; fi++)
    {
      std::cout << fi->first->vertex((fi->second+1)%4)->info() << " || ";
      std::cout << fi->first->vertex((fi->second+1)%4)->point() << std::endl;
      std::cout << fi->first->vertex((fi->second+2)%4)->info() << " || ";
      std::cout << fi->first->vertex((fi->second+2)%4)->point() << std::endl;
      std::cout << fi->first->vertex((fi->second+3)%4)->info() << " || ";
      std::cout << fi->first->vertex((fi->second+3)%4)->point() << std::endl;
      if(Base::dimension() == 3)
      {
        std::cout << "second: " << fi->first->vertex(fi->second)->info() << " ";
        std::cout << fi->first->vertex(fi->second)->point() << std::endl;
      }
    }
  }

  void print_facets() const
  {
    std::cout << "\t Star_" << index_in_star_set() << " faces :\n";
    typename Base::Finite_facets_iterator fit = this->finite_facets_begin();
    typename Base::Finite_facets_iterator fend = this->finite_facets_end();
    for(; fit != fend; ++fit)
    {
      Facet f = *fit;
      std::cout << "\t\t f(" << f.first->vertex((f.second+1)%4)->info()
                    << " " << f.first->vertex((f.second+2)%4)->info()
                    << " " << f.first->vertex((f.second+3)%4)->info()
                    << " - " << f.first->vertex(f.second)->info() << ")";
      if(is_restricted(f))
        std::cout << "restricted";
      std::cout << ",\n";
    }
  }

  void print_cells() const
  {
    std::cout << "list cells of " << m_center->info() << std::endl;
    typename Base::Finite_cells_iterator cit = this->finite_cells_begin();
    typename Base::Finite_cells_iterator cend = this->finite_cells_end();
    for(;cit!=cend;++cit)
    {
      std::cout << "Cell ";
      std::cout << cit->vertex(0)->info() << " ";
      std::cout << cit->vertex(1)->info() << " ";
      std::cout << cit->vertex(2)->info() << " ";
      std::cout << cit->vertex(3)->info() << " " << Base::is_valid(cit, true) << std::endl;
    }
    std::cout << "end" << std::endl;
  }

  void facet_indices(const Facet& f) const
  {
    std::cout << "Facet[";
    std::vector<int> indices;
    indices.push_back(f.first->vertex((f.second+1)%4)->info());
    indices.push_back(f.first->vertex((f.second+2)%4)->info());
    indices.push_back(f.first->vertex((f.second+3)%4)->info());
    std::sort(indices.begin(), indices.end());
    std::cout << indices[0] << " " << indices[1] << " " << indices[2];
    std::cout << "]";
  }

  void cell_indices(const Cell_handle& c) const
  {
    std::cout << "Cell[";
    std::vector<int> indices;
    indices.push_back(c->vertex(0)->info());
    indices.push_back(c->vertex(1)->info());
    indices.push_back(c->vertex(2)->info());
    indices.push_back(c->vertex(3)->info());
    std::sort(indices.begin(), indices.end());
    std::cout << indices[0] << " " << indices[1] << " " << indices[2] << " " << indices[3];
    std::cout << "]";
  }

// VISUALIZATION WITH OPENGL
public:
  void gl_draw_center() const
  {
    if(this->is_surface_star())
    {
      ::glPointSize(5.);
      ::glColor3f(0., 0., 0.);
    }
    else
    {
      ::glPointSize(3.);
      ::glColor3f(0.5f, 0.5f, 0.5f);
    }
    ::glBegin(GL_POINTS);
    ::glVertex3d(m_center_point.x(), m_center_point.y(), m_center_point.z());
    ::glEnd();
  }

  void gl_draw(const typename K::Plane_3& plane,
               const bool draw_edges = true) const
  {
    if(! this->is_surface_star())
      return;
    gl_draw_center();
    Facet_set_iterator fit = restricted_facets_begin();
    Facet_set_iterator fend = restricted_facets_end();
    for(; fit != fend; fit++)
    {
      const Cell_handle& cell = (*fit).first;
      const int& i = (*fit).second;

      const Point_3& pa = m_metric.inverse_transform(cell->vertex((i+1)&3)->point());
      const Point_3& pb = m_metric.inverse_transform(cell->vertex((i+2)&3)->point());
      const Point_3& pc = m_metric.inverse_transform(cell->vertex((i+3)&3)->point());
      if(is_above_plane<K>(plane, pa, pb, pc))
        gl_draw_triangle<K>(pa, pb, pc, EDGES_AND_FACES, 205., 175., 149);
    }
  }

  void gl_draw_cell(const typename K::Plane_3& plane) const
  {
    gl_draw_center();

    Cell_handle_handle ci = finite_star_cells_begin();
    Cell_handle_handle cend = finite_star_cells_end();
    for (; ci != cend; ci++)
    {
      std::vector<Point_3> cell_points(4);
      cell_points[0] = m_metric.inverse_transform((*ci)->vertex(0)->point());
      cell_points[1] = m_metric.inverse_transform((*ci)->vertex(1)->point());
      cell_points[2] = m_metric.inverse_transform((*ci)->vertex(2)->point());
      cell_points[3] = m_metric.inverse_transform((*ci)->vertex(3)->point());
      if(is_above_plane<K>(plane, cell_points[0], cell_points[1], cell_points[2], cell_points[3]))
      {
        for(int i = 0; i < 4; i++)
        {
          gl_draw_triangle<K>(cell_points[i%4],
                              cell_points[(i+1)%4],
                              cell_points[(i+2)%4],EDGES_ONLY, 255, 117, 44);
        }
      }
    }
  }

  void gl_draw_metric(const typename K::Plane_3& plane,
                      double coeff) const
  {
    if(!this->is_surface_star())
      return;
    if(!is_above_plane<K>(plane, this->center_point()))
      return;

    Point_3 p = this->center_point();
    Vector_3 vec = CGAL::NULL_VECTOR;
    double val = 0.;

    val = m_metric.get_min_eigenvalue();
    m_metric.get_min_eigenvector(vec);
    ::glColor3f(0.,0.,1.f);
    ::gl_draw_arrow<K>(p, p + coeff*vec/val);

    val = m_metric.get_max_eigenvalue();
    m_metric.get_max_eigenvector(vec);
    ::glColor3f(1.f,0.,0.);
    ::gl_draw_arrow<K>(p, p + coeff*vec/val);

    val = m_metric.get_third_eigenvalue();
    m_metric.get_third_eigenvector(vec);
    ::glColor3f(0.,1.f,0.);
    ::gl_draw_arrow<K>(p, p + coeff*vec/val);
  }

  //void gl_draw_all(const typename K::Plane_3& plane) const
  //{
  //  gl_draw_center();
  //  if(! is_surface_star())
  //    return;
  //  typename Base::Finite_facets_iterator fit2 = this->finite_facets_begin();
  //  typename Base::Finite_facets_iterator fend2 = this->finite_facets_end();
  //  for(; fit2 != fend2; fit2++)
  //  {
  //    Facet f = *fit2;
  //    const Point_3& pa = m_metric.inverse_transform(f.first->vertex((f.second+1)&3)->point());
  //    const Point_3& pb = m_metric.inverse_transform(f.first->vertex((f.second+2)&3)->point());
  //    const Point_3& pc = m_metric.inverse_transform(f.first->vertex((f.second+3)&3)->point());
  //    if(is_above_plane<K>(plane, pa, pb, pc))
  //    {
  //      if(is_in_star(f) && is_restricted(f))
  //        gl_draw_triangle<Kernel>(pa, pb, pc, EDGES_AND_FACES, 205., 175., 149);
  //      else
  //        gl_draw_triangle<Kernel>(pa, pb, pc, EDGES_ONLY/*, 139., 137., 137.*/);
  //    }
  //  }
  //}

  void gl_draw_dual(const typename K::Plane_3& plane) const
  {
    GLboolean was = (::glIsEnabled(GL_LIGHTING));
    if(was)
      ::glDisable(GL_LIGHTING);

    typename Base::Finite_facets_iterator fit = this->finite_facets_begin();
    typename Base::Finite_facets_iterator fend = this->finite_facets_end();
    for(; fit != fend; fit++)
    {
      Facet f = *fit;
      if(!is_in_star(f))
        continue;

      const Point_3& pa = m_metric.inverse_transform(f.first->vertex((f.second+1)&3)->point());
      const Point_3& pb = m_metric.inverse_transform(f.first->vertex((f.second+2)&3)->point());
      const Point_3& pc = m_metric.inverse_transform(f.first->vertex((f.second+3)&3)->point());
      if(!is_above_plane<K>(plane, pa, pb, pc))
        continue;

      if(is_restricted(f)) { ::glLineWidth(2.);}
      else                 { ::glColor3f(0.,0.,0.);  ::glLineWidth(1.);}

      typename K::Object_3 o = dual(f);
      typename K::Line_3 l;
      typename K::Ray_3 r;
      typename K::Segment_3 s;
      if(CGAL::assign(s, o))
      {
        ::glColor3f(0.,0.,1.f);
        gl_draw_segment<K>(s.source(), s.target(), true /*end points*/);
      }
      else if(CGAL::assign(r, o))
      {
        ::glColor3f(0.,1.f, 0.);
        gl_draw_segment<K>(r.source(), r.source() + r.to_vector()
            * (m_pConstrain->get_bounding_radius() * 10.0 / std::sqrt(r.to_vector() * r.to_vector())), true);
      }
      else if(CGAL::assign(l, o))
        std::cout << "gl_draw_dual : line dual\n";

      //experimental stuff
      Point_3 p;
      if(compute_dual_intersection(f, p))
      {
        ::glColor3f(1.f,0.,0.);
        ::glPointSize(5.);
        ::glBegin(GL_POINTS);
        ::glVertex3d(p.x(),p.y(),p.z());
        ::glEnd();

        //Point_3 centroid = CGAL::centroid(pa, pb, pc);
        //gl_draw_segment<K>(p, centroid);
      }
    }
    if(was)
      ::glEnable(GL_LIGHTING);
  }

  void gl_draw_surface_delaunay_balls(const typename K::Plane_3& plane) const
  {
    if(!is_above_plane<K>(plane, this->center_point()))
      return;

    Facet_set_iterator fit = restricted_facets_begin();
    Facet_set_iterator fend = restricted_facets_end();
    for(; fit != fend; fit++)
    {
      Facet f = *fit;
      Point_3 ce;
      if(!is_above_plane<K>(plane, m_metric.inverse_transform(f.first->vertex((f.second+1) % 4)->point()),
                                m_metric.inverse_transform(f.first->vertex((f.second+2) % 4)->point()),
                                m_metric.inverse_transform(f.first->vertex((f.second+3) % 4)->point())) )
        continue;
#ifdef ANISO_USE_EXACT
      this->compute_exact_dual_intersection(f, ce);
#else
      this->compute_dual_intersection(f, ce);
#endif
      TPoint_3 tce = m_metric.transform(ce);
      TPoint_3 tp2 = f.first->vertex((f.second+1) % 4)->point();
      FT sqr = m_traits->compute_squared_distance_3_object()(tce, tp2);

        //ellipsoid's a, b & c
      FT a = std::sqrt(sqr)/m_metric.get_max_eigenvalue();
      FT b = std::sqrt(sqr)/m_metric.get_min_eigenvalue();
      FT c = std::sqrt(sqr)/m_metric.get_third_eigenvalue();

        //rotation to align on the eigenvectors
      Vector_3 v1,v2,vn;
      m_metric.get_max_eigenvector(v1);
      m_metric.get_min_eigenvector(v2);
      m_metric.get_third_eigenvector(vn);

      gl_draw_ellipsoid<K>(CGAL::ORIGIN, ce, 10, 10, a, b, c, v1, v2, vn, 232, 67, 183);
    }
  }

// DEBUG
public:
  bool debug_steiner_point(const Point_3& steiner_point,
                           const Facet& f,
                           const bool more = false) const
  {
    //std::cout << "<-debug_steiner_point\n";
    bool bug = false;
    Point_3 p;
    compute_exact_dual_intersection(f,p, false, true, more);

    const TPoint_3& pf = f.first->vertex((f.second+1)&3)->point();
    const TPoint_3& center = m_metric.transform(p);
    const TPoint_3& steiner = m_metric.transform(steiner_point);

    Sphere s(center, m_traits->compute_squared_distance_3_object()(center, pf));
    if(s.has_on_unbounded_side(steiner))
    {
      facet_indices(f);
      std::cout << "\n\tSteiner point ("
        <<steiner_point<< ") outside exact surface Delaunay ball\n";
      std::cout << "\tp is " << p << std::endl;
      bug = true;
    }

    //Point_3 p2;
    //compute_dual_intersection(f,p2, false,true,more);
    //const TPoint_3& center2 = m_metric.transform(p);
    //Sphere s2(center2, m_traits->compute_squared_distance_3_object()(center2, pf));
    //if(s2.has_on_unbounded_side(steiner))
    //{
    //  facet_indices(f);
    //  std::cout << "\n\tSteiner point ("
    //    <<steiner_point<< ") outside inexact surface Delaunay ball\n";
    //  bug = true;
    //}
    //std::cout << "->\n";
    return bug;
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

  void reset(const Point_3& centerpoint,
             const int index,
             const Metric& metric_,
             const bool is_surface_star = true)
  {
    this->clear();
    m_center_point = centerpoint;
    m_metric = metric_;
    m_center = Base::insert(m_metric.transform(centerpoint));
    m_center->info() = index;
    m_is_surface_star = is_surface_star;
    this->infinite_vertex()->info() = index_of_infinite_vertex;

    invalidate_cache();
  }

  void reset(const Point_3& new_centerpoint)
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
  Stretched_Delaunay_3(const Criteria* criteria_,
                       const Constrain_surface* pconstrain_surface,
                       const bool is_surface_star = true,
                       const bool is_in_3D_mesh = false)
    :
      Base(*(m_traits = new Traits())),
      m_center_point(CGAL::ORIGIN),
      m_is_surface_star(is_surface_star),
      m_metric(),
      m_pConstrain(pconstrain_surface),
      m_criteria(new Stretched_criteria(*m_traits, criteria_)),
      m_is_in_3D_mesh(is_in_3D_mesh),
      m_is_topological_disk(false),
      m_is_valid_topo_disk(false),
      m_is_valid_bbox(false),
      m_bbox_needs_aabb_update(false),
      m_metric_needs_update(false),
      is_cache_dirty(true),
      restricted_facets_cache(),
      incident_cells_cache(),
      finite_incident_cells_cache(),
      finite_adjacent_vertices_cache(),
      finite_adjacent_restricted_vertices_cache(),
      m_active(true)
  {
    m_center = Vertex_handle();
    m_bbox = m_pConstrain->get_bbox(); // in M_euclidean
    this->infinite_vertex()->info() = index_of_infinite_vertex;
  }

  Stretched_Delaunay_3(const Criteria* criteria_,
                       const Point_3 &centerpoint,
                       const int &index,
                       const Metric &metric_,
                       const Constrain_surface* pconstrain_surface,
                       const bool is_surface_star = true,
                       const bool is_in_3D_mesh = false)
    :
      Base(*(m_traits = new Traits())),
      m_center_point(centerpoint),
      m_is_surface_star(is_surface_star),
      m_metric(metric_),
      m_pConstrain(pconstrain_surface),
      m_criteria(new Stretched_criteria(*m_traits, criteria_)),
      m_is_in_3D_mesh(is_in_3D_mesh),
      m_is_topological_disk(false),
      m_is_valid_topo_disk(false),
      m_is_valid_bbox(false),
      m_bbox_needs_aabb_update(false),
      m_metric_needs_update(false),
      is_cache_dirty(true),
      restricted_facets_cache(),
      incident_cells_cache(),
      finite_incident_cells_cache(),
      finite_adjacent_vertices_cache(),
      finite_adjacent_restricted_vertices_cache(),
      m_active(true)
  {
    m_center = Base::insert(m_metric.transform(centerpoint));
    m_center->info() = index;
    m_bbox = m_pConstrain->get_bbox(); // in M_euclidean
    this->infinite_vertex()->info() = index_of_infinite_vertex;
  }

  ~Stretched_Delaunay_3()
  {
    delete m_traits;
    delete m_criteria;
  }
};

} // Anisotropic_mesh_3
} // CGAL

#endif //CGAL_ANISOTROPIC_MESH_3_STRETCHED_DELAUNAY_3

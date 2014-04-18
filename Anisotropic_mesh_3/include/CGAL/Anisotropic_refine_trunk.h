#ifndef CGAL_ANISOTROPIC_MESH_3_REFINE_TRUNK_H
#define CGAL_ANISOTROPIC_MESH_3_REFINE_TRUNK_H

#include <CGAL/Constrain_surface_3.h>
#include <CGAL/Criteria.h>
#include <CGAL/Metric_field.h>
#include <CGAL/Stretched_Delaunay_3.h>
#include <CGAL/Star_consistency.h>

#include <CGAL/aabb_tree/aabb_tree_bbox.h>
#include <CGAL/aabb_tree/aabb_tree_bbox_primitive.h>

#include <CGAL/kd_tree/Kd_tree_for_star_set.h>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

//these classes are used to keep in memory the various conflict zones for each modified star.
//This is done for multiple reasons :
// * Avoid computing it multiple times (encroachment tests, pick_valid tests, point insertion...)
// * Keep in memory the facets/cells that will be destroyed to be able to track elements that might
//   need to be inserted in the pqueues (especially those not in modified stars)
template<typename K>
class Conflict_zone
{
  typedef Conflict_zone<K>                        Self;

public:
  typedef Stretched_Delaunay_3<K>                 Star;
  typedef typename Star::Facet                    Facet;
  typedef typename Star::Facet_vector             Facet_vector;
  typedef typename Star::Facet_handle             Facet_handle;
  typedef typename Star::Cell_handle              Cell_handle;
  typedef typename Star::Cell_handle_vector       Cell_handle_vector;
  typedef typename Star::Cell_handle_handle       Cell_handle_handle;

private:
  //Conflict zones defining vectors. Valid only until the insertion of the point!
  Facet_vector m_boundary_facets;
  Cell_handle_vector m_cells;
  Facet_vector m_internal_facets;

public:
  Facet_vector& boundary_facets() { return m_boundary_facets; }
  const Facet_vector& boundary_facets() const { return m_boundary_facets; }
  Cell_handle_vector& cells() { return m_cells; }
  const Cell_handle_vector& cells() const { return m_cells; }
  Facet_vector& internal_facets() { return m_internal_facets; }
  const Facet_vector& internal_facets() const { return m_internal_facets; }

  bool empty()
  {
    return m_boundary_facets.empty() && m_cells.empty() && m_internal_facets.empty();
  }

  //Constructors
  Conflict_zone()
    :
      m_boundary_facets(),
      m_cells(),
      m_internal_facets()
  { }
};


//TODO: something better than initialization of a point at (1e17, 1e17, 1e17)
//      CGAL::ORIGIN isn't satisfying since the insertion point could be the same.
//      Maybe use NaN?
template<typename K, typename KExact = K>
class Stars_conflict_zones
{
  typedef Stars_conflict_zones<K, KExact>                  Self;

public:
  typedef Stretched_Delaunay_3<K, KExact>                  Star;
  typedef Star*                                            Star_handle;
  typedef typename std::vector<Star_handle>                Star_vector;
  typedef CGAL::Anisotropic_mesh_3::Conflict_zone<K>       Czone;
  typedef typename Star::Index                             Index;
  typedef typename Star::Point_3                           Point;
  typedef typename KExact::Point_3                         Exact_Point;
  typedef typename Star::Vertex                            Vertex;
  typedef typename Star::Facet                             Facet;
  typedef typename Star::Facet_handle                      Facet_handle;
  typedef typename Star::Cell                              Cell;
  typedef typename Star::Cell_handle                       Cell_handle;
  typedef typename Star::Cell_handle_handle                Cell_handle_handle;
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
  const Star_vector& m_stars;
  std::map<Index, Czone> m_conflict_zones;
  Point m_conflict_p; // point causing the conflict (only needed to debug / assert)
  Index m_conflict_p_id; //id of the point

  bool m_cells_need_check; // no need to compute cells_to_check if we haven't reached cell level
  bool m_are_checks_computed; //already computed or not

  //we need to keep track of the "soon to be" destroyed facets/cells SEEN FROM UNMODIFIED STARS
  //not lose consistency in the unmodified stars and add them to the queue if needed
  std::map<Index, std::vector<Facet> > m_facets_to_check;
  std::map<Index, std::vector<Cell_handle> > m_cells_to_check;

public:
  const std::map<Index, Czone>& conflict_zones() const { return m_conflict_zones; }
  Czone& conflict_zone(Index i) { return m_conflict_zones[i]; }
  const Czone& conflict_zone(Index i) const { return m_conflict_zones[i]; }
  Point& conflicting_point() { return m_conflict_p; }
  const Point& conflicting_point() const { return m_conflict_p; }
  Index& conflicting_point_id() { return m_conflict_p_id; }
  const Index& conflicting_point_id() const { return m_conflict_p_id; }
  std::map<Index, std::vector<Facet> >& facets_to_check() { return m_facets_to_check; }
  std::map<Index, std::vector<Cell_handle> >& cells_to_check() { return m_cells_to_check; }
  bool& cells_need_checks() { return m_cells_need_check; }
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
    m_conflict_p = Point(1e17, 1e17, 1e17);
    m_conflict_p_id = -1;
    m_facets_to_check.clear();
    m_cells_to_check.clear();
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
    return m_facets_to_check.empty() && m_cells_to_check.empty();
  }

  bool is_boundary_facet_restricted_after_insertion(Star_handle star,
                                                    Facet f)
  {
    //this function predicts whether a boundary facet will still be restricted
    //after the insertion of the steiner point.

    Cell_handle c1 = f.first; // is in conflict
    Cell_handle c2 = c1->neighbor(f.second);

#ifdef ANISO_USE_INSIDE_EXACT
    bool inside2 = star->is_inside_exact(c2);
    Exact_Point tncc = star->compute_circumcenter(to_exact(c1->vertex((f.second+1)%4)->point()),
                                                  to_exact(c1->vertex((f.second+2)%4)->point()),
                                                  to_exact(c1->vertex((f.second+3)%4)->point()),
                                                  to_exact(star->metric().transform(m_conflict_p)));
    Exact_Point ncc = star->metric().inverse_transform(tncc);
    bool inside_new = (star->constrain()->side_of_constraint(back_from_exact(ncc) == ON_POSITIVE_SIDE);
#else
    bool inside2 = star->is_inside(c2);
    Point tncc = star->compute_circumcenter(c1->vertex((f.second+1)%4)->point(),
                                            c1->vertex((f.second+2)%4)->point(),
                                            c1->vertex((f.second+3)%4)->point(),
                                            star->metric().transform(m_conflict_p));
    Point ncc = star->metric().inverse_transform(tncc);
    bool inside_new = (star->constrain()->side_of_constraint(ncc) == ON_POSITIVE_SIDE);
#endif
    return (inside2 && !inside_new) || (!inside2 && inside_new);
  }

  struct Facet_compare
  {
    Star_handle star;
    Facet& f;
    Facet_ijk fijk;
    bool mirror;

    Facet_compare(Star_handle star_, Facet& f_, bool mirror_)
      :
        star(star_), f(f_), fijk(f_), mirror(mirror_)
    { }

    bool operator()(Facet& ff)
    {
      if(fijk == Facet_ijk(ff))
      {
        //if mirror is true (case of a boundary facet), we want f(c,i) to have
        //a 'c' that will not be invalidated by the point insertion. Since we
        //know that the boundary facets are such that f.first is in conflict, we
        //mirror the facet such that their .second is different
        if(mirror &&
           f.first->vertex(f.second)->info() == ff.first->vertex(ff.second)->info())
        {
          f = star->mirror_facet(f); // no need to change fijk
        }
        return true;
      }
      return false;
    }
  };

  void check_from_facet(Facet_handle fit,
                        const std::vector<Index>& indices)
  {
    //fit is an internal facet or a boundary facet that won't be restricted after insertion
    for(std::size_t j=0; j<indices.size(); ++j)
    {
      Star_handle sj = m_stars[indices[j]];

      Facet f_in_sj; //the facet f seen from sj
      if(!sj->has_facet(Facet_ijk(*fit), f_in_sj)) //only searches through restricted facets
        continue;

      // At this point, we don't know for f(c,i) if c is a conflict cell

      iterator czcit = m_conflict_zones.find(indices[j]);
      if(czcit == m_conflict_zones.end())
      {
        m_facets_to_check[indices[j]].push_back(f_in_sj);
        continue;
      }

      //TODO: might be easier to sort all the facets once, since the predicate
      //here sorts a facet (creates a Facet_ijk) every time it compares
      Czone czj = czcit->second;
      Facet_handle fhit = std::find_if(czj.internal_facets().begin(),
                                       czj.internal_facets().end(),
                                       Facet_compare(sj, f_in_sj, false/*no mirroring*/));
      //if si's f is found in sj's internals => nothing to do
      if(fhit != czj.internal_facets().end())
        continue;

      Facet_handle fhitb = std::find_if(czj.boundary_facets().begin(),
                                        czj.boundary_facets().end(),
                                        Facet_compare(sj, f_in_sj, true/*mirroring*/));

      //if si's f is found in sj's boundaries,
      //    if sj's is still restricted after insertion => add it
      //    if sj's is not restricted anymore after insertion => nothing to do
      if(fhitb != czj.boundary_facets().end())
      {
         if(is_boundary_facet_restricted_after_insertion(sj, *fhitb))
         {
           m_facets_to_check[indices[j]].push_back(f_in_sj);
         }
         continue;
      }

      //if si's f is not found in sj's internals AND is not found in
      //sj's boundaries => add it
      m_facets_to_check[indices[j]].push_back(f_in_sj);
    }
  }

  void check_from_internal_facets(iterator mit,
                                  const std::map<Facet_ijk, int>& internal_facets_counter)
  {
    Index sid = mit->first;
    Star_handle si = m_stars[sid];
    Czone& czi = mit->second;

    Facet_handle fit = czi.internal_facets().begin();
    Facet_handle fend = czi.internal_facets().end();
    for(; fit!=fend; ++fit)
    {
      if(!si->is_restricted(*fit))
        continue;

      typename std::map<Facet_ijk, int>::const_iterator fcit = internal_facets_counter.find(Facet_ijk(*fit));
      if(fcit != internal_facets_counter.end() && fcit->second == 3)
        continue;

      std::vector<Index> indices; //the other stars making up the facet
      for(int i=0; i<3; ++i)
      {
        Index id = fit->first->vertex((fit->second+i+1)%4)->info();
        if(id != sid)
          indices.push_back(id);
      }

      //Ignore facets that are not incident to the star's center point (they do
      //not exist in the queues). They will be cleaned eventually
      if(indices.size() == 3)
        continue;

      check_from_facet(fit, indices);
    }
  }

  void check_from_boundary_facets(iterator mit)
  {
    //boundary facets could stop being restricted and create inconsistencies
    Index sid = mit->first;
    Star_handle si = m_stars[sid];
    Czone& czi = mit->second;

    Facet_handle fit = czi.boundary_facets().begin();
    Facet_handle fend = czi.boundary_facets().end();
    for(; fit!=fend; ++fit)
    {
      if(!si->is_restricted(*fit))
        continue;

      std::vector<Index> indices; //the other stars making up the boundary facet
      for(int i=0; i<3; ++i)
      {
        Index id = fit->first->vertex((fit->second+i+1)%4)->info();
        if(id != sid)
          indices.push_back(id);
      }

      //Ignore facets that are not incident to the star's center point
      if(indices.size() == 3)
        continue;

      if(is_boundary_facet_restricted_after_insertion(si, *fit))
        continue;

      check_from_facet(fit, indices);
    }
  }

  void check_from_cells(iterator mit,
                        const std::map<Cell_ijkl, int>& cells_counter)
  {
    Index sid = mit->first;
    Star_handle si = m_stars[sid];
    Czone& czi = mit->second;

    //compute cells that need a check
    Cell_handle_handle cit = czi.cells().begin();
    Cell_handle_handle cend = czi.cells().end();
    for(; cit!=cend; ++cit)
    {
      Cell_handle c = *cit;

      if(!si->is_inside(c)) //checks infinity too
        continue;
      std::map<Cell_ijkl, int>::const_iterator ccit = cells_counter.find(Cell_ijkl(c));
      if(ccit != cells_counter.end() && ccit->second == 4) //cell is in conflict in all stars
        continue;

      std::vector<Index> indices; //the other stars making up the cell
      for(int i=0; i<4; ++i)
      {
        Index id = c->vertex(i)->info();
        if(id != sid)
          indices.push_back(id);
      }

      //Ignore cells that are not incident to the star's center point (they do
      //not exist in the queues). They will be cleaned eventually
      if(indices.size() == 4)
        continue;

      for(std::size_t j=0; j<indices.size(); ++j)
      {
        Star_handle sj = m_stars[indices[j]];

        Cell_handle c_in_sj; //the cell c seen from sj
        // sj does not have c => nothing to do
        if(!sj->has_cell(Cell_ijkl(c), c_in_sj) || !sj->is_inside(c_in_sj))
          continue;

        // sj has c but sj is not in conflict => add to cells_to_check
        iterator czcit = m_conflict_zones.find(indices[j]);
        if(czcit == m_conflict_zones.end())
        {
          m_cells_to_check[indices[j]].push_back(c_in_sj);
          continue;
        }

        Czone czj = czcit->second;
        Cell_handle_handle chit = std::find(czj.cells().begin(),
                                            czj.cells().end(),
                                            c_in_sj);
        // sj has c and c is NOT in the czone of cj => add to cells_to_check
        if(chit == czj.cells().end())
          m_cells_to_check[indices[j]].push_back(c_in_sj);

        // sj has c and c is in the czone of cj => nothing to do
      }
    }
  }

  void compute_elements_needing_check()
  {
//todo -----------------------
    //WHEN REFINEMENT_CONDITIONS ARE ADDED AGAIN THERE NEEDS TO BE A CHECK HERE
    //NOT TO CONSIDER USELESS FACETS/CELLS
//-----------------------

    //What is being done here:
    //we need to make sure that no inconsistency appears in stars that are not
    //in conflict with the point. This might appear when a simplex is in conflict
    //in one star, but not in the others.

    if(m_are_checks_computed)
      return;


    //Count the elements. In the general case, an internal facet(cell) is three(four)
    //times in conflict, and there is thus nothing to check.
    std::map<Facet_ijk, int> internal_facets_counter;
    std::map<Cell_ijkl, int> cells_counter;

    iterator mit = m_conflict_zones.begin();
    iterator mend = m_conflict_zones.end();
    for(; mit!=mend; ++mit)
    {
      Star_handle si = m_stars[mit->first];
      Czone& czi = mit->second;
      if(czi.empty())
        continue;

      //internal facets
      Facet_handle fit = czi.internal_facets().begin();
      Facet_handle fend = czi.internal_facets().end();
      for(; fit!=fend; ++fit)
      {
        if(!si->is_restricted(*fit) &&
           (fit->first->vertex((fit->second+1)%4)->info() == mit->first || //incident facets only
            fit->first->vertex((fit->second+2)%4)->info() == mit->first ||
            fit->first->vertex((fit->second+3)%4)->info() == mit->first))
          add_to_map(Facet_ijk(*fit), internal_facets_counter);
      }

      if(!m_cells_need_check)
        continue;

      //cells
      Cell_handle_handle cit = czi.cells().begin();
      Cell_handle_handle cend = czi.cells().end();
      for(; cit!=cend; ++cit)
      {
        Cell_handle c = *cit;
        if(si->is_inside(c) &&
           (c->vertex(0)->info() == mit->first || //incident cells only!
            c->vertex(1)->info() == mit->first ||
            c->vertex(2)->info() == mit->first ||
            c->vertex(3)->info() == mit->first))
          add_to_map(Cell_ijkl(*cit), cells_counter);
      }
    }

    mit = m_conflict_zones.begin();
    for(; mit!=mend; ++mit)
    {
      check_from_internal_facets(mit, internal_facets_counter);
      check_from_boundary_facets(mit);

      if(m_cells_need_check)
        check_from_cells(mit, cells_counter);
    }

#ifdef ANISO_DEBUG_UNMODIFIED_STARS
    if(!m_facets_to_check.empty())
    {
      std::cout << "facets_to_check is not empty. Size (n° of stars): ': " << m_facets_to_check.size() << std::endl;
      std::cout << "Detail: " << std::endl;

      typename std::map<Index, std::vector<Facet> >::iterator mvfit = m_facets_to_check.begin();
      typename std::map<Index, std::vector<Facet> >::iterator mvfend = m_facets_to_check.end();
      for(; mvfit!= mvfend; ++mvfit)
      {
        Index sid = mvfit->first;
        std::cout << "-- Star " << sid << std::endl;
        std::vector<Facet> vf = mvfit->second;

        Facet_handle fit = vf.begin();
        Facet_handle fend = vf.end();
        for(; fit!=fend; ++fit)
        {
          std::cout << fit->first->vertex((fit->second+1)%4)->info() << " ";
          std::cout << fit->first->vertex((fit->second+2)%4)->info() << " ";
          std::cout << fit->first->vertex((fit->second+3)%4)->info() << " || ";
          std::cout << fit->first->vertex(fit->second)->info() << std::endl;
        }
      }
    }

    if(!m_cells_to_check.empty())
    {
      std::cout << "Cells_to_check is not empty. Size (n° of stars): ': " << m_cells_to_check.size() << std::endl;
      std::cout << "Detail: " << std::endl;

      typename std::map<Index, std::vector<Cell_handle> >::iterator mvchit = m_cells_to_check.begin();
      typename std::map<Index, std::vector<Cell_handle> >::iterator mvchend = m_cells_to_check.end();
      for(; mvchit!= mvchend; ++mvchit)
      {
        Index sid = mvchit->first;
        std::cout << "-- Star " << sid << std::endl;
        std::vector<Cell_handle> vch = mvchit->second;

        Cell_handle_handle vchit = vch.begin();
        Cell_handle_handle vchend = vch.end();
        for(; vchit!=vchend; ++vchit)
        {
          Cell_handle c = *vchit;
          std::cout << c->vertex(0)->info() << " ";
          std::cout << c->vertex(1)->info() << " ";
          std::cout << c->vertex(2)->info() << " ";
          std::cout << c->vertex(3)->info() << std::endl;
        }
      }
    }
#endif
    m_are_checks_computed = true;
  }

//Constructors
  Stars_conflict_zones(const Star_vector& m_stars_)
    :
      m_stars(m_stars_),
      m_conflict_zones(),
      m_conflict_p(Point(1e17, 1e17, 1e17)),
      m_conflict_p_id(-1),
      m_cells_need_check(false),
      m_are_checks_computed(false),
      m_facets_to_check(),
      m_cells_to_check()
  { }
};

template<typename K>
inline std::ostream& operator<<(std::ostream& os, Stars_conflict_zones<K>& src)
{
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
    os << czi.boundary_facets().size() << " ";
    os << czi.cells().size() << " ";
    os << czi.internal_facets().size() << " ";
    os << std::endl;
    os << "-*-*-*-*-" << std::endl;
    os << "-*-*-*-*-" << std::endl;
    if(0)
      continue;

    typename Stars_conflict_zones<K>::Facet_handle fit = czi.boundary_facets().begin();
    typename Stars_conflict_zones<K>::Facet_handle fend = czi.boundary_facets().end();
    for(; fit!=fend; ++fit)
    {
      os << " " << fit->first->vertex((fit->second+1)%4)->info() << " ";
      os << fit->first->vertex((fit->second+2)%4)->info() << " ";
      os << fit->first->vertex((fit->second+3)%4)->info() << " || ";
      os << fit->first->vertex(fit->second)->info() << std::endl;
    }
    os << "-----" << std::endl;

    typename Stars_conflict_zones<K>::Cell_handle_handle chit = czi.cells().begin();
    typename Stars_conflict_zones<K>::Cell_handle_handle chend = czi.cells().end();
    for(; chit!=chend; ++chit)
    {
      typename Stars_conflict_zones<K>::Cell_handle ch = *chit;
      os << " " << ch->vertex(0)->info() << " ";
      os << ch->vertex(1)->info() << " ";
      os << ch->vertex(2)->info() << " ";
      os << ch->vertex(3)->info() << std::endl;
    }
    os << "-----" << std::endl;

    fit = czi.internal_facets().begin();
    fend = czi.internal_facets().end();
    for(; fit!=fend; ++fit)
    {
      os << " " << fit->first->vertex((fit->second+1)%4)->info() << " ";
      os << fit->first->vertex((fit->second+2)%4)->info() << " ";
      os << fit->first->vertex((fit->second+3)%4)->info() << " || ";
      os << fit->first->vertex(fit->second)->info() << std::endl;
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

  typedef CGAL::Anisotropic_mesh_3::Conflict_zone<K>               Conflict_zone;
  typedef CGAL::Anisotropic_mesh_3::Stars_conflict_zones<K>        Stars_conflict_zones;

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

  Stars_conflict_zones& m_stars_czones;

public:
  double duration(const time_t& start) const
  {
    return ((clock() - start + 0.) / ((double)CLOCKS_PER_SEC));
  }

public:
  Star_vector& stars() const { return m_stars; }
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
  Star_handle get_star(typename Stars_conflict_zones::const_iterator it) const { return m_stars[it->first]; }

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

  //checks if a facet is encroached by any existing star (debug only)
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
#if 1//def ANISO_DEBUG_ENCROACHMENT
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
    std::set<Facet_ijk> facets;
    for(unsigned int i = 0; i < m_stars.size(); i++)
    {
      Facet_set_iterator fit = m_stars[i]->begin_restricted_facets();
      Facet_set_iterator fitend = m_stars[i]->end_restricted_facets();
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

  FT star_distortion(Star_handle star) const
  {
    FT max_distortion = 0.;

    /*
    Facet_set_iterator fit = star->begin_restricted_facets();
    Facet_set_iterator fitend = star->end_restricted_facets();
    for(; fit != fitend; fit++)
    {
      Facet f = *fit;
      max_distortion = (std::max)(max_distortion, compute_distortion(f));
    }
    return max_distortion;
  */

    Star_iterator vit = star->begin_neighboring_vertices();
    Star_iterator vend = star->end_neighboring_vertices();
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

  void debug_show_distortions() const
  {
    for(std::size_t i = 0; i < m_stars.size(); ++i)
    {
      Star_handle s = m_stars[i];
      if(!s->is_surface_star())
        continue;
      std::cout << "  " << i << " : ";
      Facet_set_iterator fit = s->begin_restricted_facets();
      Facet_set_iterator fend = s->end_restricted_facets();
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

  double average_facet_distortion(bool verbose = true) const
  {
    double avg_coh_dis = 0., min_coh_dis = 1e30, max_coh_dis = -1e30;
    double avg_incoh_dis = 0., min_incoh_dis = 1e30, max_incoh_dis = -1e30;
    int nb_coh = 0, nb_incoh = 0;

    std::set<Facet_ijk> done;
    for(std::size_t i = 0; i <number_of_stars(); i++)
    {
      Star_handle star = m_stars[i];
      if(!star->is_surface_star())
        continue;

      Facet_set_iterator fit = star->begin_restricted_facets();
      Facet_set_iterator fitend = star->end_restricted_facets();
      for(; fit != fitend; fit++)
      {
        Facet f = *fit;
        std::pair<typename std::set<Facet_ijk>::iterator, bool> is_insert_successful;
        is_insert_successful = done.insert(Facet_ijk(f));
        if(!is_insert_successful.second)
          continue;

        FT facet_distortion = compute_distortion(f);

        if(is_consistent(m_stars, f))
        {
          nb_coh++;
          avg_coh_dis += facet_distortion;
          if(facet_distortion > max_coh_dis)
            max_coh_dis = facet_distortion;
          if(facet_distortion < min_coh_dis)
            min_coh_dis = facet_distortion;
        }
        else
        {
          nb_incoh++;
          avg_incoh_dis += facet_distortion;
          if(facet_distortion > max_incoh_dis)
            max_incoh_dis = facet_distortion;
          if(facet_distortion < min_incoh_dis)
            min_incoh_dis = facet_distortion;
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

  double average_cell_distortion(bool verbose = true) const
  {
    double avg_coh_dis = 0., min_coh_dis = 1e30, max_coh_dis = -1e30;
    double avg_incoh_dis = 0., min_incoh_dis = 1e30, max_incoh_dis = -1e30;
    int nb_coh = 0, nb_incoh = 0;

    std::set<Cell_ijkl> done;
    for(std::size_t i = 0; i <number_of_stars(); i++)
    {
      Star_handle star = m_stars[i];
      if(!star->is_surface_star())
        continue;

      Cell_handle_handle ci = star->begin_finite_star_cells();
      Cell_handle_handle ciend = star->end_finite_star_cells();
      for(; ci != ciend; ci++)
      {
        Cell_handle c = *ci;
        if(!star->is_inside(c))
          continue;

        std::pair<typename std::set<Cell_ijkl>::iterator, bool> is_insert_successful;
        is_insert_successful = done.insert(Cell_ijkl(c));
        if(!is_insert_successful.second)
          continue;

        FT cell_distortion = compute_distortion(c);

        if(is_consistent(m_stars, c))
        {
          nb_coh++;
          avg_coh_dis += cell_distortion;
          if(cell_distortion > max_coh_dis)
            max_coh_dis = cell_distortion;
          if(cell_distortion < min_coh_dis)
            min_coh_dis = cell_distortion;
        }
        else
        {
          nb_incoh++;
          avg_incoh_dis += cell_distortion;
          if(cell_distortion > max_incoh_dis)
            max_incoh_dis = cell_distortion;
          if(cell_distortion < min_incoh_dis)
            min_incoh_dis = cell_distortion;
        }
      }
    }

    if(verbose)
    {
      std::cout << "SUM UP OF THE DISTORTION (CELL):" << std::endl;
      std::cout << nb_coh << " coherent cells with avg: " << avg_coh_dis/(double) nb_coh << std::endl;
      std::cout << "min: " << min_coh_dis << " max: " << max_coh_dis << std::endl;
      std::cout << nb_incoh << " incoherent cells with avg " << avg_incoh_dis/(double) nb_incoh << std::endl;
      std::cout << "min: " << min_incoh_dis << " max: " << max_incoh_dis << std::endl;
      std::cout << "---------------------------------------------" << std::endl;
    }

    return avg_coh_dis/(double) nb_coh;
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

      if(is_consistent(m_stars, si))
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
    std::cout << "clean call, size: " << m_stars.size() << std::endl;
    for(std::size_t i = 0; i < m_stars.size(); i++)
      m_stars[i]->clean();
  }

//Conflict zones
  Index simulate_insert_to_stars(const Point_3& p) const
  {
    Index this_id = static_cast<Index>(m_stars.size());

    // find conflicted stars
    Star_set stars;
    finite_stars_in_conflict(p, std::inserter(stars, stars.end())); //aabb tree
    infinite_stars_in_conflict(p, std::inserter(stars, stars.end())); //convex hull

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
      {
        //add them to the map of stars in conflict (conflict zones are not computed yet)
        m_stars_czones.conflict_zone(si->index_in_star_set());
      }
    }

    return this_id;
  }

  //those functions maybe could (should) be in the Czones class. But then it means
  //that all the filters (aabb, kd, etc.) need also a ref in the czone class(es) TODO (?)
  Index compute_conflict_zones(const Point_3& p) const
  {
    if(m_stars_czones.status() == Stars_conflict_zones::CONFLICT_ZONES_ARE_KNOWN
       && p == m_stars_czones.conflicting_point()) // just to be safe
    {
      std::cout << "Zones are already konwn, can skip computations!" << std::endl;
      return m_stars_czones.conflicting_point_id();
    }

    if(m_stars_czones.status() == Stars_conflict_zones::CONFLICTS_UNKNOWN)
      m_stars_czones.conflicting_point() = p;

    if(p != m_stars_czones.conflicting_point())
    {
      std::cout << "Points should have been equal there; ";
      std::cout << "the conflict zones are not cleaned correctly..." << std::endl;
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

      si->find_conflicts(p, std::back_inserter(ci.boundary_facets()),
                            std::back_inserter(ci.cells()),
                            std::back_inserter(ci.internal_facets()));
    }

    return m_stars_czones.conflicting_point_id();
  }

  void clear_conflict_zones() const
  {
    m_stars_czones.clear();
  }

//Point insertion functions
  template<typename Key, typename Val>
  struct set_map_comp
  {
    bool operator()(Key k, const std::pair<Key, Val>& p) const
    {
      return k < p.first;
    }

    bool operator()(const std::pair<Key, Val>& p, Key k) const
    {
      return p.first < k;
    }
  };

  void create_star(const Point_3 &p,
                   int pid,
                   Star_handle& star,
                   const bool surface_star = true,
                   const bool metric_reset = true) const
  {
#ifdef USE_ANISO_TIMERS
    std::clock_t start_time = clock();
#endif
     //TODO
    if(1) //metric_reset
    {
      if(1) //surface_star
      {
        Metric m_p;

        /* TODO
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

    if(m_stars_czones.is_empty())
      std::cout << "Warning: empty conflict map at the creation of the new star" << std::endl;

    typename Stars_conflict_zones::iterator czit = m_stars_czones.begin();
    typename Stars_conflict_zones::iterator czend = m_stars_czones.end();
    for(; czit!=czend; ++czit)
    {
      Index i = czit->first;
      Star_handle star_i = get_star(i);
      star->insert_to_star(star_i->center_point(), star_i->index_in_star_set(), false);
        //"false": no condition because they should be there for consistency
    }

    typename Star::Bbox bbox = star->bbox(); // volume bbox + surface bbox when star is not a topo_disk
    Point_3 pmin(bbox.xmin(), bbox.ymin(), bbox.zmin());
    Point_3 pmax(bbox.xmax(), bbox.ymax(), bbox.zmax());
    Kd_Box_query query(pmin, pmax, /*3=dim,*/ 0./*epsilon*/, typename Kd_tree::Star_pmap(m_stars));
    std::set<Kd_point_info> indices;
    m_kd_tree.search(std::inserter(indices, indices.end()), query);

    if(indices.size() != m_stars_czones.size())// because 'indices' contains 'modified_stars'
    {
      Index_set diff;
      std::set_difference(indices.begin(), indices.end(),
                          m_stars_czones.begin(), m_stars_czones.end(),
                          std::inserter(diff, diff.end()),
                          set_map_comp<Index, Conflict_zone>());
      typename Index_set::iterator it = diff.begin();
      while(it != diff.end())
      {
        Star_handle si = get_star(it++);
        star->insert_to_star(si->center_point(), si->index_in_star_set(), true/*conditional*/);
      }
    }
  }

  Star_handle create_star(const Point_3 &p,
                          int pid)
  {
    Star_handle star = new Star(m_criteria, m_pConstrain, true /*surface star*/);
    create_star(p, pid, star, true /*surface star*/);
    return star;
  }

  Star_handle create_inside_star(const Point_3 &p,
                                 int pid) const
  {
    Star_handle star = new Star(m_criteria, m_pConstrain, false/*not surface star*/);
    create_star(p, pid, star, false/*not surface star*/);
    return star;
  }

  //this version of "perform_insertions()" uses the conflict zones previous computed
  Index perform_insertions(const Point_3& p,
                           const Index& this_id,
                           const bool infinite_stars = false)
  {
    Vertex_handle v_in_ch = Vertex_handle();

    typename Stars_conflict_zones::iterator czit = m_stars_czones.begin();
    typename Stars_conflict_zones::iterator czend = m_stars_czones.end();
    for(; czit!=czend; ++czit)
    {
      Index i = czit->first;
      Star_handle si = get_star(i);
      Vertex_handle vi;

      Conflict_zone& si_czone = m_stars_czones.conflict_zone(i);
      Facet f = *(si_czone.boundary_facets().begin());
      vi = si->insert_to_star_hole(p, this_id, si_czone.cells(), f);

      if(vi == Vertex_handle()) //no conflict (should not happen)
      {
        std::cout << "No conflict in the insertion of a point in a conflict star" << std::endl;
        continue;
      }
      else if(vi->info() < this_id) // already in star set (should not happen)
      {
        std::cout << "Warning! Insertion of p"<< this_id
                  << " (" << p << ") in S" << si->index_in_star_set()
                  << " failed. vi->info() :"<< vi->info() << std::endl;
        if(v_in_ch != Vertex_handle())
          m_ch_triangulation.remove(v_in_ch);
        remove_from_stars(this_id, m_stars_czones.begin(), ++czit);
          //should probably be (++czit)-- but it doesn't really matter since something went wrong if we're in there

        si->print_vertices(true);
        si->print_faces();
        std::cout << "Metric : \n" << si->metric().get_transformation() << std::endl;
        return vi->info();
      }
      else // inserted, standard configuration
      {
        // update the triangulation of the convex hull
        if(infinite_stars)//(si->is_infinite() && v_in_ch != Vertex_handle())
        {
          v_in_ch = m_ch_triangulation.insert(p);
          v_in_ch->info() = this_id;
        }
      }
    }
    return this_id;
  }

  //this version is called when conditional is false (insert in target_stars regardless of conflicts)
  template<typename Stars>
  Index perform_insertions(const Point_3& p,
                           const Index& this_id,
                           const Stars& target_stars,
                           const bool infinite_stars = false)
  {
    std::cout << "perform insertions with cond false. nstars: " << target_stars.size() << std::endl;

    Vertex_handle v_in_ch = Vertex_handle();
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
        //Conflict zones are not computed for the target_stars, so we need to create
        //entries in the conflict zones map (for fill_ref_queue)
        m_stars_czones.conflict_zone(i);

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
                        const bool conditional)
  {
    Index this_id = static_cast<Index>(m_stars.size());

    Index id;
    if(conditional)
    {
      id = perform_insertions(p, this_id, false/*no convex hull upgrades*/); //stars in conflict

      //TODO CHECK THAT BELOW IS RELEVANT... IT SHOULD ALSO BE OKAY FOR IT TO BE THERE
      //RATHER THAN FOR BOTH VALUES OF "conditional" SINCE IF WE TAKE ALL STARS,
      //WE HAVE THE INFINITE STARS ALREADY INCLUDED...
      Star_set target_stars;
      infinite_stars_in_conflict(p, std::inserter(target_stars, target_stars.end())); //convex hull
      id = perform_insertions(p, this_id, target_stars, true/*update convex hull*/);
    }
    else
      id = perform_insertions(p, this_id, m_stars, true/*update convex hull*/); // insert in all stars

    return id;
  }

  Index insert(const Point_3 &p,
               const bool conditional,
               const bool surface_point = false)
  {
    std::cout << "insert call. p: " << p << " ";
    std::cout << "surface_point: " << surface_point << " ||  conditional: " << conditional << std::endl;

    if(conditional && m_stars_czones.status() != Stars_conflict_zones::CONFLICT_ZONES_ARE_KNOWN)
      std::cout << "Conflict zones unknown at insertion time...(insert)" << std::endl;

    Index id = insert_to_stars(p, conditional);
    if(id < 0 || id < (int)m_stars.size())
      return id;

    Star_handle star;
    if(surface_point)
      star = create_star(p, id); //implies surface star
    else
      star = create_inside_star(p, id);

    if(star->index_in_star_set() != m_stars.size())
      std::cout << "WARNING in insert..." << std::endl;

    if(conditional)
    {
      //need to have the id of the new star in the list of ids to check in fill_ref_queue
      m_stars_czones.conflict_zone(m_stars.size());
    }
    else
    {
      //if(!conditional), the queue is filled from all stars anyway so the czones
      //have no more utility and should be cleared ASAP.
      m_stars_czones.clear();
    }

    m_stars.push_back(star);
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
                           Kd_tree& m_kd_tree_,
                           Stars_conflict_zones& m_stars_czones_)
  :
    m_stars(m_stars_),
    m_pConstrain(m_pConstrain_),
    m_criteria(m_criteria_),
    m_metric_field(m_metric_field_),
    m_ch_triangulation(m_ch_triangulation_),
    m_aabb_tree(m_aabb_tree_),
    m_kd_tree(m_kd_tree_),
    m_stars_czones(m_stars_czones_)
  { }

};  // Anisotropic_refine_trunk

}  // Anisotropic_mesh_3
}  // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_REFINE_trunk_H

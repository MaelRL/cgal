#ifndef CGAL_ANISOTROPIC_MESH_3_REFINE_TRUNK_H
#define CGAL_ANISOTROPIC_MESH_3_REFINE_TRUNK_H

//#define ANISO_DEBUG_CZONES
#define ANISO_USE_BOUNDING_BOX_VERTICES_AS_POLES

#include <CGAL/Constrain_surface_3.h>
#include <CGAL/Criteria.h>
#include <CGAL/Metric_field.h>
#include <CGAL/Stretched_Delaunay_3.h>
#include <CGAL/Starset.h>

#include <CGAL/aabb_tree/aabb_tree_bbox.h>
#include <CGAL/aabb_tree/aabb_tree_bbox_primitive.h>

#include <CGAL/kd_tree/Kd_tree_for_star_set.h>

#include <CGAL/IO/Star_set_IO.h>

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
  typedef CGAL::Anisotropic_mesh_3::Starset<K, KExact>     Starset;
  typedef Star*                                            Star_handle;
  typedef std::vector<Star_handle>                         Star_vector;
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
  const Starset& m_starset;
  std::map<Index, Czone> m_conflict_zones;
  Point m_conflict_p; // point causing the conflict (only needed to debug / assert)
  Index m_conflict_p_id; // id of the point

  bool m_cells_need_check; // no need to compute cells_to_check if we haven't reached cell level
  bool m_are_checks_computed; // already computed or not

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

  bool empty() const { return m_conflict_zones.empty(); }
  void clear()
  {
    m_conflict_zones.clear();
    m_conflict_p = Point(1e17, 1e17, 1e17); // todo change that to something better
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
    for(std::size_t j=0, is=indices.size(); j<is; ++j)
    {
      Star_handle sj = m_starset[indices[j]];

      Facet f_in_sj; //the facet f seen from sj
      if(!sj->has_facet(Facet_ijk(*fit), f_in_sj)) //only searches through restricted facets
        continue;

      // At this point, we don't know for f(c,i) if c is a conflict cell

      iterator czcit = m_conflict_zones.find(indices[j]);
      if(czcit == m_conflict_zones.end())
      {
#ifdef ANISO_DEBUG_CZONES
        std::cout << "1. need to check @ " << indices[j] << " facet: ";
        std::cout << f_in_sj.first->vertex((f_in_sj.second+1)%4)->info() << " ";
        std::cout << f_in_sj.first->vertex((f_in_sj.second+2)%4)->info() << " ";
        std::cout << f_in_sj.first->vertex((f_in_sj.second+3)%4)->info() << " ";
        std::cout << f_in_sj.first->vertex(f_in_sj.second)->info() << std::endl;
        std::cout << "Original facet: ";
        std::cout << fit->first->vertex((fit->second+1)%4)->info() << " ";
        std::cout << fit->first->vertex((fit->second+2)%4)->info() << " ";
        std::cout << fit->first->vertex((fit->second+3)%4)->info() << " ";
        std::cout << fit->first->vertex(fit->second)->info() << std::endl;
        sj->print_facets();
#endif
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
#ifdef ANISO_DEBUG_CZONES
           std::cout << "2. need to check @ " << indices[j] << " facet: ";
           std::cout << f_in_sj.first->vertex((f_in_sj.second+1)%4)->info() << " ";
           std::cout << f_in_sj.first->vertex((f_in_sj.second+2)%4)->info() << " ";
           std::cout << f_in_sj.first->vertex((f_in_sj.second+3)%4)->info() << " ";
           std::cout << f_in_sj.first->vertex(f_in_sj.second)->info() << std::endl;
           std::cout << "Original facet: ";
           std::cout << fit->first->vertex((fit->second+1)%4)->info() << " ";
           std::cout << fit->first->vertex((fit->second+2)%4)->info() << " ";
           std::cout << fit->first->vertex((fit->second+3)%4)->info() << " ";
           std::cout << fit->first->vertex(fit->second)->info() << std::endl;
#endif
           m_facets_to_check[indices[j]].push_back(f_in_sj);
         }
         continue;
      }

      //if si's f is not found in sj's internals AND is not found in
      //sj's boundaries => add it

#ifdef ANISO_DEBUG_CZONES
      std::cout << "3. need to check @ " << indices[j] << " facet: ";
      std::cout << f_in_sj.first->vertex((f_in_sj.second+1)%4)->info() << " ";
      std::cout << f_in_sj.first->vertex((f_in_sj.second+2)%4)->info() << " ";
      std::cout << f_in_sj.first->vertex((f_in_sj.second+3)%4)->info() << " ";
      std::cout << f_in_sj.first->vertex(f_in_sj.second)->info() << std::endl;
      std::cout << "Original facet: ";
      std::cout << fit->first->vertex((fit->second+1)%4)->info() << " ";
      std::cout << fit->first->vertex((fit->second+2)%4)->info() << " ";
      std::cout << fit->first->vertex((fit->second+3)%4)->info() << " ";
      std::cout << fit->first->vertex(fit->second)->info() << std::endl;
#endif
      m_facets_to_check[indices[j]].push_back(f_in_sj);
    }
  }

  void check_from_internal_facets(iterator mit,
                                  const std::map<Facet_ijk, int>& internal_facets_counter)
  {
    Index sid = mit->first;
    Star_handle si = m_starset[sid];
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

      std::vector<Index> indices; // the other stars making up the facet
      for(int i=0; i<3; ++i)
      {
        Index id = fit->first->vertex((fit->second+i+1)%4)->info();
        if(id != sid)
          indices.push_back(id);
      }

      // Ignore facets that are not incident to the star's center point (they do
      // not exist in the queues). They will be cleaned eventually
      if(indices.size() == 3)
        continue;

      check_from_facet(fit, indices);
    }
  }

  void check_from_boundary_facets(iterator mit)
  {
    //boundary facets could stop being restricted and create inconsistencies
    Index sid = mit->first;
    Star_handle si = m_starset[sid];
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
    Star_handle si = m_starset[sid];
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

      for(std::size_t j=0, is=indices.size(); j<is; ++j)
      {
        Star_handle sj = m_starset[indices[j]];

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
      Star_handle si = m_starset[mit->first];
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

    //TODO : std unique on both if not empty since we can have duplicates
    //       make sure we actually have duplicates and they are actually
    //       adjacent first..

#ifdef ANISO_DEBUG_UNMODIFIED_STARS
    if(!m_facets_to_check.empty())
    {
      std::cout << "facets_to_check is not empty. Size (n° of stars): ': "
                << m_facets_to_check.size() << std::endl;
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
      std::cout << "Cells_to_check is not empty. Size (n° of stars): ': "
                << m_cells_to_check.size() << std::endl;
      std::cout << "Detail: " << std::endl;

      typename std::map<Index, std::vector<Cell_handle> >::iterator mvchit =
                                                       m_cells_to_check.begin();
      typename std::map<Index, std::vector<Cell_handle> >::iterator mvchend =
                                                       m_cells_to_check.end();
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
  Stars_conflict_zones(const Starset& m_starset_)
    :
      m_starset(m_starset_),
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
  typedef typename Star::Facet                                     Facet;
  typedef typename Star::Facet_handle                              Facet_handle;
  typedef typename Star::Facet_set_iterator                        Facet_set_iterator;
  typedef typename Star::Cell_handle_handle                        Cell_handle_handle;
  typedef typename Star::Vector_3                                  Vector_3;
  typedef typename Star::Constrain_surface                         Constrain_surface;
  typedef typename Star::Criteria                                  Criteria;

  typedef CGAL::Anisotropic_mesh_3::Starset<K>                     Starset;

  typedef CGAL::Anisotropic_mesh_3::Metric_field<K>                Metric_field;
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
  // the following boolean is useful for different purposes:
  //    --used to indicate whether non-surfacic points are poles (set metric to identity) or
  //      ''real'' points with a metric that needs to be computed from the MF
  //    --used to indicate to a surface level that it should use 3D checks for pick_valid
  //    --used to indicate to a surface level that the stars should use volume bboxes
  mutable bool m_is_3D_level;

  Starset& m_starset;

  const Constrain_surface* m_pConstrain;
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
  const Constrain_surface* constrain_surface() const { return m_pConstrain; }
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

  Point_3 barycenter(const Facet& f) const
  {
    Point_3 p1 = get_star(f.first->vertex((f.second + 1) % 4)->info())->center_point();
    Point_3 p2 = get_star(f.first->vertex((f.second + 2) % 4)->info())->center_point();
    Point_3 p3 = get_star(f.first->vertex((f.second + 3) % 4)->info())->center_point();
    FT third = 1./3.;
    return Point_3(third * (p1.x() + p2.x() + p3.x()),
                   third * (p1.y() + p2.y() + p3.y()),
                   third * (p1.z() + p2.z() + p3.z()));
  }

  FT sq_distance_to_surface(const Facet& f, const Star_handle /*s*/) const
  {
    //Point_3 steiner;
    //s->compute_dual_intersection(f, steiner);
    //Point_3 cc = transform_from_star_point(s->compute_circumcenter(f), s);
    //return CGAL::squared_distance(steiner, cc);

    return m_pConstrain->compute_sq_approximation(barycenter(f));
  }

  bool is_surface_facet(const Facet& f) const
  {
    return get_star(f.first->vertex((f.second + 1) % 4)->info())->is_surface_star()
        && get_star(f.first->vertex((f.second + 2) % 4)->info())->is_surface_star()
        && get_star(f.first->vertex((f.second + 3) % 4)->info())->is_surface_star();
  }

  //checks if a facet is encroached by any existing star (debug only)
  bool is_encroached(Star_handle star, const Facet &facet)
  {
    Index p1 = facet.first->vertex((facet.second + 1)%4)->info();
    Index p2 = facet.first->vertex((facet.second + 2)%4)->info();
    Index p3 = facet.first->vertex((facet.second + 3)%4)->info();

    Star_iterator sit = m_starset.begin();
    Star_iterator sitend = m_starset.end();
    for(; sit!=sitend; ++sit)
    {
      Star_handle si = get_star(sit);

      Index nsi = si->index_in_star_set();
      if(nsi == p1 || nsi == p2 || nsi == p3)
        continue;

      if(star->is_facet_encroached(transform_to_star_point(si->center_point(), star), facet))
      {
#ifdef ANISO_DEBUG_ENCROACHMENT
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

//aabb functions
public:
  void update_aabb_tree(Star_handle star) const
  {
    m_aabb_tree.update_primitive(AABB_primitive(star));
  }

  void build_aabb_tree()
  {
    m_aabb_tree.rebuild(m_starset.begin(), m_starset.end());
  }

  void update_bboxes() const
  {
#ifdef DEBUG_UPDATE_AABB_TREE
    std::cout << "updating bboxes. tree : " << m_aabb_tree.size() << " " << m_aabb_tree.m_insertion_buffer_size() << std::endl;
#endif
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
  void finite_stars_in_conflict(const Point_3& p,
                                OutputIterator oit) const
  {
    update_bboxes();
#ifndef NO_USE_AABB_TREE_OF_BBOXES
    // bigger set of stars
    m_aabb_tree.all_intersected_primitives(p, oit);
#else
    //exact set
    for(std::size_t i = 0; i < number_of_stars(); i++)
    {
      Star_handle s = get_star(i);
      Cell_handle ch;
      if(s->is_conflicted(transform_to_star_point(p, s), ch))
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
  void cells_created(Star_handle star,
                     std::map<Cell_ijkl, int>& cells) const
  {
    Cell_handle_handle ci = star->finite_star_cells_begin();
    Cell_handle_handle ciend = star->finite_star_cells_end();
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
    typename std::vector<Facet>::const_iterator end = bfacets.end();
    for( ; fit!=end; fit++)
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

      Point_3 c = compute_circumcenter(this->m_starset[id1]->center_point(),
                                       this->m_starset[id2]->center_point(),
                                       this->m_starset[id3]->center_point(),
                                       p, star);
      if(is_inside_domain(c))
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
      Star_handle si = get_star(i);
      cells_created(p, p_index, si, cells);
    }

    /*
    std::cout << "cells sum up predicts : " << cells.size() << std::endl;
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
      int nmax = number_of_stars();
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
        tp0 = (c0 == nmax) ? tp : this->transform_to_star_point(this->m_starset[c0]->center_point(),
                                                                to_be_refined);
        tp1 = (c1 == nmax) ? tp : this->transform_to_star_point(this->m_starset[c1]->center_point(),
                                                                to_be_refined);
        tp2 = (c2 == nmax) ? tp : this->transform_to_star_point(this->m_starset[c2]->center_point(),
                                                                to_be_refined);
        tp3 = (c3 == nmax) ? tp : this->transform_to_star_point(this->m_starset[c3]->center_point(),
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
          star = this->m_starset[(*itc).first.vertex(i)];

        Point_3 p0, p1, p2, p3;
        p0 = (c0 == nmax) ? p : this->m_starset[c0]->center_point();
        p1 = (c1 == nmax) ? p : this->m_starset[c1]->center_point();
        p2 = (c2 == nmax) ? p : this->m_starset[c2]->center_point();
        p3 = (c3 == nmax) ? p : this->m_starset[c3]->center_point();

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

  bool is_valid_point_3D(const Point_3 &p,
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

    create_star(p, id, new_star, new_star->is_surface_star());

    bool is = check_consistency_and_sliverity(to_be_refined, new_star, sq_radiusbound);
    if(!is)
      new_star->invalidate();

    return is;
  }

//Conflict zones
  Index simulate_insert_to_stars(const Point_3& p) const
  {
    int s_in_conflict = 0.;
    std::set<Cell_ijkl> cell_set;

    Index this_id = static_cast<Index>(number_of_stars());

    // find conflicted stars
    Star_set stars;
    finite_stars_in_conflict(p, std::inserter(stars, stars.end())); //aabb tree

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

    return this_id;
  }

  //those functions maybe could (should) be in the Czones class. But then it means
  //that all the filters (aabb, kd, etc.) need also a ref in the czone class(es) TODO (?)
  Index compute_conflict_zones(const Point_3& p) const
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

// Point insertion functions
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

  void insert_from_kd_tree(Star_handle star) const
  {
    int counter = 0;
    bool grew = true;
    while(grew)
    {
      counter++;

      const typename Star::Bbox& bbox = star->bbox();
      Point_3 pmin(bbox.xmin(), bbox.ymin(), bbox.zmin());
      Point_3 pmax(bbox.xmax(), bbox.ymax(), bbox.zmax());

      Kd_Box_query query(pmin, pmax, /*3=dim,*/ 0./*epsilon*/,
                         typename Kd_tree::Star_pmap(m_starset.star_vector()));
      std::set<Kd_point_info> indices;
      m_kd_tree.search(std::inserter(indices, indices.end()), query);

#ifdef ANISO_DEBUG_KDTREE
      std::cout << "star of size: " << star->number_of_vertices() << std::endl;
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

  void create_star(const Point_3 &p,
                   int pid,
                   Star_handle& star,
                   const bool surface_star,
                   const bool reset_star = true) const
  {
#ifdef USE_ANISO_TIMERS
    std::clock_t start_time = clock();
#endif
    if(reset_star)
    {
      if(m_is_3D_level || surface_star)
      {
        Metric m_p = metric_field()->compute_metric(p);
        star->reset(p, pid, m_p, surface_star);
      }
      else
      {
        star->reset(p, pid, metric_field()->uniform_metric(p), surface_star);
      }
    }

#ifdef ANISO_USE_BOUNDING_BOX_VERTICES_AS_POLES
    //adding the 8 bounding box vertices (they are the first 8 stars)
    if(this->m_starset.size() > 8)
    {
      for(int i=0; i<8; ++i)
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
      star->insert_to_star(star_i->center_point(), star_i->index_in_star_set(),
                           false/*conditional*/);
        // "false": no condition because they should be there for consistency
    }

    insert_from_kd_tree(star);
    star->clean();

#ifdef ANISO_DEBUG_INSERT
    star->print_vertices();
    std::ofstream out_before("star_before.off");
    output_star_off(star, out_before);

    //check if we missed any other star in conflict
    typename Star_vector::iterator it = stars().begin();
    typename Star_vector::iterator iend = stars().end();
    for(; it!=iend; it++)
    {
      Star_handle star_i = get_star(it);
      if(star->has_vertex(star_i->index_in_star_set()))
        continue;

      const Point_3& cpi = star_i->center_point();
      const TPoint_3& tcpi = star->metric().transform(cpi);
      Cell_handle dummy;
      CGAL_assertion(!star->is_conflicted(tcpi, dummy));
    }

    // check if rebuilding with everything would give the expected same result
    if(m_is_3D_level)
    {
      // can't use this test if we're not using 3D bboxes & stuff like that since
      // we intentionally restrict stars when in a mesher level that is only 2D

      std::size_t mem = star->finite_adjacent_vertices_end() -
                        star->finite_adjacent_vertices_begin();
      star->reset();
      for(std::size_t i=0; i<m_starset.size(); ++i)
      {
        Star_handle star_i = m_starset[i];
        star->insert_to_star(star_i->center_point(), i, false/*no cond*/);
      }
      star->clean(true/*verbose*/);
      star->print_vertices();

      std::ofstream out_after("star_after.off");
      output_star_off(star, out_after);

      CGAL_postcondition(mem == (star->finite_adjacent_vertices_end() -
                                 star->finite_adjacent_vertices_begin()));
    }
#endif
  }

  Star_handle create_star(const Point_3 &p,
                          int pid)
  {
    Star_handle star = new Star(m_criteria, m_pConstrain, true /*surface star*/, m_is_3D_level);
    create_star(p, pid, star, true /*surface star*/);
    return star;
  }

  Star_handle create_inside_star(const Point_3 &p,
                                 int pid) const
  {
    Star_handle star = new Star(m_criteria, m_pConstrain, false/*not surface star*/, m_is_3D_level);
    create_star(p, pid, star, false/*not surface star*/);
    return star;
  }

  //this version of "perform_insertions()" uses the conflict zones previous computed
  Index perform_insertions(const Point_3& p,
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
        remove_from_stars(this_id, m_stars_czones.begin(), ++czit);
        // should probably be (++czit)-- but it doesn't really matter
        // since something went wrong if we're in there

        si->print_vertices(true);
        si->print_facets();
        std::cout << "Metric : \n" << si->metric().get_transformation() << std::endl;
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

#ifdef ANISO_DEBUG_INSERT
    // check if we have inserted the new point everywhere it's needed
    typename Star_vector::iterator it = stars().begin();
    typename Star_vector::iterator iend = stars().end();
    for(; it != iend; it++)
    {
      Star_handle star_i = get_star(it);
      Index si = star_i->index_in_star_set();

      if(m_stars_czones.conflict_zones().find(si) != m_stars_czones.conflict_zones().end())
        continue;

      const TPoint_3& tpi = star_i->metric().transform(p);
      Cell_handle dummy;
      CGAL_assertion(!star_i->is_conflicted(tpi, dummy));
    }

    // check that the stars that were in conflict for the new point are correct (after insertion)
    czit = m_stars_czones.begin();
    for(; czit!=czend; ++czit)
    {
      Index i = czit->first;
      Star_handle star_i = get_star(i);

      it = stars().begin();
      for(; it != iend; it++)
      {
        Star_handle star_j = get_star(it);

        if(star_i->has_vertex(star_j->index_in_star_set()))
          continue;

        const Point_3& cpj = star_j->center_point();
        const TPoint_3& tcpj = star_i->metric().transform(cpj);
        Cell_handle dummy;
        CGAL_assertion(!star_i->is_conflicted(tcpj, dummy));
      }
    }
#endif

    return this_id;
  }

  //this version is called when conditional is false (insert in target_stars regardless of conflicts)
  template<typename Stars>
  Index perform_insertions(const Point_3& p,
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
        si->print_facets();
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

  Index insert_to_stars(const Point_3& p,
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

  Index insert(const Point_3& p,
               const bool conditional,
               const bool surface_point = false)
  {
    if(conditional && m_stars_czones.status() != Stars_conflict_zones::CONFLICT_ZONES_ARE_KNOWN)
      std::cout << "Conflict zones unknown at insertion time...insert()" << std::endl;

    Index id = insert_to_stars(p, conditional);
    if(id < 0 || id < static_cast<int>(number_of_stars()))
      return id;

    Star_handle star;
    if(surface_point)
      star = create_star(p, id); // implies surface star
    else
      star = create_inside_star(p, id);

    if(star->index_in_star_set() != static_cast<Index>(number_of_stars()))
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
    m_kd_tree.insert(star->index_in_star_set());
#ifndef NO_USE_AABB_TREE_OF_BBOXES
    m_aabb_tree.insert(AABB_primitive(star));
#endif

    return id;
  }

private:
  //Adding the 8 points of the domain's bbox to avoid infinite cells
  void initialize_bounding_box_vertices()
  {
    std::vector<Point_3> bbox_vertices;
    Bbox_3 bbox = this->m_pConstrain->get_bbox();

    FT xmin = bbox.xmin(), xmax = bbox.xmax();
    FT ymin = bbox.ymin(), ymax = bbox.ymax();
    FT zmin = bbox.zmin(), zmax = bbox.zmax();

    bbox_vertices.push_back(Point_3(xmin, ymin, zmin));
    bbox_vertices.push_back(Point_3(xmax, ymin, zmin));
    bbox_vertices.push_back(Point_3(xmax, ymin, zmax));
    bbox_vertices.push_back(Point_3(xmin, ymin, zmax));
    bbox_vertices.push_back(Point_3(xmin, ymax, zmin));
    bbox_vertices.push_back(Point_3(xmax, ymax, zmin));
    bbox_vertices.push_back(Point_3(xmax, ymax, zmax));
    bbox_vertices.push_back(Point_3(xmin, ymax, zmax));

    typename std::vector<Point_3>::iterator it = bbox_vertices.begin();
    typename std::vector<Point_3>::iterator end = bbox_vertices.end();
    for(; it!=end; ++it)
      insert(*it, false /*no condition*/, false/*not surface*/);
  }

protected:
  //Initialization
  void initialize_medial_axis()
  {
    Point_set poles;
    this->m_pConstrain->compute_poles(poles);
#ifdef ANISO_VERBOSE
    std::cout << "Initialize medial axis..." << std::endl;
    std::cout << poles.size() << " poles before resize..." << std::endl;
#endif

    //TEMP ---------------------------------------------------------------------
    //ugly code to reduce the number of poles to the desired fixed amount
    std::vector<Point_3> poles_v;
    std::copy(poles.begin(), poles.end(), std::back_inserter(poles_v));
    std::random_shuffle(poles_v.begin(), poles_v.end());
    poles_v.resize(std::min(static_cast<std::size_t>(150), poles_v.size()));
    //-------------------------------------------------------------------------

#ifdef ANISO_VERBOSE
    //insert them in all stars
    std::cout << "Insert poles..." << std::endl;
#endif

    unsigned int i = 1;
    unsigned int done = 0;
    typename std::vector<Point_3>::iterator it = poles_v.begin();
    typename std::vector<Point_3>::iterator end = poles_v.end();
    for(; it!=end; ++it, ++i)
    {
      bool conditional = (i % 10 != 0); // 1/10 with no condition
      if(true) //m_refinement_condition(*it)) TODO
      {
        insert(*it, false /*no condition*/, false/*surface point*/); // tmp
        ++done;
      }
    }

    clean_stars();
#ifdef ANISO_VERBOSE
    std::cout << done << " points." << std::endl;
#endif
  }

  void initialize_stars(const bool are_poles_used = false,
                        const int nb = 50)
  {
#ifdef ANISO_VERBOSE
    std::cout << "Initialize " << nb << " stars..." << std::endl;
#endif
    double approx = this->m_criteria->approximation/this->m_pConstrain->get_bounding_radius();
    approx = 1e-4;

#ifdef ANISO_USE_BOUNDING_BOX_VERTICES_AS_POLES
    initialize_bounding_box_vertices();
#endif

    //The initial points need to be picked more cleverly as they completely ignore
    //the input metric field right now TODO
    typename Constrain_surface::Pointset initial_points =
      this->m_pConstrain->get_surface_points(2*nb);

    //only used in the case of a pure surface meshing process
    if(are_poles_used)
      initialize_medial_axis();

#ifdef ANISO_VERBOSE
    std::cout << "Picked " << initial_points.size() << " initial points" << std::endl;
    for(typename Constrain_surface::Pointset::iterator it=initial_points.begin(); it!=initial_points.end(); ++it)
      std::cout << it->x() << " " << it->y() << " " << it->z() << std::endl;
#endif

    typename Constrain_surface::Pointset::iterator pi = initial_points.begin();
    typename Constrain_surface::Pointset::iterator pend = initial_points.end();
    int nbdone = 0;
    for (; pi != pend && nbdone < nb; pi++)
    {
      if(nbdone > 0 && nbdone % 100 == 0)
        clean_stars();

      Index this_id = static_cast<int>(this->m_starset.size());

      //if(m_refinement_condition(*pi)) TODO
      Index id = insert(*pi, false/*under no condition*/, true/*surface point*/);

      if(this_id == id)
        nbdone++;

      /*
        //this would work to stop inserting unneeded initial points, but it's expensive
        this->clean_stars();
        if(this->count_restricted_facets() > 0)
          break;
        */
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

    //if resuming for a surface, add poles
    if(!m_is_3D_level)
    {
#ifdef ANISO_USE_BOUNDING_BOX_VERTICES_AS_POLES
      initialize_bounding_box_vertices();
#endif
      this->m_pConstrain->get_surface_points(50); // initial points are not used
      initialize_medial_axis();
      build_aabb_tree();
    }

    std::ifstream in(filename);
    std::string word;
    int useless, nv;
    FT x,y,z;

    in >> word >> useless; //MeshVersionFormatted i
    in >> word >> useless; //Dimension d
    in >> word >> nv;

    std::cout << "filename: " << filename << std::endl;
    std::cout << nv << " vertices in file" << std::endl;

    for(int i=0; i<nv; ++i)
    {
      in >> x >> y >> z >> useless;
      Point_3 p(x,y,z);

      if(m_starset.size() < 10)
        insert(p, false /*no condition*/, !m_is_3D_level/*surface*/);
      else
      {
        if(m_starset.size() == 10)
          build_aabb_tree();

        compute_conflict_zones(p);
        m_stars_czones.compute_elements_needing_check();
        insert(p, true/*conditional*/, !m_is_3D_level/*surface*/);
      }
      clear_conflict_zones();
    }
    clean_stars();

    std::cout << "consistency of the triangulation: " << m_starset.is_consistent(true, FACETS_ONLY);
    std::cout << " " << m_starset.is_consistent(true, CELLS_ONLY) << std::endl;

    std::ofstream out("resumed.mesh");
    output_medit(m_starset, out);
    std::ofstream out_facet("resumed_surf.mesh");
    output_surface_medit(m_starset, out_facet);
  }

  // it is way more efficient to use below
  void resume_from_dump_file_(const char* filename)
  {
    if(!m_starset.empty())
    {
      std::cout << "resuming with a non empty star set... ?" << std::endl;
      return;
    }

    // if resuming for a surface, add poles
    if(!m_is_3D_level)
    {
#ifdef ANISO_USE_BOUNDING_BOX_VERTICES_AS_POLES
      initialize_bounding_box_vertices();
#endif
      this->m_pConstrain->get_surface_points(50); // initial points are not used
      initialize_medial_axis();
      build_aabb_tree();
    }

    std::ifstream in(filename);
    std::size_t stars_n, v_n, id; // number of stars, neighbors in a star, neighbhor id
    FT x,y,z; // coordinates

    in >> stars_n;
    std::cout << "dump file has: " << stars_n << " stars" << std::endl;

    for(std::size_t i=0; i<stars_n; ++i)
    {
      in >> x >> y >> z;
      Point_3 p(x,y,z);

      if(i%1000 == 0)
        std::cout << i << " stars" << std::endl;

      Star_handle star = new Star(m_criteria, m_pConstrain);
      const Metric& m_p = m_metric_field->compute_metric(p);
      star->reset(p, i, m_p);
      m_starset.push_back(star);
    }

    std::cout << "built stars, now filling neighbors" << std::endl;

    for(std::size_t i=0; i<stars_n; ++i)
    {
      if(i%1000 == 0)
        std::cout << i << " neighborhood" << std::endl;


      in >> v_n;
      Star_handle star_i = m_starset.get_star(i);
      for(std::size_t j=0; j<v_n; ++j)
      {
        in >> id;
        Star_handle star_j = m_starset.get_star(id);
        star_i->insert_to_star(star_j->center_point(), id, false);
      }
    }

    build_aabb_tree();

    std::cout << "starset of size: " << m_starset.size() << " stars from " << filename << std::endl;
  }

protected:
  void switch_to_volume_bboxes()
  {
    Star_iterator sit = m_starset.begin();
    Star_iterator sitend = m_starset.end();
    for(; sit!=sitend; ++sit)
    {
      Star_handle star_i = get_star(sit);
      if(!star_i->is_in_3D_mesh())
      {
        star_i->is_in_3D_mesh() = true;
        star_i->invalidate_bbox_cache();
        star_i->update_bbox();
      }
    }
  }

protected:
  //Attempt consistency by starting from an isotropic mesh and
  //progressively increasing anisotropy
  void attempt_consistency(Consistency_check_options c_opts)
  {
    std::cout << "attempt consistency" << std::endl;

    int const pass_n = 1;
    std::vector<double> coeff(number_of_stars(), 0.);
    std::vector<double> consistent_at_last_pass(number_of_stars(), true);
    std::vector<Eigen::Matrix3d> memory(number_of_stars());
    double var = 0.1;

    for(int i=0; i<pass_n; ++i)
    {
      std::cout << "i: " << i << std::endl;

      for(std::size_t si=0; si<number_of_stars(); ++si)
      {
        Star_handle star_i = get_star(si);
        if(i==0)
        {
          memory[si] = star_i->metric().get_transformation();
          std::cout << memory[si] << std::endl;
        }

        double alpha = coeff[si];
        if(consistent_at_last_pass[si])
          alpha += var;
        else
          alpha -= var;

        alpha = std::max(alpha, 0.);
        alpha = std::min(alpha, 1.);
// -----------
        assert(alpha >= 0 && alpha <= 1);
        std::cout << "si: " << si << " alpha: " << coeff[si] << " " << alpha << " c: ";
        std::cout << consistent_at_last_pass[si] << std::endl;
// -----------
        coeff[si] = alpha;

        std::vector<std::pair<Eigen::Matrix3d, double> > w_m;
        w_m.push_back(std::make_pair(memory[si], alpha));
        w_m.push_back(std::make_pair(Eigen::Matrix3d::Identity(), 1.-alpha));
        Eigen::Matrix3d mi = linear_interpolate<K>(w_m);

        std::cout << mi << std::endl;

        Point_3 const pi = star_i->center_point();
        Metric const Mi(mi);
        bool const bi = star_i->is_surface_star();

        star_i->reset(pi, si, Mi, bi);

        for(std::size_t sj=0.; sj<number_of_stars(); ++sj)
        {
          Star_handle star_j = get_star(sj);
          star_i->insert_to_star(star_j->center_point(), sj, false/*no cond*/);
        }
        star_i->clean();
      }

      int count = 0;
      for(std::size_t si=0.; si<number_of_stars(); ++si)
      {
        // ONLY CHECKS CELL CONSISTENCY !!!

        bool is_c = true;
        Star_handle star_i = get_star(si);
        Cell_handle_handle cit = star_i->finite_star_cells_begin();
        Cell_handle_handle citend = star_i->finite_star_cells_end();
        for(; cit != citend; cit++)
        {
          if(!star_i->is_inside(*cit))
            continue;
          if(!m_starset.is_consistent(*cit, true/*verbose*/))
          {
            is_c = false;
            break;
          }
        }

        if(!is_c) //might be a quasi-cosphericity, let's move the point
        {
          Cell_handle c = *cit;
          Index i0 = c->vertex(0)->info();
          Star_handle star_0 = get_star(i0);

          int max_shake = 100;
          int shake_n = 0;
          while(max_shake > shake_n++ &&
                !m_starset.is_consistent(star_i, true/*verbose*/, c_opts, i0))
          {
            star_0->shake_center();
            std::cout << "shake : " << i0 << std::endl;
            star_0->reset();

            m_starset.rebuild();
          }
        }

        if(!is_c) ++count;
        consistent_at_last_pass[si] = is_c;
      }
      std::cout << "inconsistent stars: " << count << std::endl;

      std::ofstream out("bambimboum_c.mesh");
      output_medit(m_starset, out, false);
    }
  }

//Constructors
public:
  Anisotropic_refine_trunk(Starset& starset_,
                           const Constrain_surface* pConstrain_,
                           const Criteria* criteria_,
                           const Metric_field* metric_field_,
                           AABB_tree& aabb_tree_,
                           Kd_tree& kd_tree_,
                           Stars_conflict_zones& stars_czones_,
                           bool is_3D_level_)
  :
    m_is_3D_level(is_3D_level_),
    m_starset(starset_),
    m_pConstrain(pConstrain_),
    m_criteria(criteria_),
    m_metric_field(metric_field_),
    m_aabb_tree(aabb_tree_),
    m_kd_tree(kd_tree_),
    m_stars_czones(stars_czones_)
  { }

};  // Anisotropic_refine_trunk

}  // Anisotropic_mesh_3
}  // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_REFINE_trunk_H

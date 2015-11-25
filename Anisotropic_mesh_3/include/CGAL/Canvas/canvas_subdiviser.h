#ifndef CGAL_ANISOTROPIC_MESH_3_Canvas_subdivider_H
#define CGAL_ANISOTROPIC_MESH_3_Canvas_subdivider_H

#include <CGAL/Canvas/canvas_config.h>
#include <CGAL/Canvas/canvas_enum.h>

#include <CGAL/assertions.h>
#include <CGAL/enum.h>
#include <CGAL/Kernel/global_functions_3.h>

#include <boost/array.hpp>
#include <boost/tuple/tuple.hpp>

#include <list>
#include <map>
#include <set>
#include <utility>
#include <vector>

// This class implements barycentric subdivision of tetrahedra in the canvas context

// A criterion is used to determine whether a cell should be subdivided or not.
// Calling 'subdivide_cells' applies at most ONE subdivision to a cell.

// Cells that are neigbors to subdivided cells need a special treatment to ensure
// that we obtain a proper triangulation.

namespace CGAL
{
namespace  Anisotropic_mesh_3
{

template<class Canvas>
struct Subdomain_criterion
{
  typedef typename Canvas::Tr                                     Tr;
  typedef typename Tr::Vertex_handle                              Vertex_handle;
  typedef typename Tr::Cell_handle                                Cell_handle;

  Canvas* m_canvas;
  std::set<int> seed_indices;

  bool operator()(const Cell_handle c) const
  {
    // returns true if all the vertices's closest seed ids are in subdomain_indices

    for(int i=0; i<4; ++i)
    {
      Vertex_handle vh = c->vertex(i);
      int seed_index = (m_canvas->canvas_points[vh->info()]).closest_seed_id();
      if(seed_indices.find(seed_index) == seed_indices.end())
        return false;
    }
    return true;
  }

  Subdomain_criterion(Canvas* canvas) : m_canvas(canvas) { }
};

template<class Canvas, class Criterion>
struct Canvas_subdivider
{
  typedef typename Canvas::FT                                     FT;
  typedef typename Canvas::Tr                                     Tr;
  typedef typename Canvas::Canvas_point                           Canvas_point;

  typedef typename Tr::Point                                      Point_3;
  typedef typename Tr::Vertex_handle                              Vertex_handle;
  typedef typename Tr::Facet                                      Facet;
  typedef typename Tr::Cell_handle                                Cell_handle;
  typedef std::set<Vertex_handle>                                 Vertex_handle_set;

  typedef std::vector<Cell_handle>                                Cell_handle_vector;
  typedef typename Cell_handle_vector::iterator                   Cell_handle_handle;

  Canvas* m_canvas;
  const Criterion& m_criterion;
  std::set<int> seeds_to_reset;

  // Once the cells are split, we need to split the faces to get a pretty
  // barycentric subdivision
  // a naive way is to grab a set of triplets <int,int,int> and loop cells,
  // searching for those facets...
  // Another better way is to associate to each cell the info 'should one your
  // facet be split and if so which one ?'
  // The problem is that our cell_with_info's info() is only one int (the subdomain index),
  // so we temporarily replace that subdomain index with an index to a vector
  // that will contain additional info (it's kinda dirty, but it works)

  // To make it even dirtier we have to transform the index to that vector
  // so it is not confused with the normal subdomain index of a cell that mustn't
  // be refined. In the end the possibilities for info() are :
  // ]-inf, -2] a this cell must be split and index is in the vector is j = -(i+2)
  // -1 is an infinite cell,
  // 0 is a finite cell outside the complex that won't be split
  // [1, inf[ is a finite cell in the complex  that won't be split

  // For the tuple, the first int is the index that identifies which facet needs
  // to be split // (there is always only one)
  // the second int is the cell's real info() : the subdomain index that we must
  // to keep in memory and re assign later
  std::vector<boost::tuple<Cell_handle, int, int> > extra_info;
  std::list<Cell_handle> cells_with_a_face_to_divide;

  // subdivision functions -----------------------------------------------------

  Tr& tr() { return m_canvas->m_tr; }
  const Tr& tr() const { return m_canvas->m_tr; }

  bool should_cell_be_subdivided(const Cell_handle c) const
  {
    if(tr().is_infinite(c))
      return false;
    return m_criterion(c);
  }

  Point_3 facet_barycenter(Cell_handle c, int i)
  {
    CGAL_precondition(!tr().is_infinite(c));

    const Point_3& p1 = c->vertex((i+1)%4)->point();
    const Point_3& p2 = c->vertex((i+2)%4)->point();
    const Point_3& p3 = c->vertex((i+3)%4)->point();

    FT third = 1./3.;
    return CGAL::barycenter(p1, third, p2, third, p3, third);
  }

  Point_3 cell_barycenter(Cell_handle c)
  {
    CGAL_precondition(!tr().is_infinite(c));
    return CGAL::barycenter(c->vertex(0)->point(), 0.25,
                            c->vertex(1)->point(), 0.25,
                            c->vertex(2)->point(), 0.25,
                            c->vertex(3)->point(), 0.25);
  }

  void set_info_on_edge_incident_cells(Vertex_handle opp,
                                       Vertex_handle vh, int info)
  {
    std::vector<Cell_handle> incident_cells;
    tr().incident_cells(opp, std::back_inserter(incident_cells));

    Cell_handle_handle chit = incident_cells.begin();
    for(; chit!=incident_cells.end(); ++chit)
    {
      Cell_handle c = *chit;
      if(!c->has_vertex(vh)) // consider only the cells that have the edge opp-vh
        continue;

      c->info() = info;

      // grab the seed indices that need to be repainted
      for(int i=0; i<4; ++i)
      {
        Vertex_handle vhi = c->vertex(i);
        if(!tr().is_infinite(c->vertex(i)))
        {
          std::size_t seed_id = m_canvas->canvas_points[vhi->info()].closest_seed_id();
          if(seed_id != static_cast<std::size_t>(-1))
             seeds_to_reset.insert(seed_id);
        }
      }
    }
  }

  void set_info_on_incident_cells(Vertex_handle vh, int info)
  {
    std::vector<Cell_handle> incident_cells;
    tr().incident_cells(vh, std::back_inserter(incident_cells));

    Cell_handle_handle chit = incident_cells.begin();
    for(; chit!=incident_cells.end(); ++chit)
    {
      Cell_handle c = *chit;
      CGAL_assertion(!tr().is_infinite(c));

      int index_of_vh_in_c = -1;
      c->has_vertex(vh, index_of_vh_in_c);
      CGAL_postcondition(index_of_vh_in_c != -1);

      // mark the facet for refinement
      CGAL_precondition(info >= -1);
      extra_info.push_back(boost::make_tuple(c, index_of_vh_in_c, info));
      c->info() = -(extra_info.size() + 1); // so that the new info is in ]-inf, -2]
      cells_with_a_face_to_divide.push_back(c);

      // grab the seed indices whose Voronoi cell needs to be repainted
      for(int i=0; i<4; ++i)
      {
        Vertex_handle vhi = c->vertex(i);
        if(!tr().is_infinite(c->vertex(i)))
        {
          std::size_t seed_id = m_canvas->canvas_points[vhi->info()].closest_seed_id();
          if(seed_id != static_cast<std::size_t>(-1))
             seeds_to_reset.insert(seed_id);
        }
      }
    }
  }

  int get_common_facet_id(const Cell_handle c1, const Cell_handle c2)
  {
    // id of the common face in c2
    for(int i=0; i<4; ++i)
      if(c2->neighbor(i) == c1)
        return i;
    CGAL_assertion(false);
    return -1;
  }

  void subdivide_facet(Cell_handle& c)
  {
    // very important tricky trick here: once a facet has be refined from one
    // side, the cells incident to the new facets have recovered their correct
    // info() and the test below makes sure we don't do a double refinement
    // (once from each side).
    if(c->info() >= -1)
      return;

    CGAL_assertion(c->info() < -1 &&
                   c->info() >= static_cast<int>(-(extra_info.size()+1)));
    int real_index = -(c->info() + 2);
    const boost::tuple<Cell_handle, int, int>& tup = extra_info[real_index];
    CGAL_assertion(c == tup.template get<0>());

    int i = tup.template get<1>(); // index of the facet that must be refined

    if(i < 0)
      return;
    CGAL_assertion(i >= 0 && i < 4);
    Cell_handle c_mirror = c->neighbor(i);

    // At the moment, c_mirror->info() lives in ]inf,inf[ but since we're going to
    // refine the common face, we can re assign the correct info() to both cells
    int& c_mirror_info = c_mirror->info();
    if(c_mirror_info <= -2)
    {
      // just to be safe:
      // if c_mirror was also marked for refinement, it has to be the common face
      const boost::tuple<Cell_handle, int, int>& tup_mirror =
                                               extra_info[-(c_mirror_info + 2)];
      CGAL_assertion(tup_mirror.template get<0>() == c_mirror &&
                     c_mirror->neighbor(tup_mirror.template get<1>()) == c);
      c_mirror_info = tup_mirror.template get<2>();
    }

    // do the same for c
    c->info() = tup.template get<2>();

    // create the barycentric vertex
    const Point_3& p = facet_barycenter(c, i);
    Vertex_handle vh = tr().tds().insert_in_facet(c, i);
    vh->set_point(p);

    int point_id = m_canvas->canvas_points.size();
    vh->info() = point_id;
    typename Canvas::Canvas_point cp(p, point_id, m_canvas->mf, vh, &(tr()), m_canvas);
    m_canvas->canvas_points.push_back(cp);

    // set the info of the cells
    int f_in_c_mirror = get_common_facet_id(c, c_mirror);
    CGAL_assertion(c_mirror->neighbor(f_in_c_mirror) == c);
    Vertex_handle opposite_in_c = c->vertex(i);
    Vertex_handle opposite_in_c_mirror = c_mirror->vertex(f_in_c_mirror);

    set_info_on_edge_incident_cells(opposite_in_c, vh, c->info());
    set_info_on_edge_incident_cells(opposite_in_c_mirror, vh, c_mirror_info);
  }

  void subdivide_cell(Cell_handle c)
  {
    CGAL_precondition(!tr().is_infinite(c));
//    std::cout << "split cell : ";
//    for(int i=0; i<4; ++i)
//      std::cout << c->vertex(i)->info() << " ("
//                << m_canvas->canvas_points[c->vertex(i)->info()].closest_seed_id() << ") ";
//    std::cout << std::endl;

    int info = c->info();
    CGAL_precondition(info >= -1);

    // compute mid point
    const Point_3& p = cell_barycenter(c);

    // create vertex
    Vertex_handle vh = tr().tds().insert_in_cell(c);
    vh->set_point(p);

    int point_id = m_canvas->canvas_points.size();
    vh->info() = point_id;
    typename Canvas::Canvas_point cp(p, point_id, m_canvas->mf, vh, &(tr()), m_canvas);
    m_canvas->canvas_points.push_back(cp);

    // set the info of the cells and mark the facets that need to be split
    set_info_on_incident_cells(vh, info);
  }

  void paint_new_canvas()
  {
#if (verbosity > 20)
    std::cout << "must reset everything with color : ";
    std::set<int>::iterator strit = seeds_to_reset.begin();
    for(; strit!=seeds_to_reset.end(); ++strit)
      std::cout << *strit << " ";
    std::cout << std::endl;
#endif

    // brute force ish for now : reset some of the seeds
    for(std::size_t i=0; i<m_canvas->canvas_points.size(); ++i)
    {
      Canvas_point& cp = m_canvas->canvas_points[i];
      cp.invalidate_caches();

      int seed_id = cp.closest_seed_id();
      if(true || seeds_to_reset.find(seed_id) != seeds_to_reset.end())
        cp.reset_paint();
    }

    CGAL_assertion(!seeds_to_reset.empty());
    std::set<int>::iterator sit = seeds_to_reset.begin();
    for(; sit!=seeds_to_reset.end(); ++sit)
      m_canvas->locate_and_initialize(m_canvas->seeds[*sit], *sit);

    m_canvas->reset_counters();
    m_canvas->clear_primal();
    m_canvas->paint();
  }

  void post_subdivide()
  {
    CGAL_postcondition(tr().is_valid(true));
    CGAL_postcondition(m_canvas->canvas_points.size() == tr().number_of_vertices());

    // check the correspondence between canvas points & triangulation points
    for(std::size_t i=0; i<m_canvas->canvas_points.size(); ++i)
    {
      const Canvas_point& cp = m_canvas->canvas_points[i];
      CGAL_postcondition(cp.point() == cp.m_v->point());
    }

CGAL_expensive_assertion_code(
    // make sure everything is reachable in the new canvas
    std::vector<int> visited(m_canvas->canvas_points.size(), 0);
    // 0 = not visited, 1 = in queue, 2 = visited

    std::list<int> to_be_visited;
    std::size_t visited_counter = 0;
    to_be_visited.push_back(0);
    while(!to_be_visited.empty())
    {
      const Canvas_point& cp = m_canvas->canvas_points[to_be_visited.front()];
      cp.invalidate_caches();
      visited[cp.index()] = 2;
      visited_counter++;
      to_be_visited.pop_front();

      typename std::vector<Vertex_handle>::iterator it =
                                    cp.finite_interior_adjacent_vertices_begin();
      typename std::vector<Vertex_handle>::iterator end =
                                      cp.finite_interior_adjacent_vertices_end();
      for(; it!=end; ++it)
      {
        const Canvas_point& cq = m_canvas->canvas_points[(*it)->info()];
        if(visited[cq.index()] == 0)
        {
          visited[cq.index()] = 1;
          to_be_visited.push_back(cq.index());
        }
      }
    }
    std::cout << "total visited: " << visited_counter << " out of "
              << m_canvas->canvas_points.size() << std::endl;

    for(std::size_t i=0; i<m_canvas->canvas_points.size(); ++i)
    {
      CGAL_assertion(visited[i] == 2);
    }
); // CGAL expensive assertion

    // check that the cell info() are reasonable
    Cell_handle ch = tr().all_cells_begin();
    for(; ch!=tr().all_cells_end(); ++ch)
    {
      CGAL_assertion(ch->info() >= -1);
      if(tr().is_infinite(ch))
        continue;
      CGAL_assertion(ch->info() >= 0);
    }

    // only functions that are not debug code in this function
    paint_new_canvas();
    m_canvas->output_canvas("subdivided_canvas");

    for(std::size_t i=0; i<m_canvas->canvas_points.size(); ++i)
    {
      const Canvas_point& cp = m_canvas->canvas_points[i];
      CGAL_postcondition(cp.state() == KNOWN);
      CGAL_postcondition(cp.distance_to_closest_seed() != FT_inf);
      CGAL_postcondition(cp.closest_seed_id() != static_cast<std::size_t>(-1));
    }
  }

  std::size_t subdivide_cells()
  {
    std::size_t init_point_count = m_canvas->canvas_points.size();
#if (verbosity > 15)
    std::cout << "split cells" << std::endl;
    std::cout << "before split: " << init_point_count << " vertices" << std::endl;
#endif

#if (verbosity > 20)
    std::cout << "the ids to subdivide are: ";
    std::set<int>::iterator it = m_criterion.seed_indices.begin();
    for(; it!=m_criterion.seed_indices.end(); ++it)
      std::cout << *it << " ";
    std::cout << std::endl;
#endif

    // since the triangulation gets modified and i'm not sure if the underlying
    // container is changing the order of cells or not, i'm copying the handle
    // so I am sure I only split at most once any (input) cell
    std::vector<Cell_handle> finite_cells;
    Cell_handle cit = tr().finite_cells_begin();
    Cell_handle cend = tr().finite_cells_end();
    for(; cit!=cend; ++cit)
    {
      if(cit->info() > 0)
        finite_cells.push_back(cit);
    }

    for(std::size_t i=0; i<finite_cells.size(); ++i)
    {
      Cell_handle cit = finite_cells[i];
      if(should_cell_be_subdivided(cit))
        subdivide_cell(cit);
    }

    // now we subdivide every facet that is incident to at least one cell that was split
    CGAL_assertion(cells_with_a_face_to_divide.size() == extra_info.size());
    typename std::list<Cell_handle>::iterator cwit = cells_with_a_face_to_divide.begin();
    for(; cwit!=cells_with_a_face_to_divide.end(); ++cwit)
    {
      Cell_handle c = *cwit;
      subdivide_facet(c);
    }

    std::size_t final_point_count = m_canvas->canvas_points.size();
#if (verbosity > 15)
    std::cout << "after split: " << final_point_count << " vertices" << std::endl;
#endif
    CGAL_assertion(final_point_count >= init_point_count);
    return final_point_count - init_point_count;
  }

  void subdivide()
  {
    std::size_t new_points_count = subdivide_cells();
    if(new_points_count == 0)
      return;
    post_subdivide();
  }

  Canvas_subdivider(Canvas* canvas, const Criterion& criterion)
    :
      m_canvas(canvas),
      m_criterion(criterion),
      seeds_to_reset(),
      extra_info(),
      cells_with_a_face_to_divide()
  { }
};

} // namespace Anisotropic_mesh_3
} // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_Canvas_subdivider_H

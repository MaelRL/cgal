#ifndef CGAL_ANISOTROPIC_MESH_3_CAMPEN_C3T3_CANVAS_H
#define CGAL_ANISOTROPIC_MESH_3_CAMPEN_C3T3_CANVAS_H

#include <CGAL/Canvas/canvas_config.h>
#include <CGAL/Canvas/canvas.h>
#include <CGAL/Canvas/Campen_c3t3_point.h>
#include <CGAL/Canvas/canvas_triangulation_io.h>
#include <CGAL/Canvas/canvas_subdiviser.h>

#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/assertions.h>

#include <boost/unordered_set.hpp>

#include <cstddef>
#include <deque>
#include <istream>
#include <iterator>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <utility>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K, typename Metric_field>
class Campen_canvas :
    public Canvas<K, Campen_canvas_point<K, Metric_field>, Metric_field>
{
private:
  typedef Campen_canvas<K, Metric_field>                          Self;

public:
  typedef Campen_canvas_point<K, Metric_field>                    Canvas_point;
  typedef Canvas_point*                                           Canvas_point_handle;

  typedef int                          Vertex_Info; // index of the canvas point
  typedef int                          Cell_Info; // index of the subdomain

  typedef Canvas<K, Canvas_point, Metric_field>                   Base;

  typedef typename Base::FT                                       FT;
  typedef typename Base::Point_3                                  Point_3;
  typedef typename Base::Tetrahedron_3                            Tetrahedron_3;
  typedef typename Base::Metric                                   Metric;
  typedef typename Base::Vector3d                                 Vector3d;
  typedef typename Base::Simplex                                  Simplex;
  typedef typename Base::BEdge                                    BEdge;
  typedef typename Base::BTriangle                                BTriangle;
  typedef typename Base::BTetrahedron                             BTetrahedron;

  typedef CGAL::Triangulation_vertex_base_with_info_3<Vertex_Info, K>       Vb;
  typedef CGAL::Triangulation_cell_base_with_info_3<Cell_Info, K>           Cb;
  typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                      TDS;
  typedef CGAL::Triangulation_3<K, TDS>                                     Tr;

  typedef typename Tr::Vertex_handle                       Vertex_handle;
  typedef std::vector<Vertex_handle>                       Vertex_handle_vector;
  typedef typename Vertex_handle_vector::iterator          Vertex_handle_handle;
  typedef typename Tr::Finite_vertices_iterator            Finite_vertices_iterator;
  typedef typename Tr::Cell_handle                         Cell_handle;
  typedef std::vector<Cell_handle>                         Cell_handle_vector;
  typedef typename Cell_handle_vector::iterator            Cell_handle_handle;
  typedef typename Tr::Finite_cells_iterator               Finite_cells_iterator;

  Tr m_tr;

  std::size_t compute_precise_Voronoi_vertex_on_edge(const std::size_t p,
                                                     const std::size_t q,
                                                     int call_count = 0)
  {
    const Canvas_point& cp = this->get_point(p);
    const Canvas_point& cq = this->get_point(q);

    // dichotomy with the centroid is most likely not the optimal thing...

    // create a virtual canvas point
    Point_3 centroid = CGAL::barycenter(cp.point(), 0.5, cq.point(), 0.5);
    this->canvas_points.push_back(Canvas_point(centroid, this->canvas_points.size(),
                                               NULL/*vh*/, NULL/*tr*/, this));
    Canvas_point& virtual_cp = this->canvas_points.back();

    if(cp.closest_seed_id() == cq.closest_seed_id())
    {
      CGAL_assertion(call_count == 0);
      return virtual_cp.index();
    }

    // set distance, seed_id and ancestor for the virtual point
    virtual_cp.compute_closest_seed(cp);
    virtual_cp.compute_closest_seed(cq);
    virtual_cp.state() = KNOWN;

    if(call_count > 0) // fixme don't hardcode stuff
      return virtual_cp.index();

    if(virtual_cp.closest_seed_id() == cp.closest_seed_id())
      return compute_precise_Voronoi_vertex_on_edge(virtual_cp.index(), q, ++call_count);
    else
    {
      CGAL_assertion(virtual_cp.closest_seed_id() == cq.closest_seed_id());
      return compute_precise_Voronoi_vertex_on_edge(virtual_cp.index(), p, ++call_count);
    }
  }

  std::size_t compute_precise_Voronoi_vertex_in_tetrahedron(const std::size_t p,
                                                            const std::size_t q,
                                                            const std::size_t r,
                                                            const std::size_t s,
                                                            int call_count = 0)
  {
    const Canvas_point& cp = this->get_point(p);
    const Canvas_point& cq = this->get_point(q);
    const Canvas_point& cr = this->get_point(r);
    const Canvas_point& cs = this->get_point(s);

    // centroid is probably not the most optimal...
    Point_3 centroid = CGAL::centroid(cp.point(), cq.point(), cr.point(), cs.point());
    this->canvas_points.push_back(Canvas_point(centroid, this->canvas_points.size(),
                                               NULL/*vh*/, NULL/*tr*/, this));
    Canvas_point& virtual_cp = this->canvas_points.back();

    // set distance, seed_id and ancestor for the virtual point
    virtual_cp.compute_closest_seed(cp);
    virtual_cp.compute_closest_seed(cq);
    virtual_cp.compute_closest_seed(cr);
    virtual_cp.compute_closest_seed(cs);
    virtual_cp.state() = KNOWN;

    // another potential stop is if the (max) distance between the centroid
    // and a point of the triangle is below a given distance
    if(call_count > 0)
      return virtual_cp.index();

    // since we might have ci.closest_seed_id = cj.closest_seed_id, it kinda
    // favors the earliest point, but that shouldn't matter
    if(virtual_cp.closest_seed_id() == cp.closest_seed_id())
      return compute_precise_Voronoi_vertex_in_tetrahedron(virtual_cp.index(),
                                                           q, r, s,
                                                           ++call_count);
    else if(virtual_cp.closest_seed_id() == cq.closest_seed_id())
      return compute_precise_Voronoi_vertex_in_tetrahedron(p, virtual_cp.index(),
                                                           q, s, ++call_count);
    else if(virtual_cp.closest_seed_id() == cr.closest_seed_id())
      return compute_precise_Voronoi_vertex_in_tetrahedron(p, q, virtual_cp.index(),
                                                           s, ++call_count);
    else
    {
      CGAL_assertion(virtual_cp.closest_seed_id() == cs.closest_seed_id());
      return compute_precise_Voronoi_vertex_in_tetrahedron(p, q, r, virtual_cp.index(),
                                                           ++call_count);
    }
  }

  void initialize()
  {
#if (VERBOSITY > 5)
    std::cout << "canvas initialization" << std::endl;
#endif

    // read and create the canvas points
    std::ifstream in((this->canvas_str + ".mesh").c_str());
    if(!in)
    { 
      std::cerr << "couldn't read canvas input file : " << (this->canvas_str + ".mesh") << std::endl;
      exit(0);
    }
    else
      std::cout << "reading: " << (this->canvas_str + ".mesh") << std::endl;

    std::vector<bool> border_info_vec;
    bool is_tr_well_built = CGAL::build_triangulation_from_file(in, m_tr, border_info_vec);
    CGAL_postcondition(is_tr_well_built);
    CGAL_postcondition(border_info_vec.size() == m_tr.number_of_vertices());

    // build the canvas points
    std::size_t vertex_counter = 0;
    Finite_vertices_iterator vit = m_tr.finite_vertices_begin();
    for(; vit!=m_tr.finite_vertices_end(); ++vit)
    {
      Canvas_point cp(vit->point(), vertex_counter, vit, &m_tr, this);
      cp.border_info() = border_info_vec[vertex_counter++];
      this->canvas_points.push_back(cp);
      vit->info() = this->canvas_points.size() - 1;
      CGAL_postcondition(this->canvas_points[vit->info()].point() == vit->point());
    }
    CGAL_postcondition(this->canvas_points.size() == m_tr.number_of_vertices());

    // count the interior faces (purely debug information)
    std::size_t interior_cell_count = 0;
    Cell_handle cit = m_tr.all_cells_begin();
    for(; cit!=m_tr.all_cells_end(); ++cit)
    {
      if(cit->info() > 0)
        ++interior_cell_count;
    }
    std::cout << interior_cell_count << " interior cells" << std::endl;
    CGAL_postcondition(interior_cell_count > 0);

    Base::compute_canvas_bbox();
    this->seeds.initialize_seeds();
    Base::locate_seeds_on_canvas();

#if (VERBOSITY > 5)
    std::cout << "canvas initialized with " << this->seeds.size() << " seeds" << std::endl;
#endif
  }

  void initialize_vertex(Vertex_handle vh, const std::size_t seed_id,
                         const Point_3& t)
  {
    const Metric& seed_m = this->seeds.seeds_metrics[seed_id];

    CGAL_assertion(!m_tr.is_infinite(vh));
    Canvas_point& cp = this->canvas_points[vh->info()];
    const Metric& v_m = cp.metric();
    const Eigen::Matrix3d& f = get_interpolated_transformation(seed_m, v_m);

    Vector3d v;
    v(0) = t.x() - cp.point().x();
    v(1) = t.y() - cp.point().y();
    v(2) = t.z() - cp.point().z();
    v = f * v;
    FT d = v.norm(); // d = std::sqrt(v^t * M * v) = (f*v).norm()
    Base::initialize_canvas_point(cp, d, seed_id);
  }

  void initialize_seed_on_vertex(Cell_handle c, int li,
                                 std::size_t seed_id,
                                 const Point_3& t)
  {
    Vertex_handle vh = c->vertex(li);
    if(c->info() < 1)
    {
      // there has to exist an incident finite interior cell or the seed is rejected
      Canvas_point& cp = this->canvas_points[vh->info()];
      Cell_handle_handle b = cp.finite_interior_incident_cells_begin();
      Cell_handle_handle e = cp.finite_interior_incident_cells_end();
      std::ptrdiff_t d = e - b;

      if(d == 0)
      {
#if (VERBOSITY > 2)
        std::cout << "vertex " << c->vertex(li)->info()
                  << " is not acceptable for seed " << seed_id
                  << " 's initialization. (Seed point : " << t << ")" << std::endl;
#endif
#ifndef ANISO_GEO_FORCE_SEED_INITIALIZATION
        return;
#endif
      }
    }
    initialize_vertex(vh, seed_id, t);
  }

  void initialize_seed_on_edge(Cell_handle c, int li, int lj,
                               std::size_t seed_id,
                               const Point_3& t)
  {
    if(c->info() < 1)
    {
      // there has to be an incident finite interior cell or the seed is rejected
      typename Tr::Cell_circulator current_cell = m_tr.incident_cells(c, li, lj);
      typename Tr::Cell_circulator done = current_cell;
      bool found = false;
      do
      {
        if(current_cell->info() > 0)
        {
          found = true;
          break;
        }
        ++current_cell;
      } while (current_cell != done);
      if(!found)
      {
#if (VERBOSITY > 2)
        std::cout << "edge " << c->vertex(li)->info() << " "
                             << c->vertex(lj)->info()
                  << " is not acceptable for seed " << seed_id
                  << " 's initialization. (Seed point : " << t << ")" << std::endl;
#endif
#ifndef ANISO_GEO_FORCE_SEED_INITIALIZATION
        return;
#endif
      }
    }

    for(int i=0; i<4; ++i)
      if(i == li || i == lj)
        initialize_vertex(c->vertex(i), seed_id, t);
  }

  void initialize_seed_on_triangle(Cell_handle c, int li,
                                   std::size_t seed_id,
                                   const Point_3& t)
  {
    if(c->info() < 1)
    {
      // there has to be an incident finite interior cell or the seed is rejected
      Cell_handle c_mirror = c->neighbor(li);
      if(c_mirror->info() < 1)
      {
#if (VERBOSITY > 2)
        std::cout << "face " << c->vertex((li + 1) % 4)->info() << " "
                             << c->vertex((li + 2) % 4)->info() << " "
                             << c->vertex((li + 3) % 4)->info()
                  << " is not acceptable for seed " << seed_id
                  << " 's initialization. (Seed point : " << t << ")" << std::endl;
#endif
#ifndef ANISO_GEO_FORCE_SEED_INITIALIZATION
        return;
#endif
      }
    }

    for(int i=0; i<4; ++i)
      if(i != li)
        initialize_vertex(c->vertex(i), seed_id, t);
    return;
  }

  void initialize_seed_in_cell(Cell_handle c,
                               std::size_t seed_id,
                               const Point_3& t)
  {
    // has to be an interior cell (this checks for infinite cell)
    if(c->info() < 1)
    {
#if (VERBOSITY > 2)
        std::cout << "cell " << c->vertex(0)->info() << " "
                             << c->vertex(1)->info() << " "
                             << c->vertex(2)->info() << " "
                             << c->vertex(3)->info()
                  << " is not acceptable for seed " << seed_id
                  << " 's initialization. (Seed point : " << t << ")" << std::endl;
#endif
#ifndef ANISO_GEO_FORCE_SEED_INITIALIZATION
      return;
#endif
    }

    for(std::size_t j=0; j<4; ++j)
      initialize_vertex(c->vertex(j), seed_id, t);
    return;
  }

  void locate_and_initialize(const Point_3& t, const std::size_t seed_id)
  {
    // find the tetrahedron that contains the seed and initialize the distance at
    // these four points accordingly.

#if (VERBOSITY > 10)
        std::cout << "locating seed " << seed_id << " point: " << t  << std::endl;
#endif
    typename Tr::Locate_type lt;
    int li, lj;
    Cell_handle c = m_tr.locate(t, lt, li, lj);

    if(m_tr.is_infinite(c))
      std::cout << "WARNING: seed located in an infinite cell" << std::endl;

#if (VERBOSITY > 15)
    std::cout << "locate in triangulation : " << std::endl;
    std::cout << "locate type: " << lt << " li/j: " << li << " " << lj << std::endl;
    std::cout << "cell : " << &*c << " with info : " << c->info() << std::endl;
    for(int i=0; i<4; ++i)
    {
      Vertex_handle vh = c->vertex(i);
      CGAL_assertion(!m_tr.is_infinite(vh));
      std::cout << "point: " << c->vertex(i)->point() << std::endl;
      std::cout << vh->info() << " || " << this->canvas_points[vh->info()].point() << std::endl;
    }
#endif

    typename Tr::Tetrahedron tetra(c->vertex(0)->point(), c->vertex(1)->point(),
                                   c->vertex(2)->point(), c->vertex(3)->point());
    CGAL_assertion( m_tr.geom_traits().has_on_bounded_side_3_object()(tetra, t) ||
                    m_tr.geom_traits().has_on_boundary_3_object()(tetra, t));
    CGAL_assertion(lt < 4);

      // a degenerate case that can happen : the seed is on the border of the
      // interior but the cell found by locate() is exterior.

      // we accept a partial initialization of the cell (lt=0 -> the vertex,
      // lt=1, the vertices of the edge, etc.) ONLY if there exists an incident
      // interior finite cell
    if(lt == 0) // vertex
      initialize_seed_on_vertex(c, li, seed_id, t);
    else if(lt == 1) // edge
      initialize_seed_on_edge(c, li, lj, seed_id, t);
    else if(lt == 2) // facet
      initialize_seed_on_triangle(c, li, seed_id, t);
    else
      initialize_seed_in_cell(c, seed_id, t);
  }

  void primal_shenanigans(const Canvas_point* cp)
  {
#if (VERBOSITY > 15)
    std::cout << "primal shenanigans at : " << cp->index() << std::endl;
#endif

    CGAL_assertion(cp->state() == KNOWN);
    typedef boost::unordered_set<const Canvas_point*>   Candidates_set;

    Cell_handle_handle cit = cp->finite_interior_incident_cells_begin();
    Cell_handle_handle cend = cp->finite_interior_incident_cells_end();

#if (VERBOSITY > 25)
    std::cout << "cp: " << cp->index() << " (" << cp->point() << ") ";
    std::cout << "has " << cend - cit << " incident interior cells" << std::endl;
#endif
    for(; cit!=cend; ++cit)
    {
      const Cell_handle c = *cit;
      CGAL_assertion(c->info() > 0);

      Candidates_set candidates;
      candidates.insert(cp);

      for(std::size_t j=0; j<4; ++j)
      {
        const Vertex_handle vh = c->vertex(j);
        CGAL_assertion(!m_tr.is_infinite(vh));
        const Canvas_point& cq = this->canvas_points[vh->info()];
        if(cq.state() != KNOWN)
          continue;

        candidates.insert(&cq);
//        fill_edge_bisectors_map(cp, cq); // pseudo bisectors
      }

      Base::construct_primal_elements_from_candidates(candidates);
      //todo mark/compute Voronoi vertices
    }
  }

  void compute_local_primal_elements(const Canvas_point* cp)
  {
    primal_shenanigans(cp);
  }

  void check_canvas_density()
  {
    // verify that the canvas is dense enough for this new point to have a proper
    // Voronoi cell. If it's not the case, refine the canvas.

    typedef CGAL::Anisotropic_mesh_3::Subdomain_criterion<Self> SCriterion;
    SCriterion criterion(this);
    criterion.seed_indices = this->check_canvas_density_with_primals();
    CGAL::Anisotropic_mesh_3::Canvas_subdivider<Self, SCriterion> canvas_sub(this, criterion);
    canvas_sub.subdivide();
  }

  void compute_primal()
  {
    if(this->primal_edges.empty() && this->primal_triangles.empty() && this->primal_tetrahedra.empty())
    {
#if (VERBOSITY > 5)
    std::cout << "Primal computations" << std::endl;
#endif
      typedef boost::unordered_set<const Canvas_point*>   Candidates_set;
      Finite_cells_iterator cit = m_tr.finite_cells_begin();
      Finite_cells_iterator end = m_tr.finite_cells_end();
      for(; cit!=end; ++cit)
      {
        // ignore exterior cells to ignore triangulation neighbors that are not canvas neighbors
        if(cit->info() < 1)
          continue;

        std::size_t color_0 = this->canvas_points[cit->vertex(0)->info()].closest_seed_id();
        std::size_t color_1 = this->canvas_points[cit->vertex(1)->info()].closest_seed_id();
        std::size_t color_2 = this->canvas_points[cit->vertex(2)->info()].closest_seed_id();
        std::size_t color_3 = this->canvas_points[cit->vertex(3)->info()].closest_seed_id();

        // filter tets with only one color
        if(color_0 == color_1 && color_1 == color_2 && color_2 == color_3)
          continue;

#ifndef COMPUTE_PRIMAL_ALL_DIMENSIONS
        if(color_0 != color_1 && color_0 != color_2 && color_0 != color_3 &&
           color_1 != color_2 && color_1 != color_3 &&
           color_2 != color_3)
        { /* don't filter four different colors */ }
        else // below is either 2 or 3 colors
        {
          // filter tets with two or three different colors but none of the vertices
          // are on the border of the domain
          if(!this->canvas_points[cit->vertex(0)->info()].border_info() &&
             !this->canvas_points[cit->vertex(1)->info()].border_info() &&
             !this->canvas_points[cit->vertex(2)->info()].border_info() &&
             !this->canvas_points[cit->vertex(3)->info()].border_info())
            continue;
        }
#endif

        Candidates_set candidates;
        for(std::size_t j=0; j<4; ++j)
          candidates.insert(&(this->canvas_points[cit->vertex(j)->info()]));

        Base::construct_primal_elements_from_candidates(candidates);
      }
      // todo add mark/compute Voronoi vertices etc.
    }
    else
      std::cout << "WARNING: call to compute_primal with non-empty primal data structures..." << std::endl;

#ifdef CGAL_ANISO_REFINE_LOW_CANVAS_DENSITY
    check_canvas_density();
#else
    this->check_canvas_density_with_primals();
#endif
  }

  void check_edelsbrunner()
  {
    // todo: every structure here could be unordered_
    // the algorithm is very straightfoward and explicit rather than elegant

    std::cout << "it's Edel's time" << std::endl;

    // for each simplex, associate all the tetrahedra (vector of size_t, the size_t
    // is the index of the tet in 'tetrahedra') of the base canvas that correspond
    // to that simplex...
    // Then check that this tet set satisfies the closed ball property
    typedef std::map<Simplex, Cell_handle_vector> Dual_map;

    Dual_map dual_map;
    int failed_blob_counter = 0;

    Finite_cells_iterator cit = m_tr.finite_cells_begin();
    Finite_cells_iterator cend = m_tr.finite_cells_end();
    for(; cit!=cend; ++cit)
    {
      if(cit->info() < 1)
        continue;

      Simplex primal_simplex; // has to be a set since it is used as a key in a map
      for(std::size_t j=0; j<4; ++j)
        primal_simplex.insert(cit->vertex(j)->closest_seed_id);

      std::pair<typename Dual_map::iterator, bool> is_insert_successful;
      Cell_handle_vector vec;
      vec.push_back(cit);
      is_insert_successful = dual_map.insert(std::make_pair(primal_simplex, vec));

      if(!is_insert_successful.second) // primal_simplex already has an entry in the map
        (is_insert_successful.first)->second.push_back(cit);
    }

    // now, for each simplex, we want the tetrahedra corresponding to that simplex
    // to form a single blob

    typename Dual_map::const_iterator dmit = dual_map.begin();
    typename Dual_map::const_iterator dmend = dual_map.end();
    for(; dmit!=dmend; ++dmit)
    {
      const Simplex& simplex = dmit->first;
      std::cout << "looking at the simplex : ";
      typename Simplex::iterator sit = simplex.begin();
      typename Simplex::iterator send = simplex.end();
      for(; sit!=send; ++sit)
        std::cout << *sit << " ";
      std::cout << std::endl;

      Cell_handle_vector& tetrahedra_ids = dmit->second;

      // collect all the triangles that appear once
      std::set<BTriangle> triangles;
      for(std::size_t i=0; i<tetrahedra_ids.size(); ++i)
      {
        const Cell_handle c = tetrahedra_ids[i];
        for(int j=0; j<4; ++j)
        {
          BTriangle triangle;
          triangle[0] = c->vertex((j+1)%4);
          triangle[1] = c->vertex((j+2)%4);
          triangle[2] = c->vertex((j+3)%4);
          // gotta sort since we'll use a map with edge keys
          std::sort(triangle.begin(), triangle.end());

          std::pair<typename std::set<BTriangle>::iterator, bool> is_insert_successful;
          is_insert_successful = triangles.insert(triangle);

          // already in 'triangles', thus triangle is internal & ignore it
          // it's okay to delete now because at most a triangle appears twice!!
          if(!is_insert_successful.second)
            triangles.erase(is_insert_successful.first);
        }
      }

      std::cout << triangles.size() << " triangles" << std::endl;
      CGAL_assertion(!triangles.empty());

      // check that all the border triangles of the blob actually form a blob
      // note: 'bow-ties' configurations could appear... what to do ?
      std::map<const BTriangle*, bool> visited_status;
      std::map<BEdge, std::deque<const BTriangle*> > incident_tris; // edge to tri

      typename std::set<BTriangle>::iterator it = triangles.begin();
      typename std::set<BTriangle>::iterator iend = triangles.end();
      for(; it!=iend; ++it)
      {
        const BTriangle& tr = *it;

        std::pair<typename  std::map<BEdge, std::deque<const BTriangle*> >::iterator, bool>
                                                           is_insert_successful;

        std::deque<const BTriangle*> vec;
        vec.push_back(&tr);

        BEdge e;
        e[0] = tr[0]; e[1] = tr[1];
        is_insert_successful = incident_tris.insert(std::make_pair(e, vec));
        if(!is_insert_successful.second)
          is_insert_successful.first->second.push_back(&tr);

        e[0] = tr[0]; e[1] = tr[2];
        is_insert_successful = incident_tris.insert(std::make_pair(e, vec));
        if(!is_insert_successful.second)
          is_insert_successful.first->second.push_back(&tr);

        e[0] = tr[1]; e[1] = tr[2];
        is_insert_successful = incident_tris.insert(std::make_pair(e, vec));
        if(!is_insert_successful.second)
          is_insert_successful.first->second.push_back(&tr);

        visited_status[&tr] = false;
      }

      CGAL_assertion(!incident_tris.empty());

      std::deque<const BTriangle*> triangles_to_visit;
      triangles_to_visit.push_back(visited_status.begin()->first);
      while(!triangles_to_visit.empty())
      {
//        std::cout << "triangles to visit: " << triangles_to_visit.size() << std::endl;
        const BTriangle& current_tri = *(triangles_to_visit.front());
        triangles_to_visit.pop_front();

        visited_status[&current_tri] = true;

        // add all the triangles incident to the edges of current tri
        for(int i=0; i<3; ++i)
        {
          BEdge e;
          if(i==2) // since current_tri is sorted and we need to find e in the map
          {
            e[0] = current_tri[0];
            e[1] = current_tri[2];
          }
          else
          {
            e[0] = current_tri[i];
            e[1] = current_tri[i+1];
          }

          typename std::deque<const BTriangle*>::iterator dit = incident_tris[e].begin();
          for(; dit!=incident_tris[e].begin(); ++dit)
          {
            if(!visited_status[*dit])
              triangles_to_visit.push_back(*dit);
          }
        }
      }

      typename std::map<const BTriangle*, bool>::const_iterator mit =
                                                         visited_status.begin();
      for(; mit!=visited_status.end(); ++mit)
      {
        if(!mit->second)
        {
          std::cout << "Blob not achieved" << std::endl;
          ++failed_blob_counter;
          break;
        }
      }
    }
    std::cout << "end edelsbrunner" << std::endl;
    std::cout << failed_blob_counter << " out of " << dual_map.size()
              << " failed" << std::endl;
  }

  void output_canvas(const std::string str_base) const
  {
#define ONLY_PRINT_BISECTORS
    std::ofstream out((str_base + ".mesh").c_str());
    std::ofstream out_bb((str_base + ".bb").c_str());

    out << "MeshVersionFormatted 1" << std::endl;
    out << "Dimension 3" << std::endl;
    out << "Vertices" << std::endl;
    out << this->canvas_points.size() << std::endl;
    out_bb << "3 1 " << this->canvas_points.size() << " 2" << std::endl;

    for(std::size_t i=0; i<this->canvas_points.size(); ++i)
    {
      const Canvas_point& cp = this->canvas_points[i];
      out << cp.point() << " " << i+1 << std::endl;

#ifdef ONLY_PRINT_BISECTORS
      out_bb << cp.closest_seed_id() << std::endl;
#else
      out_bb << cp.distance_to_closest_seed() << std::endl;
#endif
    }

    /*
    out << "Triangles" << std::endl;
    out << m_tr.number_of_finite_facets() << std::endl;
    typename Tr::Finite_facets_iterator fit = m_tr.finite_facets_begin();
    typename Tr::Finite_facets_iterator fend = m_tr.finite_facets_end();
    for(; fit!=fend; ++fit)
    {
      Cell_handle c = fit->first;
      Cell_handle c_opp = c->neighbor(fit->second);
      if(c->info() < 1 && c_opp->info() < 1) // ignore purely exterior facets
        continue;

      for(std::size_t j=1; j<=3; ++j)
        out << c->vertex((fit->second+j)%4)->info() + 1 << " ";
      out << "1" << std::endl;
    }
    */

    out << "Tetrahedra" << std::endl;
    out << m_tr.number_of_finite_cells() << std::endl;
    Finite_cells_iterator cit = m_tr.finite_cells_begin();
    Finite_cells_iterator cend = m_tr.finite_cells_end();
    for(; cit!=cend; ++cit)
    {
      if(cit->info() < 1) // ignore exterior cells
        continue;

      boost::unordered_set<std::size_t> materials;
      for(std::size_t j=0; j<4; ++j)
        materials.insert(this->canvas_points[cit->vertex(j)->info()].closest_seed_id());

#ifdef ONLY_PRINT_BISECTORS
      if(materials.size() == 1)
        continue;
#endif
      for(std::size_t j=0; j<4; ++j)
        out << cit->vertex(j)->info() + 1 << " ";
      std::size_t mat = (materials.size()==1) ? (*(materials.begin())) :
                                                (this->seeds.size());
      out << mat << std::endl;
    }

    out_bb << "End" << std::endl;
    out << "End" << std::endl;
  }

  Campen_canvas(const std::string& canvas_str_,
                const std::string& seeds_str_,
                const std::size_t max_seeds_n_,
                const Metric_field* mf_)
    :
      Base(canvas_str_, seeds_str_, max_seeds_n_, mf_),
      m_tr()
  { }
};

} // namespace Anisotropic_mesh_3
} // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_CAMPEN_C3T3_CANVAS_H

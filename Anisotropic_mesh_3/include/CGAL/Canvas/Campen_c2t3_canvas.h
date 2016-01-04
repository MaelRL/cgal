#ifndef CGAL_ANISOTROPIC_MESH_3_CAMPEN_C2T3_CANVAS_H
#define CGAL_ANISOTROPIC_MESH_3_CAMPEN_C2T3_CANVAS_H

#include <CGAL/Canvas/canvas_config.h>
#include <CGAL/Canvas/canvas.h>
#include <CGAL/Canvas/Campen_c2t3_point.h>
#include <CGAL/Canvas/canvas_triangulation_io.h>
#include <CGAL/Canvas/canvas_subdiviser.h>

#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/assertions.h>

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <Eigen/Dense>

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

template<typename K, typename Metric_field, typename C3t3>
class Campen_canvas :
    public Canvas<K, Campen_canvas_point<K, Metric_field, C3t3>, Metric_field>
{
private:
  typedef Campen_canvas<K, Metric_field, C3t3>             Self;

public:
  typedef Campen_canvas_point<K, Metric_field, C3t3>       Canvas_point;
  typedef Canvas_point*                                    Canvas_point_handle;

  typedef int                          Vertex_Info; // index of the canvas point
  typedef int                          Cell_Info; // index of the subdomain

  typedef Canvas<K, Canvas_point, Metric_field>            Base;


  typedef typename Base::Metric                            Metric;
  typedef typename Base::Vector3d                          Vector3d;
  typedef typename Base::Simplex                           Simplex;
  typedef typename Base::BEdge                             BEdge;
  typedef typename Base::BTriangle                         BTriangle;
  typedef typename Base::BTetrahedron                      BTetrahedron;

  typedef typename C3t3::Triangulation                     Tr;
  typedef typename Tr::Geom_traits                         Gt;

  typedef typename Gt::FT                                  FT;
  typedef typename Gt::Bare_point                          Point_3;
  typedef typename Gt::Weighted_point                      Weighted_point_3;
  typedef typename Gt::Vector_3                            Vector_3;

  typedef typename C3t3::Vertices_in_complex_iterator      Vertices_in_complex_iterator;
  typedef typename C3t3::Facets_in_complex_iterator        Facet_iterator;

  typedef typename Tr::Vertex_handle                       Vertex_handle;
  typedef std::vector<Vertex_handle>                       Vertex_handle_vector;
  typedef typename Vertex_handle_vector::iterator          Vertex_handle_handle;
  typedef typename Tr::Facet                               Facet;
  typedef std::vector<Facet>                               Facet_vector;
  typedef typename Facet_vector::iterator                  Facet_handle;
  typedef typename Tr::Cell_handle                         Cell_handle;
  typedef typename Tr::Finite_vertices_iterator            Finite_vertices_iterator;

  typedef std::pair<std::size_t, std::size_t>                           Oriented_edge;
  typedef boost::unordered_map<Oriented_edge, FT>                       Angle_map;
  typedef typename boost::unordered_map<Oriented_edge, FT>::iterator    Angle_map_iterator;

  C3t3& m_c3t3;
  Angle_map m_angles;

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
                                               NULL/*vh*/, NULL/*c3t3*/, this));
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

  FT angle_in_radian(const Vector_3& u,
                     const Vector_3& v,
                     K k = K()) const
  {
    typename K::Compute_scalar_product_3 scalar_product =
      k.compute_scalar_product_3_object();
    typename K::Construct_cross_product_vector_3 cross_product =
      k.construct_cross_product_vector_3_object();
    typename K::Compute_squared_length_3 sq_length =
      k.compute_squared_length_3_object();

    // -------------------------------------
    // Angle between two vectors (in rad)
    // uv = |u||v| cos(u,v)
    // u^v  = w
    // |w| = |u||v| |sin(u,v)|
    // -------------------------------------
    FT product = CGAL::sqrt(sq_length(u) * sq_length(v));

    // Check
    if ( product == FT(0) )
      return FT(0);

    // Sine
    Vector_3 w = cross_product(u,v);
    FT abs_sin = CGAL::sqrt(sq_length(w)) / product;

    if ( abs_sin > FT(1) ) { abs_sin = FT(1); }
    CGAL_assertion(abs_sin >= 0.);
    CGAL_assertion(abs_sin <= 1);

    // We just need cosine sign
    FT cosine_sign = scalar_product(u,v);

    if ( cosine_sign >= FT(0) )
      return FT(std::asin(abs_sin));
    else
      return FT(CGAL_PI) - FT(std::asin(abs_sin));
  }

  FT compute_loc_angle(Finite_vertices_iterator vh0,
                       Vertex_handle vh1, Vertex_handle vh2) const
  {
    // compute the angle between v0v1 and v1v2
    Vector_3 v1(vh0->point(), vh1->point());
    Vector_3 v2(vh0->point(), vh2->point());

    return angle_in_radian(v1, v2);
  }

  void initialize_cumulative_angles()
  {
    Finite_vertices_iterator vit = m_c3t3.triangulation().finite_vertices_begin();
    Finite_vertices_iterator vend = m_c3t3.triangulation().finite_vertices_end();
    for(; vit!=vend; ++vit)
    {
      if(vit->in_dimension() > 2)
        continue; // only consider vertices on the boundary of the complex

      std::size_t id = vit->info();
      const Canvas_point& cp = this->canvas_points[id];

      FT angle = 0.;
      std::list<Angle_map_iterator> entries; // keep the new entries to normalize the angles at the end

      // compute the cumulative angles in a star (the initial one is 0)
      // the entry in the map is of the form :
      // [center C of the star, vertex V on the first ring] = angle V, C, Vnext

      // the adjacent vertices have been ordered when building the cache
      Vertex_handle_handle vhit = cp.adjacent_vertices_in_complex_begin();
      Vertex_handle_handle vhend = cp.adjacent_vertices_in_complex_end();
      for(; vhit!=vhend; ++vhit)
      {
        Vertex_handle_handle vhit_next = vhit;
        if(vhit == --(cp.adjacent_vertices_in_complex_end()))
          vhit_next = cp.adjacent_vertices_in_complex_begin();
        else
          std::advance(vhit_next, 1);

        Vertex_handle vh = *vhit;
        Vertex_handle vh_next = *vhit_next;

        // set up the angle for vh, vit, vh_next in the angle memory
        std::pair<Angle_map_iterator, bool> is_insert_successful =
                      m_angles.insert(std::make_pair(std::make_pair(vit->info(),
                                                                    vh->info()),
                                                     angle));
        CGAL_postcondition(is_insert_successful.second);
        entries.push_back(is_insert_successful.first);

        FT loc_angle = compute_loc_angle(vit, vh, vh_next);
        angle += loc_angle;
      }

      // normalize the angles (sum must be 2*pi)
      FT norm_coeff = 2. * M_PI / angle;
      typename std::list<Angle_map_iterator>::iterator lit = entries.begin();
      typename std::list<Angle_map_iterator>::iterator lend = entries.end();
      for(; lit!=lend; ++lit)
      {
        const Angle_map_iterator mit = *lit;
        mit->second *= norm_coeff;
      }
    }
  }

  void initialize()
  {
#if (VERBOSITY > 5)
    std::cout << "canvas initialization" << std::endl;
#endif

    // build the canvas points
    std::size_t vertex_counter = 0;

    typename Tr::Finite_vertices_iterator vit = m_c3t3.triangulation().finite_vertices_begin();
    typename Tr::Finite_vertices_iterator vend = m_c3t3.triangulation().finite_vertices_end();
    for(; vit!=vend; ++vit)
    {
      if(vit->in_dimension() > 2)
        continue; // keep the boundary vertices

      // note that vit->point() is a weighted_point
      Canvas_point cp(vit->point(), vertex_counter++, vit, this);

      // border info is meaningless for now since we consider closed surfaces
      // fixme if the c2t3 has (1-dimensional) borders
      cp.border_info() = false;

      this->canvas_points.push_back(cp);
      vit->info() = this->canvas_points.size() - 1;
      CGAL_postcondition(this->canvas_points[vit->info()].point() == vit->point().point());
    }

    std::cout << this->canvas_points.size() << " points on canvas" << std::endl;

    initialize_cumulative_angles();

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

#ifndef ANISO_GEO_FORCE_SEED_INITIALIZATION
    CGAL_assertion(vh->in_dimension() <= 2);
#endif
    Canvas_point& cp = this->canvas_points[vh->info()];
    const Metric& v_m = cp.metric();
    const Eigen::Matrix3d& f = get_interpolated_transformation(seed_m, v_m);

    Vector3d v;
    v(0) = t.x() - cp.point().x();
    v(1) = t.y() - cp.point().y();
    v(2) = t.z() - cp.point().z();
    v = f * v;

    // (f*v).norm() = std::sqrt(v^t * M * v)
    Base::initialize_canvas_point(cp, v.norm(), seed_id);
  }

  void initialize_seed_on_vertex(Cell_handle c, int li,
                                 std::size_t seed_id,
                                 const Point_3& t)
  {
    Vertex_handle vh = c->vertex(li);

    if(c->vertex(li)->in_dimension() > 2)
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

    initialize_vertex(vh, seed_id, t);
  }

  void initialize_seed_on_edge(Cell_handle c, int li, int lj,
                               std::size_t seed_id,
                               const Point_3& t)
  {
    if(c->vertex(li)->in_dimension() > 2 ||
       c->vertex(lj)->in_dimension() > 2)
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

    for(int i=0; i<4; ++i)
      if(i == li || i == lj)
        initialize_vertex(c->vertex(i), seed_id, t);
  }

  void initialize_seed_on_triangle(Cell_handle c, int li,
                                   std::size_t seed_id,
                                   const Point_3& t)
  {
    if(!m_c3t3.is_in_complex(c, li))
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

    for(int i=0; i<4; ++i)
      if(i != li)
        initialize_vertex(c->vertex(i), seed_id, t);
  }

  void initialize_seed_in_cell(Cell_handle c,
                               std::size_t seed_id,
                               const Point_3& t)
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

    for(std::size_t j=0; j<4; ++j)
      initialize_vertex(c->vertex(j), seed_id, t);
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
    Cell_handle c = m_c3t3.triangulation().locate(t, lt, li, lj);
    CGAL_assertion(lt < 3); // must be a vertex, an edge or a triangle

    if(!m_c3t3.is_in_complex(c))
      std::cout << "WARNING: seed located in a cell outside of the complex" << std::endl;

#if (VERBOSITY > 15)
    std::cout << "locate in triangulation : " << std::endl;
    std::cout << "locate type: " << lt << " li/j: " << li << " " << lj << std::endl;
    std::cout << "cell : " << &*c << " with info : " << c->info() << std::endl;
    for(int i=0; i<4; ++i)
    {
      Vertex_handle vh = c->vertex(i);
      std::cout << "point: " << vh->info() << std::endl;
      if(!m_c3t3.triangulation().is_infinite(vh))
      {
        std::cout << c->vertex(i)->point() << " || "
                  << this->canvas_points[vh->info()].point() << std::endl;
      }
    }
#endif

    if(lt == 0) // vertex
      initialize_seed_on_vertex(c, li, seed_id, t);
    else if(lt == 1) // edge
      initialize_seed_on_edge(c, li, lj, seed_id, t);
    else if(lt == 2) // facet
      initialize_seed_on_triangle(c, li, seed_id, t);
    else
    {
      std::cerr << "WARNING: Initialization of a seed in a cell..." << std::endl;
      initialize_seed_in_cell(c, seed_id, t);
    }
  }

  void primal_shenanigans(const Canvas_point* cp)
  {
#if (VERBOSITY > 15)
    std::cout << "primal shenanigans at : " << cp->index() << std::endl;
#endif

    CGAL_assertion(cp->state() == KNOWN);
    typedef boost::unordered_set<const Canvas_point*>   Candidates_set;

    Facet_handle fit = cp->incident_facets_in_complex_begin();
    Facet_handle end = cp->incident_facets_in_complex_end();

#if (VERBOSITY > 25)
    std::cout << "cp: " << cp->index() << " (" << cp->point() << ") ";
    std::cout << "has " << end - fit << " incident facets in complex" << std::endl;
#endif
    for(; fit!=end; ++fit)
    {
      CGAL_assertion(m_c3t3.is_in_complex(fit->first, fit->second));

      Candidates_set candidates;
      candidates.insert(cp);

      for(std::size_t j=1; j<4; ++j)
      {
        const Vertex_handle vh = fit->first->vertex((fit->second+j)%4);
        CGAL_assertion(vh->in_dimension() < 3); // has to be on the surface
        const Canvas_point& cq = this->canvas_points[vh->info()];
        if(cq.state() != KNOWN)
          continue;

        candidates.insert(&cq);
      }

      Base::construct_primal_elements_from_candidates(candidates);
      //todo mark/compute Voronoi vertices
    }
  }

  void compute_local_primal_elements(const Canvas_point* cp)
  {
    primal_shenanigans(cp);
  }

  void compute_primal()
  {
    if(this->primal_edges.empty() && this->primal_triangles.empty())
    {
#if (VERBOSITY > 5)
    std::cout << "Primal computations" << std::endl;
#endif
      typedef boost::unordered_set<const Canvas_point*>   Candidates_set;
      Facet_iterator fit = m_c3t3.facets_in_complex_begin();
      Facet_iterator end = m_c3t3.facets_in_complex_end();
      for(; fit!=end; ++fit)
      {
        Cell_handle c = fit->first;
        std::size_t sec = fit->second;

        std::size_t color_0 =
            this->canvas_points[c->vertex((sec+1)%4)->info()].closest_seed_id();
        std::size_t color_1 =
            this->canvas_points[c->vertex((sec+2)%4)->info()].closest_seed_id();
        std::size_t color_2 =
            this->canvas_points[c->vertex((sec+3)%4)->info()].closest_seed_id();

        // filter facets with only one color
        if(color_0 == color_1 && color_1 == color_2)
          continue;

        // todo: when the c2t3 has borders, you can add some filters again (see c3t3)

        Candidates_set candidates;
        for(std::size_t j=1; j<4; ++j)
          candidates.insert(&(this->canvas_points[c->vertex((sec+j)%4)->info()]));

        Base::construct_primal_elements_from_candidates(candidates);
      }
      // todo add mark/compute Voronoi vertices etc.
    }
    else
      std::cout << "WARNING: call to compute_primal with non-empty primal data structures..." << std::endl;

    this->check_canvas_density_with_primals();
  }

  void output_canvas(const std::string str_base) const
  {
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

    out << "Triangles" << std::endl;
    out << m_c3t3.number_of_facets() << std::endl;
    typename C3t3::Facets_in_complex_iterator fit = m_c3t3.facets_in_complex_begin();
    typename C3t3::Facets_in_complex_iterator fend = m_c3t3.facets_in_complex_end();
    for(; fit!=fend; ++fit)
    {
      Cell_handle c = fit->first;
      std::size_t sec = fit->second;

      boost::unordered_set<std::size_t> materials;
      for(std::size_t j=1; j<4; ++j)
        materials.insert(this->canvas_points[c->vertex((sec+j)%4)->info()].closest_seed_id());

#ifdef ONLY_PRINT_BISECTORS
      if(materials.size() == 1)
        continue;
#endif

      for(std::size_t j=1; j<=3; ++j)
        out << c->vertex((sec+j)%4)->info() + 1 << " ";

      std::size_t mat = (materials.size()==1) ? (*(materials.begin())) :
                                                (this->seeds.size());
      out << mat << std::endl;
    }

/*
    out << "Tetrahedra" << std::endl;
    out << m_c3t3.number_of_cells() << std::endl;
    typename C3t3::Cell_iterator cit = m_c3t3.cells_begin();
    typename C3t3::Cell_iterator cend = m_c3t3.cells_end();
    for(; cit!=cend; ++cit)
    {
      for(std::size_t j=0; j<4; ++j)
        out << cit->vertex(j)->info() + 1 << " ";
      out << cit->subdomain_index() << std::endl;
    }
*/

    out_bb << "End" << std::endl;
    out << "End" << std::endl;
  }

  Campen_canvas(C3t3& c3t3,
                const std::string& canvas_str_,
                const std::string& seeds_str_,
                const std::size_t max_seeds_n_,
                const Metric_field* mf_)
    :
      Base(canvas_str_, seeds_str_, max_seeds_n_, mf_),
      m_c3t3(c3t3),
      m_angles()
  { }
};

} // namespace Anisotropic_mesh_3
} // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_CAMPEN_C2T3_CANVAS_H

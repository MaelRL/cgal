#ifndef CGAL_ANISOTROPIC_MESH_3_CAMPEN_STARSET_CANVAS_H
#define CGAL_ANISOTROPIC_MESH_3_CAMPEN_STARSET_CANVAS_H

#include <CGAL/Canvas/canvas_config.h>
#include <CGAL/Canvas/canvas.h>
#include <CGAL/Canvas/Campen_starset_point.h>
#include <CGAL/Canvas/canvas_subdiviser.h>

#include <CGAL/Anisotropic_mesher_3.h>
#include <CGAL/Starset.h>
#include <CGAL/IO/Star_set_IO.h>

#include <CGAL/helpers/combinatorics_helper.h>

#include <CGAL/aabb_tree/aabb_tree_bbox.h>
#include <CGAL/aabb_tree/aabb_tree_bbox_primitive.h>

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

template<typename K, typename Domain, typename MF, typename Criteria>
class Campen_starset_canvas :
    public Canvas<K, Campen_starset_point<K, Domain, MF, Criteria>, MF>
{
private:
  typedef Campen_starset_canvas<K, Domain, MF, Criteria>          Self;

public:
  typedef Campen_starset_point<K, Domain, MF, Criteria>           Canvas_point;
  typedef Canvas_point*                                           Canvas_point_handle;

  typedef int                                  Vertex_Info; // index of the canvas point
  typedef int                                  Cell_Info; // index of the subdomain

  typedef Canvas<K, Canvas_point, MF>                             Base;

  typedef typename Base::FT                                       FT;
  typedef typename Base::Point_3                                  Point_3;
  typedef typename Base::Tetrahedron_3                            Tetrahedron_3;
  typedef typename Base::Metric                                   Metric;
  typedef typename Base::Vector3d                                 Vector3d;
  typedef typename Base::Simplex                                  Simplex;
  typedef typename Base::BEdge                                    BEdge;
  typedef typename Base::BTriangle                                BTriangle;
  typedef typename Base::BTetrahedron                             BTetrahedron;

  typedef Starset_with_info<K, Domain, MF, Criteria>              Star_set;
  typedef Stretched_Delaunay_3<K>                                 Star;
  typedef typename Star_set::Star_handle                          Star_handle;
  typedef typename Star_set::Vertex_handle                        Vertex_handle;
  typedef typename Star_set::Vertex_handle_handle                 Vertex_handle_handle;
  typedef typename Star_set::Cell_handle                          Cell_handle;
  typedef typename Star_set::Cell_handle_handle                   Cell_handle_handle;

  typedef CGAL::AABB_tree_bbox<K, Star>                           AABB_tree;

  Star_set* m_ss;
  AABB_tree m_aabb_tree;

  void build_aabb_tree()
  {
    m_aabb_tree.rebuild(m_ss->begin(), m_ss->end());
  }

  void read_dump()
  {
    std::cout << "Reading dump..." << std::endl;
    std::ifstream in("dump_wip.txt");
    m_ss->clear();

    std::size_t stars_n, v_n, id;
    FT x, y, z;

    in >> stars_n;

    for(std::size_t i=0; i<stars_n; ++i)
    {
      in >> x >> y >> z;
      Point_3 p(x,y, z);

      Star_handle star = new typename Star_set::Star(m_ss->criteria(),
                                                     m_ss->constrain_surface());
      const Metric& m_p = this->mf->compute_metric(p);
      star->reset(p, i, m_p);
      m_ss->push_back(star);
    }

    for(std::size_t i=0; i<stars_n; ++i)
    {
      in >> v_n;
      Star_handle star_i = m_ss->get_star(i);
      for(std::size_t j=0; j<v_n; ++j)
      {
        in >> id;
        Star_handle star_j = m_ss->get_star(id);
        star_i->insert_to_star(star_j->center_point(), id, false/*conditional*/);
      }
    }
    std::cout << m_ss->size() << " stars from dump" << std::endl;
    build_aabb_tree();
  }

  void initialize()
  {
#if (VERBOSITY > 5)
    std::cout << "canvas initialization" << std::endl;
#endif

    // read and create the starset
#if 1//def STARSET_FROM_DUMP
    read_dump();
#else
    std::ifstream in((this->canvas_str + ".mesh").c_str());
    build_starset<Star_set>(in);
#endif

    // build the canvas points
    std::size_t vertex_counter = 0;
    for(std::size_t i=0, sss=m_ss->size(); i<sss; ++i)
    {
      Star_handle star_i = m_ss->get_star(i);
      bool is_on_border = has_a_corner_vertex(star_i);
      Canvas_point cp(star_i->center_point(), vertex_counter++, is_on_border, m_ss, this);
      this->canvas_points.push_back(cp);
    }
    CGAL_postcondition(this->canvas_points.size() == m_ss->size());

    Base::compute_canvas_bbox();
    this->seeds.initialize_seeds();
    Base::locate_seeds_on_canvas();

#if (VERBOSITY > 5)
    std::cout << "canvas initialized with " << this->seeds.size() << " seeds" << std::endl;
#endif
  }

  void initialize_vertex(std::size_t id,
                         const std::size_t seed_id,
                         const Point_3& t)
  {
    CGAL_precondition(id != static_cast<std::size_t>(-1));
    const Metric& seed_m = this->seeds.seeds_metrics[seed_id];

    Canvas_point& cp = this->canvas_points[id];
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

  void find_and_initialize_closest_vertex(const std::size_t seed_id,
                                          const Point_3& t)
  {
    // this is pretty expensive, but only used when the aabb tree failed to find a star
    FT min_d = FT_inf;
    std::size_t min_id = -1;
    for(std::size_t i=0, sss=m_ss->size(); i<sss; ++i)
    {
      Star_handle star = m_ss->get_star(i);
      typename K::Compute_squared_distance_3 sqd =
                                        K().compute_squared_distance_3_object();
      FT d = sqd(star->center_point(), t);
      if(d < min_d)
      {
        min_id = i;
        min_d = d;
      }
    }
    CGAL_postcondition(min_id != static_cast<std::size_t>(-1));

    Canvas_point& cp = this->canvas_points[min_id];
    Base::initialize_canvas_point(cp, std::sqrt(min_d), seed_id);
  }

  void initialize_seed_in_cell(Cell_handle c, std::size_t seed_id, const Point_3& t)
  {
    for(std::size_t j=0; j<4; ++j)
      initialize_vertex(c->vertex(j)->info(), seed_id, t);
  }

  void locate_and_initialize(const Point_3& t, const std::size_t seed_id)
  {
    // find the tetrahedron that contains the seed and initialize the distance at
    // these four points accordingly.
    bool found = false;

    std::list<Star_handle> intersections;
    m_aabb_tree.all_intersected_primitives(t, std::back_inserter(intersections));

    if(intersections.empty())
    {
      std::cout << "the locate-ray didn't encounter any star..." << std::endl;
      // switch to something expensive... take the closest point...
      find_and_initialize_closest_vertex(seed_id, t);
      found = true;
    }

    typename std::list<Star_handle>::const_iterator it = intersections.begin(),
                                                    end = intersections.end();
    for(; it!=end; ++it)
    {
      Star_handle star = *it;
      Cell_handle_handle cit = star->finite_star_cells_begin();
      Cell_handle_handle cend = star->finite_star_cells_end();
      for(; cit!=cend; ++cit)
      {
        Cell_handle c = *cit;
        if(!star->is_inside(c))
          continue;

        // expensive stuff...
        Tetrahedron_3 tetra(m_ss->get_star(c->vertex(0)->info())->center_point(),
                            m_ss->get_star(c->vertex(1)->info())->center_point(),
                            m_ss->get_star(c->vertex(2)->info())->center_point(),
                            m_ss->get_star(c->vertex(3)->info())->center_point());

        typename Star::Traits traits;

        if(traits.has_on_bounded_side_3_object()(tetra, t) ||
           traits.has_on_boundary_3_object()(tetra, t))
        {
          std::cout << "initialize through star: " << star->index_in_star_set() << std::endl;
          initialize_seed_in_cell(c, seed_id, t);
          found = true;
          return; // only initialize one cell
        }
      }
    }
    std::cout << "Failed to initialize a seed" << std::endl;
    CGAL_assertion(false);
  }

  void primal_shenanigans(const Canvas_point* cp)
  {
#if (VERBOSITY > 15)
    std::cout << "primal shenanigans at : " << cp->index() << std::endl;
#endif

    CGAL_assertion(cp->state() == KNOWN);
    typedef boost::unordered_set<const Canvas_point*>   Candidates_set;

    Star_handle star = m_ss->get_star(cp->index());
    Cell_handle_handle cit = star->finite_star_cells_begin();
    Cell_handle_handle cend = star->finite_star_cells_end();

#if (VERBOSITY > 25)
    std::cout << "cp: " << cp->index() << " (" << cp->point() << ") ";
    std::cout << "has " << cend - cit << " incident interior cells" << std::endl;
#endif

    for(; cit!=cend; ++cit)
    {
      const Cell_handle c = *cit;
      if(!star->is_inside(c))
        continue;

      Candidates_set candidates;
      candidates.insert(cp);

      for(std::size_t j=0; j<4; ++j)
      {
        const Vertex_handle vh = c->vertex(j);
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

  void fix_distance_and_seed_id_at_corners()
  {
    // corners are points outside the domain lazily added by the starset so
    // it doesn't have to deal with the refinement of infinite cells
    // these points are the first four stars
    for(std::size_t i=0; i<8; ++i)
    {
      this->canvas_points[i].closest_seed_id() = -1;
      this->canvas_points[i].distance_to_closest_seed() = 0;
    }
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
    if(!this->primal_edges.empty() || !this->primal_triangles.empty() ||
       !this->primal_tetrahedra.empty())
      std::cout << "WARNING: call to compute_primal with non-empty primal data structures..." << std::endl;

#if (VERBOSITY > 5)
    std::cout << "Primal computations" << std::endl;
#endif
    typedef boost::unordered_set<const Canvas_point*>   Candidates_set;
    for(std::size_t i=0, ts=m_ss->size(); i<ts; ++i)
    {
      const Star_handle star = m_ss->get_star(i);

      Cell_handle_handle cit = star->finite_star_cells_begin();
      Cell_handle_handle cend = star->finite_star_cells_end();
      for(; cit!=cend; ++cit)
      {
        Cell_handle c = *cit;

        // ignore exterior cells to ignore triangulation neighbors that are not canvas neighbors
        if(!star->is_inside(c))
          continue;

        std::size_t color_0 = this->canvas_points[c->vertex(0)->info()].closest_seed_id();
        std::size_t color_1 = this->canvas_points[c->vertex(1)->info()].closest_seed_id();
        std::size_t color_2 = this->canvas_points[c->vertex(2)->info()].closest_seed_id();
        std::size_t color_3 = this->canvas_points[c->vertex(3)->info()].closest_seed_id();

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
          if(!this->canvas_points[c->vertex(0)->info()].border_info() &&
             !this->canvas_points[c->vertex(1)->info()].border_info() &&
             !this->canvas_points[c->vertex(2)->info()].border_info() &&
             !this->canvas_points[c->vertex(3)->info()].border_info())
            continue;
        }
#endif

        Candidates_set candidates;
        for(std::size_t j=0; j<4; ++j)
          candidates.insert(&(this->canvas_points[c->vertex(j)->info()]));

        Base::construct_primal_elements_from_candidates(candidates);
      }
    }
    // todo add mark/compute Voronoi vertices etc.

#ifdef CGAL_ANISO_REFINE_LOW_CANVAS_DENSITY
    check_canvas_density();
#else
    this->check_canvas_density_with_primals();
#endif
  }

  void output_canvas(const std::string str_base) const
  {
    std::ofstream out((str_base + ".mesh").c_str());
    std::ofstream out_bb((str_base + ".bb").c_str());

    out << "MeshVersionFormatted 1" << '\n';
    out << "Dimension 3" << '\n';
    out << "Vertices" << '\n';
    out << this->canvas_points.size() << '\n';
    out_bb << "3 1 " << this->canvas_points.size() << " 2" << '\n';

    for(std::size_t i=0; i<this->canvas_points.size(); ++i)
    {
      const Canvas_point& cp = this->canvas_points[i];
      out << cp.point() << " " << i+1 << '\n';

//      out_bb << cp.closest_seed_id < < std::endl;
      out_bb << cp.distance_to_closest_seed() << '\n';
    }

    Facet_ijk_unordered_set output_facets;
    Cell_ijkl_unordered_set output_cells;
    for(std::size_t i=0, sss=m_ss->size(); i<sss; ++i)
    {
      Star_handle star = m_ss->get_star(i);
      Cell_handle_handle ci = star->finite_star_cells_begin();
      Cell_handle_handle ciend = star->finite_star_cells_end();
      for (; ci!=ciend; ++ci)
      {
        Cell_handle c = *ci;
        if(!star->is_inside(*ci))
          continue;

        output_cells.insert(Cell_ijkl(c));
        output_facets.insert(Facet_ijk(c->vertex(0)->info(), c->vertex(1)->info(),
                                       c->vertex(2)->info()));
        output_facets.insert(Facet_ijk(c->vertex(0)->info(), c->vertex(2)->info(),
                                       c->vertex(3)->info()));
        output_facets.insert(Facet_ijk(c->vertex(0)->info(), c->vertex(1)->info(),
                                       c->vertex(3)->info()));
        output_facets.insert(Facet_ijk(c->vertex(3)->info(), c->vertex(1)->info(),
                                       c->vertex(2)->info()));
      }
    }

    out << "Tetrahedra\n";
    out << output_cells.size() << std::endl;
    typename Cell_ijkl_unordered_set::iterator cit = output_cells.begin();
    typename Cell_ijkl_unordered_set::iterator citend = output_cells.end();
    for (; cit != citend; cit++)
    {
      Cell_ijkl c = *cit;
      int n0 = c.vertices()[0];
      int n1 = c.vertices()[1];
      int n2 = c.vertices()[2];
      int n3 = c.vertices()[3];

      out << (n0+1) << " " << (n1+1) << " " << (n2+1) << " " << (n3+1);
#if 0
      Cell_handle useless;
      bool is_consistent = (m_ss->get_star(n0)->has_cell(c, useless) &&
                            m_ss->get_star(n1)->has_cell(c, useless) &&
                            m_ss->get_star(n2)->has_cell(c, useless) &&
                            m_ss->get_star(n3)->has_cell(c, useless) );
      out << " " << is_consistent << '\n';
#else
      boost::unordered_set<std::size_t> materials;
      for(std::size_t j=0; j<4; ++j)
        materials.insert(this->canvas_points[c.vertices()[j]].closest_seed_id());

      std::size_t mat = (materials.size()==1) ? (*(materials.begin())) :
                                                (this->seeds.size());
      out << " " << mat << std::endl;
#endif
    }
    out_bb << "End" << std::endl;
    out << "End" << std::endl;
  }

  Campen_starset_canvas(const std::string& canvas_str_,
                        const std::string& seeds_str_,
                        const std::size_t max_seeds_n_,
                        Star_set* ss_)
    :
      Base(canvas_str_, seeds_str_, max_seeds_n_, ss_->metric_field()),
      m_ss(ss_),
      m_aabb_tree()
  { }
};

} // namespace Anisotropic_mesh_3
} // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_CAMPEN_STARSET_CANVAS_H

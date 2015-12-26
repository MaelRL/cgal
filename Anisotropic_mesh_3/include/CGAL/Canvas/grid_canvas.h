#ifndef CGAL_ANISOTROPIC_MESH_3_GRID_CANVAS_H
#define CGAL_ANISOTROPIC_MESH_3_GRID_CANVAS_H

#include <CGAL/Canvas/canvas.h>

#include <boost/array.hpp>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <ostream>
#include <set>
#include <string>
#include <vector>

// Case where the canvas is an orthogonal grid

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K, typename Cpoint, typename MF>
class Grid_canvas :
    public Canvas<K, Cpoint, MF>
{
public:
  typedef Canvas<K, Cpoint, MF>                                 Base;

  typedef Cpoint                                                Canvas_point;
  typedef MF                                                    Metric_field;

  typedef typename Canvas_point::Neighbors                      Neighbors;

  typedef typename K::FT                                        FT;
  typedef typename K::Point_3                                   Point_3;

  typedef typename Base::Simplex                                Simplex;
  typedef typename Base::BEdge                                  BEdge;
  typedef typename Base::BTriangle                              BTriangle;
  typedef typename Base::BTetrahedron                           BTetrahedron;
  typedef typename Base::Canvas_point_handle_vector             Canvas_point_handle_vector;

  typedef typename Base::Vector3d                               Vector3d;

  // Grid geometry
  const Point_3 center;
  const FT side;
  const std::size_t n; // number of points per side

  FT offset_x; // offset is the bottom left point
  FT offset_y;
  FT offset_z;
  const FT step;
  const FT sq_n;

  void locate_and_initialize(const Point_3& p,
                             const std::size_t seed_id)
  {
    int index_x = std::floor((p.x()-offset_x)/step);
    int index_y = std::floor((p.y()-offset_y)/step);
    int index_z = std::floor((p.z()-offset_z)/step);

    Canvas_point& cp = this->canvas_points[index_z*sq_n + index_y*n + index_x];

#if (VERBOSITY > 10)
    std::cout << "looking for p: " << p << std::endl;
    std::cout << "found cp: " << cp.index() << " [" << cp.point() << "]" << std::endl;
    std::cout << "index: " << index_x << " " << index_y << " " << index_z << std::endl;
    std::cout << "offset: " << offset_x << " " << offset_y << " " << offset_z << std::endl;
#endif

    Vector3d v;
    v(0) = p.x() - cp.point().x();
    v(1) = p.y() - cp.point().y();
    v(2) = p.z() - cp.point().z();
    const Eigen::Matrix3d& m = cp.metric().get_mat();
    FT d = std::sqrt(v.transpose() * m * v);

    Base::initialize_canvas_point(cp, d, seed_id);
  }

  void initialize()
  {
#if (VERBOSITY > 5)
    std::cout << "canvas initialization" << std::endl;
#endif

    // create the canvas points
    for(std::size_t k=0; k<n; ++k)
    {
      for(std::size_t j=0; j<n; ++j)
      {
        for(std::size_t i=0; i<n; ++i)
        {
          // fill from bot left to top right
          Point_3 p(offset_x+i*step, offset_y+j*step, offset_z+k*step);
          Canvas_point cp(p, i + j*n + k*sq_n, this);
          this->canvas_points.push_back(cp);
        }
      }
    }

    // assign the neighbors
    for(std::size_t k=0; k<n; ++k)
    {
      for(std::size_t j=0; j<n; ++j)
      {
        for(std::size_t i=0; i<n; ++i)
        {
          std::size_t curr_id = i + j*n + k*sq_n;

          if(k != n-1) // there is a neighbor above
            this->canvas_points[curr_id].neighbors[0] = i + j*n + (k+1)*sq_n;
          if(i != 0) // there is a neighbor left
            this->canvas_points[curr_id].neighbors[1] = i-1 + j*n + k*sq_n;
          if(j != n-1) // there is a neighbor back
            this->canvas_points[curr_id].neighbors[2] = i + (j+1)*n + k*sq_n;
          if(i != n-1) // there is a neighbor right
            this->canvas_points[curr_id].neighbors[3] = i+1 + j*n + k*sq_n;
          if(j != 0) // there is a neighbor front
            this->canvas_points[curr_id].neighbors[4] = i + (j-1)*n + k*sq_n;
          if(k != 0) // there is a neighbor below
            this->canvas_points[curr_id].neighbors[5] = i + j*n + (k-1)*sq_n;
        }
      }
    }
#if (VERBOSITY > 6)
    std::cout << "neighbors assigned" << std::endl;
#endif

    Base::compute_canvas_bbox();
    this->seeds.initialize_seeds();
    Base::locate_seeds_on_canvas();

#if (VERBOSITY > 6)
    std::cout << "canvas initialized" << std::endl;
#endif
  }

  void add_tet_to_simplices(boost::unordered_set<BTetrahedron>& tetrahedra,
                            const Canvas_point& cp,
                            const Canvas_point* n1,
                            const Canvas_point* n2,
                            const Canvas_point* n3) const
  {
    BTetrahedron t;
    if(n1 && n2 && n3)
    {
      t[0] = cp.index();
      t[1] = n1->index();
      t[2] = n2->index();
      t[3] = n3->index();
      tetrahedra.insert(t);
    }
  }

  void add_tet_to_simplices(boost::unordered_set<BTetrahedron>& tetrahedra,
                            const Canvas_point& cp, std::size_t i1,
                            std::size_t i2, std::size_t i3) const
  {
    return add_tet_to_simplices(tetrahedra, cp, cp.neighbors[i1],
                                cp.neighbors[i2],cp.neighbors[i3]);
  }

  void primal_shenanigans(const Canvas_point* cp)
  {
    CGAL_assertion(cp->state() == KNOWN);
    // Check the neighbors of cp for points that have the state KNOWN.
    // Amongst these, consider those that have a closest_seed_id different than cp's
    // to obtain primal simplices (of dim > 1).

    typedef boost::unordered_set<const Canvas_point*>   Candidates_set;
    Candidates_set candidates;
    candidates.insert(cp);

    typename Neighbors::const_iterator it = cp->neighbors.begin(),
                                       iend = cp->neighbors.end();
    for(; it!=iend; ++it)
    {
      if(*it == static_cast<std::size_t>(-1))
        continue;

      const Canvas_point& cq = this->canvas_points[*it];
      if(cq.state() != KNOWN)
        continue;
      candidates.insert(&cq);
    }
    Base::construct_primal_elements_from_candidates(candidates);
  }

  void compute_local_primal_elements(const Canvas_point* cp)
  {
    primal_shenanigans(cp);
    typename Neighbors::const_iterator it = cp->neighbors.begin(),
                                       iend = cp->neighbors.end();
    for(; it!=iend; ++it)
    {
      if(*it == static_cast<std::size_t>(-1))
        continue;

      const Canvas_point& cq = this->canvas_points[*it];
      if(cq.state() == KNOWN)
        primal_shenanigans(&cq);
    }
  }

  void cosphericity_sweeper() const
  {
    // we define a tiny (2k+1)*(2k+1) grid ('sweeper'), and we look at the colors
    // if we have more than 5 colors at a time, we have a (quasi-) cosphericity

    typedef boost::unordered_set<std::size_t> Color_set;

    std::set<boost::array<std::size_t, 5> > quasi_cosphericities;

    const int h = (std::max)(static_cast<std::size_t>(1), n/100);
    const std::size_t sweeper_size = 2*h+1;
    if(sweeper_size > n)
    {
      std::cerr << "WARNING: sweeper too big for that canvas !" << std::endl;
      return;
    }

    // we move the center of the sweeper in the canvas
    for(std::size_t i=h; i<n-h; ++i)
    {
      for(std::size_t j=h; j<n-h; ++j)
      {
        for(std::size_t k=h; k<n-h; ++k)
        {
          //const Canvas_point& sweeper_center = points[j*n + i];

          // gather the pts
          Color_set colors;
          for(int lx=-h; lx<=h; ++lx)
          {
            for(int ly=-h; ly<=h; ++ly)
            {
              for(int lz=-h; lz<=h; ++lz)
              {
                const Canvas_point& cp = this->canvas_points[(k+lz)*sq_n + (j+ly)*n + (i+lx)];
                colors.insert(cp.closest_seed_id());
              }
            }
          }

          if(colors.size() > 4)
          {
            // we have at least 5 pts in a quasi cospherical configuration...
            // get them all with combinations!
            typedef std::vector<boost::array<std::size_t, 5> >    RType;
            RType combis = combinations<5>(colors);
            quasi_cosphericities.insert(combis.begin(), combis.end());
          }
        }
      }
    }

    std::cout << "Sweeper no sweeping ! (with size: " << h << ")" << std::endl;
    std::cout << "Found " << quasi_cosphericities.size() << " (quasi-)cosphericities" << std::endl;

    std::set<boost::array<std::size_t, 5> >::const_iterator sit =
                                                   quasi_cosphericities.begin();
    std::set<boost::array<std::size_t, 5> >::const_iterator send =
                                                     quasi_cosphericities.end();
    for(; sit!=send; ++sit)
    {
      const boost::array<std::size_t, 5>& quasi_cosphericity = *sit;
      for(std::size_t i=0; i<sit->size(); ++i)
        std::cout << quasi_cosphericity[i] << " ";
      std::cout << std::endl;
    }
    std::cout << "End of quasi cosphericities" << std::endl;
  }

  void output_canvas(const std::string str_base) const
  {
#if (VERBOSITY > 3)
    std::cout << "Output canvas" << std::endl;
#endif

    std::ofstream out((str_base + ".mesh").c_str());
    std::ofstream out_bb((str_base + ".bb").c_str());
    out << "MeshVersionFormatted 1" << std::endl;
    out << "Dimension 3" << std::endl;
    out << "Vertices" << std::endl;
    out << this->canvas_points.size() << std::endl;
    out_bb << "3 1 " << this->canvas_points.size() << " 2" << std::endl;

    boost::unordered_set<BTetrahedron> tetrahedra;
    for(std::size_t i=0; i<this->canvas_points.size(); ++i)
    {
      const Canvas_point& cp = this->canvas_points[i];
      out << cp.point() << " " << i+1 << std::endl;

//      out_bb << cp.closest_seed_id << std::endl;
      out_bb << cp.distance_to_closest_seed() << std::endl;

      // compute the tets...
      const Neighbors& ns = cp.neighbors; // just to get a shorter name...
      for(std::size_t i=0; i<ns.size(); ++i)
      {
#if 0
        // lazily : the eight corner tetrahedra (creates self intersections,
        // but it doesn't really matter except visually)
        add_tet_to_simplices(tetrahedra, cp, 0, 1, 2);
        add_tet_to_simplices(tetrahedra, cp, 0, 2, 3);
        add_tet_to_simplices(tetrahedra, cp, 0, 3, 4);
        add_tet_to_simplices(tetrahedra, cp, 0, 4, 1);
        add_tet_to_simplices(tetrahedra, cp, 5, 1, 2);
        add_tet_to_simplices(tetrahedra, cp, 5, 2, 3);
        add_tet_to_simplices(tetrahedra, cp, 5, 3, 4);
        add_tet_to_simplices(tetrahedra, cp, 5, 4, 1);
#else
        // a bit smarter : a cube is decomposed in 5 tets and we pick them
        // correctly so we have no overlap
        add_tet_to_simplices(tetrahedra, cp, 0, 2, 3);
        add_tet_to_simplices(tetrahedra, cp, 0, 1, 4);
        add_tet_to_simplices(tetrahedra, cp, 3, 4, 5);
        add_tet_to_simplices(tetrahedra, cp, 1, 2, 5);

        // the last one is a bit annoying since we have to grab neighbors of neighbors
        if(cp.neighbors[3] != static_cast<std::size_t>(-1) &&
           cp.neighbors[4] != static_cast<std::size_t>(-1))
        {
          const Canvas_point& n3 = this->canvas_points[cp.neighbors[3]];
          const Canvas_point& n4 = this->canvas_points[cp.neighbors[4]];
          add_tet_to_simplices(tetrahedra, cp, n4.neighbors[3],
                               n4.neighbors[0], n3.neighbors[0]);
        }
#endif
      }
    }

    out << "Tetrahedra" << std::endl;
    out << tetrahedra.size() << std::endl;
    for(typename boost::unordered_set<BTetrahedron>::iterator it = tetrahedra.begin();
                                                              it != tetrahedra.end(); ++it)
    {
      const BTetrahedron& tet = *it;
      boost::unordered_set<std::size_t> materials;

      for(std::size_t i=0; i<tet.size(); ++i)
      {
        out << tet[i] + 1 << " ";
        materials.insert(this->canvas_points[tet[i]].closest_seed_id());
      }
      std::size_t mat = (materials.size()==1)? (*(materials.begin())) :
                                               (this->seeds.size());
      out << mat << std::endl;
    }

    boost::unordered_set<BTriangle> triangles;
    for(typename boost::unordered_set<BTetrahedron>::iterator it = tetrahedra.begin();
                                                 it != tetrahedra.end(); ++it)
    {
      const BTetrahedron& tet = *it;
      BTriangle tr; // could use a permutation to make it look nicer, I suppose... todo
      tr[0] = tet[0]; tr[1] = tet[1]; tr[2] = tet[2]; triangles.insert(tr);
      tr[0] = tet[0]; tr[1] = tet[1]; tr[2] = tet[3]; triangles.insert(tr);
      tr[0] = tet[0]; tr[1] = tet[2]; tr[2] = tet[3]; triangles.insert(tr);
      tr[0] = tet[1]; tr[1] = tet[2]; tr[2] = tet[3]; triangles.insert(tr);
    }

    out << "Triangles" << std::endl;
    out << triangles.size() << std::endl;
    for(typename boost::unordered_set<BTriangle>::iterator it = triangles.begin();
                                                           it != triangles.end(); ++it)
    {
      const BTriangle& tr = *it;
      boost::unordered_set<std::size_t> materials;
      for(std::size_t i=0; i<tr.size(); ++i)
      {
        out << tr[i] + 1 << " ";
        materials.insert(this->canvas_points[tr[i]].closest_seed_id());
      }
      std::size_t mat = (materials.size()==1) ? (*(materials.begin())) :
                                                (this->seeds.size());
      out << mat << std::endl;
    }
  }

  void compute_primal()
  {
    // this version computes the primal by looping over the canvas points
    // and looking at the amount of colors in said point and said point's neighbors
#if (VERBOSITY > 10)
    std::cout << "primal computations" << std::endl;
#endif
    if(!this->primal_edges.empty() || !this->primal_triangles.empty() || !this->primal_tetrahedra.empty())
    {
      std::cerr << "WARNING: call to compute_primal with non-empty primal data structures..." << std::endl;
      return;
    }

    typedef boost::unordered_set<const Canvas_point*>   Candidates_set;
    for(std::size_t i=0; i<this->canvas_points.size(); ++i)
    {
      const Canvas_point& cp = this->canvas_points[i];
      Candidates_set candidates;
      candidates.insert(&cp);

      typename Canvas_point::Neighbors::const_iterator it = cp.neighbors.begin(),
                                                       iend = cp.neighbors.end();
      for(; it!=iend; ++it)
      {
        if(*it == static_cast<std::size_t>(-1))
          continue;

        const Canvas_point& cq = this->canvas_points[*it];
        candidates.insert(&cq);
      }
      Base::construct_primal_elements_from_candidates(candidates);
    }
  }

  Grid_canvas(const std::string& canvas_str_,
              const std::string& seeds_str_,
              const Point_3& center_,
              const FT side_,
              const std::size_t points_per_side,
              const std::size_t max_seeds_n_,
              const Metric_field* mf_)
    :
      Base(canvas_str_, seeds_str_, max_seeds_n_, mf_),
      center(center_),
      side(side_),
      n(points_per_side),
      offset_x(center.x() - side/2.),
      offset_y(center.y() - side/2.),
      offset_z(center.z() - side/2.),
      step(side / (n-1)),
      sq_n(n*n)
  { }
};

}  // namespace Anisotropic_mesh_3
}  // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_GRID_CANVAS_H

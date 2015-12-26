#ifndef CGAL_ANISOTROPIC_MESH_3_CANVAS_CVT_OPTIMIZER_H
#define CGAL_ANISOTROPIC_MESH_3_CANVAS_CVT_OPTIMIZER_H

#include <CGAL/Canvas/canvas_config.h>
#include <CGAL/Canvas/canvas_enum.h>

#include <CGAL/helpers/metric_helper.h>

#include <CGAL/assertions.h>

#include <boost/array.hpp>
#include <Eigen/Dense>

#include <cstddef>
#include <set>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename Canvas>
class CVT_optimizer
{
public:
  typedef typename Canvas::Kernel                                Kernel;
  typedef typename Kernel::FT                                    FT;
  typedef typename Kernel::Point_3                               Point_3;
  typedef typename Kernel::Point_3                               TPoint_3;
  typedef typename Kernel::Segment_3                             Segment_3;
  typedef typename Kernel::Triangle_3                            Triangle_3;
  typedef typename Kernel::Tetrahedron_3                         Tetrahedron_3;

  typedef typename Canvas::Canvas_point                          Canvas_point;
  typedef typename Canvas::Metric                                Metric;

  typedef std::vector<std::pair<Point_3, FT> >                   Centroid_vector;
  typedef std::vector<Centroid_vector>                           Centroid_matrix;

  Canvas& canvas;
  std::size_t max_iter;

  mutable FT duration;

  FT tet_volume_in_metric(const Canvas_point& cp, const Canvas_point& cq,
                          const Canvas_point& cr, const Canvas_point& cs) const
  {
    FT fourth = 0.25;

    // could use another interpolation...
    Eigen::Matrix3d f = fourth * (cp.metric().get_transformation(),
                                  cq.metric().get_transformation(),
                                  cr.metric().get_transformation(),
                                  cs.metric().get_transformation());

    // transform the points
    TPoint_3 tcp = transform(f, cp.point());
    CGAL_assertion(tcp == old_transform(f, cp.point()));

    TPoint_3 tcr = transform(f, cq.point());
    TPoint_3 tcq = transform(f, cr.point());
    TPoint_3 tcs = transform(f, cs.point());

    typename Kernel::Compute_volume_3 o;
    return std::abs(o(tcp, tcq, tcr, tcs));
  }

  void add_to_canvas_centroids(const Canvas_point& cp, const Canvas_point& cq,
                               const Canvas_point& cr, const Canvas_point& cs,
                               Centroid_vector& seed_centroids) const
  {
    // this function inserts the centroid of pqrs in the canvas centroid matrix
    // FOR THE SEED CORRESPONDING TO p !!

    Point_3 centroid = CGAL::centroid(cp.point(), cq.point(), cr.point(), cs.point());
    FT vol = tet_volume_in_metric(cp, cq, cr, cs);
    CGAL_postcondition(vol != 0.);

    std::size_t seed_id = cp.closest_seed_id();
    CGAL_precondition(seed_id < canvas.seeds.size());
    seed_centroids.push_back(std::make_pair(centroid, vol));
  }

  void add_sub_tets_to_canvas_centroids(const std::size_t p,
                                        const std::size_t e1,
                                        const std::size_t e2,
                                        const std::size_t e3,
                                        const std::size_t e,
                                        Centroid_matrix& centroids) const
  {
    const Canvas_point& cp = canvas.get_point(p);
    const Canvas_point& ce1 = canvas.get_point(e1);
    const Canvas_point& ce2 = canvas.get_point(e2);
    const Canvas_point& ce3 = canvas.get_point(e3);
    const Canvas_point& c = canvas.get_point(e);

    std::size_t seed_id = cp.closest_seed_id();
    Centroid_vector& seed_centroids = centroids[seed_id];

    add_to_canvas_centroids(cp, ce1, ce2, c, seed_centroids);
    add_to_canvas_centroids(cp, ce2, ce3, c, seed_centroids);
    add_to_canvas_centroids(cp, ce3, ce1, c, seed_centroids);
  }

  void split_cell_multiple_colors(const std::size_t p, const std::size_t q,
                                  const std::size_t r, const std::size_t s,
                                  Centroid_matrix& centroids) const
  {
    // What we want to do here is to split the cell in fragments and assign
    // each fragment to the proper seed.
    // To do so, we compute some kind of mid point (for the geodesic distance)
    // for all the edges and for the full tet.
    // Ideally we would only compute the minimum amount of tets needed
    // to split the cell depending on which vertices have the same color etc.
    // However, it's very tedious:
    // let's say we have two different colors on a tet, the possibilities are:
    // - 1 vertex red, 3 vertices blue
    // - 2 vertices red, 2 vertices blue (and two different configurations here)
    // etc.

    // Anyway, multiply colored tetrahedra are a minority so it's not too costly
    // to simply split he cell in 12 tets that are based on the midpoints as soon
    // as there are more than 2 colors in the same tet.
    // If an edge has the same color at both extremities, then we just take
    // the Euclidean midpoint as midpoint.

    // Each vertex has 3 incident edges, which gives us 5 points: the vertex (v),
    // 3 edge midpoints and the Voronoi vertex (c) in the middle of the cell.
    // Each tet is therefore v-c + 2 edge points.

    // No need to bother with the degenerate case of c on a face, the volume will be 0
    // and the sub-tets won't count anyway...

    // would it better to split with points on the faces ? todo

    std::clock_t start = std::clock();
    std::size_t real_points_n = canvas.canvas_points.size();

    boost::array<std::size_t, 7> mid_pts; // 6 for the edges + 1 for the tet
    mid_pts[0] = canvas.compute_precise_Voronoi_vertex_on_edge(p, q);
    mid_pts[1] = canvas.compute_precise_Voronoi_vertex_on_edge(p, r);
    mid_pts[2] = canvas.compute_precise_Voronoi_vertex_on_edge(p, s);
    mid_pts[3] = canvas.compute_precise_Voronoi_vertex_on_edge(q, r);
    mid_pts[4] = canvas.compute_precise_Voronoi_vertex_on_edge(q, s);
    mid_pts[5] = canvas.compute_precise_Voronoi_vertex_on_edge(r, s);
    mid_pts[6] = canvas.compute_precise_Voronoi_vertex_in_tetrahedron(p, q, r, s);

    add_sub_tets_to_canvas_centroids(p, mid_pts[0], mid_pts[1], mid_pts[2],
                                     mid_pts[6], centroids);
    add_sub_tets_to_canvas_centroids(q, mid_pts[0], mid_pts[3], mid_pts[4],
                                     mid_pts[6], centroids);
    add_sub_tets_to_canvas_centroids(r, mid_pts[1], mid_pts[3], mid_pts[5],
                                     mid_pts[6], centroids);
    add_sub_tets_to_canvas_centroids(s, mid_pts[2], mid_pts[4], mid_pts[5],
                                     mid_pts[6], centroids);

    CGAL_assertion(real_points_n <= canvas.canvas_points.size());
    canvas.canvas_points.resize(real_points_n);

    duration += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  }

  void compute_canvas_centroids_for_all_seeds(Centroid_matrix& centroids,
                                              const std::size_t seed_id = -1) const
  {
    // This will only compile if you use a (Campen) triangulation canvas...
    // It should be easy enough to make it work for a grid canvas (simply split
    // the cubes in tets and do something like here) but it's kinda pointless since
    // I never use grids anymore... todo

    typedef typename Canvas::Tr                              Tr;
    typedef typename Tr::Vertex_handle                       Vertex_handle;
    typedef typename Tr::Finite_cells_iterator               Finite_cells_iterator;

    int one_c = 0, m_c = 0;
    double duration_one = 0., duration_m = 0.;

    std::cout << "cells : " << canvas.m_tr.number_of_finite_cells() << std::endl;

    Finite_cells_iterator cit = canvas.m_tr.finite_cells_begin();
    Finite_cells_iterator cend = canvas.m_tr.finite_cells_end();
    for(; cit!=cend; ++cit)
    {
      // ignore finite exterior cells (that is cit->info() == 0)
      if(cit->info() < 1)
        continue;

      if(seed_id != static_cast<std::size_t>(-1) &&
         canvas.canvas_points[cit->vertex(0)->info()].closest_seed_id() != seed_id &&
         canvas.canvas_points[cit->vertex(1)->info()].closest_seed_id() != seed_id &&
         canvas.canvas_points[cit->vertex(2)->info()].closest_seed_id() != seed_id &&
         canvas.canvas_points[cit->vertex(3)->info()].closest_seed_id() != seed_id)
          continue;

      std::set<std::size_t> colors;
      for(int i=0; i<4; ++i)
      {
        Vertex_handle vi = cit->vertex(i);
        colors.insert(canvas.canvas_points[vi->info()].closest_seed_id());
      }

      if(colors.size() == 1)
      {
        std::clock_t start = std::clock();
        std::size_t seed_id = canvas.canvas_points[cit->vertex(0)->info()].closest_seed_id();
        add_to_canvas_centroids(canvas.canvas_points[cit->vertex(0)->info()],
                                canvas.canvas_points[cit->vertex(1)->info()],
                                canvas.canvas_points[cit->vertex(2)->info()],
                                canvas.canvas_points[cit->vertex(3)->info()],
                                centroids[seed_id]);
        duration_one += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        one_c++;
      }
      else // colors.size() >= 2
      {
        std::clock_t start = std::clock();
        split_cell_multiple_colors(cit->vertex(0)->info(), cit->vertex(1)->info(),
                                   cit->vertex(2)->info(), cit->vertex(3)->info(),
                                   centroids);

        duration_m += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        m_c++;
      }
    }

    std::cout << "coloration (one/multiple): " << one_c << " " << m_c << std::endl;
    std::cout << duration_one << " " << duration_m << std::endl;
  }

  Point_3 compute_centroid(std::size_t seed_id,
                           const Centroid_vector& seed_centroids) const
  {
    CGAL_precondition(seed_id < seed_centroids.size());

    FT total_area = 0.;
    FT centroid_x = 0., centroid_y = 0., centroid_z = 0.;
    for(std::size_t i=0; i<seed_centroids.size(); ++i)
    {
      const std::pair<Point_3, FT>& centroid = seed_centroids[i];
      const Point_3& c = centroid.first;
      FT area = centroid.second;

      total_area += area;
      centroid_x += area * c.x();
      centroid_y += area * c.y();
      centroid_z += area * c.z();
    }

    CGAL_assertion(total_area != 0.);
    centroid_x /= total_area;
    centroid_y /= total_area;
    centroid_z /= total_area;

    return Point_3(centroid_x, centroid_y, centroid_z);
  }

  FT optimize_seed(std::size_t seed_id,
                   const Centroid_vector& seed_centroids)
  {
#if (VERBOSITY > 15)
    std::cout << "Optimizing seed: " << seed_id << std::endl;
#endif

    const Point_3 old_seed = canvas.seeds[seed_id];
    Point_3 new_seed = compute_centroid(seed_id, seed_centroids);

    canvas.seeds[seed_id] = new_seed;
    canvas.seeds.seeds_metrics[seed_id] = canvas.mf->compute_metric(new_seed);

    typename Kernel::Compute_squared_distance_3 sqd =
                                   Kernel().compute_squared_distance_3_object();
    FT displacement = sqd(old_seed, new_seed);

#if (VERBOSITY > 17)
    std::cout << "centroid with grid triangles " << new_seed << std::endl;
    std::cout << "seed: " << seed_id << " displacement: " << displacement << std::endl;
#endif

    return displacement;
  }

  void optimize_seeds(const std::string str_base)
  {
    if(max_iter == 0)
      return;

#if (VERBOSITY > 10)
    std::cout << "Optimizing all seeds" << std::endl;
#endif

    for(std::size_t i=0; i<canvas.canvas_points.size(); ++i)
      CGAL_assertion(canvas.canvas_points[i].state() == KNOWN);

    bool is_optimized = false;
    std::size_t counter = 0;

    do
    {
#if (VERBOSITY > 12)
      std::cout << "Optimizing: iteration " << counter << std::endl;
#endif

      Centroid_matrix centroids(canvas.seeds.size());
      compute_canvas_centroids_for_all_seeds(centroids);

      std::cout << "duration: " << duration << std::endl;
      duration = 0;

      FT cumulated_displacement = 0;
      for(std::size_t i=0; i<canvas.seeds.size(); ++i)
        cumulated_displacement += optimize_seed(i, centroids[i]);

#if (VERBOSITY > 15)
       std::cout << "cumulated displacement : " << cumulated_displacement << std::endl;
#endif

       canvas.reset();
       canvas.locate_seeds_on_canvas();
       canvas.paint(false/*refining*/);

       std::ostringstream opti_out;
       opti_out << "optimized_" << str_base << "_tr_" << counter << std::ends;
       canvas.output_canvas_data_and_primal(opti_out.str().c_str());

       is_optimized = (++counter > max_iter || cumulated_displacement < 1e-10); //fixme hardcoded
    } while(!is_optimized);
  }

  CVT_optimizer() { }

  CVT_optimizer(Canvas& canvas_, std::size_t max_iter_ = 1000)
    :
      canvas(canvas_),
      max_iter(max_iter_),
      duration(0)
  { }
};

}  // namespace Anisotropic_mesh_3
}  // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_CANVAS_CVT_OPTIMIZER_H

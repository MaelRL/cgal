#ifndef CGAL_ANISOTROPIC_MESH_3_CANVAS_CVT_OPTIMIZER_SURF_H
#define CGAL_ANISOTROPIC_MESH_3_CANVAS_CVT_OPTIMIZER_SURF_H

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

template<typename Canvas, typename C3t3>
class CVT_surf_optimizer
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
  typedef typename Canvas::Vector3d                              Vector3d;

  typedef typename Canvas::Tr                              Tr;
  typedef typename Tr::Vertex_handle                       Vertex_handle;

  typedef typename C3t3::Facets_in_complex_iterator        Facets_in_complex_iterator;

  typedef std::vector<std::pair<Point_3, FT> >                   Centroid_vector;
  typedef std::vector<Centroid_vector>                           Centroid_matrix;

  typedef boost::array<Vector3d, 3>                              Basis;
  typedef typename std::vector<Basis>                            Bases;

  Canvas& canvas;
  std::size_t max_iter;

  Bases tangent_bases;

  mutable FT duration;

  void compute_tangent_bases()
  {
    std::size_t css = canvas.seeds.size();

    for(std::size_t i=0; i<css; ++i)
    {
      Point_3 s = canvas.seeds[i];

      // locate seed
      FT min_d = FT_inf;
      Vertex_handle min_v;

      typename C3t3::Triangulation::Finite_vertices_iterator vit =
        canvas.m_c3t3.triangulation().finite_vertices_begin();
      typename C3t3::Triangulation::Finite_vertices_iterator vend =
        canvas.m_c3t3.triangulation().finite_vertices_end();
      for(; vit!=vend; ++vit)
      {
        FT sqd = CGAL::squared_distance(s, vit->point());
        if(sqd < min_d)
        {
          min_d = sqd;
          min_v = vit;
        }
      }

      CGAL_assertion(min_d < 1e-10);
      CGAL_assertion((normal.norm() - 1) < 1e-10 );

      Vector3d normal = canvas.m_vertex_normals[min_v];

      // compute two opposites
      Vector3d t1, t2;
      t1(0) = normal(1); t1(1) = - normal(0); t1(2) = 0;
      t2 = normal.cross(t1);

      t1.normalize();
      t2.normalize();

      tangent_bases[i][0] = normal;
      tangent_bases[i][1] = t1;
      tangent_bases[i][2] = t2;
    }
    std::cout << "computed tangent spaces" << std::endl;
  }

  FT tri_volume_in_metric(const Canvas_point& cp, const Canvas_point& cq,
                          const Canvas_point& cr) const
  {
    FT third = 1./3.;

    // could use another interpolation...
    Eigen::Matrix3d f = third * (cp.metric().get_transformation(),
                                 cq.metric().get_transformation(),
                                 cr.metric().get_transformation());

    // transform the points
    TPoint_3 tcp = transform(f, cp.point());
    CGAL_assertion(tcp == old_transform(f, cp.point()));
    TPoint_3 tcr = transform(f, cq.point());
    TPoint_3 tcq = transform(f, cr.point());

    typename Kernel::Compute_area_3 o;
    return std::abs(o(tcp, tcq, tcr));
  }

  void add_to_canvas_centroids(const Canvas_point& cp, const Canvas_point& cq,
                               const Canvas_point& cr,
                               Centroid_vector& seed_centroids)
  {
    // this function inserts the centroid of pqr in the canvas centroid matrix
    // FOR THE SEED CORRESPONDING TO p !!

    Point_3 centroid = CGAL::centroid(cp.point(), cq.point(), cr.point());
    FT area = tri_volume_in_metric(cp, cq, cr);
    CGAL_postcondition(area != 0.);

    std::size_t seed_id = cp.closest_seed_id();
    CGAL_precondition(seed_id < canvas.seeds.size());

    const Point_3& seed = canvas.seeds[seed_id];

    // lift to the tangent plane of p
    Basis& tangent_basis = tangent_bases[seed_id];

    Vector3d orig_to_centroid;
    orig_to_centroid(0) = centroid.x() - seed.x();
    orig_to_centroid(1) = centroid.y() - seed.y();
    orig_to_centroid(2) = centroid.z() - seed.z();

    FT c1 = tangent_basis[1].dot(orig_to_centroid);
    FT c2 = tangent_basis[2].dot(orig_to_centroid);

    // position in 3D space is orign + coeff
    FT x = seed.x() + c1 * tangent_basis[1](0) + c2 * tangent_basis[2](0);
    FT y = seed.y() + c1 * tangent_basis[1](1) + c2 * tangent_basis[2](1);
    FT z = seed.z() + c1 * tangent_basis[1](2) + c2 * tangent_basis[2](2);

    Point_3 lifted_centroid(x, y, z);

    seed_centroids.push_back(std::make_pair(lifted_centroid, area));
  }

  void compute_canvas_centroids_for_all_seeds(Centroid_matrix& centroids,
                                              const std::size_t seed_id = -1)
  {
    // This will only compile if you use a (Campen) triangulation canvas...
    // It should be easy enough to make it work for a grid canvas (simply split
    // the cubes in tets and do something like here) but it's kinda pointless since
    // I never use grids anymore... todo

    int one_c = 0, m_c = 0;
    double duration_one = 0., duration_m = 0.;

    Facets_in_complex_iterator fit = canvas.m_c3t3.facets_in_complex_begin();
    Facets_in_complex_iterator fend = canvas.m_c3t3.facets_in_complex_end();
    for(; fit!=fend; ++fit)
    {
      Vertex_handle v1 = fit->first->vertex((fit->second + 1)%4);
      Vertex_handle v2 = fit->first->vertex((fit->second + 2)%4);
      Vertex_handle v3 = fit->first->vertex((fit->second + 3)%4);

      if(seed_id != static_cast<std::size_t>(-1) &&
         canvas.canvas_points[v1->info()].closest_seed_id() != seed_id &&
         canvas.canvas_points[v2->info()].closest_seed_id() != seed_id &&
         canvas.canvas_points[v3->info()].closest_seed_id() != seed_id)
          continue;

      std::set<std::size_t> colors;

      colors.insert(canvas.canvas_points[v1->info()].closest_seed_id());
      colors.insert(canvas.canvas_points[v2->info()].closest_seed_id());
      colors.insert(canvas.canvas_points[v3->info()].closest_seed_id());

      if(colors.size() == 1)
      {
        std::clock_t start = std::clock();
        std::size_t seed_id = canvas.canvas_points[v1->info()].closest_seed_id();
        add_to_canvas_centroids(canvas.canvas_points[v1->info()],
                                canvas.canvas_points[v2->info()],
                                canvas.canvas_points[v3->info()],
                                centroids[seed_id]);
        duration_one += ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        one_c++;
      }
      else // colors.size() >= 2
      {
        std::clock_t start = std::clock();
//        split_cell_multiple_colors(v1->info(), v2->info(), v3->info(), centroids);

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

  Point_3 project_new_seed(const Point_3& new_seed)
  {
    Point_3 s;
    FT min_d = FT_inf;

    typename C3t3::Triangulation::Finite_vertices_iterator vit =
      canvas.m_c3t3.triangulation().finite_vertices_begin();
    typename C3t3::Triangulation::Finite_vertices_iterator vend =
      canvas.m_c3t3.triangulation().finite_vertices_end();
    for(; vit!=vend; ++vit)
    {
      FT sqd = CGAL::squared_distance(new_seed, vit->point());
      if(sqd < min_d)
      {
        min_d = sqd;
        s = vit->point();
      }
    }
    return s;
  }

  FT optimize_seed(std::size_t seed_id,
                   const Centroid_vector& seed_centroids)
  {
#if (VERBOSITY > 15)
    std::cout << "Optimizing seed: " << seed_id << std::endl;
#endif

    const Point_3 old_seed = canvas.seeds[seed_id];
    Point_3 new_seed = compute_centroid(seed_id, seed_centroids);

    // project
    Point_3 projected_new_seed = project_new_seed(new_seed);

    canvas.seeds[seed_id] = projected_new_seed;
    canvas.seeds.seeds_metrics[seed_id] = canvas.mf->compute_metric(projected_new_seed);

    typename Kernel::Compute_squared_distance_3 sqd =
                                   Kernel().compute_squared_distance_3_object();
    FT displacement = sqd(old_seed, projected_new_seed);

#if (VERBOSITY > 17)
    std::cout << "centroid with grid triangles " << projected_new_seed << std::endl;
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

      compute_tangent_bases();

      Centroid_matrix centroids(canvas.seeds.size());
      compute_canvas_centroids_for_all_seeds(centroids);

      std::cout << "duration: " << duration << std::endl;
      duration = 0;

      FT cumulated_displacement = 0;
      for(std::size_t i=0; i<canvas.seeds.size(); ++i)
        cumulated_displacement += optimize_seed(i, centroids[i]);

       std::cout << "cumulated displacement : " << cumulated_displacement << std::endl;

       canvas.reset();
       canvas.locate_seeds_on_canvas();
       canvas.paint(false/*refining*/);

       std::ostringstream opti_out;
       opti_out << "optimized_" << str_base << "_tr_" << counter << std::ends;
       canvas.output_canvas_data_and_primal(opti_out.str().c_str());

       is_optimized = (++counter > max_iter || cumulated_displacement < 1e-10); //fixme hardcoded
    } while(!is_optimized);
  }

  CVT_surf_optimizer() { }

  CVT_surf_optimizer(Canvas& canvas_, std::size_t max_iter_ = 1000)
    :
      canvas(canvas_),
      max_iter(max_iter_),
      tangent_bases(canvas.seeds.size()),
      duration(0)
  { }
};

}  // namespace Anisotropic_mesh_3
}  // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_CANVAS_CVT_OPTIMIZER_SURF_H

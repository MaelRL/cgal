#ifndef CGAL_ANISOTROPIC_MESH_3_CANVAS_SEEDS_H
#define CGAL_ANISOTROPIC_MESH_3_CANVAS_SEEDS_H

#include <CGAL/Canvas/canvas_config.h>

#include <CGAL/Metric.h>

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <istream>
#include <string>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename Canvas>
class Canvas_seeds
{
public:
  typedef typename Canvas::Kernel                              Kernel;
  typedef typename Kernel::FT                                  FT;
  typedef typename Kernel::Point_3                             Point_3;

  typedef Metric_base<Kernel>                                  Metric;

  const std::string seeds_str;
  std::size_t max_seeds_n;

  // the seeds and their corresponding metrics
  std::vector<Point_3> seeds;
  std::vector<Metric> seeds_metrics;

  const Canvas& canvas;

  const Point_3& operator[](const int i) const { return seeds[i]; }
  Point_3& operator[] (const int i) { return seeds[i]; }
  std::size_t size() const { return seeds.size(); }

  const std::vector<Metric>& metrics() const { return seeds_metrics; }

  std::size_t insert_new_seed(const FT x, const FT y, const FT z)
  {
#ifdef ANISO_GEO_FILTER_SEEDS_OUTSIDE_CANVAS_BBOX
    if(canvas.is_point_outside_canvas(x, y, z))
    {
  #if (VERBOSITY > 1)
      std::cout << "filtered : " << x << " " << y << " " << z << std::endl;
  #endif
      return seeds.size();
    }
#endif

#if (VERBOSITY > 0)
    std::cout << "added new seed: " << x << " " << y << " " << z;
    std::cout << " (" << seeds.size() << ")" << std::endl;
#endif

    seeds.push_back(Point_3(x, y, z));
    seeds_metrics.push_back(canvas.mf->compute_metric(seeds.back()));
    return seeds.size();
  }

  std::size_t build_seeds()
  {
    std::ifstream in(seeds_str.c_str());
    std::string word;
    std::size_t useless, nv, dim;
    FT r_x, r_y, r_z;

    in >> word >> useless; // MeshVersionFormatted i
    in >> word >> dim; // Dimension d
    in >> word >> nv;
    std::cout << "seeds nv: " << nv << std::endl;
    CGAL_assertion(dim == 3);

    std::size_t min_nv = (std::min)(nv, max_seeds_n);

    seeds.reserve(min_nv);
    seeds_metrics.reserve(min_nv);

    for(std::size_t i=0; i<nv; ++i)
    {
      in >> r_x >> r_y >> r_z >> useless;
      insert_new_seed(r_x, r_y, r_z);

      if(seeds.size() >= max_seeds_n)
        break;
    }
#if (VERBOSITY > 0)
    std::cout << "seeds: " << seeds.size() << std::endl;
#endif
    return seeds.size();
  }

  void initialize_seeds()
  {
    max_seeds_n = build_seeds();
    CGAL_assertion(max_seeds_n > 0 && "No seed in domain..." );
  }

  Canvas_seeds(const Canvas& _canvas,
               const std::string _seeds_str,
               std::size_t _max_seeds_n)
    :
      seeds_str(_seeds_str),
      max_seeds_n(_max_seeds_n),
      seeds(),
      seeds_metrics(),
      canvas(_canvas)
  { }
};

}  // namespace Anisotropic_mesh_3
}  // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_CANVAS_SEEDS_H

#ifndef CGAL_ANISOTROPIC_MESH_3_CANVAS_MESHER_H
#define CGAL_ANISOTROPIC_MESH_3_CANVAS_MESHER_H

#include <CGAL/Canvas/canvas_config.h>
#include <CGAL/Canvas/refine_queue.h>

#include <CGAL/intersections.h>
#include <CGAL/Tetrahedron_3_Tetrahedron_3_do_intersect.h>
#include <CGAL/Tetrahedron_3_Segment_3_do_intersect.h>

#include <CGAL/assertions.h>

#include <algorithm>
#include <cmath>
#include <ctime>
#include <cstddef>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename Canvas>
class Canvas_mesher
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
  typedef typename Canvas::Primal_triangle                       Primal_triangle;
  typedef typename Canvas::Primal_tetrahedron                    Primal_tetrahedron;
  typedef typename Primal_triangle::BSimplex                     TrSimplex;
  typedef typename Primal_tetrahedron::BSimplex                  TSimplex;
  typedef typename Canvas::Primal_triangles_container            Primal_triangles_container;
  typedef typename Canvas::Primal_tetrahedra_container           Primal_tetrahedra_container;

  typedef Canvas_refinement_queue<Canvas, 3>                     Face_queues;
  typedef typename Face_queues::Queue_entry_iterator             F_Queue_it;
  typedef Canvas_refinement_queue<Canvas, 4>                     Tet_queues;
  typedef typename Tet_queues::Queue_entry_iterator              T_Queue_it;

  Canvas& canvas;
  Face_queues fqueues;
  Tet_queues tqueues;
  std::size_t max_refine_n;

  FT compute_volume(const Point_3& p, const Point_3& q, const Point_3& r) const
  {
    typename Kernel::Compute_area_3 o;
    return std::abs(o(p, q, r));
  }

  FT compute_volume(const Point_3& p, const Point_3& q, const Point_3& r, const Point_3& s) const
  {
    typename Kernel::Compute_volume_3 o;
    return std::abs(o(p, q, r, s));
  }

  template<typename PSimplex>
  FT compute_volume(const PSimplex& s)
  {
    if(s.size() == 3)
      return compute_volume(canvas.seeds[s[0]], canvas.seeds[s[1]],
                            canvas.seeds[s[2]]);
    else // if(s.size() == 4)
      return compute_volume(canvas.seeds[s[0]], canvas.seeds[s[1]],
                            canvas.seeds[s[2]], canvas.seeds[s[3]]);
  }

  FT quality(const Point_3& p, const Point_3& q, const Point_3& r) const
  {
    typename Kernel::Compute_squared_distance_3 sqd =
                                    Kernel().compute_squared_distance_3_object();

    FT alpha = 4.*std::sqrt(3.);
    FT A = compute_volume(p, q, r);
    FT a = sqd(p, q);
    FT b = sqd(p, r);
    FT c = sqd(q, r);
    FT quality = alpha*A/(a+b+c);
    return quality;
  }

  FT compute_quality(const TrSimplex& tr) const
  {
    const Point_3& p0 = canvas.seeds[tr[0]];
    const Point_3& p1 = canvas.seeds[tr[1]];
    const Point_3& p2 = canvas.seeds[tr[2]];

    const Metric& m0 = canvas.seeds.seeds_metrics[tr[0]];
    const Metric& m1 = canvas.seeds.seeds_metrics[tr[1]];
    const Metric& m2 = canvas.seeds.seeds_metrics[tr[2]];

    const TPoint_3& tp0_0 = m0.transform(p0);
    const TPoint_3& tp1_0 = m0.transform(p1);
    const TPoint_3& tp2_0 = m0.transform(p2);

    const TPoint_3& tp0_1 = m1.transform(p0);
    const TPoint_3& tp1_1 = m1.transform(p1);
    const TPoint_3& tp2_1 = m1.transform(p2);

    const TPoint_3& tp0_2 = m2.transform(p0);
    const TPoint_3& tp1_2 = m2.transform(p1);
    const TPoint_3& tp2_2 = m2.transform(p2);

    FT qual = (std::min)((std::min)(quality(tp0_0, tp1_0, tp2_0),
                                    quality(tp0_1, tp1_1, tp2_1)),
                                    quality(tp0_2, tp1_2, tp2_2));
    return qual;
  }

  FT quality(const Point_3& p, const Point_3& q, const Point_3& r, const Point_3& s) const
  {
    typename Kernel::Compute_squared_distance_3 sqd =
                                    Kernel().compute_squared_distance_3_object();

    FT alpha = 216.*std::sqrt(3.);
    FT V = compute_volume(p, q, r, s);
    FT a = compute_volume(p, q, r);
    FT b = compute_volume(p, q, s);
    FT c = compute_volume(p, r, s);
    FT d = compute_volume(q, r, s);
    FT denom = (a+b+c+d);
    FT quality = alpha*V*V/(denom*denom*denom);
    return quality;
  }

  FT compute_quality(const TSimplex& tet) const
  {
    const Point_3& p0 = canvas.seeds[tet[0]];
    const Point_3& p1 = canvas.seeds[tet[1]];
    const Point_3& p2 = canvas.seeds[tet[2]];
    const Point_3& p3 = canvas.seeds[tet[3]];

    const Metric& m0 = canvas.seeds.seeds_metrics[tet[0]];
    const Metric& m1 = canvas.seeds.seeds_metrics[tet[1]];
    const Metric& m2 = canvas.seeds.seeds_metrics[tet[2]];
    const Metric& m3 = canvas.seeds.seeds_metrics[tet[3]];

    const TPoint_3& tp0_0 = m0.transform(p0);
    const TPoint_3& tp1_0 = m0.transform(p1);
    const TPoint_3& tp2_0 = m0.transform(p2);
    const TPoint_3& tp3_0 = m0.transform(p3);

    const TPoint_3& tp0_1 = m1.transform(p0);
    const TPoint_3& tp1_1 = m1.transform(p1);
    const TPoint_3& tp2_1 = m1.transform(p2);
    const TPoint_3& tp3_1 = m1.transform(p3);

    const TPoint_3& tp0_2 = m2.transform(p0);
    const TPoint_3& tp1_2 = m2.transform(p1);
    const TPoint_3& tp2_2 = m2.transform(p2);
    const TPoint_3& tp3_2 = m2.transform(p3);

    const TPoint_3& tp0_3 = m3.transform(p0);
    const TPoint_3& tp1_3 = m3.transform(p1);
    const TPoint_3& tp2_3 = m3.transform(p2);
    const TPoint_3& tp3_3 = m3.transform(p3);

    FT qual = (std::min)((std::min)((std::min)(quality(tp0_0, tp1_0, tp2_0, tp3_0),
                                               quality(tp0_1, tp1_1, tp2_1, tp3_1)),
                                               quality(tp0_2, tp1_2, tp2_2, tp3_2)),
                                               quality(tp0_3, tp1_3, tp2_3, tp3_3));
    return qual;
  }

  void clean_refinement_queue()
  {
    // should be called during the refinement process when we overwrite colors
    // and some simplices disappear (thus need to be cleaned from the queues)...
    // Tedious... todo
  }

  template<typename PSimplex, typename Queue>
  void test_primal_simplex(const PSimplex& primal_simplex, Queue& queue)
  {
    // todo this whole thing can be extended to a primal simplex of any dimension
    const Canvas_point* cp = primal_simplex.dual_point();
    const PSimplex& s = primal_simplex.simplex();

    // todo criteria in a criteria class...
    FT max_distortion = 1.; // <= 1 means unused
    FT max_size = 0.01; // <= 0 means unused
    bool intersection_ref = false;
    FT min_qual = 0.; // <= 0 means unsused

    // distortion
    if(max_distortion > 1.)
    {
      FT gamma = cp->distortion_to_seed();
      if(gamma > max_distortion)
      {
        queue.push(primal_simplex, gamma, 0/*queue_id*/);
        return;
      }
    }

    // size
    if(max_size > 0.)
    {
      FT cpd = cp->distance_to_closest_seed();
      if(cpd > max_size)
      {
        queue.push(primal_simplex, cpd, 1/*queue_id*/);
        return;
      }
    }

    // intersection
    if(intersection_ref)
    {
      if(primal_simplex.m_is_intersected)
      {
        FT vol = compute_volume(s);
        queue.push(primal_simplex, vol, 2/*queue_id*/);
        return;
      }
    }

    // quality fixme quality should be radius edge ratio or something else
    if(min_qual > 0.)
    {
      FT qual = compute_quality(primal_simplex.simplex());
      if(qual < min_qual)
      {
        CGAL_assertion(qual > 1e-17);
        // 1./qual to refine the worst element first
        queue.push(primal_simplex, 1./qual, 3/*queue_id*/);
        return;
      }
    }
  }

  void fill_refinement_queue(const std::size_t i)
  {
    // version where we only consider new simplices todo
    // this would be used with a smart version of the primal shenanigans that
    // cleans the queue while spreading (while refining)
  }

  void build_refinement_queue()
  {
#if (VERBOSITY > 8)
    std::cout << "build refinement queue" << std::endl;
#endif

    if(canvas.primal_tetrahedra.empty() && canvas.primal_triangles.empty())
      canvas.compute_primal();

    fqueues.clear();
    tqueues.clear();

    // loop border triangles
    for(typename Primal_triangles_container::iterator it = canvas.primal_triangles.begin();
                                                      it != canvas.primal_triangles.end();
                                                      ++it)
    {
      const Primal_triangle& primal_tr = *it;
      test_primal_simplex(primal_tr, fqueues);
    }

    // make sure we've computed the potential intersections between tetrahedra
    canvas.detect_tetrahedra_self_intersections();

    // loop tetrahedron & test them
    for(typename Primal_tetrahedra_container::iterator it = canvas.primal_tetrahedra.begin();
                                                       it != canvas.primal_tetrahedra.end();
                                                       ++it)
    {
      const Primal_tetrahedron& primal_tet = *it;
      test_primal_simplex(primal_tet, tqueues);
    }
  }

  bool get_next_refinement_point(const Canvas_point*& new_seed) const
  {
    F_Queue_it fe;
    T_Queue_it te;
    if(!fqueues.top(fe))
    {
      if(!tqueues.top(te))
      {
        std::cerr << "Couldn't find a ref point..." << std::endl;
        return false;
      }
      else
      {
        std::cout << "selected ref point from cell queue" << std::endl;
        new_seed = te->dual_point();
      }
    }
    else
    {
      std::cout << "selected ref point from face queue" << std::endl;
      new_seed = fe->dual_point();
    }

    return true;
  }

  void refine(const Point_3& new_seed)
  {
#if (VERBOSITY > 5)
    std::cout << "refine with : " << new_seed << std::endl;
    std::cout << "now " << canvas.seeds.size()+1 << " seeds" << std::endl;
#endif

    std::size_t old_seed_size =  canvas.seeds.size();
    canvas.seeds.max_seeds_n = canvas.seeds.insert_new_seed(new_seed.x(),
                                                            new_seed.y(),
                                                            new_seed.z());
    CGAL_assertion(old_seed_size != canvas.seeds.size());

#ifdef USE_FULL_REBUILD
    canvas.reset();
    canvas.locate_seeds_on_canvas();
#else
    // we can't spread from the new seed if all states are 'KNOWN'
    canvas.refresh_canvas_point_states();
    canvas.locate_and_initialize(new_seed, canvas.seeds.size()-1);
#endif

    canvas.paint(true/*refining*/);
    if(!fqueues.pop())
      tqueues.pop();
  }

  FT refine()
  {
    std::clock_t start = std::clock();

    if(max_refine_n < 1)
      return 0;

    std::ostringstream out_init;
    out_init << "ref_" << canvas.seeds.size();
    canvas.output_canvas_data_and_primal(out_init.str());

    for(std::size_t i=1; i<=max_refine_n; ++i)
    {
      build_refinement_queue();
      fqueues.print_queues();
      tqueues.print_queues();

      const Canvas_point* new_seed = NULL;
      if(!get_next_refinement_point(new_seed))
        break;

      refine(new_seed->point());

      if(i%100 == 0)
      {
        std::ostringstream out;
        out << "ref_" << canvas.seeds.size();
        canvas.output_canvas_data_and_primal(out.str());
      }
    }

    canvas.set_points_states_to_known();

    double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cerr << "End refinement: " << duration << std::endl;
    return duration;
  }

  Canvas_mesher(Canvas& canvas_, std::size_t max_refine_n_)
    :
      canvas(canvas_),
      fqueues(),
      tqueues(),
      max_refine_n(max_refine_n_)
  { }
};

}  // namespace Anisotropic_mesh_3
}  // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_CANVAS_MESHER_H

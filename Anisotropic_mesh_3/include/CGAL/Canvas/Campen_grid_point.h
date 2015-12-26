#ifndef CGAL_ANISOTROPIC_MESH_3_CANVAS_GRID_POINT_H
#define CGAL_ANISOTROPIC_MESH_3_CANVAS_GRID_POINT_H

#include <CGAL/Canvas/canvas_point.h>

#include <CGAL/assertions.h>

#include <Eigen/Dense>
#include <boost/array.hpp>

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_3
{


template<typename K, typename CPoint, typename MF>
class Grid_canvas;

template<typename K, typename MF>
class Campen_grid_point :
    public Canvas_point<K, Grid_canvas<K, Campen_grid_point<K, MF>, MF> >
{
  typedef Campen_grid_point<K, MF>                                   Self;
public:
  typedef boost::array<std::size_t, 6>                               Neighbors;
  typedef Grid_canvas<K, Self, MF>                                   Canvas;
  typedef Canvas_point<K, Canvas>                                    Base;

  typedef typename Base::FT                                          FT;
  typedef typename Base::Point_3                                     Point_3;
  typedef typename Base::Metric                                      Metric;
  typedef typename Base::Vector3d                                    Vector3d;

  Neighbors neighbors;

  bool compute_closest_seed(std::size_t anc)
  {
    // returns true if we improved the distance
    CGAL_assertion(anc != static_cast<std::size_t>(-1) &&
                   anc < this->canvas()->canvas_points.size());
    CGAL_assertion(this->canvas()->get_point(anc).state() == KNOWN);

    const int k = 8; // depth of the ancestor edge
    FT d = FT_inf;

    // the path is 'this' to ancestor1, ancestor1 to ancestor2, etc.
    // stored as 'this', ancestor1, ancestor2, etc.
    boost::array<std::size_t, k+1> ancestor_path;
    for(int i=1; i<k+1; ++i)
      ancestor_path[i] = -1;
    ancestor_path[0] = this->index();
    ancestor_path[1] = anc;

    boost::array<Eigen::Matrix3d, k> path_metrics;
    for(int i=0; i<k; ++i)
    {
      if(i >= 1)
      {
        if(this->canvas()->get_point(ancestor_path[i]).ancestor() ==
                                                   static_cast<std::size_t>(-1))
          break;

        ancestor_path[i+1] = this->canvas()->get_point(ancestor_path[i]).ancestor();
      }

      std::size_t e0 = ancestor_path[i];
      std::size_t e1 = ancestor_path[i+1];

      const Metric& m0 = this->canvas()->get_point(e0).metric();
      const Metric& m1 = this->canvas()->get_point(e1).metric();

      path_metrics[i] = get_interpolated_transformation(m0, m1);
    }

    std::size_t curr_anc = anc;
    for(int i=1; i<=k; ++i)
    {
      const Base& curr_ancestor = this->canvas()->get_point(curr_anc);

      Vector3d ancestor_edge;
      ancestor_edge(0) = this->point().x() - curr_ancestor.point().x();
      ancestor_edge(1) = this->point().y() - curr_ancestor.point().y();
      ancestor_edge(2) = this->point().z() - curr_ancestor.point().z();
      FT ancestor_edge_length = ancestor_edge.norm();
      Vector3d normalized_anc_edge = ancestor_edge / ancestor_edge_length;

      // compute the distance for the current depth (i)
      FT dist_to_ancestor = 0.;
      for(int j=0; j<i; ++j) // we add a part for each edge in the path
      {
        // get the metric for the current edge
        const Base& e0 = this->canvas()->get_point(ancestor_path[j]);
        const Base& e1 = this->canvas()->get_point(ancestor_path[j+1]);

        Vector3d curr_edge;
        curr_edge(0) = e0.point().x() - e1.point().x();
        curr_edge(1) = e0.point().y() - e1.point().y();
        curr_edge(2) = e0.point().z() - e1.point().z();

        // interpolate between both metric and transform the normalized edge
        // then we have (transformed_edge).norm() = || e ||_M = sqrt(e^t M e)
        const Eigen::Matrix3d& f = path_metrics[j];

        Vector3d transformed_curr_edge = f*normalized_anc_edge;

        FT sp = curr_edge.dot(normalized_anc_edge);
        FT l = transformed_curr_edge.norm(); // length of the normalized anc edge in the metric

        dist_to_ancestor += sp * l;
      }
      dist_to_ancestor = (std::max)(dist_to_ancestor, 0.);

      // add ancestor edge length to the distance at that ancestor
      FT dist_at_anc = curr_ancestor.distance_to_closest_seed();
      FT new_d = dist_at_anc + dist_to_ancestor;

      if(new_d < d)
        d = new_d;

      if(!curr_ancestor.ancestor()) // can't go any farther up in the ancestor tree
        break;

      curr_anc = curr_ancestor.ancestor();
    }

    if(d < this->distance_to_closest_seed())
    {
      this->ancestor() = anc;
      this->distance_to_closest_seed() = d;
      this->closest_seed_id() = this->canvas()->get_point(anc).closest_seed_id();
      return true;
    }

    return false;
  }

  PQ_state update_neighbors_distances(std::vector<Self*>& trial_pq)
  {
    // consider all the neighbors of a KNOWN point and compute their distance to 'this'
#if (VERBOSITY > 15)
    std::cout << "update neighbors of " << this->index() << std::endl;
#endif
    CGAL_assertion(this->state() == KNOWN);

    PQ_state pqs_ret = NOTHING_TO_DO;
    typename Neighbors::const_iterator it = neighbors.begin(),
                                       iend = neighbors.end();
    for(; it!=iend; ++it)
    {
      if(*it == static_cast<std::size_t>(-1))
        continue;

      Self& cp = this->canvas()->get_point(*it);
      if(cp.state() == KNOWN)
        continue; // dual_shenanigans(cp);
      else if(cp.state() == TRIAL)
      {
        // note that we don't insert in trial_pq since it's already in
        if(cp.compute_closest_seed(this->index()))
          pqs_ret = REBUILD_TRIAL;
      }
      else // cp.state == FAR
      {
        CGAL_assertion(cp.state() == FAR);

        // note that cp.distance_to_closest_seed is not necessarily FT_inf here :
        // if we're refining, we've assigned FAR to all points after inserting a new
        // seed, therefore we must verify that compute_closest_seed is an update
        // before inserting it in the trial_queue
        if(cp.compute_closest_seed(this->index()))
        {
          cp.state() = TRIAL;
          trial_pq.push_back(&cp);
          std::push_heap(trial_pq.begin(), trial_pq.end(),
                         Canvas_point_comparer<Self>());
        }
      }
    }
    return pqs_ret;
  }

  Campen_grid_point(const Point_3& point_,
                    const std::size_t index_,
                    Canvas* canvas_)
    :
      Base(point_, index_, canvas_),
      neighbors()
  {
    for(std::size_t i=0; i<neighbors.size(); ++i)
      neighbors[i] = -1;
  }
};

} // namespace Anisotropic_mesh_3
} // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_CANVAS_GRID_POINT_H

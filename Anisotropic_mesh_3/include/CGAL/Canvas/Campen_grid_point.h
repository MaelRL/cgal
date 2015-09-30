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

template<typename K>
class Campen_canvas_point :
    public Canvas_point<K>
{
public:
  typedef boost::array<Campen_canvas_point*, 6>                      Neighbors;
  typedef Canvas_point<K>                                            Base;

  typedef typename Base::FT                                          FT;
  typedef typename Base::Point_3                                     Point_3;
  typedef typename Base::Metric                                      Metric;
  typedef typename Base::Vector3d                                    Vector3d;

  Neighbors neighbors;

  bool compute_closest_seed(const Base* anc)
  {
    // returns true if we improved the distance
    CGAL_assertion(anc->state() == KNOWN);

    const int k = 8; // depth of the ancestor edge
    FT d = FT_inf;

    // the path is 'this' to ancestor1, ancestor1 to ancestor2, etc.
    // stored as 'this', ancestor1, ancestor2, etc.
    boost::array<const Base*, k+1> ancestor_path;
    for(int i=1; i<k+1; ++i)
      ancestor_path[i] = NULL;
    ancestor_path[0] = this;

    const Base* curr_ancestor = anc;
    for(int i=1; i<=k; ++i)
    {
      // add the new segment to the ancestor path
      ancestor_path[i] = curr_ancestor;

      Vector3d ancestor_edge;
      ancestor_edge(0) = this->point().x() - curr_ancestor->point().x();
      ancestor_edge(1) = this->point().y() - curr_ancestor->point().y();
      ancestor_edge(2) = this->point().z() - curr_ancestor->point().z();
      FT ancestor_edge_length = ancestor_edge.norm();
      Vector3d normalized_anc_edge = ancestor_edge / ancestor_edge_length;

      // compute the distance for the current depth (i)
      FT dist_to_ancestor = 0.;
      for(int j=0; j<i; ++j) // we add a part for each edge in the path
      {
        // get the metric for the current edge
        const Base* e0 = ancestor_path[j];
        const Base* e1 = ancestor_path[j+1];

        CGAL_assertion(e0 && e1);

        const Metric& m0 = e0->metric();
        const Metric& m1 = e1->metric();

        Vector3d curr_edge;
        curr_edge(0) = e0->point().x() - e1->point().x();
        curr_edge(1) = e0->point().y() - e1->point().y();
        curr_edge(2) = e0->point().z() - e1->point().z();

        // interpolate between both metric and transform the normalized edge
        // then we have (transformed_edge).norm() = || e ||_M = sqrt(e^t M e)
        Eigen::Matrix3d f = get_interpolated_transformation(m0, m1);
        Vector3d transformed_curr_edge = f*normalized_anc_edge;

        FT sp = curr_edge.dot(normalized_anc_edge);
        FT l = transformed_curr_edge.norm(); // length of the normalized anc edge in the metric

        dist_to_ancestor += sp * l;
      }
      dist_to_ancestor = (std::max)(dist_to_ancestor, 0.);

      // add ancestor edge length to the distance at that ancestor
      FT dist_at_anc = curr_ancestor->distance_to_closest_seed();
      FT new_d = dist_at_anc + dist_to_ancestor;

      if(new_d < d)
        d = new_d;

      if(!curr_ancestor->ancestor()) // can't go any farther up in the ancestor tree
        break;

      curr_ancestor = curr_ancestor->ancestor();
    }

    if(d < this->distance_to_closest_seed())
    {
      this->m_ancestor = anc;
      this->distance_to_closest_seed() = d;
      this->closest_seed_id() = anc->closest_seed_id();
      return true;
    }

    return false;
  }

  PQ_state update_neighbors_distances(std::vector<Campen_canvas_point*>& trial_pq) const
  {
    // consider all the neighbors of a KNOWN point and compute their distance to 'this'
#if (verbosity > 15)
    std::cout << "update neighbors of " << this->index() << std::endl;
#endif
    CGAL_assertion(this->state() == KNOWN);

    PQ_state pqs_ret = NOTHING_TO_DO;
    typename Neighbors::const_iterator it = neighbors.begin(),
                                       iend = neighbors.end();
    for(; it!=iend; ++it)
    {
      Campen_canvas_point* cp = *it;
      if(!cp)
        continue;
      else if(cp->state() == KNOWN)
        continue; // dual_shenanigans(cp);
      else if(cp->state() == TRIAL)
      {
        // note that we don't insert in trial_pq since it's already in
        if(cp->compute_closest_seed(this))
          pqs_ret = REBUILD_TRIAL;
      }
      else // cp->state == FAR
      {
        CGAL_assertion(cp->state() == FAR);

        // note that cp->distance_to_closest_seed is not necessarily FT_inf here :
        // if we're refining, we've assigned FAR to all points after inserting a new
        // seed, therefore we must verify that compute_closest_seed is an update
        // before inserting it in the trial_queue
        if(cp->compute_closest_seed(this))
        {
          cp->state() = TRIAL;
          trial_pq.push_back(cp);
          std::push_heap(trial_pq.begin(), trial_pq.end(),
                         Canvas_point_comparer<Campen_canvas_point>());
        }
      }
    }
    return pqs_ret;
  }

  template<typename MF>
  Campen_canvas_point(const Point_3& _point,
                      const std::size_t _index,
                      const MF _mf)
    :
      Base(_point, _index, _mf),
      neighbors()
  {
    for(std::size_t i=0; i<neighbors.size(); ++i)
      neighbors[i] = NULL;
  }
};

} // namespace Anisotropic_mesh_3
} // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_CANVAS_GRID_POINT_H

#ifndef CGAL_ANISOTROPIC_MESH_3_CAMPEN_STARSET_POINT_H
#define CGAL_ANISOTROPIC_MESH_3_CAMPEN_STARSET_POINT_H

#include <CGAL/Canvas/canvas_config.h>
#include <CGAL/Canvas/canvas_enum.h>
#include <CGAL/Canvas/canvas_point.h>
#include <CGAL/Canvas/canvas_helper.h>

#include <CGAL/Anisotropic_mesher_3.h>
#include <CGAL/Starset.h>
#include <CGAL/IO/Star_set_IO.h>

#include <CGAL/Metric.h>

#include <CGAL/assertions.h>

#include <boost/array.hpp>
#include <Eigen/Dense>

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

bool is_a_corner_vertex(std::size_t id)
{
  return id < 8;
}

template<typename Vh>
bool is_a_corner_vertex(Vh vh)
{
  return vh->info() < 8;
}

template<typename Star>
bool has_a_corner_vertex(Star* star)
{
  typename Star::Vertex_handle_handle vit = star->finite_adjacent_vertices_begin();
  typename Star::Vertex_handle_handle vend = star->finite_adjacent_vertices_end();
  for(; vit != vend; vit++)
  {
    typename Star::Vertex_handle vh = *vit;
    if(is_a_corner_vertex(vh))
      return true;
  }
  return false;
}

template<typename K, typename Domain, typename MF, typename Criteria>
class Campen_starset_canvas;

template<typename K, typename Domain, typename MF, typename Criteria>
class Campen_starset_point :
    public Canvas_point<K, Campen_starset_canvas<K, Domain, MF, Criteria> >
{
private:
  typedef Campen_starset_point<K, Domain, MF, Criteria>               Self;

public:
  typedef std::vector<Self>                            Campen_starset_point_vector;
  typedef Self*                                        Campen_starset_point_handle;

  typedef int                          Vertex_Info; // index of the canvas point
  typedef int                          Cell_Info; // index of the subdomain

  typedef Campen_starset_canvas<K, Domain, MF, Criteria>              Canvas;
  typedef Canvas_point<K, Canvas>                                     Base;

  typedef typename K::FT                                              FT;
  typedef typename K::Point_3                                         Point_3;
  typedef typename Base::Metric                                       Metric;
  typedef typename Base::Vector3d                                     Vector3d;

  typedef Starset_with_info<K, Domain, MF, Criteria>                  Star_set;
  typedef typename Star_set::Star_handle                              Star_handle;
  typedef typename Star_set::Vertex_handle                            Vertex_handle;
  typedef typename Star_set::Vertex_handle_handle                     Vertex_handle_handle;
  typedef typename Star_set::Cell_handle                              Cell_handle;
  typedef typename Star_set::Cell_handle_handle                       Cell_handle_handle;

  Star_set* m_ss;

  // this function is the heart of the painter
  bool compute_closest_seed(const Self& anc,
                            const bool verb = false)
  {
    // returns true if we improved the distance

#if (VERBOSITY > 20)
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "compute closest seed for : " << this->index();
    std::cout << " (" << this->point() << ") ";
    std::cout << "curr. dist: " << this->distance_to_closest_seed();
    std::cout << " anc: " << anc.index() << " & ancdist: ";
    std::cout << anc.distance_to_closest_seed() << std::endl;
#endif
    CGAL_assertion(anc.state() == KNOWN);

    const int k = 8; // depth of the ancestor edge
    FT d = FT_inf;

    // the path is 'this' to ancestor1, ancestor1 to ancestor2, etc.
    // stored as 'this', ancestor1, ancestor2, etc.
    boost::array<std::size_t, k+1> ancestor_path;
    for(int i=1; i<k+1; ++i)
      ancestor_path[i] = -1;
    ancestor_path[0] = this->index();
    ancestor_path[1] = anc.index();

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

    std::size_t curr_anc = anc.index();
    for(int i=1; i<=k; ++i)
    {
      const Base& curr_ancestor = this->canvas()->get_point(curr_anc);
#if (VERBOSITY > 25)
      std::cout << "current ancestor: " << curr_ancestor.index() << std::endl;
#endif

      Vector3d ancestor_edge;
      ancestor_edge(0) = this->point().x() - curr_ancestor.point().x();
      ancestor_edge(1) = this->point().y() - curr_ancestor.point().y();
      ancestor_edge(2) = this->point().z() - curr_ancestor.point().z();
      FT ancestor_edge_length = ancestor_edge.norm();
      CGAL_assertion(ancestor_edge_length > 1e-16);
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

//        const Metric& m0 = e0->metric();
//        const Metric& m1 = e1->metric();
//        CGAL_assertion(f == get_interpolated_transformation(m0, m1));

        Vector3d transformed_curr_edge = f*normalized_anc_edge;

        FT sp = curr_edge.dot(normalized_anc_edge);
        FT l = transformed_curr_edge.norm(); // length of the normalized anc edge in the metric

        dist_to_ancestor += sp * l;
#if (VERBOSITY > 30)
        std::cout << "e0: " << e0.point() << " e1: " << e1.point() << std::endl;
        std::cout << "normalized_anc_edge: " << normalized_anc_edge.transpose() << std::endl;
//        std::cout << "metrics:" << std::endl << m0.get_mat() << std::endl << m1.get_mat() << std::endl;
        std::cout << "transf edge: " << transformed_curr_edge.transpose() << std::endl;
        std::cout << "dist_to_anc building: " << dist_to_ancestor << " " << sp << " " << l << std::endl;
#endif
      }
      dist_to_ancestor = (std::max)(dist_to_ancestor, 0.);

      // add ancestor edge length to the distance at that ancestor
      FT dist_at_anc = curr_ancestor.distance_to_closest_seed();
      FT new_d = dist_at_anc + dist_to_ancestor;

#if (VERBOSITY > 25)
      std::cout << "potential update: " << new_d << " " << dist_at_anc
                << " " << dist_to_ancestor << std::endl;
#endif

      if(new_d < d)
        d = new_d;

      // checks if we can't go any farther up in the ancestor tree
      if(curr_ancestor.ancestor() == static_cast<std::size_t>(-1))
        break;

      curr_anc = curr_ancestor.ancestor();
    }

#if (VERBOSITY > 20)
    std::cout << "distance with that anc: " << d << std::endl;
#endif

    if(d < this->distance_to_closest_seed())
    {
      if(verb)
      {
        std::cout << "improving distance at " << Base::index() << " "
                  << " from: " << this->distance_to_closest_seed()
                  << " to " << d << std::endl;
        if(this->ancestor())
          std::cout << "prev anc: " << this->ancestor() << std::endl;
        else
          std::cout << "no prev anc" << std::endl;
        std::cout << "new anc: " << anc.index() << std::endl;
      }

      // remove 'index' from the previous ancestor's children (if needed)
      if(this->ancestor() != static_cast<std::size_t>(-1))
        this->canvas()->get_point(this->ancestor()).remove_from_children(this->index());

      this->ancestor() = anc.index();
      this->distance_to_closest_seed() = d;
      this->closest_seed_id() = anc.closest_seed_id();

      // add 'index' to the new ancestor's children
      anc.m_children.insert(this->index());

//      if(!this->children().empty())
//      {
//        std::cout << "dealing with the " << this->children().size()
//                  << " descendant(s) of " << this->index() << std::endl;
//        std::cout << "position : " << this->point() << std::endl;
//      }

      while(!this->children().empty())
      {
        Self& cp = this->canvas()->get_point(*(this->children().begin()));
        cp.reset_descendants();
        CGAL_postcondition(cp.ancestor_path_length()); // checks for circular ancestry
      }

      return true;
    }
    return false;
  }

  PQ_state update_neighbors_distances(std::vector<Self*>& trial_pq)
  {
    // consider all the neighbors of a KNOWN point and compute their distance to 'this'
#if (VERBOSITY > 15)
    std::cout << "update neighbors of " << Base::index() << std::endl;
#endif
    CGAL_assertion(this->state() == KNOWN);

    PQ_state pqs_ret = NOTHING_TO_DO;

    Star_handle star = m_ss->get_star(this->index());
    Vertex_handle_handle it = star->finite_adjacent_vertices_begin();
    Vertex_handle_handle end = star->finite_adjacent_vertices_end();
    for(; it!=end; ++it)
    {
      Vertex_handle vh = *it;
      if(is_a_corner_vertex(vh))
        continue;

      Self& cp = this->canvas()->get_point(vh->info());
      if(cp.state() == KNOWN)
        continue;
      else if(cp.state() == TRIAL)
      {
        // note that we don't insert in trial_pq since it's already in
        if(cp.compute_closest_seed(*this))
          pqs_ret = REBUILD_TRIAL;
      }
      else // cp.state == FAR
      {
        CGAL_assertion(cp.state() == FAR);

        // note that cp.distance_to_closest_seed is not necessarily FT_inf here :
        // if we're refining, we've assigned FAR to all points after inserting a new
        // seed, therefore we must verify that compute_closest_seed is an update
        // before inserting it in the trial_queue
        if(cp.compute_closest_seed(*this))
        {
          CGAL_assertion(cp.distance_to_closest_seed() != FT_inf);
          cp.state() = TRIAL;
          trial_pq.push_back(&cp);
          std::push_heap(trial_pq.begin(), trial_pq.end(),
                         Canvas_point_comparer<Self>());
        }
      }
      CGAL_assertion(cp.distance_to_closest_seed() != FT_inf);
    }
    return pqs_ret;
  }

  Campen_starset_point() : Base(), m_ss(NULL) { }

  Campen_starset_point(const Point_3& p, const std::size_t index,
                       bool border_info,
                       Star_set* ss_, Canvas* canvas)
    :
      Base(p, index, canvas, border_info),
      m_ss(ss_)
  { }
};

} // namespace Anisotropic_mesh_3
} // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_CAMPEN_STARSET_POINT_H

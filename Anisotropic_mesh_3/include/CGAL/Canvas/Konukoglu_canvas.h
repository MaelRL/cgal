#ifndef CGAL_ANISOTROPIC_MESH_3_KONUKOGLU_CANVAS_H
#define CGAL_ANISOTROPIC_MESH_3_KONUKOGLU_CANVAS_H

#include <CGAL/Canvas/canvas_point.h>
#include <CGAL/Canvas/grid_canvas.h>
#include <CGAL/Canvas/Konukoglu_point.h>

#include <boost/array.hpp>
#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K, typename KExact, typename Metric_field>
class Konukoglu_canvas :
    public Grid_canvas<K, Konukoglu_canvas_point<K, KExact, Metric_field>, Metric_field>
{
public:
  typedef Konukoglu_canvas_point<K, KExact, Metric_field>    Canvas_point;
  typedef Grid_canvas<K, Canvas_point, Metric_field>         Base;

  typedef typename Base::Kernel                              Kernel;
  typedef typename Kernel::FT                                FT;
  typedef typename Kernel::Point_3                           Point_3;

  typedef typename Base::Canvas_point_handle_vector          Canvas_point_handle_vector;

  // a priority queue specific to this algorithm for 'KNOWN' points whose value changes
  Canvas_point_handle_vector changed_points;

  // tolerance on how large of an update the new value needs to be
  FT recursive_tolerance;

  // the exact same function exists in the base class of the canvas, but since
  // the points use the canvas in the constructor and need to know the fully
  // derived canvas, it's put here again so *this is a Konukoglu_canvas...
  // Could do CRTP and used derived() etc... but it's a hassle.
  void initialize()
  {
#if (verbosity > 5)
    std::cout << "canvas initialization" << std::endl;
#endif

    // create the canvas points
    for(unsigned int k=0; k<this->n; ++k)
    {
      for(unsigned int j=0; j<this->n; ++j)
      {
        for(unsigned int i=0; i<this->n; ++i)
        {
          Point_3 p(this->offset_x + i*this->step,
                    this->offset_y + j*this->step,
                    this->offset_z + k*this->step); // fill from bot left to top right
          Canvas_point cp(p, i + j*this->n + k*this->sq_n, this);
          this->canvas_points.push_back(cp);
        }
      }
    }

    // assign the neighbors
    for(unsigned int k=0; k<this->n; ++k)
    {
      for(unsigned int j=0; j<this->n; ++j)
      {
        for(unsigned int i=0; i<this->n; ++i)
        {
          std::size_t curr_id = i + j*this->n + k*this->sq_n;
          if(k != this->n-1) // there is a neighbor above
            this->canvas_points[curr_id].neighbors[0] = &(this->canvas_points[i + j*this->n + (k+1)*this->sq_n]);
          if(i != 0) // there is a neighbor left
            this->canvas_points[curr_id].neighbors[1] = &(this->canvas_points[i-1 + j*this->n + k*this->sq_n]);
          if(j != this->n-1) // there is a neighbor back
            this->canvas_points[curr_id].neighbors[2] = &(this->canvas_points[i + (j+1)*this->n + k*this->sq_n]);
          if(i != this->n-1) // there is a neighbor right
            this->canvas_points[curr_id].neighbors[3] = &(this->canvas_points[i+1 + j*this->n + k*this->sq_n]);
          if(j != 0) // there is a neighbor front
            this->canvas_points[curr_id].neighbors[4] = &(this->canvas_points[i + (j-1)*this->n + k*this->sq_n]);
          if(k != 0) // there is a neighbor below
            this->canvas_points[curr_id].neighbors[5] = &(this->canvas_points[i + j*this->n + (k-1)*this->sq_n]);
        }
      }
    }
#if (verbosity > 5)
    std::cout << "neighbors assigned" << std::endl;
#endif

    Base::compute_canvas_bbox();
    this->seeds.initialize_seeds();
    Base::locate_seeds_on_canvas();

#if (verbosity > 5)
    std::cout << "canvas initialized" << std::endl;
#endif
  }

  void determine_ancestors()
  {
    // (try to) assign the correct ancestor to each point (coloring the canvas)

    std::clock_t start = std::clock();

#if (verbosity > 0)
    std::cout << "Determine ancestors" << std::endl;
#endif
    CGAL_assertion(changed_points.empty() && this->trial_points.empty());

    // TODO do something smarter when spreading a new color in an already colored
    // grid... One can probably just start from the new seed WITHOUT RESETING EVERYTHING
    // and then spread and only consider stuff if one of the possible ancestors
    // has the new color something like that...

    Canvas_point_handle_vector known_points;
    for(std::size_t i=0; i<this->canvas_points.size(); ++i)
    {
      // reset the colors (but not the values!)

      // fixme, just leave the colors uninitialized in the first pass and
      // there's no need to reset them here...

      Canvas_point* cp = &(this->canvas_points[i]);
      cp->closest_seed_id() = -1;
      cp->m_ancestor = NULL;

      CGAL_assertion(cp->state() == KNOWN);

      known_points.push_back(cp);
    }
    std::make_heap(known_points.begin(), known_points.end(),
                   Canvas_point_comparer<Canvas_point>());

    // only re initialize the closest_seed_id for each point-seed
    for(std::size_t i=0; i<this->seeds.size(); ++i)
    {
      // if you use an 8 pts initilization for a seed, you need to modify here too...
      const Point_3& p = this->seeds[i];
      int index_x = std::floor((p.x() - this->offset_x) / this->step);
      int index_y = std::floor((p.y() - this->offset_y) / this->step);
      int index_z = std::floor((p.z() - this->offset_z) / this->step);
      Canvas_point* cp = &(this->canvas_points[index_z * this->sq_n + index_y * this->n + index_x]);
      cp->closest_seed_id() = i;
    }

    bool is_kp_empty = known_points.empty();
    while(!is_kp_empty)
    {
      Canvas_point* cp = known_points.front();
      std::pop_heap(known_points.begin(), known_points.end(),
                    Canvas_point_comparer<Canvas_point>());
      known_points.pop_back();

      cp->determine_ancestor();
      is_kp_empty = known_points.empty();
    }

    std::cerr << "End of determine ancestors. time: ";
    std::cerr << ( std::clock() - start ) / (double) CLOCKS_PER_SEC << std::endl;
  }

  void paint()
  {
    std::clock_t start = std::clock();
#if (verbosity > 0)
    std::cout << "main loop" << std::endl;
#endif

    bool is_cp_empty = changed_points.empty();
    bool is_t_empty = this->trial_points.empty();

    Canvas_point* cp;
    while(!is_cp_empty || !is_t_empty)
    {
#if (verbosity > 5)
      std::cout << "Queue sizes. Trial: " << this->trial_points.size() << " Changed: " << changed_points.size() << std::endl;
#endif

#if (verbosity > 20)
      std::cout << "changed heap: " << std::endl;
      typename std::vector<Canvas_point*>::iterator it = changed_points.begin();
      for (; it != changed_points.end(); ++it)
        std::cout << (*it)->index() << " " << (*it)->distance_to_closest_seed() << std::endl;
      std::cout << std::endl;

      std::cout << "trial heap: " << std::endl;
      it = this->trial_points.begin();
      for (; it != this->trial_points.end(); ++it)
        std::cout << (*it)->index() << " " << (*it)->distance_to_closest_seed() << std::endl;
      std::cout << std::endl;
#endif
      if(!is_cp_empty)
      {
        cp = changed_points.front();
        CGAL_assertion(cp && cp->state() == CHANGED);
        std::pop_heap(changed_points.begin(), changed_points.end(),
                      Canvas_point_comparer<Canvas_point>());
        changed_points.pop_back();
      }
      else // !is_t_empty
      {
        cp = this->trial_points.front();
        CGAL_assertion(cp && cp->state() == TRIAL);
        std::pop_heap(this->trial_points.begin(), this->trial_points.end(),
                      Canvas_point_comparer<Canvas_point>());
        this->trial_points.pop_back();
      }

#if (verbosity > 5)
      std::cout << "picked n° " << cp->index() << " (" << cp->point() << ")";
      std::cout << "at distance : " << cp->distance_to_closest_seed() << " from " << cp->closest_seed_id() << std::endl;
#endif

      cp->state() = KNOWN;
      PQ_state pqs = cp->update_neighbors_distances();

      if(pqs == REBUILD_TRIAL || pqs == REBUILD_BOTH)
        std::make_heap(this->trial_points.begin(), this->trial_points.end(),
                       Canvas_point_comparer<Canvas_point>());
      if(pqs == REBUILD_CHANGED || pqs == REBUILD_BOTH)
        std::make_heap(changed_points.begin(), changed_points.end(),
                       Canvas_point_comparer<Canvas_point>());

      is_cp_empty = changed_points.empty();
      is_t_empty = this->trial_points.empty();
    }

    std::cerr << "End of geo_grid_loop. time: ";
    std::cerr << ( std::clock() - start ) / (double) CLOCKS_PER_SEC << std::endl;

    determine_ancestors();
  }

  void debug()
  {
#ifdef USE_RECURSIVE_UPDATES
    CGAL_assertion(this->trial_points.empty() && changed_points.empty());

    std::cout << "BRUTE FORCE CHECK THAT WE ARE REALLY FINISHED : " << std::endl;

    // it would probably be best to consider all the possible 2D cases here...
    // but more generally, 3D being weaker than 2D is a problem to be fixed

    // detail of the problem :
    // 2D can be solved analytically while 3D cannot. Therefore we might sometimes
    // lose the exact result when we 'unlock' the 3D if the min is achieved on a face

    for(std::size_t i=0; i<this->canvas_points.size(); ++i)
    {
      Canvas_point* cp = &(this->canvas_points[i]);
      std::cout << "point " << i << " min distance is supposedly: ";
      std::cout << cp->distance_to_closest_seed() << std::endl;
      typename Base::Neighbors::const_iterator it = cp->neighbors.begin(),
                                               iend = cp->neighbors.end();
      for(; it!=iend; ++it)
      {
        const Canvas_point* cq = *it;
        if(cq)
          CGAL_assertion(!cp->compute_closest_seed(cq));
      }
    }
#endif
  }

  Konukoglu_canvas(const std::string& canvas_str_,
                   const std::string& seeds_str_,
                   const Point_3& center_,
                   const FT side_,
                   const std::size_t points_per_side,
                   const std::size_t max_seeds_n_,
                   const Metric_field& mf_,
                   const FT tolerance_)
    :
      Base(canvas_str_, seeds_str_,
           center_, side_, points_per_side,
           max_seeds_n_,
           mf_),
      changed_points(),
      recursive_tolerance(tolerance_)
  { }
};

} // namespace Anisotropic_mesh_3
} // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_KONUKOGLU_CANVAS_H

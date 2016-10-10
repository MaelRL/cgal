#ifndef CGAL_ANISOTROPIC_MESH_3_CANVAS_H
#define CGAL_ANISOTROPIC_MESH_3_CANVAS_H

#include <CGAL/Canvas/canvas_config.h>
#include <CGAL/Canvas/canvas_enum.h>
#include <CGAL/Canvas/canvas_helper.h>
#include <CGAL/Canvas/canvas_seeds.h>
#include <CGAL/Canvas/primal_simplex.h>
#include <CGAL/Canvas/canvas_tetrahedron_intersection.h>

#include <CGAL/Metric.h>
#include <CGAL/Bbox_3.h>

#include <CGAL/box_intersection_d.h>

#include <Eigen/Dense>
#include <boost/array.hpp>
#include <boost/unordered_set.hpp>

#include <algorithm>
#include <cstddef>
#include <ctime>
#include <iostream>
#include <ostream>
#include <set>
#include <string>
#include <vector>
#include <utility>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K, typename Cpoint, typename MF>
class Canvas
{
private:
  typedef Canvas<K, Cpoint, MF>                             Self;

public:
  typedef Canvas_seeds<Self>                                Seeds;
  typedef K                                                 Kernel;
  typedef Cpoint                                            Canvas_point;
  typedef MF                                                Metric_field;

  typedef typename K::FT                                    FT;
  typedef typename K::Point_3                               Point_3;
  typedef typename K::Segment_3                             Segment_3;
  typedef typename K::Triangle_3                            Triangle_3;
  typedef typename K::Tetrahedron_3                         Tetrahedron_3;

  typedef Metric_base<K>                                    Metric;

  typedef std::set<std::size_t>                             Simplex;
  typedef boost::array<std::size_t, 2>                      BEdge;
  typedef boost::array<std::size_t, 3>                      BTriangle;
  typedef boost::array<std::size_t, 4>                      BTetrahedron;

  typedef Primal_simplex<Self, 2>                           Primal_edge;
  typedef Primal_simplex<Self, 3>                           Primal_triangle;
  typedef Primal_simplex<Self, 4>                           Primal_tetrahedron;

  typedef boost::unordered_set<Primal_edge,
                               Primal_simplex_hash<Canvas, 2>,
                               Primal_simplex_comparer<Canvas, 2> >
                                                            Primal_edges_container;
  typedef typename Primal_edges_container::iterator         PEC_iterator;
  typedef boost::unordered_set<Primal_triangle,
                               Primal_simplex_hash<Canvas, 3>,
                               Primal_simplex_comparer<Canvas, 3> >
                                                            Primal_triangles_container;
  typedef typename Primal_triangles_container::iterator     PTrC_iterator;
  typedef boost::unordered_set<Primal_tetrahedron,
                               Primal_simplex_hash<Canvas, 4>,
                               Primal_simplex_comparer<Canvas, 4> >
                                                            Primal_tetrahedra_container;
  typedef typename Primal_tetrahedra_container::iterator    PTC_iterator;

  typedef std::vector<Canvas_point*>                        Canvas_point_handle_vector;

  typedef Eigen::Matrix<FT, 3, 1>                           Vector3d;

  typedef CGAL::Bbox_3                                      Bbox;
  typedef CGAL::Box_intersection_d::Box_with_handle_d<FT, 3, PTC_iterator>
                                                            Box;

  // Canvas name
  const std::string canvas_str;

  // Canvas geometry
  std::vector<Canvas_point> canvas_points;
  Canvas_point_handle_vector trial_points;

  // Seeds
  Seeds seeds;

  // Metric field
  const Metric_field* mf;

  // Primal triangulation
  Primal_edges_container primal_edges;
  Primal_triangles_container primal_triangles;
  Primal_tetrahedra_container primal_tetrahedra;
  bool are_tetrahedra_intersection_computed;

  // virtual points are used when we compute precise Voronoi vertices
  // that live on the canvas but are not actual canvas points
  std::vector<Canvas_point> virtual_points;

  // Bbox of the canvas
  CGAL::Bbox_3 canvas_bbox;

  // debug & info
  std::size_t known_count, trial_count, far_count;

  // ---------------------------------------------------------------------------
  // virtual functions
  virtual void initialize() = 0;
  virtual void output_canvas(const std::string str_base) const = 0;
  virtual void compute_primal() = 0;
  virtual void compute_local_primal_elements(const Canvas_point* cp) = 0;
  // ---------------------------------------------------------------------------

  Canvas_point& get_point(std::size_t i)
  {
    CGAL_precondition(i < canvas_points.size());
    return canvas_points[i];
  }

  const Canvas_point& get_point(std::size_t i) const
  {
    CGAL_precondition(i < canvas_points.size());
    return canvas_points[i];
  }

  bool is_point_outside_canvas(const FT x, const FT y, const FT z) const
  {
    return (x < canvas_bbox.xmin() || x > canvas_bbox.xmax() ||
            y < canvas_bbox.ymin() || y > canvas_bbox.ymax() ||
            z < canvas_bbox.zmin() || z > canvas_bbox.zmax());
  }

  bool is_point_outside_canvas(const Point_3& p) const
  {
    return is_point_outside_canvas(p.x(), p.y(), p.z());
  }

  void clear_primal()
  {
    primal_edges.clear();
    primal_triangles.clear();
    primal_tetrahedra.clear();
    are_tetrahedra_intersection_computed = false;
    virtual_points.clear();
  }

  void initialize_canvas_point(Canvas_point& cp,
                               const FT distance_from_seed,
                               const std::size_t seed_id)
  {
    if(cp.closest_seed_id() != static_cast<std::size_t>(-1))
    {
      std::cout << "WARNING: a new seed is overwriting the closest seed id";
      std::cout << " of a canvas point!" << std::endl;
      std::cout << "seed: " << seed_id << " wants to initialize point: " << cp.index() << "(" << cp.point() << ")";
      std::cout << " but seed " << cp.closest_seed_id() << " has done it already" << std::endl;
    }

    // We can't accept two seeds for one canvas point
    if(cp.state() == TRIAL)
      CGAL_assertion(false && "the canvas is not dense enough for the input seeds...");

    cp.initialize_from_point(distance_from_seed, seed_id);

    trial_points.push_back(&cp);
    std::push_heap(trial_points.begin(), trial_points.end(),
                   Canvas_point_comparer<Canvas_point>());
  }

  virtual void locate_and_initialize(const Point_3& p,
                                     const std::size_t seed_id)
  {
    // something pretty todo...
    // for now: brutally find the canvas point closest to the seed

    int cp_id = -1; // index of the closest seed

    FT min_d = FT_inf;
    for(std::size_t i=0, cps=canvas_points.size(); i<cps; ++i)
    {
      Vector3d v;
      v(0) = p.x() - canvas_points[i].point().x();
      v(1) = p.y() - canvas_points[i].point().y();
      v(2) = p.z() - canvas_points[i].point().z();
      const Eigen::Matrix3d& m = canvas_points[i].metric().get_mat();

      FT d = v.transpose() * m * v; // note that this is the squared_distance
      if(d < min_d)
      {
        min_d = d;
        cp_id = i;
      }
    }

    Canvas_point& cp = canvas_points[cp_id];

#if (VERBOSITY > 15)
    std::cout << "looking for p: " << p << std::endl;
    std::cout << "found cp: " << cp.index() << " [" << cp.point() << "] ";
    std::cout << "at distance: " << std::sqrt(min_d) << std::endl;
#endif

    initialize_canvas_point(cp, std::sqrt(min_d), seed_id);
  }

  void compute_canvas_bbox()
  {
    canvas_bbox = CGAL::Bbox_3();
    for(std::size_t i=0, cps=canvas_points.size(); i<cps; ++i)
      canvas_bbox += canvas_points[i].point().bbox();

    // std::cout << "WARNING: hardcoding bbox" << std::endl;
    // canvas_bbox = CGAL::Bbox_3(0.,0.,0.,3.,3.,3.);

    std::cout << "final bbox: "
              << "x " << canvas_bbox.xmin() << " " << canvas_bbox.xmax() << '\n'
              << "y " << canvas_bbox.ymin() << " " << canvas_bbox.ymax() << '\n'
              << "z " << canvas_bbox.zmin() << " " << canvas_bbox.zmax() << std::endl;
  }

  void locate_seeds_on_canvas()
  {
    // find the canvas vertex closest to the seed
    for(std::size_t i=0; i<seeds.size(); ++i)
    {
      const Point_3& p = seeds[i];
      locate_and_initialize(p, i);
    }
  }

  void print_states() const
  {
    std::cout << "known: " << known_count;
    std::cout << " trial: " << trial_count;
    std::cout << " far: " << far_count << std::endl;
  }

  virtual void paint(const bool refining = false)
  {
#if (VERBOSITY > 5)
    std::cout << "Paiting..." << std::endl;
#endif

    std::clock_t start = std::clock();
    clear_primal();
    Canvas_point* cp;

    bool is_t_empty = trial_points.empty();

    if(is_t_empty)
    {
      std::cout << "Trying to paint without anything in the PQ..." << std::endl;
      exit(0);
    }
    else
      std::cout << trial_points.size() << " initial points in the queue" << std::endl;

    while(!is_t_empty)
    {
#if (VERBOSITY > 10)
      if(known_count % ((std::max)(static_cast<std::size_t>(1),
                                   canvas_points.size()/100)) == 0)
        print_states();
#endif

#if (VERBOSITY > 15)
      std::cout << "Trial queue size : " << trial_points.size() << std::endl;
#endif

#if (VERBOSITY > 55)
      std::cout << "trial heap: " << std::endl;
      for(typename std::vector<Canvas_point*>::iterator it = trial_points.begin();
           it != trial_points.end(); ++it)
        std::cout << (*it)->index() << " " << (*it)->distance_to_closest_seed() << std::endl;
      std::cout << std::endl;
#endif

      cp = trial_points.front();
      CGAL_assertion(cp);
      CGAL_assertion(cp->state() == TRIAL);
      std::pop_heap(trial_points.begin(), trial_points.end(),
                    Canvas_point_comparer<Canvas_point>());
      trial_points.pop_back();

#if (VERBOSITY > 15)
      std::cout << "picked n° " << cp->index() << " (" << cp->point() << ") ";
      std::cout << "at distance : " << cp->distance_to_closest_seed()
                << " from " << cp->closest_seed_id() << std::endl;
#endif

      cp->state() = KNOWN;
      ++known_count; // todo --> use change_state()

      if(false/*!refining*/) // tmp fix primal shenanigans properly before using them again
      {
        // if we're refining, we're only painting a small part of the canvas
        // (the new voronoi cell) and maintening the primal triangulation is
        // not currently supported (todo)
        compute_local_primal_elements(cp);
      }

      PQ_state pqs = cp->update_neighbors_distances(trial_points);

      if(pqs == REBUILD_TRIAL)
        std::make_heap(trial_points.begin(), trial_points.end(),
                       Canvas_point_comparer<Canvas_point>());

      is_t_empty = trial_points.empty();
    }

    CGAL_expensive_assertion_code(debug());

#if (VERBOSITY > 15)
    std::cout << "final states after painting: " << std::endl;
    for(std::size_t i=0, cps=canvas_points.size(); i<cps; ++i)
    {
      const Canvas_point& cp = canvas_points[i];
      std::cout << cp.index() << " (" << cp.point() << ")";
      std::cout << " at distance: " << cp.distance_to_closest_seed();
      std::cout << " from " << cp.closest_seed_id() << std::endl;
      cp.print_ancestor_tree();
    }
#endif

#if (VERBOSITY > 5)
    std::cout << "End of paint. time: ";
    std::cout << ( std::clock() - start ) / (double) CLOCKS_PER_SEC << std::endl;
#endif
  }

  virtual void debug()
  {
    CGAL_assertion(trial_points.empty() );

    for(std::size_t i=0, cps=canvas_points.size(); i<cps; ++i)
    {
      const Canvas_point& cp = canvas_points[i];

      if(cp.ancestor() != static_cast<std::size_t>(-1) &&
         canvas_points[cp.ancestor()].children().find(cp.index()) ==
         canvas_points[cp.ancestor()].children().end())
      {
        std::cout << "failure in ancestor/children relationship at " << cp.index() << std::endl;
        CGAL_assertion(false);
      }

      CGAL_postcondition(cp.distance_to_closest_seed() != FT_inf);
      CGAL_postcondition(cp.closest_seed_id() < seeds.size());
    }
  }

  void reset_counters()
  {
    known_count = 0;
    trial_count = 0;
    far_count = canvas_points.size();
  }

  void refresh_canvas_point_states()
  {
    CGAL_assertion(trial_points.empty());
    for(std::size_t i=0, cps=canvas_points.size(); i<cps; ++i)
      canvas_points[i].state() = FAR;

    reset_counters();
  }

  void set_points_states_to_known()
  {
    for(std::size_t i=0, cps=canvas_points.size(); i<cps; ++i)
      canvas_points[i].state() = KNOWN;

    known_count = canvas_points.size();
    trial_count = 0;
    far_count = 0;
  }

  void reset()
  {
    // Remove everything related to "paint" but not the geometric information
    // (therefore the point position, metric, neighbors, etc. aren't reset)

    for(std::size_t i=0, cps=canvas_points.size(); i<cps; ++i)
    {
      Canvas_point& cp = canvas_points[i];
      cp.state() = FAR;
      cp.distance_to_closest_seed() = FT_inf;
      cp.closest_seed_id() = -1;
      cp.ancestor() = -1;
      cp.children().clear();
    }

    reset_counters();
    clear_primal();
  }

  // ---------------------------------------------------------------------------
  // Functions to insert in primal data structure

  // These functions are a bit heavy, but it's a tedious process... there's
  // probably something more elegant out there.
  // The purpose of these functions is : given a set of canvas points, build the
  // primal simplices associated with all the possible color combinations

  // the main issue that makes this so heavy is that the simplicial complex
  // is not necessarily pure so we've got to be careful and add at all levels

  // there might be a way to make it less lengthy with partial specilization
  // but it's not worth the hassle for now...

  template<std::size_t size>
  const Canvas_point*
  compute_primal_simplex_dual_point(const boost::array<const Canvas_point*, size>& witnesses)
  {
    // this base version simply returns the farthest witness, but a more complicated
    // version might want to compute a more precise value for the dual of a simplex
    // which is why the compute_primal_*_dual_point functions are virtual
    // (this one can't be since it's templated...)

    const Canvas_point* ret_point = NULL;
    FT max_dist = -FT_inf;
    for(std::size_t i=0; i<size; ++i)
    {
      if(witnesses[i]->distance_to_closest_seed() > max_dist)
      {
        max_dist = witnesses[i]->distance_to_closest_seed();
        ret_point = witnesses[i];
      }
    }
    return ret_point;
  }

  template<typename CPoint_Barray>
  bool check_primal_simplex(const CPoint_Barray& witnesses,
                            boost::array<std::size_t, CPoint_Barray::static_size>& bsimplex)
  {
    // the witnesses must have all different colors here

    for(std::size_t i=0; i<witnesses.size(); ++i)
      bsimplex[i] = witnesses[i]->closest_seed_id();

    std::sort(bsimplex.begin(), bsimplex.end());

    bool duplicate_found = false;
    std::size_t pos = 0;
    while(pos < bsimplex.size()-1)
    {
      if(bsimplex[pos] == bsimplex[pos+1])
      {
        duplicate_found = true;
        break;
      }
      ++pos;
    }

    return !duplicate_found;
  }

  virtual const Canvas_point* compute_primal_edge_dual_point(const boost::array<const Canvas_point*, 2>& witnesses)
  {
    return compute_primal_simplex_dual_point<2>(witnesses);
  }

  void insert_in_primal_edges(const Primal_edge& p_edge)
  {
    // There might already exist an entry for this simplex,
    // but we want to keep the primal edge whose dual point is farthest
    typedef typename Primal_edges_container::iterator           PE_container_it;

    std::pair<PE_container_it, bool> is_insert_successful;
    is_insert_successful = primal_edges.insert(p_edge);

#if (VERBOSITY > 20)
    std::cout << "built primal edge : " << p_edge[0] << " " << p_edge[1] << std::endl;
#endif

    if(!is_insert_successful.second)
    {
      // already exist, check the distances
      FT new_distance = p_edge.dual_point()->distance_to_closest_seed();
      FT old_distance = is_insert_successful.first->dual_point()->distance_to_closest_seed();
#if (VERBOSITY > 20)
      std::cout << "already in; old/new: " << old_distance << " " << new_distance << std::endl;
#endif
      if(new_distance > old_distance)
      {
        // the newer element is better, we must replace the old entry in the primal edge set
#if (VERBOSITY > 20)
        std::cout << "found a farther dual point at : " << p_edge.dual_point()->index();
        std::cout << " (" << p_edge.dual_point()->point() << ")" << std::endl;
#endif
        // fixme what hint should be taken for the unordered set insertion ?
        PE_container_it hint = primal_edges.erase(is_insert_successful.first);
        primal_edges.insert(hint, p_edge);
      }
    }
  }

  void create_and_insert_primal_edge(const boost::array<const Canvas_point*, 2>& witnesses)
  {
    BEdge edge;
    if(!check_primal_simplex(witnesses, edge))
      return;

    const Canvas_point* dual_point = compute_primal_edge_dual_point(witnesses);
    Primal_edge p_edge(edge, dual_point);
    insert_in_primal_edges(p_edge);
  }

  void add_to_primal_edges(typename std::vector<boost::array<const Canvas_point*, 2> >::const_iterator it,
                           typename std::vector<boost::array<const Canvas_point*, 2> >::const_iterator end)
  {
    for(; it!=end; ++it)
      create_and_insert_primal_edge(*it);
  }

  virtual const Canvas_point* compute_primal_triangle_dual_point(const boost::array<const Canvas_point*, 3>& witnesses)
  {
    return compute_primal_simplex_dual_point<3>(witnesses);
  }

  void insert_in_primal_triangles(const Primal_triangle& p_triangle)
  {
    // There might already exist an entry for this simplex,
    // but we want to keep the primal edge whose dual point is farthest
    typedef typename Primal_triangles_container::iterator       PT_container_it;

    std::pair<PT_container_it, bool> is_insert_successful;
    is_insert_successful = primal_triangles.insert(p_triangle);

#if (VERBOSITY > 20)
    std::cout << "built primal triangle : ";
    std::cout << p_triangle[0] << " " << p_triangle[1] << " " << p_triangle[2] << std::endl;
#endif

    if(!is_insert_successful.second)
    {
      // already exist, check the distances
      FT new_distance = p_triangle.dual_point()->distance_to_closest_seed();
      FT old_distance = is_insert_successful.first->dual_point()->distance_to_closest_seed();
#if (VERBOSITY > 20)
      std::cout << "already in; old/new: " << old_distance << " " << new_distance << std::endl;
#endif
      if(new_distance > old_distance)
      {
        // the newer element is better, we must replace the old entry in the primal triangles
#if (VERBOSITY > 20)
        std::cout << "found a farther dual point at : " << p_triangle.dual_point()->index();
        std::cout << " (" << p_triangle.dual_point()->point() << ")" << std::endl;
#endif
        // fixme what hint should be taken for the unordered set insertion ?
        PT_container_it hint = primal_triangles.erase(is_insert_successful.first);
        primal_triangles.insert(hint, p_triangle);
      }
    }
  }

  void create_and_insert_primal_triangle(const boost::array<const Canvas_point*, 3>& witnesses)
  {
#ifndef COMPUTE_PRIMAL_ALL_DIMENSIONS
    // if we're not interested in dual simplicies in all dimensions we only
    // want to keep canvas triangles with 3 different colors AND on the border
    if(!witnesses[0]->border_info() || !witnesses[1]->border_info() ||
       !witnesses[2]->border_info())
      return;
#endif

    BTriangle triangle;
    if(!check_primal_simplex(witnesses, triangle))
      return;

    const Canvas_point* dual_point = compute_primal_triangle_dual_point(witnesses);
    Primal_triangle p_triangle(triangle, dual_point);
    insert_in_primal_triangles(p_triangle);
  }

  void add_to_primal_triangles(typename std::vector<boost::array<const Canvas_point*, 3> >::const_iterator it,
                               typename std::vector<boost::array<const Canvas_point*, 3> >::const_iterator end)
  {
    for(; it!=end; ++it)
      create_and_insert_primal_triangle(*it);
  }

  virtual const Canvas_point* compute_primal_tetrahedron_dual_point(const boost::array<const Canvas_point*, 4>& witnesses)
  {
    return compute_primal_simplex_dual_point<4>(witnesses);
  }

  void insert_in_primal_tetrahedra(const Primal_tetrahedron& p_tetrahedron)
  {
    // There might already exist an entry for this simplex,
    // but we want to keep the primal edge whose dual point is farthest
    typedef typename Primal_tetrahedra_container::iterator      PT_container_it;

    std::pair<PT_container_it, bool> is_insert_successful;
    is_insert_successful = primal_tetrahedra.insert(p_tetrahedron);

#if (VERBOSITY > 20)
    std::cout << "built primal tetrahedron : ";
    std::cout << p_tetrahedron[0] << " " << p_tetrahedron[1] << " "
              << p_tetrahedron[2] << " " << p_tetrahedron[3] << std::endl;
#endif

    if(!is_insert_successful.second)
    {
      // already exist, check the distances
      FT new_distance = p_tetrahedron.dual_point()->distance_to_closest_seed();
      FT old_distance = is_insert_successful.first->dual_point()->distance_to_closest_seed();
#if (VERBOSITY > 20)
      std::cout << "already in; old/new: " << old_distance << " " << new_distance << std::endl;
#endif
      if(new_distance > old_distance)
      {
        // the newer element is better, we must replace the old entry in the primal tetrahedra
#if (VERBOSITY > 20)
        std::cout << "found a farther dual point at : " << p_tetrahedron.dual_point()->index();
        std::cout << " (" << p_tetrahedron.dual_point()->point() << ")" << std::endl;
#endif
        // fixme what hint should be taken for the unordered set insertion ?
        PT_container_it hint = primal_tetrahedra.erase(is_insert_successful.first);
        primal_tetrahedra.insert(hint, p_tetrahedron);
      }
    }
  }

  void create_and_insert_primal_tetrahedron(const boost::array<const Canvas_point*, 4>& witnesses)
  {
    BTetrahedron tetrahedron;
    if(!check_primal_simplex(witnesses, tetrahedron))
      return;

    const Canvas_point* dual_point = compute_primal_tetrahedron_dual_point(witnesses);
    Primal_tetrahedron p_tetrahedron(tetrahedron, dual_point);
    insert_in_primal_tetrahedra(p_tetrahedron);
  }

  void add_to_primal_tetrahedra(typename std::vector<boost::array<const Canvas_point*, 4> >::const_iterator it,
                                typename std::vector<boost::array<const Canvas_point*, 4> >::const_iterator end)
  {
    for(; it!=end; ++it)
      create_and_insert_primal_tetrahedron(*it);
  }

  template<typename Candidates_set>
  void construct_primal_elements_from_candidates(Candidates_set& candidates)
  {
    boost::unordered_set<std::size_t> colors;
    for(typename Candidates_set::iterator it=candidates.begin(); it!=candidates.end(); ++it)
      colors.insert((*it)->closest_seed_id());

#if (VERBOSITY > 25)
    std::cout << "construct from canditate : ";
    typename Candidates_set::iterator it = candidates.begin();
    for(; it!=candidates.end(); ++it)
      std::cout << (*it)->index() << " [color: " << (*it)->closest_seed_id() << "], ";
    std::cout << std::endl;
#endif

#ifdef COMPUTE_PRIMAL_ALL_DIMENSIONS
    if(colors.size() >= 2)
    {
      // get the 2-combinations & build primal edges
      typedef std::vector<boost::array<const Canvas_point*, 2> > RType;
      RType combis = combinations<2>(candidates);
      add_to_primal_edges(combis.begin(), combis.end());
    }
#endif

    if(colors.size() >= 3)
    {
      // get the 3-combinations & build primal triangles
      typedef std::vector<boost::array<const Canvas_point*, 3> > RType;
      RType combis = combinations<3>(candidates);
      add_to_primal_triangles(combis.begin(), combis.end());
    }

    if(colors.size() >= 4)
    {
      // get the 4-combinations & build primal tetrahedra
      typedef std::vector<boost::array<const Canvas_point*, 4> > RType;
      RType combis = combinations<4>(candidates);
      add_to_primal_tetrahedra(combis.begin(), combis.end());
    }

    if(colors.size() >= 5)
      std::cout << "WARNING: cosphericity in the candidate set" << std::endl;
  }

  void detect_tetrahedra_self_intersections()
  {
    if(are_tetrahedra_intersection_computed)
      return;

    typedef typename Intersect_tetrahedra<Canvas>::Box                     Box;

    std::vector<Box> boxes;
    boxes.reserve(primal_tetrahedra.size());
    PTC_iterator it = primal_tetrahedra.begin();
    for(; it!=primal_tetrahedra.end(); ++it)
    {
      boxes.push_back(Box(seeds[it->simplex()[0]].bbox() +
                          seeds[it->simplex()[1]].bbox() +
                          seeds[it->simplex()[2]].bbox() +
                          seeds[it->simplex()[3]].bbox(), it));
    }

    std::vector<const Box*> box_ptr;
    box_ptr.reserve(primal_tetrahedra.size());
    for (typename std::vector<Box>::iterator bit=boxes.begin();
                                             bit!=boxes.end(); ++bit)
    {
      box_ptr.push_back(&*bit);
    }

    // todo check that cutoff parameter...
    CGAL::box_self_intersection_d(box_ptr.begin(), box_ptr.end(),
                                  Intersect_tetrahedra<Self>(this),
                                  std::ptrdiff_t(2000));

    // count them
    std::size_t intersected_tetrahedron = 0;
    it = primal_tetrahedra.begin();
    for(; it!=primal_tetrahedra.end(); ++it)
      if(it->m_is_intersected)
        ++intersected_tetrahedron;
    std::cout << intersected_tetrahedron << " intersected tetrahedra" << std::endl;

    are_tetrahedra_intersection_computed = true;
  }

  std::set<int> check_canvas_density_with_primals(const std::size_t ancestor_minimum_length = 8)
  {
#if (VERBOSITY > 10)
    std::cout << "density check" << std::endl;
#endif
    detect_tetrahedra_self_intersections();
    std::set<int> subdomain_indices;

    // we have taken the farthest witness point as dual information, we check
    // how many ancestors they have, if they have very few ancestors, the canvas
    // is most likely not dense enough for the seed set.

    PEC_iterator peit = primal_edges.begin();
    for(; peit!=primal_edges.end(); ++peit)
    {
      const Primal_edge& pt = *peit;
      const Canvas_point* dual_point = pt.dual_point();
      std::size_t ap_length = dual_point->count_ancestors();
      if(ap_length < ancestor_minimum_length)
        std::cout << "WARNING : the canvas is thin (" << ap_length
                  << ") for the primal : "<< pt << std::endl;
    }

    PTrC_iterator ptrit = primal_triangles.begin();
    for(; ptrit!=primal_triangles.end(); ++ptrit)
    {
      const Primal_triangle& pt = *ptrit;
      const Canvas_point* dual_point = pt.dual_point();
      std::size_t ap_length = dual_point->count_ancestors();
      if(ap_length < ancestor_minimum_length)
        std::cout << "WARNING : the canvas is thin (" << ap_length
                  << ") for the primal : " << pt << std::endl;
    }

    PTC_iterator ptit = primal_tetrahedra.begin();
    for(; ptit!=primal_tetrahedra.end(); ++ptit)
    {
      const Primal_tetrahedron& pt = *ptit;
      const Canvas_point* dual_point = pt.dual_point();
      std::size_t ap_length = dual_point->count_ancestors();
      if(ap_length < ancestor_minimum_length)
      {
        for(int i=0; i<4; ++i)
          subdomain_indices.insert(pt[i]);
        std::cout << "WARNING : the canvas is thin (" << ap_length
                  << ") for the primal : " << pt << std::endl;
      }
    }

    return subdomain_indices;
  }

  // ---------------------------------------------------------------------------
  // output related functions

  void output_geodesic_primal(std::ostream&) const
  {
    // todo when the definition of a geodesic tet is established...
  }

  void output_straight_primal(const std::string str_base) const
  {
    if(primal_edges.empty() && primal_triangles.empty() && primal_tetrahedra.empty())
    {
      std::cout << "WARNING: empty primal data structures at output" << std::endl;
      return;
    }

#if (VERBOSITY > 6)
    std::cout << "captured: ";
    std::cout << seeds.size() << " vertices ";
    std::cout << primal_edges.size() << " edges, ";
    std::cout << primal_triangles.size() << " triangles, ";
    std::cout << primal_tetrahedra.size() << " tetrahedra ";
    std::cout << " (" << this->seeds.size() << " seeds)" << std::endl;
#endif

    std::ofstream out((str_base + "_primal.mesh").c_str());
    std::ofstream out_bb((str_base + "_primal.bb").c_str());

    out << "MeshVersionFormatted 1" << std::endl;
    out << "Dimension 3" << std::endl;
    out << "Vertices" << std::endl;
    out << seeds.size() << std::endl;
    for(std::size_t i=0; i<seeds.size(); ++i)
      out << seeds[i] << " " << i+1 << std::endl;

#ifdef COMPUTE_PRIMAL_ALL_DIMENSIONS
    std::size_t bb_entries_n = primal_triangles.size() + primal_tetrahedra.size();
#else
    std::size_t bb_entries_n = primal_tetrahedra.size();
#endif
    out_bb << "3 1 " << bb_entries_n << " 1" << std::endl;

#ifdef COMPUTE_PRIMAL_ALL_DIMENSIONS
    // only output edges & triangles when we have computed them all
    // side note : if the macro isn't defined we only know _some_ of the triangles

    out << "Edges" << std::endl;
    out << primal_edges.size() << std::endl;
    for(typename Primal_edges_container::iterator it = primal_edges.begin();
                                                  it != primal_edges.end(); ++it)
    {
      const Primal_edge& edge = *it;
      out << edge[0]+1 << " " << edge[1]+1 << " 1" << std::endl;
    }

#ifdef USE_BRUTE_FORCE_TRIANGLE_INTERSECTIONS
    // build from the primal edges, the set of segments that will be used to check
    // self intersections of triangles
    std::set<typename K::Segment_3, Segment_3_comparer<K> > segments;
    for(typename Primal_edges_container::const_iterator it = primal_edges.begin();
                                                        it != primal_edges.end();
                                                        ++it)
    {
      const Primal_edge& edge = *it;
      const typename K::Segment_3 s(seeds[edge[0]], seeds[edge[1]]);
      segments.insert(s);
    }
#endif

    out << "Triangles" << std::endl;
    out << primal_triangles.size() << std::endl;
    for(typename Primal_triangles_container::iterator it = primal_triangles.begin();
                                                      it != primal_triangles.end();
                                                      ++it)
    {
      const Primal_triangle& tr = *it;
      out << tr[0]+1 << " " << tr[1]+1 << " " << tr[2]+1 << " ";
#ifdef USE_BRUTE_FORCE_TRIANGLE_INTERSECTIONS
      const typename K::Triangle_3 triangle(seeds[tr[0]], seeds[tr[1]], seeds[tr[2]]);
      out << is_triangle_intersected<K>(triangle, segments) << std::endl;
#else
      out << tr.m_is_intersected << std::endl;
#endif
      FT gamma = tr.compute_distortion(seeds.metrics());
      out_bb << gamma << std::endl;
    }
#endif

    out << "Tetrahedra" << std::endl;
    out << primal_tetrahedra.size() << std::endl;
    for(typename Primal_tetrahedra_container::iterator it = primal_tetrahedra.begin();
                                                       it != primal_tetrahedra.end();
                                                       ++it)
    {
      const Primal_tetrahedron& tet = *it;
      for(std::size_t i=0; i<4; ++i)
        out << tet[i] + 1 << " ";
      out << tet.m_is_intersected << std::endl;

      FT gamma = tet.compute_distortion(seeds.metrics());
      out_bb << gamma << std::endl;
    }
    out_bb << "End" << std::endl;
    out << "End" << std::endl;
  }

  void output_primal(const std::string str_base)
  {
    compute_primal();
    detect_tetrahedra_self_intersections();
    output_straight_primal(str_base);
  }

  void output_canvas_data_and_primal(const std::string str_base)
  {
    compute_primal();
    detect_tetrahedra_self_intersections();

    output_canvas(str_base);
    output_primal(str_base);
  }

  Canvas(const std::string& canvas_str_,
         const std::string& seeds_str_,
         const std::size_t max_seeds_n_,
         const Metric_field* mf_)
    :
      canvas_str(canvas_str_),
      canvas_points(),
      trial_points(),
      seeds(*this, seeds_str_, max_seeds_n_),
      mf(mf_),
      primal_edges(),
      primal_triangles(),
      primal_tetrahedra(),
      are_tetrahedra_intersection_computed(false),
      canvas_bbox(),
      known_count(0),
      trial_count(0),
      far_count(0)
  { }
};

}  // namespace Anisotropic_mesh_3
}  // namespace CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_CANVAS_H

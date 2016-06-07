#ifndef CGAL_ANISOTROPIC_MESH_3_GEODESIC_DRAWING_HELPER_H
#define CGAL_ANISOTROPIC_MESH_3_GEODESIC_DRAWING_HELPER_H

#include <CGAL/Canvas/canvas_enum.h>
#include <CGAL/Canvas/canvas_helper.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <CGAL/assertions.h>

#include <boost/unordered_set.hpp>

#include <iostream>
#include <map>
#include <vector>

// #define COMPUTE_GEODESICS_WITH_GRADIENT

namespace CGAL
{
namespace Anisotropic_mesh_3
{

// Bezier stuff
double evaluate_Berstein_poly(const std::size_t i, const std::size_t n,
                              const double t)
{
  // computes (i n) (1-t)^{n-i} t^i

  std::size_t c = combi(n, i);
  double a = std::pow((1-t), n-i);
  double b = std::pow(t, i);

//    std::cout << "eval berstein poly : " << i << " " << n << " " << t << std::endl;
//    std::cout << "c: " << c << " a: " << a << " b: " << b << " cab: " << c*a*b << std::endl;

  return c * a * b;
}

template<typename Canvas>
typename Canvas::Point_3
evaluate_Bezier_curve(const std::vector<typename Canvas::Point_3>& control_points,
                     const typename Canvas::FT t)
{
  typedef typename Canvas::FT          FT;
  typedef typename Canvas::Point_3     Point_3;

  // using de casteljau's algorithm
  std::size_t number_of_control_points = control_points.size();
  std::vector<FT> cp_xs(number_of_control_points);
  std::vector<FT> cp_ys(number_of_control_points);
  std::vector<FT> cp_zs(number_of_control_points);

  for(std::size_t i=0; i<number_of_control_points; ++i)
  {
    cp_xs[i] = control_points[i].x();
    cp_ys[i] = control_points[i].y();
    cp_zs[i] = control_points[i].z();
  }

  for(std::size_t i=1; i<number_of_control_points; ++i)
  {
    for(std::size_t j=0; j<number_of_control_points - i; ++j)
    {
      cp_xs[j] = (1-t) * cp_xs[j] + t * cp_xs[j+1];
      cp_ys[j] = (1-t) * cp_ys[j] + t * cp_ys[j+1];
      cp_zs[j] = (1-t) * cp_zs[j] + t * cp_zs[j+1];
    }
  }

  return Point_3(cp_xs[0], cp_ys[0], cp_zs[0]);
}

template<typename Canvas>
void compute_control_points(Canvas& canvas,
                            const std::deque<std::pair<std::size_t,
                                                       typename Canvas::FT> >& geodesic,
                            std::vector<typename Canvas::Point_3>& control_points)
{
  // the deque is made of points and a time between [0,1] obtained by knowing
  // that the speed is constant along a geodesic

  // A Bézier curve is B(t) = \sum_i (i n) (1-t)^{n-i} t^i P_i
  // and we know B(t0), \dots, B(tn). This is thus solving a system AX=B

  // the construction is definitely not the most efficient...

  typedef typename Canvas::FT             FT;
  typedef typename Canvas::Point_3        Point_3;

  std::size_t geo_path_size = geodesic.size();
  Eigen::MatrixXd A(geo_path_size, geo_path_size);
  Eigen::VectorXd B(geo_path_size);
  Eigen::VectorXd X(geo_path_size), Y(geo_path_size), Z(geo_path_size);

  // solve for the x coefficients
  std::size_t row_id = 0.;
  for(; row_id<geo_path_size; ++row_id)
  {
    const std::pair<std::size_t, FT>& pair = geodesic[row_id];
    const Point_3& p = canvas.get_point(pair.first).point();
    const FT time = pair.second;

    B(row_id) = p.x();
//    std::cout << "set_b@x: " << row_id << " " << p.x() << std::endl;

    for(std::size_t col_id=0; col_id < geo_path_size; ++col_id)
    {
      FT e = evaluate_Berstein_poly(col_id, geo_path_size-1, time);
//      std::cout << "set_a: " << row_id << " " << col_id << " " << e << std::endl;
      A(row_id, col_id) = e;
    }
  }
  X = A.fullPivLu().solve(B);

  // solve for the y coefficients (the matrix A stays the same)
  row_id = 0.;
  for(; row_id < geo_path_size; ++row_id)
  {
    const std::pair<std::size_t, FT>& pair = geodesic[row_id];
    const Point_3& p = canvas.get_point(pair.first).point();

    B(row_id) = p.y();
//    std::cout << "set_b@y: " << row_id << " " << p.y() << std::endl;
  }
  Y = A.fullPivLu().solve(B);

  // solve for the z coefficients (the matrix A stays the same)
  row_id = 0.;
  for(; row_id < geo_path_size; ++row_id)
  {
    const std::pair<std::size_t, FT>& pair = geodesic[row_id];
    const Point_3& p = canvas.get_point(pair.first).point();

    B(row_id) = p.z();
//    std::cout << "set_b@z: " << row_id << " " << p.z() << std::endl;
  }
  Z = A.fullPivLu().solve(B);

  // combine to obtain the control points !
  CGAL_precondition(control_points.size() == geodesic.size());

  std::cout << "solutions : " << std::endl << X.transpose()
                              << std::endl << Y.transpose()
                              << std::endl << Z.transpose() << std::endl;

  for(std::size_t i=0; i<geo_path_size; ++i)
  {
    Point_3 p(X(i), Y(i), Z(i));

//    std::cout << "control point i: " << i << " is " << p << std::endl;
    control_points[i] = p;
  }

  // debug (make sure that we go through the points that we wanted to interpolate)
  for(std::size_t i=0; i<geodesic.size(); ++i)
  {
    const std::pair<std::size_t, FT>& pair = geodesic[i];
    const Point_3& p = evaluate_Bezier_curve<Canvas>(control_points,
                                                     pair.second);

//      std::cout << "i/t: " << i << " " << pair.second << std::endl;
//      std::cout << "compare : " << p << " and " << points[pair.first].point << std::endl;
    CGAL_assertion(CGAL::squared_distance(p,
                                          canvas.get_point(pair.first).point()) < 1e-8);
  }
}

template<typename Canvas>
void compute_Bezier_curve(Canvas& canvas,
                          const std::deque<std::pair<std::size_t,
                                                     typename Canvas::FT> >& geodesic,
                          std::vector<typename Canvas::Point_3>& control_points)
{
  typedef typename Canvas::FT             FT;
  typedef typename Canvas::Point_3        Point_3;

  // todo simplify some notations maybe...
  std::cout << "compute Bezier approximation of the geodesic between: ";
  std::cout << canvas.get_point(geodesic.front().first).closest_seed_id() << " and ";
  std::cout << canvas.get_point(geodesic.back().first).closest_seed_id() << std::endl;

  // Before computing the control points, we might want to split the geodesic
  // in small segments of n points to limit the degree of the interpolating
  // Bezier curve
  CGAL_precondition(geodesic.size() == control_points.size());
  CGAL_precondition(geodesic.size() >= 2);

  std::size_t size_of_geodesic_segment = 4;
  std::size_t size_of_geodesic_path = geodesic.size();
  std::size_t start_of_segment = 0;
  std::size_t remaining_points = size_of_geodesic_path - 1;

  while(start_of_segment < size_of_geodesic_path - 1)
  {
    CGAL_assertion(remaining_points > 0); // can't make a segment if there is no remaining point
    if(remaining_points < size_of_geodesic_segment - 1) // can't make a full segment
      size_of_geodesic_segment = 1 + remaining_points;

    CGAL_assertion(size_of_geodesic_segment > 1); // proper segment

    std::cout << "total length of the path : " << size_of_geodesic_path << std::endl;
    std::cout << "Segment from : " << start_of_segment << " and size: " << size_of_geodesic_segment << std::endl;
    CGAL_assertion(start_of_segment + size_of_geodesic_segment - 1 < size_of_geodesic_path);

    // We want to go through the first and last point of the segment, thus we need
    // to locally scale up the time associated with these points
    FT start_of_segment_t = geodesic[start_of_segment].second;
    FT end_of_segment_t = geodesic[start_of_segment + size_of_geodesic_segment - 1].second;
    FT diff_t = end_of_segment_t - start_of_segment_t;
    FT scaling = 1. / diff_t;

    // fill the segment (todo make it lighter without copy)
    std::deque<std::pair<std::size_t, FT> > geodesic_segment;

    std::cout << "geodesic segment : " << std::endl;
    for(std::size_t i=0; i<size_of_geodesic_segment; ++i)
    {
      // same point but the time is scaled
      geodesic_segment.push_back(std::make_pair(geodesic[start_of_segment + i].first,
          scaling * (geodesic[start_of_segment + i].second - start_of_segment_t)));
      std::cout << "init: " << canvas.get_point(geodesic[start_of_segment + i].first).point()
                << " and time: " << geodesic[start_of_segment + i].second << std::endl;
      std::cout << "local : " << canvas.get_point(geodesic_segment[i].first).point()
                << " and time: " << geodesic_segment[i].second << std::endl;
    }

    std::vector<Point_3> segment_control_points(size_of_geodesic_segment);
    compute_control_points<Canvas>(canvas, geodesic_segment, segment_control_points);

    std::cout << "geodesic control points: " << std::endl;
    for(std::size_t i=0; i<size_of_geodesic_segment; ++i)
    {
      CGAL_assertion(start_of_segment + i < control_points.size());
      control_points[start_of_segment + i] = segment_control_points[i];
      std::cout << segment_control_points[i] << std::endl;
    }

    CGAL_postcondition(CGAL::squared_distance(control_points[start_of_segment],
              canvas.get_point(geodesic[start_of_segment].first).point()) < 1e-8);
    CGAL_postcondition(CGAL::squared_distance(control_points[start_of_segment + size_of_geodesic_segment - 1],
              canvas.get_point(geodesic[start_of_segment + size_of_geodesic_segment - 1].first).point()) < 1e-8);

    // -1 because the next segment starts at the end of the previous one
    start_of_segment += size_of_geodesic_segment - 1;

    // -1 because it's not useful to count the next starting point as part of the remaining points
    remaining_points -= size_of_geodesic_segment - 1;
  }

  for(std::size_t i=0; i<size_of_geodesic_path; ++i)
    CGAL_postcondition(control_points[i] != Point_3());
}

template<typename Canvas>
void output_Bezier_curves(const std::vector<std::vector<typename Canvas::Point_3> >& control_points,
                          std::size_t offset = 0,
                          const bool draw_control_points = false)
{
  // fixme if the geodesics are approximated by multiple low level Bezier curves
  // then what is drawn should be each curve independantly, not treating the
  // list of control points of different Bezier curves as if they were the
  // control points of a single high degree Bezier curve

  typedef typename Canvas::FT               FT;
  typedef typename Canvas::Point_3          Point_3;

  std::ofstream out("Bezier_test.mesh");
  std::size_t number_of_sample_points = 1000;
  FT step = 1./static_cast<FT>(number_of_sample_points - 1.);
  std::size_t number_of_curves = control_points.size();

  std::size_t total_vertices_n = number_of_curves * number_of_sample_points;
  std::size_t total_edges_n = number_of_curves * (number_of_sample_points - 1);

  if(draw_control_points)
  {
    for(std::size_t i=0; i<number_of_curves; ++i)
    {
      total_vertices_n += control_points[i].size();
      total_edges_n += control_points[i].size()-1;
    }
  }

  out << "MeshVersionFormatted 1" << '\n';
  out << "Dimension 3" << '\n';

  out << "Vertices" << '\n';
  out << total_vertices_n << std::endl;

  for(std::size_t i=0; i<number_of_curves; ++i)
  {
    const std::vector<Point_3>& local_control_points = control_points[i];

    if(draw_control_points)
    {
      for(std::size_t j=0; j<local_control_points.size(); ++j)
        out << local_control_points[j] << " " << i << std::endl;
    }

    for(std::size_t j=0; j<number_of_sample_points; ++j)
    {
      out << evaluate_Bezier_curve<Canvas>(local_control_points,
                                           step * static_cast<FT>(j))
          << " " << i << std::endl;
    }
  }

  out << "Edges" << '\n';
  out << total_edges_n << '\n';
  for(std::size_t i=0; i<number_of_curves; ++i)
  {
    if(draw_control_points)
    {
      std::size_t number_of_local_control_points = control_points[i].size();

      // edges to draw the control points
      for(std::size_t j=0; j<number_of_local_control_points-1; ++j)
        out << offset+j+1 << " " << offset+j+2 << " " << i << '\n'; // +1 for medit

      offset += number_of_local_control_points;
    }

    // edges to draw the bezier curve
    for(std::size_t j=0; j<number_of_sample_points-1; ++j)
      out << offset+j+1 << " " << offset+j+2 << " " << i << '\n'; // +1 for medit

    offset += number_of_sample_points;
  }

  out << "End" << std::endl;
}

template<typename Canvas>
void approximate_geodesics_with_Bezier(Canvas& canvas,
           std::vector<std::deque<std::pair<std::size_t,
                                            typename Canvas::FT> > >& geodesics)
{
  typedef typename Canvas::FT               FT;
  typedef typename Canvas::Point_3          Point_3;

  std::vector<std::vector<Point_3> > control_points;

  for(std::size_t i=0, gs=geodesics.size(); i<gs; ++i)
  {
    const std::deque<std::pair<std::size_t, FT> >& geodesic = geodesics[i];
    const std::size_t geo_size = geodesic.size();

    std::vector<Point_3> single_geodesic_control_points(geo_size);
    compute_Bezier_curve(canvas, geodesic, single_geodesic_control_points);

    std::cout << "full geodesic :" << std::endl;
    for(std::size_t j=0; j<geo_size; ++j)
    {
      const Point_3& p = canvas.get_point(geodesic[j].first).point();
      std::cout << "{" << p.x() << ", " << p.y() << ", " << p.z() << "}, " << std::endl;
    }

    std::cout << "full control points: " << std::endl;
    for(std::size_t j=0; j<geo_size; ++j)
    {
      const Point_3& p = single_geodesic_control_points[j];
      std::cout << "{" << p.x() << ", " << p.y() << ", " << p.z() << "}, " << std::endl;
    }

    control_points.push_back(single_geodesic_control_points);
  }

  CGAL_postcondition(geodesics.size() == control_points.size());

  output_Bezier_curves<Canvas>(control_points,
                               0, //points.size() /*offset*/,
                               false/*draw control points*/);
}

// not Bezier stuff

template<typename Canvas>
void add_to_neighboring_map(const typename Canvas::Primal_edge& e,
                            std::map<std::size_t,
                                     boost::unordered_set<std::size_t> >& m)
{
  typedef std::map<std::size_t, boost::unordered_set<std::size_t> >       Map;

  std::size_t min = e[0], max = e[1];
  CGAL_precondition(min != max);
  if(min > max)
  {
    min = e[1];
    max = e[0];
  }
  CGAL_postcondition(min < max);

  std::pair<typename Map::iterator, bool> is_insert_successful;

  boost::unordered_set<std::size_t> s;
  s.insert(max);

  is_insert_successful = m.insert(std::make_pair(min, s));
  if(!is_insert_successful.second)
    (is_insert_successful.first)->second.insert(max);
}

template<typename Canvas>
void build_neighboring_info(Canvas& canvas,
                            std::map<std::size_t,
                                     boost::unordered_set<std::size_t> >& neighbors)
{
  typedef typename Canvas::Primal_edge          PEdge;

  typename boost::unordered_set<PEdge>::iterator it = canvas.primal_edges.begin();
  typename boost::unordered_set<PEdge>::iterator end = canvas.primal_edges.end();

  for(; it!=end; ++it)
  {
    const PEdge& e = *it;
    add_to_neighboring_map<Canvas>(e, neighbors);
  }
}

template<typename Canvas>
std::size_t ancestor_with_depth(Canvas& canvas,
                                const typename Canvas::Canvas_point* cp)
{
  CGAL_precondition(cp);
  std::size_t depth = cp->depth();

  if(depth == 0)
    return -1;

  CGAL_precondition(depth > 0);
  std::size_t anc = cp->ancestor();
  --depth;

  while(depth > 0)
  {
    CGAL_precondition(anc != static_cast<std::size_t>(-1));
    anc = canvas.get_point(anc).ancestor();
    --depth;
  }

  CGAL_postcondition(anc != static_cast<std::size_t>(-1));
  return anc;
}

template<typename Canvas>
std::deque<std::pair<std::size_t, typename Canvas::FT> >
simplify_geodesic(const std::deque<std::pair<std::size_t,
                                             typename Canvas::FT> >& geodesic_path)
{
  typedef typename Canvas::FT FT;

  std::size_t max_point_n = 50;
  FT step = 1./static_cast<FT>(max_point_n - 1);
  FT curr_time = 0.;
  std::deque<std::pair<std::size_t, FT> > simplified_geodesic_path;

  // first point is naturally in
  simplified_geodesic_path.push_back(geodesic_path[0]);

  for(std::size_t i=0; i<geodesic_path.size(); ++i)
  {
    const std::pair<std::size_t, FT>& p = geodesic_path[i];
    const FT time = p.second;

    if(time - curr_time >= step)
    {
      // this point is sufficiently far from the previous one
      simplified_geodesic_path.push_back(p);
      curr_time = time;
    }
  }

  // need to make sure we have the final point of the geodesic in
  if(simplified_geodesic_path.back() != geodesic_path.back())
    simplified_geodesic_path.push_back(geodesic_path.back());

  std::cout << "simplified from " << geodesic_path.size()
            << " to " << simplified_geodesic_path.size() << std::endl;

  return simplified_geodesic_path;
}

template<typename Canvas>
void extract_geodesic(Canvas& canvas,
                      const typename Canvas::Canvas_point* cp,
                      std::vector<std::deque<std::pair<std::size_t,
                                                       typename Canvas::FT> > >& geodesics)
{
  typedef typename Canvas::FT FT;

  std::deque<std::pair<std::size_t, FT> > geodesic_path;
  FT total_length = cp->distance_to_closest_seed();
  geodesic_path.push_front(std::make_pair(cp->index(), 1.));
  CGAL_assertion(total_length != 0.);
  FT inv_total_length = 1./total_length;

#define EXTRACT_WITH_DEPTH
#ifdef EXTRACT_WITH_DEPTH
  std::size_t anc = ancestor_with_depth(canvas, cp);
#else
  std::size_t anc = cp->ancestor();
#endif
  CGAL_precondition(anc != static_cast<std::size_t>(-1));

  while(anc != static_cast<std::size_t>(-1))
  {
    FT t = canvas.get_point(anc).distance_to_closest_seed() * inv_total_length;
    if(t < 0) t = 0.;
    CGAL_postcondition(t >= 0. && t< 1.);

    geodesic_path.push_front(std::make_pair(anc, t));
#ifdef EXTRACT_WITH_DEPTH
    anc = ancestor_with_depth(canvas, &(canvas.get_point(anc)));
#else
    anc = canvas.get_point(anc).ancestor();
#endif
  }

  std::size_t s1 = cp->closest_seed_id();
  int s2 = cp->is_seed_holder();

  CGAL_precondition(s2 >= 0);
  CGAL_precondition(s1 < static_cast<std::size_t>(s2));

  std::cout << "built geodesic of " << s1 << " " << s2 << std::endl;
  geodesics.push_back(geodesic_path);
}

template<typename Canvas>
bool get_next_trial_point(Canvas& canvas,
                          typename Canvas::Canvas_point*& cp)
{
  typedef typename Canvas::Canvas_point Canvas_point;

  while(!canvas.trial_points.empty())
  {
    cp = canvas.trial_points.front();
    std::pop_heap(canvas.trial_points.begin(), canvas.trial_points.end(),
                  Canvas_point_comparer<Canvas_point>());
    canvas.trial_points.pop_back();

    CGAL_assertion(cp);

    if(cp->state() != TRIAL)
    {
      std::cout << "WARNING : point with a state non-TRIAL in the PQ : ";
      std::cout << cp->index() << "... Ignoring it!" << std::endl;
    }

    return true;
  }
  return false;
}

template<typename Canvas>
void spread_distances_from_one_seed(Canvas& canvas,
      const std::size_t seed_id,
      std::map<std::size_t,
               boost::unordered_set<std::size_t> > neighbors, // ugly but intentional copy
      std::vector<std::deque<std::pair<std::size_t,
                                       typename Canvas::FT> > >& geodesics,
      std::vector<typename Canvas::Canvas_point>& canvas_point_memory)
{
  // Spread distances from the seed_id, until it has reached its neighboring seeds
  CGAL_precondition(canvas.trial_points.empty());

  typedef typename Canvas::Canvas_point   Canvas_point;

  boost::unordered_set<std::size_t>& seeds_to_reach = neighbors[seed_id];

#if (verbose > 20)
  std::cout << "neighbors to reach : ";
  typename boost::unordered_set<std::size_t>::iterator it = seeds_to_reach.begin();
  typename boost::unordered_set<std::size_t>::iterator end = seeds_to_reach.end();
  for(; it!=end; ++it)
    std::cout << *it << " ";
  std::cout << std::endl;
#endif

  // (really) uglily for now, restart from the seed by searching the previous seed holder...
  for(std::size_t i=0, cps=canvas.canvas_points.size(); i<cps; ++i)
  {
    Canvas_point& cp = canvas.get_point(i);

    if(static_cast<std::size_t>(cp.is_seed_holder()) == seed_id)
    {
      CGAL_assertion(cp.closest_seed_id() == seed_id);
      CGAL_assertion(cp.ancestor() == static_cast<std::size_t>(-1));
      CGAL_assertion(cp.depth() == 0);

      cp.state() = TRIAL;
      canvas.trial_points.push_back(&cp);
      std::push_heap(canvas.trial_points.begin(), canvas.trial_points.end(),
                     Canvas_point_comparer<Canvas_point>());
      break;
    }
  }

  bool is_t_empty = canvas.trial_points.empty();
  CGAL_precondition(!is_t_empty);

  canvas.known_count = 0;
  Canvas_point* cp = NULL;
  while(!is_t_empty)
  {
    if(canvas.known_count%10000 == 0)
      canvas.print_states();

#if (verbose > 10)
    std::cout << "Trial queue size : " << canvas.trial_points.size() << std::endl;
#endif

#if (verbose > 40)
    std::cout << "trial heap: " << std::endl;
    for (std::vector<Canvas_point*>::iterator it = canvas.trial_points.begin();
         it != canvas.trial_points.end(); ++it)
      std::cout << (*it)->index << " " << (*it)->distance_to_closest_seed() << std::endl;
    std::cout << std::endl;
#endif

    if(!get_next_trial_point(canvas, cp))
      break;

    CGAL_postcondition(cp);

#if (verbose > 10)
    std::cout << "-------------------------------------------------------" << std::endl;
    std::cout << "picked n° " << cp->index() << " (" << cp->point() << ") "
              << "at distance : " << cp->distance_to_closest_seed()
              << " from the seed " << cp->closest_seed_id()
              << " ancestor : " << cp->ancestor() << std::endl;
#endif

    // Check if we have reached one of the required seeds
    if(cp->is_seed_holder() >= 0) // really ugly, fixme (see comment at member def)
    {
      std::cout << "found the seed " << cp->is_seed_holder()
                << " at " << cp->index() << " (" << cp->point() << ") " << std::endl;

      // That's the center of a Voronoi cell (no ancestor)
      std::size_t reached_seed_id = cp->is_seed_holder();
      CGAL_assertion(cp->closest_seed_id() == seed_id);

      // Geodesic are symmetrical, no need to compute it both ways.
      // This also dodges the first point of the trial queue (which will be the
      // seed_holder for the seed 'seed_id' and we don't care about it
      if(reached_seed_id < seed_id)
        continue;

      boost::unordered_set<std::size_t>::iterator sit =
                                          seeds_to_reach.find(reached_seed_id);
      if(sit != seeds_to_reach.end())
      {
        // the seed that we have reached is part of the neighbors we want to reach
        extract_geodesic(canvas, cp, geodesics);
        seeds_to_reach.erase(sit);
      }
    }

#if 1//ndef COMPUTE_GEODESICS_WITH_GRADIENT
    // only spread as little as possible
    if(seeds_to_reach.empty())
      break;
#endif

    cp->state() = KNOWN;
    canvas.known_count++; // tmp --> use change_state()

    // Change the distance_to_closest_seed of points to allow for seed_id's
    // cell to spread
    typename Canvas::Vertex_handle_handle it =
                                       cp->adjacent_vertices_in_complex_begin();
    typename Canvas::Vertex_handle_handle end =
                                       cp->adjacent_vertices_in_complex_end();
    for(; it!=end; ++it)
    {
      typename Canvas::Vertex_handle v = *it;
      Canvas_point& cq = canvas.get_point(v->info());

      if(cq.state() == FAR)
      {
        canvas_point_memory.push_back(cq);
        cq.distance_to_closest_seed() = FT_inf;
      }
    }

    cp->update_neighbors_distances(canvas.trial_points);

    // always rebuild the priority queue
    std::make_heap(canvas.trial_points.begin(), canvas.trial_points.end(),
                   Canvas_point_comparer<Canvas_point>());
    is_t_empty = canvas.trial_points.empty();
  }

  std::cout << "reached all the seeds in " << canvas.known_count << " points" << std::endl;

  // make sure we've reached all the neighboring seeds
  CGAL_postcondition(seeds_to_reach.empty());
}

template<typename Canvas>
void rollback_points(Canvas& canvas,
                     const std::vector<typename Canvas::Canvas_point>& cp_memory)
{
  typedef typename Canvas::Canvas_point Canvas_point;

  std::cout << "rollback" << std::endl;
  for(std::size_t i=0, cpms=cp_memory.size(); i<cpms; ++i)
  {
    const Canvas_point& cp_m = cp_memory[i];
    std::size_t id = cp_m.index();
    CGAL_assertion(canvas.get_point(id).point() == cp_m.point());
    canvas.get_point(id) = cp_m;
  }

  canvas.trial_points.clear();
}

template<typename Canvas>
void output_geodesics(Canvas& canvas,
                      const std::vector<std::deque<std::pair<std::size_t,
                                                             typename Canvas::FT> > >& geodesics,
                      const std::string str_base)
{
  typedef typename Canvas::FT FT;

  if(geodesics.empty())
  {
    std::cout << "no geodesic to output" << std::endl;
    return;
  }

  std::ofstream out((str_base + "_geodesics.mesh").c_str());

  out << "MeshVersionFormatted 1" << '\n';
  out << "Dimension 3" << '\n';

  out << "Vertices" << '\n';
  out << canvas.canvas_points.size() << '\n';
  for(std::size_t i=0, cps=canvas.canvas_points.size(); i<cps; ++i)
    out << canvas.get_point(i).point() << " " << i+1 << '\n';

  std::size_t edges_n = 0.;
  for(std::size_t i=0, gs=geodesics.size(); i<gs; ++i)
    edges_n += geodesics[i].size() - 1;

  out << "Edges" << '\n';
  out << edges_n << '\n';
  for(std::size_t i=0, gs=geodesics.size(); i<gs; ++i)
  {
    const std::deque<std::pair<std::size_t, FT> >& geo_path = geodesics[i];
    std::size_t geo_path_size_m1 = geo_path.size()-1;
    for(std::size_t j=0; j<geo_path_size_m1; ++j)
    {
      const std::pair<std::size_t, FT>& p = geo_path[j];
      const std::pair<std::size_t, FT>& pp1 = geo_path[j+1];

      // +1 due to medit
      out << p.first + 1 << " " << pp1.first + 1 << " 5" << std::endl;
    }
  }

  out << "Triangles" << std::endl;
  out << canvas.m_c3t3.number_of_facets() << std::endl;
  typename Canvas::Facet_iterator fit = canvas.m_c3t3.facets_in_complex_begin();
  typename Canvas::Facet_iterator fend = canvas.m_c3t3.facets_in_complex_end();
  for(; fit!=fend; ++fit)
  {
    typename Canvas::Cell_handle c = fit->first;
    std::size_t sec = fit->second;

    boost::unordered_set<std::size_t> materials;
    for(std::size_t j=1; j<4; ++j)
      materials.insert(canvas.get_point(c->vertex((sec+j)%4)->info()).closest_seed_id());

#ifdef ONLY_PRINT_BISECTORS
    if(materials.size() == 1)
      continue;
#endif

    for(std::size_t j=1; j<=3; ++j)
      out << c->vertex((sec+j)%4)->info() + 1 << " ";

    std::size_t mat = (materials.size() == 1) ? 2 + (*(materials.begin())) :
                                                1;
    out << mat << std::endl;
  }

  out << "End" << std::endl;
}


template<typename Canvas>
void compute_gradient_at_facet(Canvas& canvas,
                               typename Canvas::Cell_handle c,
                               std::size_t sec,
                               std::vector<Eigen::Vector3d>& vertex_gradients)
{
  // compute the grad on facet f = (c,sec) and update each vertex with
  // 1/d(p_i, g) grad(f)

  typedef typename Canvas::FT                   FT;
  typedef typename Canvas::Point_3              Point_3;
  typedef typename Canvas::Vertex_handle        Vertex_handle;

  // get vertices
  boost::array<Vertex_handle, 3> vertices;
  for(std::size_t i=0; i<3; ++i)
    vertices[i] = c->vertex((sec + i + 1) % 4);

  // compute normal
  boost::array<Eigen::Vector3d, 3> edges; // e01, e12, e20

  for(std::size_t i=0; i<3; ++i)
  {
    edges[i](0) = vertices[(i+1)%3]->point().x() - vertices[i]->point().x();
    edges[i](1) = vertices[(i+1)%3]->point().y() - vertices[i]->point().y();
    edges[i](2) = vertices[(i+1)%3]->point().z() - vertices[i]->point().z();
  }

  Eigen::Vector3d n = edges[0].cross(-edges[2]);
  n.normalize();

  // compute area of the triangle
  typename Canvas::Kernel::Compute_area_3 o;
  FT A = std::abs(o(vertices[0]->point(),
                    vertices[1]->point(),
                    vertices[2]->point()));

  // compute the gradient
  Eigen::Vector3d grad = Eigen::Vector3d::Zero();

  for(std::size_t i=0; i<3; ++i)
  {
    // edge i is e_{i,(i+1)%3}
    Eigen::Vector3d cp = n.cross(edges[i]);

    // value at the opposite vertex (i+2)%3
    FT val = canvas.get_point(vertices[(i+2)%3]->info()).distance_to_closest_seed();
    grad += val * cp;
  }

  grad *= 0.5 / A;

  // center of the triangle
  FT third = 1./3.;
  Point_3 bg = CGAL::barycenter(vertices[0]->point(), third,
                                vertices[1]->point(), third,
                                vertices[2]->point(), third);

  // add it to the vertex grad parts at each vertex
  typename Canvas::Kernel::Compute_squared_distance_3 sqd;
  for(std::size_t i=0; i<3; ++i)
  {
    FT d = sqd(bg, vertices[i]->point());
    CGAL_assertion(d != 0.);

    std::size_t id = vertices[i]->info();

    vertex_gradients[id] += grad / A;
  }
}

template<typename Canvas, typename Tree>
void trace_geodesic(Canvas& canvas,
                    const Tree& tree,
                    const std::vector<Eigen::Vector3d>& vertex_gradients,
                    std::vector<std::deque<typename Canvas::Point_3> >& geodesics,
                    std::size_t starting_id)
{
  // tracing the geodesic from starting id to the seed
  CGAL_precondition(starting_id < canvas.canvas_points.size());

  typedef typename Canvas::FT                   FT;
  typedef typename Canvas::Point_3              Point_3;
  typedef typename Canvas::Facet_iterator       Facet_iterator;

  typedef typename Canvas::Vector3d             Vector3d;

  typedef typename Tree::Point_and_primitive_id Point_and_primitive_id;

  FT tau = 0.01; // gradient step

  bool is_seed_reached = false;
  std::size_t iterations = 0.;
  Point_3 current_point = canvas.get_point(starting_id).point();

  std::deque<Point_3> geodesic;
  FT dist_mem = 1e30;

  while(!is_seed_reached)
  {
    //find closest point on the surface
    Point_and_primitive_id pp = tree.closest_point_and_primitive(current_point);

    Point_3 closest_point = pp.first;
    Facet_iterator f = pp.second; // closest primitive id

    geodesic.push_front(closest_point);

    typename Canvas::Cell_handle c = f->first;
    std::size_t second = f->second;

    // get gradient & value at point p
    Point_3 p = c->vertex((second + 1)%4)->point();
    Point_3 q = c->vertex((second + 2)%4)->point();
    Point_3 r = c->vertex((second + 3)%4)->point();
    std::size_t p_id = c->vertex((second + 1)%4)->info();
    std::size_t q_id = c->vertex((second + 2)%4)->info();
    std::size_t r_id = c->vertex((second + 3)%4)->info();

    // compute the bary weights in the face
    Vector3d v1, v2, v3;
    v1(0) = q.x() - p.x();
    v1(1) = q.y() - p.y();
    v1(2) = q.z() - p.z();
    v2(0) = r.x() - p.x();
    v2(1) = r.y() - p.y();
    v2(2) = r.z() - p.z();

    FT lambda_p, lambda_q, lambda_r;
    v3(0) = closest_point.x() - p.x();
    v3(1) = closest_point.y() - p.y();
    v3(2) = closest_point.z() - p.z();

    FT d11 = v1.dot(v1);
    FT d12 = v1.dot(v2);
    FT d22 = v2.dot(v2);
    FT d31 = v3.dot(v1);
    FT d32 = v3.dot(v2);
    FT den = 1. / (d11 * d22 - d12 * d12);
    lambda_q = (d22 * d31 - d12 * d32) * den;
    lambda_r = (d11 * d32 - d12 * d31) * den;
    lambda_p = 1. - lambda_q - lambda_r;

    CGAL_precondition(std::abs(vertex_gradients[p_id].norm() - 1.) < 1e-10);
    CGAL_precondition(std::abs(vertex_gradients[q_id].norm() - 1.) < 1e-10);
    CGAL_precondition(std::abs(vertex_gradients[r_id].norm() - 1.) < 1e-10);

    Vector3d grad_at_closest_point = lambda_p * vertex_gradients[p_id] +
                                     lambda_q * vertex_gradients[q_id] +
                                     lambda_r * vertex_gradients[r_id];


//    std::cout << "closest point is: " << closest_point << std::endl;

    FT dist_to_seed = lambda_p * canvas.get_point(p_id).distance_to_closest_seed() +
                      lambda_q * canvas.get_point(q_id).distance_to_closest_seed() +
                      lambda_r * canvas.get_point(r_id).distance_to_closest_seed();

//    std::cout << "distance to seed : " << dist_to_seed << std::endl;

    // following doesn't solve cycles
    if(std::abs(dist_to_seed - dist_mem) < 1e-5) // hardcode fixme
      break;

    dist_mem = dist_to_seed;

    // new point
    Vector3d cur_p;
    cur_p(0) = current_point.x();
    cur_p(1) = current_point.y();
    cur_p(2) = current_point.z();

//    std::cout << "grad: " << grad_at_closest_point << std::endl;

    Vector3d new_p = cur_p - tau * grad_at_closest_point;

    // point below is not on the surface, it'll be projected next loop...
    current_point = Point_3(new_p(0), new_p(1), new_p(2));

    iterations++;
    if(iterations > 1000)
    {
      std::cout << "time to debug stuff" << std::endl;
      break;
    }
  }

  CGAL_postcondition(!geodesic.empty());
  geodesics.push_back(geodesic);
}

template<typename Canvas>
struct Base_mesh_primitive
{
  typedef typename Canvas::Kernel           K;
  typedef typename K::Point_3               Point;
  typedef typename K::Triangle_3            Datum;

  typedef typename Canvas::Facet_iterator   Id;

  Id m_id;
  Datum m_datum; // cache the datum
  Point m_p; // cache the reference point

  Base_mesh_primitive() { }
  Base_mesh_primitive(Id i, const Datum& datum, const Point& point)
    :
      m_id(i),
      m_datum(datum),
      m_p(point)
  { }

   Id id() const { return m_id; }
   Datum datum() const { return m_datum; }
   Point reference_point() const { return m_p; }
};

template<typename Canvas, typename Tree>
void build_facet_aabb_tree(Canvas& canvas,
                           Tree& tree)
{
  typedef typename Canvas::Kernel   K;
  typedef typename K::Point_3                            Point_3;
  typedef typename K::Triangle_3                         Triangle_3;

  typename Canvas::Facet_iterator fit = canvas.m_c3t3.facets_in_complex_begin();
  typename Canvas::Facet_iterator fend = canvas.m_c3t3.facets_in_complex_end();
  for(; fit!=fend; ++fit)
  {
    typename Canvas::Cell_handle c = fit->first;
    std::size_t sec = fit->second;

    Point_3 p = c->vertex((sec + 1)%4)->point();
    Triangle_3 t(p,
                 c->vertex((sec + 2)%4)->point(),
                 c->vertex((sec + 3)%4)->point());
    typename Tree::Primitive pri(fit, t, p);
    tree.insert(pri);
  }
}

template<typename Canvas>
void draw_real_geodesics(const std::vector<std::deque<typename Canvas::Point_3> >& geodesics,
                         const std::string str_base)
{
  typedef typename Canvas::Point_3          Point_3;

  if(geodesics.empty())
  {
    std::cout << "no geodesic to output" << std::endl;
    return;
  }

  // count vertices and edges
  std::size_t vertices_n = 0.;
  std::size_t edges_n = 0.;
  std::size_t offset = 0;

  for(std::size_t i=0, gs=geodesics.size(); i<gs; ++i)
  {
    vertices_n += geodesics[i].size();
    edges_n += geodesics[i].size() - 1;
  }

  std::ofstream out((str_base + "_real_geodesics.mesh").c_str());

  out << "MeshVersionFormatted 1" << '\n';
  out << "Dimension 3" << '\n';

  out << "Vertices" << '\n';
  out << vertices_n << '\n';
  for(std::size_t i=0, gss=geodesics.size(); i<gss; ++i)
  {
    const std::deque<Point_3>& geodesic = geodesics[i];
    for(std::size_t j=0, gs=geodesic.size(); j<gs; ++j)
    {
      out << geodesic[j] << " " << offset + j + 1 << '\n';
    }
    offset += geodesic.size();
  }

  offset = 0;

  out << "Triangles" << '\n';
  out << edges_n << '\n';
  for(std::size_t i=0, gss=geodesics.size(); i<gss; ++i)
  {
    const std::deque<Point_3>& geodesic = geodesics[i];
    for(std::size_t j=0, gs=geodesic.size()-1; j<gs; ++j)
    {
      out << offset + j + 1 << " " << offset + j + 2 << " " << offset + j + 1 << " " << i << '\n';
    }
    offset += geodesic.size();
  }

  out << "End" << std::endl;
}

template<typename Canvas>
void draw_gradient(const Canvas& canvas,
                   const std::vector<Eigen::Vector3d>& vertex_gradients)
{
  std::ofstream out("geodesics_gradient.mesh");

  out << "MeshVersionFormatted 1" << '\n';
  out << "Dimension 3" << '\n';

  out << "Vertices" << '\n';
  out << 2 * canvas.canvas_points.size() << '\n';
  for(std::size_t i=0, cps=canvas.canvas_points.size(); i<cps; ++i)
  {
    typename Canvas::Vector_3 v(vertex_gradients[i](0),
                                vertex_gradients[i](1),
                                vertex_gradients[i](2));

    out << canvas.get_point(i).point() << " " << 2*i+1 << '\n';
    out << canvas.get_point(i).point() + 0.01 * v << " " << 2*(i+1) << '\n';
  }

  out << "Triangles" << '\n';
  out << canvas.canvas_points.size() << '\n';

  for(std::size_t i=0, cps=canvas.canvas_points.size(); i<cps; ++i)
    out << 2*i+1 << " " << 2*(i+1) << " " << 2*i+1 << " " << i << std::endl;

  out << "End" << std::endl;
}

template<typename Canvas, typename Tree>
void trace_geodesic_with_gradient(Canvas& canvas,
        std::size_t seed_id,
        std::vector<std::deque<typename Canvas::Point_3> >& geodesics,
        std::map<std::size_t/*seed_id*/,
                 boost::unordered_set<std::size_t>/*neighbor ids*/> neighbors,
        const Tree& tree,
        const std::string str_base)
{
  std::size_t v_n = canvas.canvas_points.size();

  // compute gradient at all triangles/vertices of the canvas
  std::vector<Eigen::Vector3d> vertex_gradients(v_n, Eigen::Vector3d::Zero());

  typename Canvas::Facet_iterator fit = canvas.m_c3t3.facets_in_complex_begin();
  typename Canvas::Facet_iterator fend = canvas.m_c3t3.facets_in_complex_end();
  for(; fit!=fend; ++fit)
  {
    typename Canvas::Cell_handle c = fit->first;
    std::size_t sec = fit->second;

    compute_gradient_at_facet(canvas, c, sec, vertex_gradients);
  }

  // normalize gradients
  for(std::size_t i=0; i<v_n; ++i)
    vertex_gradients[i].normalize();

  draw_gradient(canvas, vertex_gradients);
  canvas.output_canvas_data_and_primal("test");

  // trace geodesic
  const boost::unordered_set<std::size_t>& seeds_to_reach = neighbors[seed_id];

  typename boost::unordered_set<std::size_t>::iterator it = seeds_to_reach.begin();
  typename boost::unordered_set<std::size_t>::iterator end = seeds_to_reach.end();
  for(; it!=end; ++it)
  {
    std::size_t start_seed_id = *it;

    // obviously really ugly
    std::size_t starting_id = -1;
    for(std::size_t i=0, cps=canvas.canvas_points.size(); i<cps; ++i)
    {
      if(canvas.get_point(i).is_seed_holder() == start_seed_id)
      {
        starting_id = i;
        break;
      }
    }

    CGAL_postcondition(starting_id != static_cast<std::size_t>(-1));
    CGAL_postcondition(canvas.get_point(starting_id).ancestor() == static_cast<std::size_t>(-1));

    std::cout << "trying to reach the seed at : " << canvas.get_point(seed_id).point();
    std::cout << " from seed at " << canvas.get_point(starting_id).point() << std::endl;

    trace_geodesic(canvas, tree, vertex_gradients, geodesics, starting_id);
  }
}

template<typename Canvas>
void compute_geodesics(Canvas& canvas,
                       const std::string str_base)
{
  typedef typename Canvas::Kernel K;
  typedef typename K::FT FT;
  typedef typename K::Point_3 Point_3;

  std::cout << "computing geodesics" << std::endl;

  // disgusting hack
  for(std::size_t i=0, cps=canvas.canvas_points.size(); i<cps; ++i)
    canvas.get_point(i).ignore_children = true;

  // For each seed, we must spread to the nearest neighbors to grab the geodesic
  // between both seeds
  CGAL_precondition(!canvas.primal_edges.empty());

  std::map<std::size_t/*seed_id*/,
           boost::unordered_set<std::size_t>/*neighbor ids*/> neighbors;
  build_neighboring_info(canvas, neighbors);

  std::cout << "built neighboring info : " << neighbors.size() << std::endl;

  std::vector<std::deque<std::pair<std::size_t/*id*/,
                                   FT/*time*/> > >          geodesics;

#ifdef COMPUTE_GEODESICS_WITH_GRADIENT
  std::vector<std::deque<Point_3> > real_geodesics;

  // build aabb tree of facets used to project points on the surface
  typedef Base_mesh_primitive<Canvas>             Primitive;
  typedef CGAL::AABB_traits<K, Primitive>         AABB_triangle_traits;
  typedef CGAL::AABB_tree<AABB_triangle_traits>   Tree;

  Tree tree;
  build_facet_aabb_tree(canvas, tree);
#endif

  for(std::size_t seed_id=0, ss=canvas.seeds.size(); seed_id<ss; ++seed_id)
  {
    std::cout << "computing geodesics from : " << seed_id << std::endl;

    // below is kind of ugly... but it's clearer than switching everything to
    // 'FAR' earlier and then having the memory overwrite 'KNOWN' with 'FAR'
    // every loop iteration
    for(std::size_t i=0, ptss=canvas.canvas_points.size(); i<ptss; ++i)
      canvas.get_point(i).state() = FAR;

    // since we're spreading from a seed over the other cells, we must keep
    // the previous colors in memory
    std::vector<typename Canvas::Canvas_point> canvas_point_memory;

    spread_distances_from_one_seed(canvas, seed_id, neighbors,
                                   geodesics, canvas_point_memory);

#ifdef COMPUTE_GEODESICS_WITH_GRADIENT
    trace_geodesic_with_gradient(canvas, seed_id, real_geodesics, neighbors, tree, str_base);
#endif

    rollback_points(canvas, canvas_point_memory); // reset the points we have overwritten
  }

#ifdef COMPUTE_GEODESICS_WITH_GRADIENT
  draw_real_geodesics<Canvas>(real_geodesics, str_base);
#endif

  output_geodesics(canvas, geodesics, str_base);

  // simplify geodesics
  for(std::size_t i=0; i<geodesics.size(); ++i)
  {
//      geodesics[i] = simplify_geodesic(geodesics[i]);
  }

//  approximate_geodesics_with_Bezier(canvas, geodesics);
}

} // Anisotropic_mesh_3
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_GEODESIC_DRAWING_HELPER_H

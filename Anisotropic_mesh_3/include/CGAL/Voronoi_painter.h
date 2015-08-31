#ifndef CGAL_ANISOTROPIC_MESH_3_PAINTER_HELPER_3_H
#define CGAL_ANISOTROPIC_MESH_3_PAINTER_HELPER_3_H

#include <CGAL/Starset.h>
#include <CGAL/Stretched_Delaunay_3.h>

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Implicit_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/make_mesh_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <boost/array.hpp>

#include <limits>
#include <set>
#include <vector>
#include <utility>

#define FILTER_SEEDS_OUTSIDE_CANVAS
#define verbose 2

//fixme with a proper FT
const double FT_inf = std::numeric_limits<double>::infinity();

// tiny helper functions
std::size_t fact(std::size_t n)
{
  return (n==0)?1:n*fact(n-1);
}

template<std::size_t k_max>
std::vector<boost::array<std::size_t, k_max> >
combinations(const std::set<std::size_t>& set,
             std::set<std::size_t>::const_iterator it,
             const std::size_t k)
{
  // the idea is to build the combinations recursively :
  // at one level we pick one element, and we'll then build the combination
  // as : this element added to the k-1 sized combinations

  CGAL_assertion(set.size() >= k_max && k > 0);
//    if(k_max == k)
//      combis.reserve(fact(n)/(fact(k)*fact(n-k))); // is there really a point?

  std::vector<boost::array<std::size_t, k_max> > combis; // current level combis

  std::set<std::size_t>::const_iterator sit = it;
  std::set<std::size_t>::const_iterator send = set.end();

  // we still need to pick k-1 elements afterwards so we need to keep some space
  std::advance(send, -(k-1));

  for(; sit!=send; ++sit)
  {
    if(k == 1)
    {
      boost::array<std::size_t, k_max> combi;

      //this is only useful to debug, but it should be bug free now!
//        for(std::size_t i=0; i<k_max; ++i)
//          combi[i] = static_cast<std::size_t>(-1);

      // initialize the last element
      combi[k_max-k] = *sit;
      combis.push_back(combi);
    }
    else
    {
      // initialize the k_max-k element
      std::vector<boost::array<std::size_t, k_max> > lower_level_combis =
                                               combinations<k_max>(set, ++sit, k-1);
      --sit;
      for(std::size_t i=0; i<lower_level_combis.size(); ++i)
        lower_level_combis[i][k_max-k] = *(sit);

      combis.insert(combis.end(), lower_level_combis.begin(),
                                  lower_level_combis.end());
    }
  }
  return combis;
}

template<typename K>
struct Segment_3_comparer
{
  typedef typename K::Point_3       Point;
  typedef typename K::Segment_3     Segment;

  bool operator() (const Segment left, const Segment right)
  {
    const Point& left_o = left.source();
    const Point& left_t = left.target();
    const Point& right_o = right.source();
    const Point& right_t = right.target();

    for(int i=0; i<3; ++i)
    {
      if(left_o[i] < right_o[i])
        return true;
      if(left_o[i] > right_o[i])
        return false;
    }

    for(int i=0; i<2; ++i)
    {
      if(left_t[i] < right_t[i])
        return true;
      if(left_t[i] > right_t[i])
        return false;
    }

    return left_t[2] < right_t[2];
  }
};

template<typename K, typename KExact>
bool is_triangle_intersected(const typename K::Triangle_3& triangle,
                             const std::set<typename K::Segment_3,
                                            Segment_3_comparer<K> >& segments)
{
  typedef typename K::Point_3                                  Point;
  typedef typename K::Segment_3                                Segment;

  typedef typename KExact::Point_3                             EPoint;
  typedef typename KExact::Segment_3                           ESegment;
  typedef typename KExact::Triangle_3                          ETriangle;

  typedef CGAL::Cartesian_converter<K, KExact>                 To_exact;
  To_exact to_exact;

  bool is_intersected = false;

  // need to switch to epeck for the correct intersections here...
  const Point& p0 = triangle[0];
  const Point& p1 = triangle[1];
  const Point& p2 = triangle[2];

//#pragma omp parallel shared(is_intersected, p0, p1, p2)
  {
    typename std::set<Segment>::const_iterator it = segments.begin();
    for(; it != segments.end(); ++it)
    {
#pragma omp single nowait // hack to parallelize the std::set
      {
#pragma omp flush (is_intersected)
        if(!is_intersected) // hack because we're not allowed to break (or continue)
          // inside a pragma omp for
        {
          const Segment& segment = *it;
          bool local_intersected = false;

          ESegment esegment = to_exact(segment);
          ETriangle etriangle = to_exact(triangle);
          typename CGAL::cpp11::result_of<typename KExact::Intersect_3(ESegment,
                                                                       ETriangle)>::type
              result = CGAL::intersection(esegment, etriangle);

          if (result)
          {
            if (const EPoint* p = boost::get<EPoint>(&*result))
            {
              const EPoint& ep = *p;
              local_intersected = (ep != to_exact(p0) &&
                                   ep != to_exact(p1) &&
                                   ep != to_exact(p2));
            }
            else if(const ESegment* s = boost::get<ESegment>(&*result))
            {
              const EPoint& ep0 = to_exact(p0);
              const EPoint& ep1 = to_exact(p1);
              const EPoint& ep2 = to_exact(p2);
              local_intersected = ((s->source() != ep0 &&
                                    s->source() != ep1 &&
                                    s->source() != ep2) ||
                                   (s->target() != ep0 &&
                                    s->target() != ep1 &&
                                    s->target() != ep2));
            }
          }

#pragma omp critical
          {      if(local_intersected)
              is_intersected = true;
          }
#pragma omp flush (is_intersected)
        } // end !is_intersected
      } // end single nowait
    } // end for
  } // end parallel
  return is_intersected;
}

namespace CGAL
{
// stuff to generate the canvas using Mesh_3 since Aniso_mesh_3 is a turtle.
// Need to define our own refinement criterion based on metric shenanigans

template<typename K>
typename K::FT cube_function(const typename K::Point_3& p)
{
  if( p.x() >= 0 && p.x() <= 2 &&
      p.y() >= 0 && p.y() <= 2 &&
      p.z() >= 0 && p.z() <= 2 )
    return -1.;
  return 1.;
}

template<typename K>
void create_polylines(std::list<std::vector<typename K::Point_3> >& polylines)
{
  typedef typename K::FT                      FT;
  typedef typename K::Point_3                 Point_3;
  typedef std::vector<Point_3>                Polyline_3;

  FT side_x = 2.; // fixme
  FT side_y = 2.;
  FT side_z = 2.;

  //bot 4 edges
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(0,0,0));
    polyline.push_back(Point_3(side_x,0,0));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(side_x,0,0));
    polyline.push_back(Point_3(side_x,side_y,0));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(side_x,side_y,0));
    polyline.push_back(Point_3(0,side_y,0));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(0,side_y,0));
    polyline.push_back(Point_3(0,0,0));
    polylines.push_back(polyline);
  }

  // vertical edges
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(0,0,0));
    polyline.push_back(Point_3(0,0,side_z));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(side_x,0,0));
    polyline.push_back(Point_3(side_x,0,side_z));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(0,side_y,0));
    polyline.push_back(Point_3(0,side_y,side_z));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(side_x,side_y,0));
    polyline.push_back(Point_3(side_x,side_y,side_z));
    polylines.push_back(polyline);
  }

  // top 4 edges
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(0,0,side_z));
    polyline.push_back(Point_3(side_x,0,side_z));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(side_x,0,side_z));
    polyline.push_back(Point_3(side_x,side_y,side_z));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(side_x,side_y,side_z));
    polyline.push_back(Point_3(0,side_y,side_z));
    polylines.push_back(polyline);
  }
  {
    Polyline_3 polyline;
    polyline.push_back(Point_3(0,side_y,side_z));
    polyline.push_back(Point_3(0,0,side_z));
    polylines.push_back(polyline);
  }
}

template <typename Tr>
class Mesh_edge_criteria_w_distortion_3 :
    public Mesh_edge_criteria_3<Tr>

{

};

template<typename Tr>
class Mesh_facet_criteria_w_distortion_3 :
    public Mesh_facet_criteria_3<Tr>
{

};

template<typename Tr>
class Mesh_cell_criteria_w_distortion_3 :
    public Mesh_cell_criteria_3<Tr>
{

};

template<typename K>
void generate_canvas()
{
  typedef typename K::FT                                         FT;
  typedef typename K::Point_3                                    Point_3;

  typedef std::vector<Point_3>                                   Polyline_3;
  typedef std::list<Polyline_3>                                  Polylines;

  typedef FT (Function)(const Point_3&);
  typedef Implicit_mesh_domain_3<Function, K>                    Mesh_domain;
  typedef Mesh_domain_with_polyline_features_3<Mesh_domain>      Mesh_domain_with_features;

  typedef typename Mesh_triangulation_3<Mesh_domain>::type       Tr;
  typedef Mesh_complex_3_in_triangulation_3<Tr>                  C3t3;

  typedef Mesh_criteria_3<Tr>                                    Mesh_criteria;
  typedef typename Mesh_criteria::Edge_criteria                  Edge_criteria;
  typedef typename Mesh_criteria::Facet_criteria                 Facet_criteria;
  typedef typename Mesh_criteria::Cell_criteria                  Cell_criteria;

  // Domain
  Mesh_domain_with_features domain(cube_function<K>,
                                   typename K::Sphere_3(Point_3(1., 1., 1.), 5.*5.),
                                   1e-6);

  // Polylines
  Polylines polylines;
  create_polylines<K>(polylines);
//  domain.add_features(polylines.begin(), polylines.end());

  // Set mesh criteria
  Edge_criteria edge_criteria(0.04);
  Facet_criteria facet_criteria(30, 0.04, 0.04); // angle, size, approximation
  Cell_criteria cell_criteria(2., 0.04); // radius-edge ratio, size
  Mesh_criteria criteria(edge_criteria, facet_criteria, cell_criteria);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

  std::cout << "Number of vertices: " << c3t3.triangulation().number_of_vertices() << std::endl;

  // Output
  std::ofstream medit_file("Canvas.mesh");
  c3t3.output_to_medit(medit_file);
}

namespace Anisotropic_mesh_3
{

enum FMM_state
{
  CHANGED = 0, // used only in the Konukoglu algorithm. 'CHANGED' => 'KNOWN'
  KNOWN,
  TRIAL,
  FAR
};

enum PQ_state
{
  NOTHING_TO_DO = 0,
  REBUILD_TRIAL,
  REBUILD_CHANGED, // used only in the Konukoglu algorithm
  REBUILD_BOTH // used only in the Konukoglu algorithm
};

template<typename Cp>
struct Canvas_point_comparer
{
  bool operator()(Cp const * const cp1, Cp const * const cp2)
  {
    return cp1->distance_to_closest_seed > cp2->distance_to_closest_seed;
  }
};

template<typename K>
class Canvas_point
{
public:
  typedef typename K::FT                                    FT;
  typedef typename K::Point_3                               Point_3;
  typedef Metric_base<K>                                    Metric;
  typedef typename Eigen::Matrix<FT, 3, 1>                  Vector3d;

  Point_3 point;
  std::size_t index;
  FT distance_to_closest_seed;
  std::size_t closest_seed_id;
  FMM_state state;
  Metric metric;
  const Canvas_point* ancestor;

  void change_state(FMM_state new_state,
                    FT& known_count, FT& trial_count, FT& far_count)
  {
    if(new_state == state)
      std::cerr << "WARNING: useless state change..." << std::endl;

    if(state == KNOWN)
      known_count--;
    else if(state == TRIAL)
      trial_count--;
    else if(state == FAR)
      far_count--;

    if(new_state == KNOWN)
      known_count++;
    else if(new_state == TRIAL)
      trial_count++;
    else if(new_state == FAR)
      far_count++;

    state = new_state;
  }

  void initialize_from_point(const FT d,
                             const std::size_t seed_id)
  {
    closest_seed_id = seed_id;
    distance_to_closest_seed = d;
    state = TRIAL;
    ancestor = NULL;
  }

  void print_ancestor_tree() const
  {
    std::cout << "ancestors: ";
    const Canvas_point* anc = ancestor;
    while(anc)
    {
      std::cout << anc->index << " ";
      anc = anc->ancestor;
    }
    std::cout << std::endl;
  }

  Canvas_point() { }

  template<typename MF>
  Canvas_point(const Point_3& point_, const std::size_t index_, const MF mf)
    :
      point(point_),
      index(index_),
      distance_to_closest_seed(FT_inf),
      closest_seed_id(-1),
      state(FAR),
      metric(mf->compute_metric(point)),
      ancestor(NULL)
  { }
};

template<typename K, typename KExact, typename Canvas_point, typename Metric_field>
class Canvas;

template<typename K, typename Canvas>
class Canvas_seeds
{
public:
  typedef typename K::FT                                       FT;
  typedef typename K::Point_3                                  Point_3;

  typedef Stretched_Delaunay_3<K>                              Star;
  typedef typename Star::Metric                                Metric;

  const std::string seeds_str;
  std::size_t max_seeds_n;

  // the metric field and the seeds
  std::vector<Point_3> seeds;
  std::vector<Metric> seeds_metrics;

  const Canvas& canvas;

  const Point_3& operator[](const int i) const { return seeds[i]; }
  Point_3& operator[] (const int i) { return seeds[i]; }
  std::size_t size() const { return seeds.size(); }

  std::size_t insert_new_seed(const FT x, const FT y, const FT z)
  {
#ifdef FILTER_SEEDS_OUTSIDE_CANVAS
    if(canvas.is_point_outside_canvas(x, y, z))
    {
  #if (verbose > 1)
      std::cout << "filtered : " << x << " " << y << " " << z << std::endl;
  #endif
      return seeds.size();
    }
#endif

#if (verbose > 0)
    std::cout << "added new seed: " << x << " " << y << " " << z << std::endl;
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

    in >> word >> useless; //MeshVersionFormatted i
    in >> word >> dim; //Dimension d
    in >> word >> nv;
    std::cout << "seeds nv: " << nv << std::endl;
    CGAL_assertion(dim == 3);

    std::size_t min_nv = (std::min)(nv, max_seeds_n);

    seeds.reserve(min_nv);
    seeds_metrics.reserve(min_nv);

    for(std::size_t i=0; i<nv; ++i)
    {
      in >> r_x >> r_y >> r_z >> useless;
      r_z *= 3.;
//      insert_new_seed(r_x, r_y, r_z);

      insert_new_seed(0., 0., 0.);
      insert_new_seed(1.47, 1.21, 0.47);
      insert_new_seed(0.86, 0.34, 0.22);
      insert_new_seed(0.0152104, 0.740058, 0.12102);
//      insert_new_seed(0.625, 0.898728, 0.775848);
//      insert_new_seed(0.412919, 0.014439, 0.343341);
//      insert_new_seed(0.212919, 0.014839, 0.656659);
//      insert_new_seed(0.488955, 0.73307, 0.944632);
//      insert_new_seed(0.727102, 0.268218, 0.995812);
//      insert_new_seed(0.929122, 0.5, 0.961361);
//      insert_new_seed(0.62922, 0.14, 0.0386388);
//      insert_new_seed(0.127102, 0.268218, 0.00418825);
//      insert_new_seed(0.0552183, 0.749382, 0.919194);
//      insert_new_seed(0.754483, 0.479166, 0.671616);
//      insert_new_seed(0.25, 0.911091, 0.998673);

      if(seeds.size() >= max_seeds_n)
        break;
    }
#if (verbose > 0)
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

template<typename K, typename KExact, typename Cpoint, typename MF>
class Canvas
{
private:
  typedef Canvas<K, KExact, Cpoint, MF>                     Self;

public:
  typedef Canvas_seeds<K, Self>                             Seeds;
  typedef Cpoint                                            Canvas_point;
  typedef MF                                                Metric_field;

  typedef typename K::FT                                    FT;
  typedef typename K::Point_3                               Point_3;

  typedef Stretched_Delaunay_3<K>                           Star;
  typedef typename Star::Metric                             Metric;

  typedef std::set<std::size_t>                             Simplex;
  typedef boost::array<std::size_t, 2>                      Edge;
  typedef boost::array<std::size_t, 3>                      Tri;
  typedef boost::array<std::size_t, 4>                      Tet;

  // tricky bit here: the priority queue uses the BASE canvas point type
  // rather than the derived canvas point type because, when we update distances,
  // we call a function from the base canvas point and look through the neighbors
  // while IN THE BASE CLASS (therefore base canvas points). We want to add those
  // to the priority queue so the PQ has to be pointers to base canvas points.
  // The point of the derived type is just to provide function overloads, so
  // everything works out !
  typedef std::vector<Canvas_point*>                        Canvas_point_vector;
  typedef typename Canvas_point::Neighbors                  Neighbors;
  typedef typename Eigen::Matrix<FT, 3, 1>                  Vector3d;

  // Canvas name
  const std::string canvas_str;

  // The Canvas
  std::vector<Canvas_point> points;
  std::vector<Tet> tetrahedra;
  Canvas_point_vector trial_points;

  // Seeds
  Seeds seeds;

  // Metric field
  const Metric_field mf;

  // Duality
  mutable std::set<Edge> dual_edges;
  mutable std::set<Tri> dual_triangles;
  mutable std::set<Tet> dual_tetrahedra;

  // Refinement
  std::size_t n_refine;
  mutable const Canvas_point* refinement_point;
  mutable FT refinement_distance;

  // Bbox of the canvas
  CGAL::Bbox_3 canvas_bbox;

  //debug & info
  std::size_t known_count, trial_count, far_count;

  // ---------------------------------------------------------------------------

  virtual void initialize() = 0;
  virtual void output_canvas(const std::string str_base) const = 0;
  virtual void compute_dual() const = 0;

  // ---------------------------------------------------------------------------

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

  void clear_dual()
  {
    dual_edges.clear();
    dual_triangles.clear();
    dual_tetrahedra.clear();

    refinement_point = NULL;
    refinement_distance = 0.;
  }

  void initialize_canvas_point(Canvas_point* cp,
                               const FT distance_from_seed,
                               const std::size_t seed_id)
  {
    if(cp->closest_seed_id != static_cast<std::size_t>(-1))
    {
      std::cerr << "WARNING: a new seed is overwriting the closest seed id";
      std::cerr << " of a canvas point! previous index is: "
                << cp->closest_seed_id << std::endl;
    }

    // We can't accept two seeds for one canvas point
    if(cp->state == TRIAL)
      CGAL_assertion(false && "the canvas is not dense enough for the input seeds...");

    cp->initialize_from_point(distance_from_seed, seed_id);

    trial_points.push_back(cp);
    std::push_heap(trial_points.begin(), trial_points.end(),
                   Canvas_point_comparer<Canvas_point>());
  }

  virtual void locate_and_initialize(const Point_3& p,
                                     const std::size_t seed_id)
  {
    // something pretty todo...
    // for now: brutally find the canvas point closest to the seed

    Canvas_point* cp; // this will be the closest seed

    FT min_d = FT_inf;
    for(std::size_t i=0; i<points.size(); ++i)
    {
      Vector3d v;
      v(0) = p.x() - points[i].point.x();
      v(1) = p.y() - points[i].point.y();
      v(2) = p.z() - points[i].point.z();
      const Eigen::Matrix3d& m = points[i].metric.get_mat();

      FT d = v.transpose() * m * v; // note that this is the squared_distance
      if(d < min_d)
      {
        min_d = d;
        cp = &(points[i]);
      }
    }

#if (verbose > 10)
    std::cout << "looking for p: " << p << std::endl;
    std::cout << "found cp: " << cp->index << " [" << cp->point << "] ";
    std::cout << "at distance: " << std::sqrt(min_d) << std::endl;
#endif

    initialize_canvas_point(cp, std::sqrt(min_d), seed_id);
  }

  void compute_canvas_bbox()
  {
    canvas_bbox = CGAL::Bbox_3();
    for(std::size_t i=0; i<points.size(); ++i)
    {
      canvas_bbox += CGAL::Bbox_3(points[i].point.x(), points[i].point.x(),
                                  points[i].point.y(), points[i].point.y(),
                                  points[i].point.z(), points[i].point.z());
    }

    std::cout << "final bbox: x "
              << canvas_bbox.xmin() << " " << canvas_bbox.xmax() << " y "
              << canvas_bbox.ymin() << " " << canvas_bbox.ymax() << " z "
              << canvas_bbox.zmin() << " " << canvas_bbox.zmax() << std::endl;
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

  void dual_shenanigans(const Canvas_point* cp)
  {
    CGAL_assertion(cp->state == KNOWN);

    // Check the neighbors of cp for points that have the state KNOWN.
    // In these, consider those that have a closest_seed_id different than cp's
    // those give dual simplices.

    std::size_t p_id = cp->closest_seed_id;
    std::set<std::size_t> simplex;
    simplex.insert(p_id);

    typename Neighbors::const_iterator it = cp->neighbors.begin(),
                                       iend = cp->neighbors.end();
    for(; it!=iend; ++it)
    {
      const Canvas_point* cq = *it;
      if(!cq || cq->state != KNOWN)
        continue;

      std::size_t q_id = cq->closest_seed_id;
      if(p_id != q_id) // we're on the dual of [PQ]
        simplex.insert(q_id);
    }

    // the simplex of size ? is added to the dual simplices
    add_simplex_to_triangulation(cp, simplex);
  }

  virtual void paint()
  {
    std::clock_t start = std::clock();

#if (verbose > 5)
    std::cout << "main loop" << std::endl;
#endif
    clear_dual();
    Canvas_point* cp;

    bool is_t_empty = trial_points.empty();
    while(!is_t_empty)
    {
      if(known_count%10000 == 0)
        print_states();

#if (verbose > 5)
      std::cout << "Trial queue size : " << trial_points.size() << std::endl;
#endif

#if (verbose > 55)
      std::cout << "trial heap: " << std::endl;
      for(typename std::vector<Canvas_point*>::iterator it = trial_points.begin();
           it != trial_points.end(); ++it)
        std::cout << (*it)->index << " " << (*it)->distance_to_closest_seed << std::endl;
      std::cout << std::endl;
#endif

      cp = trial_points.front();
      CGAL_assertion(cp && cp->state == TRIAL);
      std::pop_heap(trial_points.begin(), trial_points.end(),
                    Canvas_point_comparer<Canvas_point>());
      trial_points.pop_back();

#if (verbose > 10)
      std::cout << "picked n° " << cp->index << " (" << cp->point << ") ";
      std::cout << "at distance : " << cp->distance_to_closest_seed << " from " << cp->closest_seed_id << std::endl;
#endif

      cp->state = KNOWN;
      known_count++; // tmp --> use change_state()

//      dual_shenanigans(cp);
//      boost::unordered_set<Canvas_point*>::const_iterator it = cp->neighbors.begin(),
//                                                        iend = cp->neighbors.end();
//      for(; it!=iend; ++it)
//      {
//        const Canvas_point* cq = *it;
//        if(cq && cq->state == KNOWN)
//          dual_shenanigans(cq);
//      }

      PQ_state pqs = cp->update_neighbors_distances(trial_points);

      if(pqs == REBUILD_TRIAL)
        std::make_heap(trial_points.begin(), trial_points.end(),
                       Canvas_point_comparer<Canvas_point>());

      is_t_empty = trial_points.empty();
    }

#if (verbose > 20)
    std::cout << "final states after painting: " << std::endl;
    for(std::size_t i=0; i<points.size(); ++i)
    {
      const Canvas_point* cp = &(points[i]);
      std::cout << cp->index << "( " << cp->point << ")";
      std::cout << " at distance: " << cp->distance_to_closest_seed << std::endl;
      cp->print_ancestor_tree();
    }
#endif
    std::cerr << "End of paint. time: ";
    std::cerr << ( std::clock() - start ) / (double) CLOCKS_PER_SEC << std::endl;
  }

  virtual void debug()
  {
    CGAL_assertion(trial_points.empty() );

    std::cout << "BRUTE FORCE CHECK THAT WE ARE REALLY FINISHED : " << std::endl;
    for(std::size_t i=0; i<points.size(); ++i)
    {
      Canvas_point* cp = &(points[i]);

      CGAL_assertion(cp->state != FAR);

      std::cout << "point " << i << " min distance is supposedly: ";
      std::cout << cp->distance_to_closest_seed << std::endl;

#if 0
      typename Neighbors::const_iterator it = cp->neighbors.begin(),
                                         iend = cp->neighbors.end();
      for(; it!=iend; ++it)
      {
        const Canvas_point* cq = *it;
        if(cq)
          CGAL_assertion(!cp->compute_closest_seed(cq));
      }
#else
      if(cp->ancestor)
        CGAL_assertion(!cp->compute_closest_seed(static_cast<const Canvas_point*>(cp->ancestor)));
      // other option would be to crtp the base canvas point...
      // also, can't dynamic cast since there's no polymorphism
#endif
    }
  }

  void refresh_canvas_point_states()
  {
    CGAL_assertion(trial_points.empty());
    for(std::size_t i=0; i<points.size(); ++i)
      points[i].state = FAR;

    known_count = 0;
    trial_count = 0;
    far_count = points.size();
  }

  void refine_seeds(const Point_3 new_seed)
  {
#if (verbose > 2)
    std::cout << "refine with : " << new_seed << std::endl;
#endif

    std::size_t old_seed_size =  seeds.size();
    seeds.max_seeds_n = seeds.insert_new_seed(new_seed.x(), new_seed.y(), new_seed.z());
    CGAL_assertion(old_seed_size != seeds.size());
    refresh_canvas_point_states(); // we can't spread from the new seed if all states are 'KNOWN'
    locate_and_initialize(new_seed, seeds.size()-1);
    paint();
  }

  void refine_seeds_with_self_computed_ref_point()
  {
    if(!refinement_point)
    {
      std::cerr << "Couldn't find a ref point, need more initial points" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    refine_seeds(refinement_point->point);
  }

  FT refine()
  {
    if(n_refine < 1)
      return 0;

    std::clock_t start = std::clock();

    for(std::size_t i=0; i<n_refine; ++i)
    {
      std::ostringstream out;
      out << "ref_" << i;
      output_canvas_data_and_dual(out.str()); // this temporarily computes the ref pt
      refine_seeds_with_self_computed_ref_point();
    }

    double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cerr << "End refinement: " << duration << std::endl;
    return duration;
  }

  void add_triangles_edges_to_dual_edges(const Tri& tr) const
  {
    // no need to sort as long as we keep the indices in order
    Edge e;
    e[0] = tr[0]; e[1] = tr[1]; dual_edges.insert(e);
    e[0] = tr[0]; e[1] = tr[2]; dual_edges.insert(e);
    e[0] = tr[1]; e[1] = tr[2]; dual_edges.insert(e);
  }

  void add_tets_triangles_and_edges_to_dual_triangles(const Tet& tet) const
  {
    // no need to sort as long as we keep the indices in order
    Edge e;
    e[0] = tet[0]; e[1] = tet[1]; dual_edges.insert(e);
    e[0] = tet[0]; e[1] = tet[2]; dual_edges.insert(e);
    e[0] = tet[0]; e[1] = tet[3]; dual_edges.insert(e);
    e[0] = tet[1]; e[1] = tet[2]; dual_edges.insert(e);
    e[0] = tet[1]; e[1] = tet[3]; dual_edges.insert(e);
    e[0] = tet[2]; e[1] = tet[3]; dual_edges.insert(e);

    Tri tr;
    tr[0] = tet[0]; tr[1] = tet[1]; tr[2] = tet[2]; dual_triangles.insert(tr);
    tr[0] = tet[0]; tr[1] = tet[1]; tr[2] = tet[3]; dual_triangles.insert(tr);
    tr[0] = tet[0]; tr[1] = tet[2]; tr[2] = tet[3]; dual_triangles.insert(tr);
    tr[0] = tet[1]; tr[1] = tet[2]; tr[2] = tet[3]; dual_triangles.insert(tr);
  }

  template<typename Iterator>
  void add_tets_triangles_and_edges_to_dual_triangles(Iterator it, Iterator end) const
  {
    for(; it!=end; ++it)
    {
      const Tet& tet = *it;
      add_tets_triangles_and_edges_to_dual_triangles(tet);
    }
  }

  void add_simplex_to_triangulation(const std::set<std::size_t>& dual_simplex) const
  {
    typename Simplex::const_iterator it = dual_simplex.begin();
    if(dual_simplex.size() == 2) // an edge!
    {
      Edge e; e[0] = *it; e[1] = (*++it);
      dual_edges.insert(e);
    }
    else if(dual_simplex.size() == 3) // a triangle!
    {
      Tri tr; tr[0] = *it; tr[1] = (*++it); tr[2] = (*++it);
      dual_triangles.insert(tr);

      add_triangles_edges_to_dual_edges(tr);
    }
    else if(dual_simplex.size() == 4)
    {
      Tet tet; tet[0] = *it; tet[1] = (*++it); tet[2] = (*++it); tet[3] = (*++it);
      dual_tetrahedra.insert(tet);

      add_tets_triangles_and_edges_to_dual_triangles(tet);
    }
    else // dual is made of > 4 points (that's a cosphericity)
    {
      std::cerr << "WARNING: COSPHERICITY... " << dual_simplex.size() << std::endl;

      // insert all combinations of tetrahedras... (this will create illegal intersections)
      std::vector<boost::array<std::size_t, 4> > combis =
                         combinations<4>(dual_simplex, dual_simplex.begin(), 4);
      this->dual_tetrahedra.insert(combis.begin(), combis.end());
      add_tets_triangles_and_edges_to_dual_triangles(combis.begin(), combis.end());
    }
  }

  void add_simplex_to_triangulation(const Tet& canvas_tet,
                                    const std::set<std::size_t>& dual_simplex) const
  {
    if(dual_simplex.size() <= 1)
      return;

    // to get the next refinement point (farthest point in a dual...)
    // Not sure if we should allow dual of edges/triangles to be considered...
    // Maybe only allow the dual of an edge/triangle to be considered if the dual is
    // on the border of the domain ?

    typename Tet::const_iterator it = canvas_tet.begin();
    typename Tet::const_iterator iend = canvas_tet.end();
    for(; it!=iend; ++it)
    {
      const Canvas_point& cp = points[*it];
      if(cp.distance_to_closest_seed > refinement_distance &&
         !is_point_outside_canvas(cp.point))
      {
        refinement_distance = cp.distance_to_closest_seed;
        refinement_point = &cp;
      }
    }

    add_simplex_to_triangulation(dual_simplex);
  }

  void add_simplex_to_triangulation(const Canvas_point* cp,
                                    const std::set<std::size_t>& dual_simplex) const
  {
    if(dual_simplex.size() <= 1)
      return;

    if(cp->distance_to_closest_seed > refinement_distance &&
       !is_point_outside_canvas(cp->point))
    {
      refinement_distance = cp->distance_to_closest_seed;
      refinement_point = cp;
    }

    add_simplex_to_triangulation(dual_simplex);
  }

  void output_geodesic_dual(std::ostream&) const
  {
    // todo when the def of a geodesic tet is established...
  }

  void output_straight_dual(const std::string str_base) const
  {
    if(dual_edges.empty() && dual_triangles.empty() && dual_tetrahedra.empty())
      compute_dual();

#if (verbose > 5)
    std::cout << "captured: ";
    std::cout << dual_edges.size() << " edges, ";
    std::cout << dual_triangles.size() << " triangles, ";
    std::cout << dual_tetrahedra.size() << " tetrahedra" << std::endl;
#endif

    std::ofstream out((str_base + "_dual.mesh").c_str());
    out << "MeshVersionFormatted 1" << std::endl;
    out << "Dimension 3" << std::endl;
    out << "Vertices" << std::endl;
    out << seeds.size() << std::endl;
    for(std::size_t i=0; i<seeds.size(); ++i)
      out << seeds[i] << " " << i+1 << std::endl;

    out << "Edges" << std::endl;
    out << dual_edges.size() << std::endl;
    for(typename std::set<Edge>::iterator it = dual_edges.begin();
                                          it != dual_edges.end(); ++it)
    {
      const Edge& edge = *it;
      out << edge[0]+1 << " " << edge[1]+1 << " 1" << std::endl;
    }

    // build from the dual edges, the set of segments that will be used to check
    // self intersections
    std::set<typename K::Segment_3, Segment_3_comparer<K> > segments;
    for(std::set<Edge>::const_iterator it = dual_edges.begin();
                                       it != dual_edges.end(); ++it)
    {
      const Edge& edge = *it;
      const typename K::Segment_3 s(seeds[edge[0]], seeds[edge[1]]);
      segments.insert(s);
    }

    out << "Triangles" << std::endl;
    out << dual_triangles.size() << std::endl;
    for(typename std::set<Tri>::iterator it = dual_triangles.begin();
                                         it != dual_triangles.end(); ++it)
    {
      const Tri& tr = *it;
      out << tr[0]+1 << " " << tr[1]+1 << " " << tr[2]+1 << " ";

      const typename K::Triangle_3 triangle(seeds[tr[0]], seeds[tr[1]], seeds[tr[2]]);
      out << is_triangle_intersected<K,KExact>(triangle, segments) << std::endl;
    }

    out << "Tetrahedra" << std::endl;
    out << dual_tetrahedra.size() << std::endl;
    for(typename std::set<Tet>::iterator it = dual_tetrahedra.begin();
                                         it != dual_tetrahedra.end(); ++it)
    {
      const Tet& tet = *it;
      for(std::size_t i=0; i<tet.size(); ++i)
        out << tet[i] + 1 << " ";
      out << "1" << std::endl;
    }
    out << "End" << std::endl;
  }

  void output_canvas_data_and_dual(const std::string str_base) const
  {
    output_canvas(str_base);
    output_straight_dual(str_base);
  }

  Canvas(const std::string& _canvas_str,
         const std::string& _seeds_str,
         const std::size_t _max_seeds_n,
         const std::size_t _n_refine,
         const Metric_field _mf)
    :
      canvas_str(_canvas_str),
      points(),
      tetrahedra(),
      trial_points(),
      seeds(*this, _seeds_str, _max_seeds_n),
      mf(_mf),
      dual_edges(),
      dual_triangles(),
      dual_tetrahedra(),
      n_refine(_n_refine),
      refinement_point(NULL),
      refinement_distance(0.),
      canvas_bbox(),
      known_count(0),
      trial_count(0),
      far_count(0)
  { }
};


// Case where the canvas is an orthogonal grid
template<typename K, typename KExact, typename Cpoint, typename MF>
class Grid_canvas :
    public Canvas<K, KExact, Cpoint, MF>
{
public:
  typedef Canvas<K, KExact, Cpoint, MF>                         Base;

  typedef Cpoint                                                Canvas_point;
  typedef MF                                                    Metric_field;

  typedef typename K::FT                                        FT;
  typedef typename K::Point_3                                   Point_3;

  typedef typename Base::Simplex                                Simplex;
  typedef typename Base::Edge                                   Edge;
  typedef typename Base::Tri                                    Tri;
  typedef typename Base::Tet                                    Tet;

  typedef typename Base::Vector3d                               Vector3d;

  // Grid geometry
  const Point_3 center;
  const FT side;
  const std::size_t n; // number of points per side

  FT offset_x; // offset is the bottom left point
  FT offset_y;
  FT offset_z;
  const FT step;
  const FT sq_n;

  void locate_and_initialize(const Point_3& p,
                             const std::size_t seed_id)
  {
    int index_x = std::floor((p.x()-offset_x)/step);
    int index_y = std::floor((p.y()-offset_y)/step);
    int index_z = std::floor((p.z()-offset_z)/step);

    Canvas_point* cp = &(this->points[index_z*sq_n + index_y*n + index_x]);

#if (verbose > 10)
    std::cout << "looking for p: " << p << std::endl;
    std::cout << "found cp: " << cp->index << " [" << cp->point << "]" << std::endl;
    std::cout << "index: " << index_x << " " << index_y << " " << index_z << std::endl;
    std::cout << "offset: " << offset_x << " " << offset_y << " " << offset_z << std::endl;
#endif

    Vector3d v;
    v(0) = p.x() - cp->point.x();
    v(1) = p.y() - cp->point.y();
    v(2) = p.z() - cp->point.z();
    const Eigen::Matrix3d& m = cp->metric.get_mat();
    FT d = std::sqrt(v.transpose() * m * v);

    Base::initialize_canvas_point(cp, d, seed_id);
  }

  virtual void initialize()
  {
#if (verbose > 5)
    std::cout << "canvas initialization" << std::endl;
#endif

    // create the canvas points
    for(unsigned int k=0; k<n; ++k)
    {
      for(unsigned int j=0; j<n; ++j)
      {
        for(unsigned int i=0; i<n; ++i)
        {
          // fill from bot left to top right
          Point_3 p(offset_x+i*step, offset_y+j*step, offset_z+k*step);
          Canvas_point cp(p, i + j*n + k*sq_n, this->mf);
          this->points.push_back(cp);
        }
      }
    }

    // assign the neighbors
    for(unsigned int k=0; k<n; ++k)
    {
      for(unsigned int j=0; j<n; ++j)
      {
        for(unsigned int i=0; i<n; ++i)
        {
          std::size_t curr_id = i + j*n + k*sq_n;
          if(k != n-1) // there is a neighbor above
            this->points[curr_id].neighbors[0] = &(this->points[i + j*n + (k+1)*sq_n]);
          if(i != 0) // there is a neighbor left
            this->points[curr_id].neighbors[1] = &(this->points[i-1 + j*n + k*sq_n]);
          if(j != n-1) // there is a neighbor back
            this->points[curr_id].neighbors[2] = &(this->points[i + (j+1)*n + k*sq_n]);
          if(i != n-1) // there is a neighbor right
            this->points[curr_id].neighbors[3] = &(this->points[i+1 + j*n + k*sq_n]);
          if(j != 0) // there is a neighbor front
            this->points[curr_id].neighbors[4] = &(this->points[i + (j-1)*n + k*sq_n]);
          if(k != 0) // there is a neighbor below
            this->points[curr_id].neighbors[5] = &(this->points[i + j*n + (k-1)*sq_n]);
        }
      }
    }
#if (verbose > 5)
    std::cout << "neighbors assigned" << std::endl;
#endif

    Base::compute_canvas_bbox();
    this->seeds.initialize_seeds();
    Base::locate_seeds_on_canvas();

#if (verbose > 5)
    std::cout << "canvas initialized" << std::endl;
#endif
  }

  void add_tet_to_simplices(std::set<Tet>& simplices,
                            const Canvas_point& cp,
                            const Canvas_point* n1,
                            const Canvas_point* n2,
                            const Canvas_point* n3) const
  {
    Tet t;
    if(n1 && n2 && n3)
    {
      t[0] = cp.index;
      t[1] = n1->index;
      t[2] = n2->index;
      t[3] = n3->index;
      std::sort(t.begin(), t.end()); // is that really needed ?
      simplices.insert(t);
    }
  }

  void add_tet_to_simplices(std::set<Tet>& simplices, const Canvas_point& cp,
                            std::size_t i1, std::size_t i2, std::size_t i3) const
  {
    return add_tet_to_simplices(simplices, cp, cp.neighbors[i1], cp.neighbors[i2], cp.neighbors[i3]);
  }

  void output_canvas(const std::string str_base) const
  {
#if (verbose > 0)
    std::cout << "Output canvas" << std::endl;
#endif

    std::ofstream out((str_base + ".mesh").c_str());
    std::ofstream out_bb((str_base + ".bb").c_str());
    out << "MeshVersionFormatted 1" << std::endl;
    out << "Dimension 3" << std::endl;
    out << "Vertices" << std::endl;
    out << this->points.size() << std::endl;
    out_bb << "3 1 " << this->points.size() << " 2" << std::endl;

    std::set<Tet> simplices;
    for(std::size_t i=0; i<this->points.size(); ++i)
    {
      const Canvas_point& cp = this->points[i];
      out << cp.point << " " << i+1 << std::endl;

//      out_bb << cp.closest_seed_id << std::endl;
      out_bb << cp.distance_to_closest_seed << std::endl;

      // compute the tets...
      const boost::array<Canvas_point*, 6>& ns = cp.neighbors; // just to get a shorter name...
      for(std::size_t i=0; i<ns.size(); ++i)
      {
#if 0
        // lazily : the eight corner tetrahedra (creates self intersections,
        // but it doesn't really matter except visually)
        add_tet_to_simplices(simplices, cp, 0, 1, 2);
        add_tet_to_simplices(simplices, cp, 0, 2, 3);
        add_tet_to_simplices(simplices, cp, 0, 3, 4);
        add_tet_to_simplices(simplices, cp, 0, 4, 1);
        add_tet_to_simplices(simplices, cp, 5, 1, 2);
        add_tet_to_simplices(simplices, cp, 5, 2, 3);
        add_tet_to_simplices(simplices, cp, 5, 3, 4);
        add_tet_to_simplices(simplices, cp, 5, 4, 1);
#else
        // a bit smarter : a cube is decomposed in 5 tets and we pick them
        // correctly so we have no overlap
        add_tet_to_simplices(simplices, cp, 0, 2, 3);
        add_tet_to_simplices(simplices, cp, 0, 1, 4);
        add_tet_to_simplices(simplices, cp, 3, 4, 5);
        add_tet_to_simplices(simplices, cp, 1, 2, 5);

        // the last one is a bit annoying since we have to grab neighbors of neighbors
        const Canvas_point* n4 = cp.neighbors[4];
        const Canvas_point* n3 = cp.neighbors[3];
        if(n4 && n3)
          add_tet_to_simplices(simplices, cp, n4->neighbors[3],
                               n4->neighbors[0], n3->neighbors[0]);
#endif
      }
    }

    out << "Tetrahedra" << std::endl;
    out << simplices.size() << std::endl;
    for(typename std::set<Tet>::iterator it = simplices.begin();
                                             it != simplices.end(); ++it)
    {
      const Tet& tet = *it;
      std::set<std::size_t> materials;

      for(std::size_t i=0; i<tet.size(); ++i)
      {
        out << tet[i] + 1 << " ";
        materials.insert(this->points[tet[i]].closest_seed_id);
      }
      std::size_t mat = (materials.size()==1)?(*(materials.begin())):(this->seeds.size());
      out << mat << std::endl;
    }

    std::set<Tri> triangles;
    for(typename std::set<Tet>::iterator it = simplices.begin();
                                         it != simplices.end(); ++it)
    {
      const Tet& tet = *it;
      Tri tr; // could use a permutation to make it look nicer, I suppose...
      tr[0] = tet[0]; tr[1] = tet[1]; tr[2] = tet[2]; triangles.insert(tr);
      tr[0] = tet[0]; tr[1] = tet[1]; tr[2] = tet[3]; triangles.insert(tr);
      tr[0] = tet[0]; tr[1] = tet[2]; tr[2] = tet[3]; triangles.insert(tr);
      tr[0] = tet[1]; tr[1] = tet[2]; tr[2] = tet[3]; triangles.insert(tr);
    }

    out << "Triangles" << std::endl;
    out << triangles.size() << std::endl;
    for(typename std::set<Tri>::iterator it = triangles.begin();
                                         it != triangles.end(); ++it)
    {
      const Tri& tr = *it;
      std::set<std::size_t> materials;
      for(std::size_t i=0; i<tr.size(); ++i)
      {
        out << tr[i] + 1 << " ";
        materials.insert(this->points[tr[i]].closest_seed_id);
      }
      std::size_t mat = (materials.size()==1)?(*(materials.begin())):(this->seeds.size());
      out << mat << std::endl;
    }
  }

  void compute_dual() const
  {
#if (verbose > 5)
    std::cout << "dual computations" << std::endl;
#endif

    for(std::size_t i=0; i<this->points.size(); ++i)
    {
      const Canvas_point& cp = this->points[i];
      Simplex dual_simplex;
      dual_simplex.insert(cp.closest_seed_id);

      typename Canvas_point::Neighbors::const_iterator it = cp.neighbors.begin(),
                                                       iend = cp.neighbors.end();
      for(; it!=iend; ++it)
      {
        const Canvas_point* cq = *it;
        if(!cq)
          continue;

        dual_simplex.insert(cq->closest_seed_id);
      }

      Base::add_simplex_to_triangulation(&cp, dual_simplex);
    }
  }

  Grid_canvas(const std::string& _canvas_str,
              const std::string& _seeds_str,
              const Point_3& _center,
              const FT _side,
              const std::size_t points_per_side,
              const std::size_t _max_seeds_n,
              const std::size_t _n_refine,
              const Metric_field _mf)
    :
      Base(_canvas_str, _seeds_str, _max_seeds_n, _n_refine, _mf),
      center(_center),
      side(_side),
      n(points_per_side),
      offset_x(center.x() - side/2.),
      offset_y(center.y() - side/2.),
      offset_z(center.z() - side/2.),
      step(side / (n-1)),
      sq_n(n*n)
  { }
};

}  // Anisotropic_mesh_3
}  // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_PAINTER_HELPER_3_H

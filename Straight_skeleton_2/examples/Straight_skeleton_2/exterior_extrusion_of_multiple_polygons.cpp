#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/create_weighted_straight_skeleton_from_polygon_with_holes_2.h>
#include <CGAL/extrude_skeleton.h>

#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Random.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <CGAL/draw_polygon_with_holes_2.h>
#include <CGAL/draw_straight_skeleton_2.h>
#include <CGAL/draw_surface_mesh.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <vector>
#include <cassert>

namespace SSi = CGAL::CGAL_SS_i;
namespace SEi = CGAL::Straight_skeleton_extrusion::internal;
namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::FT FT;
typedef K::Point_2 Point_2;
typedef K::Point_3 Point_3;

typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes_2;

using Mesh = CGAL::Surface_mesh<Point_3>;

Polygon_2 generate_random_polygon(const int nv,
                                  CGAL::Random& rnd)
{
  typedef CGAL::Random_points_in_square_2<Point_2> Point_generator;

  Polygon_2 poly;
  CGAL::random_polygon_2(nv, std::back_inserter(poly), Point_generator(0.25, rnd));
  return poly;
}

bool exterior_extrusion_of_multiple_polygons(const Polygon_with_holes_2& pwh,
                                             std::vector<std::vector<FT> >& weights,
                                             const FT height,
                                             Mesh& out)
{
  K k;

  // This is a hacked version of extrude_skeleton.h for this particular context
  const bool verbose = true;

  SEi::Slope slope;
  bool valid_input;
  FT vertical_weight;
  std::tie(slope, valid_input, vertical_weight) = SEi::preprocess_weights<FT>(weights);

  if(!valid_input)
  {
    if(verbose)
      std::cerr << "Error: invalid input weights" << std::endl;
    return false;
  }

  if(slope != SEi::Slope::INWARD && slope != SEi::Slope::VERTICAL)
  {
    std::cerr << "Error: disallowed slope" << std::endl;
    return false;
  }

  // End of preprocessing, start the actual skeleton computation

  // build a soup, to be converted to a mesh afterwards
  std::vector<Point_3> points;
  std::vector<std::vector<std::size_t> > faces;
  points.reserve(2 * pwh.outer_boundary().size()); // just a reasonnable guess
  faces.reserve(2 * pwh.outer_boundary().size() + 2*pwh.number_of_holes());

  //
  using Straight_skeleton_2 = CGAL::Straight_skeleton_2<K>;
  using HDS = typename Straight_skeleton_2::Base;
  using HDS_Halfedge_const_handle = typename HDS::Halfedge_const_handle;

  using Offset_builder_traits = CGAL::Polygon_offset_builder_traits_2<K>;
  using Visitor = SEi::Skeleton_offset_correspondence_builder_visitor<Straight_skeleton_2, Offset_builder_traits, K>;
  using Offset_builder = CGAL::Polygon_offset_builder_2<Straight_skeleton_2,
                                                        Offset_builder_traits,
                                                        Polygon_2,
                                                        Visitor>;

  using Offset_polygons = std::vector<boost::shared_ptr<Polygon_2> >;
  using Offset_polygons_with_holes = std::vector<boost::shared_ptr<Polygon_with_holes_2> >;

  SEi::Extrusion_builder<K> builder(k);

  // bottom
  for(Polygon_with_holes_2::Hole_const_iterator it = pwh.holes_begin(); it!= pwh.holes_end(); ++it)
  {
    Polygon_with_holes_2 hwh(*it); // hole becomes a boundary
    // don't reverse faces because the hole itself is already reversed
    builder.construct_horizontal_faces(hwh, 0 /*height*/, points, faces, true/*invert faces*/);
  }

  std::cout << points.size() << " points; " << faces.size() << " faces [HORIZONTAL BOTTOM]" << std::endl;

  std::unordered_map<HDS_Halfedge_const_handle, Point_2> offset_points;

#ifndef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
# error // I have extracted the code that uses snapping
#endif

  std::map<Point_2, Point_2> snapped_positions;

  auto ss_ptr = CGAL::create_interior_weighted_straight_skeleton_2(
                 SSi::vertices_begin(pwh.outer_boundary()),
                 SSi::vertices_end(pwh.outer_boundary()),
                 pwh.holes_begin(), pwh.holes_end(),
                 std::begin(weights[0]), std::end(weights[0]),
                 std::next(std::begin(weights)), std::end(weights),
                 k);

  if(!ss_ptr)
  {
    std::cerr << "Error: encountered an error during skeleton construction" << std::endl;
    return false;
  }

  CGAL::draw(*ss_ptr);

  Visitor visitor(*ss_ptr, offset_points, vertical_weight, snapped_positions);
  Offset_builder ob(*ss_ptr, Offset_builder_traits(), visitor);

  // top
  Offset_polygons raw_output;
  ob.construct_offset_contours(height, std::back_inserter(raw_output));
  std::cout << raw_output.size() << " raw_output" << std::endl;

  if(raw_output.size() == 0)
  {
    std::cerr << "Error: failed to compute offset, maybe the frame isn't big enough" << std::endl;
    return false;
  }

  // Manually filter the offset of the outer frame
  std::swap(raw_output[0], raw_output.back());
  raw_output.pop_back();

  // convert hole offsets into polygons with holes
  Offset_polygons_with_holes output;
  for(const auto& ro : raw_output)
  {
    Polygon_with_holes_2 pwh(*ro);
    boost::shared_ptr<Polygon_with_holes_2> spwh = boost::make_shared<Polygon_with_holes_2>(pwh);
    output.push_back(spwh);
  }

  builder.construct_horizontal_faces(output, height, points, faces);
  std::cout << points.size() << " points; " << faces.size() << " faces [HORIZONTAL TOP]" << std::endl;

  // lateral
  builder.construct_lateral_faces(*ss_ptr, ob, height, points, faces, offset_points, vertical_weight, snapped_positions, true/*ignore frame faces*/, true/*invert faces*/);
  std::cout << points.size() << " points; " << faces.size() << " faces [LATERAL]" << std::endl;

  // snap
  SEi::apply_snapping<K>(points, snapped_positions);

  PMP::merge_duplicate_points_in_polygon_soup(points, faces);

  CGAL::IO::write_polygon_soup("out_soup.off", points, faces, CGAL::parameters::stream_precision(17));

  if(!PMP::is_polygon_soup_a_polygon_mesh(faces))
    PMP::orient_polygon_soup(points, faces);
  CGAL_assertion(PMP::is_polygon_soup_a_polygon_mesh(faces));

  PMP::polygon_soup_to_polygon_mesh(points, faces, out);

  return true;
}

bool exterior_extrusion_of_multiple_polygons(std::vector<Polygon_2>& polygons,
                                             std::vector<std::vector<FT> >& poly_weights,
                                             const FT height,
                                             Mesh& out)
{
  CGAL_precondition(polygons.size() == poly_weights.size());

  // turn the polygons into holes
  for(std::size_t i=0; i<polygons.size(); ++i)
    polygons[i].reverse_orientation();

  // compute an enclosing outer boundary
  FT lm = 0;

  for(std::size_t i=0; i<polygons.size(); ++i)
  {
    boost::optional<FT> margin = CGAL::compute_outer_frame_margin(SSi::vertices_begin(polygons[i]),
                                                                  SSi::vertices_begin(polygons[i]),
                                                                  poly_weights[i].begin(),
                                                                  poly_weights[i].end(),
                                                                  height);
    if(margin)
      lm = (std::max)(lm, *margin);
  }

  lm = 1000;

  CGAL::Bbox_2 bbox;
  for(const Polygon_2& poly : polygons)
    for(const Point_2& p : poly)
      bbox += p.bbox();

  const FT fxmin = FT(bbox.xmin()) - lm;
  const FT fxmax = FT(bbox.xmax()) + lm;
  const FT fymin = FT(bbox.ymin()) - lm;
  const FT fymax = FT(bbox.ymax()) + lm;

  Point_2 frame[4];
  frame[0] = Point_2(fxmin, fymin);
  frame[1] = Point_2(fxmax, fymin);
  frame[2] = Point_2(fxmax, fymax);
  frame[3] = Point_2(fxmin, fymax);

  // Create the PWH
  Polygon_2 obp(frame, frame + 4);
  Polygon_with_holes_2 pwh(obp);
  for(std::size_t i=0; i<polygons.size(); ++i)
    pwh.add_hole(polygons[i]);

  std::cout << pwh.outer_boundary().size() << " vertices in the outer boundary" << std::endl;
  std::cout << pwh.number_of_holes() << " holes" << std::endl;

  CGAL::draw(pwh);

  // put a weight large enough such that frame edges are not relevant
  FT max_weight = FT(0);
  for(const std::vector<FT>& ws : poly_weights)
  {
    for(const FT& w : ws)
    {
      if(w < FT(0))
      {
        std::cerr << "Error: no negative weights allowed in this program" << std::endl;
        return false;
      }

      max_weight = (std::max)(max_weight, w);
    }
  }

  const FT boundary_weight = 1e7 * max_weight;
  std::cout << "boundary weight = " << boundary_weight << std::endl;

  std::vector<std::vector<FT> > weights;
  weights.push_back(std::vector<FT>(4, boundary_weight));

  for(std::size_t i=0; i<polygons.size(); ++i)
  {
    std::vector<FT> hole_weights = poly_weights[i];
    std::reverse(hole_weights.begin(), hole_weights.end());

    // subtelty:
    // If w[0] pointed to v_0, then when we reverse the polygon, the last polygon is pointing to v_{n-1}
    // but it is the edge v_0 v_{n-1}, which has the weight w_0.
    std::rotate(hole_weights.rbegin(), hole_weights.rbegin()+1, hole_weights.rend());

    weights.push_back(hole_weights);
  }

  CGAL_postcondition(weights.size() == pwh.number_of_holes() + 1);

  return exterior_extrusion_of_multiple_polygons(pwh, weights, height, out);
}

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  const int argc_check = argc - 1;

  int poly_n = 2;
  int poly_nv = 10;
  FT height = FT{(std::numeric_limits<double>::max)()};

  // below is only used for random weight generation
  double min_weight = 1., max_weight = 10.;
  int seed = std::time(nullptr);

  for(int i=1; i<argc; ++i)
  {
    if(!strcmp("-h", argv[i]) || !strcmp("--help", argv[i]) || !strcmp("-?", argv[i]))
    {
      std::cout << "Usage: " << argv[0] << "[options].\n"
        "Options:\n"
        "   -n <value>:  number of polygons.\n"
        "   -v <value>:  number of vertices per polygons.\n"
        "   -t <value>: height. Must be non-zero.\n"
        "   -mw <value>: min value for the random weights (must be positive).\n"
        "   -Mw <value>:  max value for the random weights (must be positive and greater than mw).\n"
        "   -s <value>:  random seed.\n"
                << std::endl;

      return EXIT_FAILURE;
    } else if(!strcmp("-n", argv[i]) && i < argc_check) {
      poly_n = std::stoi(argv[++i]);
    } else if(!strcmp("-v", argv[i]) && i < argc_check) {
      poly_nv = std::stoi(argv[++i]);
    } else if(!strcmp("-t", argv[i]) && i < argc_check) {
      height = std::stod(argv[++i]);
    } else if(!strcmp("-mw", argv[i]) && i < argc_check) {
      min_weight = std::stod(argv[++i]);
    } else if(!strcmp("-Mw", argv[i]) && i < argc_check) {
      max_weight = std::stod(argv[++i]);
    } else if(!strcmp("-s", argv[i]) && i < argc_check) {
      seed = std::stoi(argv[++i]);
    }
  }

  if(height == FT{(std::numeric_limits<double>::max)()} || height <= FT(0))
  {
    std::cerr << "Error: you must specify a positive height" << std::endl;
    return EXIT_FAILURE;
  }

  if(poly_n <= 0 || poly_nv <= 0)
  {
    std::cerr << "Error: invalid (negative) values for #polygons and #vertices per polygon: " << poly_n << " " << poly_nv << std::endl;
    return EXIT_FAILURE;
  }

  if(min_weight >= max_weight)
  {
    std::cerr << "Error: invalid random weight bounds: " << min_weight << " " << max_weight << std::endl;
    return EXIT_FAILURE;
  }

  if(CGAL::is_zero(height))
  {
    std::cerr << "Error: height must be non-zero: " << height << std::endl;
    return EXIT_FAILURE;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Generate random polygons
  //////////////////////////////////////////////////////////////////////////////////////////////////

  std::cout << poly_n << " polygons of size " << poly_nv << std::endl;
  std::cout << "height = " << height << std::endl;

  CGAL::Random rnd(seed);
  std::cout << "seed = " << rnd.get_seed() << std::endl;

#if 0
  // generate random polygons and weights
  std::vector<Polygon_2> polygons;
  std::vector<std::vector<FT> > weights;
  polygons.reserve(poly_n);
  weights.reserve(poly_n);

  for(int pol_i=0; pol_i<poly_n; ++pol_i)
  {
    Polygon_2 rp = generate_random_polygon(poly_nv, rnd);

    // Random translation to avoid overlap of polygons
    std::vector<Point_2> poly_pts = rp.container();
    double nudge = rnd.get_double(-0.25, 0.25);
    for(Point_2& pt : poly_pts)
      pt = Point_2(pt.x() + pol_i + 1, pt.y() + nudge);
    polygons.push_back(Polygon_2(poly_pts.begin(), poly_pts.end()));

    // polygon's weights
    std::vector<FT> poly_weights;
    poly_weights.reserve(poly_nv);
    for(int i=0; i<poly_nv; ++i)
      poly_weights.push_back(rnd.get_double(min_weight, max_weight));
    weights.push_back(poly_weights);
  }
#else
  // Manual example
  std::vector<Polygon_2> polygons;
  std::vector<std::vector<FT> > weights;

  for(int i=0; i<4; ++i)
  {
    Polygon_2 poly;
    poly.push_back(Point_2(i + 3 + i*3, 3));
    poly.push_back(Point_2(i + 6 + i*3, 3));
    poly.push_back(Point_2(i + 6 + i*3, 6));
    poly.push_back(Point_2(i + 3 + i*3, 6));
    polygons.push_back(poly);

    std::vector<FT> poly_weights;
    poly_weights.push_back(rnd.get_double(min_weight, max_weight) * i);
    poly_weights.push_back(rnd.get_double(min_weight, max_weight) * i);
    poly_weights.push_back(rnd.get_double(min_weight, max_weight) * i);
    poly_weights.push_back(rnd.get_double(min_weight, max_weight) * i);
    weights.push_back(poly_weights);
  }
#endif

  Mesh sm;
  bool res = exterior_extrusion_of_multiple_polygons(polygons, weights, height, sm);

  if(!res)
  {
    std::cerr << "Error encountered during the process" << std::endl;
    return EXIT_SUCCESS;
  }

  CGAL::draw(sm);

  CGAL::IO::write_polygon_mesh("out.off", sm, CGAL::parameters::stream_precision(17));

  return EXIT_SUCCESS;
}

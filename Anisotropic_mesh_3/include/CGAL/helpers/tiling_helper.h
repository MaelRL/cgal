#ifndef TILING_HELPER_H
#define TILING_HELPER_H

// metric computation part from Tiling (cf IMR 12, Pierre Alliez)

#include <algorithm>
#include <vector>
#include <queue>

#include <CGAL/helpers/tiling_include/Largest_empty_iso_rectangle_2.h>
#include <CGAL/helpers/tiling_include/Pfacet.h>
#include <CGAL/helpers/tiling_include/algorithms.h>

#include <CGAL/bounding_box.h>

template <class Polyhedron>
class Tiling
{
  typedef typename Polyhedron::Traits K;
  typedef typename K::Iso_rectangle_2                     Iso_rectangle_2;
  typedef typename K::Point_2                             Point_2;
  typedef typename K::Point_3                             Point_3;
  typedef typename K::Triangle_2                          Triangle_2;
  typedef typename K::Triangle_3                          Triangle_3;
  typedef typename K::Vector_2                            Vector_2;
  typedef typename K::Vector_3                            Vector_3;
  typedef typename K::Segment_2                           Segment_2;
  typedef typename K::Line_2                              Line_2;
  typedef typename K::RT                                  RT;
  typedef typename K::FT                                  FT;

  typedef typename Polyhedron::Facet_handle                             Facet_handle;
  typedef typename Polyhedron::Halfedge_handle                          Halfedge_handle;
  typedef typename Polyhedron::Vertex_handle                            Vertex_handle;
  typedef typename Polyhedron::Halfedge_around_facet_circulator         HF_circulator;
  typedef typename Polyhedron::Halfedge_around_vertex_circulator        HV_circulator;
  typedef typename Polyhedron::Facet_iterator                           Facet_iterator;
  typedef typename Polyhedron::Edge_iterator                            Edge_iterator;
  typedef typename Polyhedron::Halfedge_iterator                        Halfedge_iterator;
  typedef typename Polyhedron::Vertex_iterator                          Vertex_iterator;

  typedef typename CGAL::Exact_predicates_exact_constructions_kernel    kernel_new;
  typedef typename CGAL::Plane_3<kernel_new>                            Plane_3_new;
  typedef typename CGAL::Point_3<kernel_new>                            Point_3_new;
  typedef typename CGAL::Point_2<kernel_new>                            Point_2_new;
  typedef typename CGAL::Vector_3<kernel_new>                           Vector_3_new;
  typedef typename CGAL::Largest_empty_iso_rectangle_2<K>               Largest_empty_iso_rect_2;
  typedef typename CGAL::Aff_transformation_2<K>                        Transformation;

  typedef CPfacet<Facet_handle> Pfacet;
  typedef typename std::priority_queue<Pfacet, std::vector<Pfacet>, more_pfacet<Pfacet> > PQueue_facet;

//functions related to computing the coverage region
  void tolerance_domain(Vertex_handle& start_vh,
                        std::set<Facet_handle>& facets,
                        const double tolerance)
  {
    const Vector_3 ref_vector = start_vh->normal();
    double cos_tolerance = std::cos ( tolerance );
    PQueue_facet queue;

    HV_circulator h = start_vh->vertex_begin();
    HV_circulator hend = h;
    
    do
    {
      Facet_handle fh = h->facet();
      Pfacet pf ( fh, 0.0 );
      queue.push(pf);
      h++;
    }
    while(h != hend);

    while ( !queue.empty() )
    {
      Pfacet pf = queue.top();
      queue.pop();

      Facet_handle f = pf.facet();

      if ( facets.find ( f ) != facets.end() )
        continue;

      facets.insert ( f );

      Halfedge_handle h = f->facet_begin();
      Halfedge_handle seed = h;

      do
      {
        if ( h->is_border_edge() )
        {
          h = h->next();
          continue;
        }

        Facet_handle nf = h->opposite()->facet();

        if ( facets.find ( nf ) != facets.end() )
        {
          h = h->next();
          continue;
        }

        Triangle_3 tr(nf->halfedge()->vertex()->point(),
                      nf->halfedge()->next()->vertex()->point(),
                      nf->halfedge()->next()->next()->vertex()->point());
        Vector_3 nfn = tr.supporting_plane().orthogonal_vector();
        nfn = nfn / CGAL::sqrt(nfn*nfn);

        double cos_angle = nfn * ref_vector;
        if ( cos_angle < cos_tolerance )
        {
          h = h->next();
          continue;
        }

        Pfacet pf2 ( nf, cos_angle );
        queue.push ( pf2 );

        h = h->next();
      }
      while ( h != seed );
    }
  }

//functions related to largest fitting geometry

  //rectangle
  FT largest_inscribed_rectangle ( std::list<Point_2>& rectangle,
                                   const std::list<Point_2>& points,
                                   const Point_2& query,
                                   const std::list<Point_2>& points_for_bbox)
  {
    Iso_rectangle_2 box = CGAL::bounding_box ( points_for_bbox.begin(), points_for_bbox.end() );

    Largest_empty_iso_rect_2 leir ( box );

    typename std::list<Point_2>::const_iterator pit;
    for ( pit = points.begin();
          pit != points.end();
          ++pit )
    {
      const Point_2 p2d = *pit;
      leir.insert ( p2d );
    }

    Iso_rectangle_2 b = leir.get_largest_empty_iso_rectangle ( query );

    //rectangle.clear();

    for ( int i = 0; i != 4 ; ++i )
      rectangle.push_back ( b.vertex ( i ) );

    return b.area();
  }

  void best_rectangle_axis_aligned ( std::list<Point_2>& best_rectangle_2d,
                                     const Point_2& query,
                                     const std::list<Point_2>& points,
                                     const std::list<Segment_2>& segments,
                                     const std::list<Point_2>& points_box)
  {
    FT max_area = largest_inscribed_rectangle ( best_rectangle_2d, points, query, points_box );

    Iso_rectangle_2 box = CGAL::bounding_box ( points_box.begin(), points_box.end() );

    for ( double angle_degree = 1.5; angle_degree < 90.0; angle_degree += 3.0 )
    {
      double angle_radian = degree_to_radian<FT> ( angle_degree );
      Transformation rotate ( CGAL::ROTATION, sin ( angle_radian ), cos ( angle_radian ) );
      Transformation translate ( CGAL::TRANSLATION, Vector_2 ( -query.x(), -query.y() ) );

      std::list<Point_2> new_points;
      rotate_points ( new_points, points, translate, rotate );
      std::list<Point_2> pts_rect;

      std::list<Segment_2> new_segments;
      rotate_segments(new_segments, segments, translate, rotate);

      // compute cardinal points to restrict bbox
      std::list<Point_2> new_points_bbox;
      compute_cardinal_points ( new_points_bbox, new_segments, query);

      FT area = largest_inscribed_rectangle ( pts_rect, new_points, query, new_points_bbox);

      if ( area > max_area )
      {
        max_area = area;
        Transformation inverse_rotate = rotate.inverse();
        best_rectangle_2d.clear();
        rotate_points ( best_rectangle_2d, pts_rect, translate, inverse_rotate );
      }
    }
  }

  bool best_rotated_rectangle ( const Vertex_handle& vh,
                                std::list<Point_3>& best_rectangle,
                                const std::set<Facet_handle>& facets)
  {
    const Point_3 pquery = vh->point();

    // compute mean normal
    Vector_3 normal = compute_mean_normal_facets<Facet_handle>(facets.begin(), facets.end()); //vh->normal();

    Plane_3_new plane_new(Point_3_new(0.0, 0.0, 0.0), Vector_3_new(normal.x(), normal.y(), normal.z()));

    // project boundary in a 2d plane
    std::list<Point_2> points_2d;
    std::list<Segment_2> segments;
    typename std::set<Facet_handle>::const_iterator fit;
    for (fit = facets.begin();
         fit != facets.end();
         ++fit)
    {
      Facet_handle f = *fit;
      Halfedge_handle h = f->facet_begin();
      Halfedge_handle seed = h;
      do
      {
        if (h->is_border_edge() || is_boundary_of_tolerance_domain ( h, facets ) )
        {
          const Point_3 ptarget = h->vertex()->point();
          Point_3_new ptarget_new(ptarget.x(), ptarget.y(), ptarget.z());

          Point_2_new p2d_target = plane_new.to_2d ( ptarget_new );
          points_2d.push_back ( Point_2(CGAL::to_double(p2d_target.x()), CGAL::to_double(p2d_target.y())) );

          const Point_3 psource = h->opposite()->vertex()->point();
          Point_3_new psource_new(psource.x(), psource.y(), psource.z());

          Point_2_new p2d_source = plane_new.to_2d ( psource_new );
          points_2d.push_back ( Point_2(CGAL::to_double(p2d_source.x()), CGAL::to_double(p2d_source.y())) );

          Segment_2 seg(Point_2(CGAL::to_double(p2d_source.x()), CGAL::to_double(p2d_source.y())),
                        Point_2(CGAL::to_double(p2d_target.x()), CGAL::to_double(p2d_target.y())));

          segments.push_back(seg);
        }

        h = h->next();

      } while (h != seed);
    }

    // find 2d analogous of query
    Point_3_new pquery_new(pquery.x(), pquery.y(), pquery.z());
    Point_2_new query_new = plane_new.to_2d ( pquery_new );

    Point_2 query( CGAL::to_double(query_new.x()), CGAL::to_double(query_new.y()));

    // compute cardinal points to restrict bbox
    std::list<Point_2> points_bbox;
    compute_cardinal_points ( points_bbox, segments, query);

    // compute best fit rectangle
    std::list<Point_2> best_rectangle_2d;
    best_rectangle_axis_aligned ( best_rectangle_2d, query, points_2d, segments, points_bbox);

    Point_3_new new_query_point = plane_new.to_3d ( query_new );
    Point_3 nquery(CGAL::to_double(new_query_point.x()), CGAL::to_double(new_query_point.y()), CGAL::to_double(new_query_point.z()));


    // fill best rectangle
    double x = 0.0;
    double y = 0.0;
    typename std::list<Point_2>::iterator pit;
    for ( pit = best_rectangle_2d.begin();
          pit != best_rectangle_2d.end();
          ++pit )
    {
      Point_2 p_2d = *pit;
      Point_2_new p_2d_new(p_2d.x(), p_2d.y());
      x += 0.25 * p_2d.x();
      y += 0.25 * p_2d.y();
      Point_3_new p_new = plane_new.to_3d ( p_2d_new );
      Point_3 p(CGAL::to_double(p_new.x()) + (pquery.x() - nquery.x()),
              CGAL::to_double(p_new.y()) + (pquery.y() - nquery.y()),
              CGAL::to_double(p_new.z()) + (pquery.z() - nquery.z()) );
      best_rectangle.push_back ( p );
    }

    return true;
  }

  void compute_rectangle(const Vertex_handle& vh,
                         const std::set<Facet_handle>& facets)
  {
    std::cout << "compute rectangle" << std::endl;

    std::list<Point_3> best_rectangle;
    best_rotated_rectangle ( vh, best_rectangle, facets);

    // compute metric
    vh->tile_points().insert(vh->tile_points().begin(), best_rectangle.begin(), best_rectangle.end());

    std::vector<Point_3> corners;
    corners.insert(corners.begin(), best_rectangle.begin(), best_rectangle.end());
    Vector_3 vec1 = ( corners[1] - corners[0] ) + ( corners[2] - corners[3] );
    Vector_3 vec2 = ( corners[2] - corners[1] ) + ( corners[3] - corners[0] );

    FT size_1 = 0.5 * std::sqrt ( vec1.squared_length() );
    FT size_2 = 0.5 * std::sqrt ( vec2.squared_length() );
    vh->aniso_ratio() = (std::max)(size_1/size_2, size_2/size_1);

    /*
    // save
    vh->vec1() = 0.5 * vec1;
    vh->vec2() = 0.5 * vec2;
    vh->base_1() = vec1 / std::sqrt(vec1.squared_length () );
    vh->base_2() = vec2 / std::sqrt(vec2.squared_length () );
    */
  }

  //triangle
    // TODO: largest inscribed triangle
  void best_triangle_2d ( std::list<Point_2>& best_tri_2d,
                          const Point_2& query,
                          const std::list<Point_2>& points)
  {  
    typename std::list<Point_2>::const_iterator pit1, pit2, pit3;
    double max_area = 0.0;
    for (pit1 = points.begin(); pit1 != points.end(); ++pit1)
    {
      Point_2 p1 = *pit1;
      for (pit2 = points.begin(); pit2 != points.end(); ++pit2)
      {
        if (pit2 == pit1)
          continue;

        Point_2 p2 = *pit2;
        for (pit3 = points.begin(); pit3 != points.end(); ++pit3)
        {
          if (pit3 == pit2 || pit3 == pit1)
            continue;

          Point_2 p3 = *pit3;
          Triangle_2 tri(p1, p2, p3);

          // if the triangle has negative orientation
          double area = tri.area();
          if (area < 0.0)
            continue;

          if (is_good_candidate_triangle(tri, query, points))
          {
            max_area = (area > max_area) ? area : max_area;
            best_tri_2d.clear();
            best_tri_2d.push_back(tri.vertex(0));
            best_tri_2d.push_back(tri.vertex(1));
            best_tri_2d.push_back(tri.vertex(2));
          }
        }
      }
    }
  }

  bool best_triangle ( const Vertex_handle& vh,
                       std::list<Point_3>& best_tri,
                       const std::set<Facet_handle>& facets)
  {
    const Point_3 pquery = vh->point();

    // compute mean normal
    Vector_3 normal = compute_mean_normal_facets<Facet_handle>(facets.begin(), facets.end()); //vh->normal();

    Plane_3_new plane_new(Point_3_new(0.0, 0.0, 0.0), Vector_3_new(normal.x(), normal.y(), normal.z()));

    // project boundary in a 2d plane
    std::list<Point_2> points_2d;
    std::list<Segment_2> segments;
    typename std::set<Facet_handle>::const_iterator fit;
    for (fit = facets.begin();
         fit != facets.end();
         ++fit)
    {
      Facet_handle f = *fit;
      Halfedge_handle h = f->facet_begin();
      Halfedge_handle seed = h;
      do
      {
        if (h->is_border_edge() || is_boundary_of_tolerance_domain ( h, facets ) )
        {
          const Point_3 ptarget = h->vertex()->point();
          Point_3_new ptarget_new(ptarget.x(), ptarget.y(), ptarget.z());

          Point_2_new p2d_target = plane_new.to_2d ( ptarget_new );
          points_2d.push_back ( Point_2(CGAL::to_double(p2d_target.x()), CGAL::to_double(p2d_target.y())) );

          const Point_3 psource = h->opposite()->vertex()->point();
          Point_3_new psource_new(psource.x(), psource.y(), psource.z());

          Point_2_new p2d_source = plane_new.to_2d ( psource_new );
          points_2d.push_back ( Point_2(CGAL::to_double(p2d_source.x()), CGAL::to_double(p2d_source.y())) );

          Segment_2 seg(Point_2(CGAL::to_double(p2d_source.x()), CGAL::to_double(p2d_source.y())),
                        Point_2(CGAL::to_double(p2d_target.x()), CGAL::to_double(p2d_target.y())));

          segments.push_back(seg);
        }

        h = h->next();

      } while (h != seed);
    }

    // find 2d analogous of query
    Point_3_new pquery_new(pquery.x(), pquery.y(), pquery.z());
    Point_2_new query_new = plane_new.to_2d ( pquery_new );

    Point_2 query( CGAL::to_double(query_new.x()), CGAL::to_double(query_new.y()));

    // compute best fit rectangle
    std::list<Point_2> best_tri_2d;
    best_triangle_2d ( best_tri_2d, query, points_2d);

    Point_3_new new_query_point = plane_new.to_3d ( query_new );
    Point_3 nquery(CGAL::to_double(new_query_point.x()), CGAL::to_double(new_query_point.y()), CGAL::to_double(new_query_point.z()));


    // fill best triangle
    double x = 0.0;
    double y = 0.0;
    typename std::list<Point_2>::iterator pit;
    for ( pit = best_tri_2d.begin();
          pit != best_tri_2d.end();
          ++pit )
    {
      Point_2 p_2d = *pit;
      Point_2_new p_2d_new(p_2d.x(), p_2d.y());
      x += 0.25 * p_2d.x();
      y += 0.25 * p_2d.y();
      Point_3_new p_new = plane_new.to_3d ( p_2d_new );
      Point_3 p(CGAL::to_double(p_new.x()) + (pquery.x() - nquery.x()),
              CGAL::to_double(p_new.y()) + (pquery.y() - nquery.y()),
              CGAL::to_double(p_new.z()) + (pquery.z() - nquery.z()) );
      best_tri.push_back ( p );
    }
    return true;
  }

    // TODO: largest inscribed triangle
  void compute_triangle(const Vertex_handle& vh,
                        const std::set<Facet_handle>& facets)
  {
    // compute rectangle
    std::list<Point_3> best_tri;
    best_triangle ( vh, best_tri, facets);

    //m_pts_triangle = best_tri;

    //m_facet->tri() = m_pts_triangle;
  }

    // TODO: etre moins restrictif
  bool is_good_candidate_triangle( Triangle_2 tri_candidate,
                                   const Point_2& query,
                                   const std::list<Point_2>& points)
  {
    if (!tri_candidate.has_on_positive_side(query))
      return false;

    typename std::list<Point_2>::const_iterator pit;
    for (pit = points.begin(); pit != points.end(); ++pit)
    {
      Point_2 p = *pit;
      if(p == tri_candidate.vertex(0) || p == tri_candidate.vertex(1) || p == tri_candidate.vertex(2))
        continue;

      if (tri_candidate.has_on_positive_side(p))
        return false;
    }
    return true;
  }

  //ellipsoid
  void compute_ellipse(const Vertex_handle& vh,
                       const std::set<Facet_handle>& facets)
  {
    // gather boundary vertices
    std::set<Vertex_handle> bvertices;
    compute_boundary_vertices<Polyhedron, Facet_handle, Vertex_handle>(facets, bvertices);

    std::cout << "# boundary vertices: " << bvertices.size() << std::endl;

    // project query in 2d


    // project boundary in 2d


    // test for all rotations

    //scale
  }

  //rhombus


//level 0
public:
  void compute(Vertex_handle& vh,
               const double tolerance_angle)
  {
    //get geometry
    int geometry = -1;//vh->classification();

    std::cout << "entered compute with : " << vh->point();
    std::cout << " geo : " << geometry << " normal : " << vh->normal() << std::endl;

    // gather facets in tolerance
    std::set<Facet_handle> facets;
    tolerance_domain(vh, facets, tolerance_angle);

    std::cout << facets.size() << " facets." << std::endl;

    //compute best fitting for the given geometry
    switch(geometry)
    {
      case 1:
        compute_rectangle(vh, facets);
      case 2:
        compute_triangle(vh, facets);
      default:
        compute_rectangle(vh, facets);
      //compute_ellipse(vh, facets);
      //compute_rhombus(vh, facets);
    }
  }

  void compute( Polyhedron & polyhedron,
                const double tolerance_angle,
                const bool print )
  {

  }


//functions related to drawing
  void gl_render()
  {

  }

private:
//misc
  bool is_boundary_of_tolerance_domain ( const Halfedge_handle& h,
                                         const std::set<Facet_handle>& facets)
  {
    if (h->is_border_edge())
      return true;

    Facet_handle f = h->facet();
    Facet_handle nf = h->opposite()->facet();

    return ( (facets.find(f) != facets.end() && facets.find(nf) == facets.end()) );
  }

  void compute_cardinal_points ( std::list<Point_2>& points,
                                 const std::list<Segment_2>& segments,
                                 const Point_2 query )
  {

    Point_2 point_1a, point_1b, point_2a, point_2b;
    double sqd_length_1a_min = 1e308;
    double sqd_length_1b_min = 1e308;
    double sqd_length_2a_min = 1e308;
    double sqd_length_2b_min = 1e308;

    Line_2 line1(query, Vector_2(1.0, 0.0));
    Line_2 line2(query, Vector_2(0.0, 1.0));

    typename std::list<Segment_2>::const_iterator sit;
    for (sit = segments.begin();
         sit != segments.end();
         ++sit)
    {
      Segment_2 seg = *sit;

      CGAL::Object result1;
      Point_2 ipoint1;
      Segment_2 iseg1;
      result1 = CGAL::intersection(seg, line1);
      if (CGAL::assign(ipoint1, result1)) { // handle the point intersection case.
        FT sqd_dist = CGAL::squared_distance(query, ipoint1);
        if (ipoint1.x() > query.x())
        {
          if (sqd_dist < sqd_length_1a_min)
          {
            sqd_length_1a_min = sqd_dist;
            point_1a = ipoint1;
          }
        }
        else
        {
          if (sqd_dist < sqd_length_1b_min)
          {
            sqd_length_1b_min = sqd_dist;
            point_1b = ipoint1;
          }
        }
      }
      else if (CGAL::assign(iseg1, result1)) 
      {// handle the segment intersection case.
      } else { // handle the no intersection case.
      }

      CGAL::Object result2;
      Point_2 ipoint2;
      Segment_2 iseg2;
      result2 = CGAL::intersection(seg, line2);
      if (CGAL::assign(ipoint2, result2)) { // handle the point intersection case.
        FT sqd_dist = CGAL::squared_distance(query, ipoint2);
        if (ipoint2.y() > query.y())
        {
          if (sqd_dist < sqd_length_2a_min)
          {
            sqd_length_2a_min = sqd_dist;
            point_2a = ipoint2;
          }
        }
        else
        {
          if (sqd_dist < sqd_length_2b_min)
          {
            sqd_length_2b_min = sqd_dist;
            point_2b = ipoint2;
          }
        }
      }
      else if (CGAL::assign(iseg2, result2)) 
      { // handle the segment intersection case.
      } else
      { // handle the no intersection case.
      }
    }
    points.push_back(point_1a);
    points.push_back(point_2a);
    points.push_back(point_1b);
    points.push_back(point_2b);

  }

  void rotate_points ( std::list<Point_2>& new_points, 
                       const std::list<Point_2>& old_points,
                       Transformation translate, 
                       Transformation rotate )
  {
    Transformation inverse_translate = translate.inverse();

    typename std::list<Point_2>::const_iterator pit;
    for ( pit = old_points.begin();
          pit != old_points.end();
          ++pit )
    {
      Point_2 p = *pit;
      Point_2 q = translate ( p );
      q = rotate ( q );
      q = inverse_translate ( q );
      new_points.push_back ( q );
    }

  }

  void rotate_segments ( std::list<Segment_2>& new_segments, 
                         const std::list<Segment_2>& old_segments, 
                         Transformation translate, 
                         Transformation rotate )
  {
    Transformation inverse_translate = translate.inverse();

    typename std::list<Segment_2>::const_iterator sit;
    for ( sit = old_segments.begin();
          sit != old_segments.end();
          ++sit )
    {
      Segment_2 seg = *sit;
      Point_2 pa = seg.vertex(0);
      Point_2 pb = seg.vertex(1);

      Point_2 qa = translate ( pa );
      qa = rotate ( qa );
      qa = inverse_translate ( qa );

      Point_2 qb = translate ( pb );
      qb = rotate ( qb );
      qb = inverse_translate ( qb );

      Segment_2 new_seg(qa, qb);

      new_segments.push_back ( new_seg );
    }

  }
};

#endif

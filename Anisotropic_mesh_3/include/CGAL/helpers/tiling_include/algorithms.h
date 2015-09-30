#ifndef _ALGORITHMS_
#define _ALGORITHMS_

#include <CGAL/basic.h>

#include <vector>
#include <set>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

template <typename T>
void shift_rotate_to_begin(std::vector<T>& vect, const T& value)
{
  while (vect.front() != value)
  {
    T v = vect.front();
    vect.erase(vect.begin());
    vect.push_back(v);
  }
}

template <typename T>
void compute_intersection(std::vector<T>& inter,
                          const std::vector<T>& vector_1,
                          const std::vector<T>& vector_2)
{
  typename std::vector<T>::const_iterator it_find;
  typename std::vector<T>::const_iterator it;
  for (it = vector_1.begin();
       it != vector_1.end();
       ++it)
  {
    T v = *it;
    it_find = std::find(vector_2.begin(), vector_2.end(), v);
    if (it_find != vector_2.end())
      inter.push_back(v);
  }
}

template <typename T>
void compute_intersection_sets (std::set<T>& inter,
                                const std::set<T>& set_1,
                                const std::set<T>& set_2)
{
  typename std::set<T>::const_iterator it_find;
  typename std::set<T>::const_iterator it;
  for (it = set_1.begin();
       it != set_1.end();
       ++it)
  {
    T v = *it;
    it_find = std::find(set_2.begin(), set_2.end(), v);
    if (it_find != set_2.end())
      inter.insert(v);
  }
}

template <typename T>
T degree_to_radian(const T degree)
{
  return (degree * CGAL_PI) / 180.0;
}

template <typename T>
T radian_to_degree(const T radian)
{
  return (180.0 * radian ) / CGAL_PI;
}




template <typename T>
void remove_common_elements(std::set<T>& result,
                            const std::set<T>& set_1,
                            const std::set<T>& set_2)
{
  result.insert(set_1.begin(), set_1.end());

  typename std::set<T>::const_iterator it_find;

  typename std::set<T>::const_iterator it;
  for (it = set_2.begin();
       it != set_2.end();
       ++it)
  {
    T v = *it;
    it_find = std::find(result.begin(), result.end(), v);
    if (it_find == result.end())
      result.insert(v);
    else
      result.erase(it_find);
  }
}

template <class T, class InputIterator>
Vector compute_mean_normal_facets ( InputIterator begin,
                                    InputIterator end )
{
  Vector mean = CGAL::NULL_VECTOR;
  double sum_areas = 0.0;
  for ( InputIterator fit = begin; fit != end; ++fit )
  {
    T f = *fit;

    Triangle tr(f->halfedge()->vertex()->point(),
                  f->halfedge()->next()->vertex()->point(),
                  f->halfedge()->next()->next()->vertex()->point());
    Vector fn = tr.supporting_plane().orthogonal_vector();

    FT area = std::sqrt(tr.squared_area());
    sum_areas += area;
    mean = mean + area * fn;
  }

  mean = mean / sum_areas;
  mean =  mean / std::sqrt ( mean.squared_length() );
  return mean;
}

// FIXME: value for simple cases?
template <typename F, typename V, typename M>
double conformity_energy(const V& vquery,
                         const std::set<F>& set_facets)
{
  double sum = 0.0;

  if (set_facets.empty())
    return 1.0;
  else if (set_facets.size() == 1)
    return 1.0;
  else
  {
    std::vector<F> facets;
    facets.insert(facets.end(), set_facets.begin(), set_facets.end());

    unsigned int nb = 0;
    for (unsigned int i = 0 ; i != facets.size(); ++i)
    {
      F f1 = facets[i];
      for (unsigned int j = 1; j != facets.size(); ++j)
      {
        if (i == j)
          continue;

        F f2 = facets[j];

        M metric_1(f1);
        M metric_2(f2);

        std::set<V> corners1, corners2;

        V corner1 = metric_1.nearest_corner(/*p, */vquery);
        V corner2 = metric_2.nearest_corner(/*p, */vquery);

        FT d1a = std::fabs ( ( vquery->point() - corner1->point() ) * metric_1.base_1() ) / metric_1.size_1();
        FT d1b = std::fabs ( ( vquery->point() - corner1->point() ) * metric_1.base_2() ) / metric_1.size_2();
        FT d2a = std::fabs ( ( vquery->point() - corner2->point() ) * metric_2.base_1() ) / metric_2.size_1();
        FT d2b = std::fabs ( ( vquery->point() - corner2->point() ) * metric_2.base_2() ) / metric_2.size_2();


        double value1 = (d1a > 1.0 || d1b > 1.0) ? 0.0 : ( (1.0 - d1a) * (1.0 - d1b) );
        double value2 = (d2a > 1.0 || d2b > 1.0) ? 0.0 : ( (1.0 - d1b) * (1.0 - d2b) );

        double sqd_dist = CGAL::squared_distance(corner1->point(), corner2->point());

        double mean1 = 2.0 * (f1->size_1() + f1->size_2());
        double mean2 = 2.0 * (f2->size_1() + f2->size_2());

        double partial_energy = std::sqrt(sqd_dist) * ( value1 / mean1 + value2 / mean2);
        //double partial_energy = 2.0 * std::sqrt(sqd_dist) / (mean1 + mean2);
        sum += partial_energy;
        ++nb;

      }
    }
    sum /= (double)nb;

    //return 1.0 - sum;
    return sum;
  }
}

template <class T>
T only_index_in_first_but_not_in_second ( const std::vector<T>& vect1, const std::vector<T>& vect2 )
{
  typename std::vector<T>::const_iterator it;
  for ( it = vect1.begin();
        it != vect1.end();
        ++it )
  {
    T index = *it;
    if ( std::find ( vect2.begin(), vect2.end(), index ) == vect2.end() )
      return index;
  }
  return -1;
}

template <class T, class P2>
void project_points(const Vector& normal,
                    const Point p,
                    std::vector<Point>& corners,
                    std::vector<P2>& projected_points)
{
  typedef typename CGAL::Plane_3<T> Plane_3_T;
  typedef typename CGAL::Point_3<T> Point_3_T;
  typedef typename CGAL::Point_2<T> Point_2_T;
  typedef typename CGAL::Vector_3<T> Vector_3_T;

  Plane_3_T plane_new(Point_3_T(p.x(), p.y(), p.z()), Vector_3_T(normal.x(), normal.y(), normal.z()));

  for (unsigned int i = 0; i != corners.size(); ++i)
  {
    Point p = corners[i];
    Point_3_T p3d(p.x(), p.y(), p.z());
    Point_2_T p2d = plane_new.to_2d ( p3d );
    projected_points.push_back(P2(CGAL::to_double(p2d.x()), CGAL::to_double(p2d.y())));
  }
}

template <typename kernel1, typename kernel_exact, typename F, typename V>
double estimated_area(F& f)
{
  typedef typename CGAL::Plane_3<kernel_exact> Plane_exact;

  typedef typename CGAL::Point_2<kernel_exact> Point_2;
  typedef typename CGAL::Point_3<kernel_exact> Point_3;
  typedef typename CGAL::Vector_3<kernel_exact> Vector_3;
  typedef CGAL::Polygon_2<kernel_exact> Polygon_2;

  std::vector<V> vertices;
  f->get_vertices(vertices);

  Point_3 pref(vertices[0]->point().x(), vertices[0]->point().y(), vertices[0]->point().z());
  const Vector_3 nref(f->normal().x(), f->normal().y(), f->normal().z());

  Plane_exact plane(pref, nref);

  Polygon_2 pgn;

  for (unsigned int i = 0; i != vertices.size(); ++i) {
    Point_3 p3d(vertices[i]->point().x(), vertices[i]->point().y(), vertices[i]->point().z());
    Point_2 p2d = plane.to_2d( p3d );
    pgn.push_back(p2d);
  }

  return CGAL::to_double(pgn.area());
}

template <class Polyhedron, class PointOutputIterator, class IndexOutputIterator>
void generate( Polyhedron& polyhedron, PointOutputIterator pout, IndexOutputIterator iout)
{
  typedef typename Polyhedron::Vertex_handle Vertex_handle;
  typedef typename Polyhedron::Facet_iterator Facet_iterator;
  typedef typename Polyhedron::Halfedge_around_facet_circulator HF_circulator;

  std::map<Point, unsigned int> my_point_to_ind;
  unsigned int index = 0;

  std::set<Vertex_handle> all_vertices;

  Facet_iterator f;
  for (f = polyhedron.facets_begin();
       f != polyhedron.facets_end();
       ++f)
  {
    std::list<unsigned int> facet;

    std::vector<Vertex_handle> vertices_facet;
    HF_circulator hf = f->facet_begin();
    do
    {
      vertices_facet.push_back(hf->vertex());
    } while( ++hf != f->facet_begin() );

    for (unsigned i = 0; i != vertices_facet.size(); ++i)
    {
      Vertex_handle v = vertices_facet[i];
      const Point p = v->point();

      if ( all_vertices.find(v) == all_vertices.end() )
      {
        all_vertices.insert ( v );
        *pout++ = p;
        my_point_to_ind[p] = index++;
      }

      facet.push_back ( my_point_to_ind[p] );
    }
    *iout++ = facet;
  }
}

template <class Polyhedron>
void set_normals(Polyhedron& polyhedron, const std::vector<Vector>& normals)
{
  typedef typename Polyhedron::Facet_iterator Facet_iterator;
  unsigned int i = 0;
  Facet_iterator f;
  for (f = polyhedron.facets_begin();
       f !=  polyhedron.facets_end();
       ++f)
    f->normal() = normals[i++];
  
}

template <class Polyhedron>
void compute_normals(Polyhedron& polyhedron)
{
  typedef typename Polyhedron::Facet_iterator Facet_iterator;
  typedef typename Polyhedron::Vertex_handle Vertex_handle;

  Facet_iterator f;
  for (f = polyhedron.facets_begin();
       f !=  polyhedron.facets_end();
       ++f)
  {
    std::vector<Vertex_handle> vertices;
    f->get_vertices ( vertices );

    // init normal facets
    Vector vect_0_1 = vertices[1]->point() - vertices[0]->point();
    Vector vect_0_2 = vertices[2]->point() - vertices[0]->point();
    Vector cross_product = CGAL::cross_product ( vect_0_1 / std::sqrt( vect_0_1.squared_length() ), vect_0_2 / std::sqrt( vect_0_2.squared_length() ) ) ;
    cross_product = cross_product / std::sqrt(cross_product.squared_length());
    f->normal() = cross_product;

  }
}

/*template <class Polyhedron>
find_nearest_facet ( Facet_handle& fcenter, const Point& pcenter )
{


    Facet_handle fmiddle = Facet_handle();
    FT min = 1e308;
    Facet_iterator f;
    for ( f = this->facets_begin();
            f != this->facets_end();
            ++f )
    {
        const Point p = f->centroid();
        FT sqd_dist = CGAL::squared_distance ( pcenter, p );
        if ( sqd_dist < min )
        {
            min = sqd_dist;
            fmiddle = f;
        }
    }
    return fmiddle;
}*/

template <class F, class Polyhedron, class Tree>
F nearest_facet_from_query(const Point& pcenter, 
                           Polyhedron& polyhedron,
                           Tree& tree)
{
  return tree.closest_point_and_primitive(pcenter).second;
}


template <class Polyhedron, class F, class V>
void compute_boundary_vertices(const std::set<F>& facets,
                               std::set<V>& bvertices)
{
  typename std::set<F>::const_iterator fit;
  for (fit = facets.begin();
       fit != facets.end();
       ++fit)
  {
    F f = *fit;
    typename Polyhedron::Halfedge_around_facet_circulator hf = f->facet_begin();
    do
    {
      if (hf->is_border_edge())
        bvertices.insert(hf->vertex());
      else{
        F nf = hf->opposite()->facet();
        if (facets.find(nf) == facets.end())
          bvertices.insert(hf->vertex());
      }
    } while( ++hf != f->facet_begin());
  }
}

template <class Pt>
void compute_largest_circle_that_contains_query(const Pt& query,
                                                const std::set<Pt>& pts,
                                                std::set<Pt>& tri )
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef K::Point_2 Point_2;
  typedef CGAL::Delaunay_triangulation_2<K>  Delaunay_triangulation_2;
  typedef Delaunay_triangulation_2::Face_handle Face;

  Delaunay_triangulation_2 dt2;

  typename std::set<Pt>::const_iterator pit;
  for (pit = pts.begin();
       pit != pts.end();
       ++pit)
  {
    Pt p = *pit;
    dt2.push_back(Point_2(p.x(), p.y()));
  }

  //std::cout << dt2.number_of_vertices() << std::endl;

  Face f = dt2.locate(Point_2(query.x(), query.y()));

  for (int i = 0; i != 3; ++i)
  {
    tri.insert(f->vertex(i)->point());
    std::cout << i << ": " << f->vertex(i)->point() << std::endl;
  }

  
}




#endif // _ALGORITHMS_





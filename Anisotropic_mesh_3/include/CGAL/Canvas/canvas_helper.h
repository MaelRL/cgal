#ifndef CGAL_ANISOTROPIC_MESH_3_CANVAS_HELPER_H
#define CGAL_ANISOTROPIC_MESH_3_CANVAS_HELPER_H

#include <CGAL/intersection_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <Eigen/Dense>

#include <boost/array.hpp>

#include <cstddef>
#include <iterator>
#include <set>
#include <vector>

std::size_t fact(std::size_t n)
{
  return (n==0)?1:n*fact(n-1);
}

std::size_t combi(std::size_t n, std::size_t k)
{
  // n choose k
  // this will (silently) overflow for large n...

  if (k > n)
    return 0;

  std::size_t r = 1;
  for(std::size_t i=1; i<=k; ++i)
  {
    r *= n--;
    r /= i;
  }

  return r;
}

template<typename Metric>
Eigen::Matrix3d get_interpolated_transformation(const Metric& m0, const Metric& m1)
{
  // note that we return the transformation and not the full metric !
  return 0.5*(m0.get_transformation() + m1.get_transformation());
}

template<std::size_t k_max, typename Container>
std::vector<boost::array<typename Container::value_type, k_max> >
combinations(const Container& container,
             typename Container::const_iterator it,
             const std::size_t k)
{
  // the idea is to build the combinations recursively :
  // at one level we pick one element, and we'll then build the combination
  // as : this element added to the (k-1)-sized combinations
  typedef typename Container::value_type                            Contained;
  typedef std::vector<boost::array<Contained, k_max> >              Return_type;

  if(container.size() < k_max)
    return Return_type();
  CGAL_assertion(k > 0 && k <= container.size());

  std::vector<boost::array<Contained, k_max> > combis; // current level combis

  typename Container::const_iterator sit = it;
  typename Container::const_iterator send = container.begin();

  // we still need to pick k-1 elements afterwards so we can't pick any of the
  // (k-1) last elements
  std::advance(send, container.size() - (k-1));

  for(; sit!=send; ++sit)
  {
    if(k == 1)
    {
      boost::array<Contained, k_max> combi;

      // initialize the last element
      combi[k_max - k] = *sit;
      combis.push_back(combi);
    }
    else
    {
      // initialize the k_max-k element
      typename Container::const_iterator next_sit = sit;
      ++next_sit;
      Return_type lower_level_combis = combinations<k_max>(container, next_sit, k-1);
      for(std::size_t i=0; i<lower_level_combis.size(); ++i)
        lower_level_combis[i][k_max-k] = *(sit);

      combis.insert(combis.end(), lower_level_combis.begin(),
                                  lower_level_combis.end());
    }
  }
  return combis;
}

template<std::size_t k_max, typename Container>
std::vector<boost::array<typename Container::value_type, k_max> >
combinations(const Container& container,
             const std::size_t k = k_max)
{
  return combinations<k_max, Container>(container, container.begin(), k);
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

template<typename K>
bool is_triangle_intersected(const typename K::Triangle_3& triangle,
                             const std::set<typename K::Segment_3,
                                            Segment_3_comparer<K> >& segments)
{
  // todo version without EPECK, it's not needed to explicit the intersection
  // to know the intersection but simply if they intersect !
  // check what's done in the refinement queue for tet-tet intersections

  typedef typename K::Point_3                                  Point;
  typedef typename K::Segment_3                                Segment;

  typedef CGAL::Exact_predicates_exact_constructions_kernel    KExact;
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
          {
            if(local_intersected)
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
namespace Anisotropic_mesh_3
{

template<typename Cp>
struct Canvas_point_comparer
{
  bool operator()(Cp const * const cp1, Cp const * const cp2)
  {
    return cp1->distance_to_closest_seed() > cp2->distance_to_closest_seed();
  }
};

}
}

#endif // CGAL_ANISOTROPIC_MESH_3_CANVAS_HELPER_H

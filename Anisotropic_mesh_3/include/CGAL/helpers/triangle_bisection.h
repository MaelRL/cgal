#ifndef TRIANGLE_BISECTION_H
#define TRIANGLE_BISECTION_H

#include <cstddef>
#include <list>
#include <queue>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K>
struct Triangle_with_priority
{
  typedef typename K::Triangle_3 Triangle_3;

  Triangle_3 t;
  unsigned int m_index_longest_edge;
  double m_sqlmax;

  Triangle_with_priority(const Triangle_3& tr)
    : t(tr)
  {
    double sql_01 = (t[0] - t[1])*(t[0] - t[1]);
    double sql_12 = (t[1] - t[2])*(t[1] - t[2]);
    double sql_20 = (t[0] - t[2])*(t[0] - t[2]);
    m_sqlmax = std::max(sql_01, sql_12);
    m_sqlmax = std::max(m_sqlmax, sql_20);

    if(m_sqlmax == sql_01)
      m_index_longest_edge = 2;
    else if(m_sqlmax == sql_12)
      m_index_longest_edge = 0;
    else
      m_index_longest_edge = 1;
  }

  Triangle_with_priority() {}

  bool operator < (const Triangle_with_priority& tr) const
  {
    return m_sqlmax < tr.m_sqlmax;
  }
};

template<typename K>
void split_and_push(Triangle_with_priority<K> tr,
                    std::priority_queue<Triangle_with_priority<K> >& triangles)
{
  typedef typename K::Point_3     Point_3;
  typedef typename K::Triangle_3  Triangle_3;

  int i = tr.m_index_longest_edge;
  int pi = (i-1)%3;
  int ni = (i+1)%3;

  Point_3 m = CGAL::midpoint(tr.t[pi], tr.t[ni]);
  triangles.push(Triangle_with_priority<K>(Triangle_3(tr.t[i], tr.t[pi], m)));
  triangles.push(Triangle_with_priority<K>(Triangle_3(m, tr.t[ni], tr.t[i])));
#ifdef ANISO_DEBUG_SAMPLING
  std::cout << "split in two! : " << std::endl;
  std::cout << tr.t[i] << " || " << tr.t[pi] << " || " << m << std::endl;
  std::cout << m << " || " << tr.t[ni] << " || " << tr.t[i] << std::endl;
#endif
}


template<typename K>
void divide_in_smaller_triangles(std::priority_queue<Triangle_with_priority<K> >& triangles,
                                 const std::size_t max_triangle_n)
{
  while(!triangles.empty() && triangles.size() < max_triangle_n)
  {
    Triangle_with_priority<K> tr = triangles.top();
    triangles.pop();
    split_and_push(tr, triangles);
  }
}

template<typename K>
void get_sample_points(const typename K::Triangle_3& tri,
                       const unsigned int max_triangle_n,
                       std::vector<typename K::Point_3>& sample_points,
                       std::vector<typename K::FT>& coeffs)
{
  typedef typename K::FT        FT;
  typedef typename K::Point_3   Point_3;
  typedef typename K::Vector_3  Vector_3;

#ifdef ANISO_DEBUG_SAMPLING
  std::cout << "sampling " << max_triangle_n << " points on the facet : " << std::endl;
  std::cout << tri[0] << " || " << tri[1] << " || " << tri[2] << std::endl;
#endif

  //part of the barycentric coordinates computations that don't need p
  Vector_3 v00(tri[0], tri[1]);
  Vector_3 v11(tri[0], tri[2]);
  FT d00 = v00*v00;
  FT d01 = v00*v11;
  FT d11 = v11*v11;
  FT invdenom = 1./(d00 * d11 - d01 * d01);

  std::priority_queue<Triangle_with_priority<K> > triangles;
  triangles.push(Triangle_with_priority<K>(tri));

  if(max_triangle_n > 1)
    divide_in_smaller_triangles(triangles, max_triangle_n);

#ifdef ANISO_DEBUG_SAMPLING
  std::cout << "divided into " << triangles.size() << " triangles... Computing points..." << std::endl;
#endif

  while(!triangles.empty())
  {
    Triangle_with_priority<K> tr = triangles.top();
    triangles.pop();

    typename K::FT third = 1./3.;
    Point_3 p = CGAL::barycenter(tr.t[0], third, tr.t[1], third, tr.t[2], third);
    sample_points.push_back(p);

    //rest of the barycentric coordinates computations
    Vector_3 v22(tri[0], p);
    FT d20 = v22*v00;
    FT d21 = v22*v11;
    FT v = (d11 * d20 - d01 * d21) * invdenom;
    FT w = (d00 * d21 - d01 * d20) * invdenom;
    FT u = 1.0f - v - w;
    coeffs.push_back(u); coeffs.push_back(v); coeffs.push_back(w);
  }

#ifdef ANISO_DEBUG_SAMPLING
  std::cout << "check sizes : " << sample_points.size() << " " << coeffs.size() << std::endl;
#endif

}

}
}

#endif

#ifndef CGAL_ANISOTROPIC_MESH_3_PRIMAL_SIMPLEX_H
#define CGAL_ANISOTROPIC_MESH_3_PRIMAL_SIMPLEX_H

#include <CGAL/assertions.h>

#include <boost/array.hpp>
#include <boost/functional/hash.hpp>

#include <algorithm>
#include <cstddef>

namespace CGAL
{
namespace Anisotropic_mesh_3
{
// BSimplex is a boost::array that we must keep sorted because we use it
// as key in the unordered sets "primal_*" of the canvas

template<typename Canvas, std::size_t size_>
struct Primal_simplex
{
  typedef Primal_simplex<Canvas, size_>                Self;
  typedef boost::array<std::size_t, size_>             BSimplex;

  typedef typename Canvas::FT                          FT;
  typedef typename Canvas::Metric                      Metric;
  typedef typename Canvas::Canvas_point                Canvas_point;

  BSimplex m_simplex;
  const Canvas_point* m_dual_point; // witness of the primal (possibly a virtual canvas point)
  mutable bool m_is_intersected;

  const Canvas_point* dual_point() const { return m_dual_point; }
  const BSimplex& simplex() const { return m_simplex; }

  FT compute_distortion(const std::vector<Metric>& seeds_metrics) const
  {
    // return the largest distortion between two vertices of the primal simplex
    FT dist = 1.;

    typedef std::vector<boost::array<std::size_t, 2> >     RT_type;
    RT_type combis = combinations<2>(m_simplex);

    typename RT_type::const_iterator cit = combis.begin();
    for(; cit!=combis.end(); ++cit)
    {
      const boost::array<std::size_t, 2>& ids = *cit;
      const Metric& m0 = seeds_metrics[ids[0]];
      const Metric& m1 = seeds_metrics[ids[1]];
      FT loc_dist = m0.compute_distortion(m1);
      dist = (std::max)(loc_dist, dist);
    }
    return dist;
  }

  std::size_t size() const { return size_; }

  std::size_t& operator[](std::size_t i)
  {
    CGAL_assertion(i>=0 && i<size_);
    return m_simplex[i];
  }

  std::size_t operator[](std::size_t i) const
  {
    CGAL_assertion(i>=0 && i<size_);
    return m_simplex[i];
  }

  template<typename C, std::size_t s>
  friend std::ostream& operator<<(std::ostream& stream, const Primal_simplex<C,s>& ps);

  Primal_simplex() : m_simplex(), m_dual_point(NULL), m_is_intersected(false) { }
  Primal_simplex(const BSimplex& simplex_,
                 const Canvas_point* dual_point_ = NULL)
    :
      m_simplex(simplex_),
      m_dual_point(dual_point_),
      m_is_intersected(false)
  { }
};

template<typename Canvas, std::size_t size_>
std::ostream& operator<<(std::ostream& os, const Primal_simplex<Canvas, size_>& ps)
{
  os << "Primal simplex [";
  for(std::size_t i=0; i<size_; ++i)
    os << ps[i] << " ";
  os << "] dual point : " << ps.dual_point()->index()
     << " (" << ps.dual_point()->point() << ") ";
  os << "is intersected ? " << ps.m_is_intersected;

  return os;
}

template<typename Canvas, std::size_t size_>
struct Primal_simplex_hash
{
  typedef Primal_simplex<Canvas, size_>                      Psimplex;

  bool operator()(const Psimplex& psimplex) const
  {
    return boost::hash_range(psimplex.simplex().begin(), psimplex.simplex().end());
  }
};

template<typename Canvas, std::size_t size_>
struct Primal_simplex_comparer
{
  typedef Primal_simplex<Canvas, size_>                      Psimplex;

  bool operator()(const Psimplex& left, const Psimplex& right) const
  {
    CGAL_assertion(left.simplex().size() == right.simplex().size());
    return std::equal(left.simplex().begin(), left.simplex().end(),
                      right.simplex().begin());
  }
};

}  // Anisotropic_mesh_3
}  // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_PRIMAL_SIMPLEX_H

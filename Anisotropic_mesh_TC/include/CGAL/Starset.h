#ifndef CGAL_ANISOTROPIC_MESH_TC_STARSET_H
#define CGAL_ANISOTROPIC_MESH_TC_STARSET_H

#include <CGAL/Star.h>
#include <CGAL/helpers/metric_helper.h>
#include <CGAL/Metric_field.h>

#include <Eigen/Dense>

#include <iostream>
#include <fstream>
#include <queue>
#include <utility>

namespace CGAL
{
namespace Anisotropic_mesh_TC
{

enum Consistency_check_options
{
  TODO
};

template<typename Star>
bool is_infinite_vertex(int i)
{
  return (i == Star::infinite_vertex_index());
}

template<typename Kd, typename KD>
class Starset
{
  typedef Starset<Kd, KD>                                 Self;

public:
  typedef Tangent_star<Kd, KD>                            Star;
  typedef Star*                                           Star_handle;

  typedef typename Star::Index                            Index;
  typedef typename Star::dDim                             dDim;

  typedef typename Star::Traits                           Traits;
  typedef typename Star::FT                               FT;
  typedef typename Star::Point_d                          Point_d;
  typedef typename Star::Point_D                          Point_D;
  typedef typename Star::Vertex                           Vertex;
  typedef typename Star::Vertex_handle                    Vertex_handle;
  typedef typename Star::Vertex_h_iterator                Vertex_h_iterator;
  typedef typename Star::Vertex_handle_iterator           Vertex_handle_iterator;
  typedef typename Star::Face                             Face;
  typedef typename Star::Face_iterator                    Face_iterator;
  typedef typename Star::Full_cell                        Full_cell;
  typedef typename Star::Full_cell_handle                 Full_cell_handle;
  typedef typename Star::Full_cell_handle_iterator        Full_cell_handle_iterator;

  typedef std::vector<Star_handle>                        Star_vector;
  typedef typename Star_vector::const_iterator            const_iterator;
  typedef typename Star_vector::iterator                  iterator;

  typedef typename Star::Metric                           Metric;
  typedef Metric_field<Kd>*                               MF;

protected:
  Star_vector m_stars;
  MF m_mf;

public:
  const Star_vector& stars() const { return m_stars; }

  typename Star_vector::const_iterator begin() const { return m_stars.begin(); }
  typename Star_vector::iterator begin() { return m_stars.begin(); }
  typename Star_vector::const_iterator end() const { return m_stars.end(); }
  typename Star_vector::iterator end() { return m_stars.end(); }
  Star_handle back() { return m_stars.back(); }

  void pop_back() { return m_stars.pop_back(); }
  void push_back(Star_handle s) { m_stars.push_back(s); }

  std::size_t size() const { return m_stars.size(); }
  bool empty() const { return m_stars.empty(); }

  template<typename I>
  Star_handle operator[](I i) const { return m_stars[i]; }

public:
  std::size_t total_number_of_vertices() const
  {
    std::size_t nbv = 0;
    for(unsigned int i = 0; i < m_stars.size(); i++)
      nbv += m_stars[i]->number_of_vertices();
    return nbv;
  }

  //consistency
  //Face
  bool is_consistent(const Full_cell_handle& fch,
                     std::vector<bool>& inconsistent_points,
                     const bool verbose = false) const
  {
    bool retval = true;
    bool told = false;
    for(int i=0; i<fch->maximal_dimension(); i++)
    {
      int index = fch->vertex(i)->data();

      if(!m_stars[index]->has_cell(fch))
      {
        if(verbose)
        {
          if(!told)
          {
            std::cout << " inconsistent : ";
            told = true;
          }
          std::cout << fch->vertex(0)->data() << " ";
          std::cout << fch->vertex(1)->data() << " ";
          std::cout << fch->vertex(2)->data() << " ";
          std::cout << "c not in S_" << index << ", ";
        }
        retval = false;
        inconsistent_points[i] = true;
      }
    }
    if(verbose && !retval)
      std::cout << "." << std::endl;
    return retval;
  }

  bool is_consistent(const Full_cell_handle& fch,
                     const bool verbose = false) const
  {
    std::vector<bool> not_used(dDim::value);
    return is_consistent(fch, not_used, verbose);
  }

  //star
  bool is_consistent(Star_handle star,
                     const bool verbose = false) const
  {
    typename Star::Full_cell_handle_iterator fcit = star->finite_incident_full_cells_begin();
    typename Star::Full_cell_handle_iterator fcend = star->finite_incident_full_cells_end();
    for(; fcit!=fcend; fcit++)
    {
      if(!is_consistent(*fcit, verbose))
        return false;
    }

    return true;
  }

  //star set
  bool is_consistent(const bool verbose = true) const
  {
    bool retval = true;
    int counter = 0;
    std::size_t N = m_stars.size();
    for(std::size_t i = 0; i < N; i++)
    {
      Star_handle star = m_stars[i];
      if(!is_consistent(star, verbose))
      {
        retval = false;
        counter++;
      }
    }
    std::cout << counter << " inconsistent stars" << std::endl;
    return retval;
  }

  // distortion
  FT compute_distortion(const Full_cell& fc) const
  {
    FT distortion = 1.;
    for (int i = 0; i < 3; i++)
    {
      int i2 = (i+1) % 3;
      distortion = (std::max)(distortion,
                              m_stars[fc.vertex(i)->data()]->metric().compute_distortion(
                     m_stars[fc.vertex(i2)->data()]->metric()));
    }
    return distortion;
  }

  FT star_distortion(Star_handle star) const
  {
    FT max_distortion = 0.;
    Vertex_handle_iterator vit = star->finite_adjacent_vertices_begin();
    Vertex_handle_iterator vend = star->finite_adjacent_vertices_end();
    for(; vit!=vend; ++vit)
    {
      if(is_infinite_vertex((*vit)->data()))
        continue;

      FT distortion = m_stars[(*vit)->data()]->metric().compute_distortion(star->metric());
      max_distortion = (std::max)(distortion, max_distortion);
    }

    return max_distortion;
  }

  void clear()
  {
    m_stars.clear();
  }

  void update_caches()
  {
    for(std::size_t si=0; si<m_stars.size(); ++si)
      m_stars[si]->update_star_caches();
  }

  void rebuild() // brute forcy
  {
    //todo repair the commented parts
//    for(std::size_t si=0; si<m_stars.size(); ++si)
//    {
//      Star_handle star_i = m_stars[si];
//      star_i->reset();
//    }

    for(std::size_t si=0; si<m_stars.size(); ++si)
    {
      Star_handle star_i = m_stars[si];
      for(std::size_t sj=0; sj<m_stars.size(); ++sj)
      {
        if(si == sj)
          continue;
        Star_handle star_j = m_stars[sj];
        star_i->insert_to_star(star_j->m_center_S, sj, false/*no cond*/);
      }
      //star_i->clean();
    }
  }

  void insert_in_stars(const Point_d& p)
  {
    Star_handle s = new Star(p, m_stars.size(), m_mf->compute_metric(p));
    for(std::size_t i=0; i<m_stars.size(); ++i)
    {
      m_stars[i]->insert_to_star(s->m_center_S, s->index(), false/*no cond*/);
      s->insert_to_star(m_stars[i]->m_center_S, m_stars[i]->index(), false);
    }
    m_stars.push_back(s);
  }

  void insert_in_stars(const Point_d& p, const Metric& m)
  {
    Star_handle s = new Star(p, m_stars.size(), m);
    for(std::size_t i=0; i<m_stars.size(); ++i)
    {
      m_stars[i]->insert_to_star(s->m_center_S, s->index(), false/*no cond*/);
      s->insert_to_star(m_stars[i]->m_center_S, m_stars[i]->index(), false);
    }
    m_stars.push_back(s);
  }

  // Criteria ------------------------------------------------------------------

/* FIXME: define properly r_0 before. For example, it should probably use
 * the original points in R^d. Currently (below) are points in the tangent plane
 * of the star...
  FT compute_circumradius(const Full_cell_handle fch) const
  {
    const int d = Star::dDim::value;
    if(fch->maximal_dimension() != d)
    {
      std::cout << "compute circumradius w/ fch->max_dim != d..." << std::endl;
      assert(0);
    }

    assert(fch->data().second); // make sure the cell circumcenter has been computed
    const Point_d& c = fch->data().first;

    //compute distance between m_center and c IN THE METRIC OF M_CENTER
    Eigen::Matrix<FT,d,1> v;

    //looping just to make sure the point is eq to all (in their respective metrics)
    FT dist = 0.;
    for(int i=0; i<=d; ++i)
    {
      Point_d pi = m_stars[fch->vertex(i)->data()]->m_center;
      for(int j=0; j<d; ++j)
        v(j) = c[j]-pi[j];

      v = m_stars[fch->vertex(i)->data()]->metric().get_transformation() * v;
      dist  = std::sqrt(v.transpose() * v);
      std::cout << "metric dist: " << fch->vertex(i)->data() << " " << dist << std::endl;
    }
    return dist;
  }
*/

  FT compute_volume(const Full_cell_handle fch) const
  {
    const std::size_t d = dDim::value;
    if(fch->maximal_dimension() != d)
      return 0.;

    FT den = 1./(FT) fact(d);
    Eigen::Matrix<FT,d,d> m = Eigen::Matrix<FT, d, d>::Zero();
    for(int i=0; i<d; ++i)
    {
      Point_d v0 = m_stars[fch->vertex(0)->data()]->center_point();
      for(int j=0; j<d; ++j)
      {
        Point_d vj = m_stars[fch->vertex(j+1)->data()]->center_point();
        m(i,j) = vj[i]-v0[i];
      }
    }
    return std::abs(den * m.determinant());
  }

  Starset() : m_stars() { }
  Starset(const Star_vector& stars_) : m_stars(stars_) { }

  Starset(const MF mf_) : m_stars(), m_mf(mf_) { }

private:
  Starset(const Self& src);
  Self& operator=(const Self& src);
};

} // Anisotropic_mesh_TC
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_TC_STARSET_H

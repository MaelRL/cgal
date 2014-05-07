#ifndef CGAL_ANISOTROPIC_MESH_3_STARSET_H
#define CGAL_ANISOTROPIC_MESH_3_STARSET_H

#include <CGAL/Stretched_Delaunay_3.h>
#include <CGAL/helpers/combinatorics_helper.h>

#include <Eigen/Dense>

#include <iostream>
#include <fstream>
#include <utility>


namespace CGAL
{
namespace Anisotropic_mesh_3
{

enum Consistency_check_options
{
  CELLS_ONLY, //this could be merged with the similar drawing enum TODO
  CELLS_AND_FACETS,
  FACETS_ONLY
};

template<typename Star>
bool is_infinite_vertex(int i)
{
  return (i == Star::infinite_vertex_index());
}

template<typename K, typename KExact = K>
class Starset
{
  typedef Starset<K, KExact>                            Self;

public:
  typedef Stretched_Delaunay_3<K, KExact>               Star;
  typedef Star*                                         Star_handle;
  typedef std::vector<Star_handle>                      Star_vector;

  typedef typename Star::Traits                         Traits;
  typedef typename Star::FT                             FT;
  typedef typename Star::Vertex                         Vertex;
  typedef typename Star::Vertex_handle                  Vertex_handle;
  typedef typename Star::Vertex_handle_handle           Vertex_handle_handle;
  typedef typename Star::Facet                          Facet;
  typedef typename Star::Facet_handle                   Facet_handle;
  typedef typename Star::Facet_set_iterator             Facet_set_iterator;
  typedef typename Star::Cell_handle                    Cell_handle;
  typedef typename Star::Cell_handle_handle             Cell_handle_handle;
  typedef typename Star::Index                          Index;
  typedef typename Star::Point_3                        Point_3;
  typedef typename Star::TPoint_3                       TPoint_3;
  typedef typename Star::Vector_3                       Vector_3;

  typedef typename Star_vector::const_iterator          const_iterator;
  typedef typename Star_vector::iterator                iterator;

protected:
  Star_vector m_stars;

public:
  Star_vector& star_vector() { return m_stars; }

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
  unsigned int number_of_surface_stars() const
  {
    unsigned int nb = 0;
    for(unsigned int i = 0; i < m_stars.size(); ++i)
      if(m_stars[i]->is_surface_star())
        nb++;
    return nb;
  }

  std::size_t total_number_of_vertices() const
  {
    std::size_t nbv = 0;
    for(unsigned int i = 0; i < m_stars.size(); i++)
      nbv += m_stars[i]->number_of_vertices();
    return nbv;
  }

  std::size_t count_restricted_facets() const
  {
    std::set<Facet_ijk> facets;
    for(unsigned int i = 0; i < m_stars.size(); i++)
    {
      Facet_set_iterator fit = m_stars[i]->begin_restricted_facets();
      Facet_set_iterator fitend = m_stars[i]->end_restricted_facets();
      for(; fit != fitend; fit++)
        facets.insert(Facet_ijk(*fit));
    }
    return facets.size();
  }

  int number_of_tets_in_star_set() const
  {
    int count = 0;
    for (int i = 0; i < (int)m_stars.size(); i++)
    {
      Star_handle si = m_stars[i];
      typename Star::Cell_handle_handle cit = si->begin_finite_star_cells();
      typename Star::Cell_handle_handle cend = si->end_finite_star_cells();
      for(; cit != cend; cit++)
        // if(si->is_inside(*cit))
        count++;
    }
    return count;
  }

  //consistency
  //facet
  bool is_consistent(const Facet& f,
                     std::vector<bool> &inconsistent_points,
                     const bool verbose = false) const
  {
    bool retval = true;
    bool facet_told = false;
    for(int i = 1; i < 4; i++)
    {
      int index = f.first->vertex((f.second + i) % 4)->info();
      if(!m_stars[index]->has_facet(f))
      {
        if(verbose)
        {
          if(!facet_told)
          {
            m_stars[index]->facet_indices(f);
            std::cout << " inconsistent : ";
            facet_told = true;
          }
          std::cout << "f not in S_" << index << ", ";
        }
        retval = false;
        inconsistent_points[i-1] = true;
      }
    }
    if(verbose && !retval) std::cout << "." << std::endl;
    return retval;
  }

  bool is_consistent(const Facet& f,
                     const bool verbose = false) const
  {
    std::vector<bool> not_used(3);
    return is_consistent(f, not_used, verbose);
  }

  //cell
  bool is_consistent(Cell_handle c,
                     std::vector<bool>& inconsistent_points,
                     const bool verbose = false) const
  {
    bool retval = true;
    bool cell_told = false;
    for(int i = 0; i < 4; i++)
    {
      int index = c->vertex(i)->info();
      if(is_infinite_vertex<Star>(index))
        continue;
      if(!m_stars[index]->has_cell(c))
      {
        if(verbose)
        {
          if(!cell_told)
          {
            m_stars[index]->cell_indices(c);
            std::cout << " inconsistent : ";
            cell_told = true;
          }
          std::cout << c->vertex(0)->info() << " " << c->vertex(1)->info() << " ";
          std::cout << c->vertex(2)->info() << " " << c->vertex(3)->info();
          std::cout << " not in S_" << index << ", ";
        }
        retval = false;
        inconsistent_points[i] = true;
      }
    }
    if(verbose && !retval)
      std::cout << "." << std::endl;
    return retval;
  }

  bool is_consistent(Cell_handle c,
                     const bool verbose = false) const
  {
    std::vector<bool> not_used(4);
    return is_consistent(c, not_used, verbose);
  }

  //star
  bool is_consistent(Star_handle star,
                     const bool verbose = false,
                     Consistency_check_options option = CELLS_AND_FACETS) const
  {
    if(option != FACETS_ONLY)
    {
      Cell_handle_handle cit = star->begin_finite_star_cells();
      Cell_handle_handle citend = star->end_finite_star_cells();
      for(; cit != citend; cit++)
      {
        if(!star->is_inside(*cit))
          continue;
        if(!is_consistent(*cit, verbose))
          return false;
      }
    }

    if(option != CELLS_ONLY)
    {
      Facet_set_iterator fit = star->begin_restricted_facets();
      Facet_set_iterator fitend = star->end_restricted_facets();
      for(; fit != fitend; fit++)
        if(!is_consistent(*fit, verbose))
          return false;
    }

    return true;
  }

  //star set
  bool is_consistent(const bool verbose = true,
                     Consistency_check_options option = CELLS_AND_FACETS) const
  {
    std::size_t N = m_stars.size();
    for(std::size_t i = 0; i < N; i++)
    {
      Star_handle star = m_stars[i];
      if(!is_consistent(star, verbose, option))
        return false;
    }
    return true;
  }

  // distortion
  FT compute_distortion(const Facet& f) const
  {
    FT distortion = 1.;
    int index = f.second;
    Cell_handle c = f.first;
    for (int i = 0; i < 3; i++)
    {
      int i1 = (index + i + 1) % 4;
      int i2 = (index + (i + 1) % 3 + 1) % 4;
      distortion = (std::max)(distortion,
                              m_stars[c->vertex(i1)->info()]->metric().compute_distortion(
                     m_stars[c->vertex(i2)->info()]->metric()));
    }
    return distortion;
  }

  FT compute_distortion(const Cell_handle& c) const
  {
    FT distortion = 1.0;
    for(int i = 0; i < 4; i++)
    {
      for(int j = i + 1; j < 4; j++)
      {
        distortion = (std::max)(distortion,
                                m_stars[c->vertex(i)->info()]->metric().compute_distortion(
                       m_stars[c->vertex(j)->info()]->metric()));
      }
    }
    return distortion;
  }

  FT star_distortion(Star_handle star) const
  {
    FT max_distortion = 0.;

    /*
    Facet_set_iterator fit = star->begin_restricted_facets();
    Facet_set_iterator fitend = star->end_restricted_facets();
    for(; fit != fitend; fit++)
    {
      Facet f = *fit;
      max_distortion = (std::max)(max_distortion, compute_distortion(f));
    }
    return max_distortion;
  */

    Vertex_handle_handle vit = star->begin_neighboring_vertices();
    Vertex_handle_handle vend = star->end_neighboring_vertices();
    for(; vit!=vend; ++vit)
    {
      if(is_infinite_vertex((*vit)->info()))
        continue;
      if(!(get_star((*vit)->info())->is_surface_star())) //todo this is bad
        continue;

      FT distortion = m_stars[(*vit)->info()]->metric().compute_distortion(star->metric());
      max_distortion = (std::max)(distortion, max_distortion);
    }

    return max_distortion;
  }

  void debug_show_distortions() const
  {
    for(std::size_t i = 0; i < m_stars.size(); ++i)
    {
      Star_handle s = m_stars[i];
      if(!s->is_surface_star())
        continue;
      std::cout << "  " << i << " : ";
      Facet_set_iterator fit = s->begin_restricted_facets();
      Facet_set_iterator fend = s->end_restricted_facets();
      for(; fit != fend; ++fit)
      {
        std::cout << "(";
        Facet f = *fit;
        if(is_surface_facet(f))
        {
          for (int i = 0; i < 3; i++)
          {
            int index_1 = (f.second + i + 1) % 4;
            int index_2 = (f.second + (i + 1) % 3 + 1) % 4;
            FT distortion =
                m_stars[f.first->vertex(index_1)->info()]->metric().compute_distortion(
                  m_stars[f.first->vertex(index_2)->info()]->metric());

            FT costheta = std::abs(
                            m_stars[f.first->vertex(index_1)->info()]->metric().get_vmin()
                          * m_stars[f.first->vertex(index_2)->info()]->metric().get_vmin());

            std::cout << (distortion - 1./costheta) << "\t";
          }
        }
        std::cout << ")";
      }
      std::cout << std::endl;
    }
  }

  double average_facet_distortion(bool verbose = true) const
  {
    double avg_coh_dis = 0., min_coh_dis = 1e30, max_coh_dis = -1e30;
    double avg_incoh_dis = 0., min_incoh_dis = 1e30, max_incoh_dis = -1e30;
    int nb_coh = 0, nb_incoh = 0;

    std::set<Facet_ijk> done;
    for(std::size_t i = 0; i < m_stars.size(); i++)
    {
      Star_handle star = m_stars[i];
      if(!star->is_surface_star())
        continue;

      Facet_set_iterator fit = star->begin_restricted_facets();
      Facet_set_iterator fitend = star->end_restricted_facets();
      for(; fit != fitend; fit++)
      {
        Facet f = *fit;
        std::pair<typename std::set<Facet_ijk>::iterator, bool> is_insert_successful;
        is_insert_successful = done.insert(Facet_ijk(f));
        if(!is_insert_successful.second)
          continue;

        FT facet_distortion = compute_distortion(f);

        if(is_consistent(f))
        {
          nb_coh++;
          avg_coh_dis += facet_distortion;
          if(facet_distortion > max_coh_dis)
            max_coh_dis = facet_distortion;
          if(facet_distortion < min_coh_dis)
            min_coh_dis = facet_distortion;
        }
        else
        {
          nb_incoh++;
          avg_incoh_dis += facet_distortion;
          if(facet_distortion > max_incoh_dis)
            max_incoh_dis = facet_distortion;
          if(facet_distortion < min_incoh_dis)
            min_incoh_dis = facet_distortion;
        }
      }
    }

    if(verbose)
    {
      std::cout << "SUM UP OF THE DISTORTION (FACET):" << std::endl;
      std::cout << nb_coh << " coherent facets with avg: " << avg_coh_dis/(double) nb_coh << std::endl;
      std::cout << "min: " << min_coh_dis << " max: " << max_coh_dis << std::endl;
      std::cout << nb_incoh << " incoherent facets with avg " << avg_incoh_dis/(double) nb_incoh << std::endl;
      std::cout << "min: " << min_incoh_dis << " max: " << max_incoh_dis << std::endl;
      std::cout << "---------------------------------------------" << std::endl;
    }

    return avg_coh_dis/(double) nb_coh;
  }

  double average_cell_distortion(bool verbose = true) const
  {
    double avg_coh_dis = 0., min_coh_dis = 1e30, max_coh_dis = -1e30;
    double avg_incoh_dis = 0., min_incoh_dis = 1e30, max_incoh_dis = -1e30;
    int nb_coh = 0, nb_incoh = 0;

    std::set<Cell_ijkl> done;
    for(std::size_t i = 0; i < m_stars.size(); i++)
    {
      Star_handle star = m_stars[i];
      if(!star->is_surface_star())
        continue;

      Cell_handle_handle ci = star->begin_finite_star_cells();
      Cell_handle_handle ciend = star->end_finite_star_cells();
      for(; ci != ciend; ci++)
      {
        Cell_handle c = *ci;
        if(!star->is_inside(c))
          continue;

        std::pair<typename std::set<Cell_ijkl>::iterator, bool> is_insert_successful;
        is_insert_successful = done.insert(Cell_ijkl(c));
        if(!is_insert_successful.second)
          continue;

        FT cell_distortion = compute_distortion(c);

        if(is_consistent(c))
        {
          nb_coh++;
          avg_coh_dis += cell_distortion;
          if(cell_distortion > max_coh_dis)
            max_coh_dis = cell_distortion;
          if(cell_distortion < min_coh_dis)
            min_coh_dis = cell_distortion;
        }
        else
        {
          nb_incoh++;
          avg_incoh_dis += cell_distortion;
          if(cell_distortion > max_incoh_dis)
            max_incoh_dis = cell_distortion;
          if(cell_distortion < min_incoh_dis)
            min_incoh_dis = cell_distortion;
        }
      }
    }

    if(verbose)
    {
      std::cout << "SUM UP OF THE DISTORTION (CELL):" << std::endl;
      std::cout << nb_coh << " coherent cells with avg: " << avg_coh_dis/(double) nb_coh << std::endl;
      std::cout << "min: " << min_coh_dis << " max: " << max_coh_dis << std::endl;
      std::cout << nb_incoh << " incoherent cells with avg " << avg_incoh_dis/(double) nb_incoh << std::endl;
      std::cout << "min: " << min_incoh_dis << " max: " << max_incoh_dis << std::endl;
      std::cout << "---------------------------------------------" << std::endl;
    }

    return avg_coh_dis/(double) nb_coh;
  }

  double average_star_distortion(bool verbose = true) const
  {
    double avg_coh_dis = 0., min_coh_dis = 1e30, max_coh_dis = -1e30;
    double avg_incoh_dis = 0., min_incoh_dis = 1e30, max_incoh_dis = -1e30;
    int nb_coh = 0, nb_incoh = 0;

    for(std::size_t i=0; i< m_stars.size(); ++i)
    {
      Star_handle si = m_stars[i];
      if(!(si->is_surface_star()))
        continue;

      FT max_distortion = star_distortion(si);

      if(is_consistent(si))
      {
        nb_coh++;
        avg_coh_dis += max_distortion;
        if(max_distortion > max_coh_dis)
          max_coh_dis = max_distortion;
        if(max_distortion < min_coh_dis)
          min_coh_dis = max_distortion;
      }
      else
      {
        nb_incoh++;
        avg_incoh_dis += max_distortion;
        if(max_distortion > max_incoh_dis)
          max_incoh_dis = max_distortion;
        if(max_distortion < min_incoh_dis)
          min_incoh_dis = max_distortion;
      }
    }

    if(verbose)
    {
      std::cout << "SUM UP OF THE DISTORTION (STAR) :" << std::endl;
      std::cout << nb_coh << " coherent (surface) stars with avg: " << avg_coh_dis/(double) nb_coh << std::endl;
      std::cout << "min: " << min_coh_dis << " max: " << max_coh_dis << std::endl;
      std::cout << nb_incoh << " incoherent (surface) stars with avg " << avg_incoh_dis/(double) nb_incoh << std::endl;
      std::cout << "min: " << min_incoh_dis << " max: " << max_incoh_dis << std::endl;
      std::cout << "---------------------------------------------" << std::endl;
    }

    return avg_coh_dis/(double) nb_coh;
  }

  Starset() : m_stars() { }
  Starset(const Star_vector& stars_) : m_stars(stars_) { }

private:
  Starset(const Self& src);
  Self& operator=(const Self& src);
};

template<typename K, typename Constrain_surface, typename Metric, typename Criteria, typename KExact = K>
class Starset_with_info : public Starset<K, KExact>
{
  typedef Starset_with_info<K, Constrain_surface, Metric, Criteria, KExact>  Self;
  typedef Starset<K, KExact>                                                 Base;

  const Constrain_surface* m_pConstrain;
  const Metric* m_metric;
  const Criteria* m_criteria;

public:
  const Constrain_surface* constrain_surface() const { return m_pConstrain; }
  const Metric* metric_field() const { return m_metric; }
  const Criteria* criteria() const { return m_criteria; }

public:
  void set_criteria(const Criteria* criteria_) { m_criteria = criteria_; }
  void update_stars_criteria()
  {
    for(std::size_t i = 0; i < this->m_stars.size(); i++)
      this->m_stars[i]->set_criteria(m_criteria);
  }

  Starset_with_info(const Constrain_surface* constrain_surface_,
                    const Metric* metric_,
                    const Criteria* criteria_)
    :
      Base(),
      m_pConstrain(constrain_surface_),
      m_metric(metric_),
      m_criteria(criteria_)
  { }

  Starset_with_info(const typename Base::Star_vector& stars_,
                    const Constrain_surface* constrain_surface_,
                    const Metric* metric_,
                    const Criteria* criteria_)
    :
      Base(stars_),
      m_pConstrain(constrain_surface_),
      m_metric(metric_),
      m_criteria(criteria_)
  { }

private:
  Starset_with_info(const Self& src);
  Self& operator=(const Self& src);
};

} // Anisotropic_mesh_3
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_STARSET_H

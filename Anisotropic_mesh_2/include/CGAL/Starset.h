#ifndef CGAL_ANISOTROPIC_MESH_2_STARSET_H
#define CGAL_ANISOTROPIC_MESH_2_STARSET_H

#include <CGAL/Stretched_Delaunay_2.h>

#include <Eigen/Dense>

#include <iostream>
#include <fstream>
#include <utility>

namespace CGAL
{
namespace Anisotropic_mesh_2
{

enum Consistency_check_options
{
  FACES_ONLY
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
  typedef Stretched_Delaunay_2<K, KExact>               Star;
  typedef Star*                                         Star_handle;
  typedef std::vector<Star_handle>                      Star_vector;

  typedef typename Star::Traits                         Traits;
  typedef typename Star::FT                             FT;
  typedef typename Star::Vertex                         Vertex;
  typedef typename Star::Vertex_handle                  Vertex_handle;
  typedef typename Star::Vertex_handle_handle           Vertex_handle_handle;
  typedef typename Star::Face                           Face;
  typedef typename Star::Face_handle                    Face_handle;
  typedef typename Star::Face_handle_handle             Face_handle_handle;
  typedef typename Star::Index                          Index;
  typedef typename Star::Point_2                        Point_2;
  typedef typename Star::TPoint_2                       TPoint_2;
  typedef typename Star::Vector_2                       Vector_2;

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

  //consistency
  //Face_handle
  bool is_consistent(const Face_handle& fh,
                     std::vector<bool> &inconsistent_points,
                     const bool verbose = false) const
  {
    bool retval = true;
    bool told = false;
    for(int i=0; i<3; i++)
    {
      int index = fh->vertex(i)->info();

      if(!m_stars[index]->has_face(fh))
      {
        if(verbose)
        {
          if(!told)
          {
            m_stars[index]->face_indices(fh);
            std::cout << " inconsistent : ";
            told = true;
          }
          std::cout << "f not in S_" << index << ", ";
        }
        retval = false;
        inconsistent_points[i] = true;
      }
    }
    if(verbose && !retval) std::cout << "." << std::endl;
    return retval;
  }

  bool is_consistent(const Face_handle& fh,
                     const bool verbose = false) const
  {
    std::vector<bool> not_used(3);
    return is_consistent(fh, not_used, verbose);
  }

  //star
  bool is_consistent(Star_handle star,
                     const bool verbose = false,
                     Consistency_check_options option = FACES_ONLY) const
  {
    typename Star::Face_handle_handle fit = star->finite_adjacent_faces_begin();
    typename Star::Face_handle_handle fend = star->finite_adjacent_faces_end();
    for(; fit!=fend; fit++)
      if(!is_consistent(*fit, verbose))
        return false;

    return true;
  }

  //star set
  bool is_consistent(const bool verbose = true,
                     Consistency_check_options option = FACES_ONLY) const
  {
    bool retval = true;
    int counter = 0;
    std::size_t N = m_stars.size();
    for(std::size_t i = 0; i < N; i++)
    {
      Star_handle star = m_stars[i];
      if(!is_consistent(star, verbose, option))
      {
        retval = false;
        counter++;
      }
    }
    std::cout << counter << " inconsistent stars" << std::endl;
    return retval;
  }

  // distortion
  FT compute_distortion(const Face_handle& fh) const
  {
    FT distortion = 1.;
    for (int i = 0; i < 3; i++)
    {
      int i2 = (i+1) % 3;
      distortion = (std::max)(distortion,
                              m_stars[fh->vertex(i)->info()]->metric().compute_distortion(
                     m_stars[fh->vertex(i2)->info()]->metric()));
    }
    return distortion;
  }

  FT star_distortion(Star_handle star) const
  {
    FT max_distortion = 0.;
    Vertex_handle_handle vit = star->begin_neighboring_vertices();
    Vertex_handle_handle vend = star->end_neighboring_vertices();
    for(; vit!=vend; ++vit)
    {
      if(is_infinite_vertex((*vit)->info()))
        continue;
      if(!(m_stars[(*vit)->info()]->is_surface_star())) //todo this is bad
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
      typename Star::Face_handle_handle fit = s->finite_adjacent_faces_begin();
      typename Star::Face_handle_handle fend = s->finite_adjacent_faces_end();
      for(; fit!=fend; ++fit)
      {
        std::cout << "(";

        for (int i = 0; i < 3; i++)
        {
          int i2 = (i+1)%3;
          FT distortion =
              m_stars[fit->vertex(i)->info()]->metric().compute_distortion(
                                    m_stars[fit->vertex(i2)->info()]->metric());

          FT costheta = std::abs(
                          m_stars[fit->vertex(i)->info()]->metric().get_vmin()
                        * m_stars[fit->vertex(i2)->info()]->metric().get_vmin());

          std::cout << (distortion - 1./costheta) << "\t";
        }

        std::cout << ")";
      }
      std::cout << std::endl;
    }
  }

  double average_face_distortion(bool verbose = true) const
  {
    double avg_coh_dis = 0., min_coh_dis = 1e30, max_coh_dis = -1e30;
    double avg_incoh_dis = 0., min_incoh_dis = 1e30, max_incoh_dis = -1e30;
    int nb_coh = 0, nb_incoh = 0;

    Facet_ijk_unordered_set done;
    for(std::size_t i = 0; i < m_stars.size(); i++)
    {
      Star_handle star = m_stars[i];
      if(!star->is_surface_star())
        continue;

      typename Star::Face_handle_handle fit = star->finite_adjacent_faces_begin();
      typename Star::Face_handle_handle fend = star->finite_adjacent_faces_end();
      for(; fit!=fend; fit++)
      {
        Face_handle fh = *fit;
        std::pair<typename Facet_ijk_unordered_set::iterator, bool> is_insert_successful;
        is_insert_successful = done.insert(Facet_ijk(fh));
        if(!is_insert_successful.second)
          continue;

        FT face_distortion = compute_distortion(fh);

        if(is_consistent(fh))
        {
          nb_coh++;
          avg_coh_dis += face_distortion;
          if(face_distortion > max_coh_dis)
            max_coh_dis = face_distortion;
          if(face_distortion < min_coh_dis)
            min_coh_dis = face_distortion;
        }
        else
        {
          nb_incoh++;
          avg_incoh_dis += face_distortion;
          if(face_distortion > max_incoh_dis)
            max_incoh_dis = face_distortion;
          if(face_distortion < min_incoh_dis)
            min_incoh_dis = face_distortion;
        }
      }
    }

    if(verbose)
    {
      std::cout << "SUM UP OF THE DISTORTION (Face_handle):" << std::endl;
      std::cout << nb_coh << " coherent Face_handles with avg: " << avg_coh_dis/(double) nb_coh << std::endl;
      std::cout << "min: " << min_coh_dis << " max: " << max_coh_dis << std::endl;
      std::cout << nb_incoh << " incoherent Face_handles with avg " << avg_incoh_dis/(double) nb_incoh << std::endl;
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

  void clear()
  {
    m_stars.clear();
  }

  void rebuild() // brute forcy
  {
    for(std::size_t si=0; si<m_stars.size(); ++si)
    {
      Star_handle star_i = m_stars[si];
      star_i->reset();
    }

    for(std::size_t si=0; si<m_stars.size(); ++si)
    {
      Star_handle star_i = m_stars[si];
      for(std::size_t sj=0; sj<m_stars.size(); ++sj)
      {
        Star_handle star_j = m_stars[sj];
        star_i->insert_to_star(star_j->center_point(), sj, false/*no cond*/);
      }
      star_i->clean();
    }
  }

  Starset() : m_stars() { }
  Starset(const Star_vector& stars_) : m_stars(stars_) { }

private:
  Starset(const Self& src);
  Self& operator=(const Self& src);
};

template<typename K, typename Domain, typename Metric, typename Criteria, typename KExact = K>
class Starset_with_info : public Starset<K, KExact>
{
  typedef Starset_with_info<K, Domain, Metric, Criteria, KExact>  Self;
  typedef Starset<K, KExact>                                                 Base;

  const Domain* m_pdomain;
  const Metric* m_metric;
  const Criteria* m_criteria;

public:
  const Domain* domain() const { return m_pdomain; }
  const Metric* metric_field() const { return m_metric; }
  const Criteria* criteria() const { return m_criteria; }

public:
  void set_criteria(const Criteria* criteria_) { m_criteria = criteria_; }
  void update_stars_criteria()
  {
    for(std::size_t i = 0; i < this->m_stars.size(); i++)
      this->m_stars[i]->set_criteria(m_criteria);
  }

  Starset_with_info(const Domain* domain_,
                    const Metric* metric_,
                    const Criteria* criteria_)
    :
      Base(),
      m_pdomain(domain_),
      m_metric(metric_),
      m_criteria(criteria_)
  { }

  Starset_with_info(const typename Base::Star_vector& stars_,
                    const Domain* domain_,
                    const Metric* metric_,
                    const Criteria* criteria_)
    :
      Base(stars_),
      m_pdomain(domain_),
      m_metric(metric_),
      m_criteria(criteria_)
  { }

private:
  Starset_with_info(const Self& src);
  Self& operator=(const Self& src);
};

} // Anisotropic_mesh_2
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_2_STARSET_H

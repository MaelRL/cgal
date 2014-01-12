#ifndef CGAL_ANISOTROPIC_MESH_3_STAR_CONSISTENCY_H
#define CGAL_ANISOTROPIC_MESH_3_STAR_CONSISTENCY_H

#include <cstddef>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

enum Consistency_check_options{ CELLS_ONLY,
                                CELLS_AND_FACETS,
                                FACETS_ONLY};

template<typename Star>
bool is_infinite_vertex(int i)
{
  return (i == Star::infinite_vertex_index());
}

//Basic consistency functions
//facet
template<typename Star>
bool is_consistent(const typename std::vector<Star*>& stars,
                   const typename Star::Facet& f,
                   std::vector<bool> &inconsistent_points,
                   const bool verbose = false)
{
  bool retval = true;
  bool facet_told = false;
  for(int i = 1; i < 4; i++)
  {
    int index = f.first->vertex((f.second + i) % 4)->info();
    if(!stars[index]->has_facet(f))
    {
      if(verbose)
      {
        if(!facet_told)
        {
          stars[index]->facet_indices(f);
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

template<typename Star>
bool is_consistent(const typename std::vector<Star*>& stars,
                   const typename Star::Facet& f,
                   const bool verbose = false)
{
  std::vector<bool> not_used(3);
  return is_consistent(stars, f, not_used, verbose);
}

//cell
template<typename Star>
bool is_consistent(const typename std::vector<Star*>& stars,
                   const typename Star::Cell_handle& c,
                   std::vector<bool>& inconsistent_points,
                   const bool verbose = false)
{
  bool retval = true;
  bool cell_told = false;
  for(int i = 0; i < 4; i++)
  {
    int index = c->vertex(i)->info();
    if(is_infinite_vertex<Star>(index))
      continue;
    if(!stars[index]->has_cell(c))
    {
      if(verbose)
      {
        if(!cell_told)
        {
          stars[index]->cell_indices(c);
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

template<typename Star>
bool is_consistent(const typename std::vector<Star*>& stars,
                   const typename Star::Cell_handle& c,
                   const bool verbose = false)
{
  std::vector<bool> not_used(4);
  return is_consistent(stars, c, not_used, verbose);
}

//star
template<typename Star>
bool is_consistent(const typename std::vector<Star*>& stars,
                   Star* star,
                   const bool verbose = false,
                   Consistency_check_options option = CELLS_AND_FACETS)
{
  if(option != FACETS_ONLY)
  {
    typename Star::Cell_handle_handle cit = star->begin_finite_star_cells();
    typename Star::Cell_handle_handle citend = star->end_finite_star_cells();
    for(; cit != citend; cit++)
    {
      if(!star->is_inside(*cit))
        continue;
      if(!is_consistent(stars, *cit, verbose))
        return false;
    }
  }

  if(option != CELLS_ONLY)
  {
    typename Star::Facet_set_iterator fit = star->begin_restricted_facets();
    typename Star::Facet_set_iterator fitend = star->end_restricted_facets();
    for(; fit != fitend; fit++)
      if(!is_consistent(stars, *fit, verbose))
        return false;
  }

  return true;
}

//star set
template<typename Star>
bool is_consistent(const typename std::vector<Star*>& stars,
                   const bool verbose = true,
                   Consistency_check_options option = CELLS_AND_FACETS)
{
  std::size_t N = stars.size();
  for(std::size_t i = 0; i < N; i++)
  {
    Star* star = stars[i];
    if(!is_consistent(stars, star, verbose, option))
      return false;
  }
  return true;
}


} //namespace Aniso
} //namespace CGAL

#endif

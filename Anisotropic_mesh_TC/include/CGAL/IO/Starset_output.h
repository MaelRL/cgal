#ifndef CGAL_ANISOTROPIC_MESH_TC_STAR_SET_OUTPUT_H
#define CGAL_ANISOTROPIC_MESH_TC_STAR_SET_OUTPUT_H

#include <CGAL/helpers/combinatorics_helper.h>

#include <fstream>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_TC
{

template<typename Starset>
void dump(const Starset& stars)
{
  std::ofstream fx("dump.txt");

  std::size_t ns = stars.size();
  fx << ns << std::endl;
  for(std::size_t i=0; i<ns; ++i)
    fx << stars[i]->center_point() << std::endl;

  for(std::size_t i=0; i<ns; ++i)
  {
    typename Starset::Star_handle star_i = stars[i];
    typename Starset::Vertex_handle_iterator nsi = star_i->finite_adjacent_vertices_begin();
    typename Starset::Vertex_handle_iterator nsiend = star_i->finite_adjacent_vertices_end();
    fx << nsiend - nsi;
    for(; nsi != nsiend; nsi++)
      fx << " " << (*nsi)->data();
    fx << std::endl;
  }
}

template<typename Starset, typename Cell_set>
int get_simplices_to_draw(const Starset& starset,
                          Cell_set& output_cells,
                          const bool consistent_only = false)
{
  const std::size_t d = Starset::dDim::value;
  unsigned int nb_inconsistent_cells = 0;
  unsigned int nb_restricted_cells = 0;

  typename Starset::const_iterator it = starset.begin();
  typename Starset::const_iterator itend = starset.end();
  for (; it != itend; ++it)
  {
    typename Starset::Star_handle star = *it;
//    std::cout << "Star: " << star->m_center_v->data() << std::endl;

    if(star->current_dimension() != d)
      continue;

    typename Starset::Full_cell_handle_iterator fchi = star->finite_incident_full_cells_begin();
    typename Starset::Full_cell_handle_iterator fend = star->finite_incident_full_cells_end();
    for(; fchi!=fend; ++fchi)
    {
      typename Starset::Full_cell_handle fch = *fchi;

      if(star->is_inside(fch, starset.stars())) // recomputes the dual intersection
        nb_restricted_cells++;

      bool is_consistent = true; //starset.is_consistent(fch, true);
      if(!is_consistent)
        nb_inconsistent_cells++;

      if(!consistent_only || is_consistent)
      {
        boost::array<int, d+1> ids;

        int id = 0;
        typename Starset::Vertex_h_iterator v = fch->vertices_begin();
        typename Starset::Vertex_h_iterator vend = fch->vertices_end();
        for(; v!=vend; ++v)
        {
          typename Starset::Vertex_handle vh = *v;
//          std::cout << vh->data() << " ";
          ids[id++] = vh->data();
        }
//        std::cout << std::endl;

        Ordered_simplex_base<d+1> cell(ids);
        output_cells.insert(cell);
      }
    }
  }

  std::cout << "drawing " << output_cells.size() << " (" << nb_restricted_cells << ")" << std::endl;
  return nb_inconsistent_cells;
}

template<typename Starset>
void output_off(const Starset& stars,
                std::ofstream& fx,
                const bool consistent_only = false)
{
  std::cout << "Saving as .off... @ " << stars.size() << std::endl;
  std::cerr << "Saving as .off... @ " << stars.size() << std::endl;

  const std::size_t d = Starset::dDim::value;
  typedef typename Simplex_unordered_set<d+1>::type Cell_unordered_set;
  Cell_unordered_set output_cells;
  unsigned int nb_inconsistent_cells =
      get_simplices_to_draw(stars, output_cells, consistent_only);

  if(nb_inconsistent_cells > 0)
    std::cout << "Warning OFF : there are " << nb_inconsistent_cells << " inconsistent cells in the ouput mesh.\n";

  int tri_per_simplex = (d+1)*d*(d-1)/6.; // 3 from d+1
  int tn = tri_per_simplex * output_cells.size();

  fx << "OFF" << std::endl;
  fx << stars.size() << " " << tn << " 0" << std::endl;
  for(unsigned int i=0; i<stars.size(); i++)
  {
    const typename Starset::Point_d& p = stars[i]->center_point();
    for(int j=0; j<d; ++j)
      fx << p[j] << " ";
    for(int j=0; j<(3-d); ++j) // ugly trick: fill with 0
      fx << "0 ";
    fx << std::endl;
  }

  typename Cell_unordered_set::iterator fcit = output_cells.begin();
  typename Cell_unordered_set::iterator fcend = output_cells.end();
  for (; fcit != fcend; fcit++)
  {
    CGAL::Combination_enumerator<int> combi(3, 0, d+1);
    for(; !combi.finished(); combi++)
      fx << "3 " << (*fcit)[combi[0]] << " " << (*fcit)[combi[1]] << " " << (*fcit)[combi[2]] << std::endl;
  }
}

template<typename Starset>
void output_medit(const Starset& stars,
                  std::ofstream& fx,
                  const bool consistent_only = false)
{
  std::cout << "Saving medit... @ " << stars.size() << std::endl;

  const std::size_t d = Starset::dDim::value;
  typedef typename Simplex_unordered_set<d+1>::type Cell_unordered_set;
  Cell_unordered_set output_cells;
  unsigned int nb_inconsistent_cells =
      get_simplices_to_draw(stars, output_cells, consistent_only);

  if(nb_inconsistent_cells > 0)
    std::cout << "Warning Medit : there are " << nb_inconsistent_cells << " inconsistent cells in the ouput mesh.\n";

  fx << "MeshVersionFormatted 1\n";
  fx << "Dimension " << d << std::endl;

  fx << "Vertices\n";
  fx << stars.size() << std::endl;
  for(std::size_t i=0; i<stars.size(); i++)
  {
    const typename Starset::Point_d& p = stars[i]->center_point();
    for(int j=0; j<d; ++j)
      fx << p[j] << " ";
    fx << (i+1) << std::endl;
  }

  int tri_per_simplex = (d+1)*d*(d-1)/6.; // 3 from d+1
  int tn = tri_per_simplex * output_cells.size();
  fx << "Triangles\n" << tn << std::endl;
  typename Cell_unordered_set::iterator fcit = output_cells.begin();
  typename Cell_unordered_set::iterator fcend = output_cells.end();
  for (; fcit != fcend; fcit++)
  {
    CGAL::Combination_enumerator<int> combi(3, 0, d+1);
    for(; !combi.finished(); combi++)
      fx << (*fcit)[combi[0]]+1 << " " << (*fcit)[combi[1]]+1 << " " << (*fcit)[combi[2]]+1 << " 1" << std::endl;
  }
  fx << "End" << std::endl;
}

template<typename Starset>
void output_off(const Starset &stars,
                const char* f,
                const bool consistent_only = false)
{
  std::ofstream out(f);
  return output_off(stars, out, consistent_only);
}

template<typename Starset>
void output_medit(const Starset &stars,
                  const char* f,
                  const bool consistent_only = false)
{
  std::ofstream out(f);
  return output_medit(stars, out, consistent_only);
}

}  // Anisotropic_mesh_TC
}  // CGAL

#endif

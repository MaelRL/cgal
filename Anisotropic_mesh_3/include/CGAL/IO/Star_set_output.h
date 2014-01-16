#ifndef CGAL_ANISOTROPIC_MESH_3_STAR_SET_OUTPUT_H
#define CGAL_ANISOTROPIC_MESH_3_STAR_SET_OUTPUT_H

#include <vector>

#include <CGAL/Star_consistency.h>
#include <CGAL/IO/Output_cells.h>
#include <CGAL/IO/Output_facets.h>

#include <CGAL/helpers/combinatorics_helper.h>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename Star>
void dump(const std::vector<Star*>& stars)
{
  std::ofstream fx("dump.txt");

  std::size_t ns = stars.size();
  fx << ns << std::endl;
  for(std::size_t i=0; i<ns; ++i)
    fx << stars[i]->center_point() << std::endl;

  for(std::size_t i=0; i<ns; ++i)
  {
    Star* star_i = stars.get_star(i);
    typename Star::Vertex_handle_handle nsi = star_i->begin_neighboring_vertices();
    typename Star::Vertex_handle_handle nsiend = star_i->end_neighboring_vertices();
    fx << nsiend - nsi;
    for(; nsi != nsiend; nsi++)
      fx << " " << (*nsi)->info();
    fx << std::endl;
  }
}

template<typename Star>
void output_off(const std::vector<Star*>& stars,
                std::ofstream& fx,
                const bool consistent_only = true)
{
  std::cout << "Saving as .off..." << std::endl;
  std::map<typename Star::Index, int> match_indices;//because we won't use all of them
  int off_index = 0; // the corresponding index in .off
  std::vector<typename Star::Point_3> points;
  Output_facets output_facets;
  Output_cells output_cells;
  unsigned int nb_inconsistent_stars = 0;

  typename std::vector<Star*>::const_iterator it = stars.begin();
  typename std::vector<Star*>::const_iterator itend = stars.end();
  for (; it != itend; ++it)
  {
    Star* star = *it;

    points.push_back(star->center_point());
    match_indices[star->index_in_star_set()] = off_index++;

    typename Star::Cell_handle_handle ci = star->begin_star_cells();
    typename Star::Cell_handle_handle ciend = star->end_star_cells();
    for (; ci != ciend; ++ci)
    {
      typename Star::Cell_handle c = *ci;
      if(star->is_infinite(c))
        continue;
      if(!consistent_only || is_consistent(stars, c))
      {
        output_cells.insert(c->vertex(0)->info(), c->vertex(1)->info(),
                            c->vertex(2)->info(), c->vertex(3)->info());
        output_facets.insert(c->vertex(0)->info(), c->vertex(1)->info(), c->vertex(2)->info());
        output_facets.insert(c->vertex(0)->info(), c->vertex(2)->info(), c->vertex(3)->info());
        output_facets.insert(c->vertex(0)->info(), c->vertex(1)->info(), c->vertex(3)->info());
        output_facets.insert(c->vertex(3)->info(), c->vertex(1)->info(), c->vertex(2)->info());
      }
      else
      {
        nb_inconsistent_stars++;
        break;
      }
    }
  }

  if(nb_inconsistent_stars > 0)
    std::cout << "Warning : there are " << nb_inconsistent_stars << " inconsistent stars in the ouput mesh.\n";

  fx << "OFF" << std::endl;
  fx << points.size() << " " << output_facets.size() << " " << output_cells.size() << std::endl;
  for(unsigned int i = 0; i < points.size(); i++)
    fx << points[i] << std::endl;

  Output_facets::Facet_handle ofi = output_facets.begin();
  Output_facets::Facet_handle ofiend = output_facets.end();
  for (; ofi != ofiend; ofi++)
  {
    fx << "3  " << match_indices[ofi->vertices[0] ]
        << " "   << match_indices[ofi->vertices[1] ]
        << " "   << match_indices[ofi->vertices[2] ] << std::endl;
  }

  Output_cells::Cell_handle oci = output_cells.begin();
  Output_cells::Cell_handle ociend = output_cells.end();
  for (; oci != ociend; oci++)
  {
    fx << "4  " << match_indices[oci->vertices[0] ]
        << " "   << match_indices[oci->vertices[1] ]
        << " "   << match_indices[oci->vertices[2] ]
        << " "   << match_indices[oci->vertices[3] ] << std::endl;
  }
}

template<typename Star>
void output_surface_off(const std::vector<Star*>& stars,
                        std::ofstream& fx,
                        const bool consistent_only = true)
{
  std::cout << "Saving surface off..." << std::endl;
  std::map<typename Star::Index, int> match_indices;//because we won't use all of them
  int off_index = 0; // the corresponding index in .off
  std::vector<typename Star::Point_3> points;
  Output_facets output_facets;
  unsigned int nb_inconsistent_stars = 0;

  typename std::vector<Star*>::const_iterator si = stars.begin();
  typename std::vector<Star*>::const_iterator siend = stars.end();
  for (; si != siend; ++si)
  {
    Star* star = *si;

    points.push_back(star->center_point());
    match_indices[star->index_in_star_set()] = off_index++;

    typename Star::Facet_set_iterator fi = star->begin_restricted_facets();
    typename Star::Facet_set_iterator fiend = star->end_restricted_facets();
    for (; fi != fiend; ++fi)
    {
      typename Star::Point_3 cc;
      star->compute_dual_intersection(*fi, cc);

      if(!consistent_only || is_consistent(stars, *fi))
        output_facets.insert(fi->first->vertex((fi->second + 1) % 4)->info(),
                             fi->first->vertex((fi->second + 2) % 4)->info(),
                             fi->first->vertex((fi->second + 3) % 4)->info());
      else
      {
        nb_inconsistent_stars++;
        break;
      }
    }
  }

  if(nb_inconsistent_stars > 0)
    std::cout << "Warning : there are " << nb_inconsistent_stars << " inconsistent stars in the ouput mesh.\n";

  fx << "OFF" << std::endl;
  fx << points.size() << " " << output_facets.size() << " " << 0 << std::endl;
  for(unsigned int i = 0; i < points.size(); i++)
    fx << points[i] << std::endl;

  Output_facets::Facet_handle ofi = output_facets.begin();
  Output_facets::Facet_handle ofiend = output_facets.end();
  for (; ofi != ofiend; ofi++)
  {
    fx << "3  " << match_indices[ofi->vertices[0] ]
       << " "   << match_indices[ofi->vertices[1] ]
       << " "   << match_indices[ofi->vertices[2] ] << std::endl;
  }
}

/*
  void output()
  {
    typename std::ofstream fx("mesh.volume.cgal");
    fx << 3 << std::endl;
    fx << m_points.size() << std::endl;
    for(int i = 0; i < (int)m_points.size(); i++)
      fx << m_points[i] << std::endl;

    Output_cells output_cells;
    Star_iterator it = m_stars.begin();
    Star_iterator itend = m_stars.end();
    for(; it != itend; it++)
    {
      Star_handle star = *it;
      Cell_handle_handle ci = star->begin_finite_star_cells();
      Cell_handle_handle ciend = star->end_finite_star_cells();
      for(; ci != ciend; ci++)
      {
        if(!star->is_inside(*ci))
          continue;
        bool consistent = true;
        for(int i = 0; i < 4; i++)
          if(!m_stars[(*ci)->vertex(i)->info()]->has_cell(*ci))
          {
            consistent = false;
            break;
          }
          if(consistent)
            output_cells.insert(
            (*ci)->vertex(0)->info(), (*ci)->vertex(1)->info(),
            (*ci)->vertex(2)->info(), (*ci)->vertex(3)->info());
      }
    }
    Output_cells::Cell_handle oci = output_cells.begin();
    Output_cells::Cell_handle ociend = output_cells.end();
    fx << output_cells.size() << std::endl;
    for(; oci != ociend; oci++)
    {
      fx << oci->vertices[0] + 1 << " " << oci->vertices[1] + 1 << " "
        << oci->vertices[2] + 1 << " " << oci->vertices[3] + 1 << std::endl;
    }
  }
*/

template<typename Star>
void output_surface_medit(const std::vector<Star*>& stars,
                          std::ofstream& fx)
{
  std::cout << "Saving surface medit..." << std::endl;
  fx << "MeshVersionFormatted 1\n\n";
  fx << "Dimension\n3" << std::endl;

  fx << "Vertices" << std::endl;
  fx << (stars.size()) << std::endl;
  for (int i = 0; i < (int)stars.size(); i++)
    fx << stars[i]->center_point() << " " << (i+1) << std::endl;

  fx << "Triangles" << std::endl;
  std::vector<Facet_ijk> facets;
  for(int i = 0; i < (int)stars.size(); i++)
  {
    Star* s = stars[i];
    typename Star::Facet_set_iterator fit = s->begin_restricted_facets();
    typename Star::Facet_set_iterator fend = s->end_restricted_facets();
    for(; fit != fend; ++fit)
      facets.push_back(Facet_ijk(*fit));
  }
  std::sort(facets.begin(), facets.end());
  std::unique(facets.begin(), facets.end());

  fx << (facets.size()) << std::endl;
  for(int i = 0; i < (int)facets.size(); i++)
   fx << (facets[i].vertex(0)+1) << " "
      << (facets[i].vertex(1)+1) << " "
      << (facets[i].vertex(2)+1) << " 1" << std::endl;
}

template<typename Star>
void output_medit(const std::vector<Star*>& stars,
                  std::ofstream& fx)
{
  std::cout << "Saving medit..." << std::endl;
  unsigned int nb_inconsistent_stars = 0;
  fx << "MeshVersionFormatted 1\n";
  fx << "Dimension 3\n";

  fx << "Vertices\n";
  fx << stars.size() << std::endl;
  for(std::size_t i = 0; i < stars.size(); i++)
    fx << stars[i]->center_point() << " " << (i+1) << std::endl; // warning : indices start at 1 in Medit

  Output_cells output_cells;
  typename std::vector<Star*>::const_iterator sit = stars.begin();
  typename std::vector<Star*>::const_iterator sitend = stars.end();
  for(; sit != sitend; sit++)
  {
    Star* star = *sit;
    typename Star::Cell_handle_handle ci = star->begin_finite_star_cells();
    typename Star::Cell_handle_handle ciend = star->end_finite_star_cells();
    for(; ci != ciend; ci++)
    {
      if(!star->is_inside(*ci))
        continue;
      if(is_consistent(stars, *ci))
        output_cells.insert((*ci)->vertex(0)->info(), (*ci)->vertex(1)->info(),
                            (*ci)->vertex(2)->info(), (*ci)->vertex(3)->info());
      else
      {
        nb_inconsistent_stars++;
        break;
      }
    }
  }
  fx << "Tetrahedra\n";
  fx << output_cells.size() << std::endl;
  Output_cells::Cell_handle oci = output_cells.begin();
  Output_cells::Cell_handle ociend = output_cells.end();
  for(; oci != ociend; oci++)
  {
    fx << (oci->vertices[0] + 1) << " " << (oci->vertices[1] + 1) << " ";
    fx << (oci->vertices[2] + 1) << " " << (oci->vertices[3] + 1) << " 1" << std::endl;
  }
  if(nb_inconsistent_stars > 0)
    std::cout << "Warning : there are " << nb_inconsistent_stars << " inconsistent stars in the ouput mesh.\n";
}

}  // Anisotropic_mesh_3
}  // CGAL

#endif

#ifndef CGAL_ANISOTROPIC_MESH_3_STAR_SET_OUTPUT_H
#define CGAL_ANISOTROPIC_MESH_3_STAR_SET_OUTPUT_H

#include <CGAL/Star_consistency.h>

#include <CGAL/helpers/combinatorics_helper.h>

#include <vector>
#include <set>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename Star>
typename Star::FT
compute_distortion_t(const std::vector<Star*>& stars,
                     const typename Star::Cell_handle& c)
{
  typename Star::FT distortion = 1.0;
  for(int i = 0; i < 4; i++)
  {
    for(int j = i + 1; j < 4; j++)
      distortion = (std::max)(distortion,
                              stars[c->vertex(i)->info()]->metric().compute_distortion(
                     stars[c->vertex(j)->info()]->metric()));
  }
  return distortion;
}

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
    Star* star_i = stars[i];
    typename Star::Vertex_handle_handle nsi = star_i->begin_neighboring_vertices();
    typename Star::Vertex_handle_handle nsiend = star_i->end_neighboring_vertices();
    fx << nsiend - nsi;
    for(; nsi != nsiend; nsi++)
      fx << " " << (*nsi)->info();
    fx << std::endl;
  }
}

template<typename Star>
void output_star_off(const std::vector<Star*>& stars,
                     std::ofstream& fx,
                     typename Star::Index i)
{
  Star* star = stars[i];
  std::map<typename Star::Index, int> match_indices;
  int off_index = 0;
  std::vector<typename Star::Point_3> points;
  std::set<Facet_ijk> output_facets;
  std::set<Cell_ijkl> output_cells;

  typename Star::Cell_handle_handle ci = star->begin_finite_star_cells();
  typename Star::Cell_handle_handle ciend = star->end_finite_star_cells();
  for(; ci != ciend; ci++)
  {
    typename Star::Cell_handle c = *ci;
    if(!star->is_inside(c))
      continue;

    boost::array<typename Star::Index, 4> ids;
    for(int i=0; i<4; i++)
    {
      ids[i] = c->vertex(i)->info();
      std::pair<typename std::map<typename Star::Index, int>::iterator, bool> inserted;
      inserted = match_indices.insert(std::make_pair(ids[i], off_index));
      if(inserted.second)
      {
        points.push_back(star->metric().inverse_transform(c->vertex(i)->point()));
        off_index++;
      }
    }

    output_cells.insert(Cell_ijkl(c));
    output_facets.insert(Facet_ijk(ids[0], ids[1], ids[2]));
    output_facets.insert(Facet_ijk(ids[0], ids[2], ids[3]));
    output_facets.insert(Facet_ijk(ids[0], ids[1], ids[3]));
    output_facets.insert(Facet_ijk(ids[3], ids[1], ids[2]));
  }

  fx << "OFF" << std::endl;
  fx << points.size() << " " << output_facets.size() << " " << output_cells.size() << std::endl;
  for(unsigned int i = 0; i < points.size(); i++)
    fx << points[i] << std::endl;

  typename std::set<Facet_ijk>::iterator fit = output_facets.begin();
  typename std::set<Facet_ijk>::iterator fitend = output_facets.end();
  for (; fit != fitend; fit++)
  {
    fx << "3  " << match_indices[fit->vertices()[0] ]
        << " "   << match_indices[fit->vertices()[1] ]
        << " "   << match_indices[fit->vertices()[2] ] << std::endl;
  }

  typename std::set<Cell_ijkl>::iterator  cit = output_cells.begin();
  typename std::set<Cell_ijkl>::iterator  citend = output_cells.end();
  for (; cit != citend; cit++)
  {
    fx << "4  " << match_indices[cit->vertices()[0] ]
        << " "   << match_indices[cit->vertices()[1] ]
        << " "   << match_indices[cit->vertices()[2] ]
        << " "   << match_indices[cit->vertices()[3] ] << std::endl;
  }
}

template<typename Star>
void output_off(const std::vector<Star*>& stars,
                std::ofstream& fx,
                const bool consistent_only = true)
{
  std::cout << "Saving as .off..." << std::endl;
  std::vector<typename Star::Point_3> points;
  std::set<Facet_ijk> output_facets;
  std::set<Cell_ijkl> output_cells;
  unsigned int nb_inconsistent_stars = 0;

  typename std::vector<Star*>::const_iterator it = stars.begin();
  typename std::vector<Star*>::const_iterator itend = stars.end();
  for (; it != itend; ++it)
  {
    Star* star = *it;
    points.push_back(star->center_point());

    typename Star::Cell_handle_handle ci = star->begin_finite_star_cells();
    typename Star::Cell_handle_handle ciend = star->end_finite_star_cells();
    for (; ci != ciend; ++ci)
    {
      typename Star::Cell_handle c = *ci;
      if(!star->is_inside(*ci))
        continue;
      if(!consistent_only || is_consistent(stars, c))
      {
        output_cells.insert(Cell_ijkl(c));
        output_facets.insert(Facet_ijk(c->vertex(0)->info(), c->vertex(1)->info(),
                                       c->vertex(2)->info()));
        output_facets.insert(Facet_ijk(c->vertex(0)->info(), c->vertex(2)->info(),
                                       c->vertex(3)->info()));
        output_facets.insert(Facet_ijk(c->vertex(0)->info(), c->vertex(1)->info(),
                                       c->vertex(3)->info()));
        output_facets.insert(Facet_ijk(c->vertex(3)->info(), c->vertex(1)->info(),
                                       c->vertex(2)->info()));
      }
      else
      {
        nb_inconsistent_stars++;
        break;
      }
    }
  }

  if(nb_inconsistent_stars > 0)
    std::cout << "Warning OFF : there are " << nb_inconsistent_stars << " inconsistent stars in the ouput mesh.\n";

  fx << "OFF" << std::endl;
  fx << points.size() << " " << output_facets.size() << " " << output_cells.size() << std::endl;
  for(unsigned int i = 0; i < points.size(); i++)
    fx << points[i] << std::endl;

  typename std::set<Facet_ijk>::iterator fit = output_facets.begin();
  typename std::set<Facet_ijk>::iterator fitend = output_facets.end();
  for (; fit != fitend; fit++)
  {
    fx << "3  " << fit->vertices()[0]
       << " "   << fit->vertices()[1]
       << " "   << fit->vertices()[2] << std::endl;
  }

  typename std::set<Cell_ijkl>::iterator  cit = output_cells.begin();
  typename std::set<Cell_ijkl>::iterator  citend = output_cells.end();
  for (; cit != citend; cit++)
  {
    fx << "4  " << cit->vertices()[0]
       << " "   << cit->vertices()[1]
       << " "   << cit->vertices()[2]
       << " "   << cit->vertices()[3] << std::endl;
  }
}

template<typename Star>
void output_surface_off(const std::vector<Star*>& stars,
                        std::ofstream& fx,
                        const bool consistent_only = true)
{
  std::cout << "Saving surface off..." << std::endl;
  std::vector<typename Star::Point_3> points;
  std::set<Facet_ijk> output_facets;
  unsigned int nb_inconsistent_stars = 0;

  typename std::vector<Star*>::const_iterator si = stars.begin();
  typename std::vector<Star*>::const_iterator siend = stars.end();
  for (; si != siend; ++si)
  {
    Star* star = *si;
    points.push_back(star->center_point());

    typename Star::Facet_set_iterator fi = star->begin_restricted_facets();
    typename Star::Facet_set_iterator fiend = star->end_restricted_facets();
    for (; fi != fiend; ++fi)
    {
      typename Star::Point_3 cc;
      star->compute_dual_intersection(*fi, cc);

      if(!consistent_only || is_consistent(stars, *fi))
        output_facets.insert(Facet_ijk(*fi));
      else
      {
        nb_inconsistent_stars++;
        break;
      }
    }
  }

  if(nb_inconsistent_stars > 0)
    std::cout << "Warning Surf Off: there are " << nb_inconsistent_stars << " inconsistent stars in the ouput mesh.\n";

  fx << "OFF" << std::endl;
  fx << points.size() << " " << output_facets.size() << " " << 0 << std::endl;
  for(unsigned int i = 0; i < points.size(); i++)
    fx << points[i] << std::endl;

  typename std::set<Facet_ijk>::iterator fit = output_facets.begin();
  typename std::set<Facet_ijk>::iterator fitend = output_facets.end();
  for (; fit != fitend; fit++)
  {
    fx << "3  " << fit->vertices()[0]
       << " "   << fit->vertices()[1]
       << " "   << fit->vertices()[2] << std::endl;
  }
}

/*
  template<typename Star>
  void output(const std::vector<Star*>& stars)
  {
    typename std::ofstream fx("mesh.volume.cgal");
    fx << 3 << std::endl;
    fx << m_points.size() << std::endl;
    for(int i = 0; i < (int)m_points.size(); i++)
      fx << m_points[i] << std::endl;

    std::set<Cell_ijkl> output_cells;
    typename Star::Star_iterator it = stars.begin();
    typename Star::Star_iterator itend = stars.end();
    for(; it != itend; it++)
    {
      Star* star = *it;
      typename Star::Cell_handle_handle ci = star->begin_finite_star_cells();
      typename Star::Cell_handle_handle ciend = star->end_finite_star_cells();
      for(; ci != ciend; ci++)
      {
        if(!star->is_inside(*ci))
          continue;
        bool consistent = true;
        for(int i = 0; i < 4; i++)
          if(!stars[(*ci)->vertex(i)->info()]->has_cell(*ci))
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

    fx << output_cells.size() << std::endl;
    typename std::set<Cell_ijk>  cit = output_cells.begin();
    typename std::set<Cell_ijk>  citend = output_cells.end();
    for (; cit != citend; cit++)
    {
      fx << cit->vertices[0]+1 << " " << cit->vertices[1]+1 << " "
         << cit->vertices[2]+1 << " " << cit->vertices[3]+1 << std::endl;
    }
  }
*/

template<typename Star>
void output_surface_medit(const std::vector<Star*>& stars,
                          std::ofstream& fx,
                          const bool positive_vol = false,
                          const bool fake_3D = false) // link all facets to a dummy point
{
  typename Star::Point_3 anchor_point(0., 0., 0.);

  std::cout << "Saving surface medit..." << std::endl;
  fx << "MeshVersionFormatted 1\n\n";
  fx << "Dimension\n3" << std::endl;

  fx << "Vertices" << std::endl;
  fx << (stars.size()+fake_3D) << std::endl;
  for (int i = 0; i < (int)stars.size(); i++)
    fx << stars[i]->center_point() << " " << (i+1) << std::endl;
  if(fake_3D)
    fx << anchor_point << (stars.size()+1) << std::endl;

  std::set<Facet_ijk> facets;
  for(int i = 0; i < (int)stars.size(); i++)
  {
    Star* s = stars[i];
    typename Star::Facet_set_iterator fit = s->begin_restricted_facets();
    typename Star::Facet_set_iterator fend = s->end_restricted_facets();
    for(; fit != fend; ++fit)
      facets.insert(Facet_ijk(*fit));
  }

  fx << "Triangles" << std::endl;
  fx << (facets.size()) << std::endl;
  typename std::set<Facet_ijk>::iterator fit = facets.begin();
  typename std::set<Facet_ijk>::iterator fitend = facets.end();
  for(; fit!=fitend; fit++)
  {
    Facet_ijk f = *fit;
    int n0 = f.vertices()[0];
    int n1 = f.vertices()[1];
    int n2 = f.vertices()[2];
    typename Star::Point_3 p0 = stars[n0]->center_point();
    typename Star::Point_3 p1 = stars[n1]->center_point();
    typename Star::Point_3 p2 = stars[n2]->center_point();
    typename Star::Traits::Compute_area_3 o;
    if(!positive_vol || o(p0, p1, p2) > 0)
      fx << (n0+1) << " " << (n1+1) << " " << (n2+1) << " 1" << std::endl;
    else
      fx << (n1+1) << " " << (n0+1) << " " << (n2+1) << " 1" << std::endl;
  }

  if(fake_3D)
  {
    fx << "Tetrahedra" << std::endl;
    fx << (facets.size()) << std::endl;
    fit = facets.begin();
    for(; fit!=fitend; fit++)
    {
      Facet_ijk f = *fit;
      int n0 = f.vertices()[0];
      int n1 = f.vertices()[1];
      int n2 = f.vertices()[2];
      int n3 = stars.size(); //anchor_id
      typename Star::Point_3 p0 = stars[n0]->center_point();
      typename Star::Point_3 p1 = stars[n1]->center_point();
      typename Star::Point_3 p2 = stars[n2]->center_point();
      typename Star::Traits::Compute_volume_3 o;
      if(!positive_vol || o(p0, p1, p2, anchor_point) > 0)
        fx << (n0+1) << " " << (n1+1) << " " << (n2+1) << " " << (n3+1) << " 1" << std::endl;
      else
        fx << (n1+1) << " " << (n0+1) << " " << (n2+1) << " " << (n3+1) << " 1" << std::endl;
    }
  }
}

template<typename Star>
void output_medit(const std::vector<Star*>& stars,
                  std::ofstream& fx,
                  const bool consistent_only = true,
                  const bool positive_vol = false)
{
  std::cout << "Saving medit..." << std::endl;
  unsigned int nb_inconsistent_stars = 0;

  std::map<Cell_ijkl, int> colors;
  int step_n = 100;

  fx << "MeshVersionFormatted 1\n";
  fx << "Dimension 3\n";

  fx << "Vertices\n";
  fx << stars.size() << std::endl;
  for(std::size_t i = 0; i < stars.size(); i++)
    fx << stars[i]->center_point() << " " << (i+1) << std::endl; //indices start at 1 in Medit

  std::set<Cell_ijkl> output_cells;
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
      if(!consistent_only || is_consistent(stars, *ci))
      {
        Cell_ijkl cell(*ci);
        output_cells.insert(cell);

        typename Star::FT distortion = compute_distortion_t(stars, *ci);
        colors[cell] = (distortion - 1.) * step_n;
      }
      else
      {
        nb_inconsistent_stars++;
        break;
      }
    }
  }

  fx << "Tetrahedra\n";
  fx << output_cells.size() << std::endl;
  typename std::set<Cell_ijkl>::iterator  cit = output_cells.begin();
  typename std::set<Cell_ijkl>::iterator  citend = output_cells.end();
  for (; cit != citend; cit++)
  {
    Cell_ijkl c = *cit;
    int n0 = c.vertices()[0];
    int n1 = c.vertices()[1];
    int n2 = c.vertices()[2];
    int n3 = c.vertices()[3];
    typename Star::Point_3 p0 = stars[n0]->center_point();
    typename Star::Point_3 p1 = stars[n1]->center_point();
    typename Star::Point_3 p2 = stars[n2]->center_point();
    typename Star::Point_3 p3 = stars[n3]->center_point();
    typename Star::Traits::Compute_volume_3 o;
    if(!positive_vol || o(p0, p1, p2, p3) > 0)
      fx << (n0+1) << " " << (n1+1) << " " << (n2+1) << " " << (n3+1) << " " << colors[c] << std::endl;
    else
      fx << (n1+1) << " " << (n0+1) << " " << (n2+1) << " " << (n3+1) << " " << colors[c] << std::endl;
  }
  if(nb_inconsistent_stars > 0)
    std::cout << "Warning Medit: there are " << nb_inconsistent_stars << " inconsistent stars in the ouput mesh.\n";
}

}  // Anisotropic_mesh_3
}  // CGAL

#endif

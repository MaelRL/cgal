#ifndef CGAL_ANISOTROPIC_MESH_2_STAR_SET_OUTPUT_H
#define CGAL_ANISOTROPIC_MESH_2_STAR_SET_OUTPUT_H

#include <CGAL/helpers/combinatorics_helper.h>

#include <fstream>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_2
{

template<typename Starset>
typename Starset::FT
compute_distortion_t(const Starset& stars,
                     const typename Starset::Face_handle& fh)
{
  typename Starset::FT distortion = 1.0;
  for(int i = 0; i < 3; i++)
  {
    for(int j=i+1; j<3; j++)
      distortion = (std::max)(distortion,
                              stars[fh->vertex(i)->info()]->metric().compute_distortion(
                                        stars[fh->vertex(j)->info()]->metric()));
  }
  return distortion;
}

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
    typename Starset::Vertex_handle_handle nsi = star_i->begin_neighboring_vertices();
    typename Starset::Vertex_handle_handle nsiend = star_i->end_neighboring_vertices();
    fx << nsiend - nsi;
    for(; nsi != nsiend; nsi++)
      fx << " " << (*nsi)->info();
    fx << std::endl;
  }
}

template<typename Starset>
void output_off(const Starset& stars,
                std::ofstream& fx,
                const bool consistent_only = false)
{
  std::cout << "Saving as .off... @ " << stars.size() << std::endl;
  std::vector<typename Starset::Point_2> points;
  Facet_ijk_unordered_set output_faces;
  unsigned int nb_inconsistent_stars = 0;

  typename Starset::const_iterator it = stars.begin();
  typename Starset::const_iterator itend = stars.end();
  for (; it != itend; ++it)
  {
    typename Starset::Star_handle star = *it;
    points.push_back(star->center_point());

    typename Starset::Face_handle_handle fit = star->finite_incident_faces_begin();
    typename Starset::Face_handle_handle fend = star->finite_incident_faces_end();
    for(; fit!=fend; ++fit)
    {
      typename Starset::Face_handle fh = *fit;
//      if(!star->is_inside(fh))
//        continue;
      if(!consistent_only || stars.is_consistent(fh))
      {
        Facet_ijk face(fh);
        output_faces.insert(face);
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
  fx << points.size() << " " << output_faces.size() << " 0" << std::endl;
  for(unsigned int i = 0; i < points.size(); i++)
    fx << points[i] << " 0" << std::endl;

  typename Facet_ijk_unordered_set::iterator fit = output_faces.begin();
  typename Facet_ijk_unordered_set::iterator fitend = output_faces.end();
  for (; fit != fitend; fit++)
  {
    fx << "3  " << fit->vertices()[0]
       << " "   << fit->vertices()[1]
       << " "   << fit->vertices()[2] << std::endl;
  }
}

template<typename Starset>
void output_medit(const Starset& stars,
                  std::ofstream& fx,
                  const bool consistent_only = false,
                  const std::set<Facet_ijk>& ghost_faces = std::set<Facet_ijk>())
{
  std::cout << "Saving medit... @ " << stars.size() << std::endl;
  unsigned int nb_inconsistent_stars = 0;
  std::map<Facet_ijk, bool> color;

  fx << "MeshVersionFormatted 1\n";
  fx << "Dimension 2\n";

  fx << "Vertices\n";
  fx << stars.size() << std::endl;
  for(std::size_t i = 0; i < stars.size(); i++)
    fx << stars[i]->center_point() << " " << (i+1) << std::endl; //indices start at 1 in Medit

  Facet_ijk_unordered_set output_faces;
  typename Starset::const_iterator sit = stars.begin();
  typename Starset::const_iterator sitend = stars.end();
  for(; sit != sitend; sit++)
  {
    typename Starset::Star_handle star = *sit;

    typename Starset::Face_handle_handle fit = star->finite_incident_faces_begin();
    typename Starset::Face_handle_handle fend = star->finite_incident_faces_end();
    for(; fit!=fend; ++fit)
    {
      typename Starset::Face_handle fh = *fit;
      //if(!star->is_inside(fh))
      //  continue;

      if(ghost_faces.find(Facet_ijk(fh)) != ghost_faces.end() )
        continue;

      if(!consistent_only || stars.is_consistent(fh))
      {
        Facet_ijk face(fh);
        output_faces.insert(face);
        color[face] = stars.is_consistent(fh);
      }
      else
      {
        nb_inconsistent_stars++;
        break;
      }
    }
  }

  fx << "Triangles\n";
  fx << output_faces.size() << std::endl;
  typename Facet_ijk_unordered_set::iterator fit = output_faces.begin();
  typename Facet_ijk_unordered_set::iterator fend = output_faces.end();
  for (; fit!=fend; fit++)
  {
    Facet_ijk f(*fit);
    int n0 = f.vertices()[0];
    int n1 = f.vertices()[1];
    int n2 = f.vertices()[2];
    fx << (n0+1) << " " << (n1+1) << " " << (n2+1);

    if(true) // tmp
    {
      if(color[f])
        fx << " 2" << std::endl;
      else
        fx << " 1" << std::endl;
    }
    else
      fx << " 1" << std::endl;
  }
  fx << "End" << std::endl;

  if(nb_inconsistent_stars > 0)
    std::cout << "Warning Medit: there are " << nb_inconsistent_stars << " inconsistent stars in the ouput mesh.\n";
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

}  // Anisotropic_mesh_2
}  // CGAL

#endif

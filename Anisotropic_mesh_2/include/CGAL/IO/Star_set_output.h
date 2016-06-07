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

void mesh_to_dump(const char* filename)
{
  std::cout << "mesh to dump on file : " << filename << std::endl;

  std::string input_filename = filename;
  std::string stem = input_filename.substr(0, input_filename.find_last_of('.')); ;
  std::string extension = input_filename.substr(input_filename.find_last_of('.'));

  if(extension != ".mesh")
  {
    std::cout << "must be a .mesh in input" << std::endl;
    return;
  }

  std::ifstream in(input_filename.c_str());
  std::ofstream out((stem + ".dump").c_str());

  if(!in)
    std::cout << "couldn't open the file" << std::endl;

  std::string useless_str;
  int useless_int, nv, nt;

  in >> useless_str >> useless_int; // MeshVersionFormatted 1
  in >> useless_str >> useless_int; // Dimension 2
  CGAL_assertion(useless_int == 2);

  in >> useless_str >> nv; // Vertices
  std::cout << nv << " vertices" << std::endl;

  out << nv << std::endl;
  for(std::size_t i=0; i<nv; ++i)
  {
    double x, y;
    in >> x >> y >> useless_int;
    out << x << " " << y << '\n';
  }

  in >> useless_str >> nt;
  std::cout << nt << " triangles" << std::endl;

  std::vector<boost::unordered_set<std::size_t> > neighborhoods(nv);
  for(std::size_t i=0; i<nt; ++i)
  {
    std::size_t t1, t2, t3;
    in >> t1 >> t2 >> t3 >> useless_int;
    --t1; --t2; --t3; // thanks, medit

    neighborhoods[t1].insert(t2);
    neighborhoods[t1].insert(t3);

    neighborhoods[t2].insert(t1);
    neighborhoods[t2].insert(t3);

    neighborhoods[t3].insert(t1);
    neighborhoods[t3].insert(t2);
  }

  for(std::size_t i=0; i<nv; ++i)
  {
    const boost::unordered_set<std::size_t>& neighborhood = neighborhoods[i];
    out << neighborhood.size() << " ";
    boost::unordered_set<std::size_t>::const_iterator it = neighborhood.begin();
    boost::unordered_set<std::size_t>::const_iterator end = neighborhood.end();
    for(; it!=end; ++it)
    {
      out << *it << " ";
    }
    out << '\n';
  }
  out << std::endl;
}

template<typename Starset>
void dump(const Starset& stars,
          std::ofstream& fx)
{
  std::size_t ns = stars.size();
  fx << ns << '\n';
  for(std::size_t i=0; i<ns; ++i)
    fx << stars[i]->center_point() << '\n';

  for(std::size_t i=0; i<ns; ++i)
  {
    typename Starset::Star_handle star_i = stars[i];
    typename Starset::Vertex_handle_handle vit = star_i->finite_adjacent_vertices_begin();
    typename Starset::Vertex_handle_handle vend = star_i->finite_adjacent_vertices_end();
    fx << vend - vit;
    for(; vit != vend; vit++)
      fx << " " << (*vit)->info();
    fx << '\n';
  }
  fx << std::endl;
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

  fx << "OFF" << '\n';
  fx << points.size() << " " << output_faces.size() << " 0" << '\n';
  for(unsigned int i = 0; i < points.size(); i++)
    fx << points[i] << " 0" << '\n';

  typename Facet_ijk_unordered_set::iterator fit = output_faces.begin();
  typename Facet_ijk_unordered_set::iterator fitend = output_faces.end();
  for (; fit != fitend; fit++)
  {
    fx << "3  " << fit->vertices()[0]
       << " "   << fit->vertices()[1]
       << " "   << fit->vertices()[2] << '\n';
  }
  fx << std::endl;
}

template<typename Starset>
void output_medit(const Starset& stars,
                  std::ofstream& fx,
                  const bool consistent_only = false,
                  const std::set<Facet_ijk>& ghost_faces = std::set<Facet_ijk>())
{
  std::cout << "Saving medit... @ " << stars.size() << '\n';
  unsigned int nb_inconsistent_stars = 0;
  std::map<Facet_ijk, bool> colors;
  int step_n = 100;

  fx << "MeshVersionFormatted 1\n";
  fx << "Dimension 2\n";

  fx << "Vertices\n";
  fx << stars.size() << '\n';
  for(std::size_t i = 0; i < stars.size(); i++)
    fx << stars[i]->center_point() << " " << (i+1) << '\n'; //indices start at 1 in Medit

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
//      if(!star->is_inside(fh))
//        continue;

      if(ghost_faces.find(Facet_ijk(fh)) != ghost_faces.end() )
        continue;

      if(!consistent_only || stars.is_consistent(fh))
      {
        Facet_ijk face(fh);
        output_faces.insert(face);

        typename Starset::FT distortion = compute_distortion_t(stars, fh);
        colors[face] = (distortion - 1.) * step_n;
      }
      else
      {
        nb_inconsistent_stars++;
        break;
      }
    }
  }

  fx << "Triangles\n";
  fx << output_faces.size() << '\n';
  typename Facet_ijk_unordered_set::iterator fit = output_faces.begin();
  typename Facet_ijk_unordered_set::iterator fend = output_faces.end();
  for (; fit!=fend; fit++)
  {
    Facet_ijk f(*fit);
    int n0 = f.vertices()[0];
    int n1 = f.vertices()[1];
    int n2 = f.vertices()[2];
    fx << (n0+1) << " " << (n1+1) << " " << (n2+1) << " ";

#if 1
    if(consistent_only)
    {
      // if using consistent only, the tet is necessarily consistent
      // and it's pointless to try to check for inconsistencies again
      fx << " 1" << '\n';
    }
    else
    {
      typename Starset::Face_handle dummy;
      bool is_consistent = (stars.get_star(n0)->has_face(f, dummy) &&
                            stars.get_star(n1)->has_face(f, dummy) &&
                            stars.get_star(n2)->has_face(f, dummy));
      fx << (is_consistent ? "0":"1") << '\n';
    }
#else
    fx << colors[c] << '\n';
#endif
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

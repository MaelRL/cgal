#ifndef CGAL_ANISOTROPIC_MESH_3_STAR_SET_IO_H
#define CGAL_ANISOTROPIC_MESH_3_STAR_SET_IO_H

#include <CGAL/helpers/combinatorics_helper.h>

#include <fstream>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename Starset>
typename Starset::FT
compute_distortion_t(const Starset& stars,
                     const typename Starset::Cell_handle& c)
{
  typename Starset::FT distortion = 1.0;
  for(int i = 0; i < 4; i++)
  {
    for(int j = i + 1; j < 4; j++)
      distortion = (std::max)(distortion,
                              stars[c->vertex(i)->info()]->metric().compute_distortion(
                     stars[c->vertex(j)->info()]->metric()));
  }
  return distortion;
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
    std::set<std::size_t> neighbors;
    typename Starset::Star_handle star_i = stars[i];
    typename Starset::Vertex_handle_handle vit = star_i->finite_adjacent_vertices_begin();
    typename Starset::Vertex_handle_handle vend = star_i->finite_adjacent_vertices_end();
    for(; vit!=vend; vit++)
      neighbors.insert((*vit)->info());

    fx << neighbors.size();

    std::set<std::size_t>::iterator it = neighbors.begin();
    std::set<std::size_t>::iterator end = neighbors.end();
    for(; it!=end; ++it)
      fx << " " << *it;
    fx << '\n';
  }
  fx << std::endl;
}

template<typename Starset>
void read_dump(Starset& ss,
               std::string filename)
{
  std::cout << "Reading dump..." << std::endl;
  std::ifstream in(filename.c_str());
  ss.clear();

  std::size_t stars_n, v_n, id;
  typename Starset::FT x, y, z;

  in >> stars_n;

  for(std::size_t i=0; i<stars_n; ++i)
  {
    in >> x >> y >> z;
    typename Starset::Point_3 p(x, y, z);

    const typename Starset::Criteria * c = ss.criteria();
    const typename Starset::Constrain_surface * d = ss.constrain_surface();

    typename Starset::Star_handle star = new typename Starset::Star(c, d);
    const typename Starset::Metric& m_p = ss.metric_field()->compute_metric(p);
    star->reset(p, i, m_p);
    ss.push_back(star);
  }

  for(std::size_t i=0; i<stars_n; ++i)
  {
    in >> v_n;
    typename Starset::Star_handle star_i = ss[i];
    for(std::size_t j=0; j<v_n; ++j)
    {
      in >> id;
      typename Starset::Star_handle star_j = ss[id];
      star_i->insert_to_star(star_j->center_point(), id, false);
    }
    std::cout << "built star: " << i << std::endl;
  }

  std::cout << ss.size() << " stars from dump" << std::endl;
}

template<typename Starset>
void output_surface_star_off(const Starset& stars,
                             std::ofstream& fx,
                             typename Starset::Index i)
{
  typename Starset::Star_handle star = stars[i];
  std::map<typename Starset::Index, int> match_indices;
  int off_index = 0;
  std::vector<typename Starset::Point_3> points;
  Facet_ijk_unordered_set output_facets;

  typename Starset::Facet_set_iterator rfi = star->restricted_facets_begin();
  typename Starset::Facet_set_iterator rfend = star->restricted_facets_end();
  for(; rfi!=rfend; ++rfi)
  {
    boost::array<typename Starset::Index, 3> ids;
    for(int i=0; i<3; i++)
    {
      typename Starset::Vertex_handle vi = rfi->first->vertex((rfi->second+i+1)%4);
      ids[i] = vi->info();
      std::pair<typename std::map<typename Starset::Index, int>::iterator, bool> inserted;
      inserted = match_indices.insert(std::make_pair(ids[i], off_index));
      if(inserted.second)
      {
        points.push_back(star->metric().inverse_transform(vi->point()));
        off_index++;
      }
    }
    output_facets.insert(Facet_ijk(ids[0], ids[1], ids[2]));
  }

  fx << "OFF" << '\n';
  fx << points.size() << " " << output_facets.size() << " 0" << '\n';
  for(unsigned int i = 0; i < points.size(); i++)
    fx << points[i] << '\n';

  typename Facet_ijk_unordered_set::iterator fit = output_facets.begin();
  typename Facet_ijk_unordered_set::iterator fitend = output_facets.end();
  for (; fit != fitend; fit++)
  {
    fx << "3  " << match_indices[fit->vertices()[0]]
        << " "   << match_indices[fit->vertices()[1]]
        << " "   << match_indices[fit->vertices()[2]] << '\n';
  }
  fx << std::endl;
}

template<typename Star>
void output_star_off(const Star* star,
                     std::ofstream& fx)
{
  std::map<typename Star::Index, int> match_indices;
  int off_index = 0;
  std::vector<typename Star::Point_3> points;
  Facet_ijk_unordered_set output_facets;

  typename Star::Cell_handle_handle ci = star->finite_star_cells_begin();
  typename Star::Cell_handle_handle ciend = star->finite_star_cells_end();
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

    output_facets.insert(Facet_ijk(ids[0], ids[1], ids[2]));
    output_facets.insert(Facet_ijk(ids[0], ids[2], ids[3]));
    output_facets.insert(Facet_ijk(ids[0], ids[1], ids[3]));
    output_facets.insert(Facet_ijk(ids[3], ids[1], ids[2]));
  }

  fx << "OFF" << '\n';
  fx << points.size() << " " << output_facets.size() << " 0" << '\n';
  for(unsigned int i = 0; i < points.size(); i++)
    fx << points[i] << '\n';

  typename Facet_ijk_unordered_set::iterator fit = output_facets.begin();
  typename Facet_ijk_unordered_set::iterator fitend = output_facets.end();
  for (; fit != fitend; fit++)
  {
    fx << "3  " << match_indices[fit->vertices()[0] ]
        << " "   << match_indices[fit->vertices()[1] ]
        << " "   << match_indices[fit->vertices()[2] ] << '\n';
  }
  fx << std::endl;
}

template<typename Starset>
void output_star_off(const Starset& stars,
                     std::ofstream& fx,
                     typename Starset::Index i)
{
  typename Starset::Star_handle star = stars[i];
  return output_star_off(star, fx);
}

template<typename Starset>
void output_off(const Starset& stars,
                std::ofstream& fx,
                const bool consistent_only = true)
{
  std::cout << "Saving as .off..." << std::endl;
  std::vector<typename Starset::Point_3> points;
  Facet_ijk_unordered_set output_facets;
  Cell_ijkl_unordered_set output_cells;
  unsigned int nb_inconsistent_stars = 0;

  typename Starset::const_iterator it = stars.begin();
  typename Starset::const_iterator itend = stars.end();
  for (; it != itend; ++it)
  {
    typename Starset::Star_handle star = *it;
    points.push_back(star->center_point());

    typename Starset::Cell_handle_handle ci = star->finite_star_cells_begin();
    typename Starset::Cell_handle_handle ciend = star->finite_star_cells_end();
    for (; ci != ciend; ++ci)
    {
      typename Starset::Cell_handle c = *ci;
      if(!star->is_inside(*ci))
        continue;

      if(!consistent_only || stars.is_consistent(c))
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

  fx << "OFF" << '\n';
  fx << points.size() << " " << output_facets.size() << " " << output_cells.size() << '\n';
  for(unsigned int i = 0; i < points.size(); i++)
    fx << points[i] << '\n';

  typename Facet_ijk_unordered_set::iterator fit = output_facets.begin();
  typename Facet_ijk_unordered_set::iterator fitend = output_facets.end();
  for (; fit != fitend; fit++)
  {
    fx << "3  " << fit->vertices()[0]
       << " "   << fit->vertices()[1]
       << " "   << fit->vertices()[2] << '\n';
  }

  typename Cell_ijkl_unordered_set::iterator  cit = output_cells.begin();
  typename Cell_ijkl_unordered_set::iterator  citend = output_cells.end();
  for (; cit != citend; cit++)
  {
    fx << "4  " << cit->vertices()[0]
       << " "   << cit->vertices()[1]
       << " "   << cit->vertices()[2]
       << " "   << cit->vertices()[3] << '\n';
  }
  fx << std::endl;
}

template<typename Starset>
void output_surface_off(const Starset& stars,
                        std::ofstream& fx,
                        const bool consistent_only = true)
{
  std::cout << "Saving surface off..." << std::endl;
  std::vector<typename Starset::Point_3> points;
  Facet_ijk_unordered_set output_facets;
  unsigned int nb_inconsistent_stars = 0;

  typename Starset::const_iterator si = stars.begin();
  typename Starset::const_iterator siend = stars.end();
  for (; si != siend; ++si)
  {
    typename Starset::Star_handle star = *si;
    points.push_back(star->center_point());

    typename Starset::Facet_set_iterator fi = star->restricted_facets_begin();
    typename Starset::Facet_set_iterator fiend = star->restricted_facets_end();
    for (; fi != fiend; ++fi)
    {
      typename Starset::Point_3 cc;
      star->compute_dual_intersection(*fi, cc);

      if(!consistent_only || stars.is_consistent(*fi))
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

  fx << "OFF" << '\n';
  fx << points.size() << " " << output_facets.size() << " " << 0 << '\n';
  for(unsigned int i = 0; i < points.size(); i++)
    fx << points[i] << '\n';

  typename Facet_ijk_unordered_set::iterator fit = output_facets.begin();
  typename Facet_ijk_unordered_set::iterator fitend = output_facets.end();
  for (; fit != fitend; fit++)
  {
    fx << "3  " << fit->vertices()[0]
       << " "   << fit->vertices()[1]
       << " "   << fit->vertices()[2] << '\n';
  }
  fx << std::endl;
}

template<typename Starset>
void output_surface_medit(const Starset& stars,
                          std::ofstream& fx,
                          const bool positive_vol = false,
                          const bool fake_3D = false) // link all facets to a dummy point
{
  typename Starset::Point_3 anchor_point(0., 0., 0.);

  std::cout << "Saving surface medit..." << std::endl;
  fx << "MeshVersionFormatted 1\n";
  fx << "Dimension 3" << '\n';

  fx << "Vertices" << '\n';
  fx << (stars.size() + fake_3D) << '\n';
  for (std::size_t i=0; i<stars.size(); i++)
    fx << stars[i]->center_point() << " " << (i + 1) << '\n';
  if(fake_3D)
    fx << anchor_point << (stars.size() + 1) << '\n';

  Facet_ijk_unordered_set facets;
  for(std::size_t i=0; i<stars.size(); i++)
  {
    typename Starset::Star_handle s = stars[i];
    typename Starset::Facet_set_iterator fit = s->restricted_facets_begin();
    typename Starset::Facet_set_iterator fend = s->restricted_facets_end();
    for(; fit != fend; ++fit)
      facets.insert(Facet_ijk(*fit));
  }

  fx << "Triangles" << '\n';
  fx << (facets.size()) << '\n';
  typename Facet_ijk_unordered_set::iterator fit = facets.begin();
  typename Facet_ijk_unordered_set::iterator fitend = facets.end();
  for(; fit!=fitend; fit++)
  {
    const Facet_ijk& f = *fit;
    int n0 = f.vertices()[0];
    int n1 = f.vertices()[1];
    int n2 = f.vertices()[2];

    fx << (n0+1) << " " << (n1+1) << " " << (n2+1) << " ";

#if 1
    typename Starset::Facet useless;
    bool is_consistent = (stars.get_star(n0)->has_facet(f, useless) &&
                          stars.get_star(n1)->has_facet(f, useless) &&
                          stars.get_star(n2)->has_facet(f, useless));
    fx << is_consistent << '\n';
#else
    typename Starset::Point_3 p0 = stars[n0]->center_point();
    typename Starset::Point_3 p1 = stars[n1]->center_point();
    typename Starset::Point_3 p2 = stars[n2]->center_point();
    typename Starset::Traits::Compute_volume_3 o;
    if(!positive_vol || o(p0, p1, p2, p3) > 0) //color[c] creates the artifacts. fixme. todo
      fx << colors[c] << '\n';
    else
      fx << " 1" << '\n';
#endif
  }

  if(fake_3D)
  {
    fx << "Tetrahedra" << '\n';
    fx << (facets.size()) << '\n';
    fit = facets.begin();
    for(; fit!=fitend; fit++)
    {
      const Facet_ijk& f = *fit;
      int n0 = f.vertices()[0];
      int n1 = f.vertices()[1];
      int n2 = f.vertices()[2];
      int n3 = stars.size(); //anchor_id
      typename Starset::Point_3 p0 = stars[n0]->center_point();
      typename Starset::Point_3 p1 = stars[n1]->center_point();
      typename Starset::Point_3 p2 = stars[n2]->center_point();
      typename Starset::Traits::Compute_volume_3 o;
      if(!positive_vol || o(p0, p1, p2, anchor_point) > 0)
        fx << (n0+1) << " " << (n1+1) << " " << (n2+1) << " " << (n3+1) << " 1" << '\n';
      else
        fx << (n1+1) << " " << (n0+1) << " " << (n2+1) << " " << (n3+1) << " 1" << '\n';
    }
  }
  fx << "End" << std::endl;
}

template<typename Starset>
void output_medit(const Starset& stars,
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
  fx << stars.size() << '\n';
  for(std::size_t i = 0; i < stars.size(); i++)
    fx << stars[i]->center_point() << " " << (i+1) << '\n'; //indices start at 1 in Medit

  Facet_ijk_unordered_set output_facets;
  Cell_ijkl_unordered_set output_cells;

  typename Starset::const_iterator sit = stars.begin();
  typename Starset::const_iterator sitend = stars.end();
  for(; sit != sitend; sit++)
  {
    typename Starset::Star_handle star = *sit;
    typename Starset::Cell_handle_handle ci = star->finite_star_cells_begin();
    typename Starset::Cell_handle_handle ciend = star->finite_star_cells_end();
    for(; ci != ciend; ci++)
    {
      typename Starset::Cell_handle ch = *ci;
      if(!star->is_inside(ch))
        continue;

      if(!consistent_only || stars.is_consistent(ch))
      {
        Cell_ijkl c(ch);
        output_cells.insert(c);
        output_facets.insert(Facet_ijk(ch->vertex(0)->info(), ch->vertex(1)->info(),
                                       ch->vertex(2)->info()));
        output_facets.insert(Facet_ijk(ch->vertex(0)->info(), ch->vertex(2)->info(),
                                       ch->vertex(3)->info()));
        output_facets.insert(Facet_ijk(ch->vertex(0)->info(), ch->vertex(1)->info(),
                                       ch->vertex(3)->info()));
        output_facets.insert(Facet_ijk(ch->vertex(3)->info(), ch->vertex(1)->info(),
                                       ch->vertex(2)->info()));

        typename Starset::FT distortion = compute_distortion_t(stars, ch);
        colors[c] = (distortion - 1.) * step_n;
      }
      else
      {
        nb_inconsistent_stars++;
        break;
      }
    }
  }

  std::cout << output_facets.size() << " facets & " << output_cells.size() << " cells" << std::endl;

  if(nb_inconsistent_stars > 0)
  {
    std::cout << "Warning Medit: there are "
              << nb_inconsistent_stars << " inconsistent stars in the ouput mesh." << std::endl;
  }

/*
  fx << "Triangles\n";
  fx << output_facets.size() << '\n';
  typename Facet_ijk_unordered_set::iterator fit = output_facets.begin();
  typename Facet_ijk_unordered_set::iterator fitend = output_facets.end();
  for (; fit != fitend; fit++)
  {
    const Facet_ijk& f = *fit;
    int n0 = f.vertices()[0];
    int n1 = f.vertices()[1];
    int n2 = f.vertices()[2];
    fx << (n0+1) << " " << (n1+1) << " " << (n2+1) << " ";

    if(consistent_only)
    {
      // if using consistent only, the facet is necessarily consistent
      // and it's pointless to try to check for inconsistencies again
      fx << " 1" << '\n';
    }
    else
    {
      typename Starset::Facet dummy;
      bool is_consistent = (stars.get_star(n0)->has_facet(f, dummy) &&
                            stars.get_star(n1)->has_facet(f, dummy) &&
                            stars.get_star(n2)->has_facet(f, dummy));
      fx << (is_consistent?"1":"2") << '\n';
    }
  }
*/

  fx << "Tetrahedra\n";
  fx << output_cells.size() << '\n';
  typename Cell_ijkl_unordered_set::iterator  cit = output_cells.begin();
  typename Cell_ijkl_unordered_set::iterator  citend = output_cells.end();
  for (; cit != citend; cit++)
  {
    const Cell_ijkl& c = *cit;
    int n0 = c.vertices()[0];
    int n1 = c.vertices()[1];
    int n2 = c.vertices()[2];
    int n3 = c.vertices()[3];

    fx << (n0+1) << " " << (n1+1) << " " << (n2+1) << " " << (n3+1) << " ";

#if 1
    if(consistent_only)
    {
      // if using consistent only, the tet is necessarily consistent
      // and it's pointless to try to check for inconsistencies again
      fx << " 1" << '\n';
    }
    else
    {
      typename Starset::Cell_handle dummy;
      bool is_consistent = (stars.get_star(n0)->has_cell(c, dummy) &&
                            stars.get_star(n1)->has_cell(c, dummy) &&
                            stars.get_star(n2)->has_cell(c, dummy) &&
                            stars.get_star(n3)->has_cell(c, dummy) );
      fx << (is_consistent ? "1":"5") << '\n';
    }
#else
    typename Starset::Point_3 p0 = stars[n0]->center_point();
    typename Starset::Point_3 p1 = stars[n1]->center_point();
    typename Starset::Point_3 p2 = stars[n2]->center_point();
    typename Starset::Point_3 p3 = stars[n3]->center_point();
    typename Starset::Traits::Compute_volume_3 o;
    if(!positive_vol || o(p0, p1, p2, p3) > 0) //color[c] creates the artifacts. fixme. todo
      fx << colors[c] << '\n';
    else
      fx << " 1" << '\n';
#endif
  }

  fx << "End" << std::endl;
}

template<typename Starset>
void output_surface_voronoi(const Starset& stars,
                            std::ofstream& fx)
{
  typedef typename Starset::Point_3      Point;

  std::map<Facet_ijk, int> facets;
  std::map<std::pair<int,int>, std::vector<int> > edges;
  std::vector<Point> points;

  std::size_t count = 1; // 1 for first vertex of medit
  typename Starset::const_iterator sit = stars.begin();
  typename Starset::const_iterator sitend = stars.end();
  for(; sit != sitend; sit++)
  {
    typename Starset::Star_handle si = *sit;
    typename Starset::Facet_set_iterator fit = si->restricted_facets_begin();
    typename Starset::Facet_set_iterator fend = si->restricted_facets_end();
    for(; fit != fend; ++fit)
    {
      Point p;
      si->compute_dual_intersection(*fit, p);
      points.push_back(p);
      facets[Facet_ijk(*fit)] = count++;
    }
  }

  typename std::map<Facet_ijk, int>::iterator it = facets.begin();
  typename std::map<Facet_ijk, int>::iterator itend = facets.end();
  for(;it!=itend;++it)
  {
    Facet_ijk f = it->first;
    for(int i=0; i<3; ++i)
    {
      int j = (i+1)%3;
      int p1 = f[i];
      int p2 = f[j];
      if(p2 < p1)
      {
        p1 = f[j];
        p2 = f[i];
      }
      edges[std::make_pair(p1, p2)].push_back(facets[f]);
    }
  }

  fx << "MeshVersionFormatted 1" << '\n';
  fx << "Dimension 3" << '\n';
  fx << "Vertices" << '\n';
  fx << points.size() << '\n';
  for(std::size_t i=0; i<points.size(); ++i)
    fx << points[i].x() << " " << points[i].y() << " " << points[i].z() << " 0" << '\n';

  fx << "Edges" << '\n';
  fx << edges.size() << '\n';
  typename std::map<std::pair<int,int>, std::vector<int> >::iterator eit = edges.begin();
  typename std::map<std::pair<int,int>, std::vector<int> >::iterator eend = edges.end();
  for(;eit!=eend;++eit)
  {
    std::vector<int> edge = eit->second;
    assert(edge.size() == 2);
    fx << edge.front() << " " << edge.back() << " 0" << '\n';
  }

  fx << "End" << std::endl;
}

template<typename Starset>
void output_surface_off(const Starset& stars,
                        const char* f,
                        const bool consistent_only = true)
{
  std::ofstream out(f);
  out.precision(17);
  return output_surface_off(stars, out, consistent_only);
}

template<typename Starset>
void output_off(const Starset &stars,
                const char* f,
                const bool consistent_only = true)
{
  std::ofstream out(f);
  out.precision(17);
  return output_off(stars, out, consistent_only);
}

template<typename Starset>
void output_surface_medit(const Starset &stars,
                          const char* f,
                          const bool positive_vol = false,
                          const bool fake_3D = false)
{
  std::ofstream out(f);
  out.precision(17);
  return output_surface_medit(stars, out, positive_vol, fake_3D);
}

template<typename Starset>
void output_medit(const Starset &stars,
                  const char* f,
                  const bool consistent_only = true,
                  const bool positive_vol = false)
{
  std::ofstream out(f);
  out.precision(17);
  return output_medit(stars, out, consistent_only, positive_vol);
}

template<typename Starset>
void output_surface_voronoi(const Starset &stars,
                            const char* f)
{
  std::ofstream out(f);
  out.precision(17);
  return output_surface_voronoi(stars, out);
}

}  // Anisotropic_mesh_3
}  // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_STAR_SET_IO_H

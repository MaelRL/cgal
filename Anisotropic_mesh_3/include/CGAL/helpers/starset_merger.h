#ifndef STARSET_MERGER_H
#define STARSET_MERGER_H

#define BUILD_BOUNDARIES_INDEPENDENTLY

#include <CGAL/Starset.h>
#include <CGAL/IO/Star_set_IO.h>

#include <CGAL/aabb_tree/aabb_tree_bbox.h>
#include <CGAL/aabb_tree/aabb_tree_bbox_primitive.h>
#include <CGAL/kd_tree/Kd_tree_for_star_set.h>

#include <boost/bimap/bimap.hpp>

#include <cstddef>
#include <fstream>
#include <map>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

// Ideally this eventually uses some partitioner that adapts it to
// a metric function to obtain n (proc number) sub starsets who
// are eventually stitched together.


// the very particular case
template<typename Starset>
void merge_8_small_cubes(Starset& ss)
{
  std::string base_location = "/home/mrouxell/Data/Dump/cube_parts_shock/starset_cube_";
  // kind of a naive implementation

  typedef typename Starset::FT                       FT;
  typedef typename Starset::Point_3                  Point_3;

  typedef std::map<Point_3, std::size_t >            GMap;
  typedef typename GMap::iterator                    GMap_it;

  // first vector is the different parts, second vector is the number of vertices
  // in that loc part
  typedef std::vector<std::vector<std::size_t> >     CVector;

#ifdef BUILD_BOUNDARIES_INDEPENDENTLY
  Point_3 center(1.5, 1.5, 1.5);
  FT eps = 0.1; // boundaries are defined by points with x, y, or z \in [-eps, eps]
#endif

  GMap glob_ids;
  CVector cor_vectors(8);
  std::size_t stars_to_fully_build_count = 0;

  std::size_t curr_glob_id = 0;
  for(std::size_t i=0; i<8; ++i)
  {
    std::ostringstream in;
    in << base_location << i+1 << "_dump.txt";
    std::ifstream c_in(in.str().c_str());
    CGAL_precondition(c_in);

    std::size_t loc_stars_n;
    c_in >> loc_stars_n;

    std::vector<std::size_t>& cor_vector = cor_vectors[i];
    cor_vector.resize(loc_stars_n);

    for(std::size_t j=0; j<loc_stars_n; ++j)
    {
      FT x, y, z;
      c_in >> x >> y >> z;
      typename Starset::Point_3 p(x, y, z);

#ifdef BUILD_BOUNDARIES_INDEPENDENTLY
      if((x >= center.x()-eps && x <= center.x()+eps) ||
         (y >= center.y()-eps && y <= center.y()+eps) ||
         (z >= center.z()-eps && z <= center.z()+eps))
      {
        // need to say that there's no global id currently assigned
        // to the local element (i, j)
        cor_vector[j] = -1;
        stars_to_fully_build_count++;
        std::cout << "Filtered: " << p << std::endl;
        continue;
      }
#endif

      // handle local and global ids
      std::pair<GMap_it, bool> is_insert_successful = glob_ids.insert(std::make_pair(p, curr_glob_id));
      std::size_t glob_id; // corresponding glob_id

      if(is_insert_successful.second)
        glob_id = curr_glob_id++;
      else
        glob_id = (is_insert_successful.first)->second;

      cor_vector[j] = glob_id; // insert in the ([part,loc] = glob) map

      // create star
      typename Starset::Star_handle star = new typename Starset::Star(ss.criteria(),
                                                                      ss.constrain_surface());
      star->reset(p, glob_id, ss.metric_field()->compute_metric(p));
      ss.push_back(star);
      std::cout << "Created star: " << glob_id << " at " << p << std::endl;
    }
  }

  // need to build ALL the global ids first due to tangent (possibly overlapping) stars
  // note that this isn't really optimal atm, it would be cleaner to build the star set
  // except for boundaries and then insert boundary points... Let's see.

  for(std::size_t i=0; i<8; ++i)
  {
    std::ostringstream in;
    in << base_location << i+1 << "_dump.txt";
    std::ifstream c_in(in.str().c_str());
    CGAL_precondition(c_in);

    std::size_t loc_stars_n, v_n, neigh_id;
    c_in >> loc_stars_n;

    for(std::size_t j=0; j<loc_stars_n; ++j)
    {
      FT x, y, z;
      c_in >> x >> y >> z;
      // no need for those... could probably use std::ignore isntead
    }

    for(std::size_t j=0; j<loc_stars_n; ++j)
    {
      // i is the local id, let's just redefine it to be clear
      std::size_t loc_c_id = j; // c because Center of star
      std::size_t glob_c_id = cor_vectors[i][loc_c_id];

      c_in >> v_n;

#ifdef BUILD_BOUNDARIES_INDEPENDENTLY
      // a point whose star has not been built yet

      if(glob_c_id == static_cast<std::size_t>(-1))
      {
        for(std::size_t k=0; k<v_n; ++k)
          c_in >> neigh_id; // don't care about them but still gotta read the entries
        continue;
      }
#endif

      typename Starset::Star_handle star_j = ss[glob_c_id];
      for(std::size_t k=0; k<v_n; ++k)
      {
        c_in >> neigh_id;
        // id is a local id, must translate it to global coordinates
        std::size_t glob_neigh_id = cor_vectors[i][neigh_id];

        typename Starset::Star_handle star_k = ss[glob_neigh_id];
        star_j->insert_to_star(star_k->center_point(), glob_neigh_id,
                               false /*no conflict check*/);
      }

      std::cout << "Built star: " << glob_c_id << std::endl;
    }
  }

#ifdef BUILD_BOUNDARIES_INDEPENDENTLY
  // it's time to build boundaries!
  // points that have not been included in the starset have id -1 in the
  // correspondence map.
  std::cout << stars_to_fully_build_count << " stars to fully build" << std::endl;

  // build a kd tree of the starset to help...

  for(std::size_t i=0; i<8; ++i)
  {
    std::ostringstream in;
    in << base_location << i+1 << "_dump.txt";
    std::ifstream c_in(in.str().c_str());
    CGAL_precondition(c_in);

    std::size_t loc_stars_n;
    c_in >> loc_stars_n;

    for(std::size_t j=0; j<loc_stars_n; ++j)
    {
      FT x, y, z;
      c_in >> x >> y >> z;

      if(cor_vectors[i][j] != -1)
        continue;

      // get the point
      Point_3 p(x, y, z);

      typename Starset::Star_handle star = new typename Starset::Star(ss.criteria(),
                                                                      ss.constrain_surface());
      star->reset(p, curr_glob_id++, ss.metric_field()->compute_metric(p));
      std::cout << "Created star: " << star->index_in_star_set() << " at " << p << std::endl;

      // fun times, insert ALL the existing stars in that new one...
      for(std::size_t k=0; k<ss.size(); ++k)
      {
        typename Starset::Star_handle star_k = ss[k];
        star->insert_to_star(star_k->center_point(), k, false /*no conflict check*/);

        if(k%10000 == 0)
          std::cout << "at: " << k << " out of " << ss.size() << std::endl;

        // trying to filter a bit that mess
        FT d = CGAL::squared_distance(p, star_k->center_point());
        if(d > (4*eps*eps))
          continue;
      }
      star->clean();
      ss.push_back(star);
      std::cout << "Built star: " << star->index_in_star_set() << std::endl;
    }
  }
#endif

  // just to make sure :
  std::ofstream out("merged.mesh");
  output_medit(ss, out, false/*don't filter inconsistencies*/);

  std::ofstream out_d("dump.txt");
  dump(ss, out_d);
}


} // Anisotropic_mesh_3
} // CGAL

#endif // STARSET_MERGER_H

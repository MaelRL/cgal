#ifndef CGAL_ANISOTROPIC_MESH_2_STARSET_H
#define CGAL_ANISOTROPIC_MESH_2_STARSET_H

#include <CGAL/Stretched_Delaunay_2.h>
#include <CGAL/helpers/metric_helper.h>

#include <Eigen/Dense>

#include <iostream>
#include <fstream>
#include <queue>
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
  typedef typename Star::Segment                        Segment;
  typedef typename Star::Triangle                       Triangle;

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

  template<typename I>
  Star_handle get_star(I i) const { return m_stars[i]; }

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

  // flip consistency
  bool is_flip_consistent(const Face_handle& fh,
                          const bool verbose = false) const
  {
    if(is_consistent(fh, verbose))
    {
      return true;
    }

    // face is inconsistent.
    // check if only 4 pts are involved in the inconsistency...

    std::set<std::size_t> involved_stars;
    involved_stars.insert(fh->vertex(0)->info());
    involved_stars.insert(fh->vertex(1)->info());
    involved_stars.insert(fh->vertex(2)->info());

    bool stars_were_added = true;
    while(stars_were_added)
    {
      stars_were_added = false;
      std::set<std::size_t> involved_stars2 = involved_stars; // disgusting but whatever
      std::size_t old_size = involved_stars.size();

      // search for inconsistent simplices that contain a point in involved_stars
      std::set<std::size_t>::iterator it = involved_stars2.begin();
      for(; it!=involved_stars2.end(); ++it) //looping on 2, inserting in 1
      {
        Star_handle s = get_star(*it);
        typename Star::Face_handle_handle fit = s->finite_incident_faces_begin();
        typename Star::Face_handle_handle fend = s->finite_incident_faces_end();
        for(; fit!=fend; ++fit)
        {
          Face_handle fh2 = *fit;
          if(!is_consistent(fh2))
          {
            if(fh2->vertex(0)->info() != s->index_in_star_set())
            {
              involved_stars.insert(fh2->vertex(0)->info());

              if(old_size != involved_stars.size()) // check that it didn't already exist in involved_stars
                stars_were_added = true;
            }

            if(fh2->vertex(1)->info() != s->index_in_star_set())
            {
              involved_stars.insert(fh2->vertex(1)->info());

              if(old_size != involved_stars.size())
                stars_were_added = true;
            }

            if(fh2->vertex(2)->info() != s->index_in_star_set())
            {
              involved_stars.insert(fh2->vertex(2)->info());

              if(old_size != involved_stars.size())
                stars_were_added = true;
            }
          }
        }
      }

      if(involved_stars.size() > 4)
        return false;
    }

    if(involved_stars.size() == 4)
      std::cout << "WUT WUT --------------------------------------------" << std::endl;

    return involved_stars.size() == 4;
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
    typename Star::Face_handle_handle fit = star->finite_incident_faces_begin();
    typename Star::Face_handle_handle fend = star->finite_incident_faces_end();
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
      typename Star::Face_handle_handle fit = s->finite_incident_faces_begin();
      typename Star::Face_handle_handle fend = s->finite_incident_faces_end();
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

      typename Star::Face_handle_handle fit = star->finite_incident_faces_begin();
      typename Star::Face_handle_handle fend = star->finite_incident_faces_end();
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

  struct DIST_POINT
  {
    FT dist;
    std::size_t id;

    DIST_POINT(FT dist_, std::size_t id_):dist(dist_), id(id_)
    { }
  };


  struct Dist_point_compator
  {
    bool operator()(DIST_POINT left,
                    DIST_POINT right)
    {
      return left.dist < right.dist;
    }

    Dist_point_compator(){}
  };

  bool constrain_test(const std::size_t source, // segment's source
                      const std::size_t target, // segment's target
                      const std::set<Facet_ijk>& ghost_faces,
                      std::vector<std::size_t>& constraints,
                      const bool reject = true)
  {
    Point_2 p = m_stars[source]->center_point();
    Point_2 pi = m_stars[target]->center_point();
    Segment s(p, pi);

    std::cout << "constrain testing with : " << source << " ---> " << target << std::endl;
    std::cout << "constrain testing with : " << p << " [[[[ " << pi << std::endl;

    std::size_t count = 0; // count the number of times it appears as an edge in the existing faces

    for(std::size_t sj=0; sj<m_stars.size(); ++sj)
    {
      Star_handle star_j = m_stars[sj];
      Face_handle_handle fit = star_j->finite_incident_faces_begin();
      Face_handle_handle fend = star_j->finite_incident_faces_end();
      for(; fit!=fend; fit++)
      {
        Face_handle fh = *fit;
        Facet_ijk f(fh);
        if(ghost_faces.find(f) != ghost_faces.end()) // ignoring false faces
          continue;

        if(f.has(source) && f.has(target)) // intersection is an edge
        {
          count++;
          continue;
        }

        Point_2 p0 = star_j->metric().inverse_transform(fh->vertex(0)->point());
        Point_2 p1 = star_j->metric().inverse_transform(fh->vertex(1)->point());
        Point_2 p2 = star_j->metric().inverse_transform(fh->vertex(2)->point());

        Triangle t(p0,p1,p2);
        CGAL::Object obj = CGAL::intersection(s, t);

        Point_2 p_test;
        if(CGAL::assign(p_test, obj) &&
           CGAL::squared_distance(p_test, p)>1e-10 &&
           CGAL::squared_distance(p_test, pi)>1e-10)
        {
          std::cout.precision(20);
          std::cout << "fail p test : "  << p_test << std::endl;
          return false;
        }
        Segment s_test;
        if(CGAL::assign(s_test, obj) &&
           ((CGAL::squared_distance(s_test.source(), p)>1e-10 &&
             CGAL::squared_distance(s_test.source(), pi)>1e-10) ||
            (CGAL::squared_distance(s_test.target(), p)>1e-10 &&
             CGAL::squared_distance(s_test.target(), pi)>1e-10)))
        {
          std::cout.precision(20);
          std::cout << "fail s test : " <<  sj << " " << s_test.source() << " ||| " << s_test.target() << std::endl;
          std::cout << "face points: " << p0 << " ~~ " << p1 << " ~~ " << p2 << std::endl;
          return false;
        }
      }
    }

    CGAL_assertion(count < 3);

    if(reject && count==2) // edge is found in two triangles -> ignore
      return false;

    if(reject && count == 1) // add it as a constraint
      constraints.push_back(target);

    return true;
  }

  bool constrain_test(const Star_handle star,
                      const Face_handle& fh,
                      const std::set<Facet_ijk>& ghost_faces)
  {
    std::cout << "heh" << std::endl;
    std::vector<std::size_t> extremities;
    for(int i=0;i<3;++i)
    {
      if(star->index_in_star_set() != fh->vertex(i)->info())
        extremities.push_back(fh->vertex(i)->info());
    }
    CGAL_assertion(extremities.size() == 2);
    std::vector<std::size_t> useless_list;
    return constrain_test(extremities[0], extremities[1], ghost_faces, useless_list, false);
  }

  void constrain()
  {
    for(std::size_t si=0; si<m_stars.size(); ++si)
      m_stars[si]->reset();

    //pick the star with the highest anisotropy ratio
    std::size_t next_star = 0;
    FT max_ratio = -1e30;
    for(std::size_t si=0; si<m_stars.size(); ++si)
    {
      Star_handle star_i = m_stars[si];
      FT ratio = star_i->metric().get_max_eigenvalue() / star_i->metric().get_min_eigenvalue();
      if(ratio > max_ratio)
      {
        max_ratio = ratio;
        next_star = si;
      }
    }

    //build the stars
    std::set<Facet_ijk> ghost_faces; //facets that do not (really) exist
    typename std::multiset<DIST_POINT, Dist_point_compator> next_stars;
    std::list<std::size_t> queue;
    std::vector<bool> is_star_built(m_stars.size(), false);
    std::size_t stars_left = m_stars.size();
    while(stars_left)
    {
      std::cout << stars_left << " stars left (" << m_stars.size()-stars_left << ")" << std::endl;
      std::cout << "next star: " << next_star << " cp: " << m_stars[next_star]->center_point() << std::endl;

      Star_handle star = m_stars[next_star];

      //compute points available
      std::vector<std::size_t> available_points;
      std::vector<std::size_t> constraints;
      for(std::size_t si=0; si<m_stars.size(); ++si)
      {
        if(si == next_star || is_star_built[si])
          continue;

        if(constrain_test(next_star, si, ghost_faces, constraints, true))
          available_points.push_back(si);
      }

      std::cout << available_points.size() << " available points and ";
      std::cout << constraints.size() << " constraints : ";
      for(int i=0; i<constraints.size(); ++i)
        std::cout << constraints[i] << " ";
      std::cout << std::endl;

      for(std::size_t i=0; i<available_points.size(); ++i)
      {
        Star_handle star_i = m_stars[available_points[i]];
        Point_2 pi = star_i->center_point();
        Vertex_handle vh = star->insert_to_star(pi, available_points[i], false);
        vh->info() = available_points[i];

        if(std::find(constraints.begin(), constraints.end(), vh->info()) != constraints.end())
        {
          std::cout << "found: " << vh->info() << std::endl;
          std::cout << "constraining: " << star->center()->info() << " " << vh->info() << std::endl;
          star->insert_constraint(star->center(), vh);
        }
      }

      is_star_built[next_star] = true;
      stars_left--;
      star->clean();

      //re-insert constrains after clean
      typename Star::Finite_vertices_iterator vit = star->finite_vertices_begin();
      typename Star::Finite_vertices_iterator vend = star->finite_vertices_end();
      for(; vit != vend; vit++)
      {
        if(std::find(constraints.begin(), constraints.end(), vit->info()) != constraints.end())
        {
          std::cout << "found: " << vit->info() << std::endl;
          std::cout << "constraining: " << star->center()->info() << " " << vit->info() << std::endl;
          star->insert_constraint(star->center(), vit);
        }
      }

      std::cout << "built: " << star->index_in_star_set() << std::endl;
      std::cout << "it has "<< star->number_of_vertices() << " vertices and ";
      std::cout << star->number_of_faces() << " faces" << std::endl;
      star->print_vertices();
      star->print_faces();

      std::ostringstream str;
      str << "constrain_" << m_stars.size()-stars_left << ".mesh" << std::ends;
      std::ofstream out(str.str().c_str());
      output_medit(*this, out, false, ghost_faces);

      //find the ghost faces of the new star
      Face_handle_handle fit = star->finite_incident_faces_begin();
      Face_handle_handle fend = star->finite_incident_faces_end();
      for(; fit!=fend; fit++)
      {
        Face_handle fh = *fit;
        if(!constrain_test(star, fh, ghost_faces))
          ghost_faces.insert(Facet_ijk(fh));
      }

      std::cout << ghost_faces.size() << " ghosts, spooky!" << std::endl;

      //find the next star to build
#if 0
      vit = star->finite_vertices_begin();
      for(; vit != vend; vit++)
      {
        Star_handle star_i = m_stars[vit->info()];
        FT gamma = star->metric().compute_distortion(star_i->metric());
        std::cout << "compute distortion with: " << vit->info() << " " << gamma << std::endl;
        DIST_POINT dp(gamma, vit->info());
        next_stars.insert(dp);
      }
      std::cout << next_stars.size() << " stars to pick from" << std::endl;

      bool found = false;
      while(!found && !next_stars.empty())
      {
        next_star = next_stars.begin()->id;
        next_stars.erase(next_stars.begin());
        if(!is_star_built[next_star])
          found = true;
      }
      CGAL_assertion(found == true || (next_stars.empty() && stars_left==0));
#else
      vit = star->finite_vertices_begin();
      for(; vit != vend; vit++)
      {
        queue.push_back(vit->info());
      }
      std::cout << queue.size() << " stars to pick from" << std::endl;

      bool found = false;
      while(!found && !queue.empty())
      {
        next_star = queue.front();
        queue.pop_front();
        std::cout << "next star: " << next_star << std::endl;
        std::cout << queue.size() << " stars to pick from" << std::endl;
        if(!is_star_built[next_star])
          found = true;
      }
      CGAL_assertion(found == true || (queue.empty() && stars_left==0));
#endif

    }

    //output
    std::cout << "consistency: " << is_consistent(true) << std::endl;
    std::ofstream out("constrain.mesh");
    output_medit(*this, out, false, ghost_faces);
  }

  void output_radius_of_stars()
  {
    typename Traits::Compute_squared_distance_2 sqd =
        Traits().compute_squared_distance_2_object();
    std::ofstream out("stars_radius.txt");
    out << m_stars.size() << std::endl;

    for(std::size_t si=0; si<m_stars.size(); ++si)
    {
      Star_handle star_i = m_stars[si];

      for(std::size_t sj=0; sj<m_stars.size(); ++sj)
      {
        if(si!=sj)
        {
          Star_handle star_j = m_stars[sj];
          TPoint_2 tp = star_i->metric().transform(star_j->center_point());
          FT sq_d = sqd(star_i->center()->point(), tp);
          out << sq_d << " ";
        }
        else // at sj == si, we put the max distance to the adj vertices in star_i
        {
          FT sq_r = -1.;
          Vertex_handle_handle vhit = star_i->finite_adjacent_vertices_begin();
          Vertex_handle_handle vhend = star_i->finite_adjacent_vertices_end();
          for(; vhit != vhend; vhit++)
          {
            typename Traits::Compute_squared_distance_2 sqd =
                                   Traits().compute_squared_distance_2_object();
            FT sq_d = sqd(star_i->center()->point(), (*vhit)->point());
            if(sq_d > sq_r)
              sq_r = sq_d;
          }
          out << sq_r << " ";
          std::cout << "sqr: " << sq_r << std::endl;
        }
      }
      out << std::endl;
    }
  }

  void draw_metric_vector_field()
  {
    FT e0, e1;
    Vector_2 v0, v1;

    std::ofstream min_vector_field_os("vector_field_min.polylines.cgal");
    std::ofstream max_vector_field_os("vector_field_max.polylines.cgal");

    for(std::size_t sj=0; sj<m_stars.size(); ++sj)
    {
      Star_handle star_j = m_stars[sj];
      Eigen::Matrix2d m = star_j->metric().get_mat();

      get_eigen_vecs_and_vals<K>(m, v0, v1, e0, e1);

      e0 = 1./std::sqrt(std::abs(e0));
      e1 = 1./std::sqrt(std::abs(e1));

      Point_2 pj = star_j->center_point();
      FT scale = 0.1;

      if(e0>e1)
      {
        max_vector_field_os << "2 " << pj << " 0 " << (pj+scale*e0*v0) << " 0 " << std::endl;
        min_vector_field_os << "2 " << pj << " 0 " << (pj+scale*e1*v1) << " 0 " << std::endl;
      }
      else
      {
        max_vector_field_os << "2 " << pj << " 0 " << (pj+scale*e1*v1) << " 0 " << std::endl;
        min_vector_field_os << "2 " << pj << " 0 " << (pj+scale*e0*v0) << " 0 " << std::endl;
      }
    }
  }

  template<typename MF>
  void draw_metric_vector_field_2(MF const * const  mf) const
  {
    //Fixme, it's pretty bad to use a fixed value of offset here...
    const FT offset_x = -0.55; // offset is the bottom left point
    const FT offset_y = -0.55; // todo normalize this with aniso_mesh_2's rectangle
    const unsigned int n = 70;
    const FT grid_side = 1.1; // whose offset is the center of the rectangle...
    const FT step = grid_side / n;

    FT e0, e1;
    Vector_2 v0, v1;

    std::ofstream min_vector_field_os("vector_field_min.polylines.cgal");
    std::ofstream max_vector_field_os("vector_field_max.polylines.cgal");

    for(unsigned int i=0; i<n; ++i)
    {
      for(unsigned int j=0; j<n ;++j)
      {
        Point_2 p(offset_x+j*step, offset_y+i*step);

        typename MF::Metric m_p = mf->compute_metric(p);
        Eigen::Matrix2d m = m_p.get_mat();

        get_eigen_vecs_and_vals<K>(m, v0, v1, e0, e1);

        e0 = 1./std::sqrt(std::abs(e0));
        e1 = 1./std::sqrt(std::abs(e1));

        FT scale = 0.1;

        if(e0>e1)
        {
          max_vector_field_os << "2 " << p << " 0 " << (p+scale*e0*v0) << " 0 " << std::endl;
          min_vector_field_os << "2 " << p << " 0 " << (p+scale*e1*v1) << " 0 " << std::endl;
        }
        else
        {
          max_vector_field_os << "2 " << p << " 0 " << (p+scale*e1*v1) << " 0 " << std::endl;
          min_vector_field_os << "2 " << p << " 0 " << (p+scale*e0*v0) << " 0 " << std::endl;
        }
      }
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

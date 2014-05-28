#ifndef CGAL_ANISOTROPIC_MESH_3_STARSET_DISPLAY_H
#define CGAL_ANISOTROPIC_MESH_3_STARSET_DISPLAY_H

#include <CGAL/Starset.h>
#include <CGAL/Stretched_Delaunay_3.h>

#include <CGAL/helpers/combinatorics_helper.h>
#include <CGAL/gl_draw/drawing_helper.h>

#include <set>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename Starset_with_info>
class Starset_display
{
private:
  typedef Starset_display<Starset_with_info>          Self;
public:
  typedef typename Starset_with_info::Star            Star;
  typedef Star*                                       Star_handle;
  typedef std::vector<Star_handle>                    Star_vector;
  typedef typename Star_vector::iterator              Star_iterator;
  typedef typename Star::Kernel                       K;
  typedef typename Star::FT                           FT;
  typedef typename Star::Point_3                      Point_3;
  typedef typename Star::TPoint_3                     TPoint_3;
  typedef typename Star::Facet                        Facet;
  typedef typename Star::Facet_set_iterator           Facet_set_iterator;
  typedef typename Star::Cell_handle                  Cell_handle;
  typedef typename Star::Cell_handle_handle           Cell_handle_handle;

private:
  const Starset_with_info& m_starset;

public:
  void gl_draw(const typename K::Plane_3& plane,
               const bool draw_edges,
               const int star_id = -1/*only this one*/) const
  {
    GLboolean was = (::glIsEnabled(GL_LIGHTING));
    if(!was)
      ::glEnable(GL_LIGHTING);

    bool draw_all = (star_id < 0);
    std::size_t start = draw_all ? 0 : star_id;
    std::size_t N = draw_all ? m_starset.size() : (star_id + 1);

    std::set<Facet_ijk> done;

    for(std::size_t i = start; i < N; i++)
    {
      float rc = 1.; //((float) i)/((float) m_starset.size()); (variation of colors)
      Star_handle star = m_starset[i];

      if(m_starset.criteria()->dimension == 2)
      {
        Facet_set_iterator fit = star->restricted_facets_begin();
        Facet_set_iterator fitend = star->restricted_facets_end();
        for(; fit != fitend; fit++)
        {
          Facet f = *fit;
          std::pair<typename std::set<Facet_ijk>::iterator, bool> is_insert_successful;
          is_insert_successful = done.insert(Facet_ijk(f));
          if(!is_insert_successful.second)
            continue;

          const Point_3& pa = star->metric().inverse_transform(f.first->vertex((f.second+1)%4)->point());
          const Point_3& pb = star->metric().inverse_transform(f.first->vertex((f.second+2)%4)->point());
          const Point_3& pc = star->metric().inverse_transform(f.first->vertex((f.second+3)%4)->point());
          if(is_above_plane<K>(plane, pa, pb, pc))
            gl_draw_triangle<K>(pa, pb, pc, EDGES_AND_FACES, rc*191., rc*29., 72.);
        }
      }
      else
      {
        typename Star::Base::Finite_cells_iterator cit = star->finite_cells_begin();
        typename Star::Base::Finite_cells_iterator citend = star->finite_cells_end();
        for(; cit != citend; cit++)
        {
          for(int i=0; i<4; ++i)
          {
            Facet f(cit, i);
            std::pair<typename std::set<Facet_ijk>::iterator, bool> is_insert_successful;
            is_insert_successful = done.insert(Facet_ijk(f));
            if(!is_insert_successful.second)
              continue;

            const Point_3& pa = star->metric().inverse_transform(cit->vertex((i+1)%4)->point());
            const Point_3& pb = star->metric().inverse_transform(cit->vertex((i+2)%4)->point());
            const Point_3& pc = star->metric().inverse_transform(cit->vertex((i+3)%4)->point());
            if(is_above_plane<K>(plane, pa, pb, pc))
              gl_draw_triangle<K>(pa, pb, pc, EDGES_AND_FACES, rc*191., rc*29., 72.);
          }
        }
      }
    }
    if(!was)
      ::glDisable(GL_LIGHTING);
  }

  void gl_draw_metric_colors() const
  {
    std::ifstream input("geometry_input.off");
    std::ifstream ratio_colors("metric_colors.txt");
    if(!input || !ratio_colors)
    {
      std::cout << "\nWarning : file does not exist" << std::endl;
      return;
    }

    int nv, nv2, nt;
    input >> nv >> nt;
    ratio_colors >> nv2;

    if(nv != nv2)
      std::cout << "files not adapted" << std::endl;

    std::vector<Point_3> points(nv);
    std::vector<int> colors(nv);
    FT x,y,z;
    int color;

    std::cout << "nv, nt : " << nv << " " << nt << std::endl;

    for(int i=0; i<nv; ++i)
    {
      input >> x >> y >> z;
      ratio_colors >> color;
      points[i] = Point_3(x,y,z);
      colors[i] = color;
    }

    int p1, p2, p3;
    std::cout << "pts and colors :" << std::endl;
    for(int i=0; i<nv; ++i)
      std::cout << points[i].x() << " " << points[i].y() << " " << points[i].z() << " " << colors[i] << std::endl;

    ::glDisable(GL_LIGHTING);
    for(int i=0; i<nt; ++i)
    {
      input >> p1 >> p1 >> p2 >> p3; //first p1 is "3"
      ::glShadeModel(GL_SMOOTH);
      ::glBegin(GL_TRIANGLES);
      ::glColor3d(colors[p1]/256., colors[p1]/256., colors[p1]/256.);
      ::glVertex3f(points[p1].x(), points[p1].y(), points[p1].z());    // A

      ::glColor3d(colors[p2]/256., colors[p2]/256., colors[p2]/256.);
      ::glVertex3f(points[p2].x(), points[p2].y(), points[p2].z());    // B

      ::glColor3d(colors[p3]/256., colors[p3]/256., colors[p3]/256.);
      ::glVertex3f(points[p3].x(), points[p3].y(), points[p3].z());    // C
      ::glEnd();

      ::glBegin(GL_LINES);
      ::glColor3d(0.,0.,0.);
      ::glVertex3f(points[p1].x(), points[p1].y(), points[p1].z());    // A
      ::glVertex3f(points[p2].x(), points[p2].y(), points[p2].z());    // B
      ::glVertex3f(points[p1].x(), points[p1].y(), points[p1].z());    // A
      ::glVertex3f(points[p3].x(), points[p3].y(), points[p3].z());    // C
      ::glVertex3f(points[p2].x(), points[p2].y(), points[p2].z());    // B
      ::glVertex3f(points[p3].x(), points[p3].y(), points[p3].z());    // C
      ::glEnd();
    }
  }

  void gl_draw_cell(const typename K::Plane_3& plane,
                    const int star_id = -1/*only this one*/) const
  {

    if(star_id < 0) // draw them all
      for(std::size_t i = 0; i < m_starset.size(); i++)
        m_starset[i]->gl_draw_cell(plane);
    else
      m_starset[star_id]->gl_draw_cell(plane);
  }

#if 0
  void gl_draw_picked_points(const typename K::Plane_3& plane,
                             const int point_id = -1/*only this one*/) const
  {
    //facet trying to be refined
    if(!m_pick_valid_facet.empty())
    {
      gl_draw_triangle<K>(m_pick_valid_facet[0],
          m_pick_valid_facet[1],
          m_pick_valid_facet[2],
          EDGES_AND_FACES, 84, 204, 204);
    }
    GLboolean light = (::glIsEnabled(GL_LIGHTING));
    if(light)
      ::glDisable(GL_LIGHTING);

    float point_size = 5.;
    if(!red_points.empty() || !orange_points.empty() ||
       !yellow_points.empty() || !green_points.empty())
      point_size = 0.01f;

    //all pickvalid points that were tried
    if(point_id < 0)
    {
      for(std::size_t i = 0; i<m_pick_valid_cache.size(); ++i)
      {
        Point_3 pickvalid_point = m_pick_valid_cache[i];

        //keeping all these points independant of the cut plane
        // if(!is_above_plane<K>(plane, pickvalid_point))
        //   continue;

        //if(i != m_pick_valid_cache.size()-1) //point size depends on the number of problematic facets
        //  point_size += ((pickvalid_problematic_facets[pickvalid_point]).size())/2;
        ::glPointSize(point_size);

        ::glBegin(GL_POINTS);
        if( i == m_pick_valid_cache.size()-1 )
          ::glColor3f(0.87f, 0.14f, 0.14f); //circumcenter in red
        else
          ::glColor3f(0.14f, 0.87f, 0.14f); //others in green

        ::glVertex3d(pickvalid_point.x(), pickvalid_point.y(), pickvalid_point.z());
        ::glEnd();
      }

      if(light)
        ::glEnable(GL_LIGHTING);
    }
    else //only one
    {
      assert(point_id >= 0 && point_id < m_pick_valid_cache.size());
      Point_3 picked_point = m_pick_valid_cache[point_id];
      std::vector<int> picked_point_couples = pickvalid_problematic_facets[picked_point];

      //the pickvalid point
      ::glPointSize(5.);
      ::glBegin(GL_POINTS);
      ::glColor3f(0.14f, 0.87f, 0.14f);
      ::glVertex3d(picked_point.x(), picked_point.y(), picked_point.z());
      ::glEnd();

      if(light)
        ::glEnable(GL_LIGHTING);

      //all its problematic facets
      for(std::size_t i=0; i<picked_point_couples.size(); i+=2)
      {
        gl_draw_triangle<K>(picked_point,
                            m_starset[picked_point_couples[i]]->center_point(),
            m_starset[picked_point_couples[i+1]]->center_point(),
            EDGES_AND_FACES, 166, 247, 170);
      }
    }
  }

  void gl_draw_initial_points(const typename K::Plane_3& plane) const
  {
    typename Constrain_surface::Pointset::const_iterator pi = initial_points.begin();
    typename Constrain_surface::Pointset::const_iterator pend = initial_points.end();
    GLboolean light = (::glIsEnabled(GL_LIGHTING));

    ::glPointSize(10.);
    ::glColor3f(0.706f, 0.345f, 0.878f);
    if(light)
      ::glDisable(GL_LIGHTING);
    ::glBegin(GL_POINTS);
    for (; pi != pend; pi++)
    {
      // if(!is_above_plane<K>(plane, *it))
      //   continue;

      ::glVertex3d((*pi).x(), (*pi).y(), (*pi).z());
    }

    ::glEnd();
    if(light)
      ::glEnable(GL_LIGHTING);
  }

  void gl_draw_poles(const typename K::Plane_3& plane) const
  {
    GLboolean light = (::glIsEnabled(GL_LIGHTING));
    ::glPointSize(5.);
    ::glColor3f(0.98f, 0.757f, 0.137f);
    if(light)
      ::glDisable(GL_LIGHTING);

    ::glBegin(GL_POINTS);
    typename std::set<Point_3>::const_iterator it;
    for(it = m_poles.begin(); it != m_poles.end(); ++it)
    {
      if(!is_above_plane<K>(plane, *it))
        continue;
      ::glVertex3d((*it).x(), (*it).y(), (*it).z());
    }
    ::glEnd();
    if(light)
      ::glEnable(GL_LIGHTING);
  }

  void gl_draw_metric(const typename K::Plane_3& plane,
                      double bbox_min, double eps,
                      const int star_id = -1/*only this one*/) const
  {
    double coeff = bbox_min/10.;
    double glob_min = (std::max)(eps, m_pConstrain->global_min_curvature());
    std::cout << "eps, min curv, coeff : " << eps << " " << m_pConstrain->global_min_curvature() << " " << coeff << std::endl;
    coeff *= glob_min;

    if(star_id < 0) // draw them all
      for(std::size_t i = 0; i < m_starset.size(); i++)
        m_starset[i]->gl_draw_metric(plane, coeff);
    else
      m_starset[star_id]->gl_draw_metric(plane, coeff);
  }
#endif

  void gl_draw_approx_error(const typename K::Plane_3& plane,
                            const int star_id = -1) const
  {
    std::ofstream out("approx_error.txt");

    GLboolean was = (::glIsEnabled(GL_LIGHTING));
    if(was)
      ::glDisable(GL_LIGHTING);

    bool draw_all = (star_id < 0);
    std::size_t start = draw_all ? 0 : star_id;
    std::size_t N = draw_all ? m_starset.size() : (star_id + 1);

    std::set<Facet_ijk> done;
    for(std::size_t i = start; i < N; i++)
    {
      Star_handle star = m_starset[i];
      if(!star->is_surface_star())
        continue;

      Facet_set_iterator fit = star->begin_restricted_facets();
      Facet_set_iterator fend = star->end_restricted_facets();
      for(; fit != fend; fit++)
      {
        Facet f = *fit;
        std::pair<typename std::set<Facet_ijk>::iterator, bool> is_insert_successful;
        is_insert_successful = done.insert(Facet_ijk(f));
        if(!is_insert_successful.second)
          continue;

        Point_3 cc;
        star->compute_dual_intersection(f, cc);
        if(!m_refinement_condition(cc))
          continue;

        const Point_3& pa = transform_from_star_point(f.first->vertex((f.second+1)%4)->point(), star);
        const Point_3& pb = transform_from_star_point(f.first->vertex((f.second+2)%4)->point(), star);
        const Point_3& pc = transform_from_star_point(f.first->vertex((f.second+3)%4)->point(), star);

        if(!is_above_plane<K>(plane, pa, pb, pc))
          continue;

        Star_handle star_a = m_starset[f.first->vertex((f.second+1)%4)->info()];
        Star_handle star_b = m_starset[f.first->vertex((f.second+2)%4)->info()];
        Star_handle star_c = m_starset[f.first->vertex((f.second+3)%4)->info()];

        FT sqd_a = sq_distance_to_surface(f, star_a);
        FT sqd_b = sq_distance_to_surface(f, star_b);
        FT sqd_c = sq_distance_to_surface(f, star_c);
        FT error = std::sqrt((std::max)((std::max)(sqd_a, sqd_b), sqd_c));
        out << cc.x() << " " << cc.y() << " " << cc.z() << " |||||||||| " << error << std::endl;
        error *= 1250;

        float q = (float)error;
        gl_draw_triangle<K>(pa, pb, pc, FACES_ONLY, 255.f, q, q);
      }
    }
    if(was)
      ::glEnable(GL_LIGHTING);
  }

  void gl_draw_dual(const typename K::Plane_3& plane,
                    const int star_id = -1) const
  {
    if(star_id < 0) // draw them all
      for(std::size_t i = 0; i < m_starset.size(); i++)
        m_starset[i]->gl_draw_dual(plane);
    else
      m_starset[star_id]->gl_draw_dual(plane);
  }

  void gl_draw_surface_delaunay_balls(const typename K::Plane_3& plane,
                                      const int star_id = -1) const
  {
    if(star_id < 0)
      for(std::size_t i = 0; i < m_starset.size(); i++)
        m_starset[i]->gl_draw_surface_delaunay_balls(plane);
    else
      m_starset[star_id]->gl_draw_surface_delaunay_balls(plane);
  }

  void gl_draw_distortion(const typename K::Plane_3& plane,
                          const int star_id = -1) const
  {
    GLboolean was = (::glIsEnabled(GL_LIGHTING));
    if(was)
      ::glDisable(GL_LIGHTING);

    bool draw_all = (star_id < 0);
    std::size_t start = draw_all ? 0 : star_id;
    std::size_t N = draw_all ? m_starset.size() : (star_id + 1);
    double global_max_dis = m_starset.average_facet_distortion(false /*verbose*/);

    std::set<Facet_ijk> done;
    for(std::size_t i = start; i < N; i++)
    {
      Star_handle star = m_starset[i];
      typename Star::Facet_set_iterator fit = star->restricted_facets_begin();
      typename Star::Facet_set_iterator fitend = star->restricted_facets_end();
      for(; fit != fitend; fit++)
      {
        Facet f = *fit;
        std::pair<typename std::set<Facet_ijk>::iterator, bool> is_insert_successful;
        is_insert_successful = done.insert(Facet_ijk(f));
        if(!is_insert_successful.second)
          continue;

        const Point_3& pa = star->metric().inverse_transform(f.first->vertex((f.second+1)%4)->point());
        const Point_3& pb = star->metric().inverse_transform(f.first->vertex((f.second+2)%4)->point());
        const Point_3& pc = star->metric().inverse_transform(f.first->vertex((f.second+3)%4)->point());
        if(!is_above_plane<K>(plane, pa, pb, pc))
          continue;

        FT max_distortion = 0.;
        for (int i = 0; i < 3; i++)
        {
          int index_1 = (f.second + i + 1) % 4;
          int index_2 = (f.second + (i + 1) % 3 + 1) % 4;
          FT distortion = m_starset[f.first->vertex(index_1)->info()]->metric().compute_distortion(
                            m_starset[f.first->vertex(index_2)->info()]->metric());
          max_distortion = (std::max)(distortion, max_distortion);
        }
        float rgf = static_cast<float>(max_distortion/global_max_dis*255.);
        gl_draw_triangle<K>(pa, pb, pc, FACES_ONLY, rgf, rgf, 255);
      }
    }
    if(was)
      ::glEnable(GL_LIGHTING);
  }

  void gl_draw_metric_honoring(const typename K::Plane_3& plane,
                               const int star_id = -1) const
  {
    GLboolean was = (::glIsEnabled(GL_LIGHTING));
    if(was)
      ::glDisable(GL_LIGHTING);

    bool draw_all = (star_id < 0);
    std::size_t start = draw_all ? 0 : star_id;
    std::size_t N = draw_all ? m_starset.size() : (star_id + 1);

    std::set<Facet_ijk> done;
    for(std::size_t i = start; i < N; i++)
    {
      Star_handle star = m_starset[i];
      if(!star->is_surface_star())
        continue;

      Facet_set_iterator fit = star->restricted_facets_begin();
      Facet_set_iterator fend = star->restricted_facets_end();
      for(; fit != fend; fit++)
      {
        Facet f = *fit;
        std::pair<typename std::set<Facet_ijk>::iterator, bool> is_insert_successful;
        is_insert_successful = done.insert(Facet_ijk(f));
        if(!is_insert_successful.second)
          continue;

        const Point_3& pa = star->metric().inverse_transform(f.first->vertex((f.second+1)%4)->point());
        const Point_3& pb = star->metric().inverse_transform(f.first->vertex((f.second+2)%4)->point());
        const Point_3& pc = star->metric().inverse_transform(f.first->vertex((f.second+3)%4)->point());

        if(!is_above_plane<K>(plane, pa, pb, pc))
          continue;

        Star_handle star_a = m_starset[f.first->vertex((f.second+1)%4)->info()];
        Star_handle star_b = m_starset[f.first->vertex((f.second+2)%4)->info()];
        Star_handle star_c = m_starset[f.first->vertex((f.second+3)%4)->info()];

        const TPoint_3& tpa_a = star_a->metric().transform(pa);
        const TPoint_3& tpb_a = star_a->metric().transform(pb);
        const TPoint_3& tpc_a = star_a->metric().transform(pc);

        FT qualityf_in_a = 255.*star_a->compute_element_quality(tpa_a, tpb_a, tpc_a);

        const TPoint_3& tpa_b = star_b->metric().transform(pa);
        const TPoint_3& tpb_b = star_b->metric().transform(pb);
        const TPoint_3& tpc_b = star_b->metric().transform(pc);

        FT qualityf_in_b = 255.*star_b->compute_element_quality(tpa_b, tpb_b, tpc_b);

        const TPoint_3& tpa_c = star_c->metric().transform(pa);
        const TPoint_3& tpb_c = star_c->metric().transform(pb);
        const TPoint_3& tpc_c = star_c->metric().transform(pc);

        FT qualityf_in_c = 255.*star_c->compute_element_quality(tpa_c, tpb_c, tpc_c);
        FT qualityf = (std::min)((std::min)(qualityf_in_a, qualityf_in_b),qualityf_in_c);

        float q = (float)qualityf;
        gl_draw_triangle<K>(pa, pb, pc, FACES_ONLY, q, 255.f, q);
      }
    }
    if(was)
      ::glEnable(GL_LIGHTING);
  }

  void gl_draw_inconsistent_facets(const typename K::Plane_3& plane,
                                   const int star_id = -1) const
  {
    GLboolean was = (::glIsEnabled(GL_LIGHTING));
    if(!was)
      ::glEnable(GL_LIGHTING);

    ::glPolygonOffset(1.f, 0.1f);
    ::glEnable(GL_POLYGON_OFFSET_FILL);
    ::glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    bool draw_all = (star_id < 0);
    std::size_t start = draw_all ? 0 : star_id;
    std::size_t N = draw_all ? m_starset.size() : (star_id + 1);

    for(std::size_t i = start; i < N; i++)
    {
      Star_handle star = m_starset[i];
      typename Star::Facet_set_iterator fit = star->restricted_facets_begin();
      typename Star::Facet_set_iterator fitend = star->restricted_facets_end();
      for(; fit != fitend; fit++)
      {
        Facet f = *fit;
        std::vector<bool> inconsistent_points(3, false);
        if(m_starset.is_consistent(f, inconsistent_points))
          continue;

        const Point_3& pa = star->metric().inverse_transform(f.first->vertex((f.second+1)%4)->point());
        const Point_3& pb = star->metric().inverse_transform(f.first->vertex((f.second+2)%4)->point());
        const Point_3& pc = star->metric().inverse_transform(f.first->vertex((f.second+3)%4)->point());

        if(is_above_plane<K>(plane, pa, pb, pc))
        {
          gl_draw_triangle<K>(pa,pb,pc,EDGES_AND_FACES, 227,27,27);

          FT m_x = (1./3.)*(pa.x()+pb.x()+pc.x());
          FT m_y = (1./3.)*(pa.y()+pb.y()+pc.y());
          FT m_z = (1./3.)*(pa.z()+pb.z()+pc.z());
          Point_3 m(m_x,m_y,m_z);

          ::glColor3f(0.10f, 0.90f, 0.10f);
          if(!inconsistent_points[0])
            gl_draw_segment<K>(m, pa);
          if(!inconsistent_points[1])
            gl_draw_segment<K>(m, pb);
          if(!inconsistent_points[2])
            gl_draw_segment<K>(m, pc);
        }
      }
    }
    ::glDisable(GL_POLYGON_OFFSET_FILL);
    if(!was)
      ::glDisable(GL_LIGHTING);
  }

  Starset_display(const Starset_with_info& starset_) : m_starset(starset_) { }

private:
  Starset_display(const Self&) { }
  Self& operator=(const Self&) { }
};


}  // Anisotropic_mesh_3
}  // CGAL

#endif

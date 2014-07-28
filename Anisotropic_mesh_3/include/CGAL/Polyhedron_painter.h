#ifndef CGAL_ANISOTROPIC_MESH_3_POLYHEDRON_PAINTER_H
#define CGAL_ANISOTROPIC_MESH_3_POLYHEDRON_PAINTER_H

#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>
#include <CGAL/Metric_field.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/gl_draw/drawing_helper.h>
#include <CGAL/helpers/aniso_polyhedron_items.h>
#include <CGAL/helpers/colored_polyhedron_output.h>
#include <CGAL/helpers/metric_helper.h>
#include <CGAL/helpers/mpq_helper.h>
#include <CGAL/helpers/mvpq_helper.h>

#include <CGAL/gl.h>
#include <Eigen/Dense>

#include <cstdio>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K,
         typename Constrain_surface,
         typename Metric_field,
         typename Colored_polyhedron = CGAL::Polyhedron_3<K, Aniso_items> >
class Polyhedron_painter
{
public:
  typedef typename K::FT                       FT;
  typedef typename K::Vector_3                 Vector_3;
  typedef typename K::Point_3                  Point_3;
  typedef typename K::Object_3                 Object_3;
  typedef typename K::Segment_3                Segment_3;
  typedef typename K::Triangle_3               Triangle_3;
  typedef typename K::Ray_3                    Ray_3;
  typedef typename K::Line_3                   Line_3;
  typedef typename K::Plane_3                  Plane_3;

  typedef typename Metric_field::Metric       Metric;

  typedef typename Colored_polyhedron::Face_handle                        Face_handle;
  typedef typename Colored_polyhedron::Vertex_handle                      Vertex_handle;
  typedef typename Colored_polyhedron::Facet_iterator                     Facet_iterator;
  typedef typename Colored_polyhedron::Vertex_iterator                    Vertex_iterator;
  typedef typename Colored_polyhedron::Halfedge_around_vertex_circulator  HV_circulator;

  typedef CGAL::AABB_polyhedron_triangle_primitive<K, Colored_polyhedron> Primitive;
  typedef CGAL::AABB_traits<K, Primitive>                                 Traits;
  typedef CGAL::AABB_tree<Traits>                                         Tree;
  typedef typename Tree::Object_and_primitive_id                          Object_and_primitive_id;
  typedef typename Tree::Point_and_primitive_id                           Point_and_primitive_id;
  typedef typename Tree::Primitive_id                                     Primitive_id;

private:
  mutable Colored_polyhedron m_colored_poly;
  mutable std::vector<Vertex_handle> m_colored_poly_vertex_index;
  mutable Colored_polyhedron m_colored_poly_mem;
  mutable Colored_polyhedron m_colored_poly_prev;
  mutable Tree* m_colored_poly_tree;
  const Constrain_surface* const m_surface;
  const Metric_field* m_mf;

public:
  Colored_polyhedron& colored_poly() const { return m_colored_poly; }
  Colored_polyhedron& colored_poly_mem() const { return m_colored_poly_mem; }
  const Constrain_surface* const constrain_surface() const { return m_surface; }
  const Metric_field* metric_field() const { return m_mf; }

  void build_colored_polyhedron() const
  {
    constrain_surface()->build_colored_polyhedron(m_colored_poly);
  }

  void number_colored_poly() const
  {
    std::size_t index = 0;
    std::ofstream out("coordinates_to_num.txt");
    m_colored_poly_vertex_index.reserve(m_colored_poly.size_of_vertices());
    for(Vertex_iterator v = m_colored_poly.vertices_begin(); v != m_colored_poly.vertices_end(); ++v, ++index)
    {
      v->tag() = index;
      m_colored_poly_vertex_index[index] = v;
      out << v->point() << " -- " << index << std::endl;
    }
  }

  void count_colored_elements() const
  {
    int count = 0;
    for(Vertex_iterator v = m_colored_poly.vertices_begin(); v != m_colored_poly.vertices_end(); ++v)
      if(v->is_colored())
        count++;
    std::cout << count << " / " << m_colored_poly.size_of_vertices() << " colored vertices" << std::endl;
  }

  void count_green_elements() const
  {
    int count = 0;
    for(Vertex_iterator v = m_colored_poly.vertices_begin(); v != m_colored_poly.vertices_end(); ++v)
      if(v->green() > 0)
        count++;
    std::cout << count << " / " << m_colored_poly.size_of_vertices() << " green vertices" << std::endl;
  }

  void build_colored_poly_tree() const
  {
    std::cout << "build and accelerate tree" << std::endl;
    delete m_colored_poly_tree;
    m_colored_poly_tree = new Tree(m_colored_poly.facets_begin(), m_colored_poly.facets_end());
    m_colored_poly_tree->accelerate_distance_queries();
    std::cout << "tree has size : " << m_colored_poly_tree->size() << std::endl;
  }

  void reset(bool color_reset = false) const
  {
    for(Vertex_iterator v = m_colored_poly.vertices_begin(); v != m_colored_poly.vertices_end(); ++v)
    {
      v->colored_this_pass() = false;
      if(color_reset)
       v->metric() = Eigen::Matrix3d::Zero();
    }
  }

  void color_poly(const Point_3& p, const Eigen::Matrix3d& m, const bool use_previous_colors, int origin = 0) const
  {
    Point_and_primitive_id pp = m_colored_poly_tree->closest_point_and_primitive(p);
    Point_3 closest_point = pp.first;
    Face_handle closest_f = pp.second;

    Vertex_handle va = closest_f->halfedge()->vertex();
    Vertex_handle vb = closest_f->halfedge()->next()->vertex();
    Vertex_handle vc = closest_f->halfedge()->next()->next()->vertex();

    Vertex_handle v = va;
    FT dmin = CGAL::squared_distance(closest_point, va->point());
    FT db = CGAL::squared_distance(closest_point, vb->point());
    FT dc = CGAL::squared_distance(closest_point, vc->point());
    if(db < dmin)
    {
      dmin = db;
      v = vb;
    }
    if(dc < dmin)
      v = vc;

    if(!v->colored_this_pass())
    {
      v->green()++;
      v->contributors()++;
      v->metric_origin() = origin;
      v->colored_this_pass() = true;

      if(use_previous_colors)
      {
        std::cout << "brute force search in the previous FULLY colored polyhedron" << std::endl;
      //super brute force / ugly, replace todo
        Vertex_iterator vi;
        for(vi = m_colored_poly_prev.vertices_begin(); vi != m_colored_poly_prev.vertices_end(); ++vi)
          if(vi->point() == v->point())
            break;
        v->metric() = vi->metric();
      }
      else
      {
        double n = v->contributors();
        if(n >= 2)
        {
          std::vector<std::pair<Eigen::Matrix3d, FT> > w_metrics;
          w_metrics.push_back(std::make_pair(v->metric(), (n-1.)/n));
          w_metrics.push_back(std::make_pair(m, 1./n));
          v->metric() = CGAL::Anisotropic_mesh_3::logexp_interpolate<K>(w_metrics);
        }
        else
          v->metric() = m;
      }
    }

    /*
    Eigen::Matrix3d scaled_m = scale_matrix_to_point<K>(m, closest_point, v->point());
    if(!v->contributors())
    {
      v->metric() = scaled_m;
      v->metric_origin() = origin;
    }
    else
    {
      v->metric() = matrix_intersection<K>(v->metric(), scaled_m);
    }
    v->contributors()++;
    */
  }

  void color_vertices_with_acceptable_rg_ratio(double ratio_limit = 0.8,
                                               bool reset = false) const
  {
    std::cout << "coloring vertices with acceptable ratios"<< std::endl;
    for(Vertex_iterator v = m_colored_poly.vertices_begin(); v != m_colored_poly.vertices_end(); ++v)
    {

      double v_ratio = v->ratio();
      //std::cout << v->tag() << " " << v->green() << " " << v->red() << " " << v_ratio << std::endl;
      if(v_ratio > ratio_limit)
      {
        v->colored_this_pass() = true;
        if(v->metric() == Eigen::Matrix3d::Zero())
          std::cout << "acceptable but not colored...?" << std::endl;
      }
      else
      {
        v->colored_this_pass() = false;
        v->metric() = Eigen::Matrix3d::Zero();
      }

      if(reset)
      {
        Point_3 p = v->point();
        Metric mp = m_mf->compute_metric(p);
        Eigen::Matrix3d transf = mp.get_transformation();
        v->metric() = transf.transpose()*transf;
        continue;
      }
    }
    count_colored_elements();
  }

  void color_uncolored_vertices_red() const
  {
    for(Vertex_iterator v = m_colored_poly.vertices_begin(); v != m_colored_poly.vertices_end(); ++v)
      if(!v->colored_this_pass())
        v->red()++;
  }

  void color_poly(const Point_3& p, const FT& ratio, const bool useless, int uuseless = 0) const
  {
    Point_and_primitive_id pp = m_colored_poly_tree->closest_point_and_primitive(p);
    Face_handle closest_f = pp.second; // closest primitive id
    closest_f->color() += ratio;
    closest_f->contributors()++;
  }

  void average_facet_color_contributor() const
  {
    std::cout << "computing average on each facet + stats" << std::endl;

    double strongest_color = -1.;
    int colored_facets_count = 0;

    for(Facet_iterator f = m_colored_poly.facets_begin(); f != m_colored_poly.facets_end(); ++f)
    {
      if(f->contributors() > 0)
      {
        colored_facets_count++;
        f->color() /= f->contributors();
        if(f->color() > strongest_color)
          strongest_color = f->color();
      }
    }
    std::cout << colored_facets_count << "/" << m_colored_poly.size_of_facets() << " facets colored" << std::endl;

    std::ofstream out("colored_poly_ini.off");
    File_writer_OFF writer(false);
    output_colored_polyhedron(out, m_colored_poly, writer, strongest_color);
  }

  void spread_vertices_colors() const
  {
    m_colored_poly_mem = Colored_polyhedron(m_colored_poly);
    typedef Colored_modifiable_vertex_priority_queue<Colored_polyhedron>  Cvmpq;
    Cvmpq q(m_colored_poly.size_of_vertices(), typename Cvmpq::Compare(), typename Cvmpq::ID());
    q.initialize_cmvpq(m_colored_poly);
    q.color_all_vertices();
    // reset smooth etc ////////////
    color_vertices_with_acceptable_rg_ratio(0.8, true /*reset*/);
    smooth_colors();
    // /////////////////////////////
    m_colored_poly_prev = Colored_polyhedron(m_colored_poly);
  }

  void spread_facets_colors() const
  {
    m_colored_poly_mem = Colored_polyhedron(m_colored_poly);
    typedef Colored_modifiable_priority_queue<Colored_polyhedron>  Cmpq;
    Cmpq q(m_colored_poly.size_of_facets(), typename Cmpq::Compare(), typename Cmpq::ID());
    q.initialize_cmpq(m_colored_poly);
    q.color_all_facets();
    m_colored_poly_prev = Colored_polyhedron(m_colored_poly);
  }

  void get_enough_rings(std::set<Vertex_handle>& vs,
                        Vertex_handle v) const
  {
    int nb_of_rings = 0;
    bool found_colored = false;
    std::vector<Vertex_handle> vs_prev_level;
    std::vector<Vertex_handle> vs_current_level;
    vs_prev_level.push_back(v);

    int one_more = 1;

    while(!found_colored || one_more--)
    {
      nb_of_rings++;
      while(!vs_prev_level.empty())
      {
        Vertex_handle vi = vs_prev_level.back();
        vs_prev_level.pop_back();

        HV_circulator h = vi->vertex_begin();
        HV_circulator hend = h;
        do
        {
          Vertex_handle vj = h->opposite()->vertex();
          vs_current_level.push_back(vj);
          if(vj->colored_this_pass())
            found_colored = true;
          vs.insert(vj);
          h++;
        }
        while(h != hend);
      }
      vs_prev_level = vs_current_level;
      vs_current_level.clear();
    }
    std::cout << nb_of_rings << " rings and " << vs.size() << " vertices" << std::endl;
  }

  void smooth_colors(int nb = 10) const
  {
    while(nb--)
    {
      std::cout << "smooth " << nb << std::endl;
      std::vector<Eigen::Matrix3d> new_metrics(m_colored_poly.size_of_vertices());

      typedef typename Colored_polyhedron::Vertex_iterator VI;
      for (VI vi = m_colored_poly.vertices_begin();
              vi != m_colored_poly.vertices_end(); ++vi)
      {
        std::set<Vertex_handle> vs;

        //if we want the number of rings to be enough to reach a consistent zone
        //get_enough_rings(vs, vi);

        HV_circulator h = vi->vertex_begin();
        HV_circulator hend = h;
        do
        {
          Vertex_handle vj = h->opposite()->vertex();

          /* //second ring
          HV_circulator h2 = vj->vertex_begin();
          HV_circulator h2end = h;
          do
          {
            Vertex_handle vk = h2->opposite()->vertex();
            vs.insert(vk);
            h2++;
          }
          while(h2 != h2end);
          //end sr */

          vs.insert(vj);
          h++;
        }
        while(h != hend);

        std::vector<std::pair<Eigen::Matrix3d, typename K::FT> > w_metrics;
        double det = 1./(2.*vs.size());

        typename std::set<Vertex_handle>::iterator it = vs.begin();
        typename std::set<Vertex_handle>::iterator itend = vs.end();
        for(; it!=itend; ++it)
          w_metrics.push_back(std::make_pair((*it)->metric(), det));
        w_metrics.push_back(std::make_pair(vi->metric(), 0.5));

        new_metrics[vi->tag()] = logexp_interpolate<K>(w_metrics);
      }

      for(Vertex_iterator v = m_colored_poly.vertices_begin();
                          v != m_colored_poly.vertices_end();
                          ++v)
        v->metric() = new_metrics[v->tag()];

      count_colored_elements();
    }
  }

  void get_color_from_poly(const Point_3& p, Eigen::Matrix3d& m) const
  {
    Point_and_primitive_id pp = m_colored_poly_tree->closest_point_and_primitive(p);
    Point_3 closest_point = pp.first;
    Face_handle closest_f = pp.second; // closest primitive id

    Vertex_handle va = closest_f->halfedge()->vertex();
    Vertex_handle vb = closest_f->halfedge()->next()->vertex();
    Vertex_handle vc = closest_f->halfedge()->next()->next()->vertex();

    //get bary weights
    Vector_3 v0(va->point(), vb->point());
    Vector_3 v1(va->point(), vc->point());
    Vector_3 v2(va->point(), closest_point);
    FT d00 = v0*v0;
    FT d01 = v0*v1;
    FT d11 = v1*v1;
    FT d20 = v2*v0;
    FT d21 = v2*v1;
    FT denom = d00 * d11 - d01 * d01;
    FT v = (d11 * d20 - d01 * d21) / denom;
    FT w = (d00 * d21 - d01 * d20) / denom;
    FT u = 1.0f - v - w;

    std::vector<std::pair<Eigen::Matrix3d, FT> > w_metrics;
    w_metrics.push_back(std::make_pair(va->metric(), u));
    w_metrics.push_back(std::make_pair(vb->metric(), v));
    w_metrics.push_back(std::make_pair(vc->metric(), w));
    m = logexp_interpolate<K>(w_metrics);
  }

  void get_color_from_poly(const Point_3& p, FT& ratio) const
  {
    Point_and_primitive_id pp = m_colored_poly_tree->closest_point_and_primitive(p);
    Point_3 closest_point = pp.first;
    Face_handle closest_f = pp.second; // closest primitive id
    ratio = closest_f->color();
  }

  void gl_draw_ellipsoid_with_origin_color(Vertex_handle& vi) const
  {
    Point_3 pi = vi->point();
    Vector_3 vn, v1, v2;
    FT en, e1, e2;
    get_eigen_vecs_and_vals<K>(vi->metric(), vn, v1, v2, en, e1, e2);
    if(vi->metric_origin() == 0)
      gl_draw_ellipsoid<K>(CGAL::ORIGIN, pi, 10, 10,
                           1./std::sqrt(en), 1./std::sqrt(e1), 1./std::sqrt(e2),
                           vn, v1, v2, 43, 39, 234);
    else if(vi->metric_origin() == 1)
      gl_draw_ellipsoid<K>(CGAL::ORIGIN, pi, 10, 10,
                           1./std::sqrt(en), 1./std::sqrt(e1), 1./std::sqrt(e2),
                           vn, v1, v2, 43, 239, 234);
    else if(vi->metric_origin() == 2)
      gl_draw_ellipsoid<K>(CGAL::ORIGIN, pi, 10, 10,
                           1./std::sqrt(en), 1./std::sqrt(e1), 1./std::sqrt(e2),
                           vn, v1, v2, 43, 239, 104);
    else if(vi->metric_origin() == 3)
      gl_draw_ellipsoid<K>(CGAL::ORIGIN, pi, 10, 10,
                           1./std::sqrt(en), 1./std::sqrt(e1), 1./std::sqrt(e2),
                           vn, v1, v2, 243, 239, 13);
  }

  void gl_draw_colored_polyhedron_one_vertex(const int vertex_num) const
  {
    Vertex_handle v = m_colored_poly_vertex_index[vertex_num];
    std::cout << "-------------------------------------------------------------" << std::endl;
    std::cout << "draw one : " << v->tag() << " rank is : " << v->colored_rank() << std::endl;
    std::cout << "tensor : " << std::endl << v->metric() << std::endl;
    Point_3 p = v->point();
    int rank = v->colored_rank();

    Vector_3 vn, v1, v2;
    FT en, e1, e2;

    if(v->is_colored())
    {
      get_eigen_vecs_and_vals<K>(v->metric(), vn, v1, v2, en, e1, e2);
      gl_draw_ellipsoid<K>(CGAL::ORIGIN, p, 20, 20,
                           1./std::sqrt(en), 1./std::sqrt(e1), 1./std::sqrt(e2),
                           vn, v1, v2, 243, 29, 13);
    }

    ::glPointSize(7.);
    ::glBegin(GL_POINTS);
    ::glColor3d(0.,0.,0.);
    ::glVertex3f(p.x(), p.y(), p.z());
    ::glEnd();

    int count = 0;

    HV_circulator h = v->vertex_begin();
    HV_circulator hend = h;
    do
    {
      Vertex_handle vi = h->opposite()->vertex();
      if(!vi->is_colored())
        continue;

      count++;
      std::cout << count << " || " << vi->tag() << " is colored with rank : " << vi->colored_rank() << std::endl;
      std::cout << vi->metric() << std::endl;
      Point_3 pi = vi->point();
      get_eigen_vecs_and_vals<K>(vi->metric(), vn, v1, v2, en, e1, e2);

      if(vi->colored_rank() > rank)
        gl_draw_ellipsoid<K>(CGAL::ORIGIN, pi, 40, 40,
                             1./std::sqrt(en), 1./std::sqrt(e1), 1./std::sqrt(e2),
                             vn, v1, v2, 243, 239, 13);
      else if(vi->colored_rank() == -1)
        gl_draw_ellipsoid<K>(CGAL::ORIGIN, pi, 40, 40,
                             1./std::sqrt(en), 1./std::sqrt(e1), 1./std::sqrt(e2),
                             vn, v1, v2, 43, 39, 234);
      else // vi rank < rank
        gl_draw_ellipsoid<K>(CGAL::ORIGIN, pi, 40, 40,
                             1./std::sqrt(en), 1./std::sqrt(e1), 1./std::sqrt(e2),
                             vn, v1, v2, 43, 230 + count, 224);
      h++;
    }
    while(h != hend);
  }

  void gl_draw_colored_polyhedron(Colored_polyhedron P,
                                  const Plane_3& plane,
                                  const int vertex_id = -1) const
  {
    if(vertex_id != -1)
    {
      gl_draw_colored_polyhedron_one_vertex(vertex_id);
      return;
    }

    GLboolean was = (::glIsEnabled(GL_LIGHTING));
    if(was)
      ::glDisable(GL_LIGHTING);

    /*
    typedef typename Colored_polyhedron::Facet_iterator  FI;
    double strongest_color = -1;
    for(FI fi = P.facets_begin(); fi != P.facets_end(); ++fi)
      if(fi->color() > strongest_color)
        strongest_color = fi->color();


    if(strongest_color > 0)
      for(FI fi = P.facets_begin(); fi != P.facets_end(); ++fi)
      {
        double f_color = fi->color()/strongest_color;
        if(f_color < 0.)
          continue;
        gl_draw_colored_polyhedron_facet<Colored_polyhedron>(fi, f_color);
      }
    */

    typedef typename Colored_polyhedron::Vertex_iterator VI;
    int index = 0;
    for (VI vi = P.vertices_begin(); vi != P.vertices_end(); ++vi, ++index)
    {
      Point_3 pi = vi->point();

      if(!is_above_plane<K>(plane, pi))
        continue;

      float tot = vi->red() + vi->green();
      float redc = ((float) vi->red())/tot;
      float greenc = ((float) vi->green())/tot;

      if(greenc > 0.8)
      {
        ::glPointSize(3.);
        ::glBegin(GL_POINTS);
        ::glColor3d(0.1,0.1,0.7);
        ::glVertex3f(pi.x(), pi.y(), pi.z());
        ::glEnd();
      }

      ::glPointSize(7.);
      ::glBegin(GL_POINTS);
      /*
      if(!vi->is_colored())
        ::glColor3d(1.,0.,0.);
      else
        ::glColor3d(0.,1.,0.);
      */
      ::glColor3d(redc,greenc,0.1);
      ::glVertex3f(pi.x(), pi.y(), pi.z());
      ::glEnd();

      if(!vi->is_colored() || (vi->metric_origin() != 0 && index%100 != 0))
        continue;

      //gl_draw_ellipsoid_with_origin_color(vi);
    }

    if(was)
      ::glEnable(GL_LIGHTING);
  }

public:
  Polyhedron_painter(const Constrain_surface* const surface_,
                     const Metric_field* mf_)
    :
    m_colored_poly(Colored_polyhedron ()),
    m_colored_poly_vertex_index(),
    m_colored_poly_mem(Colored_polyhedron ()),
    m_colored_poly_prev(Colored_polyhedron ()),
    m_colored_poly_tree(NULL),
    m_surface(surface_),
    m_mf(mf_)
  { }

  ~Polyhedron_painter()
  {
    std::cout << "You're deleting trees!!" << std::endl;
    delete m_colored_poly_tree;
  }

};

}
}

#endif // CGAL_ANISOTROPIC_MESH_3_POLYHEDRON_PAINTER_H

#ifndef CGAL_ANISOTROPIC_MESH_3_POLYHEDRON_PAINTER
#define CGAL_ANISOTROPIC_MESH_3_POLYHEDRON_PAINTER

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
    std::ofstream out("coordinates_to_num.txt");

    std::size_t index = 0;
    m_colored_poly_vertex_index.reserve(m_colored_poly.size_of_vertices());
    for(Vertex_iterator v = m_colored_poly.vertices_begin(); v != m_colored_poly.vertices_end(); ++v, ++index)
    {
      v->tag() = index;
      m_colored_poly_vertex_index[index] = v;
      out << v->point() << " -- " << index << std::endl;
    }

    index = 0;
    for(Facet_iterator f = m_colored_poly.facets_begin(); f != m_colored_poly.facets_end(); ++f, ++index)
      f->tag() = index;
  }

  void build_colored_poly_tree() const
  {
    std::cout << "build and accelerate tree" << std::endl;
    delete m_colored_poly_tree;
    m_colored_poly_tree = new Tree(m_colored_poly.facets_begin(), m_colored_poly.facets_end());
    m_colored_poly_tree->accelerate_distance_queries();
    std::cout << "tree has size : " << m_colored_poly_tree->size() << std::endl;
  }

  void clear_colors() const
  {
#ifdef ANISO_COLOR_POLY_FACETS
    for(Facet_iterator f = m_colored_poly.facets_begin(); f != m_colored_poly.facets_end(); ++f)
    {
      f->color() = 0.;
      f->contributors() = 0.;
    }
#else
    for(Vertex_iterator v = m_colored_poly.vertices_begin(); v != m_colored_poly.vertices_end(); ++v)
    {
      v->metric() = Eigen::Matrix3d::Zero();
      v->contributors() = 0.;
    }
#endif
  }

  void color_poly(const Point_3& p, const Eigen::Matrix3d& m, int origin = 0) const
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

    if(!v->contributors())
    {
      Metric mp = m_mf->compute_metric(p);
      Eigen::Matrix3d transf = mp.get_transformation();
      v->metric() = transf.transpose()*transf;
      v->metric_origin() = origin;
      v->contributors()++;
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

  void color_poly(const Point_3& p, const FT& ratio) const
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
  }

  void spread_facets_colors() const
  {
    m_colored_poly_mem = Colored_polyhedron(m_colored_poly);
    typedef Colored_modifiable_priority_queue<Colored_polyhedron>  Cmpq;
    Cmpq q(m_colored_poly.size_of_facets(), typename Cmpq::Compare(), typename Cmpq::ID());
    q.initialize_cmpq(m_colored_poly);
    q.color_all_facets();
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
    m = interpolate_colors<K>(w_metrics);
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

    typedef typename Colored_polyhedron::Facet_iterator  FI;
    typedef typename Colored_polyhedron::Vertex_iterator VI;

    double strongest_color = -1;
    for(FI fi = P.facets_begin(); fi != P.facets_end(); ++fi)
      if(fi->color() > strongest_color)
        strongest_color = fi->color();

    GLboolean was = (::glIsEnabled(GL_LIGHTING));
    if(was)
      ::glDisable(GL_LIGHTING);

    if(strongest_color > 0)
      for(FI fi = P.facets_begin(); fi != P.facets_end(); ++fi)
      {
        double f_color = fi->color()/strongest_color;
        if(f_color < 0.)
          continue;
        gl_draw_colored_polyhedron_facet<Colored_polyhedron>(fi, f_color);
      }

    int index = 0;
    for (VI vi = P.vertices_begin(); vi != P.vertices_end(); ++vi, ++index)
    {
      Point_3 pi = vi->point();

      if(!is_above_plane<K>(plane, pi))
        continue;

      ::glPointSize(7.);
      ::glBegin(GL_POINTS);
      if(!vi->is_colored())
        ::glColor3d(1.,0.,0.);
      else
        ::glColor3d(0.,1.,0.);
      ::glVertex3f(pi.x(), pi.y(), pi.z());
      ::glEnd();

      if(!vi->is_colored() || (vi->metric_origin() != 0 && index%100 != 0))
        continue;

      gl_draw_ellipsoid_with_origin_color(vi);
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

#endif

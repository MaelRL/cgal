#ifndef DRAWING_HELPER_H
#define DRAWING_HELPER_H

#include <CGAL/gl.h>
#include <CGAL/glu.h>
#include <cmath>
//#include <QColor>

enum Polygon_drawing_options { EDGES_ONLY, EDGES_AND_FACES, FACES_ONLY };



template<typename Kernel>
void gl_draw_triangle(const typename Kernel::Point_3& pa,
                      const typename Kernel::Point_3& pb,
                      const typename Kernel::Point_3& pc,
                      const Polygon_drawing_options& option,
                      const float r = 205.,
                      const float g = 175., 
                      const float b = 149.)
{
  //draw triangles
  if(option != EDGES_ONLY)
  {
    ::glColor3d(r / 256., g / 256., b / 256.);
    ::glBegin(GL_TRIANGLES);
    typename Kernel::Vector_3 n = CGAL::cross_product(pb - pa, pc -pa);
    n = n / CGAL::sqrt(n*n);
    ::glNormal3d(n.x(),n.y(),n.z());
    ::glVertex3d(pa.x(),pa.y(),pa.z());
    ::glVertex3d(pb.x(),pb.y(),pb.z());
    ::glVertex3d(pc.x(),pc.y(),pc.z());
    ::glEnd();
  }
  //draw edges
  if(option != FACES_ONLY)
  {
    GLboolean was = (::glIsEnabled(GL_LIGHTING));
    if(was)
      ::glDisable(GL_LIGHTING);

    ::glColor3f(0.,0.,0.);
    if(option == EDGES_ONLY)
      ::glColor3d(r / 256., g / 256., b / 256.);
    if(option == EDGES_AND_FACES) ::glLineWidth(2.f);
    else                          ::glLineWidth(1.f);
    ::glBegin(GL_LINE_LOOP);
    ::glVertex3d(pa.x(),pa.y(),pa.z());
    ::glVertex3d(pb.x(),pb.y(),pb.z());
    ::glVertex3d(pc.x(),pc.y(),pc.z());
    ::glEnd();

    if(was)
      ::glEnable(GL_LIGHTING);
  }
}

template<typename Kernel, typename Facet>
void gl_draw_facet(const Facet& f,
                   const Polygon_drawing_options& option,
                   const float r = 205.,
                   const float g = 175., 
                   const float b = 149.)
{
  gl_draw_triangle<Kernel>(f.first->vertex((f.second+1)%4)->point(),
                   f.first->vertex((f.second+2)%4)->point(),
                   f.first->vertex((f.second+3)%4)->point(),
                   option, r, g, b);
}


template<typename Kernel>
void gl_draw_segment(const typename Kernel::Point_3& pa,
                     const typename Kernel::Point_3& pb,
                     const bool draw_extremities = false)
{
  ::glBegin(GL_LINES);
  ::glVertex3d(pa.x(),pa.y(),pa.z());
  ::glVertex3d(pb.x(),pb.y(),pb.z());
  ::glEnd();

  if(draw_extremities)
  {
    ::glColor3d(93./256.,204./256.,232./256.);
    ::glPointSize(5.);
    ::glBegin(GL_POINTS);
    ::glVertex3d(pa.x(),pa.y(),pa.z());
    ::glVertex3d(pb.x(),pb.y(),pb.z());
    ::glEnd();
  }
}


template<typename Kernel>
void gl_draw_ellipsoid(const typename Kernel::Point_3& p,
                       int stacks, int slices, double a, double b, double c,
                       const float cr = 205., const float cg = 175., const float cb = 149.)
{
    double x = p.x();
    double y = p.y();
    double z = p.z();
    double st = CGAL_PI / (double) stacks;
    double sl = CGAL_PI / (double) slices;

    GLboolean was = (::glIsEnabled(GL_LIGHTING));
    if(!was)
      ::glEnable(GL_LIGHTING);

    ::glColor3d(cr / 256., cg / 256., cb / 256.);
    ::glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    for(double u = -CGAL_PI/2.; u <= (CGAL_PI/2.); u += sl) // theta
    {
        ::glBegin(GL_TRIANGLE_STRIP);
        for(double v = -CGAL_PI; v <= CGAL_PI; v += st) // phi
        {
          ::glVertex3d(x + a * std::cos(u) * std::cos(v),
                       y + b * std::cos(u) * std::sin(v),
                       z + c * std::sin(u));

          ::glVertex3d(x + a * std::cos(u + sl) * std::cos(v),
                       y + b * std::cos(u + sl) * std::sin(v),
                       z + c * std::sin(u + sl));
        }
        ::glEnd();
    }

    if(!was)
      ::glDisable(GL_LIGHTING);
}


template<typename Kernel>
void gl_draw_sphere(const typename Kernel::Sphere_3& s)
{
  typename Kernel::Point_3 c = s.center();
  const GLdouble r = std::sqrt(s.squared_radius());

  GLboolean was = (::glIsEnabled(GL_LIGHTING));
  if(!was)
    ::glEnable(GL_LIGHTING);

  GLint polygon_mode[2];
  ::glPushMatrix();

  ::glGetIntegerv(GL_POLYGON_MODE, &polygon_mode[0]);
  ::glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
  ::glTranslated(CGAL::to_double(c.x()), CGAL::to_double(c.y()), CGAL::to_double(c.z()));

  GLUquadric *quad = gluNewQuadric();
  ::gluQuadricDrawStyle(quad, GLU_FILL);
  ::glColor3f(1.0f, 1.0f, 0.0f);
  ::gluSphere(quad, r, 20, 20);
  ::gluDeleteQuadric(quad);

  ::glPopMatrix();

  if(!was)
    ::glDisable(GL_LIGHTING);
}

template<typename Kernel>
void gl_draw_cylinder(const typename Kernel::Point_3& p1, 
                      const typename Kernel::Point_3& p2)
{
  typename Kernel::Vector_3 v(p1, p2);
  double r = 0.01 * ((p1 - p2) * (p1 - p2));//SQ_R
  double angle = 0.;
  unsigned int NSIDES = 10;
  double anglestep = 360./(float)NSIDES;
  for(unsigned int i = 0; i < NSIDES; i++) 
  {
    double nextangle = angle + anglestep;
    ::glBegin(GL_QUADS);
      ::glVertex3d(p1.x() + r*std::cos(angle), p1.y() + r*std::sin(angle), p1.z());
      ::glVertex3d(p1.x() + r*std::cos(nextangle), p1.y() + r*std::sin(nextangle), p1.z());
      ::glVertex3d(p2.x() + r*std::cos(nextangle), p2.y() + r*std::sin(nextangle), p2.z());
      ::glVertex3d(p2.x() + r*std::cos(angle), p2.y() + r*std::sin(angle), p2.z());
    ::glEnd();//GL_QUADS
    angle = nextangle;
  }
}

template<typename Kernel>
void gl_draw_cone(const typename Kernel::Point_3& p1,
                  const typename Kernel::Point_3& p2) //center of basis
{  
  double h = 0.1 * std::sqrt((p2 - p1)*(p2 - p1));
  double r = 0.5 * h;//radius of basis
  typename Kernel::Vector_3 v(p1, p2);
  typename Kernel::Point_3 tip = p2 + 0.2 * v;

  double angle = 0.;
  unsigned int NSIDES = 10;
  double anglestep = 360./(float)NSIDES;
  ::glBegin(GL_TRIANGLE_FAN);
  ::glVertex3d(tip.x(), tip.y(), tip.z());
  for(unsigned i = 0; i < NSIDES; i++, angle += anglestep)
    ::glVertex3d(p2.x() + r * std::cos(angle), p2.y() + r * std::sin(angle), p2.z());
  ::glEnd(); //GL_TRIANGLE_FAN

  angle = 0.;
  ::glBegin(GL_POLYGON);
  for(unsigned int i = 0; i < NSIDES; i++, angle += anglestep)
    ::glVertex3d(p2.x() + r * std::cos(angle), p2.y() + r * std::sin(angle), p2.z());
  ::glEnd();//GL_POLYGON
}

template<typename Kernel>
void gl_draw_arrow(const typename Kernel::Point_3& p1, 
                   const typename Kernel::Point_3& p2)
{
  gl_draw_segment<Kernel>(p1,p2);
  //// place and orient arrow as a whole
  //::glPushMatrix();
  //::glRotatef(180., 1., 0., 0.); // 180 degrees around x-axis
  //::glTranslatef(-1., 0., 0.); // move arrow so point is at origin
  //  
  //// scale and orient cylinder part of arrow
  //::glPushMatrix();
  //::glRotatef(90., 0., 1., 0.); // 90 degrees around z-axis
  //::glScalef( 2., 0.5, 0.5);
  gl_draw_cylinder<Kernel>(p1, p2);
  //::glPopMatrix();

  //// now use the cone as defined without any transforms
  
  gl_draw_cone<Kernel>(p1, p2);
  //::glPopMatrix();
}

template<typename K>
bool is_facet_above_plane(const typename K::Plane_3& plane,
                    const typename K::Point_3& pa,
                    const typename K::Point_3& pb,
                    const typename K::Point_3& pc) 
{
  typedef typename K::Oriented_side Side;
  using CGAL::ON_NEGATIVE_SIDE;
  const Side sa = plane.oriented_side(pa);
  const Side sb = plane.oriented_side(pb);
  const Side sc = plane.oriented_side(pc);
  return (sa == ON_NEGATIVE_SIDE && sb == ON_NEGATIVE_SIDE && sc == ON_NEGATIVE_SIDE);
}

template<typename K>
bool is_above_plane(const typename K::Plane_3& plane,
                    const typename K::Point_3& pa,
                    const typename K::Point_3& pb,
                    const typename K::Point_3& pc,
                    const typename K::Point_3& pd) 
{
  typedef typename K::Oriented_side Side;
  using CGAL::ON_NEGATIVE_SIDE;
  const Side sd = plane.oriented_side(pd);
  return (sd == ON_NEGATIVE_SIDE && is_facet_above_plane<K>(plane,pa,pb,pc));
}

template<typename C3T3, typename Plane>
void gl_draw_c3t3(const C3T3& c3t3,
                  const Plane& plane)
{
  ::glLineWidth(2.f);
  typedef typename C3T3::Triangulation Tr;
  typedef typename Tr::Geom_traits K;
  typedef typename Tr::Facet Facet;
  typename C3T3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin();
  for(; cit != c3t3.cells_in_complex_end(); ++cit)
  {
    const typename K::Point_3& p0 = cit->vertex(0)->point();
    const typename K::Point_3& p1 = cit->vertex(1)->point();
    const typename K::Point_3& p2 = cit->vertex(2)->point();
    const typename K::Point_3& p3 = cit->vertex(3)->point();
    if(is_above_plane<K>(plane, p0, p1, p2, p3))
      for(int i = 0; i < 4; i++){
        if(cit->is_facet_on_surface(i)){
          gl_draw_facet<K, Facet>(Facet(cit,i), EDGES_AND_FACES, 44, 117, 255);
        }
        else{
          gl_draw_facet<K, Facet>(Facet(cit,i), EDGES_ONLY, 0, 0, 0);
        }
      }
  }
}

template<typename Colored_polyhedron>
void gl_draw_polyhedron_facet(const typename Colored_polyhedron::Facet_const_iterator& f,
                              double f_color)
{
  typedef typename Colored_polyhedron::Point_3 Point;
  const Point& pa = f->halfedge()->vertex()->point();
  const Point& pb = f->halfedge()->next()->vertex()->point();
  const Point& pc = f->halfedge()->next()->next()->vertex()->point();

  ::glColor3d(f_color, f_color, f_color);
  ::glBegin(GL_TRIANGLES);
  ::glVertex3d(pa.x(),pa.y(),pa.z());
  ::glVertex3d(pb.x(),pb.y(),pb.z());
  ::glVertex3d(pc.x(),pc.y(),pc.z());
  ::glEnd();
}

template<typename Colored_polyhedron, typename Plane>
void gl_draw_colored_poly(const Colored_polyhedron& P,
                          const Plane& plane)
{
  typedef typename Colored_polyhedron::Facet_const_iterator  FCI;

  double strongest_color = -1;
  for( FCI fi = P.facets_begin(); fi != P.facets_end(); ++fi)
    if(fi->color()>strongest_color)
      strongest_color = fi->color();

  GLboolean was = (::glIsEnabled(GL_LIGHTING));
  if(was)
    ::glDisable(GL_LIGHTING);

  for( FCI fi = P.facets_begin(); fi != P.facets_end(); ++fi)
  {
    double f_color = fi->color()/strongest_color*256.;
    if(f_color < 0.)
      continue;
    gl_draw_polyhedron_facet<Colored_polyhedron>(fi, f_color);
  }

  if(was)
    ::glEnable(GL_LIGHTING);
}


#endif //DRAWING_HELPER_H

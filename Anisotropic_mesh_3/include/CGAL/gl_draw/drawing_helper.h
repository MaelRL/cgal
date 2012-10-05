#ifndef DRAWING_HELPER_H
#define DRAWING_HELPER_H

#include <CGAL/gl.h>
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
    ::glColor3f(r / 256., g / 256., b / 256.);
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
    ::glColor3f(0.,0.,0.);
    ::glBegin(GL_LINE_LOOP);
    if(option == EDGES_AND_FACES) ::glLineWidth(2.f);
    else                          ::glLineWidth(1.f);
    ::glVertex3d(pa.x(),pa.y(),pa.z());
    ::glVertex3d(pb.x(),pb.y(),pb.z());
    ::glVertex3d(pc.x(),pc.y(),pc.z());
    ::glEnd();
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
                     const typename Kernel::Point_3& pb) 
{
  ::glBegin(GL_LINES);
  ::glVertex3d(pa.x(),pa.y(),pa.z());
  ::glVertex3d(pb.x(),pb.y(),pb.z());
  ::glEnd();
}

template<typename Kernel>
void gl_draw_sphere(const typename Kernel::Sphere_3& s)
{
  typename Kernel::Point_3 c = s.center();
  const GLdouble r = std::sqrt(s.squared_radius());

  GLint polygon_mode[2];
  ::glGetIntegerv(GL_POLYGON_MODE, &polygon_mode[0]);
  ::glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);

  //  CGALglcolor(Qt::red);
  ::glPolygonMode(GL_FRONT_AND_BACK, polygon_mode[0]);
  ::glPushMatrix();

  ::glTranslated(CGAL::to_double(c.x()), CGAL::to_double(c.y()), CGAL::to_double(c.z()));
  ::glScaled(r, r, r);
  ::glPopMatrix();

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
      ::glVertex3f(p1.x() + r*std::cos(angle), p1.y() + r*std::sin(angle), p1.z());
      ::glVertex3f(p1.x() + r*std::cos(nextangle), p1.y() + r*std::sin(nextangle), p1.z());
      ::glVertex3f(p2.x() + r*std::cos(nextangle), p2.y() + r*std::sin(nextangle), p2.z());
      ::glVertex3f(p2.x() + r*std::cos(angle), p2.y() + r*std::sin(angle), p2.z());
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
  ::glVertex3f(tip.x(), tip.y(), tip.z());
  for(unsigned i = 0; i < NSIDES; i++, angle += anglestep)
    ::glVertex3f(p2.x() + r * std::cos(angle), p2.y() + r * std::sin(angle), p2.z());
  ::glEnd(); //GL_TRIANGLE_FAN

  angle = 0.;
  ::glBegin(GL_POLYGON);
  for(unsigned int i = 0; i < NSIDES; i++, angle += anglestep)
    ::glVertex3f(p2.x() + r * std::cos(angle), p2.y() + r * std::sin(angle), p2.z());
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


#endif //DRAWING_HELPER_H

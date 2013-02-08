//Author: Sebastien Loriot sebastien.loriot@sophia.inria.fr

#ifndef VIEWER_H
#define VIEWER_H

#include <QGLViewer/qglviewer.h>
#include <QGLViewer/camera.h>
#include <QKeyEvent>
#include "tools.h"


#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/Cartesian_converter.h>

#include "Circular_arc_3_subsampling.h"


class Viewer : public QGLViewer
{
  Q_OBJECT
  typedef CGAL::Exact_spherical_kernel_3 SK;
  
  typedef Kernel::Point_3 Point_3;
  typedef std::list< std::list<Point_3> > Subsampled_arcs;

  Point_3 center_;
  double radius_;
  bool draw_balls;
  bool draw_inputs;
  double min_edge_size;
  
  Subsampled_arcs subsampled_arcs;
  Subsampled_arcs contour_arcs;
  std::list<Point_3> inputs;
  
  GLUquadricObj *qsphere;

  bool antialiasing;
  
  template <class Triangulation_on_sphere>
  void build_the_boundary(const Triangulation_on_sphere&);


 
  template <class Edge>
  bool keep_edge(Edge edge) const {
    return edge.first->is_inside || edge.first->neighbor(edge.second)->is_inside;
  }
  
  template <class Edge>
   bool is_boundary_edge(Edge edge) const {
     return edge.first->is_inside != edge.first->neighbor(edge.second)->is_inside;
   } 

protected :
  virtual void draw();
  virtual void init();
  virtual void keyPressEvent(QKeyEvent *e);
public:
  Viewer(QWidget* parent) : QGLViewer(parent) {}

  bool antiAliasing() const { return antialiasing; }

public slots:
  void setAntiAliasing(bool b) 
  {
    antialiasing = b;
    updateGL();
  }

public:
  template <class Triangulation_on_sphere,class Iterator>
  void open(Iterator begin, Iterator end,Triangulation_on_sphere& T,Point_3 center,double scale)
  {
    subsampled_arcs.clear();
    contour_arcs.clear();
    inputs.clear();
    center_ = center;
    radius_ = scale;
    draw_balls = true;
    draw_inputs = false;
    min_edge_size = scale/100.;
    std::copy(begin,end,std::back_inserter(inputs));
    build_the_boundary(T);
  }

  ~Viewer(){
    gluDeleteQuadric(qsphere);
  }
  
};

void Viewer::keyPressEvent(QKeyEvent *e){
  if (e->key()==Qt::Key_P){
    draw_inputs=!draw_inputs;
    updateGL();
    return;
  }

  if (e->key()==Qt::Key_S){
    draw_balls=!draw_balls;
    updateGL();
    return;
  }
  
  QGLViewer::keyPressEvent(e);
}

void Viewer::draw()
{
  glPushMatrix();
  glScalef(radius_,radius_,radius_);
  glTranslatef(-center_.x(),-center_.y(),-center_.z());
  
  if(antiAliasing())
  {
    ::glEnable(GL_BLEND);
    ::glEnable(GL_LINE_SMOOTH);
    ::glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  }
  else
  {
    ::glDisable(GL_BLEND);
    ::glDisable(GL_LINE_SMOOTH);
    ::glDisable(GL_POLYGON_SMOOTH_HINT);
    ::glBlendFunc(GL_ONE, GL_ZERO);
    ::glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);
  }

  if (draw_inputs){
    glPointSize(5);
    glBegin(GL_POINTS);
    glColor3f(0,1,0);
    for (std::list<Point_3>::const_iterator it=inputs.begin();it!=inputs.end();++it)
      //~ drawing_tools<Kernel>::draw_a_point_as_sphere(qsphere,*it,Point_3(0,0,1));
      glVertex3f(it->x(),it->y(),it->z());
    glEnd();
  }
  
  if (draw_balls)
    drawing_tools<Kernel>::drawSphere(qsphere,center_,sqrt(radius_), Point_3(1,1,1));
  
  glDisable(GL_LIGHTING);
  for (Subsampled_arcs::iterator it_arcs=subsampled_arcs.begin();it_arcs!=subsampled_arcs.end();++it_arcs){
      glLineWidth(2);
      glBegin(GL_LINE_STRIP);
      glColor3f(0,0,1);
      for (std::list<Point_3>::iterator it_pt=it_arcs->begin();it_pt!=it_arcs->end();++it_pt)
        glVertex3f(it_pt->x(),it_pt->y(),it_pt->z());
      glEnd();
  }
  
  for (Subsampled_arcs::iterator it_arcs=contour_arcs.begin();it_arcs!=contour_arcs.end();++it_arcs){
      glLineWidth(2);
      glBegin(GL_LINE_STRIP);
      glColor3f(1,0,0);
      for (std::list<Point_3>::iterator it_pt=it_arcs->begin();it_pt!=it_arcs->end();++it_pt)
        glVertex3f(it_pt->x(),it_pt->y(),it_pt->z());
      glEnd();
  }  
  glEnable(GL_LIGHTING);
  glPopMatrix();
}


template <class Triangulation_on_sphere>
void Viewer::build_the_boundary(const Triangulation_on_sphere& T)
{
	
  for (typename Triangulation_on_sphere::All_edges_iterator 
    it=T.all_edges_begin();it!=T.all_edges_end();++it)
  {
    //~ if ( it->first->is_negative() &&  it->first->neighbor(it->second)->is_negative() )
      //~ continue;
    //~ if ( it->first->vertex((it->second+1)%3)->point().index==0 || it->first->vertex((it->second+2)%3)->point().index==0)
      //~ continue;
    if (!keep_edge(*it)) continue;
    
    Point_3 source=it->first->vertex( (it->second+1)%3 )->point();
    Point_3 target=it->first->vertex( (it->second+2)%3 )->point();
    Kernel::Plane_3  plane(source,target,center_);
    Kernel::Circle_3 circle(center_,radius_,plane);    
    
    
    if (is_boundary_edge(*it)){
      contour_arcs.push_front(std::list<Kernel::Point_3>());
      subsample_circular_arc_3<Kernel>(circle,plane,source,target,std::back_inserter(*contour_arcs.begin()),min_edge_size);      
    }
    else{
      subsampled_arcs.push_front(std::list<Kernel::Point_3>());
      subsample_circular_arc_3<Kernel>(circle,plane,source,target,std::back_inserter(*subsampled_arcs.begin()),min_edge_size);
    }
  }
}


void Viewer::init()
{
  // Restore previous viewer state.
  restoreStateFromFile();
 
  qsphere=gluNewQuadric();
  gluQuadricOrientation(qsphere,GLU_OUTSIDE);

  //Jane
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHTING);
  float lpos[4] = { -.2f, .2f, .9797958971f, 0.0f };
  glLightfv(GL_LIGHT0,GL_POSITION,lpos);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0f);

  ::glEnable(GL_NORMALIZE);
  glShadeModel(GL_SMOOTH);
}

#include "Viewer.moc"
#endif


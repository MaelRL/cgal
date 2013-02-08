//Author: Sebastien Loriot sebastien.loriot@sophia.inria.fr
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Delaunay_triangulation_sphere_traits_2.h>
#include <CGAL/Delaunay_triangulation_sphere_2.h>
#include <CGAL/utility.h>

#include <CGAL/Qt/DemosMainWindow.h>
#include <QMainWindow>
#include <QFileDialog>
#include <QApplication>
#include <QProgressDialog>

#include <list>
#include <fstream>
#include <boost/progress.hpp>

#include "Real_regular_triangulation_sphere_traits_2.h"
#include "Regular_triangulation_face_base_with_bool_on_sphere_2.h"

typedef CGAL::Exact_predicates_exact_constructions_kernel         EPEC;

typedef CGAL::Exact_predicates_inexact_constructions_kernel         Kernel;

// must be included once "Kernel" is defined
#include "Viewer.h"
#include "ui_Mainwindow.h"

typedef CGAL::Regular_triangulation_euclidean_traits_3<Kernel>      Traits;

typedef Traits::Weighted_point_3                    Weighted_point_3;
typedef Kernel::Point_3                             Point_3;

template <class Point_>
class Point_with_id:public Point_{
public:

  unsigned face;
  unsigned index;
  Point_with_id():Point_(){}
  Point_with_id(const Point_& p,unsigned f,unsigned i):Point_(p),face(f),index(i){} // AF: Windows doesn't like Point
};



typedef CGAL::Real_regular_triangulation_sphere_traits_2<Traits,Point_with_id<Weighted_point_3> >       Gt;
typedef  CGAL::Triangulation_data_structure_2 <CGAL::Regular_triangulation_vertex_base_2<Gt>,
                                               CGAL::Regular_triangulation_face_base_with_bool_on_sphere_2<Gt> > TDS;
typedef CGAL::Delaunay_triangulation_sphere_2<Gt,TDS>                                     RTOS;


void read_polygons(const char* file_path,std::list<std::vector<Point_3> >& polygons){
  double Long,lat;
  std::ifstream input(file_path);
  if (!input){
    std::cerr << "Pb while opening file\n";
    exit(EXIT_FAILURE);
  }
  
  while(!input.eof()){
    input >> Long;
    if (!input) {
      input.clear();
      input.ignore(10,'\n'); 
      if ( input.eof() ) break;
      polygons.push_front(std::vector<Point_3>());
      continue;
    }
    input >> lat;
    //convert into radian
    Long=Long/180.*M_PI;
    lat=lat/180.*M_PI;
    polygons.begin()->push_back(
      Point_3(
        cos(Long) * cos (lat),
        sin(Long) * cos (lat),
        sin(lat)    
      )
    );
  }
}

typedef Point_with_id<Weighted_point_3> Tos_point;

template<class Face>
bool is_an_inside_face(Face f) {
  const Tos_point& a=f.vertex(0)->point();
  const Tos_point& b=f.vertex(1)->point();
  const Tos_point& c=f.vertex(2)->point();
  
  //handle fake points
  if (a.index==0 || c.index==0 || b.index==0 ) return false;
  
  
  if ( a.face!=b.face || b.face!=c.face ) return false;
  if (a.index < b.index) return !( b.index < c.index || c.index < a.index );
  return !( b.index < c.index &&  a.index > c.index );
}

template <class Face_handle>
double max_edge_size(Face_handle f){
  double d=CGAL::squared_distance(f->vertex(0)->point(),f->vertex(1)->point());
  double max=d;
  d=CGAL::squared_distance(f->vertex(0)->point(),f->vertex(2)->point());
  if (d>max) max=d;
  d=CGAL::squared_distance(f->vertex(1)->point(),f->vertex(2)->point());
  if (d>max) max=d;
  return max;
}


template<class Face_handle>
Tos_point barycenter(Face_handle f){
  const Point_3& a=f->vertex(0)->point().point();
  const Point_3& b=f->vertex(1)->point().point();
  const Point_3& c=f->vertex(2)->point().point();
  EPEC::Point_3 center(0,0,0);
  
  CGAL::Cartesian_converter<Kernel,EPEC> in_exact;
  CGAL::Cartesian_converter<EPEC,Kernel> in_double;
  
  EPEC::Point_3 p=center+ (in_exact(a)-center)/3. + (in_exact(b)-center)/3. +(in_exact(c)-center)/3.;
  
  return Tos_point(Weighted_point_3(in_double(p),0),0,0);
  
}

class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::MainWindow
{
  Q_OBJECT
public:
  MainWindow() {
    setupUi(this);

    connect(this->actionAntiAliasing, SIGNAL(toggled(bool)),
            viewer, SLOT(setAntiAliasing(bool)));
  }
public slots:
  void open(QString filename) {

    std::list<Traits::Weighted_point_3> lst_pt;
    std::list<std::vector<Point_3> > polygons;
  
    read_polygons(filename.toUtf8(),polygons);
  
    std::list<Tos_point> tr_pts;
  
    unsigned face=1;
    for (std::list<std::vector<Point_3> >::const_iterator it=polygons.begin();it!=polygons.end();++it){
      unsigned index=0;
      for (std::vector<Point_3>::const_iterator itpt=it->begin();itpt!=it->end();++itpt){
        //TODO sample arc here
        tr_pts.push_back(Tos_point(Weighted_point_3(*itpt,9e-10),face,++index));
      }
      ++face;
    }
    assert(polygons.size()==face-1);
    std::cerr << polygons.size() << " polygons" << std::endl;
  

    RTOS tr=RTOS( Gt(Tos_point(Weighted_point_3(Point_3(0,0,0),1),0,0)) );
  
    tr.insert(Tos_point(Weighted_point_3(Point_3(0.1,0,0),-1),0,0));
    tr.insert(Tos_point(Weighted_point_3(Point_3(0,0.1,0),-1),0,0));
    tr.insert(Tos_point(Weighted_point_3(Point_3(0,0,0.1),-1),0,0));
    tr.insert(Tos_point(Weighted_point_3(Point_3(-sqrt(0.1/3.0),-sqrt(0.1/3.0),-sqrt(0.1/3.0)),-1),0,0));
  
    tr.insert(tr_pts.begin(),tr_pts.end());

    for (RTOS::All_faces_iterator it=tr.all_faces_begin();it!=tr.all_faces_end();++it)
    {
      if ( it->is_ghost() ) it->is_inside=false;
      else it->is_inside=is_an_inside_face(*it);
    }
  
    //brute force refining
    Kernel::Compute_squared_radius_3 sqrad;
    Kernel::Construct_circumcenter_3 circumcenter;

    Point_3 center(0.,0.,0.);
    double scale=1;
  
    unsigned nb_step=3;
    std::cout << "Refining the triangulation" << std::endl;
    boost::progress_display show_progress(nb_step-1);
    QProgressDialog progress_dialog(QObject::tr("Refining the triangulation..."),
                                    QObject::tr("Abort"),
                                    0,
                                    nb_step-1,
                                    this);
    progress_dialog.setWindowModality(Qt::WindowModal);
    while (--nb_step!=0)
    {
      double max=0;
      RTOS::Face_handle f;
      for (RTOS::All_faces_iterator it=tr.all_faces_begin();it!=tr.all_faces_end();++it)
      {
        if (!it->is_inside) continue;
      
        //~ double sqr=sqrad(it->vertex(0)->point(),it->vertex(1)->point(),it->vertex(2)->point());
        double sqr=max_edge_size(it);
        if (sqr > max){
          max=sqr;
          f=it;
        }
      }
      assert(f!=RTOS::Face_handle());
      //~ Point_3 p=circumcenter(f->vertex(0)->point(),f->vertex(1)->point(),f->vertex(2)->point());
      //~ Point_3 p=center+ (f->vertex(0)->point()-center)/3.+(f->vertex(1)->point()-center)/3.+(f->vertex(2)->point()-center)/3.;
      Point_3 p=barycenter(f);
      p=center+(p-center)/sqrt(CGAL::squared_distance(p,center));

      tr.insert(Tos_point(Weighted_point_3(p,0),0,0));
      lst_pt.push_back(Weighted_point_3(p,0));
      ++show_progress;
      progress_dialog.setValue(show_progress.count());
    }
  

  
  
  // Instantiate the viewer.
    viewer->open(lst_pt.begin(),lst_pt.end(),tr,center,scale);

  }

  void on_action_Quit_triggered() {
    close();
  }

  void on_action_Open_triggered() {
    QString filename = QFileDialog::getOpenFileName(this);
    if(!filename.isNull())
      open(filename);
  }
};

int main(int argc, char** argv)
{
  // Read command lines arguments.
  QApplication application(argc,argv);

  QStringList args = QApplication::arguments();
  args.removeAt(0);

  MainWindow mainWindow;
  mainWindow.show();
  if(!args.empty())
  {
    mainWindow.open(args[0]);
  }  

  // Run main loop.
  return application.exec();
}

#include "Mainwindow.moc"

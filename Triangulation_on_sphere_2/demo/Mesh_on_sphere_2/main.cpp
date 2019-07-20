//Author: Sebastien Loriot sebastien.loriot@sophia.inria.fr
#include <boost/iterator/transform_iterator.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/Delaunay_triangulation_sphere_traits_2.h>
#include <CGAL/Projection_sphere_traits_3.h>
#include <CGAL/Delaunay_triangulation_sphere_2.h>
#include <CGAL/utility.h>
#include <CGAL/Qt/DemosMainWindow.h>
#include <QMainWindow>
#include <QFileDialog>
#include <list>
#include <fstream>
#include<CGAL/Constrained_Delaunay_triangulation_sphere_2.h>
#include<CGAL/Delaunay_mesher_2.h>
#include<CGAL/Delaunay_mesh_face_base_2.h>
#include<CGAL/Delaunay_mesh_sphere_size_criteria_2.h>
#include<CGAL/Constrained_triangulation_face_base_sphere_2.h>
#include<CGAL/Delaunay_triangulation_sphere_traits_2.h>
#include<CGAL/Delaunay_mesh_sphere_traits_2.h>
#include <qapplication.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         Kernel;

#include "Viewer.h"
#include "ui_Mainwindow.h"

typedef CGAL::Delaunay_triangulation_sphere_traits_2<Kernel>       Gt;
typedef CGAL::Delaunay_mesh_sphere_traits_2<Gt>	                   Mgt;
typedef CGAL::Constrained_triangulation_face_base_sphere_2<Mgt> Cfb;
typedef CGAL::Triangulation_vertex_base_2<Mgt> Vb;

typedef CGAL::Delaunay_mesh_face_base_2<Mgt, Cfb> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;

typedef CGAL::Constrained_Delaunay_triangulation_sphere_2<Mgt, Tds> CDTS;

typedef CGAL::Delaunay_mesh_sphere_size_criteria_2<CDTS> Criteria;

struct Cell_info{
	Kernel::Point_3* pt;
	Cell_info():pt(NULL){}
	~Cell_info(){if (pt!=NULL) delete pt;}
};

typedef Kernel::Point_3                             Point_3;
typedef CDTS::Vertex_handle      Vertex_handle;

template <class Output_iterator>
void read_points(const char* file_path,Output_iterator out){
	int nb;
	double Long, lat;
	
	std::ifstream input(file_path);
	if (!input){
		std::cerr << "Error while reading " << file_path << std::endl;
		exit(EXIT_FAILURE);
	}
	
	input >> nb;
	
	for (int i=0;i<nb;++i){
#warning tmp : must handle the sphere on which are points (parameter of the traits?)
		input >> Long;
		input >> lat;
		Long=Long/180.*M_PI;
		lat=lat/180.*M_PI;
		*out++=	Point_3(100*cos(Long) * cos (lat),	100*sin(Long) * cos (lat),100*sin(lat));    
		
		
	}
}



class MainWindow :
public CGAL::Qt::DemosMainWindow,
public Ui::MainWindow
{
	Q_OBJECT
public:
	MainWindow() {
		setupUi(this);
	}
	public slots:
	void open(QString filename) {
		CDTS cdt;
		
		std::vector<Point_3> lst_pt;
		read_points(filename.toUtf8().data(),
					std::back_inserter(lst_pt));
		
		Vertex_handle v;
		cdt.set_radius(100);
		std::vector<Vertex_handle> vertices;
		for(int i=0;i< lst_pt.size(); i++){
			Point_3 p = lst_pt.at(i);
			double x = p.x();
			double y= p.y();
			double z = p.z();
			//std::cout<<lst_pt.at(i)<<std::endl;
			v =cdt.insert(lst_pt.at(i));
			vertices.push_back(v);
		}
		int n=lst_pt.size();
		cdt.insert_constraint(vertices.at(n-1), vertices.at(0));
		
		for (int i=1; i<n; i++)
			cdt.insert_constraint(vertices.at(i), vertices.at(i-1));
		
				
		Point_3 center(0,0,0);
		double scale=1;
		CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125,3,cdt.geom_traits()),false);
		int test = cdt.number_of_vertices();
		//copy points for drawing
		std::vector<Point_3> points;
		CDTS::All_vertices_iterator vit;
		vit = cdt.vertices_begin();
		for(;vit!=cdt.vertices_end();vit++)
			points.push_back(vit->point());
		
		MainWindow mainWindow;
		mainWindow.show();
		
		// Instantiate the viewer.
		//viewer->open(lst_pt.begin(),lst_pt.end(),cdt,center,scale);
		viewer->open(points.begin(),points.end(),cdt,center,scale);
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
	if(!args.empty())
	{
		mainWindow.open(args[0]);
	}  
	mainWindow.show();
	
	// Run main loop.
	return application.exec();
}

#include "Mainwindow.moc"

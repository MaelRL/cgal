//Only under unix like machines


//#define NDEBUG
//#define DEBUG_CONFLICT_VOLUME_LISTS
//#define DEBUG_CONFLICT_LISTS
//#define DEBUG_PRINT_ALL_UMBRELLAS
//#define DEBUG_READING_FROM_FILE
//#define DEBUG_CONFLICT_RADIUS
//#define DEBUG_CONFLICTS
//#define DEBUG_VOLUME_REFINEMENT_POINT
//#define DEBUG_SURFACE
//#define DEBUG_NEIGHBORS

#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
#include <signal.h>
#include <stdlib.h>

#include <sstream>
#include <string>
#include <iostream>
#include <fstream>

#include <CGAL/basic.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Mesher_level_visitors.h>

//====================extended=========================
#include <CGAL/Stretched_traits_3.h>
#include <CGAL/Stretched_delaunay_3.h>
#include <CGAL/Umbrella_set_3.h>
#include <CGAL/Refine_levels.h>
#include <CGAL/IO/Umbrella_to_medit.h>
#include <CGAL/IO/Output_umbrellaset.h>
#include <CGAL/IO/Input_umbrellaset.h>
#include <CGAL/All_visitors.h>

//====================local files==========================
#include "Circular_metric_field_3.h"
#include "Torus_surface.h"

using namespace std;

typedef CGAL::Exact_predicates_inexact_constructions_kernel					K;
typedef K::Point_3															Point;
typedef Torus_surface<K>													Constrain_surface;
typedef Umbrella_set<K,Circular_metric_field_3<Point>,Constrain_surface >	Stretched_torus;
typedef Surface_criteria_for_umbrella_set<Stretched_torus,Constrain_surface>	  		  		Torus_surface_criteria;

void print_umbrella_to_medit(unsigned index, Stretched_torus &schetched_torus,
		double main_radius,double minor_radius,double stretch,double max_ratio,
		double max_distortion,double max_size,double max_distance);

void print_umbrella_volume_to_medit(unsigned index, Stretched_torus &stretched_torus);

void Record_data(Stretched_torus &stretched_torus,char *dirname,
		double main_radius,double minor_radius,double stretch,double max_ratio,
		double max_distortion,double max_size,double max_distance);

void Anisotropic_mesher(Stretched_torus &stretched_torus,char *dirname=NULL);
void Refine_new_mesh();
bool Recover_mesh(char *dirname);
bool init_stretched_torus_from_file(Stretched_torus &stretched_torus, char* dirname);

/** The shape of the torus*/
static double main_radius = 5.;
static double minor_radius = 1.;

/** The stretch of the circular metric*/
static double stretch = .64;

/** The quality of the surface triangulation*/
static double max_ratio = 3.;
static double max_distortion = 10.;
static double max_size = 1.;
static double max_distance = 0.06;

static bool kill_signal = false;


void mySig(int sigNo)
{
	if(sigNo == SIGINT || SIGQUIT){
		cout<<"\nKill signal received. Print result to files before quitting...\n";
		kill_signal = true;
	}
}

int main(int argc,char *argv[]){
  if(argc <= 1){
	  Refine_new_mesh();
  }else{
	  Recover_mesh(argv[1]);
  }
  return 0;
}

bool Recover_mesh(char *dirname){
	string filename(dirname);
	filename += string("/ToreMeshQuality");
	ifstream input_file(filename.c_str());
	if(input_file.fail()){
		cout<<"File: "<<filename<<" not found."<<endl;
		return false;
	}
	input_file>>main_radius>>minor_radius>>stretch>>max_ratio>>max_distortion>>max_size>>max_distance;
	if(input_file.fail()){
		cout<<"Error in reading "<<filename<<endl;
		input_file.close();
		return false;
	}
	input_file.close();

	Stretched_torus	  stretched_torus(Circular_metric_field_3<Point>(stretch),
			Constrain_surface(main_radius,minor_radius),
			max_ratio,max_distortion,max_size,max_distance,false);
	
	if(init_stretched_torus_from_file(stretched_torus,dirname)==false){	
		cout<<"Recover failure..."<<endl;
		return false;
	}

	Anisotropic_mesher(stretched_torus,dirname);
	return true;
}

	
void Refine_new_mesh(){

	Stretched_torus	  stretched_torus(Circular_metric_field_3<Point>(stretch),
			Constrain_surface(main_radius,minor_radius),
			max_ratio,max_distortion,max_size,max_distance);

	Anisotropic_mesher(stretched_torus,NULL);
}


void Anisotropic_mesher(Stretched_torus &stretched_torus,char *dirname){
	
	struct sigaction act;
	struct sigaction old;
	act.sa_handler = mySig;
	
	CGAL::Null_mesher_level		  null_level;

	typedef Surface_mesher_level<Stretched_torus>					Surface_refiner;
	typedef Volume_mesher_level<Stretched_torus,Surface_refiner>	Volume_refiner;
	Surface_refiner					surface_refiner(stretched_torus,null_level);
	Volume_refiner					volume_refiner(stretched_torus,surface_refiner);	

	Surface_conflicts_visitor<Stretched_torus>		surface_refinement_visitor(stretched_torus);
	Volume_conflicts_visitor<Stretched_torus>		volume_refinement_visitor(surface_refinement_visitor);

	bool not_yet_finish=true;
	surface_refiner.scan_triangulation();
	volume_refiner.scan_triangulation();
	cout<<"\nTriangulation initialized with "<<stretched_torus.size()<<" points."<<endl;
	cout<<"\n \
		* Starting to refine the mesh.                             *\n \
		* Press Ctrl + c for interrupt and quit.                   *\n \
		* Program will store the triangulation automatically.      *\n \
		* Default directory to store data is ..\\torus\\ .           *\n \
		* For continue the refinement please use command:          *\n\
		* ./<program_name> <dir_name>                              *"<<endl;
	sigaction(SIGINT,&act,&old);
	sigaction(SIGQUIT,&act,&old);
	
	if(dirname == NULL){
		dirname = new char[6];
		strcpy(dirname,"torus");
	}
	struct stat dir_stat;
	stat(dirname, &dir_stat);
	if(!S_ISDIR(dir_stat.st_mode)){
		string cmd("mkdir "); 
		cmd += string(dirname);
		system(cmd.c_str());
	}
	
	int total_points=stretched_torus.size();
	for(;not_yet_finish && (!kill_signal);total_points++){
		if((total_points%50)==49){
			stretched_torus.update_all_surfaces();
			print_umbrella_to_medit(total_points,stretched_torus,main_radius,
					minor_radius,stretch,max_ratio,max_distortion,max_size,max_distance);
			stretched_torus.update_all_tetrahedra();
			print_umbrella_volume_to_medit(total_points,stretched_torus);
			Record_data(stretched_torus,dirname,main_radius,minor_radius,stretch,
					max_ratio,max_distortion,max_size,max_distance);
			cout<<"\nReaching the "<<total_points+1<<"th point";
		}
		cout<<".";
		cout.flush();
		not_yet_finish = volume_refiner.try_to_insert_one_point(volume_refinement_visitor);
	}
	total_points = stretched_torus.size();
	stretched_torus.update_all_surfaces();
	print_umbrella_to_medit(total_points,stretched_torus,main_radius,
			minor_radius,stretch,max_ratio,max_distortion,max_size,max_distance);
	stretched_torus.update_all_tetrahedra();
	print_umbrella_volume_to_medit(total_points,stretched_torus);
	Record_data(stretched_torus,dirname,main_radius,minor_radius,stretch,
			max_ratio,max_distortion,max_size,max_distance);
	if(!kill_signal){
		cout<<"\n\nTotal pointss:"<<total_points<<endl;
	}
}

void print_umbrella_volume_to_medit(unsigned index, Stretched_torus &stretched_torus){
	ofstream of;
	stringstream ss;
	ss<<"data/volume_"<<index<<".mesh";
	of.open(ss.str().c_str());
	Umbrella_set_volume_to_medit(stretched_torus,of);
	of.close();
}

void print_umbrella_to_medit(unsigned index, Stretched_torus &schetched_torus,
		double main_radius,double minor_radius,double stretch,
		double max_ratio,double max_distortion,double max_size,double max_distance){
	ofstream of;
	stringstream ss;
	ss<<"data/step"<<index<<".mesh";
	of.open(ss.str().c_str());
	of<< "\n# Torus parameter: \n# main_radius = "<<main_radius
		<<"\n# minor_radius = "<<minor_radius
		<<"\n\n# Stretch = "<<stretch
		<<"\n\n# Surface quality: \n# max_ratio = "<<max_ratio
		<<"\n# max_distortion = "<<max_distortion
		<<"\n# max_size = "<<max_size
		<<"\n# max_distance = "<<max_distance<<"\n\n";
	Umbrella_set_surface_to_medit(schetched_torus,of);
	of.close();
}

bool init_stretched_torus_from_file(Stretched_torus &stretched_torus,char *dirname){
	string filename(dirname);
	filename += string("/AllUmbs");
	ifstream input_file(filename.c_str());
	if(input_file.fail()){
		cerr<<"File: "<<filename<<" not found."<<endl;
		return false;
	}
	bool init_successful = init_umbrella_set_from_file(stretched_torus,input_file);

	input_file.close();
	return init_successful;
}


void Record_data(Stretched_torus &stretched_torus,char *dirname, double main_radius,double minor_radius,
		double stretch, double max_ratio,double max_distortion,double max_size, double max_distance){
	string filename(dirname);
	filename += string("/ToreMeshQuality");
	ofstream ofile(filename.c_str());
	if(!ofile.fail()){
		ofile<<endl<<main_radius<<endl<<minor_radius<<endl<<stretch<<endl
			<<max_ratio<<endl<<max_distortion<<endl<<max_size<<endl<<max_distance;
	}
	ofile.close();
	filename = string(dirname);
	filename += string("/AllUmbs");
	ofile.open(filename.c_str());
	if(!ofile.fail())
		output_umbrella_set(stretched_torus,ofile);
	ofile.close();
}

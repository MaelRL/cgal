#ifndef OUTPUT_UMBRELLASET_H
#define OUTPUT_UMBRELLASET_H

/** Reads a Umbrella_set type structure and
 * outputs a mesh file to the assigned out stream
 */

#include <iostream>
#include <map>
#include <iterator>
#include <vector>
#include <ctime>

#include <CGAL/Conflict_surface_triangles.h>///TODO to nest it in the umbrella_set

template <typename Umbrella>
void output_umbrella(Umbrella *umbrella, std::ostream &os);

/// Convert umbrella into Medit format
///
template<class Umbrella_set>
void output_umbrella_set(Umbrella_set &umbs,std::ostream &os)
{
	typedef typename Umbrella_set::Cell						  Cell;
	typedef typename Umbrella_set::Vertex_handle			  Vertex_handle;
	typedef typename Umbrella_set::Finite_cells_iterator	  Finite_cells_iterator;
	typedef typename Umbrella_set::Finite_facets_iterator	  Finite_facets_iterator;
	typedef typename Umbrella_set::Cell_handle				  Cell_handle;
	typedef typename Umbrella_set::Point					  Point;
	typedef typename Umbrella_set::Umbrella_iterator		  Umbrella_iterator;
	//typedef typename Umbrella_set::Triangle_list			  Triangle_list;

	os<<umbs.size()<<"\n\n";

	for(unsigned i=0 ; i<umbs.size(); i++){
		Point p = umbs.get_point(i);
		os << p.x()<<" "<<p.y()<<" "<<p.z();
		if(umbs.is_on_constraint(i)){
			os<<"   "<< 0 <<"\n";
		}else{
			os<<"   "<< 1 <<"\n";
		}
	}

	os<<"\n\n";

	for(unsigned i=0; i<umbs.size(); i++){
		output_umbrella(umbs.all_umbrella[i],os);
	}		

}

template <typename Umbrella>
void output_umbrella(Umbrella *umbrella, std::ostream &os)
{
	typedef typename Umbrella::Cell						  Cell;
	typedef typename Umbrella::Vertex_handle			  Vertex_handle;
	typedef typename Umbrella::Cell_handle				  Cell_handle;
	typedef typename Umbrella::Point					  Point;
	typedef typename Umbrella::Vertex_map				  Vertex_map;
	typedef typename Umbrella::Neighbor_bounds_iterator	  Neighbor_bounds_iterator;
	
	umbrella->update_neighbor_vertices();
	os<<umbrella->neighbor_vertex_map.size()<<"    ";
	for(typename Vertex_map::iterator vit = umbrella->neighbor_vertex_map.begin();
		vit != umbrella->neighbor_vertex_map.end(); vit ++){
		os<<(vit->second)->info()<<" ";
	}
	os<<"\n";
}


#endif //OUTPUT_UMBRELLASET_H

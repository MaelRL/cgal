#ifndef UMBRELLA_TO_MEDIT_H
#define UMBRELLA_TO_MEDIT_H

/** Reads a Umbrella_set type structure and
 * outputs a mesh file to the assigned out stream
 */

#include <iostream>
#include <map>
#include <iterator>
#include <vector>
#include <ctime>

#include <CGAL/Conflict_surface_triangles.h>///TODO to nest it in the umbrella_set
#include <CGAL/Conflict_volume_tetrahedra.h>

void time_stamp(std::ostream &os);	  //add time information into the file


/// Convert umbrella into Medit format
///
template<typename Umbrella>
void umbrella_to_medit(const Umbrella *const tr, std::ostream& os)
{
	typedef typename Umbrella::Cell					  Cell;
	typedef typename Umbrella::Vertex_handle			  Vertex_handle;
	typedef typename Umbrella::Finite_cells_iterator	  Finite_cells_iterator;
	typedef typename Umbrella::Finite_facets_iterator  Finite_facets_iterator;
	typedef typename Umbrella::Cell_handle			  Cell_handle;
	typedef typename Umbrella::Point					  Point;
	typedef typename Umbrella::Finite_vertices_iterator
															  Finite_vertices_iterator;

	std::map<unsigned, Point>	  vertex_map;
	unsigned max_index=0;

	time_stamp(os);

	for(Finite_vertices_iterator vit = tr->finite_vertices_begin();
		vit != tr->finite_vertices_end();
		vit ++){
		unsigned cur_index = vit->info();
		vertex_map[cur_index] = vit->point();
		if( cur_index >max_index){
		  max_index = cur_index;
		}
	}
	
	os << "MeshVersionFormatted 1\n"
		<< "Dimension\n"
		<< "3 \n\n"
		<< "Vertices\n"
		<< max_index+1 << " \n";


	int tetrahedra_counter = 0;
	for(Finite_cells_iterator cit = tr->finite_cells_begin();
			cit != tr->finite_cells_end();
			++cit){
		if(cit->is_inside(*tr,tr->get_constrain_surface(),tr->get_traits())){//inside of the domain
			tetrahedra_counter++;
		}
	}//for

	for(unsigned i=0; i<= max_index; i++)
	{
		if(vertex_map.find(i)==vertex_map.end()){
		  os <<"0 0 0 0\n";
		}else{
		  Point p = vertex_map[i];
		  os <<  p.x() << " " << p.y() << " " << p.z() << " " << 0  << 
			" \n";
		}
	}

	int triangle_counter = 0;
	for(Finite_cells_iterator cit = tr->finite_cells_begin();
			cit != tr->finite_cells_end();
			++cit)
	{
		if(cit->is_inside(*tr,tr->get_constrain_surface(),tr->get_traits())){//inside of the domain
			for(int i=0;i<4;++i){
				if(!cit->neighbor(i)->is_inside(
						*tr,tr->get_constrain_surface(),tr->get_traits())){//inside of the domain
					triangle_counter++;
				}
			}
		}
	}

	os<<"\nTriangles\n"<<triangle_counter<<"\n";

	for(Finite_cells_iterator cit = tr->finite_cells_begin();
			cit != tr->finite_cells_end();
			++cit)
	{
		if(cit->is_inside(*tr,tr->get_constrain_surface(),tr->get_traits())){//inside the domain
			for(int i=0;i<4;++i){
				if(!cit->neighbor(i)->is_inside(*tr,tr->get_constrain_surface(),tr->get_traits())){
					int color = 0;
					for(int j=0; j<4; ++j){
						if(i!=j){
							unsigned n=cit->vertex(j)->info()+1;
							//color the surface connecting to center
							if( n == tr->center_index()+1){
								color = 1;
							}
							os<< n << " ";
						}
					}
					os<< color << " \n ";
				}
			}
		}//os << color_set << " \n"; // without color. 
	}


	os << "\nTetrahedra\n"<< tetrahedra_counter << " \n";

	for(Finite_cells_iterator cit = tr->finite_cells_begin();
			cit != tr->finite_cells_end();
			++cit)
	{
		if(cit->is_inside(*tr,tr->get_constrain_surface(),tr->get_traits())){//inside of the domain
			int color=0;
			for(int i=0; i<4; ++i){
				unsigned n=cit->vertex(i)->info()+1;
				if(n==tr->center_index()+1){
					color=1;}
				os << n << " ";
			}
			os << color << " \n"; // without color. 
		}//os << color_set << " \n"; // without color. 
	}

	os << "\n" << "End\n";
}

template<class Umbrella_set>
void Umbrella_set_volume_to_medit(Umbrella_set &umbrella_set,std::ostream &os)
{
	typedef typename Umbrella_set::Cell						  Cell;
	typedef typename Umbrella_set::Vertex_handle			  Vertex_handle;
	typedef typename Umbrella_set::Finite_cells_iterator	  Finite_cells_iterator;
	typedef typename Umbrella_set::Finite_facets_iterator	  Finite_facets_iterator;
	typedef typename Umbrella_set::Cell_handle				  Cell_handle;
	typedef typename Umbrella_set::Point					  Point;
	typedef typename Umbrella_set::Umbrella_iterator		  Umbrella_iterator;
	typedef typename Tetrahedra_list::Tetrahedra_indice_iterator
														Tetrahedra_indice_iterator;

	time_stamp( os);

	os << "MeshVersionFormatted 1\n"
		<< "Dimension\n"
		<< "3 \n\n"
		<< "Vertices\n"
		<< umbrella_set.size() << " \n";

	for(unsigned i=0;i<umbrella_set.size();i++){
		Point p = umbrella_set.get_point(i);
		os <<  p.x() << " " << p.y() << " " << p.z() << " " << 0  << " \n";
	}
	
	os << "\nTriangles\n"
		<< 4 * umbrella_set.tetrahedra_list.size() << " \n";
	
	for(Tetrahedra_indice_iterator tit = umbrella_set.tetrahedra_list.begin();
			tit != umbrella_set.tetrahedra_list.end();
			tit ++){
	
		for(int i=0; i<4; i++){
			for(int j=0; j<4; j++){
				if(i!=j){
					os << tit->index(j)+1 << " ";
				}
			}
			os << 0 << "\n";
		}
	}


	os << "\n\nTetrahedra\n"
		<< umbrella_set.tetrahedra_list.size() << " \n";
	
	for(Tetrahedra_indice_iterator tit = umbrella_set.tetrahedra_list.begin();
			tit != umbrella_set.tetrahedra_list.end();
			tit ++){
	
		for(int i=0; i<4; i++){
			os << tit->index(i)+1 << " ";
		}
		os << 0 << "\n";
	}


	os << "\n" << "End\n";

}

template<class Umbrella_set>
void Umbrella_set_surface_to_medit(Umbrella_set &umbrella_set,std::ostream &os)
{
	typedef typename Umbrella_set::Cell						  Cell;
	typedef typename Umbrella_set::Vertex_handle			  Vertex_handle;
	typedef typename Umbrella_set::Finite_cells_iterator	  Finite_cells_iterator;
	typedef typename Umbrella_set::Finite_facets_iterator	  Finite_facets_iterator;
	typedef typename Umbrella_set::Cell_handle				  Cell_handle;
	typedef typename Umbrella_set::Point					  Point;
	typedef typename Umbrella_set::Umbrella_iterator		  Umbrella_iterator;

	time_stamp( os);

	os << "MeshVersionFormatted 1\n"
		<< "Dimension\n"
		<< "3 \n\n"
		<< "Vertices\n"
		<< umbrella_set.size() << " \n";

	for(unsigned i=0;i<umbrella_set.size();i++){
		Point p = umbrella_set.get_point(i);
		os <<  p.x() << " " << p.y() << " " << p.z() << " " << 0  << " \n";
	}
	
	os << "\nTriangles\n"
		<< umbrella_set.surface_list.size() << " \n";
	
	for(typename Surface_list::iterator tit = umbrella_set.surface_list.begin();
			tit != umbrella_set.surface_list.end();
			tit ++){
	
		for(int i=0; i<3; ++i){
			os << tit->index(i)+1 << " ";
		}
		int color_coherency = 0;
		if (tit->is_coherent()){
			color_coherency = 1;
		}
		os << color_coherency << " \n"; // without color. 
	}

	os << "\n" << "End\n";

}



void time_stamp(std::ostream &os){
	time_t raw_time = time(NULL);

	os<<"\n# Generated at "<< ctime(&raw_time) ;
	os<<"\n#    by Anisotropic_mesh_3 function. \n";
	os<<"\n#############################################\n\n\n";

}



template <typename Umbrella>
void umbrella_neighbor_surface_to_medit(Umbrella *umbrella, std::ostream &os)
{
	typedef typename Umbrella::Cell						  Cell;
	typedef typename Umbrella::Vertex_handle			  Vertex_handle;
	typedef typename Umbrella::Cell_handle				  Cell_handle;
	typedef typename Umbrella::Point					  Point;
	typedef typename Umbrella::Vertex_map				  Vertex_map;
	typedef typename Umbrella::Neighbor_bounds_iterator	  Neighbor_bounds_iterator;
	
	umbrella->update_neighbor_bounds();
	umbrella->update_neighbor_vertices();
	unsigned max_index=0;
	for(typename Vertex_map::iterator vit = umbrella->neighbor_vertex_map.begin();
		vit != umbrella->neighbor_vertex_map.end(); vit ++){
		unsigned index = (vit->second)->info();
		if(index > max_index){
		  max_index = index;
		}
	}

	time_stamp(os);
	
	os << "MeshVersionFormatted 1\n"
		<< "Dimension\n"
		<< "3 \n\n"
		<< "Vertices\n"
		<< max_index+1 << " \n";


	for(unsigned i=0;i<=max_index;i++){
		if(umbrella->neighbor_vertex_map.find(i)!=umbrella->neighbor_vertex_map.end()){
		  Point p = umbrella->neighbor_vertex_map[i]->point();
		  os <<  p.x() << " " << p.y() << " " << p.z() << " " << 0  << " \n";
		}else{
		  os<< "0 0 0 0 \n";
		}
	}
	
	os << "\nTriangles\n"
		<< umbrella->neighbor_bounds.size() << " \n";
	
	for(Neighbor_bounds_iterator  nbit = umbrella->neighbor_bounds.begin();
			nbit != umbrella->neighbor_bounds.end();
			nbit ++){
		Cell_handle cell=nbit->first;
		for(int i=0; i<4; ++i){
		  if(i!=int(nbit->second)){
			os << (cell->vertex(i))->info()+1 << " ";
		  }
		}
		os << 0 << " \n"; // without color. 
	}

	os << "\n" << "End\n";
}


#endif //UMBRELLASET_TO_MEDIT_H

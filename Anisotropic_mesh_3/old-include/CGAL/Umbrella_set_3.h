
#ifndef UMBRELLA_SET_H
#define UMBRELLA_SET_H

#include <CGAL/Delaunay_triangulation_3.h>

#include <list>		  
#include <vector>	  //container of all points
#include <utility>
#include <algorithm>  //swap
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <bitset>
#include <deque>

//=========================extended================

#include <CGAL/IO/Umbrella_to_medit.h>
#include <CGAL/Constrain_surface_3.h>
#include <CGAL/Stretched_delaunay_3.h>
#include <CGAL/Conflict_zone.h>
#include <CGAL/Conflict_surface_triangles.h>
#include <CGAL/Conflict_volume_tetrahedra.h>
#include <CGAL/Surface_criteria.h>
#include <CGAL/Volume_criteria.h>

template <typename K, typename Metric_field_3, typename Constrain_surface>
class Umbrella_set{

	public:
		typedef typename K::Point_3					  Point_3;
		typedef typename K::FT						  FT;

		typedef Stretched_traits_3<K,Metric_field_3>				Traits;
		typedef Stretched_delaunay_3<K,Traits,Constrain_surface>	Umbrella;
		typedef Surface_criteria_for_umbrella_set<Umbrella_set,Constrain_surface>	Surface_criteria;
		typedef Volume_criteria_for_umbrella_set<Umbrella_set>						Volume_criteria;
		typedef typename CGAL::Delaunay_triangulation_3<K,
				CGAL::Triangulation_data_structure_3<
					CGAL::Triangulation_vertex_base_with_info_3<unsigned, K> > >
					Global_triangulation;

		typedef typename Global_triangulation::Vertex_handle  Global_vertex_handle;
		typedef typename Global_triangulation::Cell_handle	  Global_cell_handle;
		typedef typename Umbrella::Vertex_handle              Vertex_handle;
		typedef typename Umbrella::Cell                       Cell;
		typedef typename Umbrella::Cell_handle                Cell_handle;
		typedef typename Umbrella::Finite_cells_iterator      Finite_cells_iterator;
		typedef typename Umbrella::Finite_facets_iterator     Finite_facets_iterator;
		typedef typename Umbrella::Neighbor_cells_iterator	  Neighbor_cells_iterator;
		typedef typename Umbrella::Neighbor_bounds_iterator	  Neighbor_bounds_iterator;
		typedef typename Umbrella::Point                      Point;
		typedef typename Umbrella::Facet                      Facet;
		typedef typename Umbrella::Vertex_map                 Vertex_map;
		typedef typename Vertex_map::iterator				  Vertex_map_iterator;
		typedef typename std::list<Umbrella *>::iterator	  Umbrella_iterator;
		typedef Conflict_triangle_list<Umbrella_set,Constrain_surface>		Surface_conflict_list;
		typedef Conflict_tetrahedra_list<Umbrella_set>						Volume_conflict_list;
		typedef typename Surface_conflict_list::Conflict_triangle_handle	Conflict_triangle_handle;
		typedef typename Volume_conflict_list::Conflict_tetrahedra_handle	Conflict_tetrahedra_handle;

		//====================================================================================
		//	                    Constructions and initialization
		//====================================================================================
		/** Read from the constrain_surface the initial_points on the surface
		 *  to initialize the triangulation. Then initialize the border information
		 *  i.e. record for each cell if it is "inside" and if one of its facet forms 
		 *  a border
		 * 
		 *  a cell is inside if its circumcenter is "inside"
		 *
		 *  a cell has a facet on border if it is "inside" and its neighbor who share
		 *  that facet with it is not "inside"
		 */

		Umbrella_set(const Metric_field_3 &mf, const Constrain_surface &cs,
				double max_ratio,double max_distortion,double max_size,double max_distance,
				bool default_init_points = true)
			:surface_conflict_list(new Surface_conflict_list(*this)),
			volume_conflict_list(new Volume_conflict_list(*this)),
			surface_list_valid(false),
			umbrella_cache(NULL),
			mfield(new Metric_field_3(mf)),
			constrain_surface(new Constrain_surface(cs)){
				surface_criteria = new Surface_criteria(constrain_surface,
						max_ratio,max_distortion,max_size,max_distance);
				volume_criteria = new Volume_criteria(max_ratio,max_distortion,max_size);
				//max_ratio,max_distortion,max_volume

				// get the initial points
				if(default_init_points){
					std::vector<Point_3>	all_points = constrain_surface->initial_points();
					set_all_points(all_points);
				}
			}

		Umbrella_set()
			:
				constrain_surface(NULL),
				mfield(NULL),
				surface_list_valid(false),
				umbrella_cache(NULL),
				surface_criteria(NULL){
					surface_conflict_list = new Surface_conflict_list(*this);
					volume_conflict_list = new Volume_conflict_list(*this);
				}

		~Umbrella_set(){
			for(unsigned i=0;i<all_umbrella.size();i++){
				delete all_umbrella[i];
			}
			if(surface_conflict_list != NULL){
				delete surface_conflict_list;
			}
			if(volume_conflict_list != NULL){
				delete volume_conflict_list;
			}
			if(surface_criteria != NULL){
				delete surface_criteria;
			}
			if(mfield != NULL){
				delete mfield;
			}
			if(constrain_surface != NULL){
				delete constrain_surface;
			}
			if(volume_criteria != NULL){
				delete volume_criteria;
			}
			if(umbrella_cache != NULL){
				delete umbrella_cache;
			}
		}

		//====================================================================================
		//                   Initialization of parameters of umbrella_set
		//====================================================================================

		bool set_constrain_surface(const Constrain_surface &cs){
			if(constrain_surface != NULL){
				return false;
			}
			constrain_surface = new Constrain_surface(cs);
			return true;
		}

		bool set_metric_field(const Metric_field_3 &mf){
			if(mfield != NULL){
				return false;
			}
			mfield = new Metric_field_3(mf);
			return true;
		}

		bool set_surface_criteria(double max_ratio,
				double max_distortion,double max_size,double max_distance){
			if((constrain_surface == NULL) || (surface_criteria != NULL)){
				return false;
			}
			surface_criteria = new Surface_criteria(constrain_surface,max_ratio,max_distortion,max_size,max_distance);
			volume_criteria = new Volume_criteria(max_ratio,max_distortion,max_size,max_distance);
			return true;
		}

		//is_on_constraint contains information about if the point is on constraint or not
		void set_all_points(const std::vector<Point_3> &pc,bool empty_umb = false,
				const std::vector<bool> *all_on_constraint_input = NULL){
			if(all_on_constraint_input == NULL){
				for(unsigned i=0;i<pc.size();i++){
					add_point(pc[i],true);
				}
			}else{
				for(unsigned i=0;i<pc.size();i++){
					add_point(pc[i],(*all_on_constraint_input)[i]);
				}
			}
			init_all_umbrellas(empty_umb);
		}

		//====================================================================================
		//                   refinment of surface: Initialization
		//====================================================================================
	private:
		void init_all_umbrellas(bool empty_umb){
			/// initialize all umbrellas
			std::cout<<"\nStart initializing umbrellas";
			for(unsigned umb_index = 0;umb_index < all_vertices.size();umb_index++){
				Traits	trait(umb_index,all_vertices[umb_index]->point(),*mfield);
				Umbrella *new_umb = new Umbrella(*constrain_surface,trait,
						all_vertices[umb_index]->point(),umb_index);
				if(!empty_umb){
					init_umbrella(new_umb,std::deque<unsigned>(1,umb_index));
					if((umb_index%100)==99){
						std::cout<<"\nStart initializing "<<umb_index + 1<<"th umbrella";
					}
					std::cout<<".";
					std::cout.flush();
				}
				all_umbrella.push_back(new_umb);
			}
			surface_list_valid = false;
#ifdef DEBUG_PRINT_ALL_UMBRELLAS
			if(!empty_umb){
				print_all_umbrellas();
			}
#endif
		}

		//init_umbrella with ajacent vertex not assigned
		//	in case of a point not in global_triangulation, indices to start with will be the 4 vertices surrounding
		//	umb->center_point() in the global_triangulation, otherwise indices to start with should be the vertex its self.
		void init_umbrella(Umbrella *umb, const std::deque<unsigned> &indices_to_start_with){

			std::deque<unsigned>					  to_be_tried_queue(indices_to_start_with);
			std::set<unsigned>						  already_tried_set(indices_to_start_with.begin(),
					indices_to_start_with.end());

			std::vector<Global_vertex_handle>		  adjacent_vertex_cache;
			std::back_insert_iterator<std::vector<Global_vertex_handle> >	
				adjacent_vertex_cache_it(adjacent_vertex_cache);
			typename K::Compute_squared_distance_3	csd = 
				euclidean_trait.compute_squared_distance_3_object();
			Point_3	  center_point = umb->get_center()->point();
			FT conflict_range = umb->get_absolute_squared_conflict_radius();

			while(!to_be_tried_queue.empty()){
				unsigned current_index = to_be_tried_queue.front();
				Point_3 current_point=all_vertices[current_index]->point();
				//try insert the new point
				Vertex_handle new_vertex = umb->insert_if_conflict(
						all_vertices[current_index]->point(),current_index);
				//consider the points within conflict range
				if((csd(current_point,center_point) < conflict_range) ||
						(new_vertex != Vertex_handle())){
					//if insert successful, than maybe change the conflict radius.
					if(new_vertex!=Vertex_handle()){
						conflict_range = umb->get_absolute_squared_conflict_radius();
					}
					//append all vertices into consideration
					adjacent_vertex_cache.clear();
					global_triangulation.incident_vertices(
							all_vertices[current_index],adjacent_vertex_cache_it);
					for(unsigned i=0;i<adjacent_vertex_cache.size();i++){
						if(!global_triangulation.is_infinite(adjacent_vertex_cache[i])){
							unsigned new_index = adjacent_vertex_cache[i]->info();
							//if the point has not been considered yet
							if(already_tried_set.insert(new_index).second){
								to_be_tried_queue.push_back(new_index);
							}
						}
					}//end of for adjacent_vertex_cache
				}//end of if
				to_be_tried_queue.pop_front();
			}//while
		}//init_umbrella

	public:
		//init_umbrella with adjacent vertex already assigned.
		void init_umbrella(unsigned umb_index,std::vector<unsigned> &points){
			for(unsigned i=0;i<points.size();i++){
				unsigned point_index = points[i];
				assert(point_index <= all_vertices.size());
				Vertex_handle v = 
					all_umbrella[umb_index]->insert(all_vertices[point_index]->point(),point_index);
			}
#ifdef DEBUG_PRINT_ALL_UMBRELLAS
			print_umbrella_to_medit(umb_index);
			print_umbrella_surface_to_medit(umb_index);
#endif //DEBUG_PRINT_ALL_UMBRELLAS
		}//init_umbrella

		//====================================================================================
		//                   refinement of umbrella: add new point to global_triangulation
		//====================================================================================

		void add_point(const Point &p,bool on_constraint){
			unsigned index = size();
			Global_vertex_handle v = global_triangulation.insert(p);
			if(size() == index){
				std::cout<<"Warning: Inserting a point already exist:"
					<<index<<" = "<<v->info()<<std::endl;
			}else{
				v->info() = index;
				all_vertices.push_back(v);
				all_on_constraint.push_back(on_constraint);
			}
		}

		//====================================================================================
		//                   refinement of umbrella: initiate new umbrellas
		//====================================================================================


		///Indices is the container of indices(unsigned)
		/// Initiate the new umbrella with the last pont in the all_vertices
		/// initiate this umbrella with the points of the indices and then
		/// the neighbors of its indices...(BFS)
		template<typename Indices>
		void add_umbrella(Indices &indices){
			unsigned last_index = all_vertices.size()-1;
			Traits	trait(last_index,all_vertices[last_index]->point(),*mfield);
			Umbrella *new_umb = new Umbrella(*constrain_surface,trait,
					all_vertices[last_index]->point(),last_index);
			std::deque<unsigned>  indices_to_start_with;
			for(unsigned i=0;i<indices.size();i++){
				indices_to_start_with.push_front(indices.umb_index(i));
			}
			init_umbrella(new_umb,indices_to_start_with);
			all_umbrella.push_back(new_umb);

			//surface need updating next time update_surface is called
			if(new_umb->on_boundary()){
				surface_list_valid = false;
			}

#ifdef DEBUG_PRINT_ALL_UMBRELLAS
			print_umbrella_to_medit(last_index);
			print_umbrella_surface_to_medit(last_index);
#endif
		}

		void add_umbrella_cache(){
			all_umbrella.push_back(umbrella_cache);
			if(umbrella_cache->on_boundary()){
				surface_list_valid = false;
			}
			umbrella_cache = NULL;
#ifdef DEBUG_PRINT_ALL_UMBRELLAS
			unsigned last_index = all_umbrella.size()-1;
			print_umbrella_to_medit(last_index);
			print_umbrella_surface_to_medit(last_index);
#endif
		}
		//====================================================================================
		//                   refinment of surface: Calculation of Surface Conflicts
		//====================================================================================
		bool test_surface_coherency(unsigned center,unsigned i,unsigned j){
			Facet f_dummy;
			return	  all_umbrella[j]->has_surface(i,center,f_dummy)
				&& all_umbrella[i]->has_surface(j,center,f_dummy);
		}

		///\return Cell_handle with the given indices in the umbrella centered at center;
		/// if this cell_handle no longer exists than return Cell_handle()
		Cell_handle get_cell_handle(unsigned center, unsigned i,unsigned j,unsigned k)const{
			Cell_handle to_be_returned;
			all_umbrella[center]->has_tetrahedra(i,j,k,to_be_returned);
			return to_be_returned;
		}

		bool test_volume_coherency(unsigned center,unsigned i,unsigned j,unsigned k){
			Cell_handle c_dummy;
			return	all_umbrella[i]->has_tetrahedra(center,j,k,c_dummy)
				&& all_umbrella[j]->has_tetrahedra(center,i,k,c_dummy)
				&& all_umbrella[k]->has_tetrahedra(center,i,j,c_dummy);
		}

		template<typename Indices>
		void append_new_surface_conflicts(const Indices &indices){
			for(unsigned i=0;i<indices.size();i++){
				indices_cache.push_back(indices.umb_index(i));
				append_new_surface_conflicts(indices.umb_index(i));
			}//for (all the umbrellas)
			append_new_surface_conflicts(size()-1);
		}

		void append_new_volume_conflicts_from_indices_cache(){
			for(unsigned i=0;i<indices_cache.size();i++){
				append_new_volume_conflicts(indices_cache[i]);
			}
			append_new_volume_conflicts(size()-1);
			indices_cache.clear();
		}		

		template<typename Indices>
		void append_new_volume_conflicts(const Indices &indices){
			for(unsigned i=0;i<indices.size();i++){
				append_new_volume_conflicts(indices.umb_index(i));
			}
			append_new_volume_conflicts(size()-1);
		}

		void append_new_surface_conflicts(unsigned umb_index){
			unsigned index[2];				//The other two index of the triangle, knowing the center
			all_umbrella[umb_index]->update_neighbor_bounds();
			for(typename std::list<Facet>::iterator fit = 
					all_umbrella[umb_index]->neighbor_bounds.begin();
					fit!=all_umbrella[umb_index]->neighbor_bounds.end();
					fit++){					  //for all surface triangles:
				int k=0;					  //find out the two indices that are not center, store in index[2]
				for(int j=0;j<4;j++){
					index[k]=fit->first->vertex(j)->info();
					if((fit->second != j)
							&& index[k]!=umb_index){
						k++;
					}
				}
				assert(k==2);
				//test coherency
				bool coherent = test_surface_coherency(umb_index,index[0],index[1]);
				surface_conflict_list->push_if_bad(
						all_umbrella[umb_index]->get_traits_ref(),index[0],index[1],
						umb_index,	coherent);
			}//for all the facets
		}


		void append_new_volume_conflicts(unsigned umb_index){
			unsigned index[3];
			for(Neighbor_cells_iterator cit = all_umbrella[umb_index]->neighbor_cells_begin();
					cit!=all_umbrella[umb_index]->neighbor_cells_end();
					cit++){
				int k=0;
				if(all_umbrella[umb_index]->cell_is_inside(*cit)){
#ifdef DEBUG_VOLUME_REFINEMENT_POINT
					bool neighbor_cells_around_center = false;
#endif //DEBUG_VOLUME_REFINEMENT_POINT
					for(int j=0;j<4;j++){
						unsigned cur_index = (*cit)->vertex(j)->info();
						if(cur_index!=umb_index){
							index[k] = cur_index;
							k++;
						}
#ifdef DEBUG_VOLUME_REFINEMENT_POINT
						else{
							neighbor_cells_around_center = true;
						}
#endif //DEBUG_VOLUME_REFINEMENT_POINT
					}
#ifdef DEBUG_VOLUME_REFINEMENT_POINT
					assert(k==3);
					assert(neighbor_cells_around_center);
#endif //DEBUG_VOLUME_REFINEMENT_POINT
					bool coherent = test_volume_coherency(umb_index,index[0],index[1],index[2]);
					volume_conflict_list->push_if_bad(index[0],index[1],
							index[2],umb_index,coherent);
				}//if cell is inside the domain
			}
		}

		//====================================================================================
		//                   refinment of surface: Get the "worst" confict
		//====================================================================================

		Conflict_triangle_handle  next_surface_conflict(){
			return surface_conflict_list.top();
		}

		Conflict_tetrahedra_handle next_volume_conflict(){
			return volume_conflict_list.top();
		}

		//====================================================================================
		//                   refinment of surface: insert a point into the given umbrella
		//====================================================================================

		//add p into the umbrella indexed with i, with the initial posotion of cell_loc
		void refine_umbrella_surface(const Point_3 &p, Cell_handle cell_loc,unsigned i){
#ifdef DEBUG_CONFLICTS
			std::cout<<"Insert into umb"<<i<<std::endl;
#endif //DEBUG_CONFLICTS
			Vertex_handle v=all_umbrella[i]->insert(p,all_umbrella.size(),cell_loc);
		}

		void refine_umbrella_volume(const Point_3 &p,Cell_handle cell_loc,unsigned i){
			refine_umbrella_surface(p,cell_loc,i);
		}

		//====================================================================================
		//                   refinment of surface: Generate a point from a conflict
		//====================================================================================

		Point_3	refine_surface_point(const Conflict_triangle_handle & tr){
			unsigned center_index = tr.get_center();
			Point_3 new_point = all_umbrella[center_index]->refine_surface_point(
					tr.index(0),tr.index(1));
#ifdef DEBUG_VOLUME_REFINEMENT_POINT
			std::cout<<"SURFACE center";
			print_point(new_point);
			std::cout<<"center of "<<center_index<<" "<<tr.index(0)<<" "<<tr.index(1)<<std::endl;
			print_point(all_vertices[tr.index(0)]->point());
			print_point(all_vertices[tr.index(1)]->point());
			print_point(all_vertices[center_index]->point());
#endif
			return new_point;
		}

		Point_3 refine_volume_point(const Conflict_tetrahedra_handle & tet){
#ifdef DEBUG_VOLUME_REFINEMENT_POINT
			std::cout<<"\nCenter:"<<tet.get_center();
			std::cout<<"\nCell "
				<<tet.index(0)<<"  "
				<<tet.index(1)<<"  "
				<<tet.index(2)<<"  "
				<<std::endl;
#endif
			unsigned umb_index = tet.get_center();
			Cell_handle cell;
			all_umbrella[umb_index]->has_tetrahedra(
					tet.index(0),tet.index(1),tet.index(2),cell);
			return all_umbrella[umb_index]->circumcenter(cell);
		}
		//====================================================================================
		//                   update of surface: judge if p is in conflict with umbrellas
		//====================================================================================

		bool point_out_of_convex_hull(const Point_3 &p,Global_cell_handle &gch){
			gch = global_triangulation.locate(p);
			if(global_triangulation.is_infinite(gch)){
			   return true;
			}
			return false;
		}

		template<typename Conflict_zone>
		void add_all_conflicting_surface(const Point_3 &p, const Global_cell_handle &gch,
				Conflict_zone &conflict_zone){
			std::vector<Global_vertex_handle>	conflict_vertices;
			std::back_insert_iterator<std::vector<Global_vertex_handle> >
												conflict_vertices_iterator(conflict_vertices);
			global_triangulation.vertices_in_conflict(p,gch,conflict_vertices_iterator);
			for(unsigned i=0;i<conflict_vertices.size();i++){
				if(!global_triangulation.is_infinite(conflict_vertices[i])){
					conflict_zone.push(conflict_vertices[i]->info(),Cell_handle());
				}
			}
		}

		/// \returns If the point p is in conflict with the umbrella[index] or not
		bool point_in_conflict(const Point_3 &p, Cell_handle &cell_loc,unsigned index){
			return all_umbrella[index]->test_conflict(p,cell_loc);
		}

		bool check_encroachment(const Point_3 &p,unsigned index){
			Facet facet;		//To store the facet being encroached
			if(all_umbrella[index]->check_encroachment(p,facet)==true){
				unsigned ind[2];unsigned k=0;
				for(int i=0;i<4;i++){
					if(i!=facet.second && (facet.first->vertex(i)->info()!=index) ){
						ind[k] = facet.first->vertex(i)->info();
						k++;
					}
				}
				surface_conflict_list->push_if_bad(all_umbrella[index]->get_traits_ref(),
						ind[0],ind[1],index,false);
#ifdef DEBUG_VOLUME_REFINEMENT_POINT
				std::cout<<"Encroaches surface: "<<ind[0]<<" "<<ind[1]<<" "<<index<<"\n";
#endif //DEBUG_VOLUME_REFINEMENT_POINT
				return true;
			}
			return false;
		}

		/// Initialize the trial umbrella -- "umbrella_cache", centered with the given point,
		/// and supported by the other existing points. 	Test if the center of this umbrella
		/// encroache the surface in its own metric.
		/// If no encroachment happens, the umbrella_cache will be appended to all_umbrella directly
		/// in the next step.
		/// Otherwise the umbrella_cache will be abandoned and the encroached surface is required 
		/// to be refined.
		template<typename Indices>
		bool new_center_test_encroachment(const Point_3 &p,const Indices &indices){
			assert(umbrella_cache==NULL);
			unsigned umb_index = size();
			Traits	trait(umb_index,p,*mfield);
			umbrella_cache = new Umbrella(*constrain_surface,trait,p,umb_index);

			indices.umb_index(0);

			//TODO improve location complexity with indices info
			Global_cell_handle cell_loc = global_triangulation.locate(p);
			std::deque<unsigned>  indices_to_start_with;
			for(int i=0;i<4;i++){
				indices_to_start_with.push_front(cell_loc->vertex(i)->info());
			}
			init_umbrella(umbrella_cache,indices_to_start_with);

#ifdef DEBUG_PRINT_ALL_UMBRELLAS
			std::ofstream of("data/umbrella_cache.mesh");
			umbrella_to_medit<Umbrella>(umbrella_cache,of);
			of.close();
#endif
			unsigned index[3];
			if(umbrella_cache->center_encroach_surface(this,index)){
				surface_conflict_list->push_if_bad(all_umbrella[index[2]]->get_traits_ref(),
						index[0],index[1],index[2],false);
#ifdef DEBUG_VOLUME_REFINEMENT_POINT
				std::cout<<"Encrochment to existing surfaces"<<index[0]<<" "
					<<index[1]<<" "<<index[2]<<" "<<std::endl;
#endif //DEBUG_VOLUME_REFINEMENT_POINT
				return true;
			}
#ifdef DEBUG_VOLUME_REFINEMENT_POINT
			std::cout<<"No encrochment to existing surfaces...\n";
#endif //DEBUG_VOLUME_REFINEMENT_POINT
			return false;
		}

		void invalidate_umbrella_cache(){
			delete umbrella_cache;
			umbrella_cache = NULL;
		}
		//====================================================================================
		//                   update of surface: pop out inexistant conflicts
		//====================================================================================

		///After each refinement, pop out the conflicts until empty or until the next top
		/// conflict is existant.
		void pop_invalid_surface_conflicts(){
			//a conflict is fake if it does no longer existe in umbrella
			Conflict_triangle_handle  conflict_triangle = surface_conflict_list->top();
#ifdef DEBUG_CONFLICT_LISTS
			print_top_surface_conflict();
#endif
			bool fake_conflict = ! is_surface_conflict_exist(conflict_triangle);
			while(fake_conflict && (!surface_conflict_list->empty())){
#ifdef DEBUG_CONFLICT_LISTS
				std::cout<<"pop out inexistant conflict"<<std::endl;
#endif 
				surface_conflict_list->pop();
				conflict_triangle = surface_conflict_list->top();
				fake_conflict = ! is_surface_conflict_exist(conflict_triangle);
#ifdef DEBUG_CONFLICT_LISTS
				print_top_surface_conflict();
#endif 
			}
		}

		void pop_invalid_volume_conflicts(){
			Conflict_tetrahedra_handle conflict_tetrahedron = volume_conflict_list->top();
#ifdef DEBUG_CONFLICT_VOLUME_LISTS
			print_top_volume_conflict();
#endif //DEBUG_CONFLICT_LISTS
			bool fake_conflict = ! is_volume_conflict_exist(conflict_tetrahedron);
			while(fake_conflict && (!volume_conflict_list->empty())){
#ifdef DEBUG_CONFLICT_VOLUME_LISTS
				std::cout<<"pop out inexistant conflict"<<std::endl;
#endif 
				volume_conflict_list->pop();
				conflict_tetrahedron = volume_conflict_list->top();
				fake_conflict = ! is_volume_conflict_exist(conflict_tetrahedron);
#ifdef DEBUG_CONFLICT_VOLUME_LISTS
			print_top_volume_conflict();
#endif //DEBUG_CONFLICT_LISTS
			}
		}//pop_invalid_volume_conflicts


#ifdef DEBUG_CONFLICT_VOLUME_LISTS

		void print_top_volume_conflict(){	
			Conflict_tetrahedra_handle cth = volume_conflict_list->top();
			std::cout<<"Current top volume conflict:("<<
				cth.index(0)<<","<<
				cth.index(1)<<","<<
				cth.index(2)<<","<<
				cth.get_center()<<") with center"<<
				cth.get_center()<<"\nTotal conflicts:"<<
				volume_conflict_list->size()<<std::endl;
		}
#endif

#ifdef DEBUG_CONFLICT_LISTS

		void print_top_surface_conflict(){	
			Conflict_triangle_handle  conflict_triangle = surface_conflict_list->top();
			std::cout<<"Current top suface conflict:("<<conflict_triangle.index(0)<<","
				<<conflict_triangle.index(1)<<","<<conflict_triangle.get_center()<<") with center"
				<<conflict_triangle.get_center()<<"\nTotal conflicts:"<<
				surface_conflict_list->size()<<std::endl;
		}
#endif

		bool is_surface_conflict_exist(const Conflict_triangle_handle & cth){
			unsigned center = cth.get_center();
			Facet f_dummy;
			return all_umbrella[center]->has_surface(cth.index(0),cth.index(1),f_dummy);
		}

		bool is_volume_conflict_exist(const Conflict_tetrahedra_handle &cth){
			unsigned center = cth.get_center();
			Cell_handle c_dummy;
			return all_umbrella[center]->has_tetrahedra(cth.index(0),cth.index(1),cth.index(2),c_dummy);
		}

		//====================================================================================
		//       Initialization for the surface_conflict_list
		//====================================================================================
		void init_all_surface_conflict_list(){
			std::cout<<"\nInitialize surface conflict..."<<std::endl;
			for(unsigned i=0;i<all_umbrella.size();i++){
				append_new_surface_conflicts(i);
			}
		}


		//====================================================================================
		//       Initialization for the volume_conflict_list
		//====================================================================================
		void init_all_volume_conflict_list(){
			std::cout<<"\nInitialize volume conflict..."<<std::endl;
			for(unsigned i=0;i<all_umbrella.size();i++){
				append_new_volume_conflicts(i);
			}
		}//init_all_volume_conflict_list


		//====================================================================================
		//       Update for output: get all surface triangles, conflicted and coherent
		//====================================================================================
		bool is_surface(unsigned i,unsigned j,unsigned k){
			Facet f_dummy;
			return all_umbrella[i]->has_surface(j,k,f_dummy);
		}

		void update_all_surfaces(){
			// new insertion into the surface will make surface invalid
			if(surface_list_valid){
				return;
			}
			surface_list_valid = true;
			surface_list.clear();
			unsigned index[4];
			for(unsigned i=0;i<all_umbrella.size();i++){
				all_umbrella[i]->update_neighbor_bounds();
				for(typename std::list<Facet>::iterator fit = 
						all_umbrella[i]->neighbor_bounds.begin();
						fit!=all_umbrella[i]->neighbor_bounds.end();
						fit++){
					int k=0;
					for(int j=0;j<4;j++){
						unsigned cur_index = fit->first->vertex(j)->info();
						if((fit->second != j)
								&& cur_index!=i){
							index[k]=cur_index;
							k++;
						}
					}
					if(k>2){
						std::cout<<"i = "<<i<<" fit->first->vertex(fit->second)->info() = "<<
							fit->first->vertex(fit->second)->info()<<" k="<<k<<" index =";
						for(int xx=0;xx<4;xx++){
							std::cout<<index[xx]<<"  ";
						}
						std::cout<<std::endl;
						for(int xx=0;xx<4;xx++){
							std::cout<<fit->first->vertex(xx)->info()<<"  ";
						}
						std::cout<<std::endl;
					}
					assert((index[1]!=i) && (index[0]!=i));
					assert(k==2);
					bool coherent = test_surface_coherency(i,index[0],index[1]);
					surface_list.push_back(index[0],index[1],i,coherent);
				}
			}
			surface_list.sort();
		}

		void update_all_tetrahedra(){
			tetrahedra_list.clear();
			unsigned index[3];	
			for(unsigned i=0;i<all_umbrella.size();i++){
				for(Neighbor_cells_iterator cit = all_umbrella[i]->neighbor_cells_begin();
						cit != all_umbrella[i]->neighbor_cells_end();
						cit ++){
					if(all_umbrella[i]->cell_is_inside(*cit)){
						unsigned k=0;
						for(int j=0;j<4;j++){
							if((*cit)->vertex(j)->info()!=i){
								index[k]=(*cit)->vertex(j)->info();
								k++;
							}
						}
						if(test_volume_coherency(i,index[0],index[1],index[2])){
							tetrahedra_list.push_back(index[0],index[1],index[2],i);
						}
					}//end of if cell is inside of domain
				}//end of for all cells
			}//end of all vertices
			tetrahedra_list.make_unique();
		}
		//====================================================================================
		//                 Print to files
		//====================================================================================



		void print_all_umbrellas(){
			for(unsigned i=0; i< all_umbrella.size();i++){
				print_umbrella_to_medit(i);
				print_umbrella_surface_to_medit(i);
			}
		}

		void print_umbrella_to_medit(unsigned index){
			std::stringstream ss;
			ss<<"data/umb"<<index<<".mesh";
			std::ofstream of(ss.str().c_str());
			umbrella_to_medit<Umbrella>(all_umbrella[index],of);
			of.close();
		}

		void print_umbrella_surface_to_medit(unsigned index){
			std::stringstream ss;
			ss<<"data/suf"<<index<<".mesh";
			std::ofstream of(ss.str().c_str());
			umbrella_neighbor_surface_to_medit<Umbrella>(all_umbrella[index],of);
			of.close();
		}



		//====================================================================================
		//                 attributes
		//====================================================================================

		unsigned size() const {
			return global_triangulation.number_of_vertices();
		}

		Point_3 get_point(unsigned i)const{
			return all_vertices[i]->point();
		}

		bool is_on_boundary(unsigned i)const{
			return all_umbrella[i]->on_boundary();
		}

		bool is_on_constraint(unsigned i)const{
			return all_on_constraint[i];
		}

		const Traits & get_traits(unsigned i)const{
			return all_umbrella[i]->get_traits();
		}

		const Surface_criteria * get_surface_criteria() const{
			return surface_criteria;
		}

		const Volume_criteria * get_volume_criteria() const{
			return volume_criteria;
		}

		//====================================================================================
		//                 members
		//====================================================================================


		K											  euclidean_trait;
		Surface_list								  surface_list;
		Tetrahedra_list								  tetrahedra_list;
		std::vector<Umbrella *>						  all_umbrella;
		Surface_conflict_list						 *surface_conflict_list;
		Volume_conflict_list						 *volume_conflict_list;
	private:

		bool											surface_list_valid;
		Umbrella										*umbrella_cache;
		std::vector<unsigned>							indices_cache;
		Global_triangulation							global_triangulation;
		std::vector<Global_vertex_handle>				all_vertices;
		std::vector<bool>								all_on_constraint;
		Metric_field_3									*mfield;
		Constrain_surface								*constrain_surface;
		Surface_criteria								*surface_criteria;
		Volume_criteria									*volume_criteria;	
};


#endif //UMBRELLA_SET_H

#ifndef REFINE_LEVELS_H
#define REFINE_LEVELS_H

#include <CGAL/Mesher_level.h>
#include <iostream>

//==================extended=================
#include <CGAL/Conflict_zone.h>

template<typename Umbrella_set>
class Surface_mesher_level:
	public CGAL::Mesher_level<
	Umbrella_set,											//triangulation
	Surface_mesher_level<Umbrella_set>,						//Self
	typename Umbrella_set::Conflict_triangle_handle,		//Element
	CGAL::Null_mesher_level,								//Previous
	Umbrella_set_mesher_level_traits_3<Umbrella_set> >		//traits
{
	public:

		typedef CGAL::Mesher_level<Umbrella_set,						//triangulation
				Surface_mesher_level<Umbrella_set>,						//Derived
				typename Umbrella_set::Conflict_triangle_handle,		//Element
				CGAL::Null_mesher_level,								//Previous
				Umbrella_set_mesher_level_traits_3<Umbrella_set> >			  Mesher_level;
		typedef typename Umbrella_set::Conflict_triangle_handle		Element;
		typedef typename Umbrella_set::Point						Point;
		typedef typename Umbrella_set::Cell_handle					Cell_handle;
		typedef typename Umbrella_set::Global_cell_handle			Global_cell_handle;
		typedef typename Mesher_level::Vertex_handle				Vertex_handle;
		typedef typename Mesher_level::Zone							Zone;
		typedef typename CGAL::Mesher_level_conflict_status Mesher_level_conflict_status;
		typedef typename CGAL::Null_mesher_level					Null_mesher_level;

		Surface_mesher_level(Umbrella_set& _u,Null_mesher_level &null_mesher_lever):
			Mesher_level(null_mesher_lever),
			umbrella_set(_u){}
		///< Unit element to be refined.

		/** \ name implimentation functions:*/

		/** insert implimentations*/
		/// 1.Insert the point to all umbrellas in the zone. 
		///
		///	2. Update marks the triangulation in every umbrella.
		///
		///	3. Construct the new umbrella, adding new vertex using (TODO BFS searching
		///    starting with Zone.)
		///
		/// 4. Test all new triangles, append conflicts to the list.
		Vertex_handle insert_impl(const Point &p,Zone &zone){
#ifdef DEBUG_REFINE_PROCESS
			std::cout<<"Enter insertion period."<<std::endl;
#endif
			umbrella_set.add_point(p,true);
			for(unsigned i=0;i<zone.size();i++){
				umbrella_set.refine_umbrella_surface(p,zone.umb_loc(i),zone.umb_index(i));	// Complete 1. and 2.
			}
			umbrella_set.add_umbrella(zone);				// Do 3. and append the last index to zone
			umbrella_set.append_new_surface_conflicts(zone);  //Do 4.


#ifdef DEBUG_REFINE_PROCESS
			std::cout<<"Finish insertion period."<<std::endl;
#endif
			return Vertex_handle();							//a dummy return value
		}

		/**Find the conflict zone of point and element*/
		///Calculate indices of conflicted umbrellas
		Zone conflicts_zone_impl(const Point &p,Element &){
			Zone conflict_zone;
			Global_cell_handle	gch;
			if(umbrella_set.point_out_of_convex_hull(p,gch)){
				umbrella_set.add_all_conflicting_surface(p,gch,conflict_zone);
			}
			for(unsigned i=0;i<umbrella_set.size();i++){
				Cell_handle cell_loc;
				if((!conflict_zone.exist(i)) &&
						umbrella_set.point_in_conflict(p,cell_loc,i)){
					conflict_zone.push(i,cell_loc);
				}
			}
			return conflict_zone;
		}


		Umbrella_set &triangulation_ref_impl(){
			return umbrella_set;
		}

		const Umbrella_set &triangulation_ref_impl() const {
			return umbrella_set;
		}

		/** Called before the first refinement, to initialized the queue of
		 *       elements that should be refined. */
		void scan_triangulation_impl(){
			umbrella_set.init_all_surface_conflict_list();
		}

		/** \return if there are elements to return in this 
		 * level */
		bool no_longer_element_to_refine_impl(){
			return umbrella_set.surface_conflict_list->empty();
		}

		/** \return Get the next element to refine
		 * */
		Element get_next_element_impl(){
#ifdef DEBUG_REFINE_PROCESS
			umbrella_set.print_top_surface_conflict();
#endif
			return umbrella_set.surface_conflict_list->top();
		}

		/** delete the first element in the list
		 * */
		void pop_next_element_impl(){
			umbrella_set.surface_conflict_list->pop();
		}

		/** \return Calculate the refinement point for the element
		 * */
		Point refinement_point_impl(const Element & tr){
			return umbrella_set.refine_surface_point(tr);
		}

		///No action
		void before_conflicts_impl(const Element&,const Point&){}

		/** Tells if, as regards this level of the refinement process, if the
		 *  point conflicts with something, and do what is needed. The return
		 *  type is made of two booleans:
		 *   - the first one tells if the point can be inserted,
		 *   - in case of, the first one is \c false, the second one tells if
		 *  the tested element should be reconsidered latter.
		 *  */
		/// Always no conflict for surface mesher
		Mesher_level_conflict_status 
			private_test_point_conflict_impl(const Point&,const Zone&){
				return CGAL::NO_CONFLICT;
			}

		///Given a point and its conflict zone, check if in all conflict zone, 
		/// the point does not encroach any surface of constrain.
		///If encroachment happens, the encroached surface will be pushed into
		///conflict list.
		Mesher_level_conflict_status
			test_point_conflict_from_superior_impl(const Point& p, Zone& zone){
				for(unsigned i = 0;i<zone.size();i++){
					if(umbrella_set.check_encroachment(p,zone.umb_index(i))){
						return CGAL::CONFLICT_BUT_ELEMENT_CAN_BE_RECONSIDERED;
					}
				}
				return CGAL::NO_CONFLICT;
			}

		///No action needed
		void before_insertion_impl(const Element&, const Point&, Zone&){}

		///Adjust the surface conflict for the next pop
		void after_insertion_impl(const Vertex_handle &v){
			umbrella_set.pop_invalid_surface_conflicts();
		}

		///No action needed
		void after_no_insertion_impl(const Element&,const Point&,Zone&){}

	private:

		Umbrella_set	&umbrella_set;

};

template<typename Umbrella_set,typename Surface_mesher_level>
class Volume_mesher_level:
	public CGAL::Mesher_level<
	Umbrella_set,												//triangulation
	Volume_mesher_level<Umbrella_set,Surface_mesher_level>, 
	typename Umbrella_set::Conflict_tetrahedra_handle,			//Element
	Surface_mesher_level,										//Previous
	Umbrella_set_mesher_level_traits_3<Umbrella_set> >		//traits
{
	public:

		typedef CGAL::Mesher_level<
			Umbrella_set,												//triangulation
			Volume_mesher_level<Umbrella_set,Surface_mesher_level>, 
			typename Umbrella_set::Conflict_tetrahedra_handle,			//Element
			Surface_mesher_level,										//Previous
			Umbrella_set_mesher_level_traits_3<Umbrella_set> >		Mesher_level;

		typedef typename Umbrella_set::Conflict_tetrahedra_handle	Element;
		typedef typename Umbrella_set::Point						Point;
		typedef typename Umbrella_set::Cell_handle					Cell_handle;
		typedef typename Mesher_level::Vertex_handle				Vertex_handle;
		typedef typename Mesher_level::Zone							Zone;
		typedef typename CGAL::Mesher_level_conflict_status		Mesher_level_conflict_status;
		typedef typename CGAL::Null_mesher_level					Null_mesher_level;

		Volume_mesher_level(Umbrella_set& _u,Surface_mesher_level &surface_mesher_lever):
			Mesher_level(surface_mesher_lever),
			umbrella_set(_u){}

		/** \ name implimentation functions:*/

		/** insert implimentations*/
		/// 1.Insert the point to all umbrellas in the zone. 
		///
		///	2. Update marks the triangulation in every umbrella.
		///
		///	3. Add the new umbrella in cache to the all_umbrella list
		///
		/// 4. Test all new triangles, append conflicts to the list.
		Vertex_handle insert_impl(const Point &p,Zone &zone){
#ifdef DEBUG_REFINE_PROCESS
			std::cout<<"Enter insertion period."<<std::endl;
#endif
			umbrella_set.add_point(p,false);
			for(unsigned i=0;i<zone.size();i++){
				umbrella_set.refine_umbrella_volume(p,zone.umb_loc(i),zone.umb_index(i));	// Complete 1. and 2.
			}
			umbrella_set.add_umbrella_cache();				// Do 3. and append the last index to zone
			umbrella_set.append_new_volume_conflicts(zone);  //Do 4.

#ifdef DEBUG_REFINE_PROCESS
			std::cout<<"Finish insertion period."<<std::endl;
#endif
			return Vertex_handle();							//a dummy return value
		}

		/**Find the conflict zone of point and element*/
		///Calculate indices of conflicted umbrellas
		Zone conflicts_zone_impl(const Point &p,Element &){
			Zone conflict_zone;
			for(unsigned i=0;i<umbrella_set.size();i++){
				Cell_handle cell_loc;
				//If conflict, cell_loc will note the cell in conflict
				if(umbrella_set.point_in_conflict(p,cell_loc,i)){
					conflict_zone.push(i,cell_loc);
				}	
			}
			return conflict_zone;
		}


		Umbrella_set &triangulation_ref_impl(){
			return umbrella_set;
		}

		const Umbrella_set &triangulation_ref_impl() const {
			return umbrella_set;
		}

		/** Called before the first refinement, to initialized the queue of
		 *       elements that should be refined. */
		void scan_triangulation_impl(){
			umbrella_set.init_all_volume_conflict_list();
		}

		/** \return if there are elements to return in this 
		 * level */
		bool no_longer_element_to_refine_impl(){
			return umbrella_set.volume_conflict_list->empty();
		}

		/** \return Get the next element to refine
		 * */
		Element get_next_element_impl(){
			return umbrella_set.volume_conflict_list->top();
		}

		/** delete the first element in the list
		 * */
		void pop_next_element_impl(){
			umbrella_set.volume_conflict_list->pop();
		}

		/** \return Calculate the refinement point for the element
		 * */
		Point refinement_point_impl(const Element & tr){
			Point p = umbrella_set.refine_volume_point(tr);
#ifdef DEBUG_VOLUME_REFINEMENT_POINT
			std::cout<<"Refinement point is:";
			print_point(p);
#endif
			return p;
		}

		void before_conflicts_impl(const Element&,const Point&){}

		/** Tells if, as regards this level of the refinement process, if the
		 *  point conflicts with something, and do what is needed. The return
		 *  type is made of two booleans:
		 *   - the first one tells if the point can be inserted,
		 *   - in case of, the first one is \c false, the second one tells if
		 *  the tested element should be reconsidered latter.
		 *  */
		//construct the umbrella cache and test if it encroaches any point
		Mesher_level_conflict_status 
			private_test_point_conflict_impl(const Point& p,const Zone& zone){
				if(umbrella_set.new_center_test_encroachment(p,zone)){
					return CGAL::CONFLICT_BUT_ELEMENT_CAN_BE_RECONSIDERED;
				}
				return CGAL::NO_CONFLICT;
			}

		/** Tells if, as regards this level of the refinement process, if the
		 *  point conflicts with something, and do what is needed. The return
		 *  type is made of two booleans:
		 *   - the first one tells if the point can be inserted,
		 *   - in case of, the first one is \c false, the second one tells if
		 *  the tested element should be reconsidered latter.
		 *  This function is called by the superior level, if any.
		 *   */
		Mesher_level_conflict_status
			test_point_conflict_from_superior_impl(const Point&, Zone&){
				return CGAL::NO_CONFLICT;
			}

		///No action needed
		void before_insertion_impl(const Element&, const Point&, Zone&){}

		///Adjust the surface conflict for the next pop
		void after_insertion_impl(const Vertex_handle &v){
			umbrella_set.pop_invalid_volume_conflicts();
		}

		///No action needed
		///Invalidate the cache of new_umb.
		void after_no_insertion_impl(const Element&,const Point&,Zone&){
			umbrella_set.invalidate_umbrella_cache();
		}

	private:
		Umbrella_set	&umbrella_set;
		Zone			zone_cache;

		///Important!!
		// If possible, add "zone" variable to the "after_insertion" function 
		//so that the function can act accordingly.

};


#endif //REFINE_LEVELS_H

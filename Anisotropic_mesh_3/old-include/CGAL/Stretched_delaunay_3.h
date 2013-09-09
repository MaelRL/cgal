#ifndef STRETCHED_DELAUNAY_3_H
#define STRETCHED_DELAUNAY_3_H

/** Define the stretched constrained delaunay triangulation
 * structure in 3 dimension
 * */

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_circumcenter_3.h>
#include <CGAL/Object.h>
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Color.h>


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iterator>			//insert iterator
#include <utility>
#include <map>				//map the index to vertex
#include <bitset>			//mark boundary surface in cell
#include <vector>			//neighbor_cells
#include <functional>		//less()
#include <algorithm>

//================extended=================
#include <CGAL/Stretched_traits_3.h>
#include <CGAL/IO/Umbrella_to_medit.h>
#include <CGAL/Surface_criteria.h>
#include <CGAL/IO/Debug_print.h>
#include <CGAL/Triangulation_cell_base_with_domain_info_3.h>

//Stretched_traits_3<K,Circular_mf_3<typename K::Point_3> >,

template <typename K, typename Traits,typename Constrain_surface>
class Stretched_delaunay_3:
	public CGAL::Delaunay_triangulation_3<
		Traits,
		CGAL::Triangulation_data_structure_3<
			CGAL::Triangulation_vertex_base_with_info_3<
				unsigned,
				Traits>,
			CGAL::Triangulation_cell_base_with_domain_info_3<
				Traits,
				Constrain_surface>
		>
	>//default tag: No_intersection
{
	public:
	typedef CGAL::Delaunay_triangulation_3<
		Traits,
		CGAL::Triangulation_data_structure_3<
			CGAL::Triangulation_vertex_base_with_info_3<
				unsigned,
				Traits >,
			CGAL::Triangulation_cell_base_with_domain_info_3<
				Traits,
				Constrain_surface>
		>
	>													DT;
	
	typedef typename K::Point_3							Point_3;
	typedef typename K::Object_3						Object_3;
	typedef typename K::Plane_3							Plane_3;
	typedef typename K::FT								FT;
	typedef Stretched_delaunay_3						Self;
	typedef typename DT::Segment						Segment;
	typedef typename DT::Triangle						Triangle;
	typedef typename DT::Tetrahedron	  				Tetrahedron;
	typedef typename DT::Vertex							Vertex;
	typedef typename DT::Cell							Cell;
	typedef typename DT::Cell_handle					Cell_handle;
	typedef typename DT::Edge							Edge;
	typedef typename DT::Facet							Facet;
	typedef typename DT::Vertex_handle                  Vertex_handle;


	typedef typename DT::Finite_cells_iterator			Finite_cells_iterator;
	///<iterator over finite cells

	typedef typename DT::Finite_facets_iterator			Finite_facets_iterator;
	///<iterator over finite facets

	typedef typename DT::Finite_edges_iterator			Finite_edges_iterator;
	///<iterator over finite edges

	typedef typename DT::Finite_vertices_iterator		Finite_vertices_iterator;
	///<iterator over finite vertices

	typedef typename DT::Point_iterator					Point_iterator;
	///<iterator over the points corresponding to the finite vertices of the triangulation.

	typedef typename std::vector<Point_3>::iterator		Vec_point_iterator;

	typedef typename std::map<unsigned,Vertex_handle>	Vertex_map;

	typedef typename std::pair<unsigned,Vertex_handle>	Vertex_pair;

	class Neighbor_cells;

	typedef typename Neighbor_cells::Neighbor_cells_iterator	
														Neighbor_cells_iterator;

	class Neighbor_bounds;

	typedef typename Neighbor_bounds::Neighbor_bounds_iterator	Neighbor_bounds_iterator;


//=========================================================================
//           Creation and initialization functions
//=========================================================================

	///Initialize with assigned traits
	Stretched_delaunay_3(Constrain_surface &cs,
			Traits _traits
			,Point_3 c_p // center point iterator
			,unsigned c_index
			):
		DT(_traits),
		constrain_surface(cs),
		traits(_traits),
		neighbor_cells(*this),
		neighbor_bounds(*this)
	{
		max_squared_stretched_conflict_radius = (2 * constrain_surface.get_bounding_radius());
		max_squared_stretched_conflict_radius *= max_squared_stretched_conflict_radius;

		c_vertex = DT::insert(c_p);
		c_vertex->info() = c_index;
	}



	///Copy constructor
	Stretched_delaunay_3(const Stretched_delaunay_3 & sd)
	  :DT(sd.traits),
	  constrain_surface(sd.constrain_surface),
	  traits(sd.traits),
	  neighbor_cells(*this),
	  neighbor_bounds(*this),
	  c_vertex(sd.get_center())
	{
		c_vertex = insert(sd.get_center()->point(),sd.center_index());
		for(Finite_vertices_iterator vit=sd.finite_vertices_begin();
				vit != sd.finite_vertices_end();
				vit++){
			Vertex_handle v;
			if(vit->point()!=c_vertex->point()){//already inserted
				v=insert(vit->point(),vit->info());
			}
		}
		update_max_squared_stretched_conflict_radius();
	}


	Stretched_delaunay_3 & operator =(Stretched_delaunay_3 sd)
	{
		std::swap((*this),sd);
	  return *this;
	}

//==========================================================================
//			Functions for updating informations
//==========================================================================
	
	void update_max_squared_stretched_conflict_radius(){
		neighbor_cells.update();	  
		FT new_max_squared_stretched_conflict_radius = -1.;
#ifdef DEBUG_CONFLICT_RADIUS
		//std::cout<<"Calculate max_squared radius of umb"<<center_index()<<std::endl;
#endif //DEBUG_CONFLICT_RADIUS
		for(Neighbor_cells_iterator cit = neighbor_cells.begin();
			  cit != neighbor_cells.end();
			  cit ++){
			if(!DT::is_infinite(*cit)){
				FT square_circumradius = 4 * (*cit)->squared_circumradius(traits);
				if(square_circumradius > new_max_squared_stretched_conflict_radius){
					new_max_squared_stretched_conflict_radius = square_circumradius;
				}
#ifdef DEBUG_CONFLICT_RADIUS
		//		std::cout<<square_circumradius<<"  ";
#endif //DEBUG_CONFLICT_RADIUS
			}
		}
		if(new_max_squared_stretched_conflict_radius < 0){
			return;
		}
		max_squared_stretched_conflict_radius = new_max_squared_stretched_conflict_radius;
#ifdef DEBUG_CONFLICT_RADIUS
		//std::cout<<"max squared stretched"<<max_squared_stretched_conflict_radius<<std::endl;
#endif //DEBUG_CONFLICT_RADIUS
	}

	bool is_infinite_vertex(const Vertex_handle &v)const{
		return DT::is_infinite(v);
	}

	class is_infinite_cell{
		public:
		is_infinite_cell(const Stretched_delaunay_3 *_sd):sd(_sd){};

		bool operator() (const Cell_handle &c){
			if(sd->dimension()<3){
				return true;
			}
			return sd->is_infinite(c);
		}

		const Stretched_delaunay_3 *sd;
	};

	is_infinite_cell is_infinite_cell_object()const {
		return is_infinite_cell(this);
	}

	///Update neighbor_vertex_map as the neighboring index to vertex map.
	void update_neighbor_vertices(){
		if(current_neighbor_vertices_center_index == c_vertex->info()){
			return;
		}
		current_neighbor_vertices_center_index = c_vertex->info();

		neighbor_cells.update(c_vertex);
		neighbor_vertex_map.clear();
		for(Neighbor_cells_iterator	  cit=neighbor_cells.begin();
				cit != neighbor_cells.end();
				cit ++){
			int local_index = (*cit)->index(c_vertex);
			for(int i=0;i<4;i++){
				if(i!=local_index){
					neighbor_vertex_map[(*cit)->vertex(i)->info()] = (*cit)->vertex(i);
				}//if not inside
			}//for all vertices
		}//for all cells
	}//update_neighbor_vertices()

	public:
//==========================================================================
//           Insertion functions
//==========================================================================

	/** insert point p if it is in conflict with the
	 * cells around c_vertex
	 */
	Vertex_handle insert_if_conflict(const Point_3 &p, unsigned index){
		Cell_handle loc;
		if(test_conflict(p,loc,false)){//not_from_infinite = false, the point might be on env convex
			return insert(p,index,loc);
		}else{
			return Vertex_handle();
		}
	}

	Vertex_handle insert(const Point_3 &p,unsigned index,Cell_handle loc=Cell_handle()){
		neighbor_cells.invalidate();	  //invalidate the cache of neighboring info
		neighbor_bounds.invalidate();
		current_neighbor_vertices_center_index = -1;

		unsigned before_size = DT::number_of_vertices(); //Check if the point already exist
		Vertex_handle v = DT::insert(p,loc);
		unsigned after_size = DT::number_of_vertices();
		if(before_size == after_size){
		  return Vertex_handle();
		}
		v->info() = index;
		update_max_squared_stretched_conflict_radius();

		return v;
	}

//==================================================================
//           Query functions
//===================================================================
	

	//name not_from_infinite as true if we are sure the point is not on the envenlope convex
	bool test_conflict(const Point_3 &p, Cell_handle &loc,bool not_from_infinite=true){
		
		if(DT::dimension()<3){
			loc = DT::infinite_cell();
			return true;
		}

		typename Traits::Compute_squared_distance_3	csd 
			= traits.compute_squared_distance_3_object();
		if( (not_from_infinite ) && 
				(csd(p,c_vertex->point()) > max_squared_stretched_conflict_radius)){
			return false;
		}//if within the ball of radius max_squared_stretched_conflict_radius.
		neighbor_cells.update(c_vertex,true);
		for(Neighbor_cells_iterator cit=neighbor_cells.begin();
				cit != neighbor_cells.end();
				cit ++){
			//std::cout<<DT::side_of_sphere((*cit),p)<<std::endl;
			if(DT::side_of_sphere((*cit),p)!=CGAL::ON_UNBOUNDED_SIDE){
				loc = (*cit);
				return true;
			}
		}//for all adjacent cells
		return false;
	}


	///@return: true if the given cell is inside of the constrain_surface
	///i.e. if its circumcenter is inside the constrain_surface
	bool cell_is_inside(const Cell_handle &c) const {
		return c->is_inside(*this,constrain_surface,traits);
	}

	Point_3 circumcenter(const Cell_handle &c)const {
		return c->circumcenter(traits);
	}

	// precondition the last update_neighbor_vertices is called to select the neighbor of
	// the c_vertex
	bool has_surface(unsigned i,unsigned j,Facet &facet){
		
		Vertex_handle	vi,vj;
		if(!(index_to_vertex_handle(i,vi)&&index_to_vertex_handle(j,vj))){
		  return false;
		}
		
		neighbor_bounds.update();	//make sure that the neighbor_bounds 
									//are updated to the vetices around center

		for(Neighbor_bounds_iterator bit=neighbor_bounds_begin();
				bit!=neighbor_bounds_end();
				bit ++){
			if(is_the_surface(vi,vj,(*bit))){
				facet = (*bit);
				return true;
			}
		}
		return false;
	}

	bool has_tetrahedra(unsigned i,unsigned j,unsigned k,Cell_handle &c){
		update_neighbor_vertices();
		Vertex_handle v[3];
		if(index_to_vertex_handle(i,v[0])==false){
			return false;
		}
		if(index_to_vertex_handle(j,v[1])==false){
			return false;
		}
		if(index_to_vertex_handle(k,v[2])==false){
			return false;
		}
		return DT::is_cell(v[0],v[1],v[2],c_vertex,c) && cell_is_inside(c);
	}

	///Translate index into vertex_handle, if the vertex lies in the neighboring of center,
	/// 
	/// \return if the vertex of the index requested is found with the neighboring of the center,
	/// return true, otherwise return false
	bool index_to_vertex_handle(unsigned i,Vertex_handle &v){
	  update_neighbor_vertices(); //make sure that the neighbor_vertex_map 
								  //are updated to the vetices around center
	  typename Vertex_map::iterator vp = neighbor_vertex_map.find(i);
	  if(vp == neighbor_vertex_map.end()){
		return false;
	  }else{
		v = vp->second;
		return true;
	  }
	}

	///Judge if the facet f containts vi,vj,center as vertices.
	/// A facet contains the vertex ssi. its mirror facet contains the vertices as well.
	//
	///precondition: f is ajacent to center and vi and vj are not center.
	bool is_the_surface(Vertex_handle vi,Vertex_handle vj,const Facet &f){
		Cell_handle cur_cell = f.first;
		Facet m_facet=mirror_facet(f);
		Cell_handle cur_mirror_cell = m_facet.first;
		if((cur_cell->has_vertex(vi))&&(cur_cell->has_vertex(vj))&&
				(cur_mirror_cell->has_vertex(vi))&&(cur_mirror_cell->has_vertex(vj))){
			return true;
		}
		return false;
	}


	///returns if the center encroaches the existing surface of the umbrellaset
	/// If an encroachment happens, index[3] will store the indices of the facet
	template<typename Umbrella_set>
	bool center_encroach_surface(Umbrella_set *umbs,unsigned index[3]){
		if(on_boundary()){
			return false;
		}
		for(Neighbor_cells_iterator cit = neighbor_cells_begin();
				cit != neighbor_cells_end();
				cit ++){
			assert(cell_is_inside(*cit));
			int center_vertex_cell_index = (*cit)->index(c_vertex);
			unsigned k=0;
			for(int i=0;i<4;i++){
				if(i!=center_vertex_cell_index){
					index[k] = (*cit)->vertex(i)->info();
					k++;
				}
			}
			if(umbs->is_surface(index[0],index[1],index[2]) 
					&& (!is_Gabriel(*cit,center_vertex_cell_index))){
#ifdef DEBUG_VOLUME_REFINEMENT_POINT
					std::cout<<"Encrochment to existing surfaces"<<index[0]<<" "
						<<index[1]<<" "<<index[2]<<" "<<std::endl;
#endif //DEBUG_VOLUME_REFINEMENT_POINT
					return true;
			}
		}//end of for
		return false;
	}


	/// judge if p encroaches some constraint face or not.
	/// if p encroaches a triangle, name facet to the facet to be modified.
	bool check_encroachment(const Point_3 &p,Facet &facet){
		for(Neighbor_bounds_iterator bit = neighbor_bounds_begin();
				bit != neighbor_bounds_end();
				bit ++){
			if(facet_encroached_by(*bit,p)){
				facet = *bit;
				return true;
			}
		}
		return false;
	}

	bool facet_encroached_by(const Facet &facet,const Point_3 &p){
		typename Traits::Encroachment_face_3 ef = traits.encroachment_face_3_object();
		Point_3 pf[3];
		unsigned k=0;
		for(int i=0; i<4; i++){
			if(i != facet.second){
				pf[k] = facet.first->vertex(i)->point();
				k++;
			}
		}//end of for, now, pf[3] contains the three points of facet.
		return ef(pf[0],pf[1],pf[2],p);
	}


	/// \return The circumcenter of 3 points projected to the surface of the constrain_surface.
	/// Given indices are supposed to be ajacent to the center index
	Point_3 refine_surface_point(unsigned i,unsigned j){
		update_neighbor_vertices();
		if(neighbor_vertex_map.find(i) == neighbor_vertex_map.end()){
			std::cout<<"Refine_surface_point called by umb["
				<<center_index()<<"with i="<<i<<" and j = "<<j<<std::endl;
		}
		assert(neighbor_vertex_map.find(i)!=neighbor_vertex_map.end());
		assert(neighbor_vertex_map.find(j)!=neighbor_vertex_map.end());
		
		Facet facet;
		has_surface(i,j,facet);
		assert(has_surface(i,j,facet));

		Object_3 intersect_object = constrain_surface.intersection(DT::dual(facet));
		Point_3 to_be_returned = CGAL::object_cast<Point_3>(intersect_object);

		return to_be_returned;
	}

//========================================================================
//             output functions
//========================================================================	

	void output_triangulation_to_medit(std::string name){
		std::stringstream filename;
		filename<<"data/"<<name<<traits.get_index() <<".mesh";
		typename std::ofstream of(filename.str().c_str());
		if(of){
			umbrella_to_medit<Stretched_delaunay_3>(this,of);
		}
	}

	//=====================================================================
	//	           value retrieval functions
	//=====================================================================

	FT get_absolute_squared_conflict_radius()const{
		FT absolute_distortion = traits.Absolute_distortion();
		FT to_be_returned = 
			absolute_distortion * absolute_distortion * max_squared_stretched_conflict_radius;
#ifdef DEBUG_CONFLICT_RADIUS
		//std::cout<<"get_absolute_squared_conflict_radius of umb["<<
		//	center_index()<<"] returns "<<to_be_returned<<std::endl;
#endif //DEBUG_CONFLICT_RADIUS
		return to_be_returned;
	}

	unsigned center_index() const{
		return c_vertex->info();
	}

	Vertex_handle get_center()const {
		return c_vertex;
	}

	const Traits &get_traits()const{
		return traits;
	}

	const Traits &get_traits_ref(){
		return traits;
	}

	const Constrain_surface &get_constrain_surface() const {
		return constrain_surface;
	}

	bool on_boundary(){
		return neighbor_bounds.on_boundary();
	}
//=====================================================================
//                  Structure
//=====================================================================


	class Neighbor_cells{
		public:
		Neighbor_cells(const Stretched_delaunay_3 &_umb):
			umb(_umb),
			current_neighbor_cells_include_infinite(false),
			current_neighbor_cells_center_index(-1)
		{};

		void invalidate(){
			current_neighbor_cells_include_infinite = false;
			current_neighbor_cells_center_index = -1;
		}

		typedef typename std::list<Cell_handle>::iterator	  Neighbor_cells_iterator;

		Neighbor_cells_iterator begin(){
			return n_list.begin();
		}

		Neighbor_cells_iterator end(){
			return n_list.end();
		}

	///Update neighbor_cells as the neighbor finite cells around 
	/// the vertex; defalt updates the information
	/// around the center;
		void update(Vertex_handle v=Vertex_handle(),bool include_infinite=false){

			if(v==Vertex_handle()){
				v = umb.get_center();
			}
			if(current_neighbor_cells_center_index == v->info() && 
					(current_neighbor_cells_include_infinite == include_infinite)){
				return;
			}else{
				current_neighbor_cells_include_infinite = include_infinite;
				current_neighbor_cells_center_index = v->info();
			}

			n_list.clear();
			typename std::back_insert_iterator<std::list<Cell_handle> >		
				neighbor_cells_insert_iterator(n_list);

			umb.incident_cells(v,neighbor_cells_insert_iterator);
#ifdef DEBUG_NEIGHBORS
			std::cout<<"Neighbor cells update to center "<<v->info()<<" in umb "<<umb.center_index()
				<<"\n total vertex_number = "<<umb.number_of_vertices()
				<<"\n size = "<<n_list.size()<<std::endl;
#endif //DEBUG_NEIGHBORS
			if(!include_infinite){
				n_list.remove_if(umb.is_infinite_cell_object());
			}

		}


		const Stretched_delaunay_3			&umb;
		std::list<Cell_handle>				n_list;
		bool								current_neighbor_cells_include_infinite;
		unsigned							current_neighbor_cells_center_index;
															///< Note down the current calculation 
															///  center information for neighbor_cells, 
	};

	class Neighbor_bounds{
		
		const Stretched_delaunay_3			&umb;
		std::list<Facet>					n_list;
		unsigned							current_neighbor_bounds_center_index;
															///< Note down the current calculation 
															///  center information for neighbor_cells, 
		
		public:
		
		typedef typename std::list<Facet>::iterator	  Neighbor_bounds_iterator;
		
		Neighbor_bounds(const Stretched_delaunay_3 &_umb):
			umb(_umb),
			current_neighbor_bounds_center_index(-1)
		{};

		void invalidate(){
			if(n_list.size()!=0){
				current_neighbor_bounds_center_index = -1;
			}
		}

		unsigned size(){
			return n_list.size();
		}

		Neighbor_bounds_iterator begin(){
			return n_list.begin();
		}

		Neighbor_bounds_iterator end(){
			return n_list.end();
		}

		///Update neighbor_cells as the neighbor finite cells around 
		/// the vertex; defalt updates the information
		/// around the center;
		void update(Vertex_handle v=Vertex_handle()){
			umb.update_neighbor_cells(v);
			if(v==Vertex_handle()){
				v=umb.get_center();
			}
			if(current_neighbor_bounds_center_index == v->info()){
				return;
			}
			current_neighbor_bounds_center_index = v->info();
			n_list.clear();
			for(Neighbor_cells_iterator	  cit=umb.neighbor_cells_begin();
					cit != umb.neighbor_cells_end();
					cit ++){
				int local_index = (*cit)->index(umb.get_center());

				if(umb.cell_is_inside(*cit)){
					for(int i=0;i<4;i++){
						if(i!=local_index &&		//only checks the neighboring facets
								(!umb.cell_is_inside((*cit)->neighbor(i)))){
							n_list.push_back(Facet((*cit),i));
						}
					}
				}//inside constrain
			}
#ifdef DEBUG_NEIGHBORS
			std::cout<<"Neighbor bounds update to center "<<v->info()<<" in umb "<<umb.center_index()
				<<"\nsize ="<<n_list.size()<<std::endl;
#endif //DEBUG_NEIGHBORS
		}//update_neighbor_bounds

		bool on_boundary(){
			update();
			return (n_list.size()!=0);
		}
	};

	Neighbor_cells_iterator neighbor_cells_begin(Vertex_handle v=Vertex_handle())const{
		update_neighbor_cells(v);
		return neighbor_cells.begin();
	}

	Neighbor_cells_iterator neighbor_cells_end()const{
		return neighbor_cells.end();
	}


	Neighbor_bounds_iterator neighbor_bounds_begin()const{
		update_neighbor_bounds();
		return neighbor_bounds.begin();
	}

	Neighbor_bounds_iterator neighbor_bounds_end()const{
		return neighbor_bounds.end();
	}

	void update_neighbor_bounds()const{
		neighbor_bounds.update();
	}

	void update_neighbor_cells(Vertex_handle v)const{
		neighbor_cells.update(v);
	}

	//=====================================================================
	//                  members
	//=====================================================================

	public:
	Constrain_surface						&constrain_surface;		  ///<	constrain_surface
	const Traits							traits;			  ///<	traits
	mutable Neighbor_cells					neighbor_cells;
	mutable Neighbor_bounds					neighbor_bounds;
	Vertex_map								neighbor_vertex_map;  ///<	map from unsigned(index) 
																  ///	  to Vertex_handle
	
	private:
	Vertex_handle							c_vertex;		  ///<	center vertex

	FT									max_squared_stretched_conflict_radius;
															  ///<	maximum of 2 * circumradius among
															  ///surrounding finite cells

	unsigned							current_neighbor_vertices_center_index;
															///< Note down the current calculation 
															///  center information for neighbor_vertice_map, 
};// end of class Stretched_delaunay_3


#endif //STRETCHED_DELAUNAY_3_H

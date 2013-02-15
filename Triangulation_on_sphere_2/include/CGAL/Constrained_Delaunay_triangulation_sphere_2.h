

#ifndef CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_SPHERE_2_H
#define CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_SPHERE_2_H

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_sphere_2.h>
#include <CGAL/Delaunay_triangulation_sphere_2.h>
#include <CGAL/Constrained_triangulation_face_base_sphere_2.h>
#include <CGAL/Triangulation_sphere_line_face_circulator_2.h>
#include <set>
#include <CGAL/iterator.h>


namespace CGAL {
template <class Gt,
	      class Tds=Triangulation_data_structure_2 <
	                Triangulation_vertex_base_2<Gt>,
	                Constrained_triangulation_face_base_sphere_2<Gt> > >
	
	
class Constrained_Delaunay_triangulation_sphere_2
	: public Delaunay_triangulation_sphere_2<Gt, Tds>
{
public:
	typedef Constrained_Delaunay_triangulation_sphere_2<Gt, Tds>  Self;
	typedef Gt						Geom_traits;
	typedef Delaunay_triangulation_sphere_2<Gt, Tds>    DTS;
	
	typedef typename DTS::Edge Edge;
	typedef typename DTS::Vertex Vertex;
	typedef typename DTS::Vertex_handle Vertex_handle;
	typedef typename DTS::Face_handle Face_handle;
	typedef typename DTS::size_type size_type;
	typedef typename DTS::Locate_type Locate_type;
	typedef typename DTS::All_faces_iterator All_faces_iterator;
	typedef typename DTS::Solid_faces_iterator Solid_faces_iterator;
	typedef typename DTS::Solid_edges_iterator Solid_edges_iterator;
	typedef typename DTS::All_vertices_iterator   All_vertices_iterator;
	typedef typename DTS::Face_circulator Face_circulator;
	typedef typename DTS::Edge_circulator Edge_circulator;
	typedef typename DTS::Vertex_circulator Vertex_circulator;
	typedef Triangulation_sphere_line_face_circulator_2<Self> Line_face_circulator;//not excisting yet
	

	typedef typename Geom_traits::Point_2      Point;
	typedef std::pair<Point,Point>             Constraint;
	typedef std::list<Edge>                    List_edges;
	typedef std::list<Face_handle>             List_faces;
	typedef std::list<Vertex_handle>           List_vertices;

	//for meshing
	typedef typename DTS::Solid_edges_iterator  Finite_edges_iterator;
	typedef typename DTS::Solid_faces_iterator Finite_faces_iterator;
	typedef typename DTS::All_vertices_iterator   Finite_vertices_iterator;
	 typedef Tag_false                           Constraint_hierarchy_tag;
	
	
	//-----Constructions
	Point circumcenter(Face_handle  f) const; 
	Point circumcenter(const Point& p0, 
					   const Point& p1, 
					   const Point& p2) const;
	
	
	Face_handle locate(const Point & p, Face_handle start = Face_handle()) const
	{
		Locate_type lt;
		int li;
		return DTS::locate(p, lt, li, start);
	}
	
	Solid_edges_iterator finite_edges_begin()const{
		return DTS::solid_edges_begin();
	}
	
	Solid_edges_iterator finite_edges_end(){
		return DTS::solid_edges_end();
	}
	
	Solid_faces_iterator finite_faces_begin(){
		return DTS::solid_faces_begin();
	}
	
	Solid_faces_iterator finite_faces_end(){
		return DTS::solid_faces_end();
	}
	
	All_vertices_iterator finite_vertices_end(){
		return vertices_end();
	}
	All_vertices_iterator finite_vertices_begin(){
		return vertices_begin();
	}
	
	
	//is_infinite(...) is needed for meshing
	bool is_infinite(Vertex_handle v)const{
		return false;
	}
	
	bool is_infinite(Face_handle f) const{
		//return this->is_ghost(f);
		return f->is_ghost();
	}
	
	Face_handle infinite_face()const{
		return DTS::_ghost;}
	
	
#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
	using DTS::number_of_vertices;
	using DTS::cw;
	using DTS::ccw;
	using DTS::dimension;
	using DTS::all_faces_begin;
	using DTS::all_faces_end;
	using DTS::solid_faces_begin;
	using DTS::solid_faces_end;
	using DTS::vertices_begin;
	using DTS::vertices_end;
	using DTS::all_edges_begin;
	using DTS::all_edges_end;
	using DTS::power_test;
	using DTS::collinear_between;
	using DTS::incident_edges;
	using DTS::test_dim_down;
	using DTS::test_dim_up;
	using DTS::make_hole;
	using DTS::delete_vertex;
	using DTS::delete_face;
	using DTS::create_face;
	using DTS::incident_faces;
	using DTS::locate;
	using DTS::is_edge;//instead of includes_edge
	using DTS::insert_first;
	using DTS::insert_second;
	using DTS::insert_outside_affine_hull_regular;
	using DTS::insert_cocircular;
	using DTS::update_ghost_faces;
	using DTS::test_conflict;
	using DTS::geom_traits;
	using DTS::FACE;
	using DTS::VERTEX;
	using DTS::NOT_ON_SPHERE;
	using DTS::OUTSIDE_CONVEX_HULL;
	using DTS::OUTSIDE_AFFINE_HULL;
	using DTS::TOO_CLOSE;
	
	
#endif
	class Less_edge;
	typedef std::set<Edge,Less_edge> Edge_set;
	
	//CONSTRUCTORS /DESTRUCTORS
	Constrained_Delaunay_triangulation_sphere_2 (const Gt& gt = Gt()):DTS(gt){}
			
	Constrained_Delaunay_triangulation_sphere_2(const Self& cdt)
    : DTS(cdt) {}
	
    virtual ~Constrained_Delaunay_triangulation_sphere_2() {}
	
	
	
	//QUERY
	bool is_constrained(Edge e) const;
	bool is_valid(bool verbose = false, int level = 0) const;
	/*bool includes_edge(Vertex_handle va, Vertex_handle vb,
						 Vertex_handle& vbb,
						 Face_handle& fr, int & i) const;*/
		
	//FLIIPS
	bool is_flipable(Face_handle f, int i, bool perturb = true) const;
	void flip(Face_handle& f, int i);
	void propagating_flip(Face_handle& f,int i);
	void propagating_flip(List_edges & edges);
	
	
	// INSERTION-REMOVAL
	Vertex_handle insert(const Point & a, Face_handle start = Face_handle());
	Vertex_handle insert(const Point& p,
						 Locate_type lt,
						 Face_handle loc, int li );
	void insert_constraint(const Point& a, const Point& b);//not yet implemented
	void insert_constraint(Vertex_handle  vaa, Vertex_handle vbb);
	
	void insert(Point a, Point b) { insert_constraint(a, b);}
	void insert(Vertex_handle va, Vertex_handle  vb) {insert_constraint(va,vb);}
	
	bool find_intersected_faces(Vertex_handle vaa,
								Vertex_handle vbb,
								List_faces & intersected_faces,
								List_edges & list_ab, 
								List_edges & list_ba,
								Vertex_handle & vi) ;
	
	
	Face_handle locate(const Point& p,Locate_type& lt,int& li, Face_handle start) const;
	
	class Less_edge 
	:  public std::binary_function<Edge, Edge, bool>
	{
	public:
		Less_edge() {}
		bool operator() (const Edge& e1, const Edge& e2) const
		{
			int ind1=e1.second, ind2=e2.second;
			return( (&(*e1.first) < &(*e2.first))
				   || ( (&(*e1.first) == &(*e2.first)) && (ind1 < ind2)));
		} 
		
		
	};
	
protected:
	virtual void triangulate_hole(List_faces& intersected_faces,
								  List_edges& conflict_boundary_ab,
								  List_edges& conflict_boundary_ba);
	
	void triangulate_hole(List_faces& intersected_faces,
						  List_edges& conflict_boundary_ab,
						  List_edges& conflict_boundary_ba,
						  List_edges& new_edges);
	
	void triangulate_half_hole(List_edges & list_edges, 
							   List_edges & new_edges);
	
	//CONSTRAINED
public:	
	void mark_constraint(Face_handle fr, int i);
	void update_constraints_opposite(Vertex_handle va);
	void update_constraints_incident(Vertex_handle v, Vertex_handle v1, Vertex_handle v2);
	void clear_constraints_incident(Vertex_handle va);	
	
	template <class OutputItFaces, class OutputItBoundaryEdges> 
	std::pair<OutputItFaces,OutputItBoundaryEdges>
	get_conflicts_and_boundary(const Point  &p, OutputItFaces fit, OutputItBoundaryEdges eit, 
							   Face_handle start = Face_handle() ) const
	{
		
			CGAL_triangulation_precondition( dimension() == 2);
			int li;
			Locate_type lt;
			Face_handle fh = locate(p,lt,li, start);
		
					
		*fit++ = fh; 
		std::pair<OutputItFaces,OutputItBoundaryEdges>
		pit = std::make_pair(fit,eit);
		if (lt ==TOO_CLOSE || lt == VERTEX) 
			return pit;
		
		pit = propagate_conflicts(p,fh,0,pit);
		pit = propagate_conflicts(p,fh,1,pit);
		pit = propagate_conflicts(p,fh,2,pit);
		
		return std::make_pair(fit,eit);
		
	} 
	
	
	template <class OutputItFaces, class OutputItBoundaryEdges> 
	std::pair<OutputItFaces,OutputItBoundaryEdges>
	propagate_conflicts (const Point  &p,
						 Face_handle fh, 
						 int i,
						 std::pair<OutputItFaces,OutputItBoundaryEdges>  pit)  const {
		Face_handle fn = fh->neighbor(i);
		if ( fh->is_constrained(i) || ! test_conflict(p,fn)) {
			*(pit.second)++ = Edge(fn, fn->index(fh));
		} else {
			*(pit.first)++ = fn;
			int j = fn->index(fh);
			pit = propagate_conflicts(p,fn,ccw(j),pit);
			pit = propagate_conflicts(p,fn,cw(j), pit);
		}
		return pit;
	}
	
	
	
	
	
	template <class OutputItFaces>
	OutputItFaces
	propagating_flip(List_edges & edges, 
					 OutputItFaces out = Emptyset_iterator()) {
		// makes the triangulation Delaunay by flipping 
		// List edges contains an initial list of edges to be flipped
		// Precondition : the output triangulation is Delaunay if the list 
		// edges contains all edges of the input triangulation that need to be
		// flipped (plus possibly others)
		int i, ii, indf, indn;
		Face_handle ni, f,ff;
		Edge ei,eni; 
		 typename Self::Edge_set edge_set;
		 typename Self::Less_edge less_edge;
		Edge e[4];
		typename List_edges::iterator itedge=edges.begin();
		
		// initialization of the set of edges to be flip
		while (itedge != edges.end()) {
			f=(*itedge).first;
			i=(*itedge).second;
			if (is_flipable(f,i)) {
				eni=Edge(f->neighbor(i),this->mirror_index(f,i));
				if (less_edge(*itedge,eni)) edge_set.insert(*itedge);
				else edge_set.insert(eni);
			}
			++itedge;
		}
		
		// flip edges and updates the set of edges to be flipped
		while (!(edge_set.empty())) {
			f=(*(edge_set.begin())).first;
			indf=(*(edge_set.begin())).second;
			
			// erase from edge_set the 4 edges of the wing of the edge to be
			// flipped (edge_set.begin) , i.e. the edges of the faces f and
			// f->neighbor(indf) that are distinct from the edge to be flipped
			
			ni = f->neighbor(indf); 
			indn=this->mirror_index(f,indf);
			ei= Edge(f,indf);
			edge_set.erase(ei);
			e[0]= Edge(f,cw(indf));
			e[1]= Edge(f,ccw(indf));
			e[2]= Edge(ni,cw(indn));
			e[3]= Edge(ni,ccw(indn));
			
			for(i=0;i<4;i++) { 
				ff=e[i].first;
				ii=e[i].second;
				eni=Edge(ff->neighbor(ii),this->mirror_index(ff,ii));
				if (less_edge(e[i],eni)) {edge_set.erase(e[i]);}
				else { edge_set.erase(eni);} 
			} 
			
			// here is the flip 
			*out++ = f;
			*out++ = f->neighbor(indf);
			flip(f, indf); 
			
			
			//insert in edge_set the 4 edges of the wing of the edge that
			//have been flipped 
			e[0]= Edge(f,indf);
			e[1]= Edge(f,cw(indf));
			e[2]= Edge(ni,indn);
			e[3]= Edge(ni,cw(indn));
			
			for(i=0;i<4;i++) { 
				ff=e[i].first;
				ii=e[i].second;
				if (is_flipable(ff,ii)) {
					eni=Edge(ff->neighbor(ii),this->mirror_index(ff,ii));
					if (less_edge(e[i],eni)) { 
						edge_set.insert(e[i]);}
					else {
						edge_set.insert(eni);} 
				}
			} 
		}
		return out;
	}
	
	
	
	
	
	template<class EdgeIt>
	Vertex_handle star_hole( const Point& p, EdgeIt edge_begin, EdgeIt edge_end) {
		std::list<Face_handle> empty_list;
		return star_hole(p, edge_begin, edge_end, empty_list.begin(), empty_list.end());
	}
	
	template<class EdgeIt, class FaceIt>
	Vertex_handle star_hole( const Point& p, EdgeIt edge_begin, EdgeIt edge_end,
							FaceIt face_begin, FaceIt face_end)
	{
		//Vertex_handle v = this->_tds.star_hole(edge_begin, edge_end);
		Vertex_handle v =  this->_tds.star_hole( edge_begin, edge_end, 
												face_begin, face_end);
		v->set_point(p);
		// restore constraint status for new faces.
		
		int vindex;
		Face_handle fh;
		int ih;
		Face_circulator fc = incident_faces(v), done(fc);
		do {
			fc->set_in_conflict_flag(0);
			vindex = fc->index(v);
			fc->set_constraint(cw(vindex), false);
			fc->set_constraint(ccw(vindex), false);
			fh = fc->neighbor(vindex);
			ih = this->mirror_index(fc,vindex);
			fc->set_constraint(vindex, fh->is_constrained(ih));
		} while (++fc != done);
		//DTS::delete_faces(face_begin, face_end);
		DTS::update_ghost_faces(v);
		
		return v;
	}
	
	
	
public:
	template < class InputIterator >
	int insert(InputIterator first, InputIterator last)
	{
		int n = number_of_vertices();
		
		std::vector<Point> points (first, last);
		std::random_shuffle (points.begin(), points.end());
		spatial_sort (points.begin(), points.end());
		
		Face_handle hint;
		Vertex_handle v;
		for (typename std::vector<Point>::const_iterator p = points.begin(),end = points.end(); p != end; ++p){
			v = insert (*p, hint);
			if( v != Vertex_handle())
				hint = v->face();
		}
		
		return number_of_vertices() - n;
	}
	
};	//end class
	
	
	
	
	
//QUERY
	template < class Gt, class Tds>
	inline  bool 
	Constrained_Delaunay_triangulation_sphere_2<Gt,Tds>::
	is_constrained(Edge e) const
	{
		return (e.first)->is_constrained(e.second);
	}
	
	template < class Gt, class Tds >
	inline  bool 
	Constrained_Delaunay_triangulation_sphere_2<Gt,Tds>::
	is_valid(bool verbose, int level) const
	{
		bool result = DTS::is_valid(verbose,level);
		for( All_faces_iterator it = all_faces_begin(); 
			it != all_faces_end() ; it++) {
			for(int i=0; i<3; i++) {
				Face_handle n = it->neighbor(i);
				result = result && 
				it->is_constrained(i) == n->is_constrained(n->index(it));
			}
		}
		
		Solid_faces_iterator fit= solid_faces_begin();
		for (; fit != solid_faces_end(); fit++) {
			for(int i=0;i<3;i++) {
				result = result && !is_flipable(fit,i, false);
				CGAL_triangulation_assertion( result );
			}
		}
		
		return result;
	}
	
		
  
						
	
	//FLIPS
	template < class Gt, class Tds>
	bool 
	Constrained_Delaunay_triangulation_sphere_2<Gt,Tds>::
	is_flipable(Face_handle f, int i, bool perturb) const
	// determines if edge (f,i) can be flipped 
	{
		Face_handle ni = f->neighbor(i); 
		if (f->is_constrained(i)) return false;
		return (power_test(ni, f->vertex(i)->point(), perturb) 
				== ON_POSITIVE_SIDE);
	}
	
 
	
//INSERTION
	
	template < class Gt, class Tds>
	inline
	typename Constrained_Delaunay_triangulation_sphere_2<Gt,Tds>::Vertex_handle
	Constrained_Delaunay_triangulation_sphere_2<Gt,Tds>::
	insert(const Point& a, Face_handle start)
	{
		Face_handle loc;
		int li;
		Locate_type lt;
		loc = locate(a, lt, li, start);
		switch (lt){
			case NOT_ON_SPHERE: 		
				return Vertex_handle();
			case TOO_CLOSE:
				return loc->vertex(li);
			case VERTEX:{
				if(number_of_vertices()==1)
					return vertices_begin();
				return (loc->vertex(li));
			}
				
			case DTS::EDGE:
			{
				bool constrained = loc->is_constrained(li);
				if (constrained){
					Vertex_handle v,v1, v2;
					v1=loc->vertex(ccw(li)); //endpoint of the constraint
					v2=loc->vertex(cw(li)); // endpoint of the constraint
					v = this->_tds.insert_in_edge(loc, li);
					v->set_point(a);
					update_constraints_incident(v, v1, v2);
					update_constraints_opposite(v);
					return v;
				}
			}
			default:
				return insert(a, lt, loc, li);
		}
		
	}
	
	template < class Gt, class Tds >
	typename Constrained_Delaunay_triangulation_sphere_2<Gt,Tds>::Vertex_handle
	Constrained_Delaunay_triangulation_sphere_2<Gt,Tds>::
	insert(const Point& p, Locate_type lt, Face_handle loc, int li)
	{
		Vertex_handle v;
		
	switch(dimension()){
		case -2 : 
			return insert_first(p);
		case -1 :
			return insert_second(p);
		case 0 :
			return insert_outside_affine_hull_regular(p);
		case 1 : {
			if (test_dim_up(p)){
				Face_handle f = all_edges_begin()->first;
				Vertex_handle v1=f->vertex(0);
				Vertex_handle v2=f->vertex(1);
				Vertex_handle v3=f->neighbor(0)->vertex(1);
				return insert_outside_affine_hull_regular(p);
			} 
			else 
				return insert_cocircular(p,lt,loc);	
		}
			
		case 2:{
			std::vector<Face_handle> faces;
			std::vector<Edge> edges;
			faces.reserve(32);
			edges.reserve(32);
			
			get_conflicts_and_boundary(p, std::back_inserter(faces), std::back_inserter(edges), loc);
			v =star_hole(p,edges.begin(), edges.end());
			v->set_point(p);
			this->delete_faces(faces.begin(), faces.end());
			
			if( lt != FACE )
				update_ghost_faces(v);
			
			update_constraints_opposite(v);
			//clear_constraints_incident(v);
			return v;		
		}
	 }//end switch
		CGAL_assertion(false);
		return v;
 }	

	template < class Gt, class Tds >
	inline void
	Constrained_Delaunay_triangulation_sphere_2<Gt,Tds>::
	insert_constraint(const Point& a, const Point& b)
	// the algorithm first inserts a and b, 
	// and then forces the constraint [va,vb]
	{
		//Vertex_handle va= virtual_insert(a);
		//Vertex_handle vb= virtual_insert(b);
		Vertex_handle va =insert(a);
		Vertex_handle vb=insert(b);
		
		if ( va != vb && va!=Vertex_handle() && vb!=Vertex_handle())   
			insert_constraint(va,vb);
	}
	
	template <class Gt, class Tds >
	inline void
	Constrained_Delaunay_triangulation_sphere_2<Gt,Tds>::
	insert_constraint(Vertex_handle  vaa, Vertex_handle vbb)
	// forces the constrained [va,vb]
	// [va,vb] will eventually be splitted into several edges
	// if a vertex vc of t lies on segment ab
	// of if ab intersect some constrained edges
	{
		CGAL_triangulation_precondition( vaa != vbb);
		CGAL_triangulation_precondition (vaa !=Vertex_handle());
		CGAL_triangulation_precondition (vbb != Vertex_handle());
		if(dimension()<=1){
			Face_handle fh=all_edges_begin()->first;
			mark_constraint(fh,2);
			return;
		}
		Vertex_handle vi;      
		Face_handle fr;
		int i;
		
		//first step :constrained allready an edge?
		if(is_edge(vaa,vbb,fr,i)) {    
			mark_constraint(fr,i);
			return;
		}
		
		List_faces intersected_faces;
		List_edges conflict_boundary_ab, conflict_boundary_ba;
		
		
		//following bool not needed
		bool intersection  = find_intersected_faces( vaa, vbb,
													intersected_faces,
													conflict_boundary_ab,
													conflict_boundary_ba,
													vi);
		
		CGAL_triangulation_precondition(!intersection);
		//no intersection
		triangulate_hole(intersected_faces,
						 conflict_boundary_ab,
						 conflict_boundary_ba);
		
		if (vi != vbb) {
			insert_constraint(vi,vbb); 
		}
		return;
		
	}
	
	
	template <class Gt, class Tds >
	bool
	Constrained_Delaunay_triangulation_sphere_2<Gt,Tds>::
	find_intersected_faces(Vertex_handle vaa,
						   Vertex_handle vbb,
						   List_faces & intersected_faces,
						   List_edges & list_ab, 
						   List_edges & list_ba,
						   Vertex_handle & vi) 
	// vi is set to the first vertex of the triangulation on [vaa,vbb].
	
	//if no intersection:
	// intersected_faces contains the list of faces intersected by [va,vi]
	// list_ab and list_ba represents the boundary of the union
	// of the intersected faces oriented cw
	// list_ab consists of the edges from vaa to vi (i.e. on the left of a->b)
	// list_ba    "         "        from vi to vaa (i.e. on the right of a->b)
	{
		const Point& aa = vaa->point();
		const Point& bb = vbb->point();
		Line_face_circulator current_face=Line_face_circulator(vaa, this, bb);
		int ind=current_face->index(vaa);
				
		Face_handle lf= current_face->neighbor(ccw(ind)); 
		Face_handle rf= current_face->neighbor(cw(ind));
		Orientation orient;
		Face_handle previous_face;
		Vertex_handle current_vertex;	
		
		list_ab.push_back(Edge(lf, lf->index(current_face)));
		list_ba.push_front(Edge(rf, rf->index(current_face)));
		intersected_faces.push_front(current_face);
		
		// initcd
		previous_face=current_face; 
		++current_face;
		ind=current_face->index(previous_face);  
		current_vertex=current_face->vertex(ind); 
		
		
		// loop over triangles intersected by ab
		bool done = false;
		while (current_vertex != vbb && !done)  { 
			orient = this->orientation(aa,bb,current_vertex->point());
			int i1, i2;
			switch (orient) {
				case COLLINEAR :  
					done = true; // current_vertex is the new endpoint
					break;
				case LEFT_TURN :
				case RIGHT_TURN :
					if (orient == LEFT_TURN) {
						i1 = ccw(ind) ; //index of second intersected edge of current_face
						i2 = cw(ind); //index of non intersected edge of current_face
					}
					else {
						i1 = cw(ind) ; //index of second intersected edge of current_face
						i2 = ccw(ind); //index of non intersected edge of current_face
					}
					if(current_face->is_constrained(i1)) {
						return true;
					}
					else {
						lf= current_face->neighbor(i2);
						intersected_faces.push_front(current_face);
						if (orient == LEFT_TURN) 
							list_ab.push_back(Edge(lf, lf->index(current_face)));
						else // orient == RIGHT_TURN
							list_ba.push_front(Edge(lf, lf->index(current_face)));
						previous_face=current_face;
						++current_face;
						ind=current_face->index(previous_face); 
						current_vertex=current_face->vertex(ind);
					}
					break;
			}
		}
		
		// last triangle 
		vi = current_vertex;
		intersected_faces.push_front(current_face);
		lf= current_face->neighbor(cw(ind));
		list_ab.push_back(Edge(lf, lf->index(current_face))); 
		rf= current_face->neighbor(ccw(ind));
		list_ba.push_front(Edge(rf, rf->index(current_face)));
		return false;
	}
	
	template < class Gt, class Tds >
	void 
	Constrained_Delaunay_triangulation_sphere_2<Gt,Tds>::
	triangulate_hole(List_faces& intersected_faces,
					 List_edges& conflict_boundary_ab,
					 List_edges& conflict_boundary_ba)
	{
		List_edges new_edges;
		triangulate_hole(intersected_faces,
						 conflict_boundary_ab,
						 conflict_boundary_ba,
						 new_edges);
		propagating_flip(new_edges);
		return;
	}
	
	
	template < class Gt, class Tds >
	void
	Constrained_Delaunay_triangulation_sphere_2<Gt,Tds>::
	triangulate_hole(List_faces& intersected_faces,
					 List_edges& conflict_boundary_ab,
					 List_edges& conflict_boundary_ba,
					 List_edges& new_edges)
	// triangulate the hole limited by conflict_boundary_ab
	// and conflict_boundary_ba
	// insert the new edges in new-edges 
	// delete the faces of intersected_faces
	{
		if ( !conflict_boundary_ab.empty() ) {
			triangulate_half_hole(conflict_boundary_ab, new_edges);
			triangulate_half_hole(conflict_boundary_ba, new_edges);
			
			// the two faces that share edge ab are neighbors
			// their common edge ab is a constraint
			Face_handle fr,fl;
			fl=(*conflict_boundary_ab.begin()).first;
			fr=(*conflict_boundary_ba.begin()).first;
			fl->set_neighbor(2, fr);
			fr->set_neighbor(2, fl);
			fl->set_constraint(2, true);
			fr->set_constraint(2, true);
			
			// delete intersected faces
			while( ! intersected_faces.empty()) {
				fl = intersected_faces.front();
				intersected_faces.pop_front();
				delete_face(fl);
			}
		}
	}
	
	template < class Gt, class Tds >
	void
	Constrained_Delaunay_triangulation_sphere_2<Gt,Tds>::
	triangulate_half_hole(List_edges & list_edges,  List_edges & new_edges)
	// triangulates the  polygon whose boundary consists of list_edges
	// plus the edge ab joining the two endpoints of list_edges
	// the orientation of the polygon (as provided by list_edges) must
	// be cw
	// the edges of list_edges are assumed to be edges of a
	// triangulation that will be updated by the procedure
	// the edges that are created are put in list new_edges
	// takes linear time
	{
		Vertex_handle va; // first vertex of list_edges 
		Face_handle newlf;
		Face_handle n1,n2,n;
		int ind1, ind2,ind;
		Orientation orient;
		
		typename List_edges::iterator current, next, tempo;
		current=list_edges.begin();
		
		va=((*current).first)->vertex(ccw((*current).second));
		next=current; 
		++next;
		
		do
		{
			n1=(*current).first;
			ind1=(*current).second;
			// in case n1 is no longer a triangle of the new triangulation
			if ( n1->neighbor(ind1) != Face_handle() ) {
				n=n1->neighbor(ind1);
				//ind=this->mirror_index(n1,ind1); 
				// mirror_index does not work in this case
				ind = cw(n->index(n1->vertex(cw(ind1))));
				n1=n->neighbor(ind); 
				ind1= this->mirror_index(n,ind);
			}
			n2=(*next).first;
			ind2=(*next).second;
			// in case n2 is no longer a triangle of the new triangulation
			if (n2->neighbor(ind2) != Face_handle() ) {
				n=n2->neighbor(ind2); 
				// ind=this->mirror_index(n2,ind2);
				// mirror_index does not work in this case
				ind = cw(n->index(n2->vertex(cw(ind2))));
				n2=n->neighbor(ind); 
				ind2= this->mirror_index(n,ind);
			}
			
			Vertex_handle v0=n1->vertex(ccw(ind1));
			Vertex_handle v1=n1->vertex(cw(ind1));
			Vertex_handle v2=n2->vertex(cw(ind2));
			orient = this->orientation(v0->point(),v1->point(),v2->point());
			switch (orient) {
				case RIGHT_TURN : 	  		
					// creates the new triangle v0v1v2
					// updates the neighbors, the constraints 
					//and the list of new edges
					newlf = create_face(v0,v2,v1);
					new_edges.push_back(Edge(newlf,2));
					newlf->set_neighbor(1, n1);
					newlf->set_neighbor(0, n2);
					n1->set_neighbor(ind1, newlf);
					n2->set_neighbor(ind2, newlf);
					if (n1->is_constrained(ind1)) {
						newlf->set_constraint(1,true);
					}
					if (n2->is_constrained(ind2)) {
						newlf->set_constraint(0,true);
					}
					// v0, v1 or v2.face() may have been removed
					v0->set_face(newlf); 
					v1->set_face(newlf);
					v2->set_face(newlf);
					// update list_edges
					tempo=current;
					current=list_edges.insert(current, Edge(newlf,2));
					list_edges.erase(tempo);
					list_edges.erase(next);
					next=current;
					if (v0 != va) {--current;} 
					else {++next;} 
					break;
				case LEFT_TURN : 	  
					++current; ++next;
					break;
				case COLLINEAR : 
					++current; ++next;
					break;
			}
		} while (next != list_edges.end());
	}
	
			
//---------------CONSTRAINED  Setting & updating-------------------------------------------------
		
	template < class Gt, class Tds >
	inline void
	Constrained_Delaunay_triangulation_sphere_2<Gt,Tds>::
	mark_constraint(Face_handle fr, int i)
	{
		if (dimension()<=1) fr->set_constraint(2, true);
		else{
			fr->set_constraint(i,true);
			fr->neighbor(i)->set_constraint(this->mirror_index(fr,i),true);
		}
		return;
	}
	
	template < class Gt, class Tds >
	void
	Constrained_Delaunay_triangulation_sphere_2<Gt,Tds>::  
	update_constraints_opposite(Vertex_handle va)
	// update status of edges opposite to a
	// after insertion of a
	{
		CGAL_triangulation_assertion(dimension()==2); 
		Face_handle f=va->face(), start=f;
		int indf;
		do {
			indf = f->index(va);
			if (f->neighbor(indf)->is_constrained(this->mirror_index(f,indf)) ) {
				f->set_constraint(indf,true);
			}
			else {
				f->set_constraint(indf,false);
			}
			f= f->neighbor(ccw(indf)); // turns ccw around va 
		} while (f != start);
		return;
	}
	
	
	template < class Gt, class Tds >
	void
	Constrained_Delaunay_triangulation_sphere_2<Gt, Tds>::
	update_constraints_incident(Vertex_handle va, 
								Vertex_handle c1,
								Vertex_handle c2)
	// update status of edges incident to a 
	// after insertion in the  constrained edge c1c2
	{
		if (dimension() == 0) return;
		if (dimension()== 1) {
			Edge_circulator ec=this->incident_edges(va), done(ec);
			do {
				((*ec).first)->set_constraint(2,true);
			}while (++ec != done);
		}
		else{
			//dimension() ==2
			int cwi, ccwi, indf;
			Face_circulator fc=this->incident_faces(va), done(fc);  
			CGAL_triangulation_assertion(fc != 0);
			do {
				indf = fc->index(va);
				cwi=cw(indf);
				ccwi=ccw(indf); 
				if ((fc->vertex(cwi) == c1)||(fc->vertex(cwi) == c2)) {
					fc->set_constraint(ccwi,true);
					fc->set_constraint(cwi,false);
				}	
				else {
					fc->set_constraint(ccwi,false);
					fc->set_constraint(cwi,true);
				}
				++fc;
			} while (fc != done);
		}
	}
	

	
	template < class Gt, class Tds  >
	void
	Constrained_Delaunay_triangulation_sphere_2<Gt,Tds>::
	clear_constraints_incident(Vertex_handle va)
	// make the edges incident to a newly created vertex unconstrained
	{
		Edge_circulator ec=this->incident_edges(va), done(ec);
		Face_handle f;
		int indf;
		if ( ec != 0){
			do {
				f = (*ec).first ;
				indf = (*ec).second;
				f->set_constraint(indf,false);
				if (dimension() == 2) {
					f->neighbor(indf)->set_constraint(this->mirror_index(f,indf),false);
				}
			} while (++ec != done);
		}
		return;
	}
	
	
	
	//-------------FLIPS------------------
	
									 
	template < class Gt, class Tds>
	void 
	Constrained_Delaunay_triangulation_sphere_2<Gt,Tds>::
	propagating_flip(Face_handle& f,int i)
	{ 
		if (!is_flipable(f,i)) return;
		Face_handle ni = f->neighbor(i); 
		flip(f, i); 
		propagating_flip(f,i); 
		i = ni->index(f->vertex(i)); 
		propagating_flip(ni,i); 
	} 
									 
									 
									 	
	template < class Gt, class Tds > 
	void  
	Constrained_Delaunay_triangulation_sphere_2<Gt,Tds>:: 
	propagating_flip(List_edges & edges) 
    {
		 propagating_flip(edges,Emptyset_iterator());
	}
	


template < class Gt, class Tds >
void 
Constrained_Delaunay_triangulation_sphere_2<Gt,Tds>::
flip (Face_handle& f, int i)
{
	Face_handle g = f->neighbor(i);
	int j = this->mirror_index(f,i);
	
	// save wings neighbors to be able to restore contraint status
	Face_handle f1 = f->neighbor(cw(i));
	int i1 = this->mirror_index(f,cw(i));
	Face_handle f2 = f->neighbor(ccw(i));
	int i2 = this->mirror_index(f,ccw(i));
	Face_handle f3 = g->neighbor(cw(j));
	int i3 = this->mirror_index(g,cw(j));
	Face_handle f4 = g->neighbor(ccw(j));
	int i4 = this->mirror_index(g,ccw(j));
	
	// The following precondition prevents the test suit 
	// of triangulation to work on constrained Delaunay triangulation
	//CGAL_triangulation_precondition(is_flipable(f,i));
	this->_tds.flip(f, i);
	
	// restore constraint status
	f->set_constraint(f->index(g), false);
	g->set_constraint(g->index(f), false);
	f1->neighbor(i1)->set_constraint(this->mirror_index(f1,i1),
									 f1->is_constrained(i1));
	f2->neighbor(i2)->set_constraint(this->mirror_index(f2,i2),
									 f2->is_constrained(i2));
	f3->neighbor(i3)->set_constraint(this->mirror_index(f3,i3),
									 f3->is_constrained(i3));
	f4->neighbor(i4)->set_constraint(this->mirror_index(f4,i4),
									 f4->is_constrained(i4));
	return;
}

	//---------Construction-------------
	template<class Gt, class Tds>
	inline
	typename Constrained_Delaunay_triangulation_sphere_2<Gt,Tds>::Point
	Constrained_Delaunay_triangulation_sphere_2<Gt,Tds>::
	circumcenter (const Point& p0, const Point& p1, const Point&  p2) const
	{
		return 
		DTS::geom_traits().construct_circumcenter_2_object()(p0,p1,p2);
	}
	
	
	template <class Gt, class Tds >
	typename Constrained_Delaunay_triangulation_sphere_2<Gt, Tds>::Point
	Constrained_Delaunay_triangulation_sphere_2<Gt, Tds>::
	circumcenter(Face_handle  f) const
	{
		return circumcenter((f->vertex(0))->point(), 
							(f->vertex(1))->point(), 
							(f->vertex(2))->point());
	}
	
	
	//--------
	
	
	template <class Gt, class Tds >
	typename Constrained_Delaunay_triangulation_sphere_2<Gt, Tds>::Face_handle
	Constrained_Delaunay_triangulation_sphere_2<Gt,Tds>::
	locate(const Point& p,Locate_type& lt,int& li, Face_handle start) const{
		Face_handle fh= DTS::locate(p,lt,li,start);
		if( lt == TOO_CLOSE)
		/* for being able to avoid inserting a point which is too close to an existing
		 vertex during the meshing, location_typ has to be changed to VERTEX*/
			lt = VERTEX;
		return fh;
		
	}	
	
 
 }//end namespace
#endif

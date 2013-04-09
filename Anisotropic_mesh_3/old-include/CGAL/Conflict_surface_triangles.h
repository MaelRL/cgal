#ifndef CONFLICT_SURFACE_TRIANGLES_H
#define CONFLICT_SURFACE_TRIANGLES_H

#include <list>		  //Coherenct_triangle_list
#include <queue>	//Conflict_triangle_list
#include <vector>
#include <iostream>
#include <algorithm>  //sort

//==========================================================================
//                 Structure: Triangle
//==========================================================================


/**General storage of triangle with 3 index 
 * sorted to avoid duplication in sorted list.
 */
	struct Triangle {
		Triangle(unsigned _i, unsigned _j, unsigned _k)
		{
			index[0] = _i;
			index[1] = _j;
			index[2] = _k;
			std::sort(index,index+3);
		}

		bool operator< (const Triangle& right) const {
			for(int i=0;i<2;i++){
				if(index[i] < right.index[i])
					return true;
				if(index[i] > right.index[i])
					return false;
			}
			return index[2] < right.index[2];
		}

		bool operator== (const Triangle& right) const{
			for(int i=0;i<3;i++){
				if(index[i]!=right.index[i])
					return false;
			}
			return true;
		}

		unsigned index[3];
	};

//==========================================================================
//                 Structure: Surface_list store all with information 
//                 about coherency, facilitating visualization
//==========================================================================


/** Storage of Surface triangles with information about coherency
 **/
struct Surface_list {

	class Coherent_triangle_handle:
		public Triangle
	{
		public:

			Coherent_triangle_handle(unsigned i,unsigned j,unsigned k,bool c=false):
				Triangle(i,j,k),
				coherent(c){}
			Coherent_triangle_handle(const Coherent_triangle_handle &tc):
				Triangle(tc.index(0),tc.index(1),tc.index(2)),
				coherent(tc.is_coherent()){}

			void set_coherent(bool c=true){
				coherent = c;
			};

			bool is_coherent()const {
				return coherent;
			};

			unsigned index(int i) const {
				i=i%3;
				return Triangle::index[i];
			}

		private:

			bool coherent;
	};

	class Less_tr{
		public:
			bool operator() (const Coherent_triangle_handle &a, const Coherent_triangle_handle &b)const {
				return a < b;
			}
	};

	typedef std::vector<Coherent_triangle_handle>		T_list;
	typedef T_list::iterator						iterator;



	void push_back(unsigned i, unsigned j, unsigned k,bool c=false) {
		t_list.push_back(Coherent_triangle_handle(i, j, k,c));
	}

	iterator begin(){
		return t_list.begin();
	}

	iterator end(){
		return t_list.end();
	}

	void clear(){
		t_list.clear();
	}

	//does not make unique?
	void sort(){
		std::sort(begin(),end(),Less_tr());
	}

	//TODO: make it faster by binary search
	bool exist(unsigned i,unsigned j,unsigned k){
		Coherent_triangle_handle cth(i,j,k,true);
		return std::find(begin(),end(),cth) != end();
	}

	unsigned size() {
		return t_list.size();
	}

	private:
	T_list			t_list;
};//Surface_list


//====================================================================================================
///                Structure: Triangle_list store all conflict triangles
///                and sort the conflicts by the worst to the lest bad
///                pop out an element, the user will CHECK IF THE ELEMENT IS STILL VALID 
///                in the triangulation.
///
///	The first parameter must provide methods of :
///
///	FT	score(Point_3 a,Point_3 b,Point_3 c);
///		(in which the worst scores the highest)
///
/// The second parameter must provide the geometrical traits
/// 
//====================================================================================================


//===========================================================================================
///                Structure: Triangle_list store all conflict triangles
///                and sort the conflicts by the worst to the lest bad
///
///                There are maybe "repeated" triangles, i.e. the same 3 indice
///                with different "star center", and the element in the queue will never
///                be deleted once they are inserted. So we suppose that each time 
///                pop out an element, the user will CHECK IF THE ELEMENT IS STILL VALID 
///                in the triangulation.
///
///	The first parameter must provide methods of :
///  compute_distorsion, compute_squared_distance, transform_metric, construct_circumcenter
///
/// The second parameter Point_container is the container of the points to whom the indice of 
/// triangles gives the access. It must provide direct access of points by index with operator[].
/// Its elements are expected to have the type Surface_criteria::Point_3
///
//=============================================================================================
template<typename Umbrella_set,typename Constrain_surface>
class Conflict_triangle_list {

	public:

	class Conflict_triangle_handle
	{
		typedef typename Umbrella_set::Surface_criteria			Surface_criteria;
		typedef typename Umbrella_set::Traits					Traits;
		typedef typename Traits::FT								FT;
		typedef typename Traits::Point_3						Point_3;
		public:
		Conflict_triangle_handle(const Traits *_traits,
				unsigned i,unsigned j,unsigned center):
			index_center(center),
			index_a(i),
			index_b(j),
			traits(_traits),
			score_cache(-1.){}

		Conflict_triangle_handle(const Conflict_triangle_handle &tc):
			index_center(tc.get_center()),
			index_a(tc.index(0)),
			index_b(tc.index(1)),
			traits(tc.get_traits()),
			score_cache(tc._raw_score()){}

		Conflict_triangle_handle &operator= (const Conflict_triangle_handle &tc) {
			index_center = tc.get_center();
			index_a = tc.index(0);
			index_b = tc.index(1);
			traits = tc.get_traits();
			score_cache = tc._raw_score();
			return *this;
		}

		private:
		FT _raw_score()const {
			return score_cache;
		}
		public:

		~Conflict_triangle_handle(){
		};

		void set_center(unsigned c){
			index_center = c;
		};

		bool is_bad(const Umbrella_set	&umbs){
			return umbs.get_surface_criteria()->is_bad(umbs.get_point(index_a),
					umbs.get_point(index_b),umbs.get_point(index_center),*traits);
		}

		FT	score(const Umbrella_set &umbs){
			if(score_cache<0.){
				score_cache = umbs.get_surface_criteria()->score(umbs.get_point(index_a),
						umbs.get_point(index_b),umbs.get_point(index_center),*traits);
			}
			return score_cache;
		}

		const Traits  *get_traits() const {
			return traits;
		}

		unsigned get_center()const {
			return index_center;
		};

		unsigned index(int i) const {
			i=i%2;
			if(i!=0){
				return index_b;
			}else{
				return index_a;
			}
		}


		private:

		unsigned index_center;	//The index of the vertex whose metric is the center.
		unsigned index_a,index_b;  //The other 2 indices of the triangle vertex
		const Surface_criteria  *surface_criteria;
		const Traits  *traits;
		FT		  score_cache;
	};//Conflict_triangle_handle

	class Less_conflicting{
		public:
		Less_conflicting(const Umbrella_set &_umbs):umbs(_umbs){};
		bool operator() (Conflict_triangle_handle &left,
				Conflict_triangle_handle &right) const {
			return left.score(umbs)<right.score(umbs);
		}
		private:
		const Umbrella_set &umbs;
	};


	public:
	typedef typename Umbrella_set::Traits							Traits;
		typedef typename Umbrella_set::Surface_criteria					Surface_criteria;
		typedef typename Traits::FT										FT;
		typedef typename Traits::Point_3								Point_3;
		typedef typename std::priority_queue<Conflict_triangle_handle,
				std::vector<Conflict_triangle_handle>,
				Less_conflicting>										T_queue;

		//========================================================================
		//                    Conflict_triangle_list functions
		//========================================================================

		Conflict_triangle_list(const Umbrella_set& _umbs):
			umbs(_umbs),
			t_queue(Less_conflicting(umbs))
	{}

		/// \return if the push action is successful or not
		/// the element is inserted only if its has a higher score than 0,
		/// i.e. that the element "is_bad()"
		bool push_if_bad(const Traits &traits, 
				unsigned i,unsigned j,unsigned center,
				bool coherent){
			Conflict_triangle_handle cth(&traits,i,j,center);

			//test if it is not coherent or bad triangle.
			if((!coherent) || cth.is_bad(umbs)){	
				t_queue.push(cth);
				return true;
			}
			return false;
		}

		bool empty(){
			return t_queue.empty();
		}

		unsigned size(){
			return t_queue.size();
		}

		Conflict_triangle_handle top(){
			return t_queue.top();
		}

		void pop(){
			t_queue.pop();
		}

		//========================================================================
		//                    Conflict_triangle_list members
		//========================================================================

	private:
		const Umbrella_set &	umbs;
		T_queue					t_queue;
};
#endif//CONFLICT_SURFACE_TRIANGLES_H

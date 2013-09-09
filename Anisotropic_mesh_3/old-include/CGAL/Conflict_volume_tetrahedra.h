#ifndef CONFLICT_VOLUME_TETRAHEDRA
#define CONFLICT_VOLUME_TETRAHEDRA

#include <list>		  //Coherenct_tetrahedra_list
#include <queue>	//Conflict_tetrahedra_list
#include <vector>
#include <iostream>
#include <algorithm>  //sort

#include <CGAL/Surface_criteria.h>	  //extended incl

//==========================================================================
//                 Structure: Tetrahedra
//==========================================================================

/**General storage of tetrahedra with 4 index 
 * sorted to avoid duplication in sorted list.
 */
struct Tetrahedra_indice {
	Tetrahedra_indice (unsigned _i, unsigned _j, unsigned _k,unsigned _l)
	{
		ind[0] = _i;
		ind[1] = _j;
		ind[2] = _k;
		ind[3] = _l;
		std::sort(ind,ind+4);
	}

	Tetrahedra_indice (const Tetrahedra_indice & ti)
	{
		for(int i=0;i<4;i++){
			ind[i] = ti.ind[i];
		}
	}

	bool operator< (const Tetrahedra_indice & right) const {
		for(int i=0;i<3;i++){
			if(ind[i] < right.ind[i])
				return true;
			if(ind[i] > right.ind[i])
				return false;
		}
		return ind[3] < right.ind[3];
	}

	bool operator== (const Tetrahedra_indice & right) const{
		for(int i=0;i<4;i++){
			if(ind[i]!=right.ind[i])
				return false;
		}
		return true;
	}

	unsigned index(unsigned i){
		i = i%4;
		return ind[i];
	}

  private:
	unsigned ind[4];
};

//==========================================================================
//                 Structure: Tetrahedra_list store only the coherent 
//					  tetrahedra, facilitating visualization
//==========================================================================


/** Storage of tetrahedras with information about coherency
 **/
struct Tetrahedra_list {

	class Less_tetra{
		public:
		bool operator() (const Tetrahedra_indice &a, const Tetrahedra_indice &b)const {
			return a < b;
		}
	};

	class Equal_tetra{
		public:
			bool operator() (const Tetrahedra_indice &a, const Tetrahedra_indice &b)const {
				return a == b;
			}
	};

	typedef std::list<Tetrahedra_indice>			T_list;
	typedef T_list::iterator						Tetrahedra_indice_iterator;

	void push_back(unsigned i, unsigned j, unsigned k,unsigned l) {
		t_list.push_back(Tetrahedra_indice(i, j, k,l));
	}

	Tetrahedra_indice_iterator begin(){
		return t_list.begin();
	}

	Tetrahedra_indice_iterator end(){
		return t_list.end();
	}

	void clear(){
		t_list.clear();
	}

  Tetrahedra_indice pop_back() {
	Tetrahedra_indice t(t_list.back());
	t_list.pop_back();
	return t;
  }

  void make_unique(){
	  t_list.sort(Less_tetra());
	  t_list.unique(Equal_tetra());
  }

  unsigned size() {
    return t_list.size();
  }

  private:
  T_list			t_list;
};//Tetrahedra_list


//====================================================================================================
///                Volume: Conflict_tetrahedra_list store all conflict tetrahedras
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
///                Structure: Triangle_list store all conflict tetrahedras
///                and sort the conflicts by the worst to the lest bad
///
///                There are maybe "repeated" tetrahedras, i.e. the same 4 indice
///                with different "star center", and the element in the queue will never
///                be deleted once they are inserted. So we suppose that each time 
///                pop out an element, the user will CHECK IF THE ELEMENT IS STILL VALID 
///                in the triangulation.
///
///	The first parameter must provide methods of :
///  compute_distorsion, compute_squared_distance, transform_metric, construct_circumcenter
///
/// The second parameter Point_container is the container of the points to whom the indice of 
/// tetrahedras gives the access. It must provide direct access of points by index with operator[].
/// Its elements are expected to have the type Surface_criteria::Point_3
///
//=============================================================================================
template<typename Umbrella_set>
class Conflict_tetrahedra_list {
	
  public:
		  class Conflict_tetrahedra_handle
		  {
			  typedef typename Umbrella_set::Traits					Traits;
			  typedef typename Umbrella_set::Volume_criteria			Volume_criteria;
			  typedef typename Umbrella_set::Cell_handle				Cell_handle;
			  typedef typename Umbrella_set::FT						FT;
			  typedef typename Traits::Point_3						Point_3;
			  public:
			  Conflict_tetrahedra_handle(unsigned i,unsigned j,unsigned k,unsigned center):
				  score_cache_valid(false),
				  index_center(center){
					  index_other[0]=i;
					  index_other[1]=j;
					  index_other[2]=k;
				  }

			  Conflict_tetrahedra_handle(const Conflict_tetrahedra_handle &tc):
				  score_cache(tc._get_score_cache()),
				  score_cache_valid(tc._get_score_valid()),
				  index_center(tc.get_center()){
					  for(int i=0;i<3;i++){
						  index_other[i]=tc.index(i);
					  }
				  }

			  Conflict_tetrahedra_handle &operator= (const Conflict_tetrahedra_handle &tc) {
				  index_center = tc.get_center();
				  score_cache_valid = tc._get_score_valid();
				  score_cache = tc._get_score_cache();
				  for(int i=0;i<3;i++){
					  index_other[i]=tc.index(i);
				  }
				  return *this;
			  }

			  bool is_bad(const Umbrella_set &umbrella_set) const {
				  const Volume_criteria *volume_criteria = umbrella_set.get_volume_criteria();
				  const Cell_handle cell = umbrella_set.get_cell_handle(
						  index_center,index_other[0],index_other[1],index_other[2]);
				  const Traits &traits = umbrella_set.get_traits(index_center);
				  return volume_criteria->is_bad(cell,traits);
			  }

			  FT	score(const Umbrella_set &umbrella_set){
				  if(!score_cache_valid){
					  const Volume_criteria *volume_criteria = umbrella_set.get_volume_criteria();
					  const Cell_handle cell = umbrella_set.get_cell_handle(
							  index_center,index_other[0],index_other[1],index_other[2]);
					  const Traits &traits = umbrella_set.get_traits(index_center);
					  score_cache = volume_criteria->score(cell,traits);
					  score_cache_valid = true;
				  }
				  return score_cache;
			  }

			  unsigned index(unsigned i)const{
				  i %= 3;
				  return index_other[i];
			  }

			  unsigned get_center()const {
				  return index_center;
			  };

			  private:

			  bool _get_score_valid() const{
				  return score_cache_valid;
			  }

			  FT _get_score_cache() const {
				  return score_cache;
			  }

			  FT score_cache;
			  bool score_cache_valid;
			  unsigned index_other[3];			//The other 3 indices except the center
			  unsigned index_center;				//The index of the vertex whose metric is the center.
		  };

		  class Tetrahedra_less_conflicting{
			  const Umbrella_set	&umbs;
			  public:
			  Tetrahedra_less_conflicting(const Umbrella_set &_umbs)
				  :umbs(_umbs){};

			  bool operator() (Conflict_tetrahedra_handle &left,
					  Conflict_tetrahedra_handle &right) const {
				  return left.score(umbs)<right.score(umbs);
			  }
		  };

	  typedef typename Umbrella_set::Traits					Traits;
	  typedef typename Umbrella_set::Volume_criteria			Volume_criteria;
	  typedef typename Umbrella_set::Cell_handle				Cell_handle;
	  typedef typename Traits::FT								FT;
	  typedef typename Traits::Point_3						Point_3;
	  typedef typename std::priority_queue<Conflict_tetrahedra_handle,
			  std::vector<Conflict_tetrahedra_handle>,
			  Tetrahedra_less_conflicting>					T_queue;

	  //========================================================================
	  //                    Conflict_tetrahedra_list functions
	  //========================================================================

	  Conflict_tetrahedra_list(const Umbrella_set &_umbs):
		  umbs(_umbs),
		  t_queue(Tetrahedra_less_conflicting(umbs))
	{}

	  /// \return if the push action is successful or not
	  /// the element is inserted only if its has a higher score than 0,
	  /// i.e. that the element "is_bad()"
	  bool push_if_bad(unsigned i,unsigned j,unsigned k,
			unsigned center, bool coherent){
		Conflict_tetrahedra_handle cth(i,j,k,center);

		//test if it is not coherent or bad tetrahedra.
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

	Conflict_tetrahedra_handle top(){
		return t_queue.top();
	}

	void pop(){
		t_queue.pop();
	}

	//========================================================================
	//                    Conflict_tetrahedra_list members
	//========================================================================

	private:
	const Umbrella_set &umbs;
	T_queue				  t_queue;
};

#endif			//CONFLICT_VOLUME_TETRAHEDRA

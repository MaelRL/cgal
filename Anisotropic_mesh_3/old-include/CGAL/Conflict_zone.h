#ifndef CONFLICT_ZONE_H
#define CONFLICT_ZONE_H

#include <vector>

/** The conflict_zone for the whole umbrella_set 
 * is the union of all conflict_zones of each umbrella
 * where the conflict_zones are not empty.
 *
 * For update incoherence, it is enough to update the conflict zone
 *
 **/
template<typename Umbrella_set>
class Umbrella_set_mesher_level_traits_3
{

	public:
	  typedef typename Umbrella_set::Cell_handle	Cell_handle;
	/** Zone for the surface-mesher should be the union of the conflicting 
 */
	class Zone {
		typedef std::pair<unsigned,Cell_handle>			  Zone_elem;
		typedef std::vector<Zone_elem>					  Zone_list;
		public:

		Zone() {
			zone_list.reserve(32);
		}

		void push(unsigned index,Cell_handle cell_loc){
			zone_list.push_back(Zone_elem(index,cell_loc));
		}

		bool empty() const {
			return zone_list.empty();
		}

		bool exist(unsigned u) {
			for(typename Zone_list::iterator zit = zone_list.begin();
					zit != zone_list.end();
					zit ++){
				if(zit->first == u){
					return true;
				}
			}
			return false;
		}
			

		unsigned umb_index(unsigned index) const {
		  return zone_list[index].first;
		}

		Cell_handle umb_loc(unsigned index) const {
		  return zone_list[index].second;
		}

		unsigned size()  const {
			return zone_list.size();
		}

		void clear(){
			zone_list.clear();
		}

		private:
		Zone_list		  zone_list;
	};
};


#endif //CONFLICT_ZONE_H

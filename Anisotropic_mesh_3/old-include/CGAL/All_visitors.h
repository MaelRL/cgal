#ifndef ALL_VISITORS_H
#define ALL_VISITORS_H

#include <CGAL/Mesher_level_visitors.h>

template<typename Umbrella_set>
class Surface_conflicts_visitor
{
	public:

		typedef CGAL::Null_mesh_visitor Previous_visitor;

		Surface_conflicts_visitor(Umbrella_set &_umbs)
			:umbs(_umbs){};

		const Previous_visitor& previous_level() const{ 
			return pv;
		}

		template <typename E, typename P>
			void before_conflicts(E, P) const {}

		template <typename E, typename P, typename Z>
			void before_insertion(E, P, Z) const {}

		///pop out inexistant conflicts and insert new conflicts
		template <typename V>
			void after_insertion(V) const {
				umbs.pop_invalid_volume_conflicts();
				umbs.append_new_volume_conflicts_from_indices_cache();
			}

		template <typename E, typename P, typename Z>
			void after_no_insertion(E, P, Z) const {}

	private:
		Previous_visitor	pv;
		Umbrella_set		&umbs;
};

template<typename Umbrella_set>
class Volume_conflicts_visitor
{
	public:
		typedef Surface_conflicts_visitor<Umbrella_set> Previous_visitor;

		Volume_conflicts_visitor(Previous_visitor &_scv)
			:surface_conflicts_visitor(_scv)
		{}

		const Previous_visitor& previous_level() const{ 
			return surface_conflicts_visitor;
		}

		template <typename E, typename P>
			void before_conflicts(E, P) const {}

		template <typename E, typename P, typename Z>
			void before_insertion(E, P, Z) const {}

		///pop out inexistant conflicts and insert new conflicts
		template <typename V>
			void after_insertion(V) const {}

		template <typename E, typename P, typename Z>
			void after_no_insertion(E, P, Z) const {}

	private:
		const Previous_visitor &surface_conflicts_visitor;
};

#endif //ALL_VISITORS_H

#ifndef POLYHEDRON_TYPE_FWD_H
#define POLYHEDRON_TYPE_FWD_H

#include <CGAL/Enriched_items.h>

#if 1//def USE_FORWARD_DECL

#include <CGAL/Filtered_kernel_fwd.h>
#include <memory>

namespace CGAL
{
class Epick;

template<typename K>
class Robust_circumcenter_traits_3;

class Polyhedron_items_3;

#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
template < class T, class I, class A>
class HalfedgeDS_default;
#else
struct HalfedgeDS_default;
#endif

template < class PolyhedronTraits_3,
           class PolyhedronItems_3,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
           template < class T, class I, class A>
#endif
           class T_HDS,
           class Alloc
           >
class Polyhedron_3;

} // end namespace CGAL

// kernel
typedef CGAL::Robust_circumcenter_traits_3<CGAL::Epick> Kernel;

typedef CGAL::Polyhedron_3<Kernel,
                           CGAL::Anisotropic_mesh_3::Enriched_items,
                           CGAL::HalfedgeDS_default,
                           std::allocator<int> >                     Polyhedron;

#else
#include "Polyhedron_type.h"
#endif

#endif // POLYHEDRON_TYPE_FWD_H

#ifndef _TYPES_H
#define _TYPES_H

#include <CGAL/AABB_tree.h> // must be inserted before kernel
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>
#include <CGAL/Cartesian_matrix.h>

#include <CGAL/Simple_cartesian.h>
typedef CGAL::Simple_cartesian<double> Kernel; // fastest in experiments


typedef Kernel::FT FT;
typedef Kernel::Ray_3 Ray;
typedef Kernel::Line_3 Line;
typedef Kernel::Point_3 Point;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Segment_3 Segment;
typedef Kernel::Triangle_3 Triangle;

#include <CGAL/Linear_algebraCd.h>
typedef CGAL::Linear_algebraCd<FT> LA;
typedef LA::Matrix Matrix;
typedef LA::Vector Vector_Matrix;
 

// #include <CGAL/Polyhedron_3.h>
// #include <CGAL/Polyhedron_items_3.h>
// #include "include/Primitives.h"

// typedef CGAL::Polyhedron_items_3<Kernel, Enriched_items> Polyhedron;
// typedef Polyhedron::Vertex_handle Vertex_handle;
// typedef Polyhedron::Halfedge_handle Halfedge_handle;
// typedef Polyhedron::Facet_handle Facet_handle;
// 
// // iterators
// typedef Polyhedron::Facet_iterator Facet_iterator;
// typedef Polyhedron::Edge_iterator Edge_iterator;
// typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
// typedef Polyhedron::Vertex_iterator Vertex_iterator;
// 
// // circulators
// typedef Polyhedron::Halfedge_around_facet_circulator         HF_circulator;
// typedef Polyhedron::Halfedge_around_vertex_circulator        HV_circulator;


#endif // _TYPES_H

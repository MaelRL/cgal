// Copyright (c) 2011  INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Author(s) : Kan-Le Shi

#ifndef CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_SHELL_H
#define CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_SHELL_H

#include "../CGAL/Constrain_surface_3.h"

using namespace CGAL::Anisotropic_mesh_3;

#define PI 3.1415926535897932384626433832795
#define DOT(a,b)	((a).x() * (b).x() + (a).y() * (b).y() + (a).z() * (b).z())
#define LINEAR_BLEND(p0,p1,t)	Point_3((p0).x() * (1 - t) + (p1).x() * t, \
										(p0).y() * (1 - t) + (p1).y() * t, \
										(p0).z() * (1 - t) + (p1).z() * t)

template<typename K, typename Point_container = std::vector<typename K::Point_3> >
class Constrain_surface_3_shell : public Constrain_surface_3<K, Point_container> {
public:
	typedef typename K::Point_3		Point_3;
	typedef typename K::Object_3	Object_3;
	typedef typename K::Segment_3	Segment_3;
	typedef typename K::Ray_3		Ray_3;
	typedef typename K::FT			FT;
	typedef typename CGAL::Oriented_side	Oriented_side;

public:
	FT radius_0;
	FT radius_1;
	EdgeList edges;

protected:
	Object_3 intersection_of_ray(const Ray_3 &ray) const {	
		return intersection(ray.source(), ray.source() + 
			ray.to_vector() * ((radius_1 * 3.0) / sqrt(DOT(ray.to_vector(), ray.to_vector()))));
	}

	Object_3 intersection_of_segment(const Segment_3 &seg) const {
		return intersection(seg.source(), seg.target());
	}

public:
	Object_3 intersection(const Point_3 &p0, const Point_3 &p1) const {
		FT r0, r1;
		Point_3 x0, x1;
		{
			FT C = DOT(p0, p0) - radius_0 * radius_0;
			FT B = DOT(p0, (p1 - p0)) * 2.0;
			FT A = DOT(p1 - p0, p1 - p0);
			FT delta = B * B - 4.0 * A * C;
			if (delta < 0) {
				r0 = -1.0;
			} else {
				FT root_1 = (-B - sqrt(delta)) / 2.0 / A;
				FT root_2 = (-B + sqrt(delta)) / 2.0 / A;
				if ((root_1 >= 0) && (root_1 <= 1)) {
					x0 = LINEAR_BLEND(p0, p1, root_1);
					r0 = root_1;
				} else if ((root_2 >= 0) && (root_2 <= 1)) {
					x0 = LINEAR_BLEND(p0, p1, root_2);
					r0 = root_2;
				} else {
					r0 = -1.0;
				}
			}
		}
		{
			FT C = DOT(p0, p0) - radius_1 * radius_1;
			FT B = DOT(p0, (p1 - p0)) * 2.0;
			FT A = DOT(p1 - p0, p1 - p0);
			FT delta = B * B - 4.0 * A * C;
			if (delta < 0) {
				r1 = -1.0;
			} else {
				FT root_1 = (-B - sqrt(delta)) / 2.0 / A;
				FT root_2 = (-B + sqrt(delta)) / 2.0 / A;
				if ((root_1 >= 0) && (root_1 <= 1)) {
					x1 = LINEAR_BLEND(p0, p1, root_1);
					r1 = root_1;
				} else if ((root_2 >= 0) && (root_2 <= 1)) {
					x1 = LINEAR_BLEND(p0, p1, root_2);
					r1 = root_2;
				} else {
					r1 = -1.0;
				}
			}
		}
		if ((r0 < 0) && (r1 < 0))
			return Object_3();
		else if (r0 < 0)
			return make_object(x1);
		else if (r1 < 0)
			return make_object(x0);
		else if (r0 < r1)
			return make_object(x0);
		else
			return make_object(x1);
	}

	Object_3 intersection(const Object_3 &obj) const {
		Segment_3 seg;
		Ray_3 ray;
		if (CGAL::assign(seg, obj))
			return intersection_of_segment(seg);
		else if (CGAL::assign(ray, obj))
			return intersection_of_ray(ray);
		else
			return Object_3();
	}

	FT get_bounding_radius() const {
		return radius_1 * 1.05;
	}

	Oriented_side side_of_constraint(const Point_3 &p) const {
		FT r = DOT(p, p);
		if ((r < radius_1 * radius_1 - 1e-6) && (r > radius_0 * radius_0 + 1e-6))
			return CGAL::ON_POSITIVE_SIDE;
		else if ((r > radius_1 * radius_1 + 1e-6) || (r < radius_0 * radius_0 - 1e-6))
			return CGAL::ON_NEGATIVE_SIDE;
		else
			return CGAL::ON_ORIENTED_BOUNDARY;
	}

	Point_container initial_points() const {
		Point_container points;
		const int split = 12;
		FT delta = 2.0 * PI / (FT)split;
		for (FT theta = -PI / 2.0 + delta / 2.0; theta < PI / 2.0; theta += delta)
			for (FT phi = 0; phi < 2.0 * PI - delta / 2.0; phi += delta)
			{
				points.push_back(Point_3(
					cos(phi + theta / 2.0) * cos(theta) * radius_0, 
					sin(phi + theta / 2.0) * cos(theta) * radius_0,
					sin(theta) * radius_0));
				points.push_back(Point_3(
					cos(phi + theta / 2.0) * cos(theta) * radius_1, 
					sin(phi + theta / 2.0) * cos(theta) * radius_1,
					sin(theta) * radius_1));
			}
		return points;
	}

	/*Point_3 project(const Point_3 &p) const {
		FT rate = radius / sqrt(p.x() * p.x() + p.y() * p.y() + p.z() * p.z());
		return Point_3(p.x() * rate, p.y() * rate, p.z() * rate);
	}*/

	EdgeIterator edge_begin() { return edges.begin(); }
	EdgeIterator edge_end() { return edges.end(); }
	void edge_split(EdgeIterator edge, const Point_3 & c) { }

	Constrain_surface_3_shell(const FT &radius_0_, const FT &radius_1_) : 
		radius_0(radius_0_), radius_1(radius_1_), edges() { }
	~Constrain_surface_3_shell() { };
};

#undef PI
#undef DOT
#undef LINEAR_BLEND

#endif

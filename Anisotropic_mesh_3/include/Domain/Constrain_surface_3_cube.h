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

#ifndef CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_CUBE_H
#define CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_CUBE_H

#include "../CGAL/Constrain_surface_3.h"

using namespace CGAL::Anisotropic_mesh_3;

#define PI 3.1415926535897932384626433832795
#define DOT(a,b)	((a).x() * (b).x() + (a).y() * (b).y() + (a).z() * (b).z())
#define LINEAR_BLEND(p0,p1,t)	Point_3((p0).x() * (1 - t) + (p1).x() * t, \
										(p0).y() * (1 - t) + (p1).y() * t, \
										(p0).z() * (1 - t) + (p1).z() * t)
#define SQUARE_DISTANCE(a,b)	(((a).x()-(b).x()) * ((a).x()-(b).x()) +  \
								 ((a).y()-(b).y()) * ((a).y()-(b).y()) +  \
								 ((a).z()-(b).z()) * ((a).z()-(b).z()))

template<typename K, typename Point_container = std::vector<typename K::Point_3> >
class Constrain_surface_3_cube : public Constrain_surface_3<K, Point_container> {
public:
	typedef typename K::Point_3		Point_3;
	typedef typename K::Object_3	Object_3;
	typedef typename K::Segment_3	Segment_3;
	typedef typename K::Ray_3		Ray_3;
	typedef typename K::FT			FT;
	typedef typename CGAL::Oriented_side	Oriented_side;

public:
	FT radius;
	EdgeList edges;

protected:
	Object_3 intersection_of_ray(const Ray_3 &ray) const {	
		return intersection(ray.source(), ray.source() + 
			ray.to_vector() * ((radius * 3.5) / sqrt(DOT(ray.to_vector(), ray.to_vector()))));
	}

	Object_3 intersection_of_segment(const Segment_3 &seg) const {
		return intersection(seg.source(), seg.target());
	}

	Object_3 intersection(const Point_3 &p0, const Point_3 &p1) const {
		FT min_t = 2.0;
		if (p1.x() - p0.x() != 0.0) {
			FT t = (-radius - p0.x()) / (p1.x() - p0.x());
			if ((t >= 0) && (t <= 1.0)) {
				Point_3 p = LINEAR_BLEND(p0, p1, t);
				if ((p.y() >= -radius) && (p.y() <= radius) &&
					(p.z() >= -radius) && (p.z() <= radius))
					if (t < min_t)
						min_t = t;
			}
			t = (radius - p0.x()) / (p1.x() - p0.x());
			if ((t >= 0) && (t <= 1.0)) {
				Point_3 p = LINEAR_BLEND(p0, p1, t);
				if ((p.y() >= -radius) && (p.y() <= radius) &&
					(p.z() >= -radius) && (p.z() <= radius))
					if (t < min_t)
						min_t = t;
			}
		}
		if (p1.y() - p0.y() != 0.0) {
			FT t = (-radius - p0.y()) / (p1.y() - p0.y());
			if ((t >= 0) && (t <= 1.0)) {
				Point_3 p = LINEAR_BLEND(p0, p1, t);
				if ((p.x() >= -radius) && (p.x() <= radius) &&
					(p.z() >= -radius) && (p.z() <= radius))
					if (t < min_t)
						min_t = t;
			}
			t = (radius - p0.y()) / (p1.y() - p0.y());
			if ((t >= 0) && (t <= 1.0)) {
				Point_3 p = LINEAR_BLEND(p0, p1, t);
				if ((p.x() >= -radius) && (p.x() <= radius) &&
					(p.z() >= -radius) && (p.z() <= radius))
					if (t < min_t)
						min_t = t;
			}
		}
		if (p1.z() - p0.z() != 0.0) {
			FT t = (-radius - p0.z()) / (p1.z() - p0.z());
			if ((t >= 0) && (t <= 1.0)) {
				Point_3 p = LINEAR_BLEND(p0, p1, t);
				if ((p.x() >= -radius) && (p.x() <= radius) &&
					(p.y() >= -radius) && (p.y() <= radius))
					if (t < min_t)
						min_t = t;
			}
			t = (radius - p0.z()) / (p1.z() - p0.z());
			if ((t >= 0) && (t <= 1.0)) {
				Point_3 p = LINEAR_BLEND(p0, p1, t);
				if ((p.x() >= -radius) && (p.x() <= radius) &&
					(p.y() >= -radius) && (p.y() <= radius))
					if (t < min_t)
						min_t = t;
			}
		}
		if (min_t > 1.0)
			return Object_3();
		else
			return make_object(LINEAR_BLEND(p0, p1, min_t));

		/*
		Point_3 pl = p0, pr = p1;
		Oriented_side vl = side_of_constraint(pl);
		Oriented_side vr = side_of_constraint(pr);
		if (vl == CGAL::ON_ORIENTED_BOUNDARY)
			return make_object(pl);
		if (vr == CGAL::ON_ORIENTED_BOUNDARY)
			return make_object(pr);
		FT dist = SQUARE_DISTANCE(pl, pr);
		while (true) {
			Point_3 pm = LINEAR_BLEND(pl, pr, 0.5);
			Oriented_side vm = side_of_constraint(pm);
			if (vm == CGAL::ON_ORIENTED_BOUNDARY)
				return make_object(project(pm));
			if (vm == vl) {
				pl = pm;
				vl = vm;
			} else {
				pr = pm;
				vr = vm;
			}
			dist /= 4.0;
			if (dist < 1e-10)
				return make_object(project(pm));
		}*/
	}

public:
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
		return radius * 2.0;
	}

	Oriented_side side_of_constraint(const Point_3 &p) const {
		if ((p.x() > radius) || (p.x() < -radius) ||
			(p.y() > radius) || (p.y() < -radius) ||
			(p.z() > radius) || (p.z() < -radius))
			return CGAL::ON_NEGATIVE_SIDE;
		else if ((p.x() > -radius) && (p.x() < radius) &&
			(p.y() > -radius) && (p.y() < radius) &&
			(p.z() > -radius) && (p.z() < radius))
			return CGAL::ON_POSITIVE_SIDE;
		else
			return CGAL::ON_ORIENTED_BOUNDARY;
	}

	Point_container initial_points() const {
		Point_container points;
		points.push_back(Point_3(-radius, -radius, -radius));
		points.push_back(Point_3(+radius, -radius, -radius));
		points.push_back(Point_3(-radius, +radius, -radius));
		points.push_back(Point_3(+radius, +radius, -radius));
		points.push_back(Point_3(-radius, -radius, +radius));
		points.push_back(Point_3(+radius, -radius, +radius));
		points.push_back(Point_3(-radius, +radius, +radius));
		points.push_back(Point_3(+radius, +radius, +radius));
		points.push_back(Point_3(0, 0, 0));
		
		points.push_back(Point_3(-radius, 0, 0));
		points.push_back(Point_3(+radius, 0, 0));
		points.push_back(Point_3(0, -radius, 0));
		points.push_back(Point_3(0, +radius, 0));
		points.push_back(Point_3(0, 0, -radius));
		points.push_back(Point_3(0, 0, +radius));

		points.push_back(Point_3(0, -radius, -radius));
		points.push_back(Point_3(0, -radius, +radius));
		points.push_back(Point_3(0, +radius, -radius));
		points.push_back(Point_3(0, +radius, +radius));
		points.push_back(Point_3(-radius, 0, -radius));
		points.push_back(Point_3(-radius, 0, +radius));
		points.push_back(Point_3(+radius, 0, -radius));
		points.push_back(Point_3(+radius, 0, +radius));
		points.push_back(Point_3(-radius, -radius, 0));
		points.push_back(Point_3(-radius, +radius, 0));
		points.push_back(Point_3(+radius, -radius, 0));
		points.push_back(Point_3(+radius, +radius, 0));

		return points;
	}

	Point_3 project(const Point_3 &p) const {
		Point_3 q = p;
		if (q.x() < -radius) q = Point_3(-radius, q.y(), q.z());
		if (q.y() < -radius) q = Point_3(q.x(), -radius, q.z());
		if (q.z() < -radius) q = Point_3(q.x(), q.y(), -radius);
		if (q.x() > radius) q = Point_3(radius, q.y(), q.z());
		if (q.y() > radius) q = Point_3(q.x(), radius, q.z());
		if (q.z() > radius) q = Point_3(q.x(), q.y(), radius);
		FT snap_dist = DBL_MAX, sd;
		int snap_type = 0;
		if ((sd = q.x() - (-radius)) < snap_dist) { snap_dist = sd; snap_type = 0; }
		if ((sd = q.y() - (-radius)) < snap_dist) { snap_dist = sd; snap_type = 1; }
		if ((sd = q.z() - (-radius)) < snap_dist) { snap_dist = sd; snap_type = 2; }
		if ((sd = radius - q.x()) < snap_dist) { snap_dist = sd; snap_type = 3; }
		if ((sd = radius - q.y()) < snap_dist) { snap_dist = sd; snap_type = 4; }
		if ((sd = radius - q.z()) < snap_dist) { snap_dist = sd; snap_type = 5; }
		switch (snap_type) {
		case 0:	q = Point_3(-radius, q.y(), q.z()); break;
		case 1:	q = Point_3(q.x(), -radius, q.z()); break;
		case 2:	q = Point_3(q.x(), q.y(), -radius); break;
		case 3:	q = Point_3(radius, q.y(), q.z()); break;
		case 4:	q = Point_3(q.x(), radius, q.z()); break;
		case 5:	q = Point_3(q.x(), q.y(), radius); break;
		}
		return q;
	}

	EdgeIterator edge_begin() {
		return edges.begin();
	}

	EdgeIterator edge_end() {
		return edges.end();
	}

	void edge_split(EdgeIterator edge, const Point_3 & c) {
		Edge e1 = Edge(edge->first, c);
		Edge e2 = Edge(c, edge->second);
		edges.erase(edge);
		edges.push_back(e1);
		edges.push_back(e2);
	}

	Constrain_surface_3_cube(const FT &radius_) : radius(radius_), edges() {
		Point_3 p000(-radius, -radius, -radius);
		Point_3 p100(+radius, -radius, -radius);
		Point_3 p010(-radius, +radius, -radius);
		Point_3 p110(+radius, +radius, -radius);
		Point_3 p001(-radius, -radius, +radius);
		Point_3 p101(+radius, -radius, +radius);
		Point_3 p011(-radius, +radius, +radius);
		Point_3 p111(+radius, +radius, +radius);
		edges.push_back(Edge(p000, p001));
		edges.push_back(Edge(p001, p011));
		edges.push_back(Edge(p011, p010));
		edges.push_back(Edge(p010, p000));
		edges.push_back(Edge(p100, p101));
		edges.push_back(Edge(p101, p111));
		edges.push_back(Edge(p111, p110));
		edges.push_back(Edge(p110, p100));
		edges.push_back(Edge(p000, p100));
		edges.push_back(Edge(p001, p101));
		edges.push_back(Edge(p010, p110));
		edges.push_back(Edge(p011, p111));
	}
	~Constrain_surface_3_cube() { };
};

template<typename K, typename Point_container = std::vector<typename K::Point_3> >
class Constrain_surface_3_flat_cube : public Constrain_surface_3<K, Point_container> {
public:
	typedef typename K::Point_3		Point_3;
	typedef typename K::Object_3	Object_3;
	typedef typename K::Segment_3	Segment_3;
	typedef typename K::Ray_3		Ray_3;
	typedef typename K::FT			FT;
	typedef typename CGAL::Oriented_side	Oriented_side;

public:
	FT radius;
	EdgeList edges;

protected:
	Object_3 intersection_of_ray(const Ray_3 &ray) const {	
		return intersection(ray.source(), ray.source() + 
			ray.to_vector() * ((radius * 3.5) / sqrt(DOT(ray.to_vector(), ray.to_vector()))));
	}

	Object_3 intersection_of_segment(const Segment_3 &seg) const {
		return intersection(seg.source(), seg.target());
	}

	Object_3 intersection(const Point_3 &p0, const Point_3 &p1) const {
		FT min_t = 2.0;
		if (p1.x() - p0.x() != 0.0) {
			FT t = (-radius - p0.x()) / (p1.x() - p0.x());
			if ((t >= 0) && (t <= 1.0)) {
				Point_3 p = LINEAR_BLEND(p0, p1, t);
				if ((p.y() >= -radius) && (p.y() <= radius) &&
					(p.z() >= -radius) && (p.z() <= radius))
					if (t < min_t)
						min_t = t;
			}
			t = (radius - p0.x()) / (p1.x() - p0.x());
			if ((t >= 0) && (t <= 1.0)) {
				Point_3 p = LINEAR_BLEND(p0, p1, t);
				if ((p.y() >= -radius) && (p.y() <= radius) &&
					(p.z() >= -radius) && (p.z() <= radius))
					if (t < min_t)
						min_t = t;
			}
		}
		if (p1.y() - p0.y() != 0.0) {
			FT t = (-radius - p0.y()) / (p1.y() - p0.y());
			if ((t >= 0) && (t <= 1.0)) {
				Point_3 p = LINEAR_BLEND(p0, p1, t);
				if ((p.x() >= -radius) && (p.x() <= radius) &&
					(p.z() >= -radius) && (p.z() <= radius))
					if (t < min_t)
						min_t = t;
			}
			t = (radius - p0.y()) / (p1.y() - p0.y());
			if ((t >= 0) && (t <= 1.0)) {
				Point_3 p = LINEAR_BLEND(p0, p1, t);
				if ((p.x() >= -radius) && (p.x() <= radius) &&
					(p.z() >= -radius) && (p.z() <= radius))
					if (t < min_t)
						min_t = t;
			}
		}
		if (p1.z() - p0.z() != 0.0) {
			FT t = (-radius - p0.z()) / (p1.z() - p0.z());
			if ((t >= 0) && (t <= 1.0)) {
				Point_3 p = LINEAR_BLEND(p0, p1, t);
				if ((p.x() >= -radius) && (p.x() <= radius) &&
					(p.y() >= -radius) && (p.y() <= radius))
					if (t < min_t)
						min_t = t;
			}
			t = (radius - p0.z()) / (p1.z() - p0.z());
			if ((t >= 0) && (t <= 1.0)) {
				Point_3 p = LINEAR_BLEND(p0, p1, t);
				if ((p.x() >= -radius) && (p.x() <= radius) &&
					(p.y() >= -radius) && (p.y() <= radius))
					if (t < min_t)
						min_t = t;
			}
		}
		if (min_t > 1.0)
			return Object_3();
		else
			return make_object(LINEAR_BLEND(p0, p1, min_t));

		/*
		Point_3 pl = p0, pr = p1;
		Oriented_side vl = side_of_constraint(pl);
		Oriented_side vr = side_of_constraint(pr);
		if (vl == CGAL::ON_ORIENTED_BOUNDARY)
			return make_object(pl);
		if (vr == CGAL::ON_ORIENTED_BOUNDARY)
			return make_object(pr);
		FT dist = SQUARE_DISTANCE(pl, pr);
		while (true) {
			Point_3 pm = LINEAR_BLEND(pl, pr, 0.5);
			Oriented_side vm = side_of_constraint(pm);
			if (vm == CGAL::ON_ORIENTED_BOUNDARY)
				return make_object(project(pm));
			if (vm == vl) {
				pl = pm;
				vl = vm;
			} else {
				pr = pm;
				vr = vm;
			}
			dist /= 4.0;
			if (dist < 1e-10)
				return make_object(project(pm));
		}*/
	}

public:
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
		return radius * 2.0;
	}

	Oriented_side side_of_constraint(const Point_3 &p) const {
		if ((p.x() > radius) || (p.x() < -radius) ||
			(p.y() > radius) || (p.y() < -radius) ||
			(p.z() > radius / 5.0) || (p.z() < -radius / 5.0))
			return CGAL::ON_NEGATIVE_SIDE;
		else if ((p.x() > -radius) && (p.x() < radius) &&
			(p.y() > -radius) && (p.y() < radius) &&
			(p.z() > -radius / 5.0) && (p.z() < radius / 5.0))
			return CGAL::ON_POSITIVE_SIDE;
		else
			return CGAL::ON_ORIENTED_BOUNDARY;
	}

	Point_container initial_points() const {
		Point_container points;
		points.push_back(Point_3(-radius, -radius, -radius / 5.0));
		points.push_back(Point_3(+radius, -radius, -radius / 5.0));
		points.push_back(Point_3(-radius, +radius, -radius / 5.0));
		points.push_back(Point_3(+radius, +radius, -radius / 5.0));
		points.push_back(Point_3(-radius, -radius, +radius / 5.0));
		points.push_back(Point_3(+radius, -radius, +radius / 5.0));
		points.push_back(Point_3(-radius, +radius, +radius / 5.0));
		points.push_back(Point_3(+radius, +radius, +radius / 5.0));
		points.push_back(Point_3(0, 0, 0));
		
		points.push_back(Point_3(-radius, 0, 0));
		points.push_back(Point_3(+radius, 0, 0));
		points.push_back(Point_3(0, -radius, 0));
		points.push_back(Point_3(0, +radius, 0));
		points.push_back(Point_3(0, 0, -radius / 5.0));
		points.push_back(Point_3(0, 0, +radius / 5.0));

		points.push_back(Point_3(0, -radius, -radius / 5.0));
		points.push_back(Point_3(0, -radius, +radius / 5.0));
		points.push_back(Point_3(0, +radius, -radius / 5.0));
		points.push_back(Point_3(0, +radius, +radius / 5.0));
		points.push_back(Point_3(-radius, 0, -radius / 5.0));
		points.push_back(Point_3(-radius, 0, +radius / 5.0));
		points.push_back(Point_3(+radius, 0, -radius / 5.0));
		points.push_back(Point_3(+radius, 0, +radius / 5.0));
		points.push_back(Point_3(-radius, -radius, 0));
		points.push_back(Point_3(-radius, +radius, 0));
		points.push_back(Point_3(+radius, -radius, 0));
		points.push_back(Point_3(+radius, +radius, 0));

		return points;
	}

	Point_3 project(const Point_3 &p) const {
		Point_3 q = p;
		if (q.x() < -radius) q = Point_3(-radius, q.y(), q.z());
		if (q.y() < -radius) q = Point_3(q.x(), -radius, q.z());
		if (q.z() < -radius / 5.0) q = Point_3(q.x(), q.y(), -radius / 5.0);
		if (q.x() > radius) q = Point_3(radius, q.y(), q.z());
		if (q.y() > radius) q = Point_3(q.x(), radius, q.z());
		if (q.z() > radius / 5.0) q = Point_3(q.x(), q.y(), radius / 5.0);
		FT snap_dist = DBL_MAX, sd;
		int snap_type = 0;
		if ((sd = q.x() - (-radius)) < snap_dist) { snap_dist = sd; snap_type = 0; }
		if ((sd = q.y() - (-radius)) < snap_dist) { snap_dist = sd; snap_type = 1; }
		if ((sd = q.z() - (-radius / 5.0)) < snap_dist) { snap_dist = sd; snap_type = 2; }
		if ((sd = radius - q.x()) < snap_dist) { snap_dist = sd; snap_type = 3; }
		if ((sd = radius - q.y()) < snap_dist) { snap_dist = sd; snap_type = 4; }
		if ((sd = radius / 5.0 - q.z()) < snap_dist) { snap_dist = sd; snap_type = 5; }
		switch (snap_type) {
		case 0:	q = Point_3(-radius, q.y(), q.z()); break;
		case 1:	q = Point_3(q.x(), -radius, q.z()); break;
		case 2:	q = Point_3(q.x(), q.y(), -radius / 5.0); break;
		case 3:	q = Point_3(radius, q.y(), q.z()); break;
		case 4:	q = Point_3(q.x(), radius, q.z()); break;
		case 5:	q = Point_3(q.x(), q.y(), radius / 5.0); break;
		}
		return q;
	}

	EdgeIterator edge_begin() {
		return edges.begin();
	}

	EdgeIterator edge_end() {
		return edges.end();
	}

	void edge_split(EdgeIterator edge, const Point_3 & c) {
		Edge e1 = Edge(edge->first, c);
		Edge e2 = Edge(c, edge->second);
		edges.erase(edge);
		edges.push_back(e1);
		edges.push_back(e2);
	}

	Constrain_surface_3_flat_cube(const FT &radius_) : radius(radius_), edges() {
		Point_3 p000(-radius, -radius, -radius / 5.0);
		Point_3 p100(+radius, -radius, -radius / 5.0);
		Point_3 p010(-radius, +radius, -radius / 5.0);
		Point_3 p110(+radius, +radius, -radius / 5.0);
		Point_3 p001(-radius, -radius, +radius / 5.0);
		Point_3 p101(+radius, -radius, +radius / 5.0);
		Point_3 p011(-radius, +radius, +radius / 5.0);
		Point_3 p111(+radius, +radius, +radius / 5.0);
		edges.push_back(Edge(p000, p001));
		edges.push_back(Edge(p001, p011));
		edges.push_back(Edge(p011, p010));
		edges.push_back(Edge(p010, p000));
		edges.push_back(Edge(p100, p101));
		edges.push_back(Edge(p101, p111));
		edges.push_back(Edge(p111, p110));
		edges.push_back(Edge(p110, p100));
		edges.push_back(Edge(p000, p100));
		edges.push_back(Edge(p001, p101));
		edges.push_back(Edge(p010, p110));
		edges.push_back(Edge(p011, p111));
	}
	~Constrain_surface_3_flat_cube() { };
};





template<typename K, typename Point_container = std::vector<typename K::Point_3> >
class Constrain_surface_3_free_cube : public Constrain_surface_3<K, Point_container> {
public:
	typedef typename K::Point_3		Point_3;
	typedef typename K::Object_3	Object_3;
	typedef typename K::Segment_3	Segment_3;
	typedef typename K::Ray_3		Ray_3;
	typedef typename K::FT			FT;
	typedef typename CGAL::Oriented_side	Oriented_side;

public:
	FT xmin, ymin, zmin;
	FT xmax, ymax, zmax;
	EdgeList edges;

protected:
	Object_3 intersection_of_ray(const Ray_3 &ray) const {	
		throw "not implemented";
	}

	Object_3 intersection_of_segment(const Segment_3 &seg) const {
		throw "not implemented";
	}

	Object_3 intersection(const Point_3 &p0, const Point_3 &p1) const {
		throw "not implemented";
	}

public:
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
		return xmax + ymax + zmax;
	}

	Oriented_side side_of_constraint(const Point_3 &p) const {
		if ((p.x() > xmax) || (p.x() < xmin) ||
			(p.y() > ymax) || (p.y() < ymin) ||
			(p.z() > zmax) || (p.z() < zmin))
			return CGAL::ON_NEGATIVE_SIDE;
		else if ((p.x() > xmin) && (p.x() < xmax) &&
			(p.y() > ymin) && (p.y() < ymax) &&
			(p.z() > zmin) && (p.z() < zmax))
			return CGAL::ON_POSITIVE_SIDE;
		else
			return CGAL::ON_ORIENTED_BOUNDARY;
	}

	Point_container initial_points() const {
		Point_container points;
		int split = 3;
		for (int i = 0; i <= split; i++)
			for (int j = 0; j <= split; j++)
				for (int k = 0; k <= split; k++)
					points.push_back(Point_3(
						(FT)i * xmax / (FT)split + xmin,
						(FT)j * ymax / (FT)split + ymin,
						(FT)k * zmax / (FT)split + zmin));
		return points;
	}

	Point_3 project(const Point_3 &p) const {
		throw "not implemented";
	}

	EdgeIterator edge_begin() {
		return edges.begin();
	}

	EdgeIterator edge_end() {
		return edges.end();
	}

	void edge_split(EdgeIterator edge, const Point_3 & c) {
		Edge e1 = Edge(edge->first, c);
		Edge e2 = Edge(c, edge->second);
		edges.erase(edge);
		edges.push_back(e1);
		edges.push_back(e2);
	}

	Constrain_surface_3_free_cube(
		const FT &xmin_, const FT &ymin_, const FT &zmin_,
		const FT &xmax_, const FT &ymax_, const FT &zmax_) :
		xmin(xmin_), ymin(ymin_), zmin(zmin_),	
		xmax(xmax_), ymax(ymax_), zmax(zmax_), edges() {
		Point_3 p000(xmin, ymin, zmin);
		Point_3 p100(xmax, ymin, zmin);
		Point_3 p010(xmin, ymax, zmin);
		Point_3 p110(xmax, ymax, zmin);
		Point_3 p001(xmin, ymin, zmax);
		Point_3 p101(xmax, ymin, zmax);
		Point_3 p011(xmin, ymax, zmax);
		Point_3 p111(xmax, ymax, zmax);
		edges.push_back(Edge(p000, p001));
		edges.push_back(Edge(p001, p011));
		edges.push_back(Edge(p011, p010));
		edges.push_back(Edge(p010, p000));
		edges.push_back(Edge(p100, p101));
		edges.push_back(Edge(p101, p111));
		edges.push_back(Edge(p111, p110));
		edges.push_back(Edge(p110, p100));
		edges.push_back(Edge(p000, p100));
		edges.push_back(Edge(p001, p101));
		edges.push_back(Edge(p010, p110));
		edges.push_back(Edge(p011, p111));
	}
	~Constrain_surface_3_free_cube() { };
};

#undef PI
#undef DOT
#undef LINEAR_BLEND
#undef SQUARE_DISTANCE

#endif

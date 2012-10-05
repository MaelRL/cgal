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

#ifndef CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_KLEIN_H
#define CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_KLEIN_H

#include "../CGAL/Constrain_surface_3_implicit.h"

using namespace CGAL::Anisotropic_mesh_3;

template<typename K, typename Point_container = std::vector<typename K::Point_3> >
class Constrain_surface_3_klein : public Constrain_surface_3_implicit<K, Point_container> {
public:
	typedef Constrain_surface_3_implicit<K, Point_container> Base;
	typedef typename Base::FT                                FT;

public:
	FT get_bounding_radius() const { return 3.0; }

	FT evaluate(const FT x, const FT y, const FT z) const {
		return (x*x+y*y+z*z+2*y-1)*((x * x+y*y+z*z -2*y-1)*(x * x+y*y+z*z -2*y-1)-
		 8*z*z)+16*x*z*(x*x+y*y+z*z-2*y-1);
	}

	Point_container initial_points() const {
		Point_container points;
		std::vector<Point_3> seeds;
		seeds.push_back(Point_3(0, 0, 0));
		return Base::initial_points(points, seeds, 0.2);
	}

	Constrain_surface_3_klein() { }
	~Constrain_surface_3_klein() { };
};

#endif

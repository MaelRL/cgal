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

#ifndef CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_CLEBSCH_H
#define CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_CLEBSCH_H

#include "../CGAL/Constrain_surface_3_implicit.h"

using namespace CGAL::Anisotropic_mesh_3;

template<typename K, typename Point_container = std::vector<typename K::Point_3> >
class Constrain_surface_3_clebsch : public Constrain_surface_3_implicit_with_bounding_sphere<K, Point_container> {
public:
	typedef Constrain_surface_3_implicit<K, Point_container> Base;
	typedef typename Base::FT                                FT;

public:
	FT get_bounding_radius() const { return 2.5; }

	FT evaluate(const FT x, const FT y, const FT z) const {
		return -(81*(x*x*x + y*y*y + z*z*z) - 189*(x*x*y + x*x*z + y*y*x + y*y*z + z*z*x + z*z*y) + 54*(x*y*z) + 126*(x*y + x*z + y*z) - 9*(x*x + y*y + z*z) - 9*(x + y + z) + 1);
	}

	Point_container initial_points() const {
		Point_container points;
		std::vector<Point_3> seeds;
		seeds.push_back(Point_3(0, 0, 0));
		return Base::initial_points(points, seeds, 0.2);
	}

	Constrain_surface_3_clebsch() { }
	~Constrain_surface_3_clebsch() { };
};

#endif

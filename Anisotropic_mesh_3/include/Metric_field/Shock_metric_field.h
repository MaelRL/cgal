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

#ifndef CGAL_ANISOTROPIC_MESH_3_SHOCK_METRIC_FIELD
#define CGAL_ANISOTROPIC_MESH_3_SHOCK_METRIC_FIELD

#include <CGAL/Metric_field.h>

using namespace CGAL::Anisotropic_mesh_3;

template<typename K>
class Shock_metric_field : public Metric_field<K> {
public:
	typedef Metric_field<K>        Base;
	typedef typename Base::FT      FT;
	typedef typename Base::Metric  Metric;
	typedef typename Base::Point_3 Point_3;

public:
	FT sigma;
	FT lambda;

public:
	virtual void report(typename std::ofstream &fx) const {
		fx << "type:   shock" << std::endl;
		fx << "sigma:  " << sigma << std::endl;
		fx << "lambda: " << lambda << std::endl;
	}
	virtual Metric compute_metric(const Point_3 &p) const {
		return Metric(
                    Vector_3(1, 0, 0), 
                    Vector_3(0, ((1.0 / lambda) + (1.0 - 1.0 / lambda) * (1.0 - exp(-(p.x() * p.x() / sigma / sigma)))),
			0), 
                    Vector_3(0, 0, 1));
	}
	Shock_metric_field(const FT &sigma_ = 0.5, const FT &lambda_ = 3.0) : sigma(sigma_), lambda(lambda_) { }
};


#endif

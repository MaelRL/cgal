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

#ifndef CGAL_ANISOTROPIC_MESH_3_ARCTAN_METRIC_FIELD
#define CGAL_ANISOTROPIC_MESH_3_ARCTAN_METRIC_FIELD

#include <CGAL/Metric_field.h>

namespace CGAL {
  namespace Anisotropic_mesh_3 {

template<typename K>
class Arctan_metric_field : public Metric_field<K> {
public:
	typedef Metric_field<K>         Base;
	typedef typename Base::FT       FT;
	typedef typename Base::Metric   Metric;
	typedef typename Base::Point_3  Point_3;
	typedef typename Base::Vector_3 Vector_3;

public:
	FT sigma;
	FT lambda;
  FT alpha;

public:
	virtual void report(typename std::ofstream &fx) const {
		fx << "type:   arctan" << std::endl;
		fx << "sigma:  " << sigma << std::endl;
		fx << "lambda: " << lambda << std::endl;
	}

	virtual Metric compute_metric(const Point_3 &p) const 
  {
    double epsilon = 1e-6;
    double d = ((atan(p.y() / sigma) / CGAL_PI) + 0.5) * (lambda - 1.0 / lambda) + (1.0 / lambda);
    return Metric(Vector_3(1, 0, 0), Vector_3(0, 1, 0), Vector_3(0, 0, 1), 
                  1., alpha*d, 1., epsilon);
		//return Metric(Vector_3(1, 0, 0), Vector_3(0, 
		//	((atan(p.y() / sigma) / CGAL_PI) + 0.5) * (lambda - 1.0 / lambda) + (1.0 / lambda), 0), 
		//	Vector_3(0, 0, 1));
	}
	Arctan_metric_field(FT lambda_, FT sigma_, FT alpha_) : 
          lambda(lambda_), 
          sigma(sigma_),
          alpha(alpha_){ }
};
  }
}


#endif

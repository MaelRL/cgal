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

#ifndef CGAL_ANISOTROPIC_MESH_3_STRETCHED_METRIC_FIELD
#define CGAL_ANISOTROPIC_MESH_3_STRETCHED_METRIC_FIELD

#include <CGAL/Metric_field.h>

using namespace CGAL::Anisotropic_mesh_3;

template<typename K>
class Stretched_metric_field : public Metric_field<K> {
public:
	typedef Metric_field<K>        Base;
	typedef typename Base::FT      FT;
	typedef typename Base::Metric  Metric;
	typedef typename Base::Point_3 Point_3;

public:
	FT stretch_parameter;
public:
	typedef typename Metric_field<K>::Metric	Metric;
	typedef typename Metric_field<K>::Point_3	Point_3;
public:
	virtual void report(typename std::ofstream &fx) const {
		fx << "type:  stretched" << std::endl;
		fx << "sigma: " << stretch_parameter << std::endl;
	}
	virtual Metric compute_metric(const Point_3 &p) const {
		return Metric(Vector_3(1, 0, 0), Vector_3(0, 1, 0), Vector_3(0, 0, exp(p.x() / stretch_parameter)));
	}
	Stretched_metric_field(const FT &stretch_parameter_) : stretch_parameter(stretch_parameter_) { }
};


/*
template<typename K>
class Hole_metric_field : public Metric_field<K> {
public:
	typedef typename K::FT						FT;
	typedef typename Metric_field<K>::Metric	Metric;
	typedef typename Metric_field<K>::Point_3	Point_3;
	typedef typename std::vector<std::pair<Point_3, FT> >	ParticleList;

	Metric_field<K> *base_field;
public:
	virtual Metric compute_metric(const Point_3 &p) const {
		if ((p.x() - 0.5) * (p.x() - 0.5) + p.y() * p.y() + p.z() * p.z() < 0.3 * 0.3)
			return Metric();
		else
			return base_field->compute_metric(p);
	}
	Hole_metric_field(Metric_field<K> *base_field_) : base_field(base_field_) { }
};
*/

/*
template<typename K>
class Inrimage_metric_field : public Metric_field<K> {
public:
	typedef typename K::FT						FT;
	typedef typename Metric_field<K>::Metric	Metric;
	typedef typename Metric_field<K>::Point_3	Point_3;
	FT xmax, ymax, zmax;
public:
	virtual Metric compute_metric(const Point_3 &p) const {
		// compute the force
		Chorale::API::Vector3D inp;
		inp.x = p.x();
		inp.y = p.y();
		inp.z = p.z();
		Chorale::API::Vector3D outp = Chorale::API::Session::__NRI_Evaluate(inp);
		FT l = sqrt(outp.x * outp.x + outp.y * outp.y + outp.z * outp.z);
		FT x1, y1, z1;
		if (l < 1e-3) {
			x1 = 1.0;
			y1 = 0.0;
			z1 = 0.0;
		} else {
			x1 = outp.x / l / (l / 2.0 + 1.0);
			y1 = outp.y / l / (l / 2.0 + 1.0);
			z1 = outp.z / l / (l / 2.0 + 1.0);
		}
		FT x2, y2, z2, x3, y3, z3;
		if ((fabs(x1) < 1e-7) && (fabs(y1) < 1e-7))	{
			x2 = 1.0; y2 = 0.0; z2 = 0.0;
			x3 = 0.0; y3 = 1.0; z3 = 0.0;
		} else {
			FT r = sqrt(x1 * x1 + y1 * y1);
			x2 = -y1 / r;
			y2 = x1 / r;
			z2 = 0.0;
			FT R = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
			x3 = -x1 * z1 / r / R;
			y3 = -y1 * z1 / r / R;
			z3 = r / R;
		}
		return Metric(Vector_3(x1, y1, z1), Vector_3(x2, y2, z2), Vector_3(x3, y3, z3));
	}
	Inrimage_metric_field() {
		Chorale::API::Session::__NRI_Initialize();
		Chorale::API::Vector3D p = Chorale::API::Session::__NRI_GetBox();
		xmax = p.x; ymax = p.y; zmax = p.z;
	}
};
*/


#endif

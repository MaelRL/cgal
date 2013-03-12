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

#ifndef CGAL_ANISOTROPIC_MESH_3_TORUS_METRIC_FIELD
#define CGAL_ANISOTROPIC_MESH_3_TORUS_METRIC_FIELD

#include <CGAL/Metric_field.h>

using namespace CGAL::Anisotropic_mesh_3;



template<typename K>
class Torus_metric_field : public Metric_field<K> 
{
public:
    typedef Metric_field<K>        Base;
    typedef typename Base::FT      FT;
    typedef typename Base::Metric  Metric;
    typedef typename Base::Point_3 Point_3;
    typedef typename Base::Vector_3 Vector_3;

public:
    FT R, r;

public:
    virtual void report(typename std::ofstream &fx) const
  {
    fx << "type:   torus" << std::endl;
    fx << "eps:    " << this->epsilon << std::endl;
    fx << "R:      " << R << std::endl;
    fx << "r:      " << r << std::endl;
  }

    virtual Metric compute_metric(const Point_3 &p) const
  {
    FT x = p.x(), y = p.y(), z = p.z();
    FT ll = sqrt(x * x + y * y);
    FT l = R / ll;
    FT x0 = x * l, y0 = y * l;

    FT dx = x - x0, dy = y - y0, dz = z;
    FT dl = sqrt(dx * dx + dy * dy + dz * dz);
    Vector_3 g(dx / dl, dy / dl, dz / dl);
          
    Vector_3 px(-y / ll, x / ll, 0);
    Vector_3 t = CGAL::cross_product(g, px);

    //return Metric(g, 
    //  Vector_3(px.x() / ll, px.y() / ll, px.z() / ll), 
    //  Vector_3(t.x() / r,   t.y() / r,   t.z() / r));

    Vector_3 n = g;
    Vector_3 v1 = (1./std::sqrt(px*px)) * px;
    Vector_3 v2 = (1./std::sqrt(t*t)) * t;

    bool is_closer = (ll < R);
    FT Phi = std::abs(std::asin(z/r));
    if(is_closer)
      Phi = M_PI - Phi;

      //correct curvature
    //FT e1 = std::max(std::abs(std::cos(Phi)/(R+r*std::cos(Phi))), this->epsilon );
    //FT e2 = std::max(1./r, this->epsilon);

      //naive curvature
    FT e1 = (std::max)(1./R, this->epsilon);
    FT e2 = (std::max)(1./r, this->epsilon);
    
      //good looking curvature
    //FT e1 = std::max(1./(R-r*std::cos(Phi)), this->epsilon);
    //FT e2 = std::max(1./r, this->epsilon);

    FT en = (std::max)(e1, e2);
    return Metric(n, v2, v1, en, e2, e1, this->epsilon);
  }

  Torus_metric_field(FT R_ = 0.7, FT r_ = 0.3, FT epsilon_ = 1.0) 
   : Metric_field<K>(epsilon_), R(R_), r(r_) { }
};

#endif //CGAL_ANISOTROPIC_MESH_3_TORUS_METRIC_FIELD

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

#ifndef CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_POINTSET_H
#define CGAL_ANISOTROPIC_MESH_3_CONSTRAIN_SURFACE_3_POINTSET_H


#include <CGAL/Default_configuration.h>

#include <CGAL/Constrain_surface_3_implicit.h>

#include <CGAL/trace.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/compute_average_spacing.h>

#include <vector>
#include <fstream>
#include <limits>  

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K>
class Constrain_surface_3_pointset : public Constrain_surface_3_implicit<K>
{
public:
  typedef Constrain_surface_3_implicit<K>                   Base;
  typedef typename Base::FT                                 FT;
  typedef typename Base::Point_3                            Point_3;

  typedef typename K::Sphere_3                              Sphere_3;
  typedef typename CGAL::Point_with_normal_3<K>             Point_with_normal;
  typedef typename std::vector<Point_with_normal>           PointList;
  typedef typename CGAL::Poisson_reconstruction_function<K> Poisson_reconstruction_function;
  typedef typename CGAL::Implicit_surface_3<K, Poisson_reconstruction_function> Surface_3;

  typedef typename Constrain_surface_3<K>::Pointset Point_container;

protected:
  PointList points;
  Poisson_reconstruction_function *function;
  FT radius;

public:
  FT get_bounding_radius() const { return radius; }

  virtual typename CGAL::Bbox_3 get_bbox() const 
  {
    double xmin = 0.5*DBL_MAX;
    double ymin = 0.5*DBL_MAX;
    double zmin = 0.5*DBL_MAX;
    double xmax = -0.5*DBL_MAX;
    double ymax = -0.5*DBL_MAX;
    double zmax = -0.5*DBL_MAX;
    typename PointList::const_iterator it;
    for(it = points.begin(); it != points.end(); ++it)
    {
      typename K::Point_3 p = *it; 
      xmin = (std::min)(xmin, p.x());
      xmax = (std::max)(xmax, p.x());
      ymin = (std::min)(ymin, p.y());
      ymax = (std::max)(ymax, p.y());
      zmin = (std::min)(zmin, p.z());
      zmax = (std::max)(zmax, p.z());
    }
    return CGAL::Bbox_3(xmin, ymin, zmin, xmax, ymax, zmax);
  }

  FT evaluate(const FT x, const FT y, const FT z) const
  {
    //FT val = (*function)(Point_3(x, y, z));
    //std::cout << "(" << x << ", " << y << ", " << z << ") => " << val << std::endl;
    return (*function)(Point_3(x, y, z));
  }

  virtual Point_container initial_points(const int nb = 8) const
  {
    Point_container points;
    std::vector<Point_3> seeds;
    Point_3 p = function->get_inner_point();
    std::cout << "Inner point: " << p << std::endl;
    seeds.push_back(p);
    return Base::initial_points(points, seeds, 0.2);
  }

  Constrain_surface_3_pointset* clone() const //covariant return types
  {
    return new Constrain_surface_3_pointset(*this);
  }

  Constrain_surface_3_pointset(char *filename)
  {
    std::cout << "Loading point set..." << std::flush;
    std::ifstream stream(filename);
    if (!stream ||
        !CGAL::read_xyz_points_and_normals(
                              stream,
                              std::back_inserter(points),
                              CGAL::make_normal_of_point_with_normal_pmap(std::back_inserter(points))))
    {
      std::cerr << "Error: cannot read file " << filename << std::endl;
    }

    std::cout << "\nCreating implicit function..." << std::flush;
    // Creates implicit function from the read points using the default solver (TAUCS).
    // Note: this method requires an iterator over points
    // + property maps to access each point's position and normal.
    // The position property map can be omitted here as we use iterators over Point_3 elements.
    function = new Poisson_reconstruction_function(
                              points.begin(), points.end(),
                              CGAL::make_normal_of_point_with_normal_pmap(points.begin()));

    // Computes the Poisson indicator function f()
    // at each vertex of the triangulation.
    if ( ! function->compute_implicit_function() )
      std::cout << "Error in compute_implicit_function!" << std::endl;

    // Computes average spacing
    FT average_spacing = CGAL::compute_average_spacing(points.begin(), points.end(),
                                                       6 /* knn = 1 ring */);

    // Gets one point inside the implicit surface
    // and computes implicit function bounding sphere radius.
    Sphere_3 bsphere = function->bounding_sphere();
    radius = std::sqrt(bsphere.squared_radius()) * 1.1;
    std::cout << "Bounding radius: " << radius << std::endl;
    std::cout << "\nPoisson surface construction done." << std::endl << std::endl;
  }

  ~Constrain_surface_3_pointset() { }

};

} //Anisotropic_mesh_3
} //CGAL

#endif

#ifndef CGAL_ANISOTROPIC_MESH_2_DOMAIN_2_H
#define CGAL_ANISOTROPIC_MESH_2_DOMAIN_2_H

#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>

#include <CGAL/enum.h>
#include <CGAL/Object.h>
#include <vector>
#include <set>

#include <Eigen/Dense>

namespace CGAL
{
namespace Anisotropic_mesh_2
{
template<typename K,
         typename Point_container = std::vector<typename K::Point_2> >
class Domain_2
{

public:
  typedef typename K::FT                       FT;
  typedef typename K::Vector_2                 Vector_2;
  typedef typename K::Point_2                  Point_2;
  typedef typename K::Object_2                 Object_2;
  typedef typename K::Segment_2                Segment_2;
  typedef typename K::Triangle_2               Triangle_2;
  typedef typename K::Ray_2                    Ray_2;
  typedef typename K::Line_2                   Line_2;

  typedef typename CGAL::Oriented_side         Oriented_side;
  typedef int                                  Subdomain_index;
  typedef int                                  Surface_patch_index;
  typedef int                                  Index;

  typedef Point_container                      Pointset;

public:
  virtual FT get_bounding_radius() const = 0;
  virtual Oriented_side side_of_constraint(const Point_2 &p) const = 0;
  virtual typename CGAL::Bbox_2 get_bbox() const = 0;
  virtual Point_container get_boundary_points(unsigned int nb, double facet_distance_coeff = 0.01) const = 0;
  virtual std::string name() const  = 0;

  virtual void gl_draw_intermediate_mesh_2() const
  {
    std::cout << "draw intermediate mesh3: incompatible surface" << std::endl;
    return;
  }

  template<typename C2T3>
  void output_points(const C2T3& c2t3,
                     std::ofstream& fx/*.pts file*/) const
  {

  }

  Domain_2(){ }

  virtual ~Domain_2()
  { }
};

} // Anisotropic_mesh_2
} // CGAL

#endif //DOMAIN_2_H

// Copyright (c) 2004-2009  INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Laurent RINEAU, Stephane Tayeb


#ifndef CGAL_MESH_CELL_CRITERIA_3_H
#define CGAL_MESH_CELL_CRITERIA_3_H

#include <iostream>
#include <CGAL/Mesh_3/mesh_standard_cell_criteria.h>

namespace CGAL {
  
template <typename Tr, typename Visitor_ = Mesh_3::Cell_criterion_visitor<Tr> >
class Mesh_cell_criteria_3
{
  typedef Visitor_ Visitor;
  typedef Mesh_3::Criteria<Tr,Visitor> Criteria;
  
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Geom_traits::FT FT;
  
  typedef Mesh_cell_criteria_3<Tr> Self;
  
public:
  typedef typename Visitor::Cell_quality Cell_quality;
  typedef typename Visitor::Cell_badness Cell_badness;
  
  /**
   * @brief Constructor
   * @param radius_edge_bound the radius-edge bound
   * @param radius_bound the radius bound (tet sizing)
   */
  Mesh_cell_criteria_3(const FT& radius_edge_bound,
                       const FT& radius_bound)
  {
    typedef Mesh_3::Cell_radius_criterion<Tr,Visitor> Radius_criterion;
    typedef Mesh_3::Cell_radius_edge_criterion<Tr,Visitor> Radius_edge_criterion;
    
    if ( 0 != radius_bound )
      criteria_.add(new Radius_criterion(radius_bound));
    
    if ( 0 != radius_edge_bound )
      criteria_.add(new Radius_edge_criterion(radius_edge_bound));
  }
  
  /// Destructor
  ~Mesh_cell_criteria_3() { }
  
  /**
   * @brief returns the badness of cell \c cell
   * @param cell the cell
   * @return the badness of \c cell
   */
  Cell_badness operator()(const Cell_handle& cell) const
  {
    return criteria_(cell);
  }
  
private:
  Criteria criteria_;
  
};  // end class Mesh_cell_criteria_3

}  // end namespace CGAL


#endif // CGAL_MESH_CELL_CRITERIA_3_H


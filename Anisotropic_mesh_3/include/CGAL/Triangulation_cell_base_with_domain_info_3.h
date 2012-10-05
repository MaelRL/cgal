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
// Author(s) : 

#ifndef CGAL_TRIANGULATION_CELL_BASE_WITH_DOMAIN_INFO_3_H
#define CGAL_TRIANGULATION_CELL_BASE_WITH_DOMAIN_INFO_3_H

#include <CGAL/Triangulation_cell_base_with_circumcenter_3.h>
#include <CGAL/Mesh_3/Mesh_surface_cell_base_3.h>

namespace CGAL {

template < typename GT,
           typename Constrain_surface,
           typename Cell_base = Triangulation_cell_base_3<GT> >
class Triangulation_cell_base_with_domain_info_3
  : public Mesh_3::Mesh_surface_cell_base_3<GT, 
                      Constrain_surface, 
                      Triangulation_cell_base_with_circumcenter_3<GT, Cell_base> >
{
public:
  typedef typename GT::FT					      FT;
  typedef typename GT::Point_3                                        Point_3;
  typedef typename GT::K                                              K;
  typedef typename GT::Exact_kernel                                   KExact;
  typedef Triangulation_cell_base_with_circumcenter_3<GT,Cell_base>   Cb0;
  typedef typename Mesh_3::Mesh_surface_cell_base_3<GT, Constrain_surface, Cb0> Cb;

  typedef typename Cb::Vertex_handle                   Vertex_handle;
  typedef typename Cb::Cell_handle                     Cell_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Cb::template Rebind_TDS<TDS2>::Other       Cb2;
    typedef Triangulation_cell_base_with_domain_info_3<GT, Constrain_surface, Cb2>  Other;
  };


  Triangulation_cell_base_with_domain_info_3()
    : Cb(),
	  inside_info_valid(false),
	  is_score_cache_valid(false),
          circumcenter_valid(false),
	  cell_squared_circumradius(-1.) {}

  Triangulation_cell_base_with_domain_info_3(const Triangulation_cell_base_with_domain_info_3 &c)
	: Cb(c),
	  inside_info_valid(c._get_raw_inside_info_valid()),
	  is_score_cache_valid(c._get_is_score_cache_valid()),
	  inside(c._get_raw_inside()),
          circumcenter_valid(c._get_circumcenter_valid()),
	  cell_squared_circumradius(c._get_cell_squared_circumradius()){}

  Triangulation_cell_base_with_domain_info_3 & operator= (const Triangulation_cell_base_with_domain_info_3 &c){
	  Triangulation_cell_base_with_domain_info_3 tmp = c;
	  std::swap(tmp,*this);
	  return *this;
  }

  Triangulation_cell_base_with_domain_info_3(Vertex_handle v0, Vertex_handle v1,
                                      Vertex_handle v2, Vertex_handle v3)
    : Cb(v0, v1, v2, v3),
	  inside_info_valid(false),
	  is_score_cache_valid(false),
          circumcenter_valid(false),
	  cell_squared_circumradius(-1.) {}

  Triangulation_cell_base_with_domain_info_3(Vertex_handle v0, Vertex_handle v1,
                                      Vertex_handle v2, Vertex_handle v3,
                                      Cell_handle   n0, Cell_handle   n1,
                                      Cell_handle   n2, Cell_handle   n3)
    : Cb(v0, v1, v2, v3, n0, n1, n2, n3),
	  inside_info_valid(false),
	  is_score_cache_valid(false),
          circumcenter_valid(false),
	  cell_squared_circumradius(-1.) {}

  // We must override the functions that modify the vertices.
  // And if the point inside a vertex is modified, we fail,
  // but there's not much we can do for this now.
  void set_vertex(int i, const Vertex_handle& v)
  {
      invalidate_info();
      Cb::set_vertex(i, v);
  }

  void set_vertices()
  {
      invalidate_info();
      Cb::set_vertices();
  }

  void set_vertices(const Vertex_handle& v0, const Vertex_handle& v1,
                    const Vertex_handle& v2, const Vertex_handle& v3)
  {
      invalidate_info();
      Cb::set_vertices(v0, v1, v2, v3);
  }
 
  /// \return The cell is considered as inside the constrain iff the circumcenter is inside 
  /// the domain.
  template<typename Umbrella>
  bool is_inside(const Umbrella &umb,
                 const Constrain_surface& constrain_surface, 
		 const GT &traits = GT()) const
  {
    if(inside_info_valid)
      return inside;

    if(umb.dimension() < 3)
      inside = false;
    else if(umb.is_infinite_vertex(this->vertex(0))
      || umb.is_infinite_vertex(this->vertex(1))
      || umb.is_infinite_vertex(this->vertex(2))
      || umb.is_infinite_vertex(this->vertex(3)))
      inside = false;
    else
    {
      Point_3 tcc = traits.construct_circumcenter_3_object()
        (this->vertex(0)->point(), this->vertex(1)->point(),
         this->vertex(2)->point(), this->vertex(3)->point());
      Point_3 cc = umb.metric().inverse_transform(tcc);
      inside = (constrain_surface.side_of_constraint(cc) == ON_POSITIVE_SIDE);
    }
    inside_info_valid = true;
    return inside;
  }

#ifdef ANISO_USE_INSIDE_EXACT
  template<typename Umbrella>
  bool is_inside_exact(const Umbrella &umb,
                       const Constrain_surface& constrain_surface, 
		       const GT &traits = GT()) const
  {
//    if(inside_info_valid)
//      return inside;

    if(umb.dimension() < 3)
      inside = false;
    else if(umb.is_infinite_vertex(this->vertex(0))
      || umb.is_infinite_vertex(this->vertex(1))
      || umb.is_infinite_vertex(this->vertex(2))
      || umb.is_infinite_vertex(this->vertex(3)))
      inside = false;
    else
    {
      typename KExact::Point_3 tcc = traits.construct_circumcenter_3_object()
          (to_exact(this->vertex(0)->point()), 
           to_exact(this->vertex(1)->point()),
           to_exact(this->vertex(2)->point()), 
           to_exact(this->vertex(3)->point()));
      typename KExact::Point_3 cc = umb.metric().inverse_transform(tcc);
      inside = (constrain_surface.side_of_constraint(back_from_exact(cc)) == ON_POSITIVE_SIDE);
    }
    inside_info_valid = true;
    return inside;
  }
#endif

  /// The squared_circumradius of the cell
  FT squared_circumradius(const GT &traits = GT()){
	  if(cell_squared_circumradius < 0){
		  typename GT::Compute_squared_distance_3 csd = traits.compute_squared_distance_3_object();
#ifndef ANISO_NO_USE_CC_CACHE
                  cell_squared_circumradius  = csd(this->vertex(0)->point(),this->circumcenter(traits));
#else
                  cell_squared_circumradius  = csd(this->vertex(0)->point(),Cb::circumcenter(traits));                  
#endif
	  }
	  return cell_squared_circumradius;
  }

  void invalidate_info(){
	  cell_squared_circumradius = -1;
	  inside_info_valid = false;
          circumcenter_valid = false;
  }


  template<typename Volume_criteria,typename Trait>
  bool is_bad(const Volume_criteria &vc,const Trait &trait){
	  return vc.is_bad(this,trait);
  }

  template<typename Volume_criteria,typename Trait>
  FT score(const Volume_criteria &vc,const Trait &trait){
	  if(!is_score_cache_valid){
		  score_cache = vc.score(this,trait);
		  is_score_cache_valid = true;
	  }
	  return score_cache;
  }

#ifndef ANISO_NO_USE_CC_CACHE
  const Point_3& circumcenter(const GT &traits = GT()) const 
  {
    if(circumcenter_valid)
      return circumcenter_cache;
    
    circumcenter_cache = Cb::circumcenter(traits);
    circumcenter_valid = true;
    return circumcenter_cache;
  }
#endif

private:
//=================================
  bool _get_raw_inside(){
	  return inside;
  }

  bool _get_raw_inside_info_valid(){
	  return inside_info_valid;
  }

  FT _get_cell_squared_circumradius(){
	  return cell_squared_circumradius;
  }
  
  bool _get_is_score_cache_valid(){
	  return is_score_cache_valid;
  }

  bool _get_circumcenter_valid(){
      return circumcenter_valid;
  }
//=================================

  private:
    typename CGAL::Cartesian_converter<K, KExact> to_exact;
    typename CGAL::Cartesian_converter<KExact, K> back_from_exact;

//=================================  
  mutable bool	inside;
  mutable bool	inside_info_valid;
  bool	is_score_cache_valid;
  FT	cell_squared_circumradius;
  FT	score_cache;
  mutable Point_3 circumcenter_cache;
  mutable bool circumcenter_valid;
};

} //namespace CGAL

#endif // CGAL_TRIANGULATION_CELL_BASE_WITH_DOMAIN_INFO_3_H 

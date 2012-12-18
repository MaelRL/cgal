#ifndef CGAL_REGULAR_TRIANGULATION_FACE_BASE_WITH_BOOL_ON_SPHERE_2_H
#define CGAL_REGULAR_TRIANGULATION_FACE_BASE_WITH_BOOL_ON_SPHERE_2_H

#include <list>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_face_base_sphere_2.h>

namespace CGAL {


template <class Gt, class Fb = Triangulation_face_base_sphere_2<Gt> >
class Regular_triangulation_face_base_with_bool_on_sphere_2
  :  public Fb
{
  typedef Fb                                            Fbase;
  typedef typename Fbase::Triangulation_data_structure  TDS;
public:
  typedef Gt                                   Geom_traits;
  typedef TDS                                  Triangulation_data_structure;
  typedef typename TDS::Vertex_handle          Vertex_handle;
  typedef typename TDS::Face_handle            Face_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other   Fb2;
    typedef Regular_triangulation_face_base_with_bool_on_sphere_2<Gt,Fb2>             Other;
  };

  typedef std::list<Vertex_handle>             Vertex_list;
  
  bool is_inside;

protected:
  Vertex_list vlist;
  unsigned char _in_conflict_flag;

public:
  void set_in_conflict_flag(unsigned char f) { _in_conflict_flag = f; }
  unsigned char get_in_conflict_flag() const { return _in_conflict_flag; }

  Regular_triangulation_face_base_with_bool_on_sphere_2()
  : Fbase(),  is_inside(true)
  {
    
    set_in_conflict_flag(0);
  }

  Regular_triangulation_face_base_with_bool_on_sphere_2(Vertex_handle v0, 
				    Vertex_handle v1, 
				    Vertex_handle v2)
  : Fbase(v0,v1,v2),  is_inside(true)
  {
    set_in_conflict_flag(0);
  }

  Regular_triangulation_face_base_with_bool_on_sphere_2(Vertex_handle v0, 
				    Vertex_handle v1, 
				    Vertex_handle v2,
				    Face_handle n0, 
				    Face_handle n1, 
				    Face_handle n2)
  : Fbase(v0,v1,v2,n0,n1,n2),  is_inside(true)
  { 
       set_in_conflict_flag(0);
  }

  Vertex_list& vertex_list()
  {
    return vlist;
  }



};

} //namespace CGAL 

#endif // CGAL_REGULAR_TRIANGULATION_FACE_BASE_WITH_BOOL_ON_SPHERE_2_H


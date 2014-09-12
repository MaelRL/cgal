#ifndef CGAL_TRIANGULATION_FACE_BASE_WITH_DOMAIN_INFO_2_H
#define CGAL_TRIANGULATION_FACE_BASE_WITH_DOMAIN_INFO_2_H

#include <CGAL/Delaunay_mesh_face_base_2.h>

namespace CGAL
{

template < typename GT,
           typename Domain,
           typename Fb = CGAL::Delaunay_mesh_face_base_2<GT> >
class Triangulation_face_base_with_domain_info_2
    : public Fb
{
public:
  typedef typename GT::FT                                                       FT;
  typedef typename GT::Point_2                                                  Point_2;
  typedef typename GT::K                                                        K;
  typedef typename GT::Exact_kernel                                             KExact;

  typedef typename Fb::Triangulation_data_structure                             TDS;
  typedef typename TDS::Vertex_handle                                           Vertex_handle;
  typedef typename TDS::Face_handle                                             Face_handle;

  template < typename TDS2 >
  struct Rebind_TDS
  {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other     Fb2;
    typedef Triangulation_face_base_with_domain_info_2<GT, Domain, Fb2> Other;
  };

  Triangulation_face_base_with_domain_info_2()
    :
      Fb(),
      inside_info_valid(false),
      face_squared_circumradius(-1.),
      circumcenter_valid(false)
    { }

  Triangulation_face_base_with_domain_info_2(const Triangulation_face_base_with_domain_info_2 &f)
    :
      Fb(f),
      inside(f._get_raw_inside()),
      inside_info_valid(f._get_raw_inside_info_valid()),
      face_squared_circumradius(f._get_face_squared_circumradius()),
      circumcenter_valid(f._get_circumcenter_valid())
    { }

  Triangulation_face_base_with_domain_info_2 & operator= (const Triangulation_face_base_with_domain_info_2 &c)
  {
    Triangulation_face_base_with_domain_info_2 tmp = c;
    std::swap(tmp,*this);
    return *this;
  }

  Triangulation_face_base_with_domain_info_2(Vertex_handle v0, Vertex_handle v1,
                                             Vertex_handle v2)
    :
      Fb(v0, v1, v2),
      inside_info_valid(false),
      face_squared_circumradius(-1.),
      circumcenter_valid(false)
    { }

  Triangulation_face_base_with_domain_info_2(Vertex_handle v0, Vertex_handle v1,
                                             Vertex_handle v2, Face_handle n0,
                                             Face_handle n1, Face_handle n2)
    :
      Fb(v0, v1, v2, n0, n1, n2),
      inside_info_valid(false),
      face_squared_circumradius(-1.),
      circumcenter_valid(false)
    { }

  // We must override the functions that modify the vertices.
  // And if the point inside a vertex is modified, we fail,
  // but there's not much we can do for this now.
  void set_vertex(int i, const Vertex_handle& v)
  {
    invalidate_info();
    Fb::set_vertex(i, v);
  }

  void set_vertices()
  {
    invalidate_info();
    Fb::set_vertices();
  }

  void set_vertices(const Vertex_handle& v0, const Vertex_handle& v1,
                    const Vertex_handle& v2)
  {
    invalidate_info();
    Fb::set_vertices(v0, v1, v2);
  }

  template<typename Umbrella>
  bool is_inside(const Umbrella &umb,
                 const Domain& domain,
                 const GT &traits = GT()) const
  {
    if(inside_info_valid)
      return inside;

    if(umb.dimension() < 2)
      inside = false;
    else if(umb.is_infinite_vertex(this->vertex(0))
            || umb.is_infinite_vertex(this->vertex(1))
            || umb.is_infinite_vertex(this->vertex(2)))
      inside = false;

    Point_2 tcc = traits.construct_circumcenter_2_object() (this->vertex(0)->point(),
                                                            this->vertex(1)->point(),
                                                            this->vertex(2)->point());
    Point_2 cc = umb.metric().inverse_transform(tcc);
    Oriented_side side = domain.side_of_constraint(cc);

#ifdef ANISO_DEBUG_IS_INSIDE
    std::cout << "is inside? ";
    std::cout << umb.metric().inverse_transform(this->vertex(0)->point()) << " XXX ";
    std::cout << umb.metric().inverse_transform(this->vertex(1)->point()) << " XXX ";
    std::cout << umb.metric().inverse_transform(this->vertex(2)->point()) << " XXX ";
    std::cout << cc << " --------> " << side << std::endl;
#endif

    inside = (side == ON_POSITIVE_SIDE || side == ON_ORIENTED_BOUNDARY);

    inside_info_valid = true;
    return inside;
  }

#ifdef ANISO_USE_INSIDE_EXACT
  template<typename Umbrella>
  bool is_inside_exact(const Umbrella &umb,
                       const Domain& domain,
                       const GT &traits = GT()) const
  {
    if(inside_info_valid)
      return inside;

    if(umb.dimension() < 2)
      inside = false;
    else if(umb.is_infinite_vertex(this->vertex(0))
            || umb.is_infinite_vertex(this->vertex(1))
            || umb.is_infinite_vertex(this->vertex(2)))
      inside = false;
    else
    {
      typename KExact::Point_2 tcc = traits.construct_circumcenter_2_object()
                                     (to_exact(this->vertex(0)->point()),
                                      to_exact(this->vertex(1)->point()),
                                      to_exact(this->vertex(2)->point()));
      typename KExact::Point_2 cc = umb.metric().inverse_transform(tcc);
      Oriented_side side = domain.side_of_constraint(back_from_exact(cc));
      inside = (side == ON_POSITIVE_SIDE || side == ON_ORIENTED_BOUNDARY);
    }
    inside_info_valid = true;
    return inside;
  }
#endif

  FT squared_circumradius(const GT &traits = GT())
  {
    if(face_squared_circumradius < 0)
    {
      typename GT::Compute_squared_distance_2 csd =
                                     traits.compute_squared_distance_2_object();
#ifndef ANISO_NO_USE_CC_CACHE
      face_squared_circumradius = csd(this->vertex(0)->point(),
                                      this->circumcenter(traits));
#else
      face_squared_circumradius = csd(this->vertex(0)->point(),
                                      Fb::circumcenter(traits));
#endif
    }
    return face_squared_circumradius;
  }

  void invalidate_info()
  {
    face_squared_circumradius = -1;
    inside_info_valid = false;
    circumcenter_valid = false;
  }

#ifndef ANISO_NO_USE_CC_CACHE
  const Point_2& circumcenter(const GT& traits = GT()) const
  {
    if(circumcenter_valid)
      return circumcenter_cache;

    circumcenter_cache = traits.construct_circumcenter_2_object()(this->vertex(0)->point(),
                                                                  this->vertex(1)->point(),
                                                                  this->vertex(2)->point());
    circumcenter_valid = true;
    return circumcenter_cache;
  }
#endif

private:
//=================================
  bool _get_raw_inside() const { return inside; }
  bool _get_raw_inside_info_valid() const { return inside_info_valid; }
  FT _get_face_squared_circumradius() const { return face_squared_circumradius; }
  bool _get_circumcenter_valid() const { return circumcenter_valid; }
//=================================

private:
  typename CGAL::Cartesian_converter<K, KExact> to_exact;
  typename CGAL::Cartesian_converter<KExact, K> back_from_exact;

//=================================
  mutable bool  inside;
  mutable bool  inside_info_valid;
  FT    face_squared_circumradius;
  mutable Point_2 circumcenter_cache;
  mutable bool circumcenter_valid;
};

} //namespace CGAL

#endif // CGAL_TRIANGULATION_FACE_BASE_WITH_DOMAIN_INFO_2_H

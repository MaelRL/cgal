#ifndef CGAL_ANISOTROPIC_MESH_TC_STAR
#define CGAL_ANISOTROPIC_MESH_TC_STAR

#include <CGAL/Metric_field.h>

#include <CGAL/Regular_triangulation.h>
#include <CGAL/Regular_triangulation_euclidean_traits.h>
#include <CGAL/Full_cell_refine_queue.h>

#include <CGAL/helpers/combinatorics_helper.h>

#include <CGAL/Combination_enumerator.h>

namespace CGAL
{
namespace Anisotropic_mesh_TC
{

template<typename Kd, typename KD>
class Tangent_star
    : public CGAL::Regular_triangulation<
                   CGAL::Regular_triangulation_euclidean_traits<Kd>,
                   typename CGAL::Triangulation_data_structure<
                                  typename Kd::Dimension,
                                  CGAL::Triangulation_vertex<
                                        CGAL::Regular_triangulation_euclidean_traits<Kd>,
                                        std::size_t>,
                                  CGAL::Triangulation_full_cell<
                                        CGAL::Regular_triangulation_euclidean_traits<Kd>,
                                        std::pair<typename Kd::Point_d, bool> > > >
{
  typedef Tangent_star<Kd, KD>                                     Self;
  typedef Self*                                                    Star_handle;

public:
  typedef int                                                      Index;

  typedef typename Kd::Dimension                                   dDim;
  typedef typename KD::Dimension                                   DDim;
  typedef typename Kd::FT                                          FT;
  typedef typename Kd::Point_d                                     Point_d;
  typedef typename Kd::Weighted_point_d                            WPoint_d;
  typedef typename Kd::Vector_d                                    Vector_d;
  typedef typename KD::Point_d                                     Point_D;
  typedef typename KD::Weighted_point_d                            WPoint_D;
  typedef typename KD::Vector_d                                    Vector_D;
  typedef std::vector<Vector_D>                                    Tangent_space_basis;

  typedef CGAL::Regular_triangulation_euclidean_traits<Kd>         Traits;
  typedef typename CGAL::Triangulation_data_structure<typename Kd::Dimension,
          CGAL::Triangulation_vertex<Traits, std::size_t>,
          CGAL::Triangulation_full_cell<Traits, std::pair<Point_d, bool> > >
                                                                   TDS;
  typedef CGAL::Regular_triangulation<Traits, TDS>                 Base;

  typedef typename TDS::Vertex                                     Vertex;
  typedef typename TDS::Face                                       Face;
  typedef typename TDS::Full_cell                                  Full_cell;

  typedef typename TDS::Vertex_handle                              Vertex_handle;
  typedef typename TDS::Full_cell_handle                           Full_cell_handle;
  typedef typename TDS::Full_cell::Vertex_handle_iterator          Vertex_h_iterator;

  typedef typename std::vector<Face>                               Face_vector;
  typedef typename Face_vector::iterator                           Face_iterator;
  typedef typename std::vector<Vertex_handle>                      Vertex_handle_vector;
  typedef typename Vertex_handle_vector::iterator                  Vertex_handle_iterator;
  typedef typename std::vector<Full_cell_handle>                   Full_cell_handle_vector;
  typedef typename Full_cell_handle_vector::iterator               Full_cell_handle_iterator;

  typedef Metric_base<Kd>                                          Metric;

  typedef typename Metric::E_Matrix                                E_Matrix_d;
  typedef typename Metric::E_Vector                                E_Vector_d;
  typedef Eigen::Matrix<FT, DDim::value, DDim::value>              E_Matrix_D;
  typedef Eigen::Matrix<FT, DDim::value, 1>                        E_Vector_D;
  typedef Eigen::Matrix<FT, DDim::value, dDim::value>              E_Matrix_Dd;

  typedef Ordered_simplex_base<dDim::value+1>                      Simplex;

  // members
  Traits* m_traits;
  Metric m_metric; // metric at m_center

  mutable bool is_tsb_init;
  mutable Tangent_space_basis m_tsb;

  Point_d m_center; // base point
  Point_D m_center_Q; // ... transformed to Q (the paraboloid)
  WPoint_D m_center_S; // ... transformed to S (the metric surface)
  WPoint_d m_center_T; // S projected on the tangent plane
  Vertex_handle m_center_v; // base point in rt

  static const Index index_of_infinite_vertex = -10;

  mutable bool m_is_cache_dirty;
  mutable Vertex_handle_vector m_finite_adjacent_vertices_cache;
  mutable Full_cell_handle_vector m_incident_full_cells_cache;
  mutable Full_cell_handle_vector m_finite_incident_full_cells_cache;

  // get/set etc. --------------------------------------------------------------
  const int d() const { return dDim::value; }
  const int D() const { return DDim::value; }
  const Index index() const { return m_center_v->data(); }

  const Point_d& center_point() const { return m_center; }
  const Metric& metric() const { return m_metric; }

  // cache related functions ---------------------------------------------------
  void invalidate_cache()
  {
    m_is_cache_dirty = true;
  }

  void update_star_caches() const
  {
    if(!m_is_cache_dirty)
      return;

    m_finite_adjacent_vertices_cache.clear();
    m_incident_full_cells_cache.clear();
    m_finite_incident_full_cells_cache.clear();

    // update adjacent vertices todo (need to add it to TDS_d...)
//    std::back_insert_iterator<Vertex_handle_vector>
//        finite_adjacent_vertices_insertor(m_finite_adjacent_vertices_cache);
//    Base::finite_adjacent_vertices(m_center_v, finite_adjacent_vertices_insertor);

    // and cells
    std::back_insert_iterator<Full_cell_handle_vector>
        incident_cells_insertor(m_incident_full_cells_cache);
    Base::incident_full_cells(m_center_v, incident_cells_insertor);

    Full_cell_handle_iterator cit = m_incident_full_cells_cache.begin();
    Full_cell_handle_iterator cend = m_incident_full_cells_cache.end();
    for(; cit!=cend; ++cit)
    {
      if(!Base::is_infinite(*cit))
        m_finite_incident_full_cells_cache.push_back(*cit);
    }

    m_is_cache_dirty = false;
  }

  inline Vertex_handle_iterator finite_adjacent_vertices_begin() const
  {
    update_star_caches();
    return m_finite_adjacent_vertices_cache.begin();
  }

  inline Vertex_handle_iterator finite_adjacent_vertices_end() const
  {
    return m_finite_adjacent_vertices_cache.end();
  }

  inline Full_cell_handle_iterator incident_full_cells_begin() const
  {
    update_star_caches();
    return m_incident_full_cells_cache.begin();
  }

  inline  Full_cell_handle_iterator incident_full_cells_end() const
  {
    return m_incident_full_cells_cache.end();
  }

  inline Full_cell_handle_iterator finite_incident_full_cells_begin() const
  {
    update_star_caches();
    return m_finite_incident_full_cells_cache.begin();
  }

  inline  Full_cell_handle_iterator finite_incident_full_cells_end() const
  {
    return m_finite_incident_full_cells_cache.end();
  }

  // tangent space functions ---------------------------------------------------
  void compute_tangent_space_basis() const
  {
    // a basis of the tangent space is given by the partial derivatives
    // of the paraboloid's parametrization evaluated at m_center

    E_Matrix_Dd tsb_m = E_Matrix_Dd::Zero();

    // 'first' part : x y z etc. derivates
    for(int i=0; i<d(); ++i)
      for(int j=0; j<d(); ++j)
        tsb_m(i,j) = (i==j);

    // 'second' part: x² xy xz y² yz z² etc. derivatives
    int pos = d();
    typename Kd::Compute_coordinate_d coord = Kd().compute_coordinate_d_object();

    for(int i=0; i<d(); ++i)
    {
      for(int j=i; j<d(); ++j)
      {
        if(i == j)
          tsb_m(pos, i) = 2*coord(m_center, i);
        else
        {
          tsb_m(pos, i) = coord(m_center, j);
          tsb_m(pos, j) = coord(m_center, i);
        }
        pos++;
      }
    }

//GRAM SCHMIDT ------------------------------------- todo check that it is needed
    FT nu = (tsb_m.col(0)).norm();
    std::cout << "nu: " << nu << std::endl;
    FT sp = tsb_m.col(0).dot(tsb_m.col(1));
    tsb_m.col(0) /= nu; // norming first => only dividing by nu below (instead of nu²)
    tsb_m.col(1) -= sp/nu*tsb_m.col(0);

    FT nv = (tsb_m.col(1)).norm();
    tsb_m.col(1) /= nv;

    std::cout << tsb_m.col(0).norm() << " " << tsb_m.col(1).norm() << std::endl;
    std::cout << tsb_m.col(0).dot(tsb_m.col(1)) << std::endl;
// ----------------------------------------

    std::cout << "tsb_m @ : ";
    std::cout << m_center_Q[0] << " " << m_center_Q[1] << " ";
    std::cout << m_center_Q[2] << " " << m_center_Q[3] << " ";
    std::cout << m_center_Q[4] << std::endl;
    std::cout << std::endl << tsb_m << std::endl;

    typename KD::Construct_vector_d constr_vec = KD().construct_vector_d_object();
    for(int i=0; i<d(); ++i)
      m_tsb[i] = (constr_vec(D(), tsb_m.col(i).data(), tsb_m.col(i).data() + D()));
  }

  // transformations -----------------------------------------------------------
  // from R^D on Q to R^d
  Point_d from_Q(const Point_D& p_on_Q) const
  {
    return Kd().construct_point_d_object()(d(), p_on_Q.begin(),
                                                p_on_Q.begin() + d());
  }

  // from a point in R^d to R^D on Q
  Point_D to_Q(const Point_d& p) const
  {
    E_Vector_D p_on_Q;
    for(int i=0; i<d(); ++i)
      p_on_Q(i) = p[i];

    int ind = d();
    typename Kd::Compute_coordinate_d coord = Kd().compute_coordinate_d_object();

    for(int i=0; i<d(); ++i)
      for(int j=i; j<d(); ++j)
        p_on_Q(ind++) = coord(p,i) * coord(p,j);

    return KD().construct_point_d_object()(D(), p_on_Q.data(), p_on_Q.data() + D());
  }

  // from R^D on the metric surface to R^d
  Point_d from_S(const WPoint_D& p_on_S) const
  {
    E_Matrix_d m;
    typename Kd::Compute_coordinate_d coord = Kd().compute_coordinate_d_object();

    int i=0, j=0;
    for(int k=d(); k<D(); k++)
    {
      m(i,j) = coord(p_on_S, k);
      m(j,i) = coord(p_on_S, k);
      j++;
      if(j == d())
      {
        i++;
        j = i;
      }
    }

    E_Vector_d tmp;
    for(int i=0; i<d; ++i)
      tmp(i) = coord(p_on_S, i);

    tmp = m.inverse() * tmp;
    return Kd().construct_point_d_object()(d(), tmp.data(),
                                                tmp.data() + d());
  }

  // from R^d to R^D on the metric surface
  WPoint_D to_S(const Point_d& p) const
  {
//    std::cout << "p: " << p[0] << " " << p[1] << std::endl;

    const E_Matrix_d m = m_metric.get_mat();
    E_Vector_d e_p, p_bar;

    typename Kd::Compute_coordinate_d coord = Kd().compute_coordinate_d_object();
    for(int i=0; i<d(); ++i)
      e_p(i) = coord(p, i);

    p_bar = m * e_p;

    E_Vector_D p_on_S;
    for(int i=0; i<d(); ++i)
      p_on_S(i) = p_bar(i);

    int ind = d();
    for(int i=0; i<d(); ++i)
    {
      for(int j=i; j<d(); ++j)
      {
        if(j==i)
          p_on_S(ind++) = -0.5*m(i,i);
        else
          p_on_S(ind++) = -m(i,j);
      }
    }

    FT n = p_on_S.norm();
    FT w = n * n - e_p.transpose() * p_bar;
//    std::cout << "n: " << n << std::endl;
//    std::cout << "weight: " << w << std::endl;

    return KD().construct_weighted_point_d_object()
      (KD().construct_point_d_object()(D(), p_on_S.data(), p_on_S.data()+D()), w);
  }

  // from R^D on the tangent plane to R^d
  Point_d from_T(const WPoint_d& p_on_T) const
  {
    //todo
  }

  // from R^d to R^D on the tangent plane
  WPoint_d to_T(const Point_d& p) const
  {
    if(!is_tsb_init)
    {
      compute_tangent_space_basis();
      is_tsb_init = true;
    }

    WPoint_D p_on_S(to_S(p));

    std::cout << "p: " << p[0] << " " << p[1] << std::endl;
    std::cout << "mq: " << m_center_Q[0] << " " << m_center_Q[1] << " " << m_center_Q[2] << " " << m_center_Q[3] << " " << m_center_Q[4] << std::endl;
    std::cout << "ps: " << p_on_S.point()[0] << " " << p_on_S.point()[1] << " " << p_on_S.point()[2] << " " << p_on_S.point()[3] << " " << p_on_S.point()[4] << " w: " << p_on_S.weight() << std::endl;

    typename KD::Scalar_product_d inner_pdct =
      KD().scalar_product_d_object();
    typename KD::Difference_of_points_d diff_points =
      KD().difference_of_points_d_object();

    Vector_D v = diff_points(p_on_S.point(), m_center_Q);

    std::cout << "v: " << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << " " << v[4] << std::endl;

    // Ambiant-space coords of the projected point
    std::vector<FT> coords;
    coords.reserve(d());
    std::vector<FT> p_proj(m_center_Q.cartesian_begin(), m_center_Q.cartesian_end());
    std::cout << "p_proj: " << p_proj[0] << " " << p_proj[1] << " " << p_proj[2] << " " << p_proj[3] << " " << p_proj[4] << std::endl;

    for(std::size_t i=0 ;i<d() ;++i)
    {
      // Compute the inner product p * ts[i]
      FT coord = inner_pdct(v, m_tsb[i]);
      coords.push_back(coord);

      // p_proj += coord * v;
      for(int j=0 ; j<D(); ++j)
        p_proj[j] += coord * m_tsb[i][j];
    }

    std::cout << "coords: " << coords.size() << " " << coords[0] << " " << coords[1] << std::endl;

    Point_D projected_pt(D(), p_proj.begin(), p_proj.end());
    std::cout << "ppt: " << projected_pt[0] << " " << projected_pt[1] << " ";
    std::cout << projected_pt[2] << " " << projected_pt[3] << " " << projected_pt[4] << std::endl;

    typename KD::Squared_distance_d sqdist = KD().squared_distance_d_object();
    std::cout << "sqdist: " << sqdist(p_on_S.point(), projected_pt) << std::endl;

    return m_traits->construct_weighted_point_d_object()
    (
      m_traits->construct_point_d_object()(d(), coords.begin(), coords.end()),
      p_on_S.weight() - sqdist(p_on_S.point(), projected_pt)
    );
  }

  // Has -----------------------------------------------------------------------
  bool has_cell(Full_cell_handle& fch,
                const boost::array<int, dDim::value+1>& a) const
  {
    const std::size_t d = dDim::value+1;
    int cids[d], dids[d];
    for(int i=0; i<d; i++)
      cids[i] = a[i];

    Full_cell_handle_iterator ci = finite_incident_full_cells_begin();
    Full_cell_handle_iterator cend = finite_incident_full_cells_end();
    for(; ci!=cend; ci++)
    {
      if((*ci)->maximal_dimension() < (d-1))
        continue;
      for(int i=0; i<d; i++)
        dids[i] = (*ci)->vertex(i)->data();
      if(is_same_ids<d>(cids, dids))
      {
        fch = *ci;
        return true;
      }
    }
    return false;
  }

  bool has_cell(const Full_cell_handle& fch) const
  {
    const std::size_t dp1 = dDim::value+1;
    int cids[dp1], dids[dp1];
    for(int i=0; i<dp1; i++)
      cids[i] = fch->vertex(i)->data();

    Full_cell_handle_iterator ci = finite_incident_full_cells_begin();
    Full_cell_handle_iterator cend = finite_incident_full_cells_end();
    for(; ci!=cend; ci++)
    {
      if((*ci)->maximal_dimension() < dDim::value)
        continue;
      for(int i=0; i<dp1; i++)
        dids[i] = (*ci)->vertex(i)->data();
      if(is_same_ids<dp1>(cids, dids))
        return true;
    }
    return false;
  }

  bool has(const std::vector<Star_handle>& cell, const std::size_t n) const
  {
    for(std::size_t i=0; i<cell.size(); ++i)
      if(cell[i]->index() == n)
        return true;
    return false;
  }

  // Insert functions ----------------------------------------------------------
  Vertex_handle insert_if_in_star(const Point_d& p,
                                  const Index& index_)
  {
    WPoint_D pp = to_S(p);
    Vertex_handle vh = Base::insert(to_T(p));
//    Vertex_handle vh = Base::insert_if_in_star(to_T(p), m_center_v);
    if(vh != Vertex_handle())
      vh->data() = index_;
    invalidate_cache();
    return vh;
  }

  Vertex_handle insert_to_star(const Point_d& p,
                               const Index& index,
                               const bool conditional = false)
  {
    // if(conditional) // todo
    // else
    Vertex_handle v = Base::insert(to_T(p), m_center_v);
    if(v != Vertex_handle())
      v->data() = index;
    invalidate_cache();
    this->is_valid(true);

    std::cout << "INSERT : " << index << " IN STAR " << m_center_v->data() << std::endl;
    std::cout << "finite cells: " << this->number_of_finite_full_cells() << std::endl;
    std::cout << "dim rt: " << this->current_dimension() << std::endl;

    return v;
  }

  std::size_t clean()
  {
    typedef typename std::pair<Point_d, int> PPoint_d;
    std::vector<PPoint_d> backup;
    std::size_t nbv = this->number_of_vertices();

    // backup star vertices
    Vertex_handle_iterator vhit = finite_adjacent_vertices_begin();
    Vertex_handle_iterator vhend = finite_adjacent_vertices_end();
    for(; vhit != vhend; vhit++)
      backup.push_back(std::make_pair((*vhit)->point().point(), (*vhit)->data()));

    this->clear();
    this->infinite_vertex()->data() = index_of_infinite_vertex;

    // re-insert
    for(std::size_t i = 1; i < backup.size(); i++)
      insert_to_star(backup[i].first, backup[i].second, false);

    return (nbv - backup.size());
 }

  void reset()
  {
    Index index = m_center_v->data();
    this->clear();
    m_center_v = Base::insert(m_center_T);
    m_center_v->data() = index;
    this->infinite_vertex()->data() = index_of_infinite_vertex;

    invalidate_cache();
  }

  void rebuild_star(); // ~ connectivity in R^D ?

public:
  Tangent_star(const Point_d& center_point,
               const Index& index_,
               //const Criteria* criteria_,
               //const Constrain_surface* pconstrain_surface,
               const Metric& metric_
               )
    :
      Base(dDim::value, *(m_traits = new Traits())),
      m_metric(metric_),
      is_tsb_init(false),
      m_tsb(d()),
      m_center(center_point),
      m_center_Q(to_Q(center_point)),
      m_center_S(to_S(center_point)),
      m_center_T(to_T(center_point)),
      m_is_cache_dirty(true),
      m_finite_adjacent_vertices_cache(),
      m_incident_full_cells_cache(),
      m_finite_incident_full_cells_cache()
  {
    m_center_v = Base::insert(m_center_T);
    m_center_v->data() = index_;
    this->infinite_vertex()->data() = index_of_infinite_vertex;

    std::cout << "init" << std::endl;
    std::cout << "finite cells: " << this->number_of_finite_full_cells() << std::endl;
    std::cout << "dim rt: " << this-> current_dimension() << std::endl;
  }

};

} // Anisotropic_mesh_TC
} // CGAL

#endif //CGAL_ANISOTROPIC_MESH_TC_STAR

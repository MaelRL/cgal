#ifndef CGAL_ANISOTROPIC_MESH_TC_STAR
#define CGAL_ANISOTROPIC_MESH_TC_STAR

#include <CGAL/Metric_field.h>

#include <CGAL/Regular_triangulation.h>
#include <CGAL/Regular_triangulation_euclidean_traits.h>
#include <CGAL/Full_cell_refine_queue.h>

#include <CGAL/helpers/combinatorics_helper.h>
#include <CGAL/Combination_enumerator.h>

#include <Eigen/Dense>

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

  Point_d m_center; // base point in R^d
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
  const std::size_t d() const { return dDim::value; }
  const std::size_t D() const { return DDim::value; }
  const Index index() const { return m_center_v->data(); }

  const Point_d& center_point() const { return m_center; }
  const Metric& metric() const { return m_metric; }

  // cache related functions ---------------------------------------------------
  void invalidate_cache() { m_is_cache_dirty = true; }

  void update_star_caches() const
  {
    assert(is_valid(true));
    if(!m_is_cache_dirty)
      return;

    std::cout << "updating cache @ " << index() << " ... ";

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

  bool is_topological_ball() const
  {
    // 2 things todo : verify no d-2 with more than 2 (d-1) on the border

    // verify there's only one umbrella
    return false;
  }

  // tangent space functions ---------------------------------------------------
  void compute_tangent_space_basis() const
  {
    // a basis of the tangent space is given by the partial derivatives
    // of the paraboloid's parametrization evaluated at m_center

    E_Matrix_Dd tsb_m = E_Matrix_Dd::Zero();

    // 'first' part : x y z etc. derivates
    for(std::size_t i=0; i<d(); ++i)
      for(std::size_t j=0; j<d(); ++j)
        tsb_m(i,j) = (i==j);

    // 'second' part: x² xy xz y² yz z² etc. derivatives
    int pos = d();
    typename Kd::Compute_coordinate_d coord = Kd().compute_coordinate_d_object();

    for(std::size_t i=0; i<d(); ++i)
    {
      for(std::size_t j=i; j<d(); ++j)
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

//GRAM SCHMIDT ----------------------------------- todo: check that it is needed
    FT nu = (tsb_m.col(0)).norm();
    FT sp = tsb_m.col(0).dot(tsb_m.col(1));
    tsb_m.col(0) /= nu; // norming first => only dividing by nu below (instead of nu²)
    tsb_m.col(1) -= sp/nu*tsb_m.col(0);

    FT nv = (tsb_m.col(1)).norm();
    tsb_m.col(1) /= nv;

//    std::cout << tsb_m.col(0).norm() << " " << tsb_m.col(1).norm() << std::endl;
//    std::cout << tsb_m.col(0).dot(tsb_m.col(1)) << std::endl;
// ----------------------------------------

    std::cout << "tangent space @ : ";
    std::cout << m_center_Q[0] << " " << m_center_Q[1] << " ";
    std::cout << m_center_Q[2] << " " << m_center_Q[3] << " ";
    std::cout << m_center_Q[4] << std::endl;
    std::cout << "tsb: " << std::endl << tsb_m.transpose() << std::endl;

    typename KD::Construct_vector_d constr_vec = KD().construct_vector_d_object();
    for(std::size_t i=0; i<d(); ++i)
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
    for(std::size_t i=0; i<d(); ++i)
      p_on_Q(i) = p[i];

    int ind = d();
    typename Kd::Compute_coordinate_d coord = Kd().compute_coordinate_d_object();

    for(std::size_t i=0; i<d(); ++i)
      for(std::size_t j=i; j<d(); ++j)
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
    const E_Matrix_d m = m_metric.get_mat();
    E_Vector_d e_p, p_bar;

    typename Kd::Compute_coordinate_d coord = Kd().compute_coordinate_d_object();
    for(std::size_t i=0; i<d(); ++i)
      e_p(i) = coord(p, i);

    p_bar = m * e_p;

    E_Vector_D p_on_S;
    for(std::size_t i=0; i<d(); ++i)
      p_on_S(i) = p_bar(i);

    int ind = d();
    for(std::size_t i=0; i<d(); ++i)
    {
      for(std::size_t j=i; j<d(); ++j)
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
  WPoint_d to_T(const WPoint_D& p_on_S) const
  {
    if(!is_tsb_init)
    {
      compute_tangent_space_basis();
      is_tsb_init = true;
    }

    typename KD::Scalar_product_d inner_pdct = KD().scalar_product_d_object();
    typename KD::Difference_of_points_d diff_points = KD().difference_of_points_d_object();

    Vector_D v = diff_points(p_on_S.point(), m_center_Q);

    // Ambiant-space coords of the projected point
    std::vector<FT> coords;
    coords.reserve(d());
    std::vector<FT> p_proj(m_center_Q.cartesian_begin(), m_center_Q.cartesian_end());

    for(std::size_t i=0 ;i<d() ;++i)
    {
      // Compute the inner product p * ts[i]
      FT coord = inner_pdct(v, m_tsb[i]);
      coords.push_back(coord);

      // p_proj += coord * v;
      for(std::size_t j=0 ; j<D(); ++j)
        p_proj[j] += coord * m_tsb[i][j];
    }

    Point_D projected_pt(D(), p_proj.begin(), p_proj.end());
    typename KD::Squared_distance_d sqdist = KD().squared_distance_d_object();

//    std::cout << "v: " << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << " " << v[4] << std::endl;
//    std::cout << "p_proj: " << p_proj[0] << " " << p_proj[1] << " " << p_proj[2] << " " << p_proj[3] << " " << p_proj[4] << std::endl;
//    std::cout << "coords: " << coords.size() << " " << coords[0] << " " << coords[1] << std::endl;
//    std::cout << "ppt: " << projected_pt[0] << " " << projected_pt[1] << " ";
//    std::cout << projected_pt[2] << " " << projected_pt[3] << " " << projected_pt[4] << std::endl;
//    std::cout << "w: " << p_on_S.weight() << std::endl;
//    std::cout << "sqdist: " << sqdist(p_on_S.point(), projected_pt) << std::endl;

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
    const std::size_t dp1 = dDim::value+1;
    int cids[dp1], dids[dp1];
    for(std::size_t i=0; i<dp1; i++)
      cids[i] = a[i];

    Full_cell_handle_iterator ci = finite_incident_full_cells_begin();
    Full_cell_handle_iterator cend = finite_incident_full_cells_end();
    for(; ci!=cend; ci++)
    {
      if((*ci)->maximal_dimension() < d())
        continue;
      for(std::size_t i=0; i<dp1; i++)
        dids[i] = (*ci)->vertex(i)->data();
      if(is_same_ids<dp1>(cids, dids))
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
    for(std::size_t i=0; i<dp1; i++)
      cids[i] = fch->vertex(i)->data();

    Full_cell_handle_iterator ci = finite_incident_full_cells_begin();
    Full_cell_handle_iterator cend = finite_incident_full_cells_end();
    for(; ci!=cend; ci++)
    {
      if((*ci)->maximal_dimension() < d())
        continue;
      for(std::size_t i=0; i<dp1; i++)
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
  // We take "p_on_S" and not "p" because "this" cannot transform p to p_on_S (it
  // doesn't know p's metric and calling to_T(p) would wrongly use "this"'s metric
  // instead of p's !)
  Vertex_handle insert_if_in_star(const WPoint_D& p_on_S,
                                  const Index& index_)
  {
    Vertex_handle vh = Base::insert_if_in_star(to_T(p_on_S), m_center_v);
    if(vh != Vertex_handle())
      vh->data() = index_;
    invalidate_cache();
    return vh;
  }

  Vertex_handle insert_to_star(const WPoint_D& p_on_S,
                               const Index& index,
                               const bool conditional = false)
  {
    Vertex_handle vh;
    if(conditional)
      vh = insert_if_in_star(p_on_S, index);
    else
      vh = Base::insert(to_T(p_on_S), m_center_v);

    if(vh != Vertex_handle())
      vh->data() = index;

    invalidate_cache();
#ifdef CGAL_ANISO_DEBUG
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "enter insert to star" << std::endl;
    std::cout << "insert " << index << " in " << this->index() << std::endl;
    WPoint_d tp = to_T(p_on_S);
    std::cout << "pt transformed to star's tangent plane: " << tp.point()[0] << " " << tp.point()[1] << " ";
    std::cout << tp.weight() << std::endl;
    std::cout << "finite cells: " << this->number_of_finite_full_cells() << std::endl;
    std::cout << "dim rt: " << this->current_dimension() << std::endl;
    this->is_valid(true);
#endif
    std::cout << "exit insert to star" << std::endl;
    return vh;
  }

  // Refinement ----------------------------------------------------------------
  FT norm(const Point_D& p) const // squared norm, actually
  {
    FT n = 0.;
    for(std::size_t i=0; i<D(); ++i)
      n += p[i]*p[i];
    return n;
  }

  bool is_support_intersecting(const std::vector<Star_handle>& cell,
                               Point_D& sol) const
  {
    // intersection with the tangent plane at the corresponding point
    // on the paraboloid of a point of a simplex

    // system is given by
    //   d(p,p0) - w0 = d(p,pi) - wi
    //   p \in T_{P_0}

    E_Matrix_D m = E_Matrix_D::Zero();
    E_Vector_D l;
    WPoint_D p0 = cell[0]->m_center_S;
    FT d0 = norm(p0.point());
    Point_D p0_on_Q = cell[0]->m_center_Q;

//    std::cout << "w0: " << cell[0]->index() << " || ";
//    std::cout << p0.point()[0] << " " << p0.point()[1] << " ";
//    std::cout << p0.point()[2] << " " << p0.point()[3] << " ";
//    std::cout << p0.point()[4] << " w " << p0.weight() << std::endl;

//    std::cout << "base: ";
//    std::cout << p0_on_Q[0] << " " << p0_on_Q[1] << std::endl;

    //first part of the matrix (lines 0-d) (dual equation)
    for(std::size_t i=0; i<d(); ++i)
    {
      WPoint_D pi = cell[i+1]->m_center_S;
//      std::cout << "wi: " << cell[i+1]->index() << " || ";
//      std::cout << pi.point()[0] << " " << pi.point()[1] << " ";
//      std::cout << pi.point()[2] << " " << pi.point()[3] << " ";
//      std::cout << pi.point()[4] << " w " << pi.weight() << std::endl;

      FT di = norm(pi.point());
      l(i) = d0 - p0.weight() - di + pi.weight();
      for(std::size_t j=0; j<D(); ++j)
        m(i,j) = 2*(p0.point()[j] - pi.point()[j]);
    }

    //second part of the matrix (lines d-D) (tangent plane equation)
    for(std::size_t i=d(); i<D(); ++i)
    {
      l(i) = p0_on_Q[i];
      m(i,i) = -1;
    }

    int pos = d();
    for(std::size_t i=0; i<d(); ++i)
    {
      for(std::size_t j=i; j<d(); ++j)
      {
        if(i == j)
          m(pos, i) = 2*(cell[0]->m_center)[i];
        else
        {
          m(pos, i) = (cell[0]->m_center)[j];
          m(pos, j) = (cell[0]->m_center)[i];
        }
        pos++;
      }
    }

//    std::cout << "m @ supportintersection:" << std::endl << m << std::endl;
//    std::cout << "LHS " << l.transpose() << std::endl;
    FT det = m.determinant();
//    std::cout << "The determinant of m is " << det << std::endl;

    if(std::abs(det) > 1e-10)
    {
      E_Matrix_D mm1 = m.inverse();
//      std::cout << "The inverse of m is:\n" << mm1 << std::endl;
      E_Vector_D vsol = mm1*l;
      sol = KD().construct_point_d_object()(D(), vsol.data(), vsol.data()+D());
//      std::cout << "sol: " << vsol.transpose() << std::endl;

//      for(int i=0; i<(d()+1); ++i)
//      {
//        const WPoint_D& wi = cell[i]->m_center_S;
//        FT di = KD().squared_distance_d_object()(sol, wi.point()) - wi.weight();
//        std::cout << "dist² sol " << cell[i]->index() << " : "  << di << std::endl;
//      }
      return true;
    }
    else
    {
      std::cout << "null det" << std::endl;
      return false;
    }
  }

  struct H // hyperplanes describing the dual and error computation
  {
    const std::vector<Star_handle>& cell;

    const int d() const { return dDim::value; }
    const int D() const { return DDim::value; }

    E_Vector_d operator()(const E_Vector_d& s) const
    {
      Point_d sp = Kd().construct_point_d_object()(d(), s.data(), s.data()+d());
      Point_D s_on_Q = cell[0]->to_Q(sp);

//      std::cout << "H() " << s.transpose() << std::endl;
//      std::cout << "sonQ " << s_on_Q[0] << " " << s_on_Q[1] << " " << s_on_Q[2] << " " << s_on_Q[3] << " " << s_on_Q[4] << std::endl;

      E_Vector_d ret = E_Vector_d::Zero();
      E_Vector_d F = E_Vector_d::Zero();
      E_Vector_d X = E_Vector_d::Zero();

      const WPoint_D& wp0 = cell[0]->m_center_S;
      for(int i=0; i<d(); ++i)
      {
        const WPoint_D& wi = cell[i+1]->m_center_S;
        for(int j=0; j<D(); ++j)
        {
          ret(i) += 2*(wi.point()[j]-wp0.point()[j])*s_on_Q[j];
          X(i) += 2*(wi.point()[j]-wp0.point()[j])*s_on_Q[j];
          ret(i) += wp0.point()[j]*wp0.point()[j] - wi.point()[j]*wi.point()[j];
          F(i) += wp0.point()[j]*wp0.point()[j] - wi.point()[j]*wi.point()[j];
        }
        ret(i) += wi.weight() - wp0.weight();
        F(i) += wi.weight() - wp0.weight();
      }
//      std::cout << "X: " << X.transpose() << " F: " << F.transpose() << std::endl;
      return ret;
    }

    H(const std::vector<Star_handle>& cell_)
      : cell(cell_)
    { }
  };

  struct J // Jacobian
  {
    const std::vector<Star_handle>& cell;

    const int d() const { return dDim::value; }
    const int D() const { return DDim::value; }

    // the matrix is the jacobian of the system :
    //  d(p,p0) - w0 = d(p,pi) - wi
    //  p \in Q
    // that is : J(i,j) = \frac{\partial \sum_k 2*(p_k-pi_k)x_j}{\partial p_j}
    // where the x_j are x x² (if d=1), x y x² xy y² (if d=2), etc.

    E_Matrix_d operator()(const E_Vector_d& s) const
    {
      E_Matrix_d m;
      WPoint_D wp0 = cell[0]->m_center_S;
//      std::cout << "w0: " << cell[0]->index() << " || ";
//      std::cout << wp0.point()[0] << " " << wp0.point()[1] << " ";
//      std::cout << wp0.point()[2] << " " << wp0.point()[3] << " ";
//      std::cout << wp0.point()[4] << std::endl;

      //the easy coefficients : those in the part x y z of \tilde{P}
      for(int i=0; i<d(); ++i)
      {
        WPoint_D wpi = cell[i+1]->m_center_S;

//        std::cout << "wi: " << cell[i+1]->index() << " || ";
//        std::cout << wpi.point()[0] << " " << wpi.point()[1] << " ";
//        std::cout << wpi.point()[2] << " " << wpi.point()[3] << " ";
//        std::cout << wpi.point()[4] << std::endl;

        for(int j=0; j<d(); ++j) // x y z...
          m(i,j) = 2.*(wpi.point()[j] - wp0.point()[j]);
      }

      // and now for the hard part: x² xy xz y² yz z² etc...
      for(int l=0; l<d(); ++l) // lines of m (the function)
      {
        const WPoint_D wpl = cell[l+1]->m_center_S;
        for(int c=0; c<d(); ++c) // columns of m (the derivative)
        {
          //ij are the indices of the Q matrix
          //max is the number of entries at line i (max=d-i)
          //count is the number of entries visited (lexico order)
          int i=0, j=0, max=d(), count=0;
          while(max > 0)
          {
            assert(i<d() && j<d());
            FT coeff = 0.;
            if(i == c) // the line Q(i,.) : x_i*x_?
            {
              if(j == c) // the coeff Q(i,j) is x_i*x_i
                coeff = 2*s(i);
              else
                coeff = s(j); // the coeff Q(i,j) is x_i*x_j
            }
            else if(j == c) // the coeff Q(i,j) is x_j*x_? ('?' can't be i)
              coeff = s(i);

            assert(d()+count < D());
            // "+d" since we've added the easy part already (the d first)
            m(l, c) += 2.*(wpl.point()[d()+count]-wp0.point()[d()+count])*coeff;

            j++; count++; // switch from x_i*x_j to x_i*x_{j+1}
            if(j >= max) // reached the end of a line of Q
            {
              max--; i++; j=i; // switch from x_i*x_d to x_{i+1}*x_{i+1}
            }
          }
        }
      }
      return m;
    }

    J(const std::vector<Star_handle>& cell_)
      : cell(cell_)
    { }
  };

  bool is_inside_power_cell(const Point_D& sol,
                            const std::vector<Star_handle>& cell,
                            const std::vector<Star_handle>& all_stars) const
  {
    std::cout << "Check if the point is in the Pcell of the simplex" << std::endl;

    typename KD::Squared_distance_d csd = KD().squared_distance_d_object();

    const WPoint_D& wp0 = cell[0]->m_center_S;
    FT sq_d = csd(sol, wp0.point()) - wp0.weight();
    for(std::size_t i=0; i<d(); i++)
    {
      const WPoint_D& wpi = cell[i+1]->m_center_S;
      FT sq_dd = csd(sol, wpi.point()) - wpi.weight();
//      std::cout << "DIST: " << sq_d << " " << sq_dd << std::endl;
      if(std::abs(sq_d-sq_dd)/sq_d > 1e-5) // sq_d might be null, I guess... fixme
      {
        std::cout.precision(20);
        std::cout << "WARNING: dds don't agree (sol isn't in the dual) " << sq_d << " " << sq_dd << std::endl;
      }
    }

    typename std::vector<Star_handle>::const_iterator sit = all_stars.begin();
    typename std::vector<Star_handle>::const_iterator send = all_stars.end();
    for(; sit!=send; ++sit)
    {
      Star_handle si = *sit;
      if(has(cell, si->index())) // vertex belongs to the simplex 'cell'
        continue;

      WPoint_D wpi = si->m_center_S;
      FT sq_dv = csd(sol, wpi.point()) - wpi.weight();
      if(std::abs(sq_dv-sq_d)/sq_d>1e-5 && sq_dv < sq_d) //fixme /sq_d
      {
        // not in the power cell of Pi
        std::cout << "closer : " << si->index() << " " << sq_dv;
        std::cout << " ref: " << sq_d << std::endl;
        return false;
      }
    }
    std::cout << "inside the power cells!" << std::endl;
    return true;
  }

  void tinker_jacobian(E_Matrix_d& m, const E_Vector_d& h) const
  {
    //we add a diagonal matrix of the form delta*f_i(x) to J to avoid |J| ~ 0
    E_Matrix_d m2 = E_Matrix_d::Zero();
    for(std::size_t i=0; i<d(); ++i)
    {
      m2(i,i) = m(i,i)/std::abs(m(i,i)); // sign
      m2(i,i) *= h(i);
    }

    m += m2;
  }

  bool newton(Point_D& sol_on_Q,
              const std::vector<Star_handle>& cell,
              const std::vector<Star_handle>& all_stars) const
  {
    // we are looking for the intersection of the dual of to_S(fch) with Q
    // we start with the approximation intersection of dual(to_S(fch)) with T
    if(!is_support_intersecting(cell, sol_on_Q)) // get the approximation
      return false;

    std::cout << "approximation: " << sol_on_Q[0] << " " << sol_on_Q[1] << std::endl;
    E_Vector_d vsol;
    for(std::size_t i=0; i<d(); ++i)
      vsol(i) = sol_on_Q[i];  // ugly 'projection' on the paraboloid
                              // of the intersection with the tangent plane

    FT err = 1e-4;
    H h(cell); // compute the error
    J j(cell); // compute the jacobian

    int max = 1e3, count = 0;
    E_Vector_d h_err = h(vsol);
    E_Matrix_d jac = j(vsol);

    tinker_jacobian(jac, h_err);

    E_Vector_d previous_vsol = vsol;
    FT step = 1e30;

    while(h_err.norm() > err && ++count < max)
    {
      if(std::abs(jac.determinant())<1e-10)
      {
        std::cout << "jacobian det is null" << std::endl;
        assert(0);
        return false; // aborting the refinement of this cell...
      }

      previous_vsol = vsol;

      vsol = vsol - jac.inverse()*h(vsol);
      jac = j(vsol);
      h_err = h(vsol);
//      std::cout << "J det (no tinkerino): " << jac.determinant() << " ";
      tinker_jacobian(jac, h_err);

      step = (vsol-previous_vsol).norm();

//      std::cout << "jac: " << std::endl << jac << std::endl;
//      std::cout << "J det: " << jac.determinant() << " ";
//      std::cout << "vsol: " << vsol.transpose() << " ";
//      std::cout << "h: " << h_err.transpose() << " (" << h_err.norm() << ")" << std::endl;

//      std::cout << std::min(1e3, jac.determinant()) << " ";
//      std::cout << std::min(1e3, h_err.norm()) << std::endl;
    }

    if(count==max)
    {
      std::cout << "didn't converge in " << max << " iterations. ";
      std::cout << "h: " << h_err.transpose() << std::endl;
//      assert(0);
      return false;
    }
    else
      std::cout << "converged to: " << vsol.transpose() << std::endl;

    Point_d vsolp = Kd().construct_point_d_object()(d(), vsol.data(), vsol.data()+d());
    sol_on_Q  = to_Q(vsolp);
    std::cout << "final solution : " << vsol.transpose() << " in " << count << " iterations" << std::endl;
    std::cout << "final jac/h: " << jac.determinant() << " " << h_err.norm() << std::endl;

    return is_inside_power_cell(sol_on_Q, cell, all_stars);
  }

  bool compute_dual_intersection(Point_D& sol_on_Q,
                                 const std::vector<Star_handle>& cell,
                                 const std::vector<Star_handle>& all_stars) const
  {
    // Newton's method to (try to) find precisely the root on Q
    if(newton(sol_on_Q, cell, all_stars)) // this checks if we're in the power cell too
     return true;
    std::cout << "dual does not intersect...? " << std::endl;
    return false;
  }

  bool compute_dual(const Full_cell_handle fch,
                    const std::vector<Star_handle>& all_stars)
  {
    std::vector<Star_handle> cell;
    for(std::size_t i=0; i<=d(); ++i)
      cell.push_back(all_stars[fch->vertex(i)->data()]);

    Point_D ponQ;
    if(compute_dual_intersection(ponQ, cell, all_stars))
    {
      Point_d p = from_Q(ponQ);
      fch->data().first = Kd().construct_point_d_object()(d(), p.begin(), p.end());
      fch->data().second = true;
      return true;
    }
    fch->data().second = false;
    return false;
  }

/* COMMENTED CODE DOES NOT WORK : it is pointless to compute the intersection
 * of the dual with the tangent plane since we can't actually project it back
 * in the global space afterwards. A possible hack is to project the intersection
 * with the tangent plane on the paraboloid and then to use this as a global
 * point but it's ugly.
 *
 * If you one day desire to use the function below anyway, you also need to
 * recreate the metric_helper.hpp that wasn't commited and got lost (or just
 * the logxexp matrix interpolation).

  bool compute_dual_from_tangent_plane(const Full_cell_handle fch,
                                       const std::vector<Star_handle>& all_stars)
  {
    std::cout << "compute dual from tangent plane" << std::endl;

    //interpolation (if you change it here, change it in the new metric computation too!)
    std::vector<std::pair<E_Matrix_d, FT> > metric_transfs;
    for(std::size_t i=0; i<d()+1; ++i)
    {
      const E_Matrix_d& mi = all_stars[fch->vertex(i)->data()]->metric().get_mat();
      std::cout << "mi : " << fch->vertex(i)->data() << std::endl << mi << std::endl;
      metric_transfs.push_back(std::make_pair(mi, 1./(d()+1)));
    }
    E_Matrix_d m = logexp_interpolate<Kd, E_Matrix_d>(metric_transfs);

    //compute the weighted circumcenter on the tangent plane
    typename Traits::Power_center_d pc = m_traits->power_center_d_object();
    typename Traits::Compute_coordinate_d coord = m_traits->compute_coordinate_d_object();
    typename KD::Compute_coordinate_d coordD = KD().compute_coordinate_d_object();

    std::vector<WPoint_d> pts;
    for(int i=0; i<d()+1; ++i)
      pts.push_back(fch->vertex(i)->point());

    WPoint_d c = pc(pts.begin(), pts.end());

    //solve linear system to get the ambiant space to find p in R^d
    //we solve ((hat p).n_i) = c_i (finding a point in R^D whose coordinates on the tangent planes are c_i)
    E_Matrix_d A = E_Matrix_d::Zero();
    E_Vector_d B = E_Vector_d::Zero();

    //fill A & B
    std::cout << "fill A&B: " << std::endl;
    for(std::size_t i=0; i<d(); ++i)
    {
      for(std::size_t j=0; j<d(); ++j) // filling A(i,j) and B(i)
      {
        for(std::size_t k=0; k<d(); ++k) // the M*p coefficients of hat p are in A
          A(i,j) += m(k,j)*m_tsb[i][k];

        B(i) = coord(c, i);
        int pos = d();
        for(std::size_t k=0; k<d();++k) // the Q coefficients of hat p are removed from B
        {
          for(std::size_t l=k; l<d(); ++l)
          {
            if(k==l)
              B(i) += 0.5*m(k,k)*m_tsb[i][pos++];
            else
              B(i) += m(k,l)*m_tsb[i][pos++];
          }
        }
      }

      for(std::size_t j=0; j<D(); ++j)
        B(i) += coordD(m_center_Q, j) * m_tsb[i][j];
    }

    E_Vector_d x = A.colPivHouseholderQr().solve(B);

#ifdef ANISO_TC_DEBUG
    typename Traits::Power_distance_d pd = m_traits->power_distance_d_object();

    std::cout << "interpolated metric is : " << std::endl << m << std::endl;
    std::cout << "computed c " << coord(c,0) << " " << coord(c,1) << " " << c.weight() << std::endl;
    std::cout << "distance check for c: ";
    for(std::size_t i=0; i<d()+1; ++i)
      std::cout << pd(c, fch->vertex(i)->point()) << " ";
    std::cout << std::endl << "x: " << x.transpose() << std::endl;

    E_Vector_d tilde_p = m * x;
    E_Vector_D hat_p;
    for(std::size_t i=0; i<d(); ++i)
      hat_p(i) = tilde_p(i);

    int ind = d();
    for(std::size_t i=0; i<d(); ++i)
    {
      for(std::size_t j=i; j<d(); ++j)
      {
        if(j==i)
          hat_p(ind++) = -0.5*m(i,i);
        else
          hat_p(ind++) = -m(i,j);
      }
    }
    FT n = hat_p.norm();
    FT w = n * n - x.transpose() * tilde_p;
    WPoint_D x_to_S = KD().construct_weighted_point_d_object()(
                        KD().construct_point_d_object()(D(), hat_p.data(), hat_p.data() + D()),
                        w);

    std::cout << "computed x_to_S check at: " << hat_p.transpose() << " " << w << std::endl;
    std::cout << "one of the following has to be equal to c: " << std::endl;
    for(std::size_t i=0; i<d()+1; ++i)
    {
      std::cout << fch->vertex(i)->data() << " || ";
      for(std::size_t j=0; j<d(); ++j)
        std::cout << coord(all_stars[fch->vertex(i)->data()]->to_T(x_to_S), j) << " ";
      std::cout << std::endl;
      std::cout << "local coordinates of fch in the star above : ";
      for(std::size_t j=0; j<d(); ++j)
        std::cout << coord(fch->vertex(i)->point(),j) << " ";
      std::cout << fch->vertex(i)->point().weight() << std::endl;
    }
#endif

    fch->data().first = Kd().construct_point_d_object()(d(),x.data(), x.data()+d());
    fch->data().second = true;
    return true;
  }
*/

  bool is_inside(const Point_d& p) // fixme this should check if we're in the domain
  {
    return (p[0] < -5 || p[0] > 5 || p[1] < -5 || p[1] > 5);
  }

  bool is_inside(const Full_cell_handle fch,
                 const std::vector<Star_handle>& all_stars)
  {
    std::cout << "is inside @ " << index() << " (";
    std::cout << fch->vertex(0)->data() << " ";
    std::cout << fch->vertex(1)->data() << " ";
    std::cout << fch->vertex(2)->data() << ") " << std::endl;

    if(!compute_dual(fch, all_stars) || !is_inside(fch->data().first))
    {
      fch->data().second = false;
      return false;
    }
    return true;
  }

  // Maintenance ---------------------------------------------------------------
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
      Base::insert(backup[i].first, backup[i].second, false);

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

  // Misc stuff ----------------------------------------------------------------
  void output_underlying_rt(std::ofstream& out) const
  {
    out.precision(17);
    out << d() << " " << Base::number_of_vertices() << std::endl;
    typename Base::Finite_vertex_const_iterator it = Base::finite_vertices_begin();
    typename Base::Finite_vertex_const_iterator end = Base::finite_vertices_end();
    for(; it!=end; ++it)
    {
      out << it->point().point()[0] << " ";
      out << it->point().point()[1] << " ";
      out << it->point().weight() << std::endl;
    }
  }

  void output_full_underlying_rt(std::ofstream& out,
                                 const std::vector<Star_handle>& all_stars) const
  {
    out.precision(17);
    out << d() << " " << all_stars.size() << std::endl;
    typename std::vector<Star_handle>::const_iterator it = all_stars.begin();
    typename std::vector<Star_handle>::const_iterator end = all_stars.end();
    for(; it!=end; ++it)
    {
      WPoint_d p = to_T((*it)->m_center_S);
      out << p.point()[0] << " " << p.point()[1] << " ";
      out << p.weight() << std::endl;
    }
  }

public:
  Tangent_star(const Point_d& center_point,
               const Index& index_,
               //const Criteria* criteria_, // todo
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
      m_center_T(to_T(m_center_S)),
      m_is_cache_dirty(true),
      m_finite_adjacent_vertices_cache(),
      m_incident_full_cells_cache(),
      m_finite_incident_full_cells_cache()
  {
    m_center_v = Base::insert(m_center_T);
    m_center_v->data() = index_;
    this->infinite_vertex()->data() = index_of_infinite_vertex;

#ifdef ANISO_TC_DEBUG
    typename Kd::Compute_coordinate_d coord_d = Kd().compute_coordinate_d_object();
    typename KD::Compute_coordinate_d coord_D = KD().compute_coordinate_d_object();

    std::cout << "New star constructor " << index_ << std::endl;
    std::cout << "center_p: " << coord_d(center_point, 0) << " " << coord_d(center_point, 1) << std::endl;
    std::cout << "transformation : " << std::endl << m_metric.get_transformation() << std::endl;

    std::cout << "center_S: ";
    for(int i=0; i<5; ++i)
      std::cout << coord_D(m_center_S.point(), i) << " ";
    std::cout << m_center_S.weight() << std::endl;

    std::cout << "center_Q: ";
    for(int i=0; i<5; ++i)
      std::cout << coord_D(m_center_Q, i) << " ";
    std::cout << std::endl;
#endif
  }

};

} // Anisotropic_mesh_TC
} // CGAL

#endif //CGAL_ANISOTROPIC_MESH_TC_STAR

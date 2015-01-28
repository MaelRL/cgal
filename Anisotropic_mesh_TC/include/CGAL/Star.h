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
    Vertex_handle vh = Base::insert_if_in_star(to_T(p), m_center_v);
    if(vh != Vertex_handle())
      vh->data() = index_;
    invalidate_cache();
    return vh;
  }

  Vertex_handle insert_to_star(const Point_d& p,
                               const Index& index,
                               const bool conditional = false)
  {
    Vertex_handle vh;
    if(conditional)
      vh = insert_if_in_star(p, index);
    else
      vh = Base::insert(to_T(p), m_center_v);
    if(vh != Vertex_handle())
      vh->data() = index;

    invalidate_cache();
    std::cout << "INSERT : " << index << " IN STAR " << m_center_v->data() << std::endl;
    std::cout << "finite cells: " << this->number_of_finite_full_cells() << std::endl;
    std::cout << "dim rt: " << this->current_dimension() << std::endl;
    this->is_valid(true);

    return vh;
  }

  // Refinement ----------------------------------------------------------------
  FT norm(const Point_D& p) const // squared norm, actually
  {
    FT n = 0.;
    for(int i=0; i<D(); ++i)
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
    for(int i=0; i<d(); ++i)
    {
      WPoint_D pi = cell[i+1]->m_center_S;
//      std::cout << "wi: " << cell[i+1]->index() << " || ";
//      std::cout << pi.point()[0] << " " << pi.point()[1] << " ";
//      std::cout << pi.point()[2] << " " << pi.point()[3] << " ";
//      std::cout << pi.point()[4] << " w " << pi.weight() << std::endl;

      FT di = norm(pi.point());
      l(i) = d0 - p0.weight() - di + pi.weight();
      for(int j=0; j<D(); ++j)
        m(i,j) = 2*(p0.point()[j] - pi.point()[j]);
    }

    //second part of the matrix (lines d-D) (tangent plane equation)
    for(int i=d(); i<D(); ++i)
    {
      l(i) = p0_on_Q[i];
      m(i,i) = -1;
    }

    int pos = d();
    for(int i=0; i<d(); ++i)
    {
      for(int j=i; j<d(); ++j)
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
    for(int i=0; i<d(); i++)
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
    for(int i=0; i<d(); ++i)
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
    if(!is_support_intersecting(cell, sol_on_Q)) // get the approximation
      return false;

    std::cout << "approximation: " << sol_on_Q[0] << " " << sol_on_Q[1] << std::endl;
    E_Vector_d vsol;
    for(int i=0; i<d(); ++i)
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
      std::cout << "didn't converge in " << max << " iterations.";
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
    // we are looking for the intersection of the dual of to_S(fch) with Q
    // we start with the approximation intersection of dual(to_S(fch)) with T
    if(newton(sol_on_Q, cell, all_stars))
    {
      // Newton's method to (try to) find precisely the root on Q
//      if(is_intersection_point_in_power_cell(sol_on_Q))
          return true;
    }
    std::cout << "dual does not intersect...? " << std::endl;
    return false;
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

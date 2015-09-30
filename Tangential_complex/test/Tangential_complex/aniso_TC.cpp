#include <CGAL/assertions_behaviour.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>

#include <CGAL/Tangential_complex.h>
#include "../../test/Tangential_complex/test_utilities.h"

#include <CGAL/Metric_field.h> // Anisotropic metrics
#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <cstdlib>
#include <fstream>
#include <math.h>

const int k = 2; // intrinsic
const int d = 5; // ambiant

typedef CGAL::Epick_d< CGAL::Dimension_tag<k> >                 Kk;
typedef CGAL::Epick_d< CGAL::Dimension_tag<d> >                 Kd;
typedef Kk::FT                                                  FT;
typedef Kk::Point_d                                             Point_k;
typedef Kd::Point_d                                             Point_d;
typedef Kd::Weighted_point_d                                    Weighted_point_d;
typedef Kd::Vector_d                                            Vector_d;
typedef CGAL::Tangential_complex<Kd,
                                 CGAL::Dimension_tag<k>,
                                 CGAL::Parallel_tag>            TC;

typedef CGAL::Tangential_complex_::Basis<Kd>                    Basis;
typedef CGAL::Anisotropic_mesh_TC::Metric_base<Kk>              Metric;
typedef CGAL::Anisotropic_mesh_TC::Metric_field<Kk>             Metric_field;

typedef typename Metric::E_Matrix                               E_Matrix_k;
typedef typename Metric::E_Vector                               E_Vector_k;
typedef Eigen::Matrix<FT, d, 1>                                 E_Vector_d;
typedef Eigen::Matrix<FT, d, k>                                 E_Matrix_dk;

typedef CGAL::Anisotropic_mesh_TC::Euclidean_metric_field<Kk>   Euclidean_mf;
typedef CGAL::Anisotropic_mesh_TC::Custom_metric_field<Kk>      Custom_mf;

void read_points(std::vector<Point_k>& points,
                 const Kk kerk = Kk())
{
  std::ifstream in("../../../../Anisotropic_mesh_2/examples/Anisotropic_mesh_2/build/bambimboum.mesh");
//  std::ifstream in("aniso_regular.mesh");
  std::string word;
  int useless, nv, dd;
  FT x;

  in >> word >> useless; //MeshVersionFormatted i
  in >> word >> dd; //Dimension d
  in >> word >> nv;

  assert(dd == Kk::Dimension::value);

  for(int i=0; i<nv; ++i)
  {
    std::vector<FT> ids;
    for(int j=0; j<k; ++j)
    {
      in >> x;
      ids.push_back(x);
    }
    in >> useless;
    points.push_back(kerk.construct_point_d_object()(ids.begin(), ids.begin() + k));

    if(points.size() == 200)
      break;
  }

  std::cout << points.size() << " points" << std::endl;
}

// compute the corresponding point on the paraboloid (these points will be the
// points at which we compute tangent planes)
Point_d to_Q(const Point_k& p)
{
  typename Kk::Compute_coordinate_d k_coord = Kk().compute_coordinate_d_object();

  E_Vector_d p_on_Q;
  for(int i=0; i<k; ++i)
    p_on_Q(i) = k_coord(p, i);

  int ind = k;
  for(int i=0; i<k; ++i)
    for(int j=i; j<k; ++j)
      p_on_Q(ind++) = k_coord(p,i) * k_coord(p,j);

  return Kd().construct_point_d_object()(d, p_on_Q.data(), p_on_Q.data() + d);
}

// compute the corresponding point on the metric manifold (these points will be
// the seeds of the global power diagram)
Weighted_point_d to_S(const Point_k& p,
                      const Metric& met,
                      const Kd kerd = Kd())
{
  const E_Matrix_k m = met.get_mat();
  E_Vector_k e_p, p_bar;

  typename Kk::Compute_coordinate_d k_coord = Kk().compute_coordinate_d_object();
  for(int i=0; i<k; ++i)
    e_p(i) = k_coord(p, i);

  p_bar = m * e_p;

  E_Vector_d p_on_S;
  for(int i=0; i<k; ++i)
    p_on_S(i) = p_bar(i);

  int ind = k;
  for(int i=0; i<k; ++i)
  {
    for(int j=i; j<k; ++j)
    {
      if(j==i)
        p_on_S(ind++) = -0.5*m(i,i);
      else
        p_on_S(ind++) = -m(i,j);
    }
  }

  FT n = p_on_S.norm();
  FT w = n * n - e_p.transpose() * p_bar;

  return kerd.construct_weighted_point_d_object()
      (kerd.construct_point_d_object()(d, p_on_S.data(),p_on_S.data()+d), w);
}

void compute_points_in_ambiant_space(std::vector<Point_k> const & points_k,
                                     Metric_field const * const mf,
                                     std::vector<Point_d>& points_d,
                                     std::vector<FT>& weights,
                                     Kd const kerd = Kd())
{
  typename Kd::Point_drop_weight_d k_drop_w = kerd.point_drop_weight_d_object();
  typename Kd::Point_weight_d k_point_weight = kerd.point_weight_d_object();

  for(std::size_t i=0; i<points_k.size(); ++i)
  {
    Point_k const & p = points_k[i];
    Metric met = mf->compute_metric(p);
    Weighted_point_d wp = to_S(p, met, kerd);
    points_d.push_back(k_drop_w(wp));
    weights.push_back(k_point_weight(wp));
  }
}

void compute_and_set_tangent_planes(TC& tc,
                                    const std::vector<Point_k>& points_k,
                                    Kd const kerd = Kd())
{
  typedef TC::TS_container TS_container;
  typedef TC::OS_container OS_container;

  typename Kd::Compute_coordinate_d coord = kerd.compute_coordinate_d_object();
  typename Kd::Construct_vector_d constr_vec = kerd.construct_vector_d_object();
  TS_container tsc;
  OS_container osc;
  std::size_t n = tc.number_of_vertices();
  CGAL_assertion(n == points_k.size());

  for(std::size_t c=0; c<n; ++c)
  {
// origin of the basis is the point moved the paraboloid
    Point_d origin = to_Q(points_k[c]);

// compute the tsc; the vectors are given by the partial derivatives
// of the parametrization of the paraboloid
    Basis ts(origin);

    // filling all the vectors at the same time through a (d,k)-matrix
    E_Matrix_dk ts_m = E_Matrix_dk::Zero();
    // 'first' part: x y z etc. derivates
    for(int i=0; i<k; ++i)
      for(int j=0; j<k; ++j)
        ts_m(i,j) = (i==j);

    // 'second' part: x² xy xz y² yz z² etc. derivatives
    int pos = k;
    for(int i=0; i<k; ++i)
    {
      for(int j=i; j<k; ++j)
      {
        if(i == j)
          ts_m(pos, i) = 2*coord(origin, i);
        else
        {
          ts_m(pos, i) = coord(origin, j);
          ts_m(pos, j) = coord(origin, i);
        }
        pos++;
      }
    }

    for(int i=0; i<k; ++i)
    {
      Vector_d v = constr_vec(d, ts_m.col(i).data(), ts_m.col(i).data() + d);
      ts.push_back(v);
    }
    tsc.push_back(CGAL::Tangential_complex_::compute_gram_schmidt_basis(ts, kerd));

#ifdef CGAL_ALPHA_TC
// compute the osc
    Basis os(origin);

    // todo
    osc.push_back(CGAL::Tangential_complex_::compute_gram_schmidt_basis(os, kerd));
#endif
  }

  tc.set_tangent_planes(tsc
#ifdef CGAL_ALPHA_TC
                      , osc
#endif
                        );
}

std::ostream& export_simplices_to_medit(const TC& tc,
                                        const TC::Simplicial_complex &complex,
                                        std::ostream & os,
                                        std::size_t &num_simplices,
                                        const std::vector<Point_k>& points_k)
{
  typedef TC::Simplicial_complex::Simplex                     Simplex;
  typedef TC::Simplicial_complex::Simplex_range               Simplex_range;

  num_simplices = 0;
  typename Simplex_range::const_iterator it_s = complex.simplex_range().begin();
  typename Simplex_range::const_iterator it_s_end = complex.simplex_range().end();
  // For each triangulation
  for ( ; it_s != it_s_end ; ++it_s)
  {
    Simplex c = *it_s;

    bool color_simplex = !tc.is_simplex_consistent(c);

    // Gather the triangles here
    typedef std::vector<Simplex> Triangles;
    Triangles triangles;

    std::size_t num_vertices = c.size();
    // Do not export smaller dimension simplices
    if (num_vertices < tc.m_intrinsic_dimension + 1)
      continue;

    if (num_vertices <= 3)
    {
      triangles.push_back(c);
    }
    else
    {
      // num_vertices >= 4: decompose the simplex in triangles
      std::vector<bool> booleans(num_vertices, false);
      std::fill(booleans.begin() + num_vertices - 3, booleans.end(), true);
      do
      {
        std::set<std::size_t> triangle;
        std::set<std::size_t>::iterator it = c.begin();
        for (int i = 0; it != c.end() ; ++i, ++it)
        {
          if (booleans[i])
            triangle.insert(*it);
        }
        triangles.push_back(triangle);
      }
      while (std::next_permutation(booleans.begin(), booleans.end()));
    }

    // For each cell
    Triangles::const_iterator it_tri = triangles.begin();
    Triangles::const_iterator it_tri_end = triangles.end();
    for ( ; it_tri != it_tri_end ; ++it_tri)
    {
      // Don't export infinite cells
      if (*it_tri->rbegin() == std::numeric_limits<std::size_t>::max())
        continue;

//      bool avoid_triangle = false;
//      std::set<std::size_t>::const_iterator it_point_idx0 = it_tri->begin();
//      for ( ; it_point_idx0 != it_tri->end() ; ++it_point_idx0)
//      {
//        typename Kk::Squared_distance_d sqd = Kk().squared_distance_d_object();

//        std::vector<FT> null_coords(k,0.);
//        Point_k ori = Kk().construct_point_d_object()(k, null_coords.begin(), null_coords.end());
//        std::cout << "sqd is: " << sqd(points_k[*it_point_idx0], ori) << std::endl;
//        if(sqd(points_k[*it_point_idx0], ori) > 2)
//        {
//          avoid_triangle = true;
//          break;
//        }
//      }
//      if(avoid_triangle)
//        continue;

      std::set<std::size_t>::const_iterator it_point_idx = it_tri->begin();
      for ( ; it_point_idx != it_tri->end() ; ++it_point_idx)
        os << *it_point_idx+1 << " ";

      if (color_simplex)
        os << "1" << std::endl;
      else
        os << "2" << std::endl;

      ++num_simplices;
    }
  }

  return os;
}

// translate the connectivity of the stars back to the original points
void export_complex_in_origin_space(const TC& tc,
                                    const TC::Simplicial_complex& complex,
                                    const std::vector<Point_k>& points_k,
                                    std::ofstream& out)
{
  std::stringstream output;
  std::size_t num_simplices;
  typename Kk::Compute_coordinate_d k_coord = Kk().compute_coordinate_d_object();

  out << "MeshVersionFormatted 1" << std::endl;
  out << "Dimension " << tc.m_intrinsic_dimension << std::endl;
  out << "Vertices" << std::endl;
  out << tc.number_of_vertices() << std::endl;
  for(std::size_t i=0; i<points_k.size(); ++i)
  {
    const Point_k& p = points_k[i];
    out << k_coord(p, 0) << " " << k_coord(p, 1) << " 0" << std::endl;
  }

  export_simplices_to_medit(tc, complex, output, num_simplices, points_k);
  out << "Triangles" << std::endl;
  out << num_simplices << std::endl;
  out << output.str();
}

void export_paraboloid(const TC& tc,
                       const TC::Simplicial_complex& complex,
                       std::ofstream& out)
{
  std::stringstream output;
  std::size_t num_simplices;
  typename Kd::Compute_coordinate_d d_coord = Kd().compute_coordinate_d_object();

  for(std::size_t i=0; i<tc.number_of_vertices(); ++i)
  {
    const Point_d& p = tc.m_tangent_spaces[i].origin();
    output << d_coord(p, 2) << " " << d_coord(p, 3) << " " << d_coord(p, 4) << std::endl;
  }

  tc.export_simplices_to_off(complex, output, num_simplices);

  out << "OFF" << std::endl;
  out << tc.number_of_vertices() << " " << num_simplices << " 0" << std::endl;
  out << output.str();
}

void make_tc(std::vector<Point_k>& points_k,
             Metric_field const * const mf)
{
  Kd ker_d;

  std::vector<Point_d> points_d;
  std::vector<FT> weights;
  compute_points_in_ambiant_space(points_k, mf, points_d, weights, ker_d);

  TC tc(points_d.begin(), points_d.end(), 0./*sparsity*/, k /*intr dim*/, ker_d);
  tc.set_weights(weights);
  compute_and_set_tangent_planes(tc, points_k, ker_d);

  tc.compute_tangential_complex();

  TC::Simplicial_complex complex;
  int max_dim = tc.export_TC(complex, false);
  complex.display_stats();
  {
    std::ofstream off_stream("aniso_complex.off");
    tc.export_to_off(complex, off_stream);
  }

  for(int i=0; i<tc.number_of_vertices(); ++i)
  {
    std::ostringstream os;
    os << "star_" << i << ".mesh" << std::ends;
    std::ofstream out(os.str().c_str());
    TC::Simplicial_complex star;
    tc.export_TC(star, false, i);
    export_complex_in_origin_space(tc, star, points_k, out);
  }

  // Collapse
  complex.collapse(max_dim);
  {
    std::ofstream off_stream("aniso_after_collapse.off");
    export_complex_in_origin_space(tc, complex, points_k, off_stream);
  }

  std::size_t num_wrong_dim_simplices, num_wrong_number_of_cofaces;
  bool pure_manifold = complex.is_pure_manifold(k, false, 1,
                                                &num_wrong_dim_simplices,
                                                &num_wrong_number_of_cofaces);
  complex.display_stats();

  std::ofstream out("aniso_TC.mesh");
  export_complex_in_origin_space(tc, complex, points_k, out);

  std::ofstream outp("aniso_TC_paraboloid.off");
  export_paraboloid(tc, complex, outp);

  for(std::size_t i=0; i<tc.number_of_vertices(); ++i)
    tc.compute_projected_points(i);
}

int main()
{
  std::freopen("log_TC.txt", "w", stdout); // redirect std::cout to log.txt

  CGAL::default_random = CGAL::Random();
  std::vector<Point_k> points_k;

  //Custom_mf* mf = new Custom_mf();
  Euclidean_mf* mf = new Euclidean_mf();

  read_points(points_k);
  make_tc(points_k, mf);

  return 0;
}

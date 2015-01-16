#ifndef CGAL_ANISOTROPIC_MESH_TC_METRIC_H
#define CGAL_ANISOTROPIC_MESH_TC_METRIC_H

#include <Eigen/Dense>

#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_TC
{

template<typename Kd>
class Metric_base
{
public:
  typedef typename Kd::FT                                            FT;
  typedef typename Kd::Point_d                                       Point_d;
  typedef typename Kd::Vector_d                                      Vector_d;

  typedef typename Kd::Dimension                                     Dim;

  typedef typename Eigen::Matrix<FT, Dim::value, 1>                  E_Vector;
  typedef typename Eigen::Matrix<FT, Dim::value, Dim::value>         E_Matrix;

private:
  E_Matrix eigen_transformation, eigen_inverse_transformation;

  // no sense in dim d
  std::vector<FT> evals;
  std::vector<Vector_d> evecs;

public:
  const E_Matrix& get_transformation() const { return eigen_transformation; }
  const E_Matrix& get_inverse_transformation() const { return eigen_inverse_transformation; }
  E_Matrix get_mat() const { return eigen_transformation.transpose()*eigen_transformation;}

  const Vector_d& get_evec(const int i) const { return evecs[i]; }
  const FT& get_eval(const int i) const { return evals[i]; }

public:
  FT compute_distortion(const Metric_base &m) const
  {
    FT eigen_dis1 = (m.eigen_transformation * eigen_inverse_transformation).operatorNorm();
    FT eigen_dis2 = (eigen_transformation * m.eigen_inverse_transformation).operatorNorm();

#ifdef ANISO_DEBUG_METRIC
    FT max_dist = 1.001;
    if(max(eigen_dis1, eigen_dis2) > max_dist)
    {
      std::cout << "--------" << std::endl;
      std::cout << "--THIS-ET---" << std::endl;
      std::cout << eigen_transformation << std::endl;
      std::cout << "--THIS-EIT---" << std::endl;
      std::cout << eigen_inverse_transformation << std::endl;
      std::cout << "--M-ET---" << std::endl;
      std::cout << m.eigen_transformation << std::endl;
      std::cout << "--M-EIT---" << std::endl;
      std::cout << m.eigen_inverse_transformation << std::endl;
      std::cout << "------" << std::endl;
      std::cout << "--M-ET-*-THIS-EIT-" << std::endl;
      std::cout << m.eigen_transformation * eigen_inverse_transformation << std::endl;
      std::cout << "------" << std::endl;
      std::cout << "--THIS-ET-*-M-EIT" << std::endl;
      std::cout << eigen_transformation * m.eigen_inverse_transformation << std::endl;
      std::cout << "------" << std::endl;
      std::cout << "vals dist1,dist2 : " << eigen_dis1 << " " << eigen_dis2 << std::endl;
    }
#endif

    return max(eigen_dis1, eigen_dis2);
  }

  void construct(const std::vector<Vector_d>& evecs,
                 const std::vector<FT>& vals,
                 const FT& epsilon)
  {
    std::vector<FT> evals(vals.size());
    for(std::size_t i=0; i<vals.size(); ++i)
      evals[i] = (std::max)(epsilon, std::abs(vals[i]));

    typename Kd::Compute_coordinate_d coord = Kd().compute_coordinate_d_object();

    E_Matrix eigen_m;

    for(int i=0; i<Dim::value; ++i)
      for(int j=0; j<Dim::value; ++j)
        eigen_m(i,j) = coord(evecs[j], i);

    E_Matrix eigen_diag = E_Matrix::Zero();
    for(int i=0; i<Dim::value; ++i)
      eigen_diag(i, i) = (std::max)(epsilon, std::abs(evals[i]));

    E_Matrix eigen_mtransp = eigen_m.transpose();
    eigen_transformation = eigen_m * eigen_diag * eigen_mtransp;

#ifdef ANISO_DEBUG_METRIC
    std::cout << "--- eigen_m ---" << std::endl;
    std::cout << eigen_m << std::endl;
    std::cout << "--- eigen_diag ---" << std::endl;
    std::cout << eigen_diag << std::endl;
    std::cout << "-- eigen_mtransp ---" << std::endl;
    std::cout << eigen_mtransp << std::endl;
    std::cout << "-- eigen_transformation ---" << std::endl;
    std::cout << eigen_transformation << std::endl;
#endif

    for(int i=0; i<Dim::value; ++i)
      eigen_diag(i, i) = 1./evals[i];

    eigen_inverse_transformation = eigen_m * eigen_diag * eigen_mtransp;

#ifdef ANISO_DEBUG_METRIC
    std::cout << "-- eigen_inverse_transformation ---" << std::endl;
    std::cout << eigen_inverse_transformation << std::endl;
#endif
  }

public:
  friend
  std::ostream& operator<<(std::ostream& out, const Metric_base& x)
  {
    out << "M  = " << x.eigen_transformation << std::endl;
    return out;
  }

  Metric_base()
  {
    eigen_transformation = E_Matrix::Identity();
    eigen_inverse_transformation = E_Matrix::Identity();
    evals = std::vector<FT>(Dim::value, 1.);

    typename Kd::Construct_vector_d constr_vec = Kd().construct_vector_d_object();

    evecs.clear();
    for(int i=0; i<Dim::value; ++i)
    {
      E_Vector vi = E_Vector::Zero();
      vi(i) = 1.;
      evecs.push_back(constr_vec(Dim::value, vi.data(), vi.data() + Dim::value));
    }
  }

  Metric_base(E_Matrix eigen_transformation_)
  {
    //std::cout << "building a metric from a matrix3d with input : " << std::endl;
    //std::cout << eigen_transformation_ << std::endl;

    eigen_transformation = eigen_transformation_;

    bool invertible;
    FT determinant;
    eigen_transformation.computeInverseAndDetWithCheck(eigen_inverse_transformation, determinant, invertible);
    /*
    if(invertible)
      std::cout << "It is invertible, and its inverse is:" << std::endl << eigen_inverse_transformation << std::endl;
    else
      std::cout << "Not invertible and determinant was : " << determinant << std::endl;
*/
  }

  Metric_base(const std::vector<Vector_d>& evecs,
              const std::vector<FT>& evals,
              const FT& epsilon)
  {
#ifdef ANISO_DEBUG_METRIC
    typename Kd::Scalar_product_d inner_pdct = Kd().scalar_product_d_object();

    std::cout << "produits scalaires : " << std::endl;
    for(int i=0; i<Dim::value; ++i)
    {
      for(int j=i+1; j<Dim::value; ++j)
      {
        FT sp = inner_pdct(evecs[i], evecs[j]);
        std::cout << "ij: " << i << " " << j << " " << sp << std::endl;
        if(sp > 1e-4)
          std::cout << "Warning : This dot product should be 0" << std::endl;
      }
    }

    std::cout << "sq_norms           : " << std::endl;
    for(int i=0; i<Dim::value; ++i)
    {
      FT sql = evecs[i].squared_length();
      std::cout << "i: " << i << " sq_n: " << sql << std::endl;
      if(std::abs(sql-1)>1e-4)
        std::cout << "Warning : This sqnorm should be 1" << std::endl;
    }

    std::cout << "valeurs propres    : " << std::endl;
    for(int i=0; i<Dim::value; ++i)
      std::cout << "i: " << i << " eval: " << evals[i] << std::endl;
#endif
    construct(evecs, evals, epsilon);
  }

  ~Metric_base() { }
};

} // Anisotropic_mesh_TC
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_TC_METRIC_H

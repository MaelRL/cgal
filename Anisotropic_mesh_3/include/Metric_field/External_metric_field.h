#ifndef CGAL_ANISOTROPIC_MESH_3_EXTERNAL_METRIC_FIELD_H
#define CGAL_ANISOTROPIC_MESH_3_EXTERNAL_METRIC_FIELD_H

#include <CGAL/Constrain_surface_3_polyhedral.h>
#include <CGAL/Metric_field.h>

#include <CGAL/helpers/metric_helper.h>

#include <fstream>
#include <string>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K>
class External_metric_field : public Metric_field<K>
{
public:
  typedef Metric_field<K>                                       Base;
  typedef typename Base::FT                                     FT;
  typedef typename Base::Metric                                 Metric;
  typedef typename Base::Point_3                                Point_3;
  typedef typename Base::Vector_3                               Vector_3;

  typedef Constrain_surface_3_polyhedral<K>                     Constrain_surface;
  typedef typename Constrain_surface::Vertex_handle             Vertex_handle;
  typedef typename Constrain_surface::Facet_handle              Facet_handle;

  typedef typename Constrain_surface::Primitive                 Primitive;
  typedef typename Constrain_surface::Traits                    Traits;
  typedef typename Constrain_surface::Tree                      Tree;
  typedef typename Tree::Object_and_primitive_id                Object_and_primitive_id;
  typedef typename Tree::Point_and_primitive_id                 Point_and_primitive_id;

private:
  const Constrain_surface& m_pConstrain;
  std::vector<Eigen::Matrix3d> m_metrics; // !! THIS IS MP NOT FP !!

private:
  void fetch_metric(const char* mf_filename)
  {
    std::cout << "Reading metric file: " << mf_filename << std::endl;
    std::ifstream mf_in(mf_filename);
    int size, format;
    //size: number of vertices
    //format: 6 for upper half matrix a00 a01 a02 a11 a12 a22
    //        9 for full matrix a00 a01 a02 a10 etc.
    //        3 for diagonal a00 a11 a22
    //        1 for function at the vertex TODO

    mf_in >> size >> format;
    std::cout << "MF File for a mesh of size: " << size << ". ";
    std::cout << "Format: " << format << std::endl;

    for(int j=0; j<size; ++j)
    {
      Eigen::Matrix3d met = Eigen::Matrix3d::Zero();
      if(format == 3 && mf_in >> met(0,0) >> met(1,1) >> met(2,2))
      { }
      else if(format == 6 &&
              mf_in >> met(0,0) >> met(0,1) >> met(0,2)
                    >> met(1,1) >> met(1,2) >> met(2,2))
      {
        met(1,0) = met(0,1); met(2,0) = met(0,2); met(2,1) = met(1,2);
      }
      else if(format == 9 &&
              mf_in >> met(0,0) >> met(0,1) >> met(0,2) >> met(1,0) >> met(1,1)
                    >> met(1,2) >> met(2,0) >> met(2,1) >> met(2,2))
      { }
      else if(format == 1)
      {
        std::cout << "todo!" << std::endl;
        //compute Hessian and stuff. Call another function
        //std::vector<double> values;
        //compute_metric(values, m_metrics, etc.);
      }
      else
        break;
      m_metrics.push_back(met);
    }

    if(m_metrics.size() != size)
    {
      std::cout << "problem while reading metric field input" << std::endl;
      std::cout << "File says: " << size << " but only read " << m_metrics.size() << std::endl;
    }

    std::cout << "End of read metric field" << std::endl;
  }

  void tensor_frame_polyhedral_surface(const Point_3 &p,
                                       Vector_3& v0, Vector_3& v1, Vector_3& v2,
                                       double& e0, double& e1, double& e2) const
  {
    //Compute with interpolations & barycentric coordinates
    Facet_handle facet = m_pConstrain.find_nearest_facet(p);
    Vertex_handle va = facet->halfedge()->vertex();
    Vertex_handle vb = facet->halfedge()->next()->vertex();
    Vertex_handle vc = facet->halfedge()->next()->next()->vertex();

    //get bary weights
    Vector_3 v00(va->point(), vb->point());
    Vector_3 v11(va->point(), vc->point());
    Vector_3 v22(va->point(), p);
    FT d00 = v00*v00;
    FT d01 = v00*v11;
    FT d11 = v11*v11;
    FT d20 = v22*v00;
    FT d21 = v22*v11;
    FT denom = d00 * d11 - d01 * d01;
    FT v = (d11 * d20 - d01 * d21) / denom;
    FT w = (d00 * d21 - d01 * d20) / denom;
    FT u = 1.0f - v - w;

    std::vector<std::pair<Eigen::Matrix3d, FT> > w_metrics;
    w_metrics.push_back(std::make_pair(m_metrics[va->tag()], u));
    w_metrics.push_back(std::make_pair(m_metrics[vb->tag()], v));
    w_metrics.push_back(std::make_pair(m_metrics[vc->tag()], w));
    Eigen::Matrix3d m = Eigen::Matrix3d::Zero();
    m = logexp_interpolate<K>(w_metrics);

    get_eigen_vecs_and_vals<K>(m, v0, v1, v2, e0, e1, e2);
  }

  void tensor_frame_3D(const Point_3 &p,
                       Vector_3& v0, Vector_3& v1, Vector_3& v2,
                       double& e0, double& e1, double& e2) const
  {
/* TODO
    Point_and_primitive_id pp = m_pConstrain.m_tet_tree.closest_point_and_primitive(p);
    Cell_ijkl_with_coords_access<K>* closest_tet = pp.second; // closest primitive id

    int pid0 = closest_tet->m_cell.vertex(0);
    int pid1 = closest_tet->m_cell.vertex(1);
    int pid2 = closest_tet->m_cell.vertex(2);
    int pid3 = closest_tet->m_cell.vertex(3);

    Point_3 p0 = m_coords[pid0];
    Point_3 p1 = m_coords[pid1];
    Point_3 p2 = m_coords[pid2];
    Point_3 p3 = m_coords[pid3];

//    std::cout << "locate point p : " << p << std::endl;
//    std::cout << pid0 << " " << p0 << std::endl;
//    std::cout << pid1 << " " << p1 << std::endl;
//    std::cout << pid2 << " " << p2 << std::endl;
//    std::cout << pid3 << " " << p3 << std::endl;
//    Point_3 closest_point = pp.first;
//    std::cout << "closest is (should be p) : " << closest_point << std::endl;

    //get bary weights
    Vector_3 v30(p3, p0);
    Vector_3 v31(p3, p1);
    Vector_3 v32(p3, p2);
    Vector_3 v3b(p3, p);

    FT lambda_0 = CGAL::determinant(v3b.x(), v31.x(), v32.x(),
                                    v3b.y(), v31.y(), v32.y(),
                                    v3b.z(), v31.z(), v32.z());
    FT lambda_1 = CGAL::determinant(v30.x(), v3b.x(), v32.x(),
                                    v30.y(), v3b.y(), v32.y(),
                                    v30.z(), v3b.z(), v32.z());
    FT lambda_2 = CGAL::determinant(v30.x(), v31.x(), v3b.x(),
                                    v30.y(), v31.y(), v3b.y(),
                                    v30.z(), v31.z(), v3b.z());

    FT denom = CGAL::determinant(v30.x(), v31.x(), v32.x(),
                                 v30.y(), v31.y(), v32.y(),
                                 v30.z(), v31.z(), v32.z());

    lambda_0 /= denom;
    lambda_1 /= denom;
    lambda_2 /= denom;
    FT lambda_3 = 1 - lambda_0 - lambda_1 - lambda_2;

    //interpolate the metrics
    std::vector<std::pair<Eigen::Matrix3d, FT> > w_metrics;
    w_metrics.push_back(std::make_pair(m_metrics[pid0], lambda_0));
    w_metrics.push_back(std::make_pair(m_metrics[pid1], lambda_1));
    w_metrics.push_back(std::make_pair(m_metrics[pid2], lambda_2));
    w_metrics.push_back(std::make_pair(m_metrics[pid3], lambda_3));
    Eigen::Matrix3d m = logexp_interpolate<K>(w_metrics);

    get_eigen_vecs_and_vals<K>(m, v0, v1, v2, e0, e1, e2);
*/
  }

public:
  Metric compute_metric(const Point_3 &p) const
  {
    Vector_3 v0, v1, v2;
    double e0, e1, e2;

    if(1)
      tensor_frame_polyhedral_surface(p, v0, v1, v2, e0, e1, e2);
    else
      tensor_frame_3D(p, v0, v1, v2, e0, e1, e2);

    return this->build_metric(v0, v1, v2, std::sqrt(e0), std::sqrt(e1), std::sqrt(e2));

  }

  void report(typename std::ofstream &fx) const { }

  External_metric_field(const Constrain_surface& surface_,
                        const char* mf_filename = "metric_field.txt",
                        FT epsilon_ = 1e-6,
                        const double& en_factor = 0.999)
  :
    Metric_field<K>(epsilon_, en_factor),
    m_pConstrain(surface_),
    m_metrics()
  {
    fetch_metric(mf_filename);
  }
};

} // Anisotropic_mesh_3
} // CGAL

#endif

#ifndef CGAL_ANISOTROPIC_MESH_3_POLYHEDRAL_CURVATURE_METRIC_FIELD_H
#define CGAL_ANISOTROPIC_MESH_3_POLYHEDRAL_CURVATURE_METRIC_FIELD_H

#include <CGAL/Constrain_surface_3_polyhedral.h>
#include <CGAL/Polyhedral_metric_field.h>
#include <CGAL/Monge_via_jet_fitting.h>

#include <CGAL/helpers/metric_helper.h>

#include <fstream>
#include <iostream>
#include <utility>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K>
class Polyhedral_curvature_metric_field : public Polyhedral_metric_field<K>
{
public:
  typedef Polyhedral_metric_field<K>                            Base;
  typedef typename K::FT                                        FT;
  typedef typename Base::Metric                                 Metric;
  typedef typename Base::Point_3                                Point_3;
  typedef typename Base::Vector_3                               Vector_3;

  typedef Constrain_surface_3_polyhedral<K>                     Constrain_surface;
  typedef typename Constrain_surface::Vertex_handle             Vertex_handle;

  typedef typename Constrain_surface::C3t3                      C3t3;

  typedef CGAL::Monge_via_jet_fitting<K>                        Monge_via_jet_fitting;
  typedef typename Monge_via_jet_fitting::Monge_form            Monge_form;

public:
  double global_max_sq_eigenvalue() const
  {
    if(this->m_cache_max_sq_eigenvalue)
      return this->m_max_sq_eigenvalue;

    this->m_max_sq_eigenvalue = 0.;
    this->m_min_sq_eigenvalue = DBL_MAX;
    typename C3t3::Triangulation::Finite_vertices_iterator v;
    for(v = this->m_pConstrain.c3t3().triangulation().finite_vertices_begin();
        v != this->m_pConstrain.c3t3().triangulation().finite_vertices_end();
        ++v)
    {
      double minc, maxc;
      min_max_sq_eigenvalues(v->point(), minc, maxc);
      this->m_max_sq_eigenvalue = (std::max)(this->m_max_sq_eigenvalue, maxc);
      this->m_min_sq_eigenvalue = (std::min)(this->m_min_sq_eigenvalue, minc);
    }
    this->m_cache_max_sq_eigenvalue = true;
    this->m_cache_min_sq_eigenvalue = true;
    return this->m_max_sq_eigenvalue;
  }

  double global_min_sq_eigenvalue() const
  {
    if(this->m_cache_min_sq_eigenvalue)
      return this->m_min_sq_eigenvalue;

    global_max_sq_eigenvalue(); //computes both max and min

    return this->m_min_sq_eigenvalue;
  }

  void min_max_sq_eigenvalues(const Point_3& p,
                          double& minc,
                          double& maxc) const
  {
    Vector_3 v0, v1, v2;
    double e0, e1, e2;
    this->tensor_frame_polyhedral_surface(p, v0, v1, v2, e0, e1, e2);
    minc = (std::min)(e0, (std::min)(e1, e2));
    maxc = (std::max)(e0, (std::max)(e1, e2));
  }

private:
  void tensor_frame_on_input_point(Vertex_handle v, Vector_3 &v0/*n*/,
                                   Vector_3 &v1/*vmax*/, Vector_3 &v2/*vmin*/,
                                   FT& c1/*cmax*/, FT& c2/*cmin*/)
  {
    std::vector<Vertex_handle> points;
    this->m_pConstrain.find_nearest_vertices(v, points, 36/*, max_sqd=0 by default*/);

    std::vector<Point_3> points_geo;
    for(std::size_t i=0; i<points.size(); ++i)
      points_geo.push_back(points[i]->point());

    Monge_form monge_form;
    Monge_via_jet_fitting monge_fit;
    monge_form = monge_fit(points_geo.begin(), points_geo.end(), 4, 4);

    FT maxc = monge_form.principal_curvatures(0);
    FT minc = monge_form.principal_curvatures(1);
    v0 = monge_form.normal_direction();

    if(minc >= 0.)
    {
      c1 = (std::max)(this->epsilon, fabs(maxc));
      c2 = (std::max)(this->epsilon, fabs(minc));
      v1 = monge_form.maximal_principal_direction();
      v2 = monge_form.minimal_principal_direction();
    }
    else if(maxc <= 0.)
    {
      c2 = (std::max)(this->epsilon, fabs(maxc));
      c1 = (std::max)(this->epsilon, fabs(minc));
      v2 = monge_form.maximal_principal_direction();
      v1 = monge_form.minimal_principal_direction();
    }
    else //minc < 0 && maxc > 0
    {
      FT abs_min = fabs(minc);
      if(abs_min < maxc)
      {
        c1 = (std::max)(this->epsilon, fabs(maxc));
        c2 = (std::max)(this->epsilon, fabs(minc));
        v1 = monge_form.maximal_principal_direction();
        v2 = monge_form.minimal_principal_direction();
      }
      else
      {
        c2 = (std::max)(this->epsilon, fabs(maxc));
        c1 = (std::max)(this->epsilon, fabs(minc));
        v2 = monge_form.maximal_principal_direction();
        v1 = monge_form.minimal_principal_direction();
      }
    }
  }

  void compute_local_metric()
  {
    std::cout << "\nComputing local metric..." << std::endl;
    this->m_max_sq_eigenvalue = 0.;
    this->m_min_sq_eigenvalue = DBL_MAX;

    std::ofstream out("poly_curvatures.txt");
    out << this->m_pConstrain.m_vertices.size() << std::endl;

    for(std::size_t i=0; i<this->m_pConstrain.m_vertices.size(); ++i)
    {
      if (i % 1000 == 0)
        std::cout << i << " out of " << this->m_pConstrain.m_vertices.size() << std::endl;

      Vector_3 v0, v1, v2;
      FT e0, e1, e2;
      tensor_frame_on_input_point(this->m_pConstrain.m_vertices[i],
                                  v0/*normal*/, v1/*vmax*/, v2/*vmin*/,
                                  e1/*cmax*/, e2/*cmin*/);

      e0 = this->en_factor*(std::max)(e1, e2);
      this->m_max_sq_eigenvalue = (std::max)(this->m_max_sq_eigenvalue, e0);
      this->m_min_sq_eigenvalue = (std::min)(this->m_min_sq_eigenvalue, (std::min)(e1, e2));

      e1 = std::sqrt(e1);
      e2 = std::sqrt(e2);
      e0 = this->en_factor*(std::max)(e1, e2);

      const Point_3& pi = this->m_pConstrain.m_vertices[i]->point();
//    tensor_frame_on_point(pi, vn, v1, v2, e1, e2);

      out << pi.x() << " " << pi.y() << " " << pi.z() << std::endl;
      out << e0 << " " << e1 << " " << e2 << "     ";
      out << v1.x() << " " << v1.y() << " " << v1.z() << "     ";
      out << v2.x() << " " << v2.y() << " " << v2.z() << "     ";
      out << v0.x() << " " << v0.y() << " " << v0.z() << std::endl << std::endl;

      std::size_t index = this->m_pConstrain.m_vertices[i]->tag();
      if(i != index)
        std::cerr << "Error1 in indices in compute_local_metric!" << std::endl;

      this->m_pConstrain.m_vertices[i]->normal() = v0;

      Metric_base<K, K> m(v0, v1, v2, e0, e1, e2, this->epsilon);
      Eigen::Matrix3d transf = m.get_transformation();
      this->m_metrics.push_back(transf.transpose()*transf);

      //if(this->m_metrics.size() != i+1 )
      //  std::cerr << "Error3 in indices in compute_local_metric!" << std::endl;
    }

    this->m_cache_max_sq_eigenvalue = true;
    this->m_cache_min_sq_eigenvalue = true;
    std::cout << "done." << std::endl;
  }

public:
  void report(typename std::ofstream &fx) const { }

  virtual Metric compute_metric(const Point_3 &p) const
  {
    return Base::compute_metric(p);
  }

  Polyhedral_curvature_metric_field(const Constrain_surface& surface_,
                                    const FT epsilon_ = 1e-6,
                                    const FT& en_factor_ = 0.999)
    :
      Base(surface_, epsilon_, en_factor_)
  {
    compute_local_metric();
//    this->global_smooth();
    this->output_mf_polylines();
  }
};

} // Anisotropic_mesh_3
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_POLYHEDRAL_CURVATURE_METRIC_FIELD_H

#ifndef CGAL_ANISOTROPIC_MESH_3_POLYHEDRAL_CURVATURE_METRIC_FIELD_H
#define CGAL_ANISOTROPIC_MESH_3_POLYHEDRAL_CURVATURE_METRIC_FIELD_H

#include <CGAL/Constrain_surface_3_polyhedral.h>
#include <CGAL/Metric_field.h>
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
class Polyhedral_curvature_metric_field : public Metric_field<K>
{
public:
  typedef Metric_field<K>                                       Base;
  typedef typename K::FT                                        FT;
  typedef typename Base::Metric                                 Metric;
  typedef typename Base::Point_3                                Point_3;
  typedef typename Base::Vector_3                               Vector_3;

  typedef Constrain_surface_3_polyhedral<K>                     Constrain_surface;
  typedef typename Constrain_surface::Vertex_handle             Vertex_handle;
  typedef typename Constrain_surface::Facet_handle              Facet_handle;

  typedef typename Constrain_surface::C3t3                      C3t3;

  typedef typename CGAL::Monge_via_jet_fitting<K>               Monge_via_jet_fitting;
  typedef typename Monge_via_jet_fitting::Monge_form            Monge_form;

public:
  const Constrain_surface& m_pConstrain;
  std::vector<typename Eigen::Matrix3d> m_metrics;

  double global_max_sq_eigenvalue() const
  {
    if(this->m_cache_max_sq_eigenvalue)
      return this->m_max_sq_eigenvalue;

    this->m_max_sq_eigenvalue = 0.;
    this->m_min_sq_eigenvalue = DBL_MAX;
    typename C3t3::Triangulation::Finite_vertices_iterator v;
    for(v = m_pConstrain.c3t3().triangulation().finite_vertices_begin();
        v != m_pConstrain.c3t3().triangulation().finite_vertices_end();
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
    tensor_frame(p, v0, v1, v2, e0, e1, e2);
    minc = (std::min)(e0, (std::min)(e1, e2));
    maxc = (std::max)(e0, (std::max)(e1, e2));
  }

  virtual void report(typename std::ofstream &fx) const
  {
    fx << "type: polyhedral surface" << std::endl;
    fx << "epsilon: " << this->epsilon << std::endl;
  }

  void tensor_frame_on_input_point(Vertex_handle v, Vector_3 &v0/*n*/,
                                   Vector_3 &v1/*vmax*/, Vector_3 &v2/*vmin*/,
                                   FT& c1/*cmax*/, FT& c2/*cmin*/)
  {
    std::vector<Vertex_handle> points;
    m_pConstrain.find_nearest_vertices(v, points, 36/*, max_sqd=0 by default*/);
//          nearest_start_try_radius*nearest_start_try_radius);

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

  void tensor_frame(const Point_3 &p,
                    Vector_3 &v0,     //unit normal
                    Vector_3 &v1,     //unit eigenvector
                    Vector_3 &v2,     //unit eigenvector
                    double& e0,       //eigenvalue corresponding to v0
                    double& e1,       //eigenvalue corresponding to v1
                    double& e2) const //eigenvalue corresponding to v2
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

  void compute_local_metric()
  {
    std::cout << "\nComputing local metric..." << std::endl;
    this->m_max_sq_eigenvalue = 0.;
    this->m_min_sq_eigenvalue = DBL_MAX;

    std::ofstream out("poly_curvatures.txt");
    out << m_pConstrain.m_vertices.size() << std::endl;

    for (std::size_t i = 0; i < m_pConstrain.m_vertices.size(); ++i)
    {
      if (i % 100 == 0)
        std::cout << ".";

      Vector_3 v0, v1, v2;
      FT e0, e1, e2;
      tensor_frame_on_input_point(m_pConstrain.m_vertices[i], v0/*normal*/, v1/*vmax*/, v2/*vmin*/,
                                  e1/*cmax*/, e2/*cmin*/);

      e0 = this->en_factor*(std::max)(e1, e2);
      this->m_max_sq_eigenvalue = (std::max)(this->m_max_sq_eigenvalue, e0);
      this->m_min_sq_eigenvalue = (std::min)(this->m_min_sq_eigenvalue, (std::min)(e1, e2));

      e1 = std::sqrt(e1);
      e2 = std::sqrt(e2);
      e0 = this->en_factor*(std::max)(e1, e2);

      const Point_3& pi = m_pConstrain.m_vertices[i]->point();
//          tensor_frame_on_point(pi, vn, v1, v2, e1, e2);

      out << pi.x() << " " << pi.y() << " " << pi.z() << std::endl;
      out << e0 << " " << e1 << " " << e2 << "     ";
      out << v1.x() << " " << v1.y() << " " << v1.z() << "     ";
      out << v2.x() << " " << v2.y() << " " << v2.z() << "     ";
      out << v0.x() << " " << v0.y() << " " << v0.z() << std::endl;

#ifdef CGAL_DEBUG_OUTPUT_VECTOR_FIELD
      FT f = m_pConstrain->get_bounding_radius()*0.02;
      (*(m_pConstrain->max_vector_field_os)) << "2 " << pi << " " << (pi+f*v1) << std::endl;
      (*(m_pConstrain->min_vector_field_os)) << "2 " << pi << " " << (pi+f*v2) << std::endl;
      (*(m_pConstrain->normals_field_os)) << "2 " << pi << " " << (pi+f*v0) << std::endl;
#endif

      std::size_t index = m_pConstrain.m_vertices[i]->tag();
      if(i != index)
        std::cerr << "Error1 in indices in compute_local_metric!" << std::endl;

      m_pConstrain.m_vertices[i]->normal() = v0;

      Metric_base<K, K> m(v0, v1, v2, e0, e1, e2, this->epsilon);
      Eigen::Matrix3d transf = m.get_transformation();
      m_metrics.push_back(transf.transpose()*transf);

      //if(m_metrics.size() != i+1 )
      //  std::cerr << "Error3 in indices in compute_local_metric!" << std::endl;
    }

    this->m_cache_max_sq_eigenvalue = true;
    this->m_cache_min_sq_eigenvalue = true;
    std::cout << "done." << std::endl;
  }

  void heat_kernel_smoothing()
  {
    std::cout << "Smoothing the metric field" << std::endl;
    int smooth_step_n = 3;

    while(smooth_step_n--)
    {
      std::cout << "computing new metrics... steps left : " << smooth_step_n << std::endl;

      std::vector<typename Eigen::Matrix3d> temp_m_metrics;
      for (std::size_t i = 0; i < m_pConstrain.m_vertices.size(); ++i)
      {
        Vertex_handle vi = m_pConstrain.m_vertices[i];

        std::vector<Vertex_handle> points;
        std::set<Vertex_handle> visited;
        std::map<FT,Vertex_handle> ring;
        std::vector<FT> sq_distances;

        //find first ring neighbours
        visited.insert(vi);
        m_pConstrain.find_ring_vertices(vi->point(), vi, visited, ring);
        typename std::map<FT, Vertex_handle>::iterator it = ring.begin();
        while(it != ring.end())
          points.push_back((*it++).second);

        points.push_back(vi);

        FT dist_sum = 0.;
        for(std::size_t j=0; j<points.size(); ++j)
        {
          sq_distances.push_back(CGAL::squared_distance(vi->point(), points[j]->point()));
          dist_sum += std::sqrt(sq_distances[j]);
        }
        FT sigma = dist_sum/((FT) points.size());
        FT coeff = -1./(2*sigma*sigma);

        //coeff W_sigma
        std::vector<FT> wsigmas(points.size());
        FT w_sigmas_sum = 0.;

        for(std::size_t j=0; j<points.size(); ++j)
        {
          FT sq_dist = sq_distances[j];
          FT w = std::exp(coeff * sq_dist);
          w_sigmas_sum += w;
          wsigmas[j] = w;
        }

        for(std::size_t j=0; j<points.size(); ++j)
          wsigmas[j] /= w_sigmas_sum;

        //compute diffusion for each coeff of the tensor
        Eigen::Matrix3d new_tensor_at_vi = Eigen::Matrix3d::Zero();

        for(int j = 0; j < 3; j++)
          for(int k = 0; k < 3; k++)
            for(std::size_t l=0; l<points.size(); ++l)
              new_tensor_at_vi(j,k) += wsigmas[l] * m_metrics[points[l]->tag()](j,k); //neighbours

        //compute new principal directions & curvature values
        double e0,e1,e2;
        Vector_3 v0,v1,v2;
        get_eigen_vecs_and_vals<K>(new_tensor_at_vi, v0, v1, v2, e0, e1, e2);

        Vector_3 ni = vi->normal();
        FT sp0 = std::abs(v0*ni);
        FT sp1 = std::abs(v1*ni);
        FT sp2 = std::abs(v2*ni);
        FT max_sp = (std::max)(sp0,(std::max)(sp1,sp2));

        FT temp_e, temp_sp;
        Vector_3 temp_v;

        if(sp1 == max_sp)
        {
          temp_e = e0;
          e0 = e1;
          e1 = e2;
          e2 = temp_e;

          temp_v = v0;
          v0 = v1;
          v1 = v2;
          v2 = temp_v;

          temp_sp = sp0;
          sp0 = sp1;
          sp1 = sp2;
          sp2 = temp_sp;
        }
        else if(sp2 == max_sp)
        {
          temp_e = e0;
          e0 = e2;
          e2 = e1;
          e1 = temp_e;

          temp_v = v0;
          v0 = v2;
          v2 = v1;
          v1 = temp_v;

          temp_sp = sp0;
          sp0 = sp2;
          sp2 = sp1;
          sp1 = temp_sp;
        }

        //projection on the tangent plane
        v1 = v1-(v1*ni)*ni;
        v2 = v2-(v2*ni)*ni;

        v1 = v1/CGAL::sqrt(v1*v1);
        v2 = v2/CGAL::sqrt(v2*v2);

        //orthogonalize v2 to v1
        v2 = v2 - (v2*v1)*v1;
        v2 = v2/CGAL::sqrt(v2*v2);

#ifdef ANISO_DEBUG_HEAT_KERNEL
        std::cout << "point is : " << vi->tag() << " || " << vi->point() << std::endl;

        std::cout << points.size() << " neighbours" << std::endl;
        for(std::size_t j=0; j<points.size(); ++j)
            std::cout << " " << points[j]->tag() << " ";
        std::cout << std::endl;

        std::cout << "sigma sum : " << w_sigmas_sum << std::endl;
        for(std::size_t j=0; j<points.size(); ++j)
            std::cout << wsigmas[j] << " ";
        std::cout << std::endl;

        std::cout << "old tensor : " << std::endl << m_metrics[i] << std::endl;
        std::cout << "new tensor : " << std::endl << new_tensor_at_vi << std::endl;

        std::cout << "ni : " << ni << std::endl;
        std::cout << "evs & vs: " << std::endl;
        std::cout << e0 << " || " << v0 << " || " << sp0 << std::endl; //should be the normal in new basis
        std::cout << e1 << " || " << v1 << " || " << sp1 << std::endl;
        std::cout << e2 << " || " << v2 << " || " << sp2 << std::endl;
#endif

        //new metric
        Metric_base<K, K> M(ni, v1, v2, std::sqrt(e0), std::sqrt(e1), std::sqrt(e2), this->epsilon);
        Eigen::Matrix3d transf = M.get_transformation();
        temp_m_metrics.push_back(transf.transpose()*transf);

#ifdef CGAL_DEBUG_OUTPUT_VECTOR_FIELD
        if(smooth_step_n == 0)
        {
          Point_3 pi = vi->point();
          FT f = m_pConstrain->get_bounding_radius()*0.02;
          (*(m_pConstrain->normals_field_os)) << "2 " << pi << " " << (pi+f*v0) << std::endl;
          if(e1 > e2)
          {
            (*(m_pConstrain->max_vector_field_os)) << "2 " << pi << " " << (pi+f*v1) << std::endl;
            (*(m_pConstrain->min_vector_field_os)) << "2 " << pi << " " << (pi+f*v2) << std::endl;
          }
          else
          {
            (*(m_pConstrain->min_vector_field_os)) << "2 " << pi << " " << (pi+f*v1) << std::endl;
            (*(m_pConstrain->max_vector_field_os)) << "2 " << pi << " " << (pi+f*v2) << std::endl;
          }
        }
#endif
      }

      //update all the metrics
      m_metrics.swap(temp_m_metrics);
    }

    std::cout << "smoothed the metric field!" << std::endl;
  }

  virtual Metric compute_metric(const Point_3 &p) const
  {
    Vector_3 v0, v1, v2;
    double e0, e1, e2;
    tensor_frame(p, v0, v1, v2, e0, e1, e2);

    return this->build_metric(v0, v1, v2, std::sqrt(e0), std::sqrt(e1), std::sqrt(e2));
  }

  Polyhedral_curvature_metric_field(const Constrain_surface& surface_,
                                    const FT epsilon_ = 1.0,
                                    const double& en_factor_ = 0.999)
    :
      Metric_field<K>(epsilon_, en_factor_),
      m_pConstrain(surface_)
  {
    m_metrics.reserve(m_pConstrain.m_vertices.size());

    compute_local_metric();
    heat_kernel_smoothing();
  }
};

} // Anisotropic_mesh_3
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_POLYHEDRAL_CURVATURE_METRIC_FIELD_H

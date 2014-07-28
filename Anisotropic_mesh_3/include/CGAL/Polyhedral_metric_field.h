#ifndef CGAL_ANISOTROPIC_MESH_3_POLYHEDRAL_METRIC_FIELD_H
#define CGAL_ANISOTROPIC_MESH_3_POLYHEDRAL_METRIC_FIELD_H

#include <CGAL/Constrain_surface_3_polyhedral.h>
#include <CGAL/Metric_field.h>

#include <CGAL/compute_normal.h>
#include <CGAL/helpers/metric_helper.h>

#include <iostream>
#include <fstream>
#include <utility>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K>
class Polyhedral_metric_field : public Metric_field<K>
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

protected:
  const Constrain_surface& m_pConstrain;
  std::vector<Eigen::Matrix3d> m_metrics; // !! THIS IS MP NOT FP !!

  std::vector<Eigen::Matrix3d> m_old_metrics;

/*
  template<typename Starset>
  void facet_edge_length_histogram_with_original_metric()
  {
    int histogram_size = 1000;
    double step_size = 2.0 / (double) histogram_size;
    std::vector<int> histogram(histogram_size, 0);
    typename Starset::FT val, min_value = 1e30, max_value = -1e30;

    std::size_t N = stars.size();
    for(std::size_t i=0; i<N; ++i)
    {
      typename Starset::Star_handle si = stars[i];
      if(!si->is_surface_star())
        continue;

      std::set<int> done;
      typename Starset::Traits::Compute_squared_distance_3 csd =
                                si->traits()->compute_squared_distance_3_object();

      //compute the old Fp
      double e_0, e_1, e_2;
      Vector_3 v_0, v_1, v_2;
      get_eigen_vecs_and_vals<K>(m_old_metrics[i], v_0, v_1, v_2, e_0, e_1, e_2);

      Eigen::Matrix3d eigen_m;
      eigen_m(0,0) = v_0.x();  eigen_m(0,1) = v_1.x();  eigen_m(0,2) = v_2.x();
      eigen_m(1,0) = v_0.y();  eigen_m(1,1) = v_1.y();  eigen_m(1,2) = v_2.y();
      eigen_m(2,0) = v_0.z();  eigen_m(2,1) = v_1.z();  eigen_m(2,2) = v_2.z();

      Eigen::Matrix3d eigen_diag = Eigen::Matrix3d::Zero();
      eigen_diag(0,0) = e_0;
      eigen_diag(1,1) = e_1;
      eigen_diag(2,2) = e_2;

      Eigen::Matrix3d eigen_mtransp = eigen_m.transpose();
      Eigen::Matrix3d eigen_transformation = eigen_m * eigen_diag * eigen_mtransp;

      typename Starset::Point_3 c = si->center_point();
      typename Starset::TPoint_3 tc = transform(eigen_transformation, c);

      std::cout << "c: " << si->index_in_star_set() << " || " << c << std::endl;
      std::cout << "real tc: " << si->center()->point() << std::endl;
      std::cout << "old tc: " << tc << " || " << c << std::endl;

      typename Starset::Facet_set_iterator fit = si->restricted_facets_begin();
      typename Starset::Facet_set_iterator fend = si->restricted_facets_end();
      for(; fit != fend; ++fit)
      {
        for(int i=1; i<4; ++i)
        {
          typename Starset::Vertex_handle vi = fit->first->vertex((fit->second+i)%4);
          if(vi->info() == si->index_in_star_set())
            continue;

          std::pair<typename std::set<int>::iterator, bool> is_insert_successful;
          is_insert_successful = done.insert(vi->info());
          if(!is_insert_successful.second)
            continue;

          typename Starset::Point_3 pi = si->inverse_transform(vi->point());
          typename Starset::TPoint_3 tpi = transform(eigen_transformation, pi);
          typename Starset::FT d_tc_tpi = csd(tc, tpi);

          val = d_tc_tpi / criteria->squared_facet_circumradius;
          val = std::sqrt(val);

          std::cout << "------- tpi: " << vi->info() << " || " << tpi <<  std::endl;
          std::cout << "real one for comparison: " << vi->point() << std::endl;
          std::cout << "+++++++ " << d_tc_tpi << " " << criteria->squared_facet_circumradius;
          std::cout << " ¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹¹ " << val << std::endl;

          min_value = (std::min)(min_value, val);
          max_value = (std::max)(max_value, val);

          for(int j=0; j<histogram_size; ++j)
            if(j*step_size <= val && val < (j+1)*step_size)
              histogram[j]++;

          if(val >= (histogram_size)*step_size)
            histogram[histogram_size-1]++;
        }
      }
    }

    if(verbose)
    {
      std::cout << "facet edge length histogram with original metric" << std::endl;
      std::cout << "min, max: " << min_value << " " << max_value << std::endl;
    }

  //  output_histogram(histogram, 0., 2., "histogram_facet_edge_length_orig.cvs");
    output_histogram(histogram, 0., 2.*criteria->facet_circumradius, "histogram_facet_edge_length_orig.cvs");
  }
*/

  //todo : coefficient using distances in the metric
  void global_smooth()
  {
    if(m_metrics.empty())
      return;
    m_old_metrics = m_metrics;

    std::cout << "Smoothing the metric field" << std::endl;
    int smooth_step_n = 3;

    double gamma_0 = 3.0;

    while(smooth_step_n--)
    {
      std::cout << "computing new metrics... steps left : " << smooth_step_n << std::endl;
      std::vector<typename Eigen::Matrix3d> temp_m_metrics;

      for (std::size_t i = 0; i < m_pConstrain.m_vertices.size(); ++i)
      {
        Vertex_handle vi = m_pConstrain.m_vertices[i];

        double max_gamma = -1;

        std::vector<Vertex_handle> points;
        std::set<Vertex_handle> visited;
        std::map<FT,Vertex_handle> ring;
        std::vector<FT> sq_distances;

//find first ring neighbours
        visited.insert(vi);
        m_pConstrain.find_ring_vertices(vi->point(), vi, visited, ring);
        typename std::map<FT, Vertex_handle>::iterator it = ring.begin();
        while(it != ring.end())
        {
          points.push_back((*it++).second);

//#define ANISO_LOCAL_SMOOTHING
#ifdef ANISO_LOCAL_SMOOTHING
          max_gamma = (std::max)(max_gamma,
                                 compute_distortion<K>(m_metrics[i], // EXTREMELY UNEFFICIENT
                                                       m_metrics[points.back()->tag()]));
#endif
        }
//add center point (or not)
        points.push_back(vi);
        std::size_t size_of_first_ring = points.size();

//skip or not depending on the maximum distortion in the first ring
#ifdef ANISO_LOCAL_SMOOTHING
        std::cout << "gamma: " << max_gamma << " at " << i << std::endl;
        if(max_gamma < gamma_0)
        {
          temp_m_metrics.push_back(m_metrics[i]);
          continue;
        }
        std::cout << "you're in for a little bit of smoothing boy" << std::endl;
#endif

//compute distances for the first ring
        FT dist_sum = 0.;
        for(std::size_t j=0; j<points.size(); ++j)
        {
          sq_distances.push_back(CGAL::squared_distance(vi->point(), points[j]->point()));
          dist_sum += std::sqrt(sq_distances[j]);
        }

//add a second ring with sq_distance = sq_dist_p_p1 + sq_dist_p1_p2
#if 0
        for(std::size_t j=0; j<size_of_first_ring; ++j)
        {
          ring.clear();
          m_pConstrain.find_ring_vertices(vi->point(), points[j], visited, ring);
          typename std::map<FT, Vertex_handle>::iterator it = ring.begin();
          while(it != ring.end())
          {
            points.push_back((*it++).second);
            sq_distances.push_back(sq_distances[j] +
                                   CGAL::squared_distance(points.back()->point(),
                                                          points[j]->point()));
            dist_sum += std::sqrt(sq_distances.back());
          }
        }

        std::cout << "check points for : " << vi->tag() << std::endl;
        for(std::size_t j=0; j<points.size(); ++j)
          std::cout << points[j]->tag() << " || " << points[j]->point() << std::endl;
#endif

//compute coefficients
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

//compute the new tensor
        Eigen::Matrix3d new_tensor_at_vi = Eigen::Matrix3d::Zero();

#if 1
//compute diffusion for each coeff of the tensor
        for(int j = 0; j < 3; j++)
          for(int k = 0; k < 3; k++)
            for(std::size_t l=0; l<points.size(); ++l)
              new_tensor_at_vi(j,k) += wsigmas[l] * m_metrics[points[l]->tag()](j,k);
#else
//logexp
        std::vector<std::pair<Eigen::Matrix3d, FT> > w_metrics;
        for(std::size_t l=0; l<points.size(); ++l)
          w_metrics.push_back(std::make_pair(m_metrics[points[l]->tag()], wsigmas[l]));

        new_tensor_at_vi = logexp_interpolate<K>(w_metrics);
#endif
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

#ifdef ANISO_DEBUG_SMOOTHING
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
      }

//update all the metrics
      std::cout << temp_m_metrics.size() << " " << m_metrics.size() << std::endl;
      assert(temp_m_metrics.size() == m_metrics.size());
      m_metrics.swap(temp_m_metrics);
    }

    std::cout << "smoothed the metric field!" << std::endl;
  }

  void output_mf_polylines()
  {
    std::cout << "Ouputing polyline metric field etc." << std::endl;

    std::ofstream min_vector_field_os("vector_field_min.polylines.cgal");
    std::ofstream max_vector_field_os("vector_field_max.polylines.cgal");
    std::ofstream normals_field_os("vector_field_normals.polylines.cgal");

    for (std::size_t i = 0; i < m_pConstrain.m_vertices.size(); ++i)
    {
      FT e0, e1, e2;
      Vector_3 v0, v1, v2;

      Point_3 pi = m_pConstrain.m_vertices[i]->point();
      Eigen::Matrix3d m = m_metrics[i];

      get_eigen_vecs_and_vals<K>(m, v0, v1, v2, e0, e1, e2);

      e1 = 1./std::sqrt(std::abs(e1));
      e2 = 1./std::sqrt(std::abs(e2));
      e0 = 1./std::sqrt(std::abs(e0));

      m_pConstrain.m_vertices[i]->normal() = v0;

      FT f = m_pConstrain.get_bounding_radius()*0.02;
      normals_field_os << "2 " << pi << " " << (pi+e0*v0) << std::endl;
      if(e1>e2)
      {
        max_vector_field_os << "2 " << pi << " " << (pi+e1*v1) << std::endl;
        min_vector_field_os << "2 " << pi << " " << (pi+e2*v2) << std::endl;
      }
      else
      {
        min_vector_field_os << "2 " << pi << " " << (pi+e1*v1) << std::endl;
        max_vector_field_os << "2 " << pi << " " << (pi+e2*v2) << std::endl;
      }
    }
  }

  void output_metric_at_input_vertices()
  {
    std::ofstream out("metric_at_input.txt");
    out.precision(20);

    out << m_metrics.size() << " 6" << std::endl; // todo other modes (1, 3 and 9)
    for(std::size_t i=0; i<m_metrics.size(); ++i)
    {
      Eigen::Matrix3d m = m_metrics[i];
      out << m(0,0) << " " << m(0,1) << " " << m(0,2) << " " << m(1,1) << " ";
      out << m(1,2) << " " << m(2,2) << std::endl;
    }
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

  virtual void report(typename std::ofstream &fx) const { }

  Polyhedral_metric_field(const Constrain_surface& surface_,
                        FT epsilon_ = 1e-6,
                        FT en_factor = 0.999)
  :
    Base(epsilon_, en_factor),
    m_pConstrain(surface_),
    m_metrics()
  {
    this->m_metrics.reserve(this->m_pConstrain.m_vertices.size());
  }
};

} // Anisotropic_mesh_3
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_3_POLYHEDRAL_METRIC_FIELD_H

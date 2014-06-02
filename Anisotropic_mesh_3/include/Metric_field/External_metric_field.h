#ifndef CGAL_ANISOTROPIC_MESH_3_EXTERNAL_METRIC_FIELD_H
#define CGAL_ANISOTROPIC_MESH_3_EXTERNAL_METRIC_FIELD_H

#include <CGAL/Constrain_surface_3_polyhedral.h>
#include <CGAL/Polyhedral_metric_field.h>

#include <CGAL/helpers/metric_helper.h>

#include <boost/array.hpp>

#include <fstream>
#include <string>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

template<typename K>
class External_metric_field : public Polyhedral_metric_field<K>
{
public:
  typedef Polyhedral_metric_field<K>                            Base;

  typedef typename Base::Metric                                 Metric;
  typedef typename Base::Constrain_surface                      Constrain_surface;

  typedef typename Base::FT                                     FT;
  typedef typename Base::Point_3                                Point_3;
  typedef typename Base::Vector_3                               Vector_3;


private:
  void fix_normal_eigenvalues()
  {
    for (std::size_t i = 0; i < this->m_pConstrain.m_vertices.size(); ++i)
    {
      boost::array<FT, 3> evs;
      boost::array<Vector_3, 3> vs;

      Eigen::Matrix3d m = this->m_metrics[i];

      get_eigen_vecs_and_vals<K>(m, vs[0], vs[1], vs[2], evs[0], evs[1], evs[2]);

      Vector_3 normal = compute_vertex_normal<typename Constrain_surface::Polyhedron_type, K>(this->m_pConstrain.m_vertices[i]);
      int id_of_vn = -1;
      FT max_sp = -1e30;

      //std::cout << "normal: " << normal << std::endl;
      //std::cout << "evs: " << evs[0] << " " << evs[1] << " " << evs[2] << std::endl;

      for(int j=0; j<3; ++j)
      {
        FT spi = std::abs(vs[j]*normal);
        //std::cout << "compared to: " << vs[j] << " " << spi << std::endl;
        if(spi > max_sp)
        {
          id_of_vn = j;
          max_sp = spi;
        }
      }

      //std::cout << "id of vn: " << id_of_vn << std::endl;

      evs[id_of_vn] = (std::max)(evs[(id_of_vn+1)%3], evs[(id_of_vn+2)%3]); //en_factor applied later

      this->m_pConstrain.m_vertices[i]->normal() = vs[id_of_vn];

      //todo clean that, building m_metrics[i] doesn't need all these steps
      Metric M = this->build_metric(vs[id_of_vn], vs[(id_of_vn+1)%3], vs[(id_of_vn+2)%3],
                                    std::sqrt(evs[id_of_vn]),
                                    std::sqrt(evs[(id_of_vn+1)%3]),
                                    std::sqrt(evs[(id_of_vn+2)%3]));
      Eigen::Matrix3d transf = M.get_transformation();
      this->m_metrics[i] = transf.transpose()*transf;
    }
  }

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
      this->m_metrics.push_back(met);
    }

    if(this->m_metrics.size() != size)
    {
      std::cout << "problem while reading metric field input" << std::endl;
      std::cout << "File says: " << size << " but only read " << this->m_metrics.size() << std::endl;
    }

    std::cout << "End of read metric field" << std::endl;

    fix_normal_eigenvalues();
    this->global_smooth();

    this->output_mf_polylines();
    this->output_metric_at_input_vertices();
  }

public:
  Metric compute_metric(const Point_3 &p) const
  {
    return Base::compute_metric(p);
  }

  void report(typename std::ofstream &fx) const { }

  External_metric_field(const Constrain_surface& surface_,
                        const char* mf_filename = "metric_field.txt",
                        FT epsilon_ = 1e-6,
                        FT en_factor = 0.999)
  :
    Base(surface_, epsilon_, en_factor)
  {
    fetch_metric(mf_filename);
  }
};

} // Anisotropic_mesh_3
} // CGAL

#endif

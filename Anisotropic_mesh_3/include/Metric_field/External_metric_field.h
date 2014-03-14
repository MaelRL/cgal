#ifndef CGAL_ANISOTROPIC_MESH_3_EXTERNAL_METRIC_FIELD_H
#define CGAL_ANISOTROPIC_MESH_3_EXTERNAL_METRIC_FIELD_H

#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/aabb_tree/aabb_tree_tetrahedron_primitive.h>

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
  typedef Metric_field<K>           Base;
  typedef typename Base::FT         FT;
  typedef typename Base::Metric     Metric;
  typedef typename Base::Point_3    Point_3;
  typedef typename Base::Vector_3   Vector_3;

private:
  typedef AABB_tetrahedron_primitive<K>           Tet_primitive;
  typedef CGAL::AABB_traits<K, Tet_primitive>     Traits;
  typedef CGAL::AABB_tree<Traits>                 Tree;

  typedef typename Tree::Object_and_primitive_id  Object_and_primitive_id;
  typedef typename Tree::Point_and_primitive_id   Point_and_primitive_id;
  typedef typename Tree::Primitive_id             Primitive_id;

  Tree m_tet_tree;

  //mesh related
  std::vector<Point_3> m_coords;
  std::vector<Eigen::Matrix3d> m_metrics; // !! THIS IS MP NOT FP !!
  std::vector<Cell_ijkl_with_coords_access<K> > m_tets;

private:
  void fetch_mesh(const char* mesh_filename)
  {
    std::cout << "Reading mesh file: " << mesh_filename << std::endl;
    std::ifstream mesh_in(mesh_filename);
    std::string word;
    int ne, ref, e1, e2, e3, e4;
    FT x,y,z;

    mesh_in >> word >> ne; //assuming ascii format
    mesh_in >> word >> ne; //assuming dimension 3
    while(mesh_in >> word && word != "End")
    {
      mesh_in >> ne;
      std::cout << ne << " Elements" << std::endl;
      for(int i=0; i<ne; ++i)
      {
        if(word == "Vertices")
        {
          mesh_in >> x >> y >> z >> ref;
          m_coords.push_back(Point_3(x,y,z));
        }
        else if(word == "Triangles") //unused
        {
          mesh_in >> e1 >> e2 >> e3 >> ref;
        }
        else if(word == "Tetrahedra")
        {
          mesh_in >> e1 >> e2 >> e3  >> e4 >> ref;
          m_tets.push_back(Cell_ijkl_with_coords_access<K>(e1-1, e2-1, e3-1, e4-1,
                                                           &m_coords));
        }
      }
    }
    std::cout << "Fetched mesh: ";
    std::cout << m_coords.size() << " points. ";
    std::cout << m_tets.size() << " tets" << std::endl;

    build_aabb();
    m_tet_tree.accelerate_distance_queries();
    std::cout << "the tree has size : " << m_tet_tree.size() << std::endl;
  }

  void build_aabb()
  {
    typename std::vector<Cell_ijkl_with_coords_access<K> >::iterator it = m_tets.begin();
    typename std::vector<Cell_ijkl_with_coords_access<K> >::iterator itend = m_tets.end();
    for(;it!=itend;++it)
      m_tet_tree.insert(Tet_primitive(&(*it)));
  }

  void fetch_metric(const char* mf_filename)
  {
    std::cout << "Reading metric file: " << mf_filename << std::endl;
    std::ifstream mf_in(mf_filename);
    int size, format;
    //size: number of vertices
    //format: 6 for lower half matrix a00 a10 a11 a20 a21 a22
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
              mf_in >> met(0,0) >> met(1,0) >> met(1,1)
                    >> met(2,0) >> met(2,1) >> met(2,2))
      {
        met(0,1) = met(1,0); met(0,2) = met(2,0); met(1,2) = met(2,1);
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
    if(m_metrics.size() != m_coords.size())
      std::cout << "N of vertices != N of metric matrices" << std::endl;
    std::cout << "End of read metric field" << std::endl;
  }

public:
  Metric compute_metric(const Point_3 &p) const
  {
    Point_and_primitive_id pp = m_tet_tree.closest_point_and_primitive(p);
    Cell_ijkl_with_coords_access<K>* closest_tet = pp.second; // closest primitive id

    int pid0 = closest_tet->m_cell.vertex(0);
    int pid1 = closest_tet->m_cell.vertex(1);
    int pid2 = closest_tet->m_cell.vertex(2);
    int pid3 = closest_tet->m_cell.vertex(3);

    Point_3 p0 = m_coords[pid0];
    Point_3 p1 = m_coords[pid1];
    Point_3 p2 = m_coords[pid2];
    Point_3 p3 = m_coords[pid3];

    /*
    std::cout << "locate point p : " << p << std::endl;
    std::cout << pid0 << " " << p0 << std::endl;
    std::cout << pid1 << " " << p1 << std::endl;
    std::cout << pid2 << " " << p2 << std::endl;
    std::cout << pid3 << " " << p3 << std::endl;
    Point_3 closest_point = pp.first;
    std::cout << "closest is (should be p) : " << closest_point << std::endl;
    */

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

    //build metric
    Vector_3 vn, v1, v2;
    double en, e1, e2;
    get_eigen_vecs_and_vals<K>(m, vn, v1, v2, en, e1, e2);
    Metric M = this->build_metric(vn, v1, v2, std::sqrt(en), std::sqrt(e1), std::sqrt(e2));

    return M;
  }

  void report(typename std::ofstream &fx) const { }

  External_metric_field(const char* mesh_filename = "mesh.off",
                        const char* mf_filename = "metric_field.txt",
                        FT epsilon_ = 1e-6,
                        const double& en_factor = 0.999)
  :
    Metric_field<K>(epsilon_, en_factor),
    m_tet_tree(),
    m_coords(),
    m_metrics(),
    m_tets()
  {
    fetch_mesh(mesh_filename);
    fetch_metric(mf_filename);
  }
};

} // Anisotropic_mesh_3
} // CGAL

#endif

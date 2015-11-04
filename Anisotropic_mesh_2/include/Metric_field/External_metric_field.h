#ifndef CGAL_ANISOTROPIC_MESH_2_EXTERNAL_METRIC_FIELD_H
#define CGAL_ANISOTROPIC_MESH_2_EXTERNAL_METRIC_FIELD_H

#include <Domain/Rectangle_domain.h>
#include <CGAL/Metric_field.h>

#include <CGAL/helpers/metric_helper.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <Eigen/Dense>
#include <boost/array.hpp>

#include <fstream>
#include <string>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_2
{

template<typename K>
class External_metric_field
    : public Metric_field<K>
{
public:
  typedef Metric_field<K>                                    Base;
  typedef typename Base::Metric                              Metric;

  typedef typename K::FT                                     FT;
  typedef typename K::Point_2                                Point_2;
  typedef typename K::Vector_2                               Vector_2;
  typedef typename K::Triangle_2                             Triangle_2;
  typedef typename K::Point_3                                Point_3;
  typedef typename K::Segment_3                              Segment_3;
  typedef typename K::Triangle_3                             Triangle_3;

  typedef boost::array<std::size_t, 3>                       Triangle;

  struct Base_mesh_primitive
  {
    typedef int                            Id;

    typedef typename K::Point_3            Point;
    typedef typename K::Triangle_3         Datum;

    Id m_it;
    Datum m_datum; // cache the datum
    Point m_p; // cache the reference point

    Base_mesh_primitive() { }
    Base_mesh_primitive(int i, const Datum& datum, const Point& point)
      :
        m_it(i),
        m_datum(datum),
        m_p(point)
    { }

     Id id() const { return m_it; }
     Datum datum() const { return m_datum; }
     Point reference_point() const { return m_p; }
  };

  typedef Base_mesh_primitive                                Primitive;
  typedef typename Primitive::Id                             Primitive_Id;
  typedef CGAL::AABB_traits<K, Primitive>                    AABB_triangle_traits;
  typedef CGAL::AABB_tree<AABB_triangle_traits>              Tree;

  std::vector<Eigen::Matrix2d> m_metrics;

  std::vector<Point_2> m_points;
  std::vector<Triangle> m_triangles;

  Tree tree;

private:
  void build_aabb_tree()
  {
    std::cout << "build aabb tree" << std::endl;
    tree.clear();

    // create an aabb tree of triangles to fasten it
    for(std::size_t i=0; i<m_triangles.size(); ++i)
    {
      const Triangle& t = m_triangles[i];
      const Point_2& p0 = m_points[t[0]];
      const Point_2& p1 = m_points[t[1]];
      const Point_2& p2 = m_points[t[2]];
      Point_3 p(p0.x(), p0.y(), 0.);
      Triangle_3 tr_3(p,
                      Point_3(p1.x(), p1.y(), 0.),
                      Point_3(p2.x(), p2.y(), 0.));
      Primitive pri(i, tr_3, p);
      tree.insert(pri);
    }
  }

  void fetch_geometry(const char* filename)
  {
    std::cout << "Reading geometry file: " << filename << std::endl;
    std::ifstream in(filename);
    if(!in)
    {
      std::cout << "Couldn't open geometry file" << std::endl;
      return;
    }

    // read and create the base mesh points
    std::string word;
    std::size_t useless, nv, nt, dim;
    FT r_x, r_y;

    in >> word >> useless; // MeshVersionFormatted i
    in >> word >> dim; // Dimension d
    CGAL_assertion(dim == 2);

    while(in >> word && word != "End")
    {
      if(word == "Vertices")
      {
        in >> nv;
        for(std::size_t i=0; i<nv; ++i)
        {
          in >> r_x >> r_y >> useless;
          Point_2 p(r_x, r_y);
          m_points.push_back(p);
        }
      }
      else if(word == "Edges")
      {
        std::size_t ne;
        in >> ne;
        for(std::size_t i=0; i<ne; ++i)
          in >> useless >> useless >> useless;
      }
      else if(word == "Triangles")
      {
        // read the triangles and assign the neighbors
        in >> nt;
        for(std::size_t i=0; i<nt; ++i)
        {
          Triangle tr;
          in >> tr[0] >> tr[1] >> tr[2] >> useless;
          --tr[0]; --tr[1]; --tr[2]; // because we're reading medit data

          m_triangles.push_back(tr);
        }
      }
      else
      {
        std::cout << "weird words in your input...: " << word << std::endl;
        return;
      }
    }

    std::cout << "base mesh initialized" << std::endl;
    std::cout << nv << " vertices" << std::endl;
    std::cout << nt << " triangles" << std::endl;

    build_aabb_tree();
  }

  void fetch_metric(const char* filename)
  {
    std::cout << "Reading metric file: " << filename << std::endl;
    std::ifstream in(filename);

    if(!in)
    {
      std::cout << "Couldn't open metric file" << std::endl;
      return;
    }

    std::size_t size, format;
    //size: number of vertices
    //format: 4 for full matrix a00 a01 a10 a11
    //        3 for upper half matrix a00 a01 a11
    //        2 for diagonal a00 a11
    //        1 for function at the vertex and we compute the Hessian yadayada TODO

    in >> size >> format;
    std::cout << "MF File for a mesh of size: " << size << ". ";
    std::cout << "Format: " << format << std::endl;

    for(std::size_t j=0; j<size; ++j)
    {
      Eigen::Matrix2d met = Eigen::Matrix2d::Zero();
      if(format == 4 && in >> met(0,0) >> met(0,1) >> met(1,0) >> met(1,1))
      { }
      else if(format == 3 && in >> met(0,0) >> met(0,1) >> met(1,1))
      {
        met(1,0) = met(0,1);
      }
      else if(format == 3 && in >> met(0,0) >> met(1,1))
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
      std::cout << "File says: " << size << " but only read " << this->m_metrics.size() << std::endl;
    }

    std::cout << "End of read metric field" << std::endl;

//    this->smooth();
    this->draw(m_points);
  }

public:
  void barycentric_coordinates(const Point_2& pa, const Point_2& pb, const Point_2& pc,
                               const Point_2& p, FT& a, FT& b, FT& c) const
  {
    Vector_2 v0(pa, pb);
    Vector_2 v1(pa, pc);
    Vector_2 v2(pa, p);
    FT d00 = v0*v0;
    FT d01 = v0*v1;
    FT d11 = v1*v1;
    FT d20 = v2*v0;
    FT d21 = v2*v1;
    FT denom = 1./ (d00 * d11 - d01 * d01);
    b = (d11 * d20 - d01 * d21) * denom;
    c = (d00 * d21 - d01 * d20) * denom;
    a = 1.0f - b - c;
  }

  Metric compute_metric(const Point_2 &p) const
  {
    Segment_3 query(Point_3(p.x(), p.y(), -1.), Point_3(p.x(), p.y(), 1.));
    std::list<Primitive_Id> intersections;
    tree.all_intersected_primitives(query, std::back_inserter(intersections));

    typename std::list<Primitive_Id>::const_iterator it = intersections.begin(),
                                                     end = intersections.end();
    for(; it!=end; ++it)
    {
      const Triangle& tr = m_triangles[*it];
      const Triangle_2 tr_2(m_points[tr[0]],
                            m_points[tr[1]],
                            m_points[tr[2]]);

      if(K().has_on_bounded_side_2_object()(tr_2, p) ||
         K().has_on_boundary_2_object()(tr_2, p))
      {
#if (verbose > 10)
        std::cout << "locating point " << seed_id
                  << " point: " << p.x() << " " << p.y() << std::endl;
        std::cout << "found triangle: " << *it << std::endl;
        std::cout << tr[0] << " [" << m_points[tr[0]].point << "] " << std::endl;
        std::cout << tr[1] << " [" << m_points[tr[1]].point << "] " << std::endl;
        std::cout << tr[2] << " [" << m_points[tr[2]].point << "] " << std::endl;
#endif

        // we're inside! compute the barycentric coordinates
        FT a, b, c;
        barycentric_coordinates(m_points[tr[0]], m_points[tr[1]], m_points[tr[2]],
                                p, a, b, c);

        std::vector<std::pair<Eigen::Matrix2d, FT> > w_metrics;
        w_metrics.push_back(std::make_pair(m_metrics[tr[0]], a));
        w_metrics.push_back(std::make_pair(m_metrics[tr[1]], b));
        w_metrics.push_back(std::make_pair(m_metrics[tr[2]], c));
        Eigen::Matrix2d m = Eigen::Matrix2d::Zero();
        m = logexp_interpolate<K>(w_metrics);

        Vector_2 v0, v1;
        FT e0, e1;
        get_eigen_vecs_and_vals<K>(m, v0, v1, e0, e1);
        return this->build_metric(v0, v1, std::sqrt(e0), std::sqrt(e1));
      }
    }
    CGAL_postcondition(false);
    return Metric();
  }

  void report(typename std::ofstream& ) const { }

  External_metric_field(const char* geometry_filename,
                        const char* mf_filename,
                        FT epsilon_ = 1e-6)
  :
    Base(epsilon_)
  {
    fetch_geometry(geometry_filename);
    fetch_metric(mf_filename);
  }
};

} // Anisotropic_mesh_2
} // CGAL

#endif // CGAL_ANISOTROPIC_MESH_2_EXTERNAL_METRIC_FIELD_H

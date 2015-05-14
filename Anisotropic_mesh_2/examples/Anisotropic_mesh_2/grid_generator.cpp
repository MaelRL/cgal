#include <CGAL/Starset.h>
#include <CGAL/Stretched_Delaunay_2.h>

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <Eigen/Dense>

#include <fstream>
#include <iostream>
#include <vector>

using namespace CGAL::Anisotropic_mesh_2;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef typename K::FT                                       FT;
typedef Stretched_Delaunay_2<K>                              Star;
typedef typename Star::Point_2                               Point_2;
typedef typename Star::TPoint_2                              TPoint_2;
typedef typename Star::Metric                                Metric;
typedef typename Star::Traits                                Traits;

typedef typename Eigen::Matrix<double, 5, 1>                 Vector5d;
typedef typename Eigen::Matrix<double, 2, 1>                 Vector2d;

const int vertices_nv = 1010;
//#define R2
#ifdef R2
const FT grid_side = 4;
FT offset_x = -grid_side/2.; // offset is the bottom left point
FT offset_y = -grid_side/2.; // todo normalize this with aniso_mesh_2's rectangle whose offset is the center of the rectangle...
#else
const FT grid_side = 4;
FT offset_x = -grid_side/2.;//-2.75; // offset is the bottom left point
FT offset_y = -grid_side/2.;//0.9;
#endif

Vector5d compute_hat(const Point_2& p,
                     const Metric& m)
{
  Eigen::Vector2d p_base;
  Vector5d p_hat;
  Eigen::Vector2d p_circle;
  Eigen::Vector3d p_wave;
  Eigen::Matrix2d mat = m.get_mat();
#if 1
  p_base[0] = p.x(); p_base[1] = p.y();
  p_circle = mat*p_base;
  p_wave[0] = -0.5*mat(0,0); p_wave[1] = -mat(0,1); p_wave[2] = -0.5*mat(1,1);
  p_hat[0] = p_circle[0]; p_hat[1] = p_circle[1];
  p_hat[2] = p_wave[0]; p_hat[3] = p_wave[1]; p_hat[4] = p_wave[2];
#elif 1
  p_base[0] = p.x(); p_base[1] = p.y();
  p_circle = mat*p_base;
  p_wave[0] = mat(0,0); p_wave[1] = mat(0,1); p_wave[2] = mat(1,1);
  p_hat[0] = p_circle[0]; p_hat[1] = p_circle[1];
  p_hat[2] = p_wave[0]; p_hat[3] = p_wave[1]; p_hat[4] = p_wave[2];
#else
  p_hat[0] = p.x(); p_hat[1] = p.y();
  p_hat[2] = mat(0,0); p_hat[3] = mat(0,1); p_hat[4] = mat(1,1);
#endif

  return p_hat;
}

template<typename Metric_field>
void build_seeds(std::vector<Point_2>& seeds,
                 std::vector<Vector5d>& R5seeds,
                 std::vector<Metric>& seeds_m,
                 std::vector<FT>& ws,
                 const Metric_field& mf)
{
  srand (0);

#if 1
  std::ifstream in("bambimboum.mesh");
//  std::ifstream in("/home/mrouxel/cgal/Anisotropic_mesh_TC/examples/Anisotropic_mesh_TC/build/aniso_TC.mesh");
  std::string word;
  std::size_t useless, nv, dim;
  FT r_x, r_y;

  in >> word >> useless; //MeshVersionFormatted i
  in >> word >> dim; //Dimension d
  in >> word >> nv;
  std::cout << "nv: " << nv << std::endl;
  assert(dim == 2);
  if(vertices_nv < nv)
    nv = vertices_nv;

  seeds.resize(nv); R5seeds.resize(nv); seeds_m.resize(nv); ws.resize(nv);

  for(std::size_t i=0; i<nv; ++i)
  {
    in >> r_x >> r_y >> useless;
#else
  for(std::size_t i=0; i<seeds.size(); ++i)
  {
    double r_x = offset_x + grid_side * ((double) rand() / (RAND_MAX));
    double r_y = offset_y + grid_side * ((double) rand() / (RAND_MAX));
#endif
    seeds[i] = Point_2(r_x, r_y);
    std::cout << "Seeds i: " << i << " " << seeds[i].x() << " " << seeds[i].y() << std::endl;

    seeds_m[i] = mf->compute_metric(seeds[i]);
    R5seeds[i] = compute_hat(seeds[i], seeds_m[i]);
    ws[i] = R5seeds[i].norm()*R5seeds[i].norm() - seeds[i].x()*R5seeds[i](0)
                                                  - seeds[i].y()*R5seeds[i](1);

//    std::cout << "check " << std::endl;
//    std::cout << seeds[i].x() << " " << seeds[i].y() << std::endl;
//    std::cout << seeds_m[i].get_mat() << std::endl;
//    std::cout << R5seeds[i](0) << " " << R5seeds[i](1) << std::endl;
//    std::cout << ws[i] << std::endl;
  }
}

void read_tangent_plane_seeds(std::vector<Vector2d>& R2seeds,
                              std::vector<FT>& r2ws)
{
//  std::ifstream in("/home/mrouxel/cgal/Anisotropic_mesh_TC/examples/Anisotropic_mesh_TC/build/projected_points.txt");
  std::ifstream in("/home/mrouxel/cgal/Tangential_complex/test/Tangential_complex/build/projected_points_120.txt");
  std::size_t nv, dim;
  FT x, y, w;

  in >> dim >> nv;
  assert(dim == 2); //tmp
  if(vertices_nv < nv)
    nv = vertices_nv;
  R2seeds.resize(nv); r2ws.resize(nv);

  for(std::size_t i=0; i<nv; ++i)
  {
    in >> x >> y >> w;
    R2seeds[i] = (Vector2d() << x, y ).finished();
    r2ws[i] = w;
  }
}

int main(int, char**)
{
  std::srand(0);

  double n = 1000.; // number of points per side
  double a = grid_side; // length of the side
  double step = a/n;
  double sq_n = n*n;
  double tri_n = 2*(n-1)*(n-1);
  std::size_t counter = 0;
  Traits* traits = new Traits();
  typename Star::Traits::Compute_squared_distance_2 csd =
      traits->compute_squared_distance_2_object();

  //Euclidean_metric_field<K>* metric_field = new Euclidean_metric_field<K>(16,1);
  Custom_metric_field<K>* metric_field = new Custom_metric_field<K>();

  std::ofstream out_bb("grid.bb");
  out_bb << "2 1 " << (int) sq_n << " 2" << std::endl;

  std::ofstream out("grid.mesh");
  out << "MeshVersionFormatted 1" << std::endl;
  out << "Dimension 2" << std::endl;
  out << "Vertices" << std::endl;
  out << (int) sq_n << std::endl;

  Point_2 center(0.5, 0.5);
  Metric m_c = metric_field->compute_metric(center);
  Vector5d m_hat(compute_hat(center, m_c));

  std::vector<Point_2> seeds; // p
  std::vector<Metric> seeds_m; // p metrics
  std::vector<Vector5d> R5seeds; // p_hat
  std::vector<FT> ws; // p_hat weights

  std::vector<Vector2d> R2seeds; // Pi_H(p_hat) on some tangent plane H
  std::vector<FT> r2ws; // corresponding weights
#ifdef R2
  read_tangent_plane_seeds(R2seeds, r2ws);
#else
  build_seeds(seeds, R5seeds, seeds_m, ws, metric_field);
#endif

  for(unsigned int i=0; i<n; ++i)
  {
    for(unsigned int j=0; j<n ;++j)
    {
      //filling from bot left to top right
      Point_2 p(offset_x+j*step, offset_y+i*step);
      Vector2d pv = (Vector2d() << p.x(), p.y()).finished();
      out << p.x() << " " << p.y() << " " << ++counter << std::endl;

      Metric m_p = metric_field->compute_metric(p);

// FIRST MODE : COMPUTE THE DISTANCE IN THE METRIC BETWEEN THE POINTS IN R^2
      TPoint_2 tp(m_p.transform(p));
      TPoint_2 tc(m_p.transform(center));
//      out_bb << csd(tp, tc) << std::endl;

      Vector5d p_hat(compute_hat(p, m_p));
      // SECOND MODE : COMPUTE THE EUCLIDEAN DISTANCE BETWEEN THE POINTS IN R^5
//      out_bb << (p_hat - m_hat).norm() << std::endl;

// THIRD MODE : BUILD A VORONOI
      FT min = 1e30;
      unsigned int min_id = 0;

#ifdef R2
      for(std::size_t k=0; k<R2seeds.size(); ++k)
#else
      for(std::size_t k=4; k<seeds.size(); ++k) // !!! IGNORING THE FIRST FOUR SEEDS !!!
#endif
      {
#ifndef R2
        // EUCLIDEAN DISTANCE IN R^5
//        double sq_d = (p_hat - R5seeds[k]).norm();

        // WEIGHTED DISTANCE IN R^5
//        double sq_d = std::pow((p_hat - R5seeds[k]).norm(),2.) - ws[k];

        // DU WANG : dx(p,x) < dx(q,x)
//        double sq_d = csd(tp, m_p.transform(seeds[k]));

        // LABELLE SHEWCHUK : dp(p,x) < dq(q,x)
        TPoint_2 tp2 = (seeds_m[k]).transform(p);
        TPoint_2 ts = (seeds_m[k]).transform(seeds[k]);
        double sq_d = csd(tp2, ts);
#else
        // WEIGHTED DISTANCE IN THE TANGENT PLANE (seeds are the projected p_hat)
        double sq_d = std::pow((pv - R2seeds[k]).norm(), 2.) - r2ws[k];
#endif
        if(sq_d < min)
        {
          if(std::abs(sq_d-min)<0.05*min)
            min_id = 1.01*vertices_nv; // drawing borders
          else
            min_id = k;
          min = sq_d;
        }
      }


      out_bb << min_id << std::endl;

//      std::cout << "min id: " << min_id << " " << min << std::endl;
    }
  }

  out << "Triangles" << std::endl;
  out << (int) tri_n << std::endl;

  for(int i=1; i<n; ++i)
  {
    for(int j=0; j<(n-1); ++j)
    {
      int id1 = j*n + i;
      int id2 = j*n + (i+1);
      int id3 = (j+1)*n + i;
      int id4 = (j+1)*n + (i+1);
      out << id1 << " " << id2 << " " << id4 << " 1" << std::endl;
      out << id1 << " " << id3 << " " << id4 << " 1" << std::endl;
    }
  }

  out << "End" << std::endl;
}

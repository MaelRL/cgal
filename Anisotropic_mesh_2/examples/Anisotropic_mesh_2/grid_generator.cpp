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

const FT offset_x = -0.55; // offset is the bottom left point
const FT offset_y = -0.55; // todo normalize this with aniso_mesh_2's rectangle
const FT grid_side = 1.1; // whose offset is the center of the rectangle...

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
                 std::vector<FT>& sqws,
                 const Metric_field& mf)
{
  srand (0);

  for(unsigned int i=0; i<seeds.size(); ++i)
  {
    double r_x = offset_x + grid_side * ((double) rand() / (RAND_MAX));
    double r_y = offset_y + grid_side * ((double) rand() / (RAND_MAX));
    seeds[i] = Point_2(r_x, r_y);
    std::cout << "Seeds i: " << i << " " << seeds[i].x() << " " << seeds[i].y() << std::endl;

    seeds_m[i] = mf->compute_metric(seeds[i]);
    R5seeds[i] = compute_hat(seeds[i], seeds_m[i]);
    sqws[i] = R5seeds[i].norm()*R5seeds[i].norm() - seeds[i].x()*R5seeds[i](0)
                                                  - seeds[i].y()*R5seeds[i](1);

//    std::cout << "check " << std::endl;
//    std::cout << seeds[i].x() << " " << seeds[i].y() << std::endl;
//    std::cout << seeds_m[i].get_mat() << std::endl;
//    std::cout << R5seeds[i](0) << " " << R5seeds[i](1) << std::endl;
//    std::cout << sqws[i] << std::endl;
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
  unsigned int seed_n = 20;
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

  std::vector<Point_2> seeds(seed_n);
  std::vector<Vector5d> R5seeds(seed_n);
  std::vector<Metric> seeds_m(seed_n);
  std::vector<FT> sqws(seed_n);

  build_seeds(seeds, R5seeds, seeds_m, sqws, metric_field);

  for(unsigned int i=0; i<n; ++i)
  {
    for(unsigned int j=0; j<n ;++j)
    {
      //filling from bot left to top right
      Point_2 p(offset_x+j*step, offset_y+i*step);
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
      double min = 1e30;
      unsigned int min_id = 0;
      for(unsigned int i=0; i<seed_n; ++i)
      {
//        double sq_d = (p_hat - R5seeds[i]).norm();             // EUCLIDEAN DISTANCE IN R^5
//        double sq_d = std::pow((p_hat - R5seeds[i]).norm(),2.) - sqws[i];  // WEIGHTED DISTANCE IN R^5
        double sq_d = csd(tp, m_p.transform(seeds[i]));          // dx(p,x) < dx(q,x) (DU WANG)

        TPoint_2 tp2 = (seeds_m[i]).transform(p);
        TPoint_2 ts = (seeds_m[i]).transform(seeds[i]);
//        double sq_d = csd(tp2, ts);                             // dp(p,x) < dq(q,x) (LABELLE SHEWCHUK)
        if(sq_d < min)
        {
          min = sq_d;
          min_id = i;
        }
      }

      if(min < 10*step*step) // closer to a seed than another point of the grid (sort of)
        min_id = seed_n;

      out_bb << min_id*min_id * ((double) min_id+1) / (double) seed_n  << std::endl;
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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Constrain_surface_3_polyhedral.h>
#include <CGAL/Metric.h>
#include <CGAL/Starset.h>

#include <CGAL/Implicit_curvature_metric_field.h>
#include <CGAL/Polyhedral_curvature_metric_field.h>
#include <Metric_field/Custom_metric_field.h>
#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/External_metric_field.h>
#include <Metric_field/Hyperbolic_shock_metric_field.h>

#include <ostream>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
//typedef typename CGAL::Exact_predicates_exact_constructions_kernel KExact;
typedef K                                                    KExact;

typedef typename K::FT                                       FT;
typedef typename K::Point_3                                  Point_3;
typedef typename K::Point_3                                  TPoint_3;
typedef typename K::Vector_3                                 Vector_3;

typedef CGAL::Anisotropic_mesh_3::Metric_base<K>             Metric;
typedef CGAL::Anisotropic_mesh_3::Starset<K>                 Starset;
typedef typename Starset::Star_handle                        Star_handle;
typedef typename Starset::Traits                             Star_traits;

typedef CGAL::Anisotropic_mesh_3::Constrain_surface_3_polyhedral<K> Constrain_surface;

typedef CGAL::Anisotropic_mesh_3::Euclidean_metric_field<K>             Euclidean_metric_field;
typedef CGAL::Anisotropic_mesh_3::Hyperbolic_shock_metric_field<K>      Hyperbolic_shock_metric_field;
typedef CGAL::Anisotropic_mesh_3::Custom_metric_field<K>                Custom_metric_field;
typedef CGAL::Anisotropic_mesh_3::External_metric_field<K>              External_metric_field;
typedef CGAL::Anisotropic_mesh_3::Polyhedral_curvature_metric_field<K>  Polyhedral_curvature_metric_field;
typedef CGAL::Anisotropic_mesh_3::Implicit_curvature_metric_field<K>    Implicit_curvature_metric_field;

//mesh fetching
void fetch_from_off(std::ifstream& in,
                    std::vector<Point_3>& points,
                    std::vector<int>& facets)
{
  if(!in)
    std::cerr << "Error: cannot read file " << std::endl;

  int nv, nf;
  std::string useless;

  if(in >> useless && useless != "OFF")
    return;

  in >> nv >> nf >> useless;

  std::cout << "nvnf: " << nv << " " << nf << std::endl;

  for(int i=0; i<nv; ++i)
  {
    double x,y,z;
    in >> x >> y >> z;
    points.push_back(Point_3(x,y,z));
  }

  for(int i=0; i<nf; ++i)
  {
    int n1, n2, n3;
    in >> useless >> n1 >> n2  >> n3;
    facets.push_back(n1); facets.push_back(n2); facets.push_back(n3);
  }
}

void fetch_from_mesh(std::ifstream& in,
                     std::vector<Point_3>& points,
                     std::vector<int>& facets,
                     std::vector<int>& cells)
{
  if(!in)
    std::cerr << "Error: cannot read file " << std::endl;

  int nv, nf, nc;
  std::string useless, word;

  in >> useless >> useless; // MeshVersionFormatted 1
  in >> useless >> useless; // Dimension 3

  while(in >> word && word != "End")
  {
    if(word == "Vertices")
    {
      in >> nv;
      for(int i=0; i<nv; ++i)
      {
        double x,y,z;
        in >> x >> y >> z >> useless;
        points.push_back(Point_3(x,y,z));
      }
    }

    if(word == "Triangles")
    {
      in >> nf;
      for(int i=0; i<nf; ++i)
      {
        int n1, n2, n3;
        in >> n1 >> n2 >> n3 >> useless;
        facets.push_back(n1-1); facets.push_back(n2-1); facets.push_back(n3-1);
      }
    }

    if(word == "Tetrahedra")
    {
      in >> nc;
      for(int i=0; i<nc; ++i)
      {
        int n1, n2, n3, n4;
        in >> n1 >> n2 >> n3 >> n4 >> useless;
        cells.push_back(n1-1); cells.push_back(n2-1);
        cells.push_back(n3-1); cells.push_back(n4-1);
      }
    }
  }
}

void fetch_mesh(const char* filename,
                std::vector<Point_3>& points,
                std::vector<int>& facets,
                std::vector<int>& cells)
{
  std::cout << "fetching...: " << filename << std::endl;

  std::string input_filename = filename;
  std::string extension = input_filename.substr(input_filename.find_last_of('.'));

  if (extension == ".off" || extension == ".OFF")
  {
    std::ifstream stream(filename);
    if(!stream)
    {
      std::cerr << "Error: cannot read file " << std::endl;
      return;
    }
    fetch_from_off(stream, points, facets);
  }
  else if (extension == ".mesh" || extension == ".MESH")
  {
    std::ifstream stream(filename);
    if(!stream)
    {
      std::cerr << "Error: cannot read file " << std::endl;
      return;
    }
    fetch_from_mesh(stream, points, facets, cells);
  }

  std::cout << "Fetched mesh: " << std::endl;
  std::cout << points.size() << " points" << std::endl;
  std::cout << facets.size()/3 << " facets" << std::endl;
  std::cout << cells.size()/4 << " cells" << std::endl;
}

// histograms
void output_histogram(const std::vector<int>& histogram,
                      FT min, FT max,
                      const char* filename)
{
  std::cout << "output: " << filename << std::endl;
  std::ofstream out(filename);
  std::size_t histo_n = histogram.size();
  std::cout << "histo_n: " << histo_n << std::endl;
  for(std::size_t i=0; i<histo_n; ++i)
  {
    FT val = min + (max-min)*((FT) i)/((FT) histo_n);
    out << i << "," << val << "," << histogram[i] << std::endl;
  }
}

void output_histogram(const std::vector<FT>& values,
                      const char* filename)
{
#if 0//def FORCE_MIN_MAX_VALUES
  FT min_value = 0.30512;
  FT max_value = 2.57682;
#else
  FT min_value = *(std::min_element(values.begin(), values.end()));
  FT max_value = *(std::max_element(values.begin(), values.end()));
#endif
  std::cout << "Outputing: " << values.size() << " " << min_value << " " << max_value << std::endl;

  int histogram_size = 1000;
  std::vector<int> histogram(histogram_size, 0);
  FT limit_val = histogram_size - 1.;
  FT step_size = (max_value - min_value) / (FT) histogram_size;

  for(std::size_t i=0; i<values.size(); ++i)
    histogram[(std::min)(limit_val, std::floor((values[i]-min_value)/step_size))]++;

  output_histogram(histogram, min_value, max_value, filename);
}

template<typename Metric_field>
void compute_metrics(const std::vector<Point_3>& points,
                     const Metric_field* mf,
                     std::vector<Metric>& metrics)
{
  metrics.resize(points.size());
  for(std::size_t i=0; i<points.size(); ++i)
    metrics[i] = mf->compute_metric(points[i]);
}

void facet_distortion_histogram(const std::vector<Point_3>& points,
                                const std::vector<int>& facets,
                                const std::vector<Metric>& metrics)
{
  if(facets.empty())
    return;
  std::cout << "facet distortion external histo with size: " << facets.size() << std::endl;

  std::ofstream out_bb("max_distortion.bb");
  std::vector<FT> max_dists(points.size(), 1.);
  std::vector<FT> values;

  for(std::size_t i=0; i<facets.size();)
  {
    std::vector<int> ns(3);

    for(int j=0; j<3; ++j)
      ns[j] = facets[i++];

    for(int j=0; j<3; ++j)
    {
      const Metric& m = metrics[ns[j]];
      FT d1 = m.compute_distortion(metrics[ns[(j+1)%3]]);
      values.push_back(d1);
      FT d2 = m.compute_distortion(metrics[ns[(j+2)%3]]);
      values.push_back(d2);

      max_dists[ns[j]] = (std::max)((std::max)(max_dists[ns[j]], d1), d2);
    }
  }
  output_histogram(values, "histogram_facet_distortion_external.cvs");

  out_bb << "3 1 " << points.size() << " 2" << std::endl;
  for(std::size_t i=0; i<max_dists.size(); ++i)
    out_bb << max_dists[i] << std::endl;
}

void facet_edge_length_histogram(const std::vector<Point_3>& points,
                                 const std::vector<int>& facets,
                                 const std::vector<Metric>& metrics)
{
  if(facets.empty())
    return;
  std::cout << "facet edge external histo with size: " << facets.size() << std::endl;

  std::ofstream out_bb("max_edge_ratio.bb");
  std::vector<FT> max_edge_r(points.size(), 1.);
  std::vector<FT> values;

  Star_traits star_traits;
  typename Star_traits::Compute_squared_distance_3 csd = star_traits.compute_squared_distance_3_object();

  for(std::size_t i=0; i<facets.size();)
  {
    std::vector<int> ns(3);
    std::vector<Point_3> ps(3);

    for(int j=0; j<3; ++j)
    {
      ns[j] = facets[i++];
      ps[j] = points[ns[j]];
    }

    for(int j=0; j<3; ++j)
    {
      const Metric& m = metrics[ns[j]];
      std::vector<TPoint_3> tps(3);
      for(int k=0; k<3; ++k)
        tps[k] = m.transform(ps[k]);

#if 1 // adjacent only
      FT l1 = CGAL::sqrt(csd(tps[j], tps[(j+1)%3]));
      values.push_back(l1);
      FT l2 = CGAL::sqrt(csd(tps[j], tps[(j+2)%3]));
      values.push_back(l2);

      if(l1<1.) l1 = 1./l1;
      if(l2<1.) l2 = 1./l2;

      max_edge_r[ns[j]] = (std::max)((std::max)(max_edge_r[ns[j]], l1), l2);
#else
      values.push_back(CGAL::sqrt(csd(tps[0], tps[1])));
      values.push_back(CGAL::sqrt(csd(tps[0], tps[2])));
      values.push_back(CGAL::sqrt(csd(tps[1], tps[2])));
#endif
    }
  }
  output_histogram(values, "histogram_facet_edge_length_external.cvs");

  out_bb << "3 1 " << max_edge_r.size() << " 2" << std::endl;
  for(std::size_t i=0; i<max_edge_r.size(); ++i)
    out_bb << max_edge_r[i] << std::endl;
}

void facet_angle_histogram(const std::vector<Point_3>& points,
                           const std::vector<int>& facets,
                           const std::vector<Metric>& metrics)
{
  if(facets.empty())
    return;

  std::cout << "facet angle histo with size: " << facets.size() << std::endl;

  std::vector<FT> values;

  for(std::size_t i=0; i<facets.size();)
  {
    std::vector<int> ns(3);
    std::vector<Point_3> ps(3);

    for(int j=0; j<3; ++j)
    {
      ns[j] = facets[i++];
      ps[j] = points[ns[j]];
    }

    for(int j=0; j<3; ++j)
    {
      const Metric& m = metrics[ns[j]];
      std::vector<TPoint_3> tps(3);

      for(int k=0; k<3; ++k)
        tps[k] = m.transform(ps[(j+k)%3]);

      const Point_3& a = tps[0];
      const Point_3& b = tps[1];
      const Point_3& c = tps[2];

      Vector_3 v1 = b - a;
      Vector_3 v2 = c - a;
      FT angle_cos = v1 * v2 / CGAL::sqrt(v1*v1) / CGAL::sqrt(v2 * v2);
      FT angle = std::acos(angle_cos);

      // make sure it's within 0 pi
      if(angle < 0)
        angle += CGAL_PI;

      std::cout << "angle: " << angle << std::endl;

      CGAL_assertion(angle >=0 && angle <= CGAL_PI);

      values.push_back(angle);
    }
  }
  std::cout << "n of values: " << values.size() << std::endl;

  output_histogram(values, "histogram_face_angles_external.cvs");
}

void cell_distortion_histogram(const std::vector<Point_3>& points,
                               const std::vector<int>& cells,
                               const std::vector<Metric>& metrics)
{
  if(cells.empty())
    return;
  std::cout << "cell distortion external histo with size: " << cells.size() << std::endl;

  std::ofstream out_bb("max_distortion.bb");
  std::vector<FT> max_dists(points.size(), 1.);
  std::vector<FT> values;

  for(std::size_t i=0; i<cells.size();)
  {
    std::vector<int> ns(4);

    for(int j=0; j<4; ++j)
      ns[j] = cells[i++];

    for(int j=0; j<4; ++j)
    {
      const Metric& m = metrics[ns[j]];
      FT d1 = m.compute_distortion(metrics[ns[(j+1)%4]]);
      values.push_back(d1);
      FT d2 = m.compute_distortion(metrics[ns[(j+2)%4]]);
      values.push_back(d2);
      FT d3 = m.compute_distortion(metrics[ns[(j+3)%4]]);
      values.push_back(d3);

      max_dists[ns[j]] = (std::max)((std::max)((std::max)(max_dists[ns[j]], d1), d2), d3);
    }
  }
  output_histogram(values, "histogram_cell_distortion_external.cvs");

  out_bb << "3 1 " << points.size() << " 2" << std::endl;
  for(std::size_t i=0; i<max_dists.size(); ++i)
    out_bb << max_dists[i] << std::endl;
}

void cell_edge_length_histogram(const std::vector<Point_3>& points,
                                const std::vector<int>& cells,
                                const std::vector<Metric>& metrics)
{
  if(cells.empty())
    return;
  std::cout << "cell edge length external histo with size: " << cells.size() << std::endl;

  std::ofstream out_bb("max_edge_ratio.bb");
  std::vector<FT> max_edge_r(points.size(), 1.);
  std::vector<FT> values;

  Star_traits star_traits;
  typename Star_traits::Compute_squared_distance_3 csd = star_traits.compute_squared_distance_3_object();

  for(std::size_t i=0; i<cells.size();)
  {
    std::vector<int> ns(4);
    std::vector<Point_3> ps(4);

    for(int j=0; j<4; ++j)
    {
      ns[j] = cells[i++];
      ps[j] = points[ns[j]];
    }

    for(int j=0; j<4; ++j)
    {
      const Metric& m = metrics[ns[j]];
      std::vector<TPoint_3> tps(4);
      for(int k=0; k<4; ++k)
        tps[k] = m.transform(ps[k]);

#if 1 // adjacent only
      FT l1 = CGAL::sqrt(csd(tps[j], tps[(j+1)%4]));
      values.push_back(l1);
      FT l2 = CGAL::sqrt(csd(tps[j], tps[(j+2)%4]));
      values.push_back(l2);
      FT l3 = CGAL::sqrt(csd(tps[j], tps[(j+3)%4]));
      values.push_back(l3);

      if(l1<1.) l1 = 1./l1;
      if(l2<1.) l2 = 1./l2;
      if(l3<1.) l3 = 1./l3;

      max_edge_r[ns[j]] = (std::max)((std::max)((std::max)(max_edge_r[ns[j]], l1), l2), l3);
#else
      values.push_back(CGAL::sqrt(csd(tps[0], tps[1])));
      values.push_back(CGAL::sqrt(csd(tps[0], tps[2])));
      values.push_back(CGAL::sqrt(csd(tps[0], tps[3])));
      values.push_back(CGAL::sqrt(csd(tps[1], tps[2])));
      values.push_back(CGAL::sqrt(csd(tps[1], tps[3])));
      values.push_back(CGAL::sqrt(csd(tps[2], tps[3])));
#endif
    }
  }
  output_histogram(values, "histogram_cell_edge_length_external.cvs");

  out_bb << "3 1 " << max_edge_r.size() << " 2" << std::endl;
  for(std::size_t i=0; i<max_edge_r.size(); ++i)
    out_bb << max_edge_r[i] << std::endl;
}

template<typename Metric_field>
void cell_edge_length_histogram_midpoint_metric(const std::vector<Point_3>& points,
                                                const std::vector<int>& cells,
                                                const Metric_field* mf)
{
  if(cells.empty())
    return;
  std::cout << "cell edge length external midpoint metric histo with size: " << cells.size() << std::endl;

  std::vector<FT> values;

  Star_traits star_traits;
  typename Star_traits::Compute_squared_distance_3 csd = star_traits.compute_squared_distance_3_object();

  for(std::size_t i=0; i<cells.size();)
  {
    std::vector<int> ns(4);
    std::vector<Point_3> ps(4);

    for(int j=0; j<4; ++j)
    {
      ns[j] = cells[i++];
      ps[j] = points[ns[j]];
    }

    for(int j=0; j<3; ++j)
    {
      for(int k=j+1; k<4; ++k)
      {
        Point_3 mid_point(0.5*(ps[j].x() + ps[k].x()),
                          0.5*(ps[j].y() + ps[k].y()),
                          0.5*(ps[j].z() + ps[k].z()));
        Metric m = mf->compute_metric(mid_point);

        TPoint_3 tpj = m.transform(ps[j]);
        TPoint_3 tpk = m.transform(ps[k]);

        values.push_back(CGAL::sqrt(csd(tpj, tpk)));
      }
    }
  }
  output_histogram(values, "histogram_cell_edge_length_external_midpoint_metric.cvs");
}

template<typename Metric_field>
void cell_edge_length_histogram_simplex_metric(const std::vector<Point_3>& points,
                                                const std::vector<int>& cells,
                                                const Metric_field* mf)
{
  if(cells.empty())
    return;
  std::cout << "cell edge length external midpoint metric histo with size: " << cells.size() << std::endl;

  std::vector<FT> values;

  Star_traits star_traits;
  typename Star_traits::Compute_squared_distance_3 csd = star_traits.compute_squared_distance_3_object();

  for(std::size_t i=0; i<cells.size();)
  {
    std::vector<int> ns(4);
    std::vector<Point_3> ps(4);

    for(int j=0; j<4; ++j)
    {
      ns[j] = cells[i++];
      ps[j] = points[ns[j]];
    }

    Point_3 simplex_center = CGAL::barycenter(ps[0], 0.25, ps[1], 0.25, ps[2], 0.25, ps[3], 025);
    Metric simplex_metric = mf->compute_metric(simplex_center);

    std::vector<Point_3> tps(4);
    for(int j=0; j<4; ++j)
      tps[j] = simplex_metric.transform(ps[j]);

    for(int j=0; j<3; ++j)
      for(int k=j+1; k<4; ++k)
        values.push_back(CGAL::sqrt(csd(tps[j], tps[k])));
  }
  output_histogram(values, "histogram_cell_edge_length_external_simplex_metric.cvs");
}

int main(int, char**)
{
  std::freopen("log.txt", "w", stdout); //all output is written in "log.txt"

  std::vector<Point_3> points;
  std::vector<int> facets, cells;
  std::vector<Metric> metrics;

//  Constrain_surface* pdomain = new Constrain_surface("bambimboum.mesh");
  Constrain_surface* pdomain = new Constrain_surface("../../data/Anisotropy_CMP/3DSurface/Fandisk.off");

  //----------- pick a metric field! ----
//  Euclidean_metric_field* metric_field = new Euclidean_metric_field();
//  Hyperbolic_shock_metric_field* metric_field = new Hyperbolic_shock_metric_field(0.6);
  Custom_metric_field* metric_field = new Custom_metric_field();
//  External_metric_field* metric_field = new External_metric_field(*pdomain, "../../data/Anisotropy_CMP/3DSurface/Fandisk_Metric.txt");
//  External_metric_field* metric_field = new External_metric_field(*pdomain, "metric_at_input.txt");
//  Polyhedral_curvature_metric_field* metric_field = new Polyhedral_curvature_metric_field(*pdomain);
//  Implicit_curvature_metric_field* metric_field = new Implicit_curvature_metric_field(*pdomain);

//  const char* mesh_filename = "/home/mrouxell/cgal/Anisotropic_mesh_3/data/Anisotropy_CMP/3DSurface/Fandisk.off";
//  const char* mesh_filename = "../../data/Anisotropy_CMP/Ours_Results/Tet/Spherical.mesh";
//  const char* mesh_filename = "./experiments/Yang Liu/cube/smoothed/5690/smoothed_cube_5690.mesh";
//  const char* mesh_filename = "../../data/Anisotropy_CMP/Ours_Results/Fandisk_Ours/our_metric/fandisk_7270_feature.mesh";
  const char* mesh_filename = "bambimboum_surf.mesh";
//  const char* mesh_filename = "/home/mrouxell/vrac/vrac.off";
//  const char* mesh_filename = "./experiments/Yang Liu/fandisk/iso/3D/Aniso_3/fandisk_100968.mesh";
//  const char* mesh_filename = "./experiments/Yang Liu/Sphere/sphere_interrupted.mesh";
//  const char* mesh_filename = "/home/mrouxell/Reviews/SIGASIA_Local_Convex_Aniso/AnisoMeshData/VOLUME/fig14_sine.mesh";
//  const char* mesh_filename = "/home/mrouxell/Reviews/SIGASIA_Local_Convex_Aniso/AnisoMeshData/SURFACE/block.off";

  fetch_mesh(mesh_filename, points, facets, cells);
  compute_metrics(points, metric_field, metrics);

//  facet_distortion_histogram(points, facets, metrics);
//  facet_edge_length_histogram(points, facets, metrics);
  facet_angle_histogram(points, facets, metrics);
//  cell_distortion_histogram(points, cells, metrics);
//  cell_edge_length_histogram(points, cells, metrics);
//  cell_edge_length_histogram_midpoint_metric(points, cells, metric_field);
//  cell_edge_length_histogram_simplex_metric(points, cells, metric_field);

//  delete pdomain;
  delete metric_field;

  return 0;
}


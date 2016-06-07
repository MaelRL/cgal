#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Metric.h>
#include <CGAL/Starset.h>

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <ostream>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
//typedef CGAL::Exact_predicates_exact_constructions_kernel KExact;
typedef K                                                    KExact;

typedef typename K::FT                                       FT;
typedef typename K::Point_2                                  Point_2;
typedef typename K::Point_2                                  TPoint_2;
typedef typename K::Vector_2                                 Vector_2;

typedef CGAL::Anisotropic_mesh_2::Metric_base<K>             Metric;
typedef CGAL::Anisotropic_mesh_2::Starset<K>                 Starset;
typedef typename Starset::Star_handle                        Star_handle;
typedef typename Starset::Traits                             Star_traits;

typedef CGAL::Anisotropic_mesh_2::Euclidean_metric_field<K>  Euclidean_metric_field;
typedef CGAL::Anisotropic_mesh_2::Custom_metric_field<K>     Custom_metric_field;

//mesh fetching
void fetch_from_off(std::ifstream& in,
                    std::vector<Point_2>& points,
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
    double x,y;
    in >> x >> y;
    points.push_back(Point_2(x,y));
  }

  for(int i=0; i<nf; ++i)
  {
    int n1, n2, n3;
    in >> useless >> n1 >> n2  >> n3;
    facets.push_back(n1); facets.push_back(n2); facets.push_back(n3);
  }
}

void fetch_from_mesh(std::ifstream& in,
                     std::vector<Point_2>& points,
                     std::vector<int>& facets)
{
  if(!in)
    std::cerr << "Error: cannot read file " << std::endl;

  int nv, nf;
  std::string useless, word;

  in >> useless >> useless; // MeshVersionFormatted 1
  in >> useless >> useless; // Dimension 2

  while(in >> word && word != "End")
  {
    if(word == "Vertices")
    {
      in >> nv;
      for(int i=0; i<nv; ++i)
      {
        double x,y;
        in >> x >> y >> useless;
        points.push_back(Point_2(x,y));
      }
    }

    if(word == "Triangles")
    {
      in >> nf;
      for(int i=0; i<nf; ++i)
      {
        int n1, n2, n3;
        in >> n1 >> n2 >> n3 >> useless;

#if 0 // artificial filter
        FT minv = 1, maxv = 1;
        Point_2 p(0,0,0);
        if(CGAL::squared_distance(points[n1-1], p) > maxv ||
           CGAL::squared_distance(points[n1-1], p) < minv ||
           CGAL::squared_distance(points[n2-1], p) > maxv ||
           CGAL::squared_distance(points[n2-1], p) < minv ||
           CGAL::squared_distance(points[n3-1], p) > maxv ||
           CGAL::squared_distance(points[n3-1], p) < minv)
          continue;
#endif
        facets.push_back(n1-1); facets.push_back(n2-1); facets.push_back(n3-1);
      }
    }

    if(word == "Edges")
    {
      std::cerr << "can't handle edges!" << std::endl;
      return;
    }
  }
}

void fetch_mesh(const std::string filename,
                std::vector<Point_2>& points,
                std::vector<int>& facets)
{
  std::cout << "fetching...: " << filename << std::endl;

  std::string extension = filename.substr(filename.find_last_of('.'));

  if (extension == ".off" || extension == ".OFF")
  {
    std::ifstream stream(filename.c_str());
    if(!stream)
    {
      std::cerr << "Error: cannot read file " << std::endl;
      return;
    }
    fetch_from_off(stream, points, facets);
  }
  else if (extension == ".mesh" || extension == ".MESH")
  {
    std::ifstream stream(filename.c_str());
    if(!stream)
    {
      std::cerr << "Error: cannot read file " << std::endl;
      return;
    }
    fetch_from_mesh(stream, points, facets);
  }

  std::cout << "Fetched mesh: " << std::endl;
  std::cout << points.size() << " points" << std::endl;
  std::cout << facets.size()/3 << " facets" << std::endl;
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
    out << val << "," << histogram[i] << std::endl;
  }
}

void output_histogram(std::vector<FT>& values,
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

  int histogram_size = 100;
  std::vector<int> histogram(histogram_size, 0);
  FT limit_val = histogram_size - 1.;
  FT step_size = (max_value - min_value) / (FT) histogram_size;

  for(std::size_t i=0; i<values.size(); ++i)
    histogram[(std::min)(limit_val, std::floor((values[i]-min_value)/step_size))]++;

  output_histogram(histogram, min_value, max_value, filename);
}

template<typename Metric_field>
void compute_metrics(const std::vector<Point_2>& points,
                     const Metric_field* mf,
                     std::vector<Metric>& metrics)
{
  metrics.resize(points.size());
  for(std::size_t i=0; i<points.size(); ++i)
    metrics[i] = mf->compute_metric(points[i]);
}

void face_distortion_histogram(const std::vector<Point_2>& points,
                               const std::vector<int>& faces,
                               const std::vector<Metric>& metrics,
                               const std::string filename_core)
{
  if(faces.empty())
    return;

  std::cout << "face distortion external histo with size: " << faces.size() << std::endl;

  std::ofstream out_bb((filename_core + "_distortion.bb").c_str());
  std::vector<FT> max_dists(points.size(), 1.);
  std::vector<FT> values;

  for(std::size_t i=0; i<faces.size();)
  {
    std::vector<int> ns(3);

    for(int j=0; j<3; ++j)
      ns[j] = faces[i++];

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
  output_histogram(values, "histogram_face_distortion_external.csv");

  out_bb << "2 1 " << points.size() << " 2" << std::endl;
  for(std::size_t i=0; i<max_dists.size(); ++i)
    out_bb << max_dists[i] << std::endl;
}

FT element_quality(const Metric& m,
                   const Point_2& p1,
                   const Point_2& p2,
                   const Point_2& p3)
{
  TPoint_2 tp1 = m.transform(p1);
  TPoint_2 tp2 = m.transform(p2);
  TPoint_2 tp3 = m.transform(p3);

  Star_traits star_traits;
  typename Star_traits::Compute_squared_distance_2 sqd =
      star_traits.compute_squared_distance_2_object();
  FT alpha = 4.*std::sqrt(3.);
  FT A = std::abs(K::Compute_area_2()(tp1, tp2, tp3));
  FT sq_a = sqd(tp1, tp2);
  FT sq_b = sqd(tp1, tp3);
  FT sq_c = sqd(tp2, tp3);

//Frey
  //FT quality = alpha*A/(sq_a+sq_b+sq_c);
//Zhong
  FT a = CGAL::sqrt(sq_a);
  FT b = CGAL::sqrt(sq_b);
  FT c = CGAL::sqrt(sq_c);
  FT h = (std::max)((std::max)(a,b),c);
  FT p = a+b+c;
  FT quality = alpha*A/(p*h);

  return quality;
}

void face_quality_histogram(const std::vector<Point_2>& points,
                            const std::vector<int>& faces,
                            const std::vector<Metric>& metrics,
                            const std::string filename_core)
{
  if(faces.empty())
    return;

  std::cout << "face quality external histo with size: " << faces.size() << std::endl;

  std::ofstream out_bb((filename_core + "_quality.bb").c_str());
  std::vector<FT> min_quals(faces.size()/3, 1.);
  std::vector<FT> values;

  for(std::size_t i=0; i<faces.size();)
  {
    int face_n = i/3;
    std::vector<int> ns(3);

    for(int j=0; j<3; ++j)
      ns[j] = faces[i++];

    if(points[ns[0]].x() > 2. || points[ns[0]].x() < -2. || points[ns[0]].y() > 2. || points[ns[0]].y() < -2. ||
       points[ns[1]].x() > 2. || points[ns[1]].x() < -2. || points[ns[1]].y() > 2. || points[ns[1]].y() < -2. ||
       points[ns[2]].x() > 2. || points[ns[2]].x() < -2. || points[ns[2]].y() > 2. || points[ns[2]].y() < -2.)
      continue;

    for(int j=0; j<3; ++j)
    {
      const Metric& m = metrics[ns[j]];
      FT q = element_quality(m, points[ns[0]], points[ns[1]], points[ns[2]]);
      values.push_back(q);
      min_quals[face_n] = (std::min)(min_quals[face_n], q);
    }
  }
  output_histogram(values, "histogram_face_quality_external.csv");

  out_bb << "2 1 " << faces.size() << " 1" << std::endl;
  for(std::size_t i=0; i<min_quals.size(); ++i)
    out_bb << min_quals[i] << std::endl;
}

void face_edge_length_histogram(const std::vector<Point_2>& points,
                                const std::vector<int>& faces,
                                const std::vector<Metric>& metrics,
                                const std::string filename_core)
{
  if(faces.empty())
    return;

  std::cout << "face edge external histo with size: " << faces.size() << std::endl;

  std::ofstream out_bb((filename_core + "_edge_length.bb").c_str());
  std::vector<FT> max_edge_r(points.size(), 1.);
  std::vector<FT> values;

  Star_traits star_traits;
  typename Star_traits::Compute_squared_distance_2 csd =
      star_traits.compute_squared_distance_2_object();

  for(std::size_t i=0; i<faces.size();)
  {
    std::vector<int> ns(3);
    std::vector<Point_2> ps(3);

    for(int j=0; j<3; ++j)
    {
      ns[j] = faces[i++];
      ps[j] = points[ns[j]];
    }

    for(int j=0; j<3; ++j)
    {
      const Metric& m = metrics[ns[j]];
      std::vector<TPoint_2> tps(3);
      for(int k=0; k<3; ++k)
        tps[k] = m.transform(ps[k]);

#if 1 // adjacent only
      FT l1 = CGAL::sqrt(csd(tps[j], tps[(j+1)%3]));
      FT l2 = CGAL::sqrt(csd(tps[j], tps[(j+2)%3]));

      FT r0 = 1.0;

      if(l1<r0) l1 = r0/l1;
      if(l2<r0) l2 = r0/l2;
      if(l1>10) l1 = 10;
      if(l2>10) l2 = 10;

      values.push_back(l1);
      values.push_back(l2);

      max_edge_r[ns[j]] = (std::max)((std::max)(max_edge_r[ns[j]], l1), l2);
#else
      values.push_back(CGAL::sqrt(csd(tps[0], tps[1])));
      values.push_back(CGAL::sqrt(csd(tps[0], tps[2])));
      values.push_back(CGAL::sqrt(csd(tps[1], tps[2])));
#endif
    }
  }
  output_histogram(values, "histogram_face_edge_length_external.csv");

  out_bb << "2 1 " << points.size() << " 2" << std::endl;
  for(std::size_t i=0; i<max_edge_r.size(); ++i)
    out_bb << max_edge_r[i] << std::endl;
}

void face_angle_histogram(const std::vector<Point_2>& points,
                                const std::vector<int>& faces,
                                const std::vector<Metric>& metrics,
                                const std::string filename_core)
{
  if(faces.empty())
    return;

  std::cout << "face edge external histo with size: " << faces.size() << std::endl;

  std::ofstream out_bb((filename_core + "_angles.bb").c_str());
  std::vector<FT> worst_face_angles(faces.size()/3);

  std::vector<FT> values;

  for(std::size_t i=0; i<faces.size();)
  {
    std::vector<int> ns(3);
    std::vector<Point_2> ps(3);

    FT worst_angle = 1e30;
    std::size_t entry = i/3;
//    std::cout << "e: " << entry << std::endl;

    for(int j=0; j<3; ++j)
    {
      ns[j] = faces[i++];
      ps[j] = points[ns[j]];
    }

    for(int j=0; j<3; ++j)
    {
      const Metric& m = metrics[ns[j]];
      std::vector<TPoint_2> tps(3);

      for(int k=0; k<3; ++k)
        tps[k] = m.transform(ps[(j+k)%3]);

      const Point_2& a = tps[0];
      const Point_2& b = tps[1];
      const Point_2& c = tps[2];

      Vector_2 v1 = b - a;
      Vector_2 v2 = c - a;
      FT angle_cos = v1 * v2 / CGAL::sqrt(v1 * v1) / CGAL::sqrt(v2 * v2);
      FT angle = std::acos(angle_cos);

      // make sure it's within 0 pi
      if(angle < 0)
        angle += CGAL_PI;

//      std::cout << v1 << " " << v2 << std::endl;
//      std::cout << "angle: " << angle << std::endl;

      CGAL_assertion(angle >=0 && angle <= CGAL_PI);
      worst_angle = (std::min)(worst_angle, (std::min)(angle, CGAL_PI - angle));

      values.push_back(angle);
    }

//    std::cout << "worst angle: " << worst_angle << std::endl;
    worst_face_angles[entry] = worst_angle;
  }

  std::cout << "n of values: " << values.size() << std::endl;
  output_histogram(values, "histogram_face_angles_external.csv");

  out_bb << "2 1 " << faces.size()/3 << " 1" << std::endl;
  for(std::size_t i=0; i<faces.size(); ++i)
    out_bb << worst_face_angles[i] << std::endl;
}

int main(int, char**)
{
//  std::freopen("wut.txt", "w", stdout); //all output is written in "wut.txt"

  std::vector<Point_2> points;
  std::vector<int> faces;
  std::vector<Metric> metrics;

  //----------- pick a metric field! ----
  //  Euclidean_metric_field* metric_field = new Euclidean_metric_field();
  Custom_metric_field* metric_field = new Custom_metric_field();

  const std::string filename_core = "super_dense_base_mesh_tr_dual";
  // ------------ pick OFF or MESH ------
  const std::string mesh_filename = filename_core + ".mesh";
  fetch_mesh(mesh_filename, points, faces);
//  const std::string off_filemane = filename_core + ".off";
//  fetch_off(mesh_filename, points, faces);

  compute_metrics(points, metric_field, metrics);

  face_distortion_histogram(points, faces, metrics, filename_core);
  face_quality_histogram(points, faces, metrics, filename_core);
  face_edge_length_histogram(points, faces, metrics, filename_core);
  face_angle_histogram(points, faces, metrics, filename_core);

  delete metric_field;

  return 0;
}


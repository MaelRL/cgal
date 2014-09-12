#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Metric.h>
#include <CGAL/Starset.h>

#include <Metric_field/Euclidean_metric_field.h>
#include <Metric_field/Custom_metric_field.h>

#include <ostream>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
//typedef typename CGAL::Exact_predicates_exact_constructions_kernel KExact;
typedef K                                                    KExact;

typedef typename K::FT                                       FT;
typedef typename K::Point_2                                  Point_2;
typedef typename K::Point_2                                  TPoint_2;

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

  int nv, nf, nc;
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
  }
}

void fetch_mesh(const char* filename,
                std::vector<Point_2>& points,
                std::vector<int>& facets)
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
void compute_metrics(const std::vector<Point_2>& points,
                     const Metric_field* mf,
                     std::vector<Metric>& metrics)
{
  metrics.resize(points.size());
  for(std::size_t i=0; i<points.size(); ++i)
    metrics[i] = mf->compute_metric(points[i]);
}

void facet_distortion_histogram(const std::vector<Point_2>& points,
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

void facet_edge_length_histogram(const std::vector<Point_2>& points,
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
  typename Star_traits::Compute_squared_distance_2 csd = star_traits.compute_squared_distance_2_object();

  for(std::size_t i=0; i<facets.size();)
  {
    std::vector<int> ns(3);
    std::vector<Point_2> ps(3);

    for(int j=0; j<3; ++j)
    {
      ns[j] = facets[i++];
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

int main(int, char**)
{
  std::freopen("wutt.txt", "w", stdout); //all output is written in "wut.txt"

  std::vector<Point_2> points;
  std::vector<int> facets;
  std::vector<Metric> metrics;

  //----------- pick a metric field! ----
//  Euclidean_metric_field* metric_field = new Euclidean_metric_field();
  Custom_metric_field* metric_field = new Custom_metric_field();

  const char* mesh_filename = "bambimboum.mesh";

  fetch_mesh(mesh_filename, points, facets);
  compute_metrics(points, metric_field, metrics);

  facet_distortion_histogram(points, facets, metrics);
  facet_edge_length_histogram(points, facets, metrics);

  delete metric_field;

  return 0;
}


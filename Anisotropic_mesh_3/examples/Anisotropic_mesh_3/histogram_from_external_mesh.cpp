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
#ifdef FORCE_MIN_MAX_VALUES
  FT min_value = 0.30512;
  FT max_value = 2.57682;
#else
  FT min_value = *(std::min_element(values.begin(), values.end()));
  FT max_value = *(std::max_element(values.begin(), values.end()));
  FT avg = std::accumulate( values.begin(), values.end(), 0.0) / values.size();
#endif
  std::cout << "Outputing: " << values.size() << " values." << std::endl;
  std::cout << "min: " << min_value 
            << " max: " << max_value
            << " avg: " << avg << std::endl;

  int histogram_size = 100;
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
  std::cout << "\n -- Facet distortion histo with size: " << facets.size() << std::endl;

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

FT element_quality(const Metric& m,
                   const Point_3& p1,
                   const Point_3& p2,
                   const Point_3& p3)
{
  TPoint_3 tp1 = m.transform(p1);
  TPoint_3 tp2 = m.transform(p2);
  TPoint_3 tp3 = m.transform(p3);

  Star_traits star_traits;
  typename Star_traits::Compute_squared_distance_3 sqd =
      star_traits.compute_squared_distance_3_object();
  FT alpha = 4.*std::sqrt(3.);
  FT A = std::abs(K::Compute_area_3()(tp1, tp2, tp3));
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

void facet_quality_histogram(const std::vector<Point_3>& points,
                             const std::vector<int>& facets,
                             const std::vector<Metric>& metrics)
{
  if(facets.empty())
    return;

  std::cout << "\n -- Facet Quality histo with size: " << facets.size() << std::endl;

  std::ofstream out_bb("facet_quality.bb");
  std::vector<FT> min_quals(facets.size()/3, 1.);
  std::vector<FT> values;

  for(std::size_t i=0; i<facets.size();)
  {
    int face_n = i/3;
    std::vector<int> ns(3);

    for(int j=0; j<3; ++j)
      ns[j] = facets[i++];

    for(int j=0; j<3; ++j)
    {
      const Metric& m = metrics[ns[j]];
      FT q = element_quality(m, points[ns[0]], points[ns[1]], points[ns[2]]);
      values.push_back(q);
      min_quals[face_n] = (std::min)(min_quals[face_n], q);
    }
  }
  output_histogram(values, "histogram_facet_quality_external.csv");

  out_bb << "2 1 " << facets.size() << " 1" << std::endl;
  for(std::size_t i=0; i<min_quals.size(); ++i)
    out_bb << min_quals[i] << std::endl;
}


void facet_edge_length_histogram(const std::vector<Point_3>& points,
                                 const std::vector<int>& facets,
                                 const std::vector<Metric>& metrics)
{
  if(facets.empty())
    return;
  std::cout << "\n -- Edge length histo with size: " << facets.size() << std::endl;

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

  std::cout << "\n --  Facet angle histo with size: " << facets.size() << std::endl;

  std::vector<FT> values;

  for(std::size_t i=0; i<facets.size();)
  {
    std::vector<int> ns(3);
    std::vector<Point_3> ps(3);
    FT worst_angle = 1e30;

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

      worst_angle = (std::min)(worst_angle, (std::min)(angle, CGAL_PI - angle));

//      std::cout << "angle: " << angle << std::endl;


      CGAL_assertion(angle >=0 && angle <= CGAL_PI);
    }

    // convert to radius
    worst_angle = worst_angle / CGAL_PI * 180.;

    values.push_back(worst_angle);
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
  std::cout << "\n -- Cell distortion histo with size: " << cells.size() << std::endl;

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
  std::cout << "\n -- Cell edge length histo with size: " << cells.size() << std::endl;

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
  std::cout << "\n --  Cell edge length midpoint metric histo with size: " << cells.size() << std::endl;

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
  std::cout << "\n -- Cell edge length simplex metric histo with size: " << cells.size() << std::endl;

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

FT minimum_dihedral_angle(const Metric& metric,
                          const Point_3& p0,
                          const Point_3& p1,
                          const Point_3& p2,
                          const Point_3& p3)
{
  const Point_3& tp0 = metric.transform(p0);
  const Point_3& tp1 = metric.transform(p1);
  const Point_3& tp2 = metric.transform(p2);
  const Point_3& tp3 = metric.transform(p3);

  Star_traits star_traits;
  typename Star_traits::Compute_squared_distance_3 sqd =
      star_traits.compute_squared_distance_3_object();

  typename Star_traits::Compute_determinant_3 determinant =
      star_traits.compute_determinant_3_object();
  typename Star_traits::Construct_cross_product_vector_3 cp =
      star_traits.construct_cross_product_vector_3_object();
  typename Star_traits::Compute_scalar_product_3 sp =
      star_traits.compute_scalar_product_3_object();

  Vector_3 v01 = tp1-tp0;
  Vector_3 v02 = tp2-tp0;
  Vector_3 v03 = tp3-tp0;
  Vector_3 v12 = tp2-tp1;
  Vector_3 v13 = tp3-tp1;
  Vector_3 v23 = tp3-tp2;

  Vector_3 v_01_02 = cp(v01,v02);
  FT a_012 = v_01_02*v_01_02;

  Vector_3 v_01_03 = cp(v01,v03);
  FT a_013 = v_01_03*v_01_03;

  Vector_3 v_12_13 = cp(v12,v13);
  FT a_123 = v_12_13*v_12_13;

  Vector_3 v_02_03 = cp(v02,v03);
  FT a_023 = v_02_03*v_02_03;

  FT min_quotient = sp(v01,v01) / (a_012 * a_013);
  min_quotient = (CGAL::min)(min_quotient,
                             sp(v02,v02) / (a_012 * a_023));
  min_quotient = (CGAL::min)(min_quotient,
                             (v03*v03) / (a_013 * a_023));
  min_quotient = (CGAL::min)(min_quotient,
                             sp(v12,v12) / (a_012 * a_123));
  min_quotient = (CGAL::min)(min_quotient,
                             sp(v13,v13) / (a_013 * a_123));
  min_quotient = (CGAL::min)(min_quotient,
                             sp(v23,v23) / (a_023 * a_123));
  min_quotient =  sqrt(min_quotient);

  return CGAL::abs(std::asin(determinant(v01, v02, v03) * min_quotient)
                   * FT(180) / FT(CGAL_PI));
}

FT compute_squared_shortest_edge(const Point_3& p,
                                 const Point_3& q,
                                 const Point_3& r,
                                 const Point_3& s)
{
  typename Star_traits::Compute_squared_distance_3 o;
  return min(min(min(min(min(o(p, q), o(p, r)), o(p, s)), o(q, r)), o(q, s)), o(r, s));
}

FT compute_sliverity(const Metric& metric,
                     const Point_3& p0,
                     const Point_3& p1,
                     const Point_3& p2,
                     const Point_3& p3)
{
  // not very efficient, but heh...
  const Point_3& p = metric.transform(p0);
  const Point_3& q = metric.transform(p1);
  const Point_3& r = metric.transform(p2);
  const Point_3& s = metric.transform(p3);

  typename Star_traits::Compute_volume_3 o;

  FT V_third = std::pow(std::abs(o(p, q, r, s)), 1.0 / 3.0);
  FT h = std::sqrt(compute_squared_shortest_edge(p, q, r, s));

  return V_third / h;
}

template<typename Metric_field>
void cell_dihedral_angle_histogram(const std::vector<Point_3>& points,
                                   const std::vector<int>& cells,
                                   const Metric_field* mf)
{
  if(cells.empty())
    return;
  std::cout << "\n -- Cell dihedral angle histo with size: " << cells.size() << std::endl;

  std::vector<FT> values;

  for(std::size_t i=0; i<cells.size();)
  {
    std::vector<int> ns(4);
    std::vector<Point_3> ps(4);

    FT smallest_dihedral_angle = 1e30;

    for(int j=0; j<4; ++j)
    {
      ns[j] = cells[i++];
      ps[j] = points[ns[j]];
    }

    for(int j=0; j<4; ++j)
    {
      Metric m = mf->compute_metric(ps[j]);
      FT diha = minimum_dihedral_angle(m, ps[0], ps[1], ps[2], ps[3]);

/*
      if(diha < 10)
      {
        std::cout << "small dihedral angle. Sliverity: ";
        std::cout << compute_sliverity(m, ps[0], ps[1], ps[2], ps[3]) << std::endl;
      }
*/

      if(diha < smallest_dihedral_angle)
        smallest_dihedral_angle = diha;
    }

    values.push_back(smallest_dihedral_angle);
  }

  output_histogram(values, "histogram_cell_dihedral_angle_external_simplex_metric.cvs");
}

FT element_quality_3D(const Point_3& p, const Point_3& q,
                      const Point_3& r, const Point_3& s)
{
  Star_traits star_traits;
  typename Star_traits::Compute_area_3 ar;
  typename Star_traits::Compute_volume_3 vo;

  FT alpha = 216.*std::sqrt(3.);
  FT V = std::abs(vo(p, q, r, s));
  FT a = std::abs(ar(p, q, r));
  FT b = std::abs(ar(p, q, s));
  FT c = std::abs(ar(p, r, s));
  FT d = std::abs(ar(q, r, s));

  FT denom = (a + b + c + d);
  FT quality = alpha * V * V / (denom * denom * denom);

  if(quality > 1e30)
    std::cout << V << " " << a << " " << b << " " << c << " " << d << std::endl;

  return quality;
}


template<typename Metric_field>
void cell_quality_histogram(const std::vector<Point_3>& points,
                            const std::vector<int>& cells,
                            const Metric_field* mf)
{
  if(cells.empty())
    return;
  std::cout << "\n -- Cell quality histo with size: " << cells.size() << std::endl;

  std::vector<FT> values;

  for(std::size_t i=0; i<cells.size();)
  {
    std::vector<int> ns(4);
    std::vector<Point_3> ps(4);

    FT min_quality = 1e30;

    for(int j=0; j<4; ++j)
    {
      ns[j] = cells[i++];
      ps[j] = points[ns[j]];
    }

    for(int j=0; j<4; ++j)
    {
      Metric m = mf->compute_metric(ps[j]);
      std::vector<TPoint_3> tps(4);
      for(int k=0; k<4; ++k)
        tps[k] = m.transform(ps[k]);

      FT quality = element_quality_3D(tps[0], tps[1], tps[2], tps[3]);
      values.push_back(quality);

      if(quality < min_quality)
        min_quality = quality;
    }

//    values.push_back(min_quality);
  }

  output_histogram(values, "histogram_cell_quality_external_simplex_metric.cvs");
}

int main(int, char**)
{
//  std::freopen("log.txt", "w", stdout); //all output is written in "log.txt"

  std::vector<Point_3> points;
  std::vector<int> facets, cells;
  std::vector<Metric> metrics;

//  Constrain_surface* pdomain =
//    new Constrain_surface("bambimboum.mesh");
  Constrain_surface* pdomain =
    new Constrain_surface("/home/mrouxell/Data/OFF/fertility.off");

  //----------- pick a metric field! ----
//  Euclidean_metric_field* metric_field = new Euclidean_metric_field();
//  Hyperbolic_shock_metric_field* metric_field = new Hyperbolic_shock_metric_field(0.6);
  Custom_metric_field* metric_field = new Custom_metric_field();
//  External_metric_field* metric_field = new External_metric_field(*pdomain, "../../data/Anisotropy_CMP/3DSurface/Fandisk_Metric.txt");
//  External_metric_field* metric_field = new External_metric_field(*pdomain, "metric_at_input.txt");
//  Polyhedral_curvature_metric_field* metric_field = new Polyhedral_curvature_metric_field(*pdomain);
//  Implicit_curvature_metric_field* metric_field = new Implicit_curvature_metric_field(*pdomain);

//  const char* mesh_filename = "c2t3_fertility_tr_primal.mesh";
  const char* mesh_filename =
    "/home/mrouxell/anisomeshes/Thesis_Mael/results/ref_1350_primal.mesh";
//  const char* mesh_filename =
//    "/home/mrouxell/anisomeshes/research_report/Results/fertility.off";

  fetch_mesh(mesh_filename, points, facets, cells);
  compute_metrics(points, metric_field, metrics);

  facet_distortion_histogram(points, facets, metrics);
  facet_edge_length_histogram(points, facets, metrics);
  facet_angle_histogram(points, facets, metrics);
  facet_quality_histogram(points, facets, metrics);
  cell_distortion_histogram(points, cells, metrics);
  cell_edge_length_histogram(points, cells, metrics);
  cell_edge_length_histogram_midpoint_metric(points, cells, metric_field);
  cell_edge_length_histogram_simplex_metric(points, cells, metric_field);
  cell_dihedral_angle_histogram(points, cells, metric_field);
  cell_quality_histogram(points, cells, metric_field);

//  delete pdomain;
//  delete metric_field;

  return 0;
}


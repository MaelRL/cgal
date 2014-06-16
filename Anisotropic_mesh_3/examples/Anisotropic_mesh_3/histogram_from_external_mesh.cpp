#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Constrain_surface_3_polyhedral.h>
#include <CGAL/Metric.h>
#include <CGAL/Starset.h>

#include <CGAL/Euclidean_metric_field.h>
#include <CGAL/Implicit_curvature_metric_field.h>
#include <CGAL/Polyhedral_curvature_metric_field.h>
#include <Metric_field/Custom_metric_field.h>
#include <Metric_field/Hyperbolic_shock_metric_field.h>
#include <Metric_field/External_metric_field.h>

#include <ostream>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
//typedef typename CGAL::Exact_predicates_exact_constructions_kernel KExact;
typedef K                                                    KExact;

typedef typename K::FT                                       FT;
typedef typename K::Point_3                                  Point_3;

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
        facets.push_back(n1); facets.push_back(n2); facets.push_back(n3);
      }
    }

    if(word == "Tetrahedra")
    {
      in >> nc;
      for(int i=0; i<nc; ++i)
      {
        int n1, n2, n3, n4;
        in >> n1 >> n2 >> n3 >> n4 >> useless;
        cells.push_back(n1); cells.push_back(n2); cells.push_back(n3); cells.push_back(n4);
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
  std::cout << "lel: " << input_filename << std::endl;

  std::string extension = input_filename.substr(input_filename.find_last_of('.'));

  if (extension == ".off" || extension == ".OFF")
  {
    std::ifstream stream(filename);
    if(!stream)
      std::cerr << "Error: cannot read file " << std::endl;
    fetch_from_off(stream, points, facets);
  }
  else if (extension == ".mesh" || extension == ".MESH")
  {
    std::ifstream stream(filename);
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
                      const char* filename = "histogram.cvs")
{
  std::cout << "output: " << filename << std::endl;
  std::ofstream out(filename);
  std::size_t histo_n = histogram.size();
  std::cout << "histon: " << histo_n << std::endl;
  for(std::size_t i=0; i<histo_n; ++i)
  {
    FT val = min + (max-min)*((FT) i)/((FT) histo_n);
    std::cout << "valhisto: " << val << std::endl;
    out << i << "," << val << "," << histogram[i] << std::endl;
  }
}

template<typename Metric_field>
void facet_edge_length_histogram(std::vector<Point_3>& points,
                                 std::vector<int>& facets,
                                 const Metric_field* mf)
{
  if(facets.empty())
    return;
  std::cout << "facet edge external histo with size: " << facets.size() << std::endl;

  int histogram_size = 1000;
  FT limit_val = limit_val;
  FT facet_circumradius = 1.0;
  FT upper_bound = 2.0 * facet_circumradius;
  FT step_size = upper_bound / (FT) histogram_size;
  std::vector<int> histogram(histogram_size, 0);
  FT val, min_value = 1e30, max_value = -1e30;

  Star_traits star_traits;
  typename Star_traits::Compute_squared_distance_3 csd = star_traits.compute_squared_distance_3_object();

  for(int i=0; i<facets.size()/3;)
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
      Metric m = mf->compute_metric(ps[j]); // not very pretty...
      std::vector<Point_3> tps(3);

      for(int k=0; k<3; ++k)
        tps[k] = m.transform(ps[k]);


      val = CGAL::sqrt( csd(tps[0], tps[1]) );
      std::cout << "val: " << val << " " << val/step_size << std::endl;
      min_value = (std::min)(min_value, val);
      max_value = (std::max)(max_value, val);
      histogram[(std::min)(limit_val, std::floor(val/step_size))]++;

      val = CGAL::sqrt( csd(tps[0], tps[2]) );
      min_value = (std::min)(min_value, val);
      max_value = (std::max)(max_value, val);
      histogram[(std::min)(limit_val, std::floor(val/step_size))]++;

      val = CGAL::sqrt( csd(tps[1], tps[2]) );
      min_value = (std::min)(min_value, val);
      max_value = (std::max)(max_value, val);
      histogram[(std::min)(limit_val, std::floor(val/step_size))]++;
    }
  }

  output_histogram(histogram, 0., facet_circumradius, "histogram_facet_edge_length_external.cvs");
}

template<typename Metric_field>
void cell_edge_length_histogram(std::vector<Point_3>& points,
                                std::vector<int>& cells,
                                const Metric_field* mf)
{
  if(cells.empty())
    return;
  std::cout << "cell edge external histo with size: " << cells.size() << std::endl;

  int histogram_size = 1000;
  FT limit_val = histogram_size - 1.;
  FT cell_circumradius = 1.0;
  FT upper_bound = 2.0 * cell_circumradius;
  FT step_size = upper_bound / (FT) histogram_size;
  std::vector<int> histogram(histogram_size, 0);
  FT val, min_value = 1e30, max_value = -1e30;

  Star_traits star_traits;
  typename Star_traits::Compute_squared_distance_3 csd = star_traits.compute_squared_distance_3_object();

  for(int i=0; i<cells.size()/4;)
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
      Metric m = mf->compute_metric(ps[j]); // not very pretty...
      std::vector<Point_3> tps(4);

      for(int k=0; k<4; ++k)
        tps[k] = m.transform(ps[k]);


      val = CGAL::sqrt( csd(tps[0], tps[1]) );
      std::cout << "val: " << val << " " << val/step_size << std::endl;
      min_value = (std::min)(min_value, val);
      max_value = (std::max)(max_value, val);
      histogram[(std::min)(limit_val, std::floor(val/step_size))]++;

      val = CGAL::sqrt( csd(tps[0], tps[2]) );
      min_value = (std::min)(min_value, val);
      max_value = (std::max)(max_value, val);
      histogram[(std::min)(limit_val, std::floor(val/step_size))]++;

      val = CGAL::sqrt( csd(tps[0], tps[3]) );
      min_value = (std::min)(min_value, val);
      max_value = (std::max)(max_value, val);
      histogram[(std::min)(limit_val, std::floor(val/step_size))]++;

      val = CGAL::sqrt( csd(tps[1], tps[2]) );
      min_value = (std::min)(min_value, val);
      max_value = (std::max)(max_value, val);
      histogram[(std::min)(limit_val, std::floor(val/step_size))]++;

      val = CGAL::sqrt( csd(tps[1], tps[3]) );
      min_value = (std::min)(min_value, val);
      max_value = (std::max)(max_value, val);
      histogram[(std::min)(limit_val, std::floor(val/step_size))]++;

      val = CGAL::sqrt( csd(tps[2], tps[3]) );
      min_value = (std::min)(min_value, val);
      max_value = (std::max)(max_value, val);
      histogram[(std::min)(limit_val, std::floor(val/step_size))]++;
    }
  }
  std::cout << "done, outputing:" << std::endl;

  output_histogram(histogram, 0., cell_circumradius, "histogram_cell_edge_length_external.cvs");
}

int main(int argc, char** argv)
{
  std::vector<Point_3> points;
  std::vector<int> facets, cells;

  Constrain_surface* pdomain = new Constrain_surface("bambimboum.mesh");
//  Constrain_surface* pdomain = new Constrain_surface("../../data/Anisotropy_CMP/3DSurface/Fandisk.off");

  //----------- pick a metric field! ----
  Euclidean_metric_field* metric_field = new Euclidean_metric_field();
//  Hyperbolic_shock_metric_field* metric_field = new Hyperbolic_shock_metric_field(0.6);
//  Custom_metric_field* metric_field = new Custom_metric_field();

//  External_metric_field* metric_field = new External_metric_field(*pdomain, "../../data/Anisotropy_CMP/3DSurface/Fandisk_Metric.txt");
//  External_metric_field* metric_field = new External_metric_field(*pdomain, "metric_at_input.txt");
//  Polyhedral_curvature_metric_field* metric_field = new Polyhedral_curvature_metric_field(pdomain);
//  Implicit_curvature_metric_field* metric_field = new Implicit_curvature_metric_field(pdomain);

//  const char* mesh_filename = (argc>1)?argv[1]:"/home/mrouxel/cgal/Anisotropic_mesh_3/data/Anisotropy_CMP/3DSurface/Fandisk.off";
  const char* mesh_filename = (argc>1)?argv[1]:"/home/mrouxel/cgal/Anisotropic_mesh_3/examples/Anisotropic_mesh_3/bambimboum.mesh";

  fetch_mesh(mesh_filename, points, facets, cells);

  facet_edge_length_histogram(points, facets, metric_field);
  cell_edge_length_histogram(points, cells, metric_field);

  delete pdomain;
  delete metric_field;

  return 0;
}


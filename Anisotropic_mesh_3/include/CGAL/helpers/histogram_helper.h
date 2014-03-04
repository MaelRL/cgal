#ifndef CGAL_ANISOTROPIC_MESH_3_HISTOGRAM_HELPER_H
#define CGAL_ANISOTROPIC_MESH_3_HISTOGRAM_HELPER_H

#include <CGAL/helpers/combinatorics_helper.h>

#include <cmath>
#include <complex>
#include <cstdlib>
#include <ostream>
#include <set>
#include <vector>

namespace CGAL
{
namespace Anisotropic_mesh_3
{

enum Facet_histogram_type{APPROXIMATION,
                          FACET_DISTORTION,
                          FACET_QUALITY,
                          FACET_SIZE,
                          FACET_RATIO,
                          FACET_ANGLE};

enum Cell_histogram_type{CELL_DISTORTION,
                         CELL_QUALITY,
                         CELL_SIZE,
                         CELL_RATIO,
                         CELL_ANGLE};

template<typename FT>
void output_histogram(const std::vector<int>& histogram,
                      FT min, FT max,
                      const char* filename = "histogram.cvs")
{
  std::ofstream out(filename);
  std::size_t histo_n = histogram.size();
  for(std::size_t i=0; i<histo_n; ++i)
  {
    FT val = min + (max-min)*((FT) i)/((FT) histo_n);
    out << i << "," << val << "," << histogram[i] << std::endl;
  }
}

template<typename Star, typename Constrain_surface, typename Criteria>
void facet_histogram(const typename std::vector<Star*>& stars,
                     const Constrain_surface* const m_pConstrain,
                     const Criteria* const m_criteria,
                     const Facet_histogram_type hist_type = APPROXIMATION,
                     const bool verbose = false)
{
  typedef typename Star::FT FT;

  int divisions = 100;
  std::vector<int> histogram(divisions + 1, 0);
  int count = 0;
  FT max = -1e10, min = 1e10, sum = 0.;
  std::set<Facet_ijk> done;

  FT max_val = 0.;
  if(hist_type == APPROXIMATION)
  {
    std::cout << "Approximation Histo" << std::endl;
    max_val = m_criteria->approximation;
  }
  else if(hist_type == FACET_DISTORTION)
  {
    std::cout << "Facet Distortion Histo" << std::endl;
    max_val = m_criteria->distortion - 1.;
  }
  else if(hist_type == FACET_QUALITY)
  {
    std::cout << "Facet Quality Histo" << std::endl;
    max_val = 1.;
  }
  else if(hist_type == FACET_SIZE)
  {
    std::cout << "Facet Size Histo" << std::endl;
    max_val = m_criteria->facet_circumradius;
  }
  else if(hist_type == FACET_RATIO)
  {
    std::cout << "Facet Ratio Histo" << std::endl;
    max_val = m_criteria->facet_radius_edge_ratio;
  }
  else if(hist_type == FACET_ANGLE)
  {
    std::cout << "Facet Angle Histo" << std::endl;
    max_val = 0.; //todo
  }

  std::size_t N = stars.size();
  for(std::size_t i=0; i<N; ++i)
  {
    Star si = stars[i];

    if(!si->is_surface_star())
      continue;

    typename Star::Facet_set_iterator fit = si->begin_restricted_facets();
    typename Star::Facet_set_iterator fend = si->end_restricted_facets();
    for(; fit != fend; ++fit)
    {
      typename Star::Facet f = *fit;

      std::pair<typename std::set<Facet_ijk>::iterator, bool> is_insert_successful;
      is_insert_successful = done.insert(Facet_ijk(f));
      if(!is_insert_successful.second)
        continue;

      FT val = 0.;

      if(hist_type == APPROXIMATION)
      {
        typename Star::Point_3 cc;
        si->compute_dual_intersection(f, cc);
        FT sqd = m_pConstrain->compute_sq_approximation(barycenter(f));
        val = CGAL::sqrt(sqd);
      }
      else if(hist_type == FACET_DISTORTION)
      {
        int index = f.second;
        typename Star::Cell_handle c = f.first;
        for (int i = 0; i < 3; i++)
        {
          int i1 = (index + i + 1) % 4;
          int i2 = (index + (i + 1) % 3 + 1) % 4;
          val = (std::max)(val,
                           stars[c->vertex(i1)->info()]->metric().compute_distortion(
                                              stars[c->vertex(i2)->info()]->metric()));
        }
      }
      else if(hist_type == FACET_QUALITY)
      {
        const typename Star::Point_3& pa = si->metric().inverse_transform(f.first->vertex((f.second+1)%4)->point());
        const typename Star::Point_3& pb = si->metric().inverse_transform(f.first->vertex((f.second+2)%4)->point());
        const typename Star::Point_3& pc = si->metric().inverse_transform(f.first->vertex((f.second+3)%4)->point());

        Star* star_a = stars[f.first->vertex((f.second+1)%4)->info()];
        Star* star_b = stars[f.first->vertex((f.second+2)%4)->info()];
        Star* star_c = stars[f.first->vertex((f.second+3)%4)->info()];


        const typename Star::TPoint_3& tpa_a = star_a->metric().transform(pa);
        const typename Star::TPoint_3& tpb_a = star_a->metric().transform(pb);
        const typename Star::TPoint_3& tpc_a = star_a->metric().transform(pc);
        FT quality_in_a = star_a->compute_element_quality(tpa_a, tpb_a, tpc_a);

        const typename Star::TPoint_3& tpa_b = star_b->metric().transform(pa);
        const typename Star::TPoint_3& tpb_b = star_b->metric().transform(pb);
        const typename Star::TPoint_3& tpc_b = star_b->metric().transform(pc);
        FT quality_in_b = star_b->compute_element_quality(tpa_b, tpb_b, tpc_b);

        const typename Star::TPoint_3& tpa_c = star_c->metric().transform(pa);
        const typename Star::TPoint_3& tpb_c = star_c->metric().transform(pb);
        const typename Star::TPoint_3& tpc_c = star_c->metric().transform(pc);
        FT quality_in_c = star_c->compute_element_quality(tpa_c, tpb_c, tpc_c);

        val = (std::min)((std::min)(quality_in_a, quality_in_b),quality_in_c);
      }
      else if(hist_type == FACET_SIZE)
      {
        FT sqr = si->compute_squared_circumradius(f);
        val = CGAL::sqrt(sqr);
      }
      else if(hist_type == FACET_RATIO)
      {
        FT sqer = si->compute_squared_radius_edge_ratio(f);
        val = CGAL::sqrt(sqer);
      }
      //else if(hist_type == ANGLE) //todo

      if(val < min)
        min = val;
      if(val > max)
        max = val;

      count++;
      sum += val;

      if(hist_type == FACET_DISTORTION)
        val -= 1.;

      if(std::floor(((double) divisions)*val/max_val) >= divisions)
        histogram[divisions]++;
      else
        histogram[std::floor(((double) divisions)*val/max_val)]++;
    }
  }

  if(verbose)
  {
    for(unsigned int i = 0; i<divisions; i++)
    {
      double inf = ((double) i)/((double) divisions)*max_val;
      double sup = ((double) i+1)/((double) divisions)*max_val;

      if(hist_type == FACET_DISTORTION)
      {
        inf += 1.;
        sup += 1.;
      }

      std::cout << "\t" << inf << " " << sup << " : " << histogram[i];
      if(i%4)
        std::cout << "    |||    ";
      else
        std::cout << std::endl;
    }
    std::cout << std::endl << "min, max: " << min << " " << max << std::endl;
    std::cout << "average : " << sum / count << " (" << count << ")" << std::endl;
    std::cout << "\tabove : " << histogram[divisions] << std::endl;
  }

  if(hist_type == APPROXIMATION)
    output_histogram(histogram, 0., max_val, "histogram_approximation.cvs");
  else if(hist_type == FACET_DISTORTION)
    output_histogram(histogram, 1., max_val + 1., "histogram_facet_distortion.cvs");
  else if(hist_type == FACET_QUALITY)
    output_histogram(histogram, 0., 1., "histogram_facet_quality.cvs");
  else if(hist_type == FACET_SIZE)
    output_histogram(histogram, 0., max_val, "histogram_facet_size.cvs");
  else if(hist_type == FACET_RATIO)
    output_histogram(histogram, 0., max_val, "histogram_facet_ratio.cvs");
  else// if(hist_type == FACET_ANGLE)
    output_histogram(histogram, 0., max_val, "histogram_facet_angle.cvs");
}

template<typename Star, typename Criteria>
void cell_histogram(const typename std::vector<Star*>& stars,
                               const Criteria* const m_criteria,
                               const Cell_histogram_type hist_type = CELL_DISTORTION,
                               const bool verbose = false)
{
  typedef typename Star::FT FT;

  int divisions = 100;
  std::vector<int> histogram(divisions + 1, 0);
  int count = 0;
  FT max = -1e10, min = 1e10, sum = 0.;
  std::set<Cell_ijkl> done;

  FT max_val = 0.;
  if(hist_type == CELL_DISTORTION)
  {
    std::cout << "Cell Distortion Histo" << std::endl;
    max_val = m_criteria->distortion - 1.;
  }
  else if(hist_type == CELL_QUALITY)
  {
    std::cout << "Cell Quality Histo" << std::endl;
    max_val = 1.;
  }
  else if(hist_type == CELL_SIZE)
  {
    std::cout << "Cell Size Histo" << std::endl;
    max_val = m_criteria->facet_circumradius;
  }
  else if(hist_type == CELL_RATIO)
  {
    std::cout << "Cell Ratio Histo" << std::endl;
    max_val = m_criteria->facet_radius_edge_ratio;
  }
  else if(hist_type == CELL_ANGLE)
  {
    std::cout << "Cell Angle Histo" << std::endl;
    max_val = 0.; //todo
  }

  std::size_t N = stars.size();
  for(std::size_t i=0; i<N; ++i)
  {
    Star* si = stars[i];

    typename Star::Cell_handle_handle cit = si->begin_finite_star_cells();
    typename Star::Cell_handle_handle cend = si->end_finite_star_cells();
    for(; cit != cend; ++cit)
    {
      typename Star::Cell_handle c = *cit;
      if(!si->is_inside(c))
        continue;

      std::pair<typename std::set<Cell_ijkl>::iterator, bool> is_insert_successful;
      is_insert_successful = done.insert(Cell_ijkl(c));
      if(!is_insert_successful.second)
        continue;

      FT val = 0.;

      if(hist_type == CELL_DISTORTION)
      {
        for(int i = 0; i < 4; i++)
          for(int j = i + 1; j < 4; j++)
            val = (std::max)(val,
                             stars[c->vertex(i)->info()]->metric().compute_distortion(
                    stars[c->vertex(j)->info()]->metric()));
      }
      else if(hist_type == CELL_QUALITY)
      {

        const typename Star::Point_3& pa = si->metric().inverse_transform(c->vertex(0)->point());
        const typename Star::Point_3& pb = si->metric().inverse_transform(c->vertex(1)->point());
        const typename Star::Point_3& pc = si->metric().inverse_transform(c->vertex(2)->point());
        const typename Star::Point_3& pd = si->metric().inverse_transform(c->vertex(3)->point());

        Star* star_a = stars[c->vertex(0)->info()];
        Star* star_b = stars[c->vertex(1)->info()];
        Star* star_c = stars[c->vertex(2)->info()];
        Star* star_d = stars[c->vertex(3)->info()];

        const typename Star::TPoint_3& tpa_a = star_a->metric().transform(pa);
        const typename Star::TPoint_3& tpb_a = star_a->metric().transform(pb);
        const typename Star::TPoint_3& tpc_a = star_a->metric().transform(pc);
        const typename Star::TPoint_3& tpd_a = star_a->metric().transform(pd);
        FT quality_in_a = star_a->compute_element_quality(tpa_a, tpb_a, tpc_a, tpd_a);

        const typename Star::TPoint_3& tpa_b = star_b->metric().transform(pa);
        const typename Star::TPoint_3& tpb_b = star_b->metric().transform(pb);
        const typename Star::TPoint_3& tpc_b = star_b->metric().transform(pc);
        const typename Star::TPoint_3& tpd_b = star_b->metric().transform(pd);
        FT quality_in_b = star_b->compute_element_quality(tpa_b, tpb_b, tpc_b, tpd_b);

        const typename Star::TPoint_3& tpa_c = star_c->metric().transform(pa);
        const typename Star::TPoint_3& tpb_c = star_c->metric().transform(pb);
        const typename Star::TPoint_3& tpc_c = star_c->metric().transform(pc);
        const typename Star::TPoint_3& tpd_c = star_c->metric().transform(pd);
        FT quality_in_c = star_c->compute_element_quality(tpa_c, tpb_c, tpc_c, tpd_c);

        const typename Star::TPoint_3& tpa_d = star_d->metric().transform(pa);
        const typename Star::TPoint_3& tpb_d = star_d->metric().transform(pb);
        const typename Star::TPoint_3& tpc_d = star_d->metric().transform(pc);
        const typename Star::TPoint_3& tpd_d = star_d->metric().transform(pd);
        FT quality_in_d = star_c->compute_element_quality(tpa_d, tpb_d, tpc_d, tpd_d);

        val = (std::min)((std::min)((std::min)(quality_in_a, quality_in_b),quality_in_c),quality_in_d);
      }
      else if(hist_type == CELL_SIZE)
      {
        FT sqr = si->compute_squared_circumradius(c);
        val = CGAL::sqrt(sqr);
      }
      else if(hist_type == CELL_RATIO)
      {
        FT sqer = si->compute_squared_radius_edge_ratio(c);
        val = CGAL::sqrt(sqer);
      }

      if(val < min)
        min = val;
      if(val > max)
        max = val;

      count++;
      sum += val;

      if(hist_type == CELL_DISTORTION)
        val -= 1.;

      if(std::floor(((double) divisions)*val/max_val) >= divisions)
        histogram[divisions]++;
      else
        histogram[std::floor(((double) divisions)*val/max_val)]++;
    }
  }

  if(verbose)
  {
    for(unsigned int i = 0; i<divisions; i++)
    {
      double inf = ((double) i)/((double) divisions) * max_val;
      double sup = ((double) i+1)/((double) divisions) * max_val;

      if(hist_type == CELL_DISTORTION)
      {
        inf += 1.;
        sup += 1.;
      }

      std::cout << "\t" << inf << " " << sup << " : " << histogram[i];
      if(i%4)
        std::cout << "    |||    ";
      else
        std::cout << std::endl;
    }
    std::cout << std::endl << "min, max: " << min << " " << max << std::endl;
    std::cout << "average : " << sum / count << " (" << count << ")" << std::endl;
    std::cout << "above : " << histogram[divisions] << std::endl;
    std::cout << done.size() << " cells" << std::endl;
  }

  if(hist_type == CELL_DISTORTION)
    output_histogram(histogram, 1., max_val + 1., "histogram_cell_distortion.cvs");
  else if(hist_type == CELL_QUALITY)
    output_histogram(histogram, 0., 1., "histogram_cell_quality.cvs");
  else if(hist_type == CELL_SIZE)
    output_histogram(histogram, 0., max_val, "histogram_cell_size.cvs");
  else if(hist_type == CELL_RATIO)
    output_histogram(histogram, 0., max_val, "histogram_cell_ratio.cvs");
  else //if(hist_type == CELL_ANGLE)
    output_histogram(histogram, 0., max_val, "histogram _cell_angle.cvs");

}

template<typename Star, typename Constrain_surface, typename Criteria>
void all_facet_histograms(const typename std::vector<Star*>& stars,
                          const Constrain_surface* const m_pConstrain,
                          const Criteria* const m_criteria,
                          const bool verbose = false)
{
  facet_histogram(stars, m_pConstrain, m_criteria, APPROXIMATION, verbose);
  facet_histogram(stars, m_pConstrain, m_criteria, FACET_DISTORTION, verbose);
  facet_histogram(stars, m_pConstrain, m_criteria, FACET_QUALITY, verbose);
  facet_histogram(stars, m_pConstrain, m_criteria, FACET_SIZE, verbose);
  facet_histogram(stars, m_pConstrain, m_criteria, FACET_RATIO, verbose);
  //facet_histogram(stars, m_pConstrain, m_criteria, FACET_ANGLE, verbose);
}

template<typename Star, typename Criteria>
void all_cell_histograms(const typename std::vector<Star*>& stars,
                          const Criteria* const m_criteria,
                          const bool verbose = false)
{
  cell_histogram(stars, m_criteria, CELL_DISTORTION, verbose);
  cell_histogram(stars, m_criteria, CELL_QUALITY, verbose);
  cell_histogram(stars, m_criteria, CELL_SIZE, verbose);
  cell_histogram(stars, m_criteria, CELL_RATIO, verbose);
  //cell_histogram(stars, m_criteria, CELL_ANGLE, verbose);
}

} //namespace Aniso
} //namespace CGAL

#endif

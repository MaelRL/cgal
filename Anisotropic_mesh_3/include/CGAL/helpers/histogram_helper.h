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

enum Facet_histogram_type
{
  APPROXIMATION,
  FACET_DISTORTION,
  FACET_QUALITY,
  FACET_SIZE,
  FACET_RATIO,
  FACET_ANGLE
};

enum Cell_histogram_type
{
  CELL_DISTORTION,
  CELL_QUALITY,
  CELL_SIZE,
  CELL_RATIO,
  CELL_ANGLE
};

template <typename Star>
typename Star::FT
minimum_dihedral_angle(Star* star,
                       const typename Star::Point_3& p0,
                       const typename Star::Point_3& p1,
                       const typename Star::Point_3& p2,
                       const typename Star::Point_3& p3)
{
  typedef typename Star::Kernel             K;
  typedef typename K::FT                    FT;
  typedef typename Star::TPoint_3           TPoint;
  K k = K();

  const TPoint& tp0 = star->metric().transform(p0);
  const TPoint& tp1 = star->metric().transform(p1);
  const TPoint& tp2 = star->metric().transform(p2);
  const TPoint& tp3 = star->metric().transform(p3);

  typename K::Compute_determinant_3 determinant =
      k.compute_determinant_3_object();
  typename K::Construct_cross_product_vector_3 cp =
      k.construct_cross_product_vector_3_object();

  typename K::Compute_scalar_product_3 sp =
      k.compute_scalar_product_3_object();

  typename K::Vector_3 v01 = tp1-tp0;
  typename K::Vector_3 v02 = tp2-tp0;
  typename K::Vector_3 v03 = tp3-tp0;
  typename K::Vector_3 v12 = tp2-tp1;
  typename K::Vector_3 v13 = tp3-tp1;
  typename K::Vector_3 v23 = tp3-tp2;

  typename K::Vector_3 v_01_02 = cp(v01,v02);
  FT a_012 = v_01_02*v_01_02;

  typename K::Vector_3 v_01_03 = cp(v01,v03);
  FT a_013 = v_01_03*v_01_03;

  typename K::Vector_3 v_12_13 = cp(v12,v13);
  FT a_123 = v_12_13*v_12_13;

  typename K::Vector_3 v_02_03 = cp(v02,v03);
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

template <typename Star>
typename Star::FT
dihedral_angle(Star* star,
               const typename Star::Point_3& a,
               const typename Star::Point_3& b,
               const typename Star::Point_3& c,
               const typename Star::Point_3& d)
{
  typedef typename Star::Kernel           K;
  typedef typename Star::TPoint_3         TPoint;
  typedef typename K::Vector_3            Vector_3;
  typedef typename K::FT                  FT;
  K k = K();

  const TPoint& ta = star->metric().transform(a);
  const TPoint& tb = star->metric().transform(b);
  const TPoint& tc = star->metric().transform(c);
  const TPoint& td = star->metric().transform(d);

  typename K::Construct_vector_3 vector = k.construct_vector_3_object();
  typename K::Construct_cross_product_vector_3 cross_product =
      k.construct_cross_product_vector_3_object();
  typename K::Compute_squared_distance_3 sq_distance =
      k.compute_squared_distance_3_object();
  typename K::Compute_scalar_product_3 scalar_product =
      k.compute_scalar_product_3_object();

  const Vector_3 ab = vector(ta, tb);
  const Vector_3 ac = vector(ta, tc);
  const Vector_3 ad = vector(ta, td);

  const Vector_3 abad = cross_product(ab, ad);
  const double x = CGAL::to_double(scalar_product(cross_product(ab, ac), abad));
  const double l_ab = CGAL::sqrt(CGAL::to_double(sq_distance(a, b)));
  const double y = l_ab * CGAL::to_double(scalar_product(ac, abad));

  return FT(std::atan2(y, x) * 180 / CGAL_PI );
}

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

template<typename Starset, typename Constrain_surface, typename Criteria>
void facet_histogram(const Starset& stars,
                     const Constrain_surface* const m_pConstrain,
                     const Criteria* const m_criteria,
                     const Facet_histogram_type hist_type = APPROXIMATION,
                     const bool verbose = false)
{
  typedef typename Starset::FT FT;
  typedef typename Starset::Point_3 Point_3;
  typedef typename Starset::TPoint_3 TPoint_3;
  typedef typename Starset::Facet Facet;
  typedef typename Starset::Star_handle Star_handle;

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
    max_val = 3.0; //m_criteria->distortion - 1.;
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
    Star_handle si = stars[i];

    if(!si->is_surface_star())
      continue;

    typename Starset::Facet_set_iterator fit = si->restricted_facets_begin();
    typename Starset::Facet_set_iterator fend = si->restricted_facets_end();
    for(; fit != fend; ++fit)
    {
      Facet f = *fit;

      std::pair<typename std::set<Facet_ijk>::iterator, bool> is_insert_successful;
      is_insert_successful = done.insert(Facet_ijk(f));
      if(!is_insert_successful.second)
        continue;

      Star_handle star_a = stars[f.first->vertex((f.second+1)%4)->info()];
      Star_handle star_b = stars[f.first->vertex((f.second+2)%4)->info()];
      Star_handle star_c = stars[f.first->vertex((f.second+3)%4)->info()];

      const Point_3& pa = si->metric().inverse_transform(f.first->vertex((f.second+1)%4)->point());
      const Point_3& pb = si->metric().inverse_transform(f.first->vertex((f.second+2)%4)->point());
      const Point_3& pc = si->metric().inverse_transform(f.first->vertex((f.second+3)%4)->point());

      const TPoint_3& tpa_a = star_a->metric().transform(pa);
      const TPoint_3& tpb_a = star_a->metric().transform(pb);
      const TPoint_3& tpc_a = star_a->metric().transform(pc);

      const TPoint_3& tpa_b = star_b->metric().transform(pa);
      const TPoint_3& tpb_b = star_b->metric().transform(pb);
      const TPoint_3& tpc_b = star_b->metric().transform(pc);

      const TPoint_3& tpa_c = star_c->metric().transform(pa);
      const TPoint_3& tpb_c = star_c->metric().transform(pb);
      const TPoint_3& tpc_c = star_c->metric().transform(pc);

      FT val = 0.;

      if(hist_type == APPROXIMATION)
      {
        typename Starset::Point_3 cc, bf;
        si->compute_dual_intersection(f, cc);
        bf = CGAL::barycenter(stars[f.first->vertex((f.second+1)%4)->info()]->center_point(), 1./3.,
                              stars[f.first->vertex((f.second+2)%4)->info()]->center_point(), 1./3.,
                              stars[f.first->vertex((f.second+3)%4)->info()]->center_point(), 1./3.);
        FT sqd = m_pConstrain->compute_sq_approximation(bf);
        val = CGAL::sqrt(sqd);
      }
      else if(hist_type == FACET_DISTORTION)
      {
        int index = f.second;
        typename Starset::Cell_handle c = f.first;
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
        FT quality_in_a = star_a->compute_element_quality(tpa_a, tpb_a, tpc_a);
        FT quality_in_b = star_b->compute_element_quality(tpa_b, tpb_b, tpc_b);
        FT quality_in_c = star_c->compute_element_quality(tpa_c, tpb_c, tpc_c);

        val = (std::min)((std::min)(quality_in_a, quality_in_b),quality_in_c);
      }
      else if(hist_type == FACET_SIZE)
      {
        FT sqr_in_a = star_a->compute_squared_circumradius(tpa_a, tpb_a, tpc_a);
        FT sqr_in_b = star_b->compute_squared_circumradius(tpa_b, tpb_b, tpc_b);
        FT sqr_in_c = star_c->compute_squared_circumradius(tpa_c, tpb_c, tpc_c);

        val = CGAL::sqrt( (std::max)((std::max)(sqr_in_a, sqr_in_b), sqr_in_c) );
      }
      else if(hist_type == FACET_RATIO)
      {
        FT sqer_in_a = star_a->compute_squared_radius_edge_ratio(tpa_a, tpb_a, tpc_a);
        FT sqer_in_b = star_b->compute_squared_radius_edge_ratio(tpa_b, tpb_b, tpc_b);
        FT sqer_in_c = star_c->compute_squared_radius_edge_ratio(tpa_c, tpb_c, tpc_c);

        val = CGAL::sqrt( (std::max)((std::max)(sqer_in_a, sqer_in_b), sqer_in_c) );
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
    std::cout << "\tabove : " << histogram[divisions] << std::endl;
  }
  std::cout << std::endl << "min, max: " << min << " " << max << std::endl;
  std::cout << "average : " << sum / count << " (" << count << ")" << std::endl;
  std::cout << done.size() << " facets" << std::endl;

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

template<typename Starset, typename Criteria>
void cell_histogram(const Starset& stars,
                    const Criteria* const m_criteria,
                    const Cell_histogram_type hist_type = CELL_DISTORTION,
                    const bool verbose = false)
{
  typedef typename Starset::FT FT;
  typedef typename Starset::Point_3 Point_3;
  typedef typename Starset::TPoint_3 TPoint_3;
  typedef typename Starset::Star_handle Star_handle;
  typedef typename Starset::Cell_handle Cell_handle;

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
    max_val = 180.;
  }

  std::size_t N = stars.size();
  for(std::size_t i=0; i<N; ++i)
  {
    Star_handle si = stars[i];

    typename Starset::Cell_handle_handle cit = si->finite_star_cells_begin();
    typename Starset::Cell_handle_handle cend = si->finite_star_cells_end();
    for(; cit != cend; ++cit)
    {
      Cell_handle c = *cit;
      if(!si->is_inside(c))
        continue;

      std::pair<typename std::set<Cell_ijkl>::iterator, bool> is_insert_successful;
      is_insert_successful = done.insert(Cell_ijkl(c));
      if(!is_insert_successful.second)
        continue;

      const Point_3& pa = si->metric().inverse_transform(c->vertex(0)->point());
      const Point_3& pb = si->metric().inverse_transform(c->vertex(1)->point());
      const Point_3& pc = si->metric().inverse_transform(c->vertex(2)->point());
      const Point_3& pd = si->metric().inverse_transform(c->vertex(3)->point());

      Star_handle star_a = stars[c->vertex(0)->info()];
      Star_handle star_b = stars[c->vertex(1)->info()];
      Star_handle star_c = stars[c->vertex(2)->info()];
      Star_handle star_d = stars[c->vertex(3)->info()];

      const TPoint_3& tpa_a = star_a->metric().transform(pa);
      const TPoint_3& tpb_a = star_a->metric().transform(pb);
      const TPoint_3& tpc_a = star_a->metric().transform(pc);
      const TPoint_3& tpd_a = star_a->metric().transform(pd);

      const TPoint_3& tpa_b = star_b->metric().transform(pa);
      const TPoint_3& tpb_b = star_b->metric().transform(pb);
      const TPoint_3& tpc_b = star_b->metric().transform(pc);
      const TPoint_3& tpd_b = star_b->metric().transform(pd);

      const TPoint_3& tpa_c = star_c->metric().transform(pa);
      const TPoint_3& tpb_c = star_c->metric().transform(pb);
      const TPoint_3& tpc_c = star_c->metric().transform(pc);
      const TPoint_3& tpd_c = star_c->metric().transform(pd);

      const TPoint_3& tpa_d = star_d->metric().transform(pa);
      const TPoint_3& tpb_d = star_d->metric().transform(pb);
      const TPoint_3& tpc_d = star_d->metric().transform(pc);
      const TPoint_3& tpd_d = star_d->metric().transform(pd);

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
        FT quality_in_a = star_a->compute_element_quality(tpa_a, tpb_a, tpc_a, tpd_a);
        FT quality_in_b = star_b->compute_element_quality(tpa_b, tpb_b, tpc_b, tpd_b);
        FT quality_in_c = star_c->compute_element_quality(tpa_c, tpb_c, tpc_c, tpd_c);
        FT quality_in_d = star_c->compute_element_quality(tpa_d, tpb_d, tpc_d, tpd_d);

        val = (std::min)((std::min)((std::min)(quality_in_a, quality_in_b),quality_in_c),quality_in_d);
      }
      else if(hist_type == CELL_SIZE)
      {
        FT sqr_in_a = star_a->compute_squared_circumradius(tpa_a, tpb_a, tpc_a, tpd_a);
        FT sqr_in_b = star_b->compute_squared_circumradius(tpa_b, tpb_b, tpc_b, tpd_b);
        FT sqr_in_c = star_c->compute_squared_circumradius(tpa_c, tpb_c, tpc_c, tpd_c);
        FT sqr_in_d = star_c->compute_squared_circumradius(tpa_d, tpb_d, tpc_d, tpd_d);

        val = CGAL::sqrt( (std::max)((std::max)((std::min)(sqr_in_a, sqr_in_b), sqr_in_c), sqr_in_d) );
      }
      else if(hist_type == CELL_RATIO)
      {
        FT sqer_in_a = star_a->compute_squared_radius_edge_ratio(tpa_a, tpb_a, tpc_a, tpd_a);
        FT sqer_in_b = star_b->compute_squared_radius_edge_ratio(tpa_b, tpb_b, tpc_b, tpd_b);
        FT sqer_in_c = star_c->compute_squared_radius_edge_ratio(tpa_c, tpb_c, tpc_c, tpd_c);
        FT sqer_in_d = star_c->compute_squared_radius_edge_ratio(tpa_d, tpb_d, tpc_d, tpd_d);

        val = CGAL::sqrt( (std::max)((std::max)((std::max)(sqer_in_a, sqer_in_b), sqer_in_c), sqer_in_d) );
      }
      else if(hist_type == CELL_ANGLE)
      {
        FT min_angle_a = minimum_dihedral_angle(star_a, pa, pb, pc, pd);
        FT min_angle_b = minimum_dihedral_angle(star_b, pa, pb, pc, pd);
        FT min_angle_c = minimum_dihedral_angle(star_c, pa, pb, pc, pd);
        FT min_angle_d = minimum_dihedral_angle(star_d, pa, pb, pc, pd);

        val = (std::min)((std::min)((std::min)(min_angle_a, min_angle_b), min_angle_c), min_angle_d);
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
    std::cout << "above : " << histogram[divisions] << std::endl;
  }
  std::cout << "min, max: " << min << " " << max << std::endl;
  std::cout << "average : " << sum / count << " (" << count << ")" << std::endl;
  std::cout << done.size() << " cells" << std::endl;

  if(hist_type == CELL_DISTORTION)
    output_histogram(histogram, 1., max_val + 1., "histogram_cell_distortion.cvs");
  else if(hist_type == CELL_QUALITY)
    output_histogram(histogram, 0., 1., "histogram_cell_quality.cvs");
  else if(hist_type == CELL_SIZE)
    output_histogram(histogram, 0., max_val, "histogram_cell_size.cvs");
  else if(hist_type == CELL_RATIO)
    output_histogram(histogram, 0., max_val, "histogram_cell_ratio.cvs");
  else //if(hist_type == CELL_ANGLE)
    output_histogram(histogram, 0., max_val, "histogram_cell_angle.cvs");

}

template<typename Starset>
void dihedral_angle_histogram(const Starset& stars,
                              const bool verbose = false)
{
  std::set<Cell_ijkl> done;
  std::vector<int> histogram(181, 0);
  typename Starset::FT val, min_value = 1e30, max_value = -1e30;

  std::size_t N = stars.size();
  for(std::size_t i=0; i<N; ++i)
  {
    typename Starset::Star_handle si = stars[i];

    typename Starset::Cell_handle_handle cit = si->finite_star_cells_begin();
    typename Starset::Cell_handle_handle cend = si->finite_star_cells_end();
    for(; cit != cend; ++cit)
    {
      typename Starset::Cell_handle c = *cit;
      if(!si->is_inside(c))
        continue;

      std::pair<typename std::set<Cell_ijkl>::iterator, bool> is_insert_successful;
      is_insert_successful = done.insert(Cell_ijkl(c));
      if(!is_insert_successful.second)
        continue;

      const typename Starset::Point_3& pa = si->metric().inverse_transform(c->vertex(0)->point());
      const typename Starset::Point_3& pb = si->metric().inverse_transform(c->vertex(1)->point());
      const typename Starset::Point_3& pc = si->metric().inverse_transform(c->vertex(2)->point());
      const typename Starset::Point_3& pd = si->metric().inverse_transform(c->vertex(3)->point());

      for(int i=0; i<4; ++i)
      {
        typename Starset::Star_handle star = stars[c->vertex(i)->info()];
        val = dihedral_angle(star, pa, pb, pc, pd);
        histogram[std::abs(std::floor(val))]++;
        min_value = (std::min)(min_value, val);
        max_value = (std::max)(max_value, val);

        val = dihedral_angle(star, pa, pc, pb, pd);
        histogram[std::abs(std::floor(val))]++;
        min_value = (std::min)(min_value, val);
        max_value = (std::max)(max_value, val);

        val = dihedral_angle(star, pa, pd, pb, pc);
        histogram[std::abs(std::floor(val))]++;
        min_value = (std::min)(min_value, val);
        max_value = (std::max)(max_value, val);

        val = dihedral_angle(star, pb, pc, pa, pd);
        histogram[std::abs(std::floor(val))]++;
        min_value = (std::min)(min_value, val);
        max_value = (std::max)(max_value, val);

        val = dihedral_angle(star, pb, pd, pa, pc);
        histogram[std::abs(std::floor(val))]++;
        min_value = (std::min)(min_value, val);
        max_value = (std::max)(max_value, val);

        val = dihedral_angle(star, pc, pd, pa, pb);
        histogram[std::abs(std::floor(val))]++;
        min_value = (std::min)(min_value, val);
        max_value = (std::max)(max_value, val);
      }
    }
  }

  if(verbose)
  {
    std::cout << "full dihedral angle histogram" << std::endl;
    std::cout << "min, max: " << min_value << " " << max_value << std::endl;
    std::cout << done.size() << " cells" << std::endl;
  }
  output_histogram(histogram, 0., 180., "histogram_dihedral_angle.cvs");
}

template<typename Starset>
typename Starset::FT edge_length(const typename Starset::Star_handle s,
                                 const typename Starset::Point_3& p1,
                                 const typename Starset::Point_3& p2)
{
  typename Starset::Traits::Compute_squared_distance_3 csd =
                            s->traits()->compute_squared_distance_3_object();

  typename Starset::TPoint_3 tp1 = s->metric().transform(p1);
  typename Starset::TPoint_3 tp2 = s->metric().transform(p2);

  typename Starset::FT d_tp1_tp2 = csd(tp1, tp2);

  return CGAL::sqrt(d_tp1_tp2);
}

template<typename Starset, typename Criteria>
void facet_edge_length_histogram(const Starset& stars,
                                 const Criteria* const criteria,
                                 const bool verbose = false)
{
  typedef typename Starset::FT FT;

  int histogram_size = 1000;
  FT limit_val = histogram_size - 1.;
  FT upper_bound = 2.0 * criteria->facet_circumradius; // will obv fail if r_0 = 0
  FT step_size = upper_bound / (FT) histogram_size;
  std::vector<int> histogram(histogram_size, 0);
  FT val, min_value = 1e30, max_value = -1e30;

  std::size_t N = stars.size();
  for(std::size_t i=0; i<N; ++i)
  {
    typename Starset::Star_handle si = stars[i];
    if(!si->is_surface_star())
      continue;

    typename Starset::Facet_set_iterator fit = si->restricted_facets_begin();
    typename Starset::Facet_set_iterator fend = si->restricted_facets_end();
    for(; fit != fend; ++fit)
    {
      for(int i=1; i<4; ++i)
      {
        typename Starset::Star_handle star = stars[fit->first->vertex((fit->second+i)%4)->info()];

        typename Starset::Traits::Compute_squared_distance_3 csd =
                                  star->traits()->compute_squared_distance_3_object();

        typename Starset::TPoint_3 tp1 = star->metric().transform(fit->first->vertex((fit->second+1)%4)->point());
        typename Starset::TPoint_3 tp2 = star->metric().transform(fit->first->vertex((fit->second+2)%4)->point());
        typename Starset::TPoint_3 tp3 = star->metric().transform(fit->first->vertex((fit->second+3)%4)->point());

        val = CGAL::sqrt( csd(tp1, tp2) );
        min_value = (std::min)(min_value, val);
        max_value = (std::max)(max_value, val);
        histogram[(std::min)(limit_val, std::floor(val/step_size))]++;

        val = CGAL::sqrt( csd(tp2, tp3) );
        min_value = (std::min)(min_value, val);
        max_value = (std::max)(max_value, val);
        histogram[(std::min)(limit_val, std::floor(val/step_size))]++;

        val = CGAL::sqrt( csd(tp1, tp3) );
        min_value = (std::min)(min_value, val);
        max_value = (std::max)(max_value, val);
        histogram[(std::min)(limit_val, std::floor(val/step_size))]++;
      }
    }
  }

  std::cout << "facet edge length histogram" << std::endl;
  std::cout << "min, max: " << min_value << " " << max_value << std::endl;

  output_histogram(histogram, 0., 2.0 * criteria->facet_circumradius, "histogram_facet_edge_length.cvs");
}

template<typename Starset, typename Criteria>
void cell_edge_length_histogram(const Starset& stars,
                                const Criteria* const criteria,
                                const bool verbose = false)
{
  typedef typename Starset::FT FT;

  int histogram_size = 1000;
  FT limit_val = histogram_size - 1.;
  FT upper_bound = 2.0 * criteria->cell_circumradius; // will obv fail if r_0 = 0
  FT step_size = upper_bound / (FT) histogram_size;
  std::vector<int> histogram(histogram_size, 0);
  FT val, min_value = 1e30, max_value = -1e30;

  std::size_t N = stars.size();
  std::cout << "N: " << N << std::endl;
  for(std::size_t i=0; i<N; ++i)
  {
    typename Starset::Star_handle si = stars[i];

    typename Starset::Cell_handle_handle chit = si->finite_star_cells_begin();
    typename Starset::Cell_handle_handle chend = si->finite_star_cells_end();
    for(; chit!=chend; ++chit)
    {
      typename Starset::Cell_handle cit = *chit;

      if(!si->is_inside(cit))
        continue;

      for(int i=0; i<4; ++i)
      {
        typename Starset::Star_handle star = stars[cit->vertex(i)->info()];

        typename Starset::Traits::Compute_squared_distance_3 csd =
                                  star->traits()->compute_squared_distance_3_object();

        typename Starset::TPoint_3 tp0 = star->metric().transform(cit->vertex(0)->point());
        typename Starset::TPoint_3 tp1 = star->metric().transform(cit->vertex(1)->point());
        typename Starset::TPoint_3 tp2 = star->metric().transform(cit->vertex(2)->point());
        typename Starset::TPoint_3 tp3 = star->metric().transform(cit->vertex(3)->point());

        val = CGAL::sqrt( csd(tp0, tp1) );
        min_value = (std::min)(min_value, val);
        max_value = (std::max)(max_value, val);
        histogram[(std::min)(limit_val, std::floor(val/step_size))]++;

        val = CGAL::sqrt( csd(tp0, tp2) );
        min_value = (std::min)(min_value, val);
        max_value = (std::max)(max_value, val);
        histogram[(std::min)(limit_val, std::floor(val/step_size))]++;

        val = CGAL::sqrt( csd(tp0, tp3) );
        min_value = (std::min)(min_value, val);
        max_value = (std::max)(max_value, val);
        histogram[(std::min)(limit_val, std::floor(val/step_size))]++;

        val = CGAL::sqrt( csd(tp1, tp2) );
        min_value = (std::min)(min_value, val);
        max_value = (std::max)(max_value, val);
        histogram[(std::min)(limit_val, std::floor(val/step_size))]++;

        val = CGAL::sqrt( csd(tp1, tp3) );
        min_value = (std::min)(min_value, val);
        max_value = (std::max)(max_value, val);
        histogram[(std::min)(limit_val, std::floor(val/step_size))]++;

        val = CGAL::sqrt( csd(tp2, tp3) );
        min_value = (std::min)(min_value, val);
        max_value = (std::max)(max_value, val);
        histogram[(std::min)(limit_val, std::floor(val/step_size))]++;
      }
    }
  }

  std::cout << "cell edge length histogram" << std::endl;
  std::cout << "min, max: " << min_value << " " << max_value << std::endl;

  output_histogram(histogram, 0., 2.0 * criteria->cell_circumradius, "histogram_cell_edge_length.cvs");
}


template<typename Starset, typename Constrain_surface, typename Criteria>
void all_facet_histograms(const Starset& stars,
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

  facet_edge_length_histogram(stars, m_criteria, verbose);
}

template<typename Starset, typename Criteria>
void all_cell_histograms(const Starset& stars,
                         const Criteria* const m_criteria,
                         const bool verbose = false)
{
  cell_histogram(stars, m_criteria, CELL_DISTORTION, verbose);
  cell_histogram(stars, m_criteria, CELL_QUALITY, verbose);
  cell_histogram(stars, m_criteria, CELL_SIZE, verbose);
  cell_histogram(stars, m_criteria, CELL_RATIO, verbose);
  cell_histogram(stars, m_criteria, CELL_ANGLE, verbose);
  dihedral_angle_histogram(stars, verbose);

  cell_edge_length_histogram(stars, m_criteria, verbose);
}

} //namespace Aniso
} //namespace CGAL

#endif

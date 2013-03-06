#include <vector>


template<typename C3T3, typename OutputIterator>
void compute_triangulation_poles(const C3T3& c3t3,
                                 OutputIterator oit)
{   
  typedef typename C3T3::Triangulation Triangulation;
  typedef typename Triangulation::Geom_traits::Point_3 Point_3;
  typedef typename Triangulation::Geom_traits::Vector_3 Vector_3;
  typedef typename Triangulation::Cell_handle Cell_handle;

  const Triangulation& tr = c3t3.triangulation();
  typename Triangulation::Finite_vertices_iterator v;
  for(v = tr.finite_vertices_begin();
      v != tr.finite_vertices_end();
      ++v)
  {
    Point_3 pi = v->point();

    if(c3t3.in_dimension(v) == 3)
    {
      std::cerr << "Warning : in compute_triangulation_poles ("<< pi <<") is inside" << std::endl;
      continue;
    }
                    
    std::vector<Point_3> candidates;
    Point_3 pole;
    // loop (1/2) : collect circumcenters and one-side poles
    std::vector<Cell_handle> cells;
    tr.incident_cells(v, std::back_inserter(cells));
    double max_sqd = 0.;
    for(unsigned int i = 0; i < cells.size(); ++i)
    {
      Cell_handle c = cells[i];
      if(tr.is_infinite(c))
        continue;
      Point_3 cc = c->circumcenter();
      candidates.push_back(cc);
      double sqd = CGAL::squared_distance(pi, cc);
      if(sqd > max_sqd)
      {
        pole = cc;
        max_sqd = sqd;
      }
    }
    if(max_sqd == 0.) //no pole could be found
      continue;
    else
      *oit++ = pole;

    // loop (2/2) : collect, among cc's, the other side poles 
    Vector_3 vec(pi, pole);
    double min_scal_prod = 0.; // will be < 0
    for(std::size_t j = 0; j < candidates.size(); ++j)
    {
      Point_3 cc = candidates[j];
      Vector_3 vj(pi, cc);
      double scal_prod = vec * vj;
      if(scal_prod < min_scal_prod)
      {
        pole = cc;
        min_scal_prod = scal_prod;
      }
    }
    if(min_scal_prod < 0.)
      *oit++ = pole;
  }
}

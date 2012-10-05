#ifndef STATISTICS_HELPER_H
#define STATISTICS_HELPER_H

#include <map>
#include <fstream>

template<typename StarSet>
void histogram_vertices_per_star(StarSet& stars)
{
  typedef typename StarSet::Star_handle Star_handle;
  typedef typename StarSet::Star_iterator Star_iterator;

  std::map<unsigned int/*nb of vertices*/, unsigned int/*nb of stars for this number*/> countmap;
  unsigned int maxnb = 0;
  Star_iterator si = stars.m_stars.begin();
  Star_iterator siend = stars.m_stars.end();
  for (; si != siend; si++) 
  {
    Star_handle star = *si;
    unsigned int nb = star->number_of_vertices(); 
    maxnb = (std::max)(maxnb, nb);

    if(countmap.find(nb) == countmap.end())
      countmap.insert(std::pair<int,int>(nb, 1));
    else
      countmap[nb]++;
  }

  typename std::ofstream fx("vertices_per_star.cgal.txt");
  typename std::map<unsigned int, unsigned int>::iterator it = countmap.begin();
  unsigned int i;
  for(i = 0; i < maxnb; i++)
  {
//    CGAL_HISTOGRAM_PROFILER("[nbv / star]", it->first);
    if(it->first == i)
    {
      fx << i << " " << it->second << "\n";
      it++; //o.w. it stays at this rank and i grows
    }
    else
      fx << i << " " << 0 << "\n";
  }
  fx.close();
}


#endif //STATISTICS_HELPER_H
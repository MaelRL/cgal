#ifndef CGAL_KD_TREE_STAR_SET_H
#define CGAL_KD_TREE_STAR_SET_H

#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>

#include <CGAL/Fuzzy_iso_box.h>

#include <boost/iterator/counting_iterator.hpp>
#include <vector>

namespace CGAL
{
//definition of the property map and get
//function as friend function to have access to
//private member
template<typename K, typename StarHandle>
class Star_property_map
{
public:
  const std::vector<StarHandle>& stars;
public:
  typedef typename K::Point_3   value_type;
  typedef const value_type&     reference;
  //typedef std::size_t           Point;
  typedef int                   Point;
  typedef Point                 key_type;
  typedef boost::readable_property_map_tag category;  

  Star_property_map(const std::vector<StarHandle>& s) : stars(s){}

  friend reference get(const Star_property_map& ppmap, key_type i) 
  { return ppmap.stars[i]->center_point(); }

}; // end class Star_property_map


template<typename K, typename StarHandle>
class Kd_tree_for_star_set : public CGAL::Kd_tree<
    CGAL::Search_traits_adapter< typename Star_property_map<K, StarHandle>::Point,
                                 Star_property_map<K, StarHandle>,
                                 CGAL::Search_traits_3<K> > >
{

private:
  typedef typename CGAL::Kd_tree<
    CGAL::Search_traits_adapter< typename Star_property_map<K, StarHandle>::Point,
                                 Star_property_map<K, StarHandle>,
                                 CGAL::Search_traits_3<K> > >                Base;
public:
  typedef Star_property_map<K, StarHandle> Star_pmap;
  typedef typename Star_pmap::Point                 Point;
  typedef Point                            Point_d; // as in Kd_tree
  typedef CGAL::Search_traits_3<K>         Traits_base;
  typedef CGAL::Search_traits_adapter<Point, Star_pmap, Traits_base> Traits;
  typedef CGAL::Kd_tree<Traits>                                      Self;
  typedef typename Base::Splitter Splitter;

public:
  typedef CGAL::Fuzzy_iso_box<Traits> Box_query;
  typedef typename Star_pmap::key_type key_type;
  
private:
  std::size_t m_insertion_buffer_max_size;
  std::vector<Point_d> m_insertion_buffer;

public:
  Kd_tree_for_star_set(const std::vector<StarHandle>& stars, 
                       const std::size_t bs = 10) :
    Base(boost::counting_iterator<std::size_t>(0),
         boost::counting_iterator<std::size_t>(stars.size()),
         typename Base::Splitter(),
         Traits(Star_pmap(stars))), 
    m_insertion_buffer_max_size(bs),
    m_insertion_buffer() {}


public:
  bool empty()
  {
    return Base::empty() && m_insertion_buffer.empty();
  }
  void clear()
  {
    Base::clear();
    m_insertion_buffer.clear();
  }
  
  void insert(const Point_d& p)
  {
    if(m_insertion_buffer.size() >= m_insertion_buffer_max_size)
      insert_buffer();
    m_insertion_buffer.push_back(p);
    
    //if(m_insertion_buffer.size() < m_insertion_buffer_max_size)
    //  m_insertion_buffer.push_back(p);
    //else
    //{
    //  insert_buffer();
    //  Base::insert(p);
    //}
  }

  void remove_last() 
  {
    m_insertion_buffer.pop_back();
    // note : the last inserted point always is in the insertion buffer
    // (because of the way insert(p) is written)
  }

  void insert_buffer()
  {
    std::size_t N = m_insertion_buffer.size();
    std::size_t i;
    for(i = 0; i < N; i++)
      Base::insert(m_insertion_buffer[i]/*Point_d*/);
    m_insertion_buffer.clear();
  }

  template<typename ConstPointIterator>
  void insert(ConstPointIterator first,
              ConstPointIterator beyond)
  {
    while(first != beyond)
    {
      insert(*first); //update_bbox done by insert(BboxPrimitive)
      ++first;
    }
  }

  template<class OutputIterator, class FuzzyQueryItem>
  OutputIterator
  search(OutputIterator oit, const FuzzyQueryItem& q) const
  {
    std::size_t N = m_insertion_buffer.size();
    for(std::size_t i = 0; i < N; i++)
      if(q.contains(m_insertion_buffer[i]))
        *oit++ = m_insertion_buffer[i];

    return Base::search(oit,q);
  }
  
}; //end class Kd_tree_for_star_set

}//end namespace CGAL

#endif //CGAL_KD_TREE_STAR_SET_H

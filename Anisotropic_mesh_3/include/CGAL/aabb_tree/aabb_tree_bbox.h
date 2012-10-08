#ifndef CGAL_AABB_TREE_BBOX_H
#define CGAL_AABB_TREE_BBOX_H

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

#include <CGAL/aabb_tree/aabb_tree_bbox_primitive.h>

namespace CGAL
{
  template<typename GeomTraits, typename Star>
  class AABB_tree_bbox :
    public CGAL::AABB_tree< CGAL::AABB_traits< GeomTraits,
                                        CGAL::AABB_bbox_primitive<Star> > >
  {
    typedef AABB_tree_bbox<GeomTraits, Star>            Self;

    typedef typename CGAL::AABB_bbox_primitive<Star>             BboxPrimitive;
    typedef typename CGAL::AABB_traits<GeomTraits,BboxPrimitive> AABBtraits; 
    typedef typename CGAL::AABB_tree<AABBtraits>                 Base;
    typedef Star*                                                Star_handle;

  private:
    unsigned int m_insertion_buffer_max_size;
    std::vector<Star_handle> m_insertion_buffer;

#ifdef USE_ANISO_TIMERS
    mutable double m_rebuild_timer;
    mutable double m_query_timer;
    mutable int m_rebuild_counter;
    mutable bool m_rebuild_needed;
#endif
    
  public:
    AABB_tree_bbox() :
        m_insertion_buffer_max_size(20),
        m_insertion_buffer()
#ifdef USE_ANISO_TIMERS
        , m_rebuild_timer(0.)
        , m_query_timer(0.)
        , m_rebuild_counter(0)
        , m_rebuild_needed(false)
#endif
        { };

    AABB_tree_bbox(const unsigned int bs) :
        m_insertion_buffer_max_size(bs),
        m_insertion_buffer()
#ifdef USE_ANISO_TIMERS
        , m_rebuild_timer(0.)
        , m_query_timer(0.)
        , m_rebuild_counter(0)
        , m_rebuild_needed(false)
#endif
        { };
        
  /// Insertion
  public:
    void insert(const BboxPrimitive& p)
    {
      Star_handle pstar = p.id();
      pstar->update_bbox();

      if(m_insertion_buffer.size() >= m_insertion_buffer_max_size)
      {
        insert_buffer();
#ifdef USE_ANISO_TIMERS
        m_rebuild_needed = true;
#endif
      }
      m_insertion_buffer.push_back(pstar);
      
//      if(m_insertion_buffer.size() < m_insertion_buffer_max_size)
//        m_insertion_buffer.push_back(pstar);
//      else
//      {
//        insert_buffer();
//        Base::insert(BboxPrimitive(pstar));
//#ifdef USE_ANISO_TIMERS
//        m_rebuild_needed = true;
//#endif
//      }
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
      for(i = 0; i < N; i++)//update_bbox done when inserted to buffer
        Base::insert(m_insertion_buffer[i]/*Star_handle*/);

      m_insertion_buffer.clear();
    }

    template<typename ConstPrimitiveIterator>
    void insert(ConstPrimitiveIterator first, //Star_handle_iterator
                ConstPrimitiveIterator beyond)
    {
      while(first != beyond)
      {
        insert(BboxPrimitive(*first)); //update_bbox done by insert(BboxPrimitive)
        ++first;
      }
    }
    
    template<typename ConstPrimitiveIterator>
    void rebuild(ConstPrimitiveIterator first, //Star_handle_iterator
                 ConstPrimitiveIterator beyond)
    {
#ifdef USE_ANISO_TIMERS
     std::clock_t start = clock();
#endif

     Base::clear();
     m_insertion_buffer.clear(); //before things are added by insert(f,b)
     this->insert(first, beyond);
     insert_buffer(); // rebuild only with first-beyond data
     this->build();

#ifdef USE_ANISO_TIMERS
      m_rebuild_timer += (clock()-start+0.) / ((double)CLOCKS_PER_SEC);
      m_rebuild_counter++;
#endif
    }

    
    void clear()
    {
      Base::clear();
      m_insertion_buffer.clear();
    }

    /// Traversal
  private:
    template <class Query, class Traversal_traits>
    void traversal(const Query& query,
                   Traversal_traits& traits) const
    {
      if(this->empty() && m_insertion_buffer.empty())
        return;
      
      // first traverse the buffer
      std::size_t N = m_insertion_buffer.size();
      std::size_t i = 0;
      while(traits.go_further() && i < N)
      {
        CGAL_PROFILER("[intersections Bbox-Pt in buffer]");
        traits.intersection(query, m_insertion_buffer[i]);
        i++;      
      }
      // then use Base::traversal if not already ended in the buffer
      if(!this->empty() && traits.go_further())
      {
        CGAL_PROFILER("[Base traversals]");
        Base::traversal(query, traits);

#ifdef USE_ANISO_TIMERS
        if(m_rebuild_needed)
        {
          ++m_rebuild_counter;
          m_rebuild_needed = false; //done in traversal
        }
#endif
      }
    }

#ifdef USE_ANISO_TIMERS
public:
  void report_timers() const
  {
    std::cout << "\tAABBtree traversal: " << m_query_timer << " sec." << std::endl;
    std::cout << "\tAABBtree rebuild  : " //<< m_rebuild_timer << " sec"
      << "done " << m_rebuild_counter << " times."<< std::endl;
  }
#endif

  public:
    template<typename Query, typename OutputIterator>
    OutputIterator all_intersected_primitives(const Query& query,
		                              OutputIterator out) const
    {
#ifdef USE_ANISO_TIMERS
     std::clock_t start = clock();
#endif
      using namespace CGAL::internal::AABB_tree;
      Listing_primitive_traits<AABBtraits, Query, OutputIterator> traversal_traits(out);
      this->traversal(query, traversal_traits);

#ifdef USE_ANISO_TIMERS
      m_query_timer += (clock()-start+0.) / ((double)CLOCKS_PER_SEC);
#endif
      return out;
    }
    
  }; // end class AABB_tree_bbox

} // end namespace CGAL

#endif // CGAL_AABB_TREE_BBOX_H

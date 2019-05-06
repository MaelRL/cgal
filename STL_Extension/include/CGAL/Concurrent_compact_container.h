// Copyright (c) 2012  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Clement Jamin

#ifdef CGAL_LINKED_WITH_TBB

#ifndef CGAL_CONCURRENT_COMPACT_CONTAINER_H
#define CGAL_CONCURRENT_COMPACT_CONTAINER_H

#include <CGAL/disable_warnings.h>

#include <CGAL/basic.h>
#include <CGAL/Default.h>

#include <iterator>
#include <algorithm>
#include <vector>
#include <cstring>
#include <cstddef>

#include <CGAL/Compact_container.h>

#include <CGAL/memory.h>
#include <CGAL/iterator.h>
#include <CGAL/CC_safe_handle.h>
#include <CGAL/Time_stamper.h>
#include <CGAL/atomic.h>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/queuing_mutex.h>

#include <boost/mpl/if.hpp>

namespace CGAL {

#define CGAL_GENERATE_MEMBER_DETECTOR(X)                                           \
template<typename T> class has_##X {                                          \
    struct Fallback { int X; };                                               \
    struct Derived : T, Fallback { };                                         \
                                                                              \
    template<typename U, U> struct Check;                                     \
                                                                              \
    typedef char ArrayOfOne[1];                                               \
    typedef char ArrayOfTwo[2];                                               \
                                                                              \
    template<typename U> static ArrayOfOne & func(                            \
                                            Check<int Fallback::*, &U::X> *); \
    template<typename U> static ArrayOfTwo & func(...);                       \
  public:                                                                     \
    typedef has_##X type;                                                     \
    enum { value = sizeof(func<Derived>(0)) == 2 };                           \
} // semicolon is after the macro call

#define CGAL_INIT_CONCURRENT_COMPACT_CONTAINER_BLOCK_SIZE 14
#define CGAL_INCREMENT_CONCURRENT_COMPACT_CONTAINER_BLOCK_SIZE 16

// The traits class describes the way to access the pointer.
// It can be specialized.
template < class T >
struct Concurrent_compact_container_traits {
  static void *   pointer(const T &t) { return t.for_compact_container(); }
  static void * & pointer(T &t)       { return t.for_compact_container(); }
};

namespace CCC_internal {
  CGAL_GENERATE_MEMBER_DETECTOR(increment_erase_counter);
  
  // A basic "no erase counter" strategy
  template <bool Has_erase_counter_tag>
  class Erase_counter_strategy {
  public:
    // Do nothing
    template <typename Element>
    static unsigned int erase_counter(const Element &) { return 0; }
    template <typename Element>
    static void set_erase_counter(Element &, unsigned int) {}
    template <typename Element>
    static void increment_erase_counter(Element &) {}
  };


  // A strategy managing an internal counter
  template <>
  class Erase_counter_strategy<true>
  {
  public:
    template <typename Element>
    static unsigned int erase_counter(const Element &e)
    {
      return e.erase_counter();
    }

    template <typename Element>
    static void set_erase_counter(Element &e, unsigned int c)
    {
      e.set_erase_counter(c);
    }

    template <typename Element>
    static void increment_erase_counter(Element &e)
    {
      e.increment_erase_counter();
    }
  };
}

// Free list (head and size)
template< typename pointer, typename size_type, typename CCC >
class Free_list {
public:
  Free_list() : m_head(NULL), m_size(0) {}

  void init()                { m_head = NULL; m_size = 0; }
  pointer head() const       { return m_head; }
  void set_head(pointer p)   { m_head = p; }
  size_type size() const     { return m_size; }
  void set_size(size_type s) { m_size = s; }
  void inc_size()            { ++m_size; }
  void dec_size()            { --m_size; }
  bool empty()               { return size() == 0; }
  // Warning: copy the pointer, not the data!
  Free_list& operator= (const Free_list& other)
  {
    m_head = other.m_head;
    m_size = other.m_size;
    return *this;
  }

  void merge(Free_list &other)
  {
    if (m_head == NULL) {
      *this = other;
    }
    else if (!other.empty())
    {
      pointer p = m_head;
      while (CCC::clean_pointee(p) != NULL)
        p = CCC::clean_pointee(p);
      CCC::set_type(p, other.m_head, CCC::FREE);
      m_size += other.m_size;
    }
    other.init(); // clear other
  }

protected:
  pointer   m_head;  // the free list head pointer
  size_type m_size;  // the free list size
};

// Class Concurrent_compact_container
//
// Safe concurrent "insert" and "erase".
// Do not parse the container while others are modifying it.
//
template < class T, class Allocator_ = Default >
class Concurrent_compact_container
{
  typedef Allocator_                                                Al;
  typedef typename Default::Get<Al, CGAL_ALLOCATOR(T) >::type       Allocator;
  typedef Concurrent_compact_container <T, Al>                      Self;
  typedef Concurrent_compact_container_traits <T>                   Traits;

public:
  typedef CGAL::Time_stamper_impl<T>                Time_stamper_impl;

  typedef T                                         value_type;
  typedef Allocator                                 allocator_type;

  typedef value_type&                               reference;
  typedef const value_type&                         const_reference;

#ifdef CGAL_CXX11
  typedef typename std::allocator_traits<Allocator>::pointer               pointer;
  typedef typename std::allocator_traits<Allocator>::const_pointer         const_pointer;
  typedef typename std::allocator_traits<Allocator>::size_type             size_type;
  typedef typename std::allocator_traits<Allocator>::difference_type       difference_type;
#else
  typedef typename Allocator::pointer               pointer;
  typedef typename Allocator::const_pointer         const_pointer;
  typedef typename Allocator::size_type             size_type;
  typedef typename Allocator::difference_type       difference_type;
#endif

  typedef internal::CC_iterator<Self, false>        iterator;
  typedef internal::CC_iterator<Self, true>         const_iterator;
  typedef std::reverse_iterator<iterator>           reverse_iterator;
  typedef std::reverse_iterator<const_iterator>     const_reverse_iterator;

private:
  typedef Free_list<pointer, size_type, Self>       FreeList;
  typedef tbb::enumerable_thread_specific<FreeList> Free_lists;

  // FreeList can access our private function (clean_pointee...)
  friend class Free_list<pointer, size_type, Self>;

public:
  friend class internal::CC_iterator<Self, false>;
  friend class internal::CC_iterator<Self, true>;

  explicit Concurrent_compact_container(const Allocator &a = Allocator())
  : m_alloc(a)
  , m_time_stamper(new Time_stamper_impl())
  {
    init ();
  }

  template < class InputIterator >
  Concurrent_compact_container(InputIterator first, InputIterator last,
                    const Allocator & a = Allocator())
  : m_alloc(a)
  , m_time_stamper(new Time_stamper_impl())
  {
    init();
    std::copy(first, last, CGAL::inserter(*this));
  }

  // The copy constructor and assignment operator preserve the iterator order
  Concurrent_compact_container(const Concurrent_compact_container &c)
  : m_alloc(c.get_allocator())
  , m_time_stamper(new Time_stamper_impl())
  {
    init();
    m_block_size = c.m_block_size;
    std::copy(c.begin(), c.end(), CGAL::inserter(*this));
  }

  Concurrent_compact_container & operator=(const Concurrent_compact_container &c)
  {
    if (&c != this) {
      Self tmp(c);
      swap(tmp);
    }
    return *this;
  }

  ~Concurrent_compact_container()
  {
    clear();
    delete m_time_stamper;
  }

  bool is_used(const_iterator ptr) const
  {
    return (type(&*ptr)==USED);
  }

  void swap(Self &c)
  {
    std::swap(m_alloc, c.m_alloc);
    std::swap(m_capacity, c.m_capacity);
    { // non-atomic swap
      size_type other_size = c.m_size;
      c.m_size = size_type(m_size);
      m_size = other_size;
    }
    std::swap(m_block_size, c.m_block_size);
    std::swap(m_first_item, c.m_first_item);
    std::swap(m_last_item, c.m_last_item);
    std::swap(m_free_lists, c.m_free_lists);
    m_all_items.swap(c.m_all_items);
    std::swap(m_time_stamper, c.m_time_stamper);
  }

  iterator begin() { return iterator(m_first_item, 0, 0); }
  iterator end()   { return iterator(m_last_item, 0); }

  const_iterator begin() const { return const_iterator(m_first_item, 0, 0); }
  const_iterator end()   const { return const_iterator(m_last_item, 0); }

  reverse_iterator rbegin() { return reverse_iterator(end()); }
  reverse_iterator rend()   { return reverse_iterator(begin()); }

  const_reverse_iterator
  rbegin() const { return const_reverse_iterator(end()); }
  const_reverse_iterator
  rend()   const { return const_reverse_iterator(begin()); }

  // Boost.Intrusive interface
  iterator iterator_to(reference value) const {
    return iterator(&value, 0);
  }
  const_iterator iterator_to(const_reference value) const {
    return const_iterator(&value, 0);
  }
  static iterator s_iterator_to(reference value) {
    return iterator(&value, 0);
  }
  static const_iterator s_iterator_to(const_reference value) {
    return const_iterator(&value, 0);
  }

  // Special insert methods that construct the objects in place
  // (just forward the arguments to the constructor, to optimize a copy).
  template < typename... Args >
  iterator
  emplace(const Args&... args)
  {
    FreeList * fl = get_free_list();
    pointer ret = init_insert(fl);
    new (ret) value_type(args...);
    return finalize_insert(ret, fl);
  }

  iterator insert(const T &t)
  {
    FreeList * fl = get_free_list();
    pointer ret = init_insert(fl);
    std::allocator_traits<allocator_type>::construct(m_alloc, ret, t);
    return finalize_insert(ret, fl);
  }

  template < class InputIterator >
  void insert(InputIterator first, InputIterator last)
  {
    for (; first != last; ++first)
      insert(*first);
  }

  template < class InputIterator >
  void assign(InputIterator first, InputIterator last)
  {
    clear(); // erase(begin(), end()); // ?
    insert(first, last);
  }

private:
  void erase(iterator x, FreeList * fl)
  {
    typedef CCC_internal::Erase_counter_strategy<
      CCC_internal::has_increment_erase_counter<T>::value> EraseCounterStrategy;

    CGAL_precondition(type(x) == USED);
    EraseCounterStrategy::increment_erase_counter(*x);

    std::allocator_traits<allocator_type>::destroy(m_alloc, &*x);

/* WE DON'T DO THAT BECAUSE OF THE ERASE COUNTER
#ifndef CGAL_NO_ASSERTIONS
    std::memset(&*x, 0, sizeof(T));
#endif*/
    --m_size;
    put_on_free_list(&*x, fl);
  }
public:

  void erase(iterator x)
  {
    erase(x, get_free_list());
  }

  void erase(iterator first, iterator last) {
    while (first != last)
      erase(first++);
  }

  void clear();

  // Merge the content of d into *this.  d gets cleared.
  // The complexity is O(size(free list = capacity-size)).
  void merge(Self &d);

  // If `CGAL_NO_ATOMIC` is defined, do not call this function while others
  // are inserting/erasing elements
  size_type size() const
  {
#ifdef CGAL_NO_ATOMIC
    size_type size = m_capacity;
    for( typename Free_lists::iterator it_free_list = m_free_lists.begin() ;
         it_free_list != m_free_lists.end() ;
         ++it_free_list )
    {
      size -= it_free_list->size();
    }
    return size;
#else // atomic can be used
    return m_size;
#endif
  }

  size_type max_size() const
  {
#ifdef CGAL_CXX11
    return std::allocator_traits<allocator_type>::max_size(m_alloc);
#else
    return m_alloc.max_size();
#endif
  }

  size_type capacity() const
  {
    return m_capacity;
  }

  // void resize(size_type sz, T c = T()); // TODO  makes sense ???

  bool empty() const
  {
    return size() == 0;
  }

  allocator_type get_allocator() const
  {
    return m_alloc;
  }

  // Returns whether the iterator "cit" is in the range [begin(), end()].
  // Complexity : O(#blocks) = O(sqrt(capacity())).
  // This function is mostly useful for purposes of efficient debugging at
  // higher levels.
  bool owns(const_iterator cit) const
  {
    // We use the block structure to provide an efficient version :
    // we check if the address is in the range of each block,
    // and then test whether it is valid (not a free element).

    if (cit == end())
      return true;

    const_pointer c = &*cit;

    Mutex::scoped_lock lock(m_mutex);

    for (typename All_items::const_iterator it = m_all_items.begin(), itend = m_all_items.end();
         it != itend; ++it) {
      const_pointer p = it->first;
      size_type s = it->second;

      // Are we in the address range of this block (excluding first and last
      // elements) ?
      if (c <= p || (p+s-1) <= c)
        continue;

      CGAL_assertion_msg( (c-p)+p == c, "wrong alignment of iterator");

      return type(c) == USED;
    }
    return false;
  }

  bool owns_dereferencable(const_iterator cit) const
  {
    return cit != end() && owns(cit);
  }

  /** Reserve method to ensure that the capacity of the Concurrent_compact_container be
   * greater or equal than a given value n.
   */
  // TODO?
  //void reserve(size_type n)
  //{
    // Does it really make sense: it will reserve size for the current
    // thread only!
    /*Mutex::scoped_lock lock;
    if ( m_capacity >= n ) return;
    size_type tmp = m_block_size;
    // TODO: use a tmpBlockSize instead of m_block_size
    m_block_size = (std::max)( n - m_capacity, m_block_size );
    allocate_new_block(free_list());
    m_block_size = tmp + CGAL_INCREMENT_CONCURRENT_COMPACT_CONTAINER_BLOCK_SIZE;*/
  //}

private:

  FreeList*       get_free_list()       { return & m_free_lists.local(); }
  const FreeList* get_free_list() const { return & m_free_lists.local(); }

  // Two helper functions for the emplace() methods

  // allocate new space if needed get the pointer from
  // the free list and then clean it
  pointer init_insert(FreeList * fl)
  {
    pointer fl2 = fl->head();
    if (fl2 == NULL) {
      allocate_new_block(fl);
      fl2 = fl->head();
    }
    pointer ret = fl2;
    fl->set_head(clean_pointee(ret));
    return ret;
  }

  // get verify the return pointer increment size and
  // return as iterator
  iterator finalize_insert(pointer ret, FreeList * fl)
  {
    CGAL_assertion(type(ret) == USED);
    fl->dec_size();
    ++m_size;
    m_time_stamper->set_time_stamp(ret);
    return iterator(ret, 0);
  }

  void allocate_new_block(FreeList *fl);

  void put_on_free_list(pointer x, FreeList * fl)
  {
    set_type(x, fl->head(), FREE);
    fl->set_head(x);
    fl->inc_size();
  }

  // Definition of the bit squatting :
  // =================================
  // ptr is composed of a pointer part and the last 2 bits.
  // Here is the meaning of each of the 8 cases.
  //
  //                          value of the last 2 bits as "Type"
  // pointer part     0              1                2              3
  //         NULL     user elt       unused           free_list end  start/end
  //      != NULL     user elt       block boundary   free elt       unused
  //
  // meaning of ptr : user stuff     next/prev block  free_list      unused

  enum Type { USED = 0, BLOCK_BOUNDARY = 1, FREE = 2, START_END = 3 };

  // The bit squatting is implemented by casting pointers to (char *), then
  // subtracting to NULL, doing bit manipulations on the resulting integer,
  // and converting back.

  static char * clean_pointer(char * p)
  {
    return reinterpret_cast<char*>(reinterpret_cast<std::ptrdiff_t>(p) &
                                   ~ (std::ptrdiff_t) START_END);
  }

  // Returns the pointee, cleaned up from the squatted bits.
  static pointer clean_pointee(const_pointer ptr)
  {
    return (pointer) clean_pointer((char *) Traits::pointer(*ptr));
  }

  // Get the type of the pointee.
  static Type type(const_pointer ptr)
  {
    char * p = (char *) Traits::pointer(*ptr);
    return (Type) (reinterpret_cast<std::ptrdiff_t>(p) -
                   reinterpret_cast<std::ptrdiff_t>(clean_pointer(p)));
  }

  static Type type(const_iterator it)
  {
    return type(it.operator->());
  }

  // Sets the pointer part and the type of the pointee.
  static void set_type(pointer ptr, void * p, Type t)
  {
    // This out of range compare is always true and causes lots of
    // unnecessary warnings.
    // CGAL_precondition(0 <= t && t < 4);
    Traits::pointer(*ptr) = reinterpret_cast<void *>
      (reinterpret_cast<std::ptrdiff_t>(clean_pointer((char *) p)) + (int) t);
  }

  typedef tbb::queuing_mutex Mutex;

  // We store a vector of pointers to all allocated blocks and their sizes.
  // Knowing all pointers, we don't have to walk to the end of a block to reach
  // the pointer to the next block.
  // Knowing the sizes allows to deallocate() without having to compute the size
  // by walking through the block till its end.
  // This opens up the possibility for the compiler to optimize the clear()
  // function considerably when has_trivial_destructor<T>.
  typedef std::vector<std::pair<pointer, size_type> >  All_items;


  void init()
  {
    m_block_size = CGAL_INIT_CONCURRENT_COMPACT_CONTAINER_BLOCK_SIZE;
    m_capacity  = 0;
    for( typename Free_lists::iterator it_free_list = m_free_lists.begin() ;
         it_free_list != m_free_lists.end() ;
         ++it_free_list )
    {
      it_free_list->set_head(0);
      it_free_list->set_size(0);
    }
    m_first_item = NULL;
    m_last_item  = NULL;
    m_all_items  = All_items();
    m_size = 0;
    m_time_stamper->reset();
  }

  allocator_type    m_alloc;
  size_type         m_capacity;
  size_type         m_block_size;
  Free_lists        m_free_lists;
  pointer           m_first_item;
  pointer           m_last_item;
  All_items         m_all_items;
  mutable Mutex     m_mutex;
#ifdef CGAL_NO_ATOMIC
  size_type         m_size;
#else
  CGAL::cpp11::atomic<size_type> m_size;
#endif

  // This is a pointer, so that the definition of Compact_container does
  // not require a complete type `T`.
  Time_stamper_impl* m_time_stamper;
};

template < class T, class Allocator >
void Concurrent_compact_container<T, Allocator>::merge(Self &d)
{
  m_size += d.m_size;
  CGAL_precondition(&d != this);

  // Allocators must be "compatible" :
  CGAL_precondition(get_allocator() == d.get_allocator());

  // Concatenate the free_lists.
  // Iterates over TLS free lists of "d". Note that the number of TLS freelists
  // may be different.
  typename Free_lists::iterator it_free_list = m_free_lists.begin();
  if (it_free_list == m_free_lists.end())
  {
    // No free list at all? Create our local one... empty.
    get_free_list()->set_head(0);
    get_free_list()->set_size(0);
    // Now there is one TLS free list: ours!
    it_free_list = m_free_lists.begin();
  }
  for( typename Free_lists::iterator it_free_list_d = d.m_free_lists.begin() ;
       it_free_list_d != d.m_free_lists.end() ;
       ++it_free_list_d, ++it_free_list )
  {
    // If we run out of TLS free lists in *this, let's start again from "begin"
    if (it_free_list == m_free_lists.end())
      it_free_list = m_free_lists.begin();

    it_free_list->merge(*it_free_list_d);
  }
  // Concatenate the blocks.
  if (m_last_item == NULL) { // empty...
    m_first_item = d.m_first_item;
    m_last_item  = d.m_last_item;
  } else if (d.m_last_item != NULL) {
    set_type(m_last_item, d.m_first_item, BLOCK_BOUNDARY);
    set_type(d.m_first_item, m_last_item, BLOCK_BOUNDARY);
    m_last_item = d.m_last_item;
  }
  m_all_items.insert(m_all_items.end(), d.m_all_items.begin(), d.m_all_items.end());
  // Add the capacities.
  m_capacity += d.m_capacity;
  // It seems reasonnable to take the max of the block sizes.
  m_block_size = (std::max)(m_block_size, d.m_block_size);
  // Clear d.
  d.init();
}

template < class T, class Allocator >
void Concurrent_compact_container<T, Allocator>::clear()
{
  for (typename All_items::iterator it = m_all_items.begin(), itend = m_all_items.end();
       it != itend; ++it) {
    pointer p = it->first;
    size_type s = it->second;
    for (pointer pp = p + 1; pp != p + s - 1; ++pp) {
      if (type(pp) == USED)
        m_alloc.destroy(pp);
    }
    m_alloc.deallocate(p, s);
  }
  init();
}

template < class T, class Allocator >
void Concurrent_compact_container<T, Allocator>::
  allocate_new_block(FreeList * fl)
{
  typedef CCC_internal::Erase_counter_strategy<
    CCC_internal::has_increment_erase_counter<T>::value> EraseCounterStrategy;

  size_type old_block_size;
  pointer new_block;

  {
    Mutex::scoped_lock lock(m_mutex);
    old_block_size = m_block_size;
    new_block = m_alloc.allocate(old_block_size + 2);
    m_all_items.push_back(std::make_pair(new_block, old_block_size + 2));
    m_capacity += old_block_size;

    // We insert this new block at the end.
    if (m_last_item == NULL) // First time
    {
        m_first_item = new_block;
        m_last_item  = new_block + old_block_size + 1;
        set_type(m_first_item, NULL, START_END);
    }
    else
    {
        set_type(m_last_item, new_block, BLOCK_BOUNDARY);
        set_type(new_block, m_last_item, BLOCK_BOUNDARY);
        m_last_item = new_block + old_block_size + 1;
    }
    set_type(m_last_item, NULL, START_END);
    // Increase the m_block_size for the next time.
    m_block_size += CGAL_INCREMENT_CONCURRENT_COMPACT_CONTAINER_BLOCK_SIZE;
  }

  // We don't touch the first and the last one.
  // We mark them free in reverse order, so that the insertion order
  // will correspond to the iterator order...
  for (size_type i = old_block_size; i >= 1; --i)
  {
    EraseCounterStrategy::set_erase_counter(*(new_block + i), 0);
    m_time_stamper->initialize_time_stamp(new_block + i);
    put_on_free_list(new_block + i, fl);
  }
}

template < class T, class Allocator >
inline
bool operator==(const Concurrent_compact_container<T, Allocator> &lhs,
                const Concurrent_compact_container<T, Allocator> &rhs)
{
  return lhs.size() == rhs.size() &&
    std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

template < class T, class Allocator >
inline
bool operator!=(const Concurrent_compact_container<T, Allocator> &lhs,
                const Concurrent_compact_container<T, Allocator> &rhs)
{
  return ! (lhs == rhs);
}

template < class T, class Allocator >
inline
bool operator< (const Concurrent_compact_container<T, Allocator> &lhs,
                const Concurrent_compact_container<T, Allocator> &rhs)
{
  return std::lexicographical_compare(lhs.begin(), lhs.end(),
                                      rhs.begin(), rhs.end());
}

template < class T, class Allocator >
inline
bool operator> (const Concurrent_compact_container<T, Allocator> &lhs,
                const Concurrent_compact_container<T, Allocator> &rhs)
{
  return rhs < lhs;
}

template < class T, class Allocator >
inline
bool operator<=(const Concurrent_compact_container<T, Allocator> &lhs,
                const Concurrent_compact_container<T, Allocator> &rhs)
{
  return ! (lhs > rhs);
}

template < class T, class Allocator >
inline
bool operator>=(const Concurrent_compact_container<T, Allocator> &lhs,
                const Concurrent_compact_container<T, Allocator> &rhs)
{
  return ! (lhs < rhs);
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_CONCURRENT_COMPACT_CONTAINER_H

#endif // CGAL_LINKED_WITH_TBB

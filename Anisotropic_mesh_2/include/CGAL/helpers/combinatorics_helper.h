#ifndef COMBINATORICS_HELPER_H
#define COMBINATORICS_HELPER_H

#include <boost/array.hpp>
#include <boost/functional/hash.hpp>
#include <boost/unordered_set.hpp>

// helper for combinatorics

class Edge_ij;
class Facet_ijk;

template<int count>
inline bool is_same_ids(int *cids, int *dids)
{
  int t;
  for (int i = count - 1; i >= 0; i--)
  {
    for (int j = 1; j <= i; j++) //bubble sort + check on top bubble
    {
      if (cids[j] < cids[j - 1])
      {
        t = cids[j];
        cids[j] = cids[j - 1];
        cids[j - 1] = t;
      }
      if (dids[j] < dids[j - 1])
      {
        t = dids[j];
        dids[j] = dids[j - 1];
        dids[j - 1] = t;
      }
    }
    if (cids[i] != dids[i])
      return false;
  }
  return true;
}

//template<typename Cell_handle>
//inline bool is_same(const Cell_handle &c, int *vertices)
//{
//  int cids[4];
//  for (int i = 0; i < 4; i++)
//    cids[i] = c->vertex(i)->info();
//  return is_same_ids<4>(cids, vertices);
//}

template<typename Element, typename OneMap>
void add_to_map(const Element& e, OneMap& m)
{
  std::pair<typename OneMap::iterator, bool> is_insert_successful;
  is_insert_successful = m.insert(std::make_pair(e,1));
  if(!is_insert_successful.second)
    (is_insert_successful.first)->second += 1; // m[e] += 1
}

template<std::size_t ssize>
class Ordered_simplex_base
{
  typedef Ordered_simplex_base<ssize> Self;

  boost::array<int, ssize> m_simplex;

public :
  std::size_t size() const { return ssize; }
  void sort() { std::sort(m_simplex.begin(), m_simplex.end()); }
  boost::array<int, ssize>& vertices() { return m_simplex; }
  const boost::array<int, ssize>& vertices() const { return m_simplex; }

  int vertex(const std::size_t index) const
  {
    if(index >= ssize || index < 0)
    {
      std::cerr << "Ordered_simplex does not have vertex number " << index << "!\n";
      return -1;
    }
    else
      return m_simplex[index];
  }

  bool has(const int index) const
  {
    for(std::size_t i=0; i<ssize; ++i)
      if(vertex(i) == index)
        return true;
    return false;
  }

  bool has(const int index, boost::array<int, ssize-1>& others) const
  {
    for(std::size_t i=0; i<ssize; ++i)
      if(vertex(i) == index)
      {
        for(std::size_t j=1; j<ssize; ++j)
          others[j-1] = vertex((i+j)%ssize);
        return true;
      }
    return false;
  }

  int operator[](const std::size_t index) const
  {
    if(index >= ssize || index < 0)
    {
      std::cerr << "Ordered_simplex does not have vertex number " << index << "!\n";
      return -1;
    }
    else
      return m_simplex[index];
  }

  bool operator==(const Self& os) const
  {
    if(ssize != os.size())
      return false;
    for(std::size_t i=0; i<ssize; ++i)
      if(vertex(i) != os.vertex(i))
        return false;
    return true;
  }

  bool operator<(const Self& os) const
  {
    for(std::size_t i=0; i<ssize; ++i)
      if (vertex(i) != os.vertex(i))
        return (vertex(i) < os.vertex(i));
    return false;
  }

  Self& operator=(const Self& os)
  {
    for(std::size_t i=0; i<ssize; ++i)
      m_simplex[i] = os.vertex(i);
    return *this;
  }

  bool is_infinite() const
  {
    return has(-10);
  }

  bool is_valid()
  {
    return !has(-1);
  }

  Ordered_simplex_base() { std::fill(m_simplex.begin(), m_simplex.end(),-1); }

  Ordered_simplex_base(const boost::array<int, ssize>& ids)
  {
    std::copy(ids.begin(), ids.end(), m_simplex.begin());
    sort();
  }
};

class Edge_ij : public Ordered_simplex_base<2>
{
public:
  Edge_ij(){}
  Edge_ij(int a, int b)
  {
    vertices()[0] = (std::min)(a,b);
    vertices()[1] = (std::max)(a,b);
  }

  template<typename Tr_simplex>
  Edge_ij(const Tr_simplex &tr_s)
  {
    int i = tr_s.first->vertex(tr_s.second)->info();
    int j = tr_s.first->vertex(tr_s.third)->info();
    vertices()[0] = (std::min)(i,j);
    vertices()[1] = (std::max)(i,j);
  }
};

class Facet_ijk : public Ordered_simplex_base<3>
{
public:

  Facet_ijk(){}
  Facet_ijk(int a, int b, int c)
  {
    vertices()[0] = a;
    vertices()[1] = b;
    vertices()[2] = c;
    sort();
  }

  template<typename Tr_simplex>
  Facet_ijk(const Tr_simplex &tr_s)
  {
    vertices()[0] = tr_s->vertex(0)->info();
    vertices()[1] = tr_s->vertex(1)->info();
    vertices()[2] = tr_s->vertex(2)->info();
    sort();
  }
};

template<std::size_t ssize>
class Ordered_simplex_base_hash
{
public:
  typedef Ordered_simplex_base<ssize>           Simplex;

  Ordered_simplex_base_hash(){}

  std::size_t operator()(const Simplex& s) const
  {
    return boost::hash_range(s.vertices().begin(),
                             s.vertices().end());
  }
};

template<std::size_t ssize>
class Ordered_simplex_base_comparer
{
public:
  typedef Ordered_simplex_base<ssize>           Simplex;

  Ordered_simplex_base_comparer() { }
  bool operator() (const Simplex& left, const Simplex& right) const
  {
    return std::equal(left.vertices().begin(), left.vertices().end(),
                      right.vertices().begin());
  }
};

typedef typename boost::unordered_set<Edge_ij,
                                      Ordered_simplex_base_hash<2>,
                                      Ordered_simplex_base_comparer<2> >
                                                    Edge_ij_unordered_set;
typedef typename boost::unordered_set<Facet_ijk,
                                      Ordered_simplex_base_hash<3>,
                                      Ordered_simplex_base_comparer<3> >
                                                    Facet_ijk_unordered_set;

#endif

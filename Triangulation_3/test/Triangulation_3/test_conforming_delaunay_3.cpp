bool del = false;

#include <CGAL/_test_types.h>
#include <CGAL/_test_cls_delaunay_3.h>
#include <CGAL/_test_cls_conforming_triangulation_3.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Conforming_Delaunay_triangulation_3.h>

#include <cassert>

template <typename Gt>
using Vb = CGAL::Conforming_Delaunay_triangulation_vertex_base_3<Gt>;

template <typename Gt>
using Cb = CGAL::Delaunay_triangulation_cell_base_3<Gt>;

template <typename Gt>
using Tds = CGAL::Triangulation_data_structure_3<Vb<Gt>, Cb<Gt> >;

template <typename Gt>
using Dt3 = CGAL::Delaunay_triangulation_3<Gt, Tds<Gt> >;

// Explicit instantiation of the whole class :
template class CGAL::Conforming_Delaunay_triangulation_3<Dt3<EPIC> >;

template <typename K>
void test_kernel()
{
  using Cls3 = CGAL::Conforming_Delaunay_triangulation_3<Dt3<K> >;

  // conforming as Delaunay
  _test_cls_delaunay_3(Cls3());

  _test_cls_conforming_delaunay_3(Cls3());

  // Test operator== between triangulations of different Tds types.
  using Vb_bis = CGAL::Triangulation_vertex_base_with_info_3<int, K>;
  using Cb_bis = CGAL::Triangulation_cell_base_3<K>;
  using Tds_bis = CGAL::Triangulation_data_structure_3<Vb_bis, Cb_bis>;

  typedef CGAL::Triangulation_3<K, Tds_bis> Cls3_2;
  assert(Cls3() == Cls3_2());
  std::cout << "done" << std::endl;
}

int main()
{
  test_kernel<EPIC>();
  test_kernel<EPEC>();

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
}

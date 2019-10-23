// Copyright (c) 2001,2002  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Sylvain Pion
 

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

int main()
{
  CGAL::Exact_predicates_inexact_constructions_kernel k;
  CGAL::Epick_without_intervals kiwi;

  std::cout << "typeid k: " << typeid(k).name() << std::endl;
  CGAL::Exact_predicates_inexact_constructions_kernel::Point_3 p, q, r, s;
  k.orientation_3_object()(p, q, r, s);

  std::cout << "typeid kiwi: " << typeid(kiwi).name() << std::endl;
  CGAL::Epick_without_intervals::Point_3 pi, qi, ri, si;
  kiwi.orientation_3_object()(pi, qi, ri, si);
}

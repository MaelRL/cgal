// Copyright (c) 2018 GeometryFactory Sarl (France).
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
// Author(s)     : Simon Giraudot

#ifndef CGAL_THREAD_H
#define CGAL_THREAD_H

#include <CGAL/config.h>

/*
  This file defines the following:
  - CGAL::cpp11::thread
  - CGAL::cpp11::this_thread::sleep_for
  - CGAL::cpp11::atomic
  - CGAL::cpp11::chrono::seconds

  It uses either TBB or STD depending on what's available: as TBB can
  quite often override `std::thread`, it is possible that TBB will be
  used instead of STD even if the real CXX11 `std::thread` is
  available.
*/

#if defined(CGAL_LINKED_WITH_TBB)
#  include <tbb/tbb_config.h>
#  if TBB_IMPLEMENT_CPP0X
#    include <tbb/compat/thread>
#    include <tbb/atomic.h>
#    include <tbb/tick_count.h>
#    define CGAL_USE_TBB_THREADS 1
#  else
#    define CGAL_USE_TBB_THREADS 0
#  endif
#else
#  define CGAL_USE_TBB_THREADS 0
#endif

#if !CGAL_USE_TBB_THREADS
#  if !defined(CGAL_CFG_NO_STD_THREAD)
#    include <thread>
#    include <atomic>
#    include <chrono>
#    define CGAL_USE_STD_THREADS 1
#  else
#    define CGAL_USE_STD_THREADS 0
#  endif
#else
#  define CGAL_USE_STD_THREADS 0
#endif

#if CGAL_USE_STD_THREADS || CGAL_USE_TBB_THREADS
#  define CGAL_HAS_STD_THREADS // Useful define
#endif

namespace CGAL {

namespace cpp11 {

#if CGAL_USE_TBB_THREADS
  using std::thread; // std::thread is declared by TBB if TBB_IMPLEMENT_CPP0X == 1
  namespace this_thread
  {
    using std::this_thread::sleep_for;  // std::this_thread::sleep_for is declared by TBB if TBB_IMPLEMENT_CPP0X == 1
  }
  using tbb::atomic;
  namespace chrono
  {
    typedef tbb::tick_count::interval_t seconds;
  }
#elif CGAL_USE_STD_THREADS
  using std::thread;
  namespace this_thread
  {
    using std::this_thread::sleep_for;
  }
  using std::atomic;
  namespace chrono
  {
    typedef std::chrono::duration<double, std::ratio<1> > seconds;
  }
#endif

} // cpp11

} //namespace CGAL

#undef CGAL_USE_STD_THREADS
#undef CGAL_USE_TBB_THREADS

#endif // CGAL_THREAD_H

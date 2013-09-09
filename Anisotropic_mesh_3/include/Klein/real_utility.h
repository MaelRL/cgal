// Copyright (c) 2011  INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Author(s) : Kan-Le Shi

#ifndef KLEIN_KERNEL_MATHS_REAL_UTILITY_INCLUDE
#define KLEIN_KERNEL_MATHS_REAL_UTILITY_INCLUDE

#include <cmath>
#include <cfloat>

namespace Klein
{

namespace Kernel
{

namespace Maths
{

#ifndef PI
#define PI			3.1415926535897932384626433832795
#endif

//#ifndef min
//#define min(a,b)	(((a)>(b))?(b):(a))
//#endif

//#ifndef max
//#define max(a,b)	(((a)>(b))?(a):(b))
//#endif

/**
 * @brief The type representing a real value.
 */
typedef double Real;

/**
 * @brief The utility functions of the Real type.
 */
class RealUtility
{
public:

	/**
	 * @brief Compare whether the two Real values are equal within the specified tolerance.
	 * @param x,y Real values to be compared.
	 * @param tolerance User specified tolerance. 
	 * @return true if equal.
	 */
	static inline bool areEqual(Real x, Real y, Real tolerance);

	/**
	 * @brief Compare a Real value with setZero within the specified tolerance.
	 * @param x Real value to be compared.
	 * @param tolerance User specified tolerance. 
	 * @return true if x equals to setZero.
	 */
	static inline bool isZero(Real x, Real tolerance);

	/**
	 * @brief Get the maximum value that can be represented by Real.
	 * @return Maximum Real value.
	 */
	static inline Real maxReal();
	/**
	 * @brief Get the minimum value that can be represented by Real.
	 * @return Minimum Real value.
	 */
	static inline Real minReal();
};


bool RealUtility::areEqual(Real x, Real y, Real tolerance)
{
	return fabs(x - y) <= tolerance;
}

bool RealUtility::isZero(Real x, Real tolerance)
{
	return (x < tolerance) && (x > -tolerance);
}

Real RealUtility::maxReal()
{
	return DBL_MAX;
}

Real RealUtility::minReal()
{
	return DBL_MIN;
}

}	// Maths

}	// Kernel

}	// Klein

#endif
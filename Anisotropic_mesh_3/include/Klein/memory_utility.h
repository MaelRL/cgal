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

#ifndef KLEIN_KERNEL_UTILITY_MEMORY_UTILITY_INCLUDE
#define KLEIN_KERNEL_UTILITY_MEMORY_UTILITY_INCLUDE

namespace Klein
{

namespace Kernel
{

namespace Utility
{

#ifndef NULL
#define NULL	0
#endif

/**
 * @brief Memory management utilities.
 */
class MemoryUtility
{
public:

	/**
	 * @brief Memory block copy.
	 * @param dest the pointer of the destination memory block.
	 * @param destSize the size, in bytes, of the destination memory block.
	 * @param src the pointer of the source memory block.
	 * @param srcSize the size, in bytes, of the source memory block.
	 * @return true if succeeds.
	 */
	static bool copy(void *dest, int destSize, const void *src, int srcSize) {
		char *dest_p = (char*)dest;
		char *src_p = (char*)src;
		for (int i = (std::min)(destSize, srcSize); i > 0; i--)
			*(dest_p++) = *(src_p++);
		return true;
	}

	/**
	 * @brief Set all the data in the specified memory block setZero.
	 * @param dest the pointer of the destination memory block.
	 * @param destSize the size, in bytes, of the destination memory block.
	 * @return true if succeeds.
	 */
	static void setZero(void *dest, int destSize) {
		char *p = (char *)dest;
		for (int i = 0; i < destSize; i++)
			*(p++) = 0;
	}
};

}	// Utility

}	// Kernel

}	// Klein

#endif	// KLEIN_KERNEL_UTILITY_MEMORY_UTILITY_INCLUDE

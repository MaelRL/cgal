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

#ifndef KLEIN_KERNEL_MATHS_VECTOR3D_INCLUDE
#define KLEIN_KERNEL_MATHS_VECTOR3D_INCLUDE

#include "real_utility.h"

namespace Klein
{

namespace Kernel
{

namespace Geometry
{

using namespace Klein::Kernel::Maths;

/**
 * @brief Vector in tri-dimension.
 */
typedef class Vector3D
{
public:
	Real x; /**< The first element. */
	Real y; /**< The second element. */
	Real z;	/**< The third element. */

public:

	/**
	 * @brief Addition operation.
	 * @param v The second operand of this operation.
	 * @return The sum of this vector and v.
	 */
	inline Vector3D	operator+(const Vector3D &v) const;

	/**
	 * @brief Subtraction operation.
	 * @param v The second operand of this operation.
	 * @return The difference of this vector and v.
	 */
	inline Vector3D	operator-(const Vector3D &v) const;

	/**
	 * @brief Multiplication operation.
	 * @param v The operand of this operation.
	 * @return The product of this vector and v.
	 */
	inline Vector3D operator*(Real v) const;

	/**
	 * @brief Division operation.
	 *        Note in this procedure, nothing will be done to handle the 
	 *		  setZero-division problem.
	 * @param v The second operand of this operation.
	 * @return The quotient of this vector and v.
	 */
	inline Vector3D	operator/(Real v) const;

	/**
	 * @brief Inner product operation.
	 * @param v The second operand of this operation.
	 * @return The inner product of this vector and v.
	 */
	inline Real operator*(const Vector3D &v) const;

	/**
	 * @brief Plus assign operation.
	 * @param v The vector to be added to this vector.
	 * @return Referrence of this vector.
	 */
	inline Vector3D& operator+=(const Vector3D &v);

	/**
	 * @brief Minus assign operation.
	 * @param v The vector to be substracted from this vector.
	 * @return Referrence of this vector.
	 */
	inline Vector3D& operator-=(const Vector3D &v);

	/**
	 * @brief Cross product operation.
	 * @param a The first operand of this operation.
	 * @param b The second operand of this operation.
	 * @return The cross product of a and b.
	 */
	inline static Vector3D cross(const Vector3D &a, const Vector3D &b);

	/**
	 * @brief Cross product operation.
	 * @param v The second operand of this operation.
	 * @return The cross product of this vector and v.
	 */
	inline Vector3D	cross(const Vector3D &v) const;

	/**
	 * @brief Get the length of this vector.
	 * @return The length of this vector.
	 */
	inline Real getLength() const;

	/**
	 * @brief Get the unit vector with the same direction of this vector.
	 *        Note in this procedure, nothing will be done to handle the 
	 *		  setZero-division problem.
	 * @return The unit vector.
	 */
	inline Vector3D getUnit() const;

	/**
	 * @brief Set the vector to a unit length.
	 *        Note in this procedure, nothing will be done to handle the 
	 *		  setZero-division problem.
	 */
	inline void setUnit();

	/**
	 * @brief Set this vector to setZero.
	 *        Note in this procedure, nothing will be done to handle the 
	 *		  setZero-division problem.
	 */
	inline void setZero();

	/**
	 * @brief Default constructor, set both elements to setZero.
	 */
	Vector3D() : x(0.0), y(0.0), z(0.0) { }

	/**
	 * @brief Constructor with specified two elements.
	 * @param paramX The first element value.
	 * @param paramY The second element value.
	 */
	Vector3D(Real paramX, Real paramY, Real paramZ) : x(paramX), y(paramY), z(paramZ) { }

} *Vector3DPtr, *Vector3DArray;

/**
 * @typedef Vector3DPtr
 * @brief Pointer to Vector3D representing reference to a single Vector3D object.
 */

/**
 * @typedef Vector3DArray
 * @brief Pointer to Vector3D representing an array of Vector3D object.
 */

/**
 * @brief Compute the distance from a point to a plane determined by a
 *        triangle.
 * @param p The point. Its distance to the plane will be computed.
 * @param tri1 The first point to determine the plane.
 * @param tri2 The second point to determine the plane.
 * @param tri3 The third point to determine the plane.
 * @return The distance.
 */
Real distance(const Vector3D &p, const Vector3D &tri1,
	 const Vector3D &tri2, const Vector3D &tri3);

 
Vector3D Vector3D::operator+(const Vector3D &v) const
{
	return Vector3D(x + v.x, y + v.y, z + v.z);
}

Vector3D Vector3D::operator-(const Vector3D &v) const
{
	return Vector3D(x - v.x, y - v.y, z - v.z);
}

Vector3D Vector3D::operator*(Real v) const
{
	return Vector3D(x * v, y * v, z * v);
}

Vector3D Vector3D::operator/(Real v) const
{
	return Vector3D(x / v, y / v, z / v);
}

Real Vector3D::operator*(const Vector3D &v) const
{
	return x * v.x + y * v.y + z * v.z;
}

Vector3D& Vector3D::operator+=(const Vector3D &v)
{
	x += v.x;
	y += v.y;
	z += v.z;
	return *this;
}

Vector3D& Vector3D::operator-=(const Vector3D &v)
{
	x -= v.x;
	y -= v.y;
	z -= v.z;
	return *this;
}

Vector3D Vector3D::cross(const Vector3D &a, const Vector3D &b)
{
	return Vector3D(
		a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x);
}

Vector3D Vector3D::cross(const Vector3D &v) const
{
	return Vector3D(
		y * v.z - z * v.y,
		z * v.x - x * v.z,
		x * v.y - y * v.x);
}

Real Vector3D::getLength() const
{
	return sqrt(x * x + y * y + z * z);
}

Vector3D Vector3D::getUnit() const
{
	double r = sqrt(x * x + y * y + z * z);
	return Vector3D(x / r, y / r, z / r);
}

void Vector3D::setUnit()
{
	double r = sqrt(x * x + y * y + z * z);
	x /= r;
	y /= r;
	z /= r;
}

void Vector3D::setZero()
{
	x = 0.0;
	y = 0.0;
	z = 0.0;
}

}	// Geometry

}	// Kernel

}	// Klein

#endif
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

#ifndef KLEIN_KERNEL_MATHS_MATRIX3X3_INCLUDE
#define KLEIN_KERNEL_MATHS_MATRIX3X3_INCLUDE

#include "real_utility.h"
#include "vector3d.h"
#include "memory_utility.h"

namespace Klein
{

namespace Kernel
{

namespace Maths
{

using namespace Klein::Kernel::Geometry;
using namespace Klein::Kernel::Utility;

/**
* @brief 3x3-dimension matrix.
*/
typedef class Matrix3x3
{
public:
	Real data_[9];	/**< The data, stored column by column. That means the storage order is 
					(0, 0), (1, 0), (2, 0), ... */

public:
	/**
	 * @brief Retrieve one column of the matrix.
	 * @param ID of the column to get.
	 * @return A pointer to the Real data in the column.
	 */
	//@{
	inline Real* operator[](int column);
	inline const Real* operator[](int column) const;
	//@}

	/**
	* @brief Matrix3x3 assignment.
	* @param m The original matrix.
	* @return The reference of the assigned matrix.
	*/
	inline Matrix3x3& operator=(const Matrix3x3 &m);

	/**
	* @brief Matrix3x3 element-wise addition.
	* @param v The right operand of the matrix addition.
	* @return The result of the matrix.
	*/
	inline Matrix3x3 operator+(const Matrix3x3 &v) const;

	/**
	* @brief Matrix3x3 element-wise subtraction.
	* @param v The right operand of the matrix addition.
	* @return The result of the matrix.
	*/
	inline Matrix3x3 operator-(const Matrix3x3 &v) const;

	/**
	* @brief Matrix3x3 multiplication.
	* @param v The right operand of the matrix multiplication.
	* @return The result of the matrix.
	*/
	inline Matrix3x3 operator*(const Matrix3x3 &v) const;

	/**
	* @brief A Matrix4x4 multiply a Vector4D.
	* @param v The right operand of the multiplication.
	* @return The result of the multiplication.
	*/
	inline Vector3D operator*(const Vector3D &v) const;

	/**
	* @brief Element-wise data addition.
	* @param res The matrix destination.
	* @param m1 The first source matrix.
	* @param m2 The second source matrix.
	*/
	inline static void add(Matrix3x3 &res, const Matrix3x3 &m1, const Matrix3x3 &m2);

	/**
	* @brief Element-wise data subtraction.
	* @param res The matrix destination.
	* @param m1 The first source matrix.
	* @param m2 The second source matrix.
	*/
	inline static void substract(Matrix3x3 &res, const Matrix3x3 &m1, const Matrix3x3 &m2);

	/**
	* @brief Matrix3x3 multiplication.
	* @param res The matrix destination.
	* @param m1 The first source matrix.
	* @param m2 The second source matrix.
	*/
	inline static void multiply(Matrix3x3 &res, const Matrix3x3 &m1, const Matrix3x3 &m2);

	/**
	* @brief A Matrix3x3 multiply a Vector3D.
	* @param res The result of the multiplication.
	* @param m1 The source matrix.
	* @param m2 The source vector.
	*/
	inline static void multiply(Vector3D &res, const Matrix3x3 &m1, const Vector3D &m2);

	/**
	* @brief Matrix3x3 inverse.
	* @param res The matrix destination.
	* @param src The source matrix. Make sure the source matrix is invertible.
	*/
	inline static bool inverse(Matrix3x3 &res, const Matrix3x3 &src);

	/**
	* @brief Matrix3x3 inverse.
	* @param res The matrix destination.
	*/
	inline bool inverse(Matrix3x3 &res);
	
	/**
	* @brief Matrix3x3 determinant.
	* @param m The matrix source.
	* @return The determinant of the matrix3x3.
	*/
	inline static Real determinant(const Matrix3x3 &m);

	/**
	* @brief Matrix3x3 determinant.
	* @return The determinant of the matrix3x3.
	*/
	inline Real determinant() const;

	/**
	* @brief Check whether the Matrix3x3 is invertible.
	* @param m The matrix source.
	* @param tolerance The tolerance used to check whether the determinant
	*		 is equal to zero.
	* @return true if the Matrix3x3 is invertible.
	*/
	inline static bool isInvertible(const Matrix3x3 &m, Real tolerance);

	/**
	* @brief Check whether the Matrix3x3 is invertible.
	* @param tolerance The tolerance used to check whether the determinant
	*		 is equal to zero.
	* @return true if the Matrix3x3 is invertible.
	*/
	inline bool isInvertible(Real tolerance) const;

	/**
	* @brief Element-wise data copy from the specified source matrix to the matrix.
	* @param dest The matrix destination.
	* @param src The source matrix.
	*/
	inline static void copy(Matrix3x3 &dest, const Matrix3x3 &src);

	/**
	* @brief Element-wise copy the specified data.
	*		  Note this function will not check the length of the specified data.
	* @return The pointer to the data memory block.
	*/
	inline void copyData(const Real *data);

	/**
	* @brief Copy the specified matrix to the local one.
	*		  Note if the size is not the same, it will copy the minimum subset.
	*/
	inline void setMatrix(const Matrix3x3 &m);

	/**
	* @brief Set all data to zero.
	*/
	inline void setZero();

	/**
	* @brief Set the matrix to an diagonal matrix, with the values in the diagonal position.
	* @param value The diagonal value of each element.
	*/
	inline void setDiagonal(Real value);

	/**
	* @brief Set the matrix to an diagonal matrix, with the values in the diagonal position.
	* @param values The diagonal values of each element.
	* @param count The count of the specified diagonal values.
	*/
	inline void setDiagonal(Real *values, int count);

	/**
	* @brief Set the matrix to an identity matrix.
	*/
	inline void setIdentity();

	/**
	* @brief Create a 3x3-identity matrix.
	* @return The created 3x3 matrix.
	*/
	inline static Matrix3x3 identityMatrix();

	/**
	* @brief Solve a linear Equation mx=v
	*		  Note that the matrix m is invertible.
	* @param m The specified coefficient of the linear equation
	* @param v The right coefficient of the linear equation
	* @param res The result of the linear Equation.
	* @return true if the matrix is invertible
	*/
	inline static bool linearEquationSolver(const Matrix3x3 &m, const Vector3D &v, Vector3D &res);

	/**
	* @brief Solve a linear Equation mx=v
	*		  Note that the matrix m is invertible.
	* @param v The right coefficient of the linear equation
	* @param res The result of the linear Equation.
	* @return true if the matrix is invertible
	*/
	inline bool linearEquationSolver(const Vector3D &v, Vector3D &res);
public:

	/**
	* @brief Constructor with specified matrix.
	* @param m The specified matrix.
	*/
	inline Matrix3x3(const Matrix3x3 &m);

	/**
	* @brief Constructor with the dimensions.
	* @param row The initial row count.
	* @param column The initial column count.
	*/
	inline Matrix3x3();

	inline Matrix3x3(const Real m00, const Real m01, const Real m02,
			  const Real m10, const Real m11, const Real m12,
			  const Real m20, const Real m21, const Real m22) {
		data_[0] = m00;
		data_[1] = m01;
		data_[2] = m02;
		data_[3] = m10;
		data_[4] = m11;
		data_[5] = m12;
		data_[6] = m20;
		data_[7] = m21;
		data_[8] = m22;
	}

	/**
	* @brief Destructor.
	*/
	inline ~Matrix3x3();

private:
	/**
	* @brief Calculate an algebraic complement of a matrix3x3
	* @param row The row number of the element
	* @param column The column number of the element
	* @param m The specified matrix.
	* @return The algebraic complement.
	*/
	inline static Real algebraicComplement(int row, int column, const Matrix3x3 &m);

} *Matrix3x3Ptr, *Matrix3x3Array;

/**
* @typedef Matrix3x3Ptr
* @brief Pointer to Matrix3x3 representing reference to a single Matrix3x3 object.
*/

/**
* @typedef Matrix3x3Array
* @brief Pointer to Matrix3x3 representing an array of Matrix3x3 object.
*/

inline Real* Matrix3x3::operator[](int column)
{
	return data_ + column * 3;
}

inline const Real* Matrix3x3::operator[](int column) const
{
	return data_ + column * 3;
}

/**
* @brief Calculate the eigen value and vector of a matrix3x3
* @param input the matrix3x3 to be caculated
* @param v the matrix of eigen vector
* @param d the diag matrix of eigen value
* @return the rank of eigen vector
*/
inline int eigenValue(const Matrix3x3 &input, Matrix3x3 &v, Matrix3x3 &d);




const double TOL = 1e-6;
inline Matrix3x3::Matrix3x3()
{
	setZero();
}

inline Matrix3x3::~Matrix3x3()
{
}

inline void Matrix3x3::add(Matrix3x3 &res, const Matrix3x3 &m1, const Matrix3x3 &m2)
{
	Real *ps = res.data_;
	for (int i = 0; i < 9; i++)
		*(ps++) = m1.data_[i] + m2.data_[i];
}

inline void Matrix3x3::substract(Matrix3x3 &res, const Matrix3x3 &m1, const Matrix3x3 &m2)
{
	Real *ps = res.data_;
	for (int i = 0; i < 9; i++)
		*(ps++) = m1.data_[i] - m2.data_[i];
}

inline void Matrix3x3::setZero()
{
	MemoryUtility::setZero(data_, 9 * sizeof(Real));
}

inline void Matrix3x3::setDiagonal(Real value)
{
	MemoryUtility::setZero(data_, 9 * sizeof(Real));
	for (int i = 2; i >= 0; i--)
		data_[i * 3 + i] = value;
}

inline void Matrix3x3::setIdentity()
{
	setDiagonal(1.0);
}

inline void Matrix3x3::copy(Matrix3x3 &dest, const Matrix3x3 &src)
{
	int size = 9 * sizeof(Real);
	MemoryUtility::copy(dest.data_, size, src.data_, size);
}

inline void Matrix3x3::multiply(Matrix3x3 &res, const Matrix3x3 &m1, const Matrix3x3 &m2)
{
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
		{
			Real sum = 0.0;
			for (int k = 0; k < 3; k++)
				sum += m1.data_[i + 3 * k] * m2.data_[k + 3 * j];
			res.data_[j * 3 + i] = sum;
		}
}

inline void Matrix3x3::multiply(Vector3D &res, const Matrix3x3 &m1, const Vector3D &m2)
{
	res.x = m1.data_[0] * m2.x + m1.data_[3] * m2.y + m1.data_[6] * m2.z;
	res.y = m1.data_[1] * m2.x + m1.data_[4] * m2.y + m1.data_[7] * m2.z;
	res.z = m1.data_[2] * m2.x + m1.data_[5] * m2.y + m1.data_[8] * m2.z;
}

inline Matrix3x3& Matrix3x3::operator=(const Matrix3x3 &m)
{
	copy(*this, m);
	return *this;
}

inline Matrix3x3 Matrix3x3::operator+(const Matrix3x3 &v) const
{
	Matrix3x3 retval;
	add(retval, *this, v);
	return retval;
}

inline Matrix3x3 Matrix3x3::operator-(const Matrix3x3 &v) const
{
	Matrix3x3 retval;
	substract(retval, *this, v);
	return retval;
}

inline Matrix3x3 Matrix3x3::operator*(const Matrix3x3 &v) const
{
	Matrix3x3 retval;
	multiply(retval, *this, v);
	return retval;
}

inline Vector3D Matrix3x3::operator*(const Vector3D &v) const
{
	Vector3D retval;
	multiply(retval, *this, v);
	return retval;
}

inline Matrix3x3::Matrix3x3(const Matrix3x3 &m)
{
	int size = 9 * sizeof(Real);
	MemoryUtility::copy(data_, size, m.data_, size);
}

inline void Matrix3x3::setMatrix(const Matrix3x3 &m)
{
	int size = 9 * sizeof(Real);
	MemoryUtility::copy(data_, size, m.data_, size);
}

inline bool Matrix3x3::inverse(Matrix3x3 &res, const Matrix3x3 &src)
{
	Real der = src.determinant();
	if (RealUtility::isZero(der, 1E-12))
		return false;
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			res.data_[i + j * 3] = (algebraicComplement(j, i, src) / der);
	return true;
}

inline bool Matrix3x3::inverse(Matrix3x3 &res)
{
	return inverse(res, *this);
}

inline Real Matrix3x3::determinant(const Matrix3x3 &m)
{
	Real retval = 0.0;
	for (int i = 0; i < 3; ++i)
		retval += m.data_[3 * i] * algebraicComplement(0, i, m);
	return retval;
}

inline Real Matrix3x3::determinant() const
{
	return determinant(*this);
}

inline bool Matrix3x3::isInvertible(const Matrix3x3 &m, Real tolerance)
{
	return RealUtility::isZero(m.determinant(), tolerance);
}

inline bool Matrix3x3::isInvertible(Real tolerance) const
{
	return RealUtility::isZero(determinant(), tolerance);
}

inline Matrix3x3 Matrix3x3::identityMatrix()
{
	Matrix3x3 m;
	m.setDiagonal(1.0);
	return m;
}

inline void Matrix3x3::copyData(const Real *data)
{
	int size = 9 * sizeof(Real);
	MemoryUtility::copy(data_, size, data, size);
}

inline Real Matrix3x3::algebraicComplement(int row, int column, const Matrix3x3 &m)
{
	Real t[4];
	int k = 0;
	for (int j = 0; j < 3; ++j)
		for (int i = 0; i < 3; ++i)
			if (i != row && j != column)
				t[k++] = m.data_[i + j * 3];
	Real retval = t[0] * t[3] - t[1] * t[2];
	if ((row + column) % 2 != 0)
		retval = -retval;
	return retval;
}

inline bool Matrix3x3::linearEquationSolver(const Matrix3x3 &m, const Vector3D &v, Vector3D &res)
{
	Matrix3x3 inverseMat;
	if (!inverse(inverseMat, m))
		return false;
	res =  inverseMat * v;
	return true;
}

inline bool Matrix3x3::linearEquationSolver(const Vector3D &v, Vector3D &res)
{
	return linearEquationSolver(*this, v, res);
}

inline static int getEigenVector(const Matrix3x3 &input, Real eig, Real* result)
{
	Matrix3x3 K = input;
	for (int i = 0; i < 3; i ++)
		K.data_[i * 4] -= eig;

	Vector3D bases[3];
	for (int i = 0; i < 3; i++)
		bases[i] = Vector3D(K.data_[i * 3], K.data_[i * 3 + 1], K.data_[i * 3 + 2]);

	Vector3D sum(0.0, 0.0, 0.0);
	int sumc = 0;
	for (int i = 0; i < 3; i++)
	{
		Vector3D n = Vector3D::cross(bases[i], bases[(i + 1) % 3]);
		Real nlen = n.getLength();
		if (nlen > 1e-7)
		{
			n = n / nlen;
			sum += n;
			if (sum.getLength() <= 1e-7) {
				sum = n;
				sumc = 1;
			} else {
				sumc++;
			}
		}
	}
	if (sumc == 0)
	{
		bool found = false;
		for (int i = 0; i < 3; i++)
		{
			Real basel = bases[i].getLength();
			if (basel > 1e-7)
			{
				sum = bases[i] / basel;
				found = true;
				break;
			}
		}
		if (!found)
			sum = Vector3D(1, 0, 0);
	}
	else
	{
		sum = sum.getUnit();
	}

	result[0] = sum.x;
	result[1] = sum.y;
	result[2] = sum.z;

	return 3 - sumc;

	/*
	Real _max = K.data_[0];
	int row = 0;
	for (int i = 0; i < 3; i ++)
	{
		if (abs(K.data_[i]) > abs(_max))
		{
			_max = abs(K.data_[i]);
			row = i;
		}
	}
	
	if (_max == 0)
	{
		result[0] = 1; result[1] = result[2] = 0;
		return 2;
	}
	else
	{
		Real temp[3] = {K.data_[0], K.data_[3], K.data_[6]};
		for (int i = 0; i < 3; i ++)
		{
			K.data_[i * 3] = K.data_[i * 3 + row];
			K.data_[i * 3 + row] = temp[i];

		}
		for (int i = 0; i < 3; i ++)
			K.data_[i * 3] = K.data_[i * 3] / _max;
		for (int i = 1; i < 3; i ++)
		{
			Real temp = K.data_[i];
			for (int j = 0; j < 3; j ++)
				K.data_[i + j * 3] = K.data_[i + j * 3] - K.data_[j * 3] * temp;
		}

		if (abs(K.data_[8] > abs(K.data_[4]))
		|| (K.data_[4] == 0 && K.data_[7] == 0))
		{
			for (int i = 1; i < 3; i ++)
			{
				Real temp = K.data_[1 + i * 3];
				K.data_[1 + i * 3] = K.data_[2 + i * 3];
				K.data_[2 + i * 3] = temp;
			}
		}
		if (K.data_[4] != 0)
		{
			Real temp = K.data_[4];
			for(int i = 1; i < 3; i ++)
				K.data_[1 + i * 3] = K.data_[1 + i * 3] / temp;
			K.data_[8] = K.data_[8] - K.data_[7] * K.data_[5];
			K.data_[5] = 0;
		}
		if (K.data_[4]  == 0.0)
		{
			result[2] = 1.0;
			result[1] = 1.0;
		}
		else
		{
			result[2] = 1.0;
			result[1] = -1.0 * K.data_[7];
		}
		result[0] = -(result[1] * K.data_[3] + result[2] * K.data_[6]);
		Real mol = 0;
		for (int i = 0; i < 3; i ++)
			mol += pow(result[i], 2.0);
		mol = sqrt(mol);
		for (int i = 0; i < 3; i ++)
			result[i] /= mol;
		if (K.data_[4] == 0.0)
			return 1;
		else
			return 0;	
	}
	*/
}

inline int eigenValue(const Matrix3x3 &input, Matrix3x3 &v, Matrix3x3 &d)
{
	Real m = (input.data_[0] + input.data_[4] + input.data_[8])/3.0;
	Matrix3x3 K = input;
	for (int i = 0; i < 3; i ++)
		K.data_[i * 4] = K.data_[i * 4] - m;
	Real q = K.determinant() / 2;
	Real p = 0.0;

	for (int i = 0; i < 3; i ++)
		for (int j = 0; j < 3; j ++)
		{
			p += pow(K.data_[i + j * 3],2.0);
		}
	p /= 6.0;

	Real acosv = q / pow(p, 1.5);
	if (q == 0) acosv = 0.0;
	if (acosv > 1.0) acosv = 1.0;
	if (acosv < -1.0) acosv = -1.0;
	Real phi = 1.0 / 3.0 * acos(acosv);
	if (phi < 0)
		phi += PI / 3.0;

	Real eig1 = m + 2 * sqrt(p) * cos(phi);
	Real eig2 = m - sqrt(p) * (cos(phi) + sqrt(3.0) * sin(phi));
	Real eig3 = m - sqrt(p) * (cos(phi) - sqrt(3.0) * sin(phi));

	Real eigs[3] = {eig1,eig2,eig3};
	for (int i = 0; i < 3; i ++)
		for (int j = i; j < 3; j ++)
		{
			if (eigs[i] > eigs[j])
			{
				Real temp = eigs[i];
				eigs[i] = eigs[j];
				eigs[j] = temp;
			}
	
		}
	d.setZero();
	for (int i = 0; i < 3; i ++)
		d.data_[i * 4] = eigs[i];

	Real temp[3];

	if (abs(eigs[0] - eigs[2]) < TOL)
	{
		v.setIdentity();
		return 2;
	}
	
	else if (abs(eigs[0] - eigs[1]) < TOL)
	{
		getEigenVector(input, eigs[0], temp);
		Vector3D l[3];
		l[1] = Vector3D(temp[0], temp[1], temp[2]);
		getEigenVector(input, eigs[2], temp);
		l[2] = Vector3D(temp[0], temp[1], temp[2]);
		l[0] = l[2].cross(l[1]);
		for (int i = 0; i < 3; i ++)
		{
				v.data_[i * 3] = l[i].x;
				v.data_[i * 3 + 1] = l[i].y;
				v.data_[i * 3 + 2] = l[i].z;
		}
		return 1;
	}
	else if (abs(eigs[1] - eigs[2]) < TOL)
	{
		getEigenVector(input, eigs[0], temp);
		Vector3D l[3];
		l[0] = Vector3D(temp[0], temp[1], temp[2]);
		getEigenVector(input, eigs[2], temp);
		l[2] = Vector3D(temp[0], temp[1], temp[2]);
		l[1] = (l[2].cross(l[0]))*-1;
		for (int i = 0; i < 3; i ++)
		{
				v.data_[i * 3] = l[i].x;
				v.data_[i * 3 + 1] = l[i].y;
				v.data_[i * 3 + 2] = l[i].z;
		}
		return 1;
	}
	else
	{
		for (int i = 0; i < 3; i ++)
		{
			getEigenVector(input, eigs[i], temp);
			for (int j = 0; j < 3; j ++)
				v.data_[j + i * 3] = temp[j];
		}
		if (v.determinant() < 0)
		{
			v.data_[0] *= -1;v.data_[1] *= -1;v.data_[2] *= -1;
		}
		return 0;
	}
}







}	// Maths

}	// Kernel

}	// Klein

#endif // KLEIN_KERNEL_MATHS_MATRIX3X3_INCLUDE
/* ================================================================
 *
 * PyCA Project
 *
 * Copyright (c) J. Samuel Preston, Linh K. Ha, Sarang C. Joshi. All
 * rights reserved.  See Copyright.txt or for details.
 *
 * This software is distributed WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the above copyright notice for more information.
 *
 * ================================================================ */

#include "Mat3D.h"

//#include <cstring>
#include <iostream>
//#include <algorithm>

namespace PyCA {

// note: while we have float and double specialization, we only use
// double version.  Internally float arrays are copied to
// fortran-order double arrays for use in lapack functions

// ################################################################
// Mat3D
// ################################################################


template <typename T>
Mat3D<T>
::Mat3D() 
{
  eye();
}

template <typename T>
inline T&
Mat3D<T>
::operator()(unsigned int elementIndex)
{
  return a[elementIndex];
}

template <typename T>
inline const T&
Mat3D<T>
::operator()(unsigned int elementIndex) const
{
  return a[elementIndex];
}

template <typename T>
inline T&
Mat3D<T>
::operator()(unsigned int rowIndex,
	     unsigned int columnIndex)
{
  return a[rowIndex * 3 + columnIndex];
}

template <typename T>
inline const T&
Mat3D<T>
::operator()(unsigned int rowIndex,
	     unsigned int columnIndex) const
{
  return a[rowIndex * 3 + columnIndex];
}

template <typename T>
inline const T&
Mat3D<T>::
get(const unsigned int rowIndex,
      const unsigned int columnIndex)
{
  if(rowIndex > 2 || columnIndex > 2){
    throw PyCAException(__FILE__,__LINE__,"Error, index out of bounds");
  }
  return a[rowIndex * 3 + columnIndex];
}

template <typename T>
void
Mat3D<T>::
set(const unsigned int rowIndex,
    const unsigned int columnIndex,
    const T& val)
{
  if(rowIndex > 2 || columnIndex > 2){
    throw PyCAException(__FILE__,__LINE__,"Error, index out of bounds");
  }
  a[rowIndex * 3 + columnIndex] = val;
}


template<class T> 
bool 
Mat3D<T>::
operator==( const Mat3D<T> &other ) const
{
  return (this->a[0] == other.a[0] &&
	  this->a[1] == other.a[1] &&
	  this->a[2] == other.a[2] &&
	  this->a[3] == other.a[3] &&
	  this->a[4] == other.a[4] &&
	  this->a[5] == other.a[5] &&
	  this->a[6] == other.a[6] &&
	  this->a[7] == other.a[7] &&
	  this->a[8] == other.a[8]);
}

template <typename T>
void 
Mat3D<T>
::eye()
{
  a[0] = a[4] = a[8] = 1;
  a[1] = a[2] = a[3] = a[5] = a[6] = a[7] = 0;
}

template <typename T>
Mat3D<T>
Mat3D<T>:: 
Transpose()
{
  Mat3D<T> transpose;
  _transposeMatrixData(this->a, transpose.a);
  return transpose;
}

template <typename T>
std::ostream& 
Mat3D<T>
::writeASCII(std::ostream& output) const
{
  output 
    << a[0] << " " << a[1] << " " << a[2] << '\n'  // Don't want to use endl
    << a[3] << " " << a[4] << " " << a[5] << '\n'  // because that flushes the
    << a[6] << " " << a[7] << " " << a[8];         // buffer.
  return output;
}

template <typename T>
std::ostream& 
operator<<(std::ostream& output, const Mat3D<T>& matrix)
{
  return matrix.writeASCII(output);
}

//
// compute the inverse of the matrix a, store the result in ainv.
// true is returned if the inverse computation was successful, false
// if there was a problem (e.g. a is singular).  a and ainv may point
// to the same memory, both must be initialized.  If false is
// returned, ainv is unchanged.
//
// uses the clapack functions dgetrf_ and dgetri_ to compute the LU
// decomposition and then the inverse.
//
// bcd 2004
//
template <typename T>
bool
Mat3D<T>
::computeInverse(const T* const a, T* const ainv)
{
  // manual computation of inverse, specific to 3x3 -JDH 2013
  ainv[0] = (a[4]*a[8]-a[5]*a[7]);
  ainv[3] = -(a[3]*a[8]-a[5]*a[6]); 
  ainv[6] = (a[3]*a[7]-a[4]*a[6]);
  ainv[1] = -(a[1]*a[8]-a[2]*a[7]);
  ainv[4] = (a[0]*a[8]-a[2]*a[6]);
  ainv[7] = -(a[0]*a[7]-a[1]*a[6]);
  ainv[2] = (a[1]*a[5]-a[2]*a[4]);
  ainv[5] = -(a[0]*a[5]-a[2]*a[3]);
  ainv[8] = (a[0]*a[4]-a[1]*a[3]);

  // determinant, using precomputed pieces
  T deta = a[0]*ainv[0] + a[1]*ainv[3] + a[2]*ainv[6];

  // divide by determinant
  for (size_t i=0;i < 9;++i)
      ainv[i] /= deta;

  return true;
}

// Compute determinant directly
template <typename T>
T
Mat3D<T>
::det()
{
  return a[0]*(a[4]*a[8]-a[5]*a[7]) 
       - a[1]*(a[3]*a[8]-a[5]*a[6]) 
       + a[2]*(a[3]*a[7]-a[4]*a[6]);
}

// If matrix is singular, returns false and leaves matrix unchanged.
template <typename T>
bool
Mat3D<T>
::invert()
{
  return Mat3D<T>::computeInverse(this->a, this->a);
}

template <typename T, typename U>
Vec3D<T> 
operator*(const Mat3D<T>& m, const Vec3D<U>& v)
{
   return Vec3D<T>(v.x * m.a[0] +
		   v.y * m.a[1] +
		   v.z * m.a[2],
		   v.x * m.a[3] +
		   v.y * m.a[4] +
		   v.z * m.a[5],
		   v.x * m.a[6] +
		   v.y * m.a[7] +
		   v.z * m.a[8]);
}

template <typename T, typename U>
Mat3D<T>& 
operator*=(Mat3D<T>& lhs, const Mat3D<U>& rhs)
{
  double tmp1, tmp2;
  tmp1     = lhs.a[0] * rhs.a[0] + lhs.a[1] * rhs.a[3] + lhs.a[2] * rhs.a[6];
  tmp2     = lhs.a[0] * rhs.a[1] + lhs.a[1] * rhs.a[4] + lhs.a[2] * rhs.a[7];
  lhs.a[2] = lhs.a[0] * rhs.a[2] + lhs.a[1] * rhs.a[5] + lhs.a[2] * rhs.a[8];
  lhs.a[0] = tmp1;
  lhs.a[1] = tmp2;

  tmp1     = lhs.a[3] * rhs.a[0] + lhs.a[4] * rhs.a[3] + lhs.a[5] * rhs.a[6];
  tmp2     = lhs.a[3] * rhs.a[1] + lhs.a[4] * rhs.a[4] + lhs.a[5] * rhs.a[7];
  lhs.a[5] = lhs.a[3] * rhs.a[2] + lhs.a[4] * rhs.a[5] + lhs.a[5] * rhs.a[8];
  lhs.a[3] = tmp1;
  lhs.a[4] = tmp2;

  tmp1     = lhs.a[6] * rhs.a[0] + lhs.a[7] * rhs.a[3] + lhs.a[8] * rhs.a[6];
  tmp2     = lhs.a[6] * rhs.a[1] + lhs.a[7] * rhs.a[4] + lhs.a[8] * rhs.a[7];
  lhs.a[8] = lhs.a[6] * rhs.a[2] + lhs.a[7] * rhs.a[5] + lhs.a[8] * rhs.a[8];
  lhs.a[6] = tmp1;
  lhs.a[7] = tmp2;

  return lhs;
}

template <typename T, typename U>
Mat3D<T>
operator*(const Mat3D<T>& lhs, const Mat3D<U>& rhs)
{
  Mat3D<T> tmp(lhs);
  tmp *= rhs;
  return tmp;
}

// template instantiation 

template class Mat3D<double>;
template class Mat3D<float>;

template Vec3D<float>
operator*(const Mat3D<float>& m, const Vec3D<float>& v);
template Vec3D<float>
operator*(const Mat3D<float>& m, const Vec3D<double>& v);
template Vec3D<double>
operator*(const Mat3D<double>& m, const Vec3D<float>& v);
template Vec3D<double>
operator*(const Mat3D<double>& m, const Vec3D<double>& v);

template Mat3D<float>&
operator*=(Mat3D<float>& m, const Mat3D<float>& v);
template Mat3D<float>&
operator*=(Mat3D<float>& m, const Mat3D<double>& v);
template Mat3D<double>&
operator*=(Mat3D<double>& m, const Mat3D<float>& v);
template Mat3D<double>&
operator*=(Mat3D<double>& m, const Mat3D<double>& v);

template Mat3D<float>
operator*(const Mat3D<float>& m, const Mat3D<float>& v);
template Mat3D<float>
operator*(const Mat3D<float>& m, const Mat3D<double>& v);
template Mat3D<double>
operator*(const Mat3D<double>& m, const Mat3D<float>& v);
template Mat3D<double>
operator*(const Mat3D<double>& m, const Mat3D<double>& v);

} // end namespace PyCA

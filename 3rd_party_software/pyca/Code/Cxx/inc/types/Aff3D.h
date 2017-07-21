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

// Mark Foskey 5/04

#ifndef AFF3D_H
#define AFF3D_H

#ifndef SWIG

#include <stdexcept>

#include "Mat3D.h"
#include "Vec3D.h"

#endif // !SWIG

namespace PyCA {

template< class T >
class Aff3D
{
public:

  typedef T CoordinateType;
  typedef Aff3D< CoordinateType > Self;
  typedef Mat3D< CoordinateType > MatrixType;
  typedef Vec3D< CoordinateType > VectorType; 

  MatrixType matrix;
  VectorType vector;

  Aff3D() : matrix(), vector() {} // Identity.
  Aff3D(const MatrixType& m, const VectorType& v ) 
    : matrix(m), vector(v) {}
  explicit Aff3D( const MatrixType& m ) : matrix(m), vector() {}
  explicit Aff3D( const VectorType& v ) : matrix(), vector(v) {}

  template< class U >
  Aff3D( const Aff3D<U>& t )
    : matrix( t.matrix ), vector( t.vector ) {}

#ifndef SWIG
  // Do not wrap assignment operator, python doesn't support it
  template< class U >
  Self& operator=( const Aff3D<U>& rhs )
  { 
    this->matrix = rhs.matrix; 
    this->vector = rhs.vector; 
    return *this; 
  }
#endif // !SWIG

  void eye() { *this = Aff3D<T>(); } // Sets to identity
  bool invert();
  Self& operator*=( const Self& rhs );
    
  Self operator*( const Self& other ) const;

template< class U >
  void ApplyTransform(const Aff3D<U>& rhs)
  {
	this->matrix = this->matrix * rhs.matrix;
	this->vector = this->vector + rhs.vector;
  }

  VectorType operator*( const VectorType& v ) const;
  void transformVector( const VectorType& vIn, VectorType& vOut ) const;
  void transformVector( VectorType& v ) const;

  template <class U, class V>
  void transformCoordinates( const U& xIn, const U& yIn, const U& zIn, 
                             V& xOut, V& yOut, V& zOut ) const
  {
    xOut =(V)(xIn*matrix.a[0] + yIn*matrix.a[1] + zIn*matrix.a[2] + vector[0]);
    yOut =(V)(xIn*matrix.a[3] + yIn*matrix.a[4] + zIn*matrix.a[5] + vector[1]);
    zOut =(V)(xIn*matrix.a[6] + yIn*matrix.a[7] + zIn*matrix.a[8] + vector[2]);
  }

  bool operator==( const Self t ) const;

  // No other operators; e.g., '+' isn't really a sensible op.

  void writePLUNCStyle(const char* filename) const;
  void writePLUNCStyle(const std::string& filename) const
    { writePLUNCStyle(filename.c_str()); }
  void readPLUNCStyle(const char* filename);
  void readPLUNCStyle(const std::string& filename)
    { readPLUNCStyle(filename.c_str()); }

double determinant()
{
double a = this->matrix.a[0];
double b = this->matrix.a[1];
double c = this->matrix.a[2];
double d = this->matrix.a[3];
double e = this->matrix.a[4];
double f = this->matrix.a[5];
double g = this->matrix.a[6];
double h = this->matrix.a[7];
double i = this->matrix.a[8];
double tp1 = (a*e*i) + (b*f*g) + (c*d*h);
double tp2 = (a*f*h) + (b*d*i) + (c*e*g);
double det = tp1 - tp2;
return det;
}
private:

  Aff3D( const Self& first, const Self& second );

};

template <typename T>
std::ostream& 
operator<<(std::ostream& output, const Aff3D<T>& t)
{
  return output << '\n' << t.matrix << "\n\n" << t.vector[0] << " "
                << t.vector[1] << " " << t.vector[2];
}

typedef Aff3D<float> Aff3Df;
typedef Aff3D<double> Aff3Dd;

} // end namespace PyCA

#endif  // ndef AFFINETRANSFORM_H

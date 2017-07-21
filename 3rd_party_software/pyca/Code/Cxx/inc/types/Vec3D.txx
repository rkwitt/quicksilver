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

#ifndef __VEC3D_TXX__
#define __VEC3D_TXX__

template <class T>
inline
Vec3D<T>
::Vec3D()
  : x(0), y(0), z(0)
{}

template <class T>
inline
Vec3D<T>
::Vec3D(const T& xIn, 
	   const T& yIn, 
	   const T& zIn)
  :x(xIn),y(yIn), z(zIn)
{}

template <class T>
inline
const T& 
Vec3D<T>::
get(unsigned int index) const
{
  if(index > 2){
    throw PyCAException(__FILE__,__LINE__,"Error, index out of bounds");
  }
  return (*this)[index];
}

template <class T>
inline
void
Vec3D<T>::
set(unsigned int index, const T& val)
{
  if(index > 2){
    throw PyCAException(__FILE__,__LINE__,"Error, index out of bounds");
  }
  (*this)[index] = val;
}

template <class T>
inline
T&
Vec3D<T>
::operator[](unsigned int index) 
{
  return (&x)[index];
}
  
template <class T>
inline
const T&
Vec3D<T>
::operator[](unsigned int index) const
{
  return (&x)[index];
}

template <class T>
inline
double
Vec3D<T>
::length() const
{
  return sqrt(x*x + y*y + z*z);
}

template <class T>
inline
double
Vec3D<T>
::lengthSquared() const
{
  return x*x + y*y + z*z;
}

template <class T>
inline
double
Vec3D<T>
::normL2() const
{
  return sqrt(x*x + y*y + z*z);
}

template <class T>
inline
double
Vec3D<T>
::normL1() const
{
  return (x > 0 ? x : -x) + (y > 0 ? y : -y) + (z > 0 ? z : -z);
}

template <class T>
inline
typename Vec3D<T>::acc_type
Vec3D<T>
::prod() const
{
  return x * y * z;
}

template <class T>
inline
typename Vec3D<T>::acc_type
Vec3D<T>
::sum() const
{
  return x + y + z;
}

template <class T>
inline
T
Vec3D<T>
::maxElement() const
{
  return std::max(x, std::max(y, z));
}

template <class T>
inline
T
Vec3D<T>
::minElement() const
{
  return std::min(x, std::min(y, z));
}

template <class T>
inline
void
Vec3D<T>
::set(const T& x,
      const T& y,
      const T& z)
{ 
  this->x = x;
  this->y = y;
  this->z = z;
}

template <class T>
inline
void
Vec3D<T>
::scale(const T& sx,
	const T& sy,
	const T& sz)
{
  x *= sx;
  y *= sy;
  z *= sz;
}

template <class T>
inline
void
Vec3D<T>
::translate(const T& tx,
	    const T& ty,
	    const T& tz)
{
  x += tx;
  y += ty;
  z += tz;
}

template <class T>
inline
void
Vec3D<T>
::negate()
{
  x = -x;
  y = -y;
  z = -z;
}  

template <class T>
inline
void
Vec3D<T>
::invert()
{
   throw PyCAException("Vec3D::invert not supported by this template type");
}  

// float specialization
template <>
inline
void
Vec3D<float>
::invert()
{
   x = x == 0.f ? 0.f : 1.0/x;
   y = y == 0.f ? 0.f : 1.0/y;
   z = z == 0.f ? 0.f : 1.0/z;
}  

template <class T>
inline
Vec3D<T>
Vec3D<T>::
inverse() const
{
   throw PyCAException("Vec3D::inverse not supported by this template type");
}  

// float specialization
template <>
inline
Vec3D<float>
Vec3D<float>::
inverse() const
{
   Vec3D<float> v = *this;
   v.invert();
   return v;
}  

template <class T>
inline
void
Vec3D<T>
::normalize()
{
   throw PyCAException("Vec3D::normalize not supported by this template type");
}  

// float specialization
template <>
inline
void
Vec3D<float>
::normalize()	
{
   
  double len = length();
  if (length() > 0){
     len = 1.0/len;
     x *= len;
     y *= len;
     z *= len;
  }else{
     x = 0.f;
     y = 0.f;
     z = 0.f;
  }   
}

template <class T>
inline
double
Vec3D<T>
::distance(const Vec3D<T>& rhs) const
{
  return sqrt(distanceSquared(rhs));
}

template <class T>
inline
double
Vec3D<T>
::distanceSquared(const Vec3D<T>& rhs) const
{
  return 
    (x - rhs.x) * (x - rhs.x) +
    (y - rhs.y) * (y - rhs.y) +
    (z - rhs.z) * (z - rhs.z);
}

template <class T>
inline
typename Vec3D<T>::acc_type
Vec3D<T>
::dot(const Vec3D<T>& rhs) const
{
  return (x * rhs.x + y * rhs.y + z * rhs.z);
}

template <class T>
inline
Vec3D<T>
Vec3D<T>
::cross(const Vec3D<T>& rhs) const
{
  return Vec3D<T>(y * rhs.z - z * rhs.y,
		     z * rhs.x - x * rhs.z,
		     x * rhs.y - y * rhs.x);
}

template <class T>
inline
::std::ostream& 
Vec3D<T>
::writeASCII(::std::ostream& output) const
{
  output << "(" << x << ", " << y << ", " << z << ")";
  return output;  
}

template <class T>
inline
::std::istream& 
Vec3D<T>
::readASCII(::std::istream& input)
{
          char paren, comma;
  input >> paren
	>> x >> comma
	>> y >> comma
	>> z >> paren;
  return input;
} 

template <class T>
inline
::std::ostream& 
Vec3D<T>
::writeBinary(::std::ostream& output) const
{
  output.write((char*)(&x), 3 * sizeof(T));
  return output;
}
  
template <class T>
inline
::std::istream& 
Vec3D<T>
::readBinary(::std::istream& input)
{
  input.read((char*)(&x), 3 * sizeof(T));
  return input;
}

/**
 * invert operator (unary minus)
 */
template <class T>
inline
Vec3D<T> 
operator-(const Vec3D<T>& rhs)
{
  return Vec3D<T>(-rhs.x, -rhs.y, -rhs.z);
}

/**
 * Element-wise addition
 */
template <class T, class U>
inline
Vec3D<T>& 
operator+=(Vec3D<T>& lhs, const Vec3D<U>& rhs)
{
  lhs.x += rhs.x;
  lhs.y += rhs.y;
  lhs.z += rhs.z;
  return lhs;
}

/**
 * Element-wise subtraction
 */
template <class T, class U>
inline
Vec3D<T>& 
operator-=(Vec3D<T>& lhs, const Vec3D<U>& rhs)
{
  lhs.x -= rhs.x;
  lhs.y -= rhs.y;
  lhs.z -= rhs.z;
  return lhs;
}

/**
 * Element-wise multiplication
 */
template <class T, class U>
inline
Vec3D<T>& 
operator*=(Vec3D<T>& lhs, const Vec3D<U>& rhs)
{
  lhs.x *= rhs.x;
  lhs.y *= rhs.y;
  lhs.z *= rhs.z;
  return lhs;
}

/**
 * Element-wise division
 */
template <class T, class U>
inline
Vec3D<T>& 
operator/=(Vec3D<T>& lhs, const Vec3D<U>& rhs)
{
  lhs.x /= rhs.x;
  lhs.y /= rhs.y;
  lhs.z /= rhs.z;
  return lhs;
}

/**
 * Element-wise addition
 */
template <class T, class U>
inline
const Vec3D<T> 
operator+(const Vec3D<T>& lhs, const Vec3D<U>& rhs)
{
  return Vec3D<T>(static_cast<T>(lhs.x + rhs.x), 
		     static_cast<T>(lhs.y + rhs.y), 
		     static_cast<T>(lhs.z + rhs.z));
}

/**
 * Element-wise subtraction
 */
template <class T, class U>
inline
const Vec3D<T> 
operator-(const Vec3D<T>& lhs, const Vec3D<U>& rhs)
{
  return Vec3D<T>(static_cast<T>(lhs.x - rhs.x), 
		     static_cast<T>(lhs.y - rhs.y), 
		     static_cast<T>(lhs.z - rhs.z));
}

/**
 * Element-wise multiply
 */
template <class T, class U>
inline
const Vec3D<T> 
operator*(const Vec3D<T>& lhs, const Vec3D<U>& rhs)
{
  return Vec3D<T>(static_cast<T>(lhs.x * rhs.x), 
		     static_cast<T>(lhs.y * rhs.y), 
		     static_cast<T>(lhs.z * rhs.z));
}

/**
 * Element-wise division
 */
template <class T, class U>
inline
const Vec3D<T> 
operator/(const Vec3D<T>& lhs, const Vec3D<U>& rhs)
{
 
  return Vec3D<T>(static_cast<T>(lhs.x / rhs.x), 
  		     static_cast<T>(lhs.y / rhs.y), 
		     static_cast<T>(lhs.z / rhs.z));
}

/**
 * Add constant to each element
 */
template <class T, class U>
inline
Vec3D<T>& 
operator+=(Vec3D<T>& v, const U &u)
{
  v.x += u;
  v.y += u;
  v.z += u;
  return v;
}

/**
 * Subtract constant from each element
 */
template <class T, class U>
inline
Vec3D<T>& 
operator-=(Vec3D<T>& v, const U &u)
{
  v.x -= u;
  v.y -= u;
  v.z -= u;
  return v;
}

/**
 * Scale by a constant
 */
template <class T, class U>
inline
Vec3D<T>& 
operator*=(Vec3D<T>& v, const U &u)
{
  v.x *= u;
  v.y *= u;
  v.z *= u;
  return v;
}

/**
 * Divide by a constant
 */
template <class T, class U>
inline
Vec3D<T>& 
operator/=(Vec3D<T>& v, const U &u)
{
  v.x /= u;
  v.y /= u;
  v.z /= u;
  return v;
}

/**
 * Add a constant to each element
 */
template <class T, class U>
inline
const Vec3D<T> 
operator+(const Vec3D<T>& v, const U &u)
{
  return Vec3D<T>(static_cast<T>(v.x + u), 
		     static_cast<T>(v.y + u), 
		     static_cast<T>(v.z + u)); 
}

/**
 * Add constant to each element
 */
template <class T, class U>
inline
const Vec3D<T> 
operator+(const U &u, const Vec3D<T>& v)
{
  return Vec3D<T>(static_cast<T>(v.x + u),
		     static_cast<T>(v.y + u), 
		     static_cast<T>(v.z + u)); 
}

/**
 * Subtract constant from each element
 */
template <class T, class U>
inline
const Vec3D<T> 
operator-(const Vec3D<T>& v, const U &u)
{
  return Vec3D<T>(static_cast<T>(v.x - u), 
		     static_cast<T>(v.y - u), 
		     static_cast<T>(v.z - u)); 
}

/**
 * Negate vector and add constant to each element
 */
template <class T, class U>
inline
const Vec3D<T> 
operator-(const U &u, const Vec3D<T>& v)
{
  return Vec3D<T>(static_cast<T>(u - v.x), 
		     static_cast<T>(u - v.y),
		     static_cast<T>(u - v.z)); 
}

/**
 * Scale by a constant
 */
template <class T, class U>
inline
const Vec3D<T> 
operator*(const Vec3D<T>& v, const U &u)
{
  return Vec3D<T>(static_cast<T>(v.x * u), 
		     static_cast<T>(v.y * u),
		     static_cast<T>(v.z * u)); 
}

/**
 * Scale by a constant
 */
template <class T, class U>
inline
const Vec3D<T> 
operator*(const U &u, const Vec3D<T>& v)
{
  return Vec3D<T>(static_cast<T>(v.x * u), 
		     static_cast<T>(v.y * u), 
		     static_cast<T>(v.z * u)); 
}

/**
 * Divide by a constant
 */
template <class T, class U>
inline
const Vec3D<T> 
operator/(const Vec3D<T>& v, const U &u)
{
     return Vec3D<T>(v.x / u, 
        	         v.y / u, 
	                 v.z / u); 
}

/**
 * Divide constant by a vector (element-wise)
 */
template <class T, class U>
inline
const Vec3D<T> 
operator/(const U &u, const Vec3D<T>& v)
{
  return Vec3D<T>(static_cast<T>(u / v.x), 
		     static_cast<T>(u / v.y), 
		     static_cast<T>(u / v.z)); 
}

//
// ======== input/output ========
//
template <class T>
::std::ostream& 
operator<<(::std::ostream& output, const Vec3D<T>& v)
{
  return v.writeASCII(output);
}

template <class T>
::std::istream& 
operator>>(::std::istream& input, Vec3D<T>& v)
{
  return v.readASCII(input);
}
  

template <class T>
inline const Vec3D<T> componentMax(const Vec3D<T>& a, const Vec3D<T>& b)
{
    return Vec3D<T>(::std::max(a.x, b.x),::std::max(a.y, b.y),::std::max(a.z, b.z));
}

template <class T>
inline const Vec3D<T> componentMin(const Vec3D<T>& a, const Vec3D<T>& b)
{
    return Vec3D<T>(::std::min(a.x, b.x),::std::min(a.y, b.y),::std::min(a.z, b.z));
}

#endif // __VEC3D_TXX__

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

#ifndef __VEC3D_H
#define __VEC3D_H

#ifndef SWIG
#include <iosfwd>
#include <AccumType.h>
#include <iostream>
#include <cmath>
#include "PyCAException.h"
#endif // SWIG

namespace PyCA {

//
// This struct is meant to be a fast 3D vector class.
// It is intended for use with built-in types (int, float, etc.)
// 
// Optimization is important for this struct
// - almost everything is inlined
// - readability and maintainability tradeoffs have been made 
//   for speed gains
//
// bcd 2003
//

template <class T>
class Vec3D 
{
public:
    typedef T    value_type;
    typedef typename AccumType<T>::AccT acc_type;
    
    //
    // IMPORTANT: DO NOT change the order of x, y, and z.  
    // DO NOT add other data members between x, y, and z.
    // Array style access depends on these being adjacent
    // in memory.  This is guarenteed for standards compliant
    // compilers only if they are declared as they are here.
    // DO NOT TOUCH THIS LINE
    //
    T x, y, z;
    // IMPORTANT DO NOT TOUCH THE ABOVE LINE //

    // 
    // constructors
    //

    Vec3D();
    explicit Vec3D(const T& xIn, const T& yIn, const T& zIn);
    explicit Vec3D(const T& v):x(v), y(v), z(v)   {}

    // copy constructor
    Vec3D(const Vec3D<T>& rhs):x(rhs.x), y(rhs.y), z(rhs.z){};
    

    // template constructor
    template <class U>
    Vec3D(const Vec3D<U>& rhs)
        :x((T)rhs.x), y((T)rhs.y), z((T)rhs.z)
        {}


   const T& get(unsigned int index) const;
   void set(unsigned int index, const T& val);

#ifndef SWIG
   //
   // Don't wrap these operators, they don't map to python.
   //

    // copy assigment
    template <class U>
    Vec3D<T>& operator=(const Vec3D<U>& rhs)
        {
            if (this != &rhs) {
                x = (T)rhs.x;
                y = (T)rhs.y;
                z = (T)rhs.z;
            }
            return *this;
        }

    // template assignment
    Vec3D<T>& operator=(const Vec3D<T>& rhs){
        if (this != &rhs) {
            x = rhs.x;
            y = rhs.y;
            z = rhs.z;
        }
        return *this;
    }

    
    //
    // comparison operators
    //

    bool operator==(const Vec3D<T>& v) const {
        return ((x == v.x) && (y == v.y) && (z == v.z));
    }

    bool operator!=(const Vec3D<T>& v) const{
        return !(*this == v);
    }

    bool operator<(const Vec3D<T>& rhs) const {
        if (x < rhs.x) return true;
        if (x > rhs.x) return false;
        if (y < rhs.y) return true;
        if (y > rhs.y) return false;
        return (z < rhs.z);
    }

    //
    // array style coordinate access
    //
    // with optomization turned on, this is just as fast
    // as .x, .y, .z access
    //
    T& operator[](unsigned int index);
    const T& operator[](unsigned int index) const;

#endif // SWIG

    //
    // the euclidean length of this vector
    //
    double length() const;
    double lengthSquared() const;

    //
    // L1 and L2 norms, normL2 is the same as length()
    //
    double normL1() const;
    double normL2() const;

    //
    // element accumulation
    //
    acc_type prod() const;
    acc_type sum() const;

    //
    // min and max coordinate values
    //
    T maxElement() const;
    T minElement() const;

    //
    // manipulation methods
    // 
    void set(const T& x, const T& y, const T& z);
    void scale(const T& sx, const T& sy, const T& sz);
    void translate(const T& tx, const T& ty, const T& tz);
    void negate();
    void invert();
   Vec3D<T> inverse() const;
    void normalize();

    //
    // the distance between the tip of this and the tip of rhs
    //
    double distance(const Vec3D<T>& rhs) const;
    double distanceSquared(const Vec3D<T>& rhs) const;

    //
    // inner (dot) product
    //
    acc_type dot(const Vec3D<T>& rhs) const;

    //
    // vector (cross) product
    //
    Vec3D<T> cross(const Vec3D<T>& rhs) const;

    //
    // output
    //
    std::ostream& writeASCII(std::ostream& output) const;
    std::istream& readASCII(std::istream& input);
    std::ostream& writeBinary(std::ostream& output) const;
    std::istream& readBinary(std::istream& input);

};

//
// comparison operators
//
template <class T> 
bool operator<(const Vec3D<T>& lhs, const Vec3D<T>& rhs);
template <class T>
bool operator==(const Vec3D<T>& lhs, const Vec3D<T>& rhs);
template <class T>
bool operator!=(const Vec3D<T>& lhs, const Vec3D<T>& rhs);

//
// invert operator (unary minus)
//
template <class T>
Vec3D<T> operator-(const Vec3D<T>& rhs);

//
// elementwise operators with Vec3D
//
template <class T, class U>
Vec3D<T>& operator+=(Vec3D<T>& lhs, const Vec3D<U>& rhs);
template <class T, class U>
Vec3D<T>& operator-=(Vec3D<T>& lhs, const Vec3D<U>& rhs);
template <class T, class U>
Vec3D<T>& operator*=(Vec3D<T>& lhs, const Vec3D<U>& rhs);
template <class T, class U>
Vec3D<T>& operator/=(Vec3D<T>& lhs, const Vec3D<U>& rhs);
template <class T, class U>
const Vec3D<T> operator+(const Vec3D<T>& lhs, const Vec3D<U>& rhs);
template <class T, class U>
const Vec3D<T> operator-(const Vec3D<T>& lhs, const Vec3D<U>& rhs);
template <class T, class U>
const Vec3D<T> operator*(const Vec3D<T>& lhs, const Vec3D<U>& rhs);
template <class T, class U>
const Vec3D<T> operator/(const Vec3D<T>& lhs, const Vec3D<U>& rhs);

//
// element-wise operators with a double
//
template <class T, class U>
Vec3D<T>& operator+=(Vec3D<T>& v, const U &u);
template <class T, class U>
Vec3D<T>& operator-=(Vec3D<T>& v, const U &u);
template <class T, class U>
Vec3D<T>& operator*=(Vec3D<T>& v, const U &u);
template <class T, class U>
Vec3D<T>& operator/=(Vec3D<T>& v, const U &u);

template <class T, class U>
const Vec3D<T> operator+(const Vec3D<T>& v, const U &u);
template <class T, class U>
const Vec3D<T> operator+(const U &u, const Vec3D<T>& v);
template <class T, class U>
const Vec3D<T> operator-(const Vec3D<T>& v, const U &u);
template <class T, class U>
const Vec3D<T> operator-(const U &u, const Vec3D<T>& v);
template <class T, class U>
const Vec3D<T> operator*(const Vec3D<T>& v, const U &u);
template <class T, class U>
const Vec3D<T> operator*(const U &u, const Vec3D<T>& v);
template <class T, class U>
const Vec3D<T> operator/(const Vec3D<T>& v, const U &u);
template <class T, class U>
const Vec3D<T> operator/(const U &u, const Vec3D<T>& v);

//
// input/output
//
template <class T>
std::ostream& operator<<(std::ostream& output, const Vec3D<T>& v);
template <class T>
std::istream& operator>>(std::istream& input, Vec3D<T>& v);

// component-wise max/min
template <class T>
const Vec3D<T> componentMax(const Vec3D<T>& a, const Vec3D<T>& b);

template <class T>
const Vec3D<T> componentMin(const Vec3D<T>& a, const Vec3D<T>& b);

#ifndef SWIG
#include "Vec3D.txx"
#endif // SWIG

// convenient types
typedef class Vec3D<int>          Vec3Di;
typedef class Vec3D<unsigned int> Vec3Du;
typedef class Vec3D<float>        Vec3Df;
typedef class Vec3D<double>       Vec3Dd;

} // end namespace PyCA

#endif

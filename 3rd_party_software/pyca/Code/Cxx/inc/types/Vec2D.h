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

#ifndef __VEC2D_H
#define __VEC2D_H

#include <AccumType.h>
#include <iostream>

namespace PyCA {

#ifdef __GNUC__
#define PYCA_ALIGNED __attribute__((__aligned__(8)))
#else
#define PYCA_ALIGNED ;
#endif

template<class T>
class Vec2D{
public:
    typedef T    value_type;
    typedef typename AccumType<T>::AccT acc_type;
    T x,y;

    explicit Vec2D(const T& x=0, const T& y=0):x(x), y(y){};

    Vec2D(const Vec2D<T>& rhs):x(rhs.x),y(rhs.y){}

    template<typename U>
    Vec2D(const Vec2D<U>& rhs):x(rhs.x),y(rhs.y){}

    Vec2D<T>& operator=(const Vec2D<T>& v) {
        x = v.x;
        y = v.y;
        return *this;
    }
    
    template<typename U>
    Vec2D<T>& operator=(const Vec2D<U>& v) {
        x = (T)v.x;
        y = (T)v.y;
        return *this;
    }

    bool operator==(const Vec2D& v) {
        return ((x==v.x) && (y==v.y));
    }

    bool operator!=(const Vec2D& v) {
        return ((x!=v.x) || (y!=v.y));
    }

    
    T& operator[](unsigned int index) {
        return (&x)[index];
    }
    
    const T& operator[](unsigned int index) const {
        return (&x)[index];
    }

    template<typename U> Vec2D<T>& operator +=(const Vec2D<U>& rhs) {
        x+=rhs.x; y+=rhs.y;
        return *this;
    }

    template<typename U>
    Vec2D<T>& operator -=(const Vec2D<U>& rhs) {
        x-=rhs.x; y-=rhs.y;
        return *this;
    }

    template<typename U> Vec2D<T>& operator *=(const Vec2D<U>& rhs) {
        x*=rhs.x; y*=rhs.y;
        return *this;
    }

    template<typename U> Vec2D<T>& operator /=(const Vec2D<U>& rhs) {
        x/=rhs.x; y/=rhs.y;
        return *this;
    }

    template<typename U> Vec2D<T>& operator +=(const U v) {
        x+=v; y+=v;
        return *this;
    }
    template<typename U> Vec2D<T>& operator -=(const U v) {
        x-=v; y-=v;
        return *this;
    }
    template<typename U> Vec2D<T>& operator *=(const U v) {
        x*=v; y*=v;
        return *this;
    }
    template<typename U> Vec2D<T>& operator /=(const U v) {
        x/=v; y/=v;
        return *this;
    }

    double length() const {
        return sqrt( x*x + y*y );
    }
    
} PYCA_ALIGNED;

template<typename T>
inline double length(const Vec2D<T> v) {
    return v.length();
}

template<typename T, typename U>
inline const Vec2D<T> operator +(const Vec2D<T>& lhs, const Vec2D<U>& rhs) {
    return Vec2D<T>(static_cast<T>(lhs.x + rhs.x), static_cast<T>(lhs.y + rhs.y));
}

template<typename T, typename U>
inline const Vec2D<T> operator *(const Vec2D<T>& lhs, const Vec2D<U>& rhs) {
    return Vec2D<T>(static_cast<T>(lhs.x * rhs.x), static_cast<T>(lhs.y * rhs.y));
}

template<typename T, typename U>
inline const Vec2D<T> operator *(const Vec2D<T>& lhs, const U& s) {
    return Vec2D<T>(static_cast<T>(lhs.x * s), static_cast<T>(lhs.y * s));
}

template<typename T, typename U>
inline const Vec2D<T> operator /(const Vec2D<T>& lhs, const Vec2D<U>& rhs) {
    return Vec2D<T>(static_cast<T>(lhs.x / rhs.x), static_cast<T>(lhs.y / rhs.y));
}

template<typename T, typename U>
inline const Vec2D<T> operator /(const Vec2D<T>& lhs, const U& s) {
    return Vec2D<T>(static_cast<T>(lhs.x / s), static_cast<T>(lhs.y / s));
}

template<typename T, typename U>
inline const Vec2D<T> operator -(const Vec2D<T>& lhs, const Vec2D<U>& rhs) {
    return Vec2D<T>(static_cast<T>(lhs.x - rhs.x), static_cast<T>(lhs.y - rhs.y));
}

template <class T>
inline std::ostream& operator<<(std::ostream& oss, const Vec2D<T>& v) {
    oss << "(" << v.x << ", " << v.y << ")";
    return oss;
}
template <class T>
inline std::istream& operator>>(std::istream& iss, Vec2D<T>& v) {
  char paren, comma;
  iss >> paren
	>> v.x >> comma
	>> v.y >> paren;
  return iss;
}
// convenient types
typedef Vec2D<int>          Vec2Di;
typedef Vec2D<unsigned int> Vec2Du;
typedef Vec2D<float>        Vec2Df;
typedef Vec2D<double>       Vec2Dd;

} // end namespace PyCA

#endif

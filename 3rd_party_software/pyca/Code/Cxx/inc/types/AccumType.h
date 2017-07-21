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

#ifndef __ACCUM_TYPE_H__
#define __ACCUM_TYPE_H__

namespace PyCA {

template<typename T>
class AccumType;

template<>
class AccumType<char>{
public:
    typedef int AccT;
};

template<>
class AccumType<unsigned char>{
public:
    typedef unsigned int AccT;
};

template<>
class AccumType<short>{
public:
    typedef int AccT;
};

template<>
class AccumType<unsigned short>{
public:
    typedef unsigned int AccT;
};

template<>
class AccumType<int>{
public:
    typedef long int AccT;
};

template<>
class AccumType<unsigned int>{
public:
    typedef unsigned long int AccT;
};

template<>
class AccumType<float>{
public:
    typedef double AccT;
};

template<>
class AccumType<double>{
public:
    typedef double AccT;
};

} // end namespace PyCA

#endif // __ACCUM_TYPE_H__

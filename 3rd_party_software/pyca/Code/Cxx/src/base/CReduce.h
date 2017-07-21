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

#ifndef __CREDUCE_H
#define __CREDUCE_H

#include <MOpers.h>
#include <estream.h>

namespace PyCA {

template<typename T>
class Vec2D;

class CReduce{
private:
    template<typename T, MATH_OPS op>
    void reduce(T* h_o, const T* h_i, size_t n, bool update, StreamT stream);

    template <typename T, MATH_OPS op, MATH_OPS op1>
    void compReduce(T* h_o, const T* h_i, size_t n, bool update, StreamT stream);

    template <typename T, MATH_OPS op, MATH_OPS op1>
    void bireduce(T* h_o, const T* h_i, size_t n, bool update, StreamT stream);
    
    template <typename T, MATH_OPS op, MATH_OPS op1>
    void product(T* h_o, const T* h_i0, const T* h_i1, size_t n, bool update, StreamT stream);
public:
    CReduce();
    template<typename T>
    void Max(T& h_o, const T* h_i, size_t n, bool update=false, StreamT stream=NULL);
    
    template<typename T>
    void Min(T& h_o, const T* h_i, size_t n, bool update=false, StreamT stream=NULL);
    template<typename T>
    void Sum(T& h_o, const T* h_i, size_t n, bool update=false, StreamT stream=NULL);

    template<typename T>
    void LInf(T& h_o, const T* h_i, size_t n, bool update=false, StreamT stream=NULL);
    template<typename T>
    void L1(T& h_o, const T* h_i, size_t n, bool update=false, StreamT stream=NULL);
    template<typename T>
    void Sum2(T& h_o, const T* h_i, size_t n, bool update=false, StreamT stream=NULL);

    template<class T>
    void MaxMin(Vec2D<T>& h_o, const T* h_i, size_t n, bool update=false, StreamT stream=NULL);

    template<typename T>
    void Dot(T& h_o, const T* h_i, const T* h_i1,  size_t n, bool update=false, StreamT stream=NULL);

    bool selfTest(size_t n);
};


//static __thread CReduce* __global_cpu_reduce;

} // end namespace PyCA

#endif


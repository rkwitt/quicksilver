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

#ifndef __GREDUCE_STREAM_H
#define __GREDUCE_STREAM_H

#include <MOpers.h>
#include <estream.h>

namespace PyCA {

template<typename T>
class Vec2D;


class GReduceStream{
private:
    int* mdBuf; // reduce buffer

    template<typename T, MATH_OPS op>
    void reduce(T* d_o, const T* d_i, size_t n, bool update, StreamT stream);

    template <typename T, MATH_OPS op, MATH_OPS op1>
    void compReduce(T* d_o, const T* d_i, size_t n, bool update, StreamT stream);

    template <typename T, MATH_OPS op, MATH_OPS op1>
    void bireduce(T* d_o, const T* d_i, size_t n, bool update, StreamT stream);
    
    template <typename T, MATH_OPS op, MATH_OPS op1>
    void product(T* d_o, const T* d_i0, const T* d_i1, size_t n, bool update, StreamT stream);

public:
    GReduceStream();
    ~GReduceStream();

    template<typename T>
    void Max(T& d_o, const T* d_i, size_t n, bool update=false, StreamT stream=NULL);
    
    template<typename T>
    void Min(T& d_o, const T* d_i, size_t n, bool update=false, StreamT stream=NULL);
    template<typename T>
    void Sum(T& d_o, const T* d_i, size_t n, bool update=false, StreamT stream=NULL);

    template<typename T>
    void LInf(T& d_o, const T* d_i, size_t n, bool update=false, StreamT stream=NULL);
    template<typename T>
    void L1(T& d_o, const T* d_i, size_t n, bool update=false, StreamT stream=NULL);
    template<typename T>
    void Sum2(T& d_o, const T* d_i, size_t n, bool update=false, StreamT stream=NULL);

    template<class T>
    void MaxMin(Vec2D<T>& d_o, const T* d_i, size_t n, bool update=false, StreamT stream=NULL);

    template<typename T>
    void Dot(T& d_o, const T* d_i, const T* d_i1,  size_t n, bool update=false, StreamT stream=NULL);

    bool selfTest(size_t n);
};

} // end namespace PyCA

#endif

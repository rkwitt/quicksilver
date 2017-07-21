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

#ifndef __REDUCE_H__
#define __REDUCE_H__

#include <Vec2D.h>
#include <pycaConst.h>
#include <estream.h>
#include <Selector.h>

#include <PyCAThread.h>

namespace PyCA {

class Field3D;

class CReduce;
class GReduce;
class GReduceStream;

template<int mode>
class Reduce{
public:
    typedef typename Selector<mode, CReduce, GReduce, GReduceStream>::Result ReduceExec;
    
    template<typename T>
    static void Max(T& a_o, const T* a_i, size_t n, bool update=false, StreamT stream=NULL);
    
    template<typename T>
    static void Min(T& a_o, const T* a_i, size_t n, bool update=false, StreamT stream=NULL);
    template<typename T>
    static void Sum(T& a_o, const T* a_i, size_t n, bool update=false, StreamT stream=NULL);

    template<typename T>
    static void LInf(T& a_o, const T* a_i, size_t n, bool update=false, StreamT stream=NULL);
    template<typename T>
    static void L1(T& a_o, const T* a_i, size_t n, bool update=false, StreamT stream=NULL);
    
    template<typename T>
    static void Sum2(T& a_o, const T* a_i, size_t n, bool update=false, StreamT stream=NULL);

    template<class T>
    static void MaxMin(Vec2D<T>& a_o, const T* a_i, size_t n, bool update=false, StreamT stream=NULL);

    template<typename T>
    static void Dot(T& a_o, const T* a_i, const T* a_i1,  size_t n, bool update=false, StreamT stream=NULL);


    static void Dot(float& a_o, const Field3D& a_i0, const Field3D& a_i1,
                    size_t n, bool update=false, StreamT stream=NULL);

    template<typename T>
    static void Max(T& a_o, const T** a_Arrays, size_t nA,
                    size_t n, bool update=false, StreamT stream=NULL);

    template<typename T>
    static void MaxMin(Vec2D<T>& a_o, const T** a_Arrays, size_t nA,
                       size_t n, bool update=false, StreamT stream=NULL);
    
    static bool selfTest(size_t n);

    static ReduceExec& GetReduceObj() {
	if (mImpl == NULL) {
	    mImpl = new ReduceExec();
	}
	return *mImpl;
    }
    static void Clear() {
	if (mImpl != NULL)
	    delete mImpl;
	mImpl = NULL;
    }
private:
    static __thread ReduceExec* mImpl;
};

} // end namespace PyCA

#endif

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

#include <Reduce.h>
#include "CReduce.h"
#include "GReduce.h"
#include "GReduceStream.h"
#include <Vec2D.h>
#include <MOpers.h>
#include <FieldOpers.h>
#include <Field3D.h>

namespace PyCA {

// template<int mode>
// typename Reduce<mode>::ReduceExec& Reduce<mode>::GetReduceObj() {
   
// }

template<int mode>
template<typename T>
void Reduce<mode>::Max(T& a_o, const T* a_i, size_t n, bool update, StreamT stream){
    GetReduceObj().Max(a_o, a_i, n, update, stream);
};

template<int mode>
template<typename T>
void Reduce<mode>::Min(T& a_o, const T* a_i, size_t n, bool update, StreamT stream){
    GetReduceObj().Min(a_o, a_i, n, update, stream);
};

template<int mode>
template<typename T>
void Reduce<mode>::Sum(T& a_o, const T* a_i, size_t n, bool update, StreamT stream){
    GetReduceObj().Sum(a_o, a_i, n, update, stream);
};

template<int mode>
template<typename T>
void Reduce<mode>::LInf(T& a_o, const T* a_i, size_t n, bool update, StreamT stream){
    GetReduceObj().LInf(a_o, a_i, n, update,stream);
};

template<int mode>
template<typename T>
void Reduce<mode>::L1(T& a_o, const T* a_i, size_t n, bool update, StreamT stream){
    GetReduceObj().L1(a_o, a_i, n, update,stream);
};

template<int mode>
template<typename T>
void Reduce<mode>::Sum2(T& a_o, const T* a_i, size_t n, bool update, StreamT stream){
    GetReduceObj().Sum2(a_o, a_i, n, update,stream);
};


template<int mode>
template<class T>
void Reduce<mode>::MaxMin(Vec2D<T>& a_o, const T* a_i, size_t n, bool update, StreamT stream){
    GetReduceObj().MaxMin(a_o, a_i, n, update,stream);
};


template<int mode>
template<typename T>
void Reduce<mode>::Dot(T& a_o, const T* a_i, const T* a_i1,  size_t n, bool update, StreamT stream){
    GetReduceObj().Dot(a_o, a_i, a_i1, n, update,stream);
};

template<int mode>
bool Reduce<mode>::selfTest(size_t n){
    return GetReduceObj().selfTest(n);
}

template<int mode>
void Reduce<mode>::
Dot(float& a_o, const Field3D& a_i0, const Field3D& a_i1,
    size_t n, bool update, StreamT stream){
    
    if (Is1D(a_i0, a_i1, n)){
        Dot(a_o, a_i0.x, a_i1.x, 3 * n, update, stream);
    } else {
        Dot(a_o, a_i0.x, a_i1.x, n, update, stream);
        Dot(a_o, a_i0.y, a_i1.y, n, true, stream);
        Dot(a_o, a_i0.z, a_i1.z, n, true, stream);
    }    
}

template<int mode>
template<typename T>
void Reduce<mode>::
Max(T& a_o, const T** a_Arrays, size_t nA,
    size_t n, bool update, StreamT stream){
   PRECONDITION(nA >= 1, "Number of elements = 1 for reduce call");
    
    Max(a_o, a_Arrays[0], n, update, stream);
    for (size_t i=1; i< nA; ++i)
        Max(a_o, a_Arrays[i], n, true, stream);
}

template<int mode>
template<typename T>
void Reduce<mode>::
MaxMin(Vec2D<T>& a_o, const T** a_Arrays, size_t nA,
       size_t n, bool update, StreamT stream){
   PRECONDITION(nA >= 1, "Number of elements = 1 for reduce call");
    MaxMin(a_o, a_Arrays[0], n, update, stream);
    for (size_t i=1; i< nA; ++i)
        MaxMin(a_o, a_Arrays[i], n, true, stream);
}

template<int mode>
__thread typename Reduce<mode>::ReduceExec* Reduce<mode>::mImpl=NULL;

// template instantiation
#include "Reduce_inst.cxx"

} // end namespace PyCA

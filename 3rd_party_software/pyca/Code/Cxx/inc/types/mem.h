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

#ifndef __MEM_H
#define __MEM_H

#ifndef SWIG
#include <stddef.h>
#include <memory>
#include <estream.h>
#include <conditionMacro.h>
#include <PyCAException.h>

// TEST -- make sure filed including boost aren't leaking into
// nvcc-compiled code
#if defined(PYCA_BOOSTTEST)
#if defined(__CUDACC__)
int bla[-1];
#endif
#endif
// END TEST

// needed for templated function args to AUTO_EXEC
#define COMMA ,

#ifdef CUDA_ENABLED
// execute correct version of a function based on memory type
#define AUTO_EXEC(mem_type, the_class, the_function)	                \
    {									\
	if(mem_type == MEM_DEVICE ||				        \
	   mem_type == MEM_HOST_PINNED){				\
	    the_class<EXEC_GPU>::the_function;			        \
	}else{								\
	    the_class<EXEC_CPU>::the_function;			        \
	}								\
    }
#else // CUDA_ENABLED

// execute correct version of a function based on memory type
#define AUTO_EXEC(mem_type, the_class, the_function)	                \
    {									\
	if(mem_type == MEM_DEVICE ||				        \
	   mem_type == MEM_HOST_PINNED){				\
	    throw PyCAException(__FILE__, __LINE__,                     \
                                "Error, GPU code not compiled");        \
	}else{								\
	    the_class<EXEC_CPU>::the_function;			        \
	}								\
    }

#endif // CUDA_ENABLED

#endif // SWIG

namespace PyCA {

// SWIG doesn't like the 'aligned' macro
#ifdef SWIG
#define __attribute__(bla) ;
#define __aligned__(bla) ;
#endif // SWIG

// MSVC doesn't like __attribute__, which is a gcc feature
#ifdef _MSC_VER
#define __attribute__(bla) ;
#endif


template<class T>
struct ComplexT {
    T x;
    T y;
} __attribute__((__aligned__(8)));

typedef ComplexT<float> Complex;

// Three different type of memory
enum MemoryType {MEM_UNINITIALIZED, 
		 MEM_HOST, 
		 MEM_HOST_PINNED, 
		 MEM_DEVICE};


const char* MemTypeStr(MemoryType type);

/*
 * Allocation/deallocation for CPU/GPU memory
 */
template<typename T>
void dmemAlloc(T*& d_ptr, size_t n);

template<typename T>
void dmemFree(T*& d_ptr);

template<typename T>
void phmemAlloc(T*& ph_ptr, size_t n);

template<typename T>
void phmemFree(T*& ph_ptr);

template<MemoryType mType, typename T>
void memAlloc(T*& p, size_t n);

template<MemoryType mType, typename T>
void memFree(T*& p);

template<typename T>
void memAlloc(T*& p, size_t n, MemoryType mType);

template<typename T>
void memFree(T*& p, MemoryType mType);

template<MemoryType mType, typename T>
T* memAlloc(size_t n);

template<typename T>
T* memAlloc(size_t n, MemoryType mType);

template<MemoryType mType, typename T>
void ArrayAllocator(T*& p, size_t n);

template<MemoryType mType, typename T>
void ArrayDeleter(T*& p);

/*
 * Memory copy function
 */
template<class T>
void cpyArrayH2H(T* a, const T* b, size_t n);

template<class T>
void cpyArrayH2D(T* d_a, const T* h_b, size_t n);

template<class T>
void cpyArrayD2H(T* h_a, const T* d_b, size_t n);

template<class T>
void cpyArrayD2D(T* d_a, const T* d_b, size_t n);

template<class T>
void cpyArrayH2H(T* a, const T* b, size_t n, StreamT s);

template<class T>
void acpyArrayH2D(T* d_a, const T* h_b, size_t n, StreamT s);

template<class T>
void acpyArrayD2H(T* h_a, const T* d_b, size_t n, StreamT s);

template<class T>
void acpyArrayD2D(T* d_a, const T* d_b, size_t n, StreamT s);

template<typename T>
void memCopy(T* a_o, const T* a_i, size_t n, MemoryType oType, MemoryType iType);

template<typename T>
void memCopy(T* a_o, const T* a_i, size_t n, MemoryType oType, MemoryType iType, StreamT stream);

/*
 * Copy function for constant or symbol
 */
template<class T>
void cpyH2C(const char* name, const T* h_i, size_t cnt);

template<class T>
void cpyC2H(T* h_i, const char* name, size_t cnt);

template<class T>
void cpyD2C(const char* name, const T* d_i, size_t cnt);

template<class T>
void cpyC2D(T* d_i, const char* name, size_t cnt);

template<class T>
void acpyH2C(const char* name, const T* h_i, size_t cnt, StreamT s);

template<class T>
void acpyC2H(T* h_i, const char* name, size_t cnt, StreamT s);

template<class T>
void acpyD2C(const char* name, const T* d_i, size_t cnt, StreamT s);

template<class T>
void acpyC2D(T* d_i, const char* name, size_t cnt, StreamT s);

} // end namespace PyCA

#ifndef SWIG
#include "mem.txx"
#endif // SWIG

#endif

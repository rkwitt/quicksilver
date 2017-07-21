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

#ifndef __MEM_INC
#define __MEM_INC

#include <cstring>

#include "CudaUtils.h"

namespace PyCA {

#ifdef CUDA_ENABLED

template<typename T>
inline void dmemAlloc(T*& d_ptr, size_t n)
{
    cudaMalloc((void**)&d_ptr, n * sizeof(T));
}

template<typename T>
inline void dmemFree(T*& d_ptr)
{
    if (d_ptr) {
        cudaFree(d_ptr);
        d_ptr = NULL;
    }
}

template<typename T>
inline void phmemAlloc(T*& ph_ptr, size_t n)
{
    cudaHostAlloc((void**)&ph_ptr, n * sizeof(T), cudaHostAllocPortable);
}

template<typename T>
inline void phmemFree(T*& ph_ptr)
{
    if (ph_ptr) {
        cudaFreeHost(ph_ptr);
        ph_ptr = NULL;
    }
}

template<MemoryType mType, typename T>
inline void memAlloc(T*& p, size_t n){
    if (mType == MEM_HOST)
        p = new T [n];
    else if (mType == MEM_HOST_PINNED)
        cudaHostAlloc((void**)&p, n * sizeof(T), cudaHostAllocPortable);
    else
        cudaMalloc((void**)&p, n * sizeof(T));
}

template<MemoryType mType, typename T>
inline void memFree(T*& p){
    if (p) {
        if (mType == MEM_HOST)
            delete []p;
        else if (mType == MEM_HOST_PINNED)
            cudaFreeHost(p);
        else
            cudaFree(p);
        p = NULL;
    }
}

template<typename T>
inline void memAlloc(T*& p, size_t n, MemoryType mType){
    if (mType == MEM_HOST)
        p = new T [n];
    else if (mType == MEM_HOST_PINNED)
        cudaHostAlloc((void**)&p, n * sizeof(T), cudaHostAllocPortable);
    else
        cudaMalloc((void**)&p, n * sizeof(T));
}

template<typename T>
inline void memFree(T*& p, MemoryType mType){
    if (p) {
        if (mType == MEM_HOST)
            delete []p;
        else if (mType == MEM_HOST_PINNED)
            cudaFreeHost(p);
        else
            cudaFree(p);
        p = NULL;
    }
}

#else // CUDA_ENABLED

template<typename T>
inline void dmemAlloc(T*& d_ptr, size_t n)
{
   throw PyCAException(__FILE__, __LINE__, "CUDA support not compiled");
}

template<typename T>
inline void dmemFree(T*& d_ptr)
{
   throw PyCAException(__FILE__, __LINE__, "CUDA support not compiled");
}

template<typename T>
inline void phmemAlloc(T*& ph_ptr, size_t n)
{
   throw PyCAException(__FILE__, __LINE__, "CUDA support not compiled");
}

template<typename T>
inline void phmemFree(T*& ph_ptr)
{
   throw PyCAException(__FILE__, __LINE__, "CUDA support not compiled");
}

template<MemoryType mType, typename T>
inline void memAlloc(T*& p, size_t n){
    if (mType == MEM_HOST)
        p = new T [n];
    else
       throw PyCAException(__FILE__, __LINE__, "CUDA support not compiled");
}

template<MemoryType mType, typename T>
inline void memFree(T*& p){
    if (p) {
        if (mType == MEM_HOST)
            delete []p;
        else
	   throw PyCAException(__FILE__, __LINE__, "CUDA support not compiled");
        p = NULL;
    }
}

template<typename T>
inline void memAlloc(T*& p, size_t n, MemoryType mType){
    if (mType == MEM_HOST)
        p = new T [n];
    else
       throw PyCAException(__FILE__, __LINE__, "CUDA support not compiled");
}

template<typename T>
inline void memFree(T*& p, MemoryType mType){
    if (p) {
        if (mType == MEM_HOST)
            delete []p;
        else
	   throw PyCAException(__FILE__, __LINE__, "CUDA support not compiled");
        p = NULL;
    }
}

#endif // CUDA_ENABLED

template<MemoryType mType, typename T>
inline T* memAlloc(size_t n){
    T* p;
    memAlloc<mType, T>(p, n);
    return p;
}

template<typename T>
inline T* memAlloc(size_t n, MemoryType mType){
     T* p;  
     memAlloc<T>(p, n, mType);
     return p;
}

template<MemoryType mType, typename T>
inline T* ArrayAllocator(size_t n){
       return memAlloc<mType, T>(n);
}

template<MemoryType mType, typename T>
inline void ArrayDeleter(T*& p){
    memFree<mType, T>(p);
}

/*
 * Memory copy function
 */
template<class T>
inline void cpyArrayH2H(T* a, const T* b, size_t n) {
    memcpy(a, b, sizeof(T) * n);
}

template<class T>
inline void acpyArrayH2H(T* a, const T* b, size_t n, StreamT s) {
    memcpy(a, b, sizeof(T) * n);
}

#ifdef CUDA_ENABLED

template<class T>
inline void cpyArrayH2D(T* d_a, const T* h_b, size_t n){
    cudaMemcpy(d_a, h_b, n * sizeof(T), cudaMemcpyHostToDevice);
    CudaUtils::CheckCUDAError(__FILE__,__LINE__, 
			      "synchronous copy host to device");
}

template<class T>
inline void cpyArrayD2H(T* h_a, const T* d_b, size_t n){
    cudaMemcpy(h_a, d_b, n * sizeof(T), cudaMemcpyDeviceToHost);
    CudaUtils::CheckCUDAError(__FILE__,__LINE__,
			      "synchronous copy device to host");
}

template<class T>
inline void cpyArrayD2D(T* d_a, const T* d_b, size_t n){
    cudaMemcpy(d_a, d_b, n * sizeof(T), cudaMemcpyDeviceToDevice);
    CudaUtils::CheckCUDAError(__FILE__,__LINE__,
			      "synchronous copy device to device");
}

template<class T>
inline void acpyArrayH2D(T* d_a, const T* h_b, size_t n, StreamT s){
    cudaMemcpyAsync(d_a, h_b, n * sizeof(T), cudaMemcpyHostToDevice, s);
    CudaUtils::CheckCUDAError(__FILE__,__LINE__,
			      "asynchronous copy host to device");
}

template<class T>
inline void acpyArrayD2H(T* h_a, const T* d_b, size_t n, StreamT s){
    cudaMemcpyAsync(h_a, d_b, n * sizeof(T), cudaMemcpyDeviceToHost, s);
    CudaUtils::CheckCUDAError(__FILE__,__LINE__,
			      "asynchronous copy device to host");
}

template<class T>
inline void acpyArrayD2D(T* d_a, const T* d_b, size_t n, StreamT s){
    cudaMemcpyAsync(d_a, d_b, n * sizeof(T), cudaMemcpyDeviceToDevice, s);
    CudaUtils::CheckCUDAError(__FILE__,__LINE__,
			      "asynchronous copy device to device");
}

#else // CUDA_ENABLED

template<class T>
inline void cpyArrayH2D(T* d_a, const T* h_b, size_t n){
   throw PyCAException(__FILE__, __LINE__, "CUDA support not compiled");
}

template<class T>
inline void cpyArrayD2H(T* h_a, const T* d_b, size_t n){
   throw PyCAException(__FILE__, __LINE__, "CUDA support not compiled");
}

template<class T>
inline void cpyArrayD2D(T* d_a, const T* d_b, size_t n){
   throw PyCAException(__FILE__, __LINE__, "CUDA support not compiled");
}

template<class T>
inline void acpyArrayH2D(T* d_a, const T* h_b, size_t n, StreamT s){
   throw PyCAException(__FILE__, __LINE__, "CUDA support not compiled");
}

template<class T>
inline void acpyArrayD2H(T* h_a, const T* d_b, size_t n, StreamT s){
   throw PyCAException(__FILE__, __LINE__, "CUDA support not compiled");
}

template<class T>
inline void acpyArrayD2D(T* d_a, const T* d_b, size_t n, StreamT s){
   throw PyCAException(__FILE__, __LINE__, "CUDA support not compiled");
}

#endif // CUDA_ENABLED

template<typename T>
inline void memCopy(T* a_o, const T* a_i, size_t n, MemoryType oType, MemoryType iType) {
    if (oType == MEM_DEVICE && iType == MEM_DEVICE) {
        cpyArrayD2D(a_o, a_i, n);
        return;
    }

    if (oType != MEM_DEVICE && iType != MEM_DEVICE) {
        cpyArrayH2H(a_o, a_i, n);
        return;
    }
    if (oType == MEM_DEVICE)
        cpyArrayH2D(a_o, a_i, n);
    else
        cpyArrayD2H(a_o, a_i, n);
}


template<typename T>
inline void memCopy(T* a_o, const T* a_i, size_t n, MemoryType oType, MemoryType iType, StreamT stream) {
    if (stream == NULL) {
        memCopy(a_o, a_i, n, oType, iType);
        return;
    }
    
    if ((oType == MEM_HOST && iType == MEM_DEVICE) || (oType == MEM_DEVICE && iType == MEM_HOST) )
        throw "Asynchronous memory coppy required host pinned memory ";

    if (oType == MEM_DEVICE && iType == MEM_DEVICE) {
        acpyArrayD2D(a_o, a_i, n, stream);
        return;
    }

    if (oType != MEM_DEVICE && iType != MEM_DEVICE) {
        cpyArrayH2H(a_o, a_i, n);
        return;
    }

    if (oType == MEM_DEVICE)
        acpyArrayH2D(a_o, a_i, n, stream);
    else //oType == MEM_HOST_PINNED
        acpyArrayD2H(a_o, a_i, n, stream);
}

/*
 *  Memory copy from to constant memory
 */

#ifdef CUDA_ENABLED

template<class T>
inline void cpyH2C(const char* name, const T* h_i, size_t cnt){
    cudaMemcpyToSymbol(name, h_i, cnt * sizeof(T), 0, cudaMemcpyHostToDevice);
    CudaUtils::CheckCUDAError(__FILE__,__LINE__,
			      "Host to constant");
}

template<class T>
inline void cpyC2H(T* h_o, const char* name, size_t cnt){
    cudaMemcpyFromSymbol(h_o, name, cnt * sizeof(T), 0, cudaMemcpyDeviceToHost);
    CudaUtils::CheckCUDAError(__FILE__,__LINE__,
			      "Constant to Host");
}

template<class T>
inline void cpyD2C(const char* name, const T* d_i, size_t cnt){
    cudaMemcpyToSymbol(name, d_i, cnt * sizeof(T), 0, cudaMemcpyDeviceToDevice);
    CudaUtils::CheckCUDAError(__FILE__,__LINE__,
			      "Device to constant");
}

template<class T>
inline void cpyC2D(T* d_o, const char* name, size_t cnt){
    cudaMemcpyFromSymbol(d_o, name, cnt * sizeof(T), 0, cudaMemcpyDeviceToDevice);
    CudaUtils::CheckCUDAError(__FILE__,__LINE__,
			      "Constant to device");
}


template<class T>
inline void acpyH2C(const char* name, const T* h_i, size_t cnt, StreamT s){
    cudaMemcpyToSymbolAsync(name, h_i, cnt * sizeof(T), 0, cudaMemcpyHostToDevice, s);
    CudaUtils::CheckCUDAError(__FILE__,__LINE__,
			      "Host to constant - host should be pinned");
}

template<class T>
inline void acpyC2H(T* h_o, const char* name, size_t cnt, StreamT s){
    cudaMemcpyFromSymbolAsync(h_o, name, cnt * sizeof(T), 0, cudaMemcpyDeviceToHost,s);
    CudaUtils::CheckCUDAError(__FILE__,__LINE__,
			      "Constant to Host - host should be pinned");
}

template<class T>
inline void acpyD2C(const char* name, const T* d_i, size_t cnt, StreamT s){
    cudaMemcpyToSymbolAsync(name, d_i, cnt * sizeof(T), 0, cudaMemcpyDeviceToDevice,s);
    CudaUtils::CheckCUDAError(__FILE__,__LINE__,
			      "Device to constant");
}

template<class T>
inline void acpyC2D(T* d_o, const char* name, size_t cnt, StreamT s){
    cudaMemcpyFromSymbolAsync(d_o, name, cnt * sizeof(T), 0, cudaMemcpyDeviceToDevice,s);
    CudaUtils::CheckCUDAError(__FILE__,__LINE__,
			      "Constant to device");
}

#endif // CUDA_ENABLED

} // end namespace PyCA

#endif // __MEM_INC

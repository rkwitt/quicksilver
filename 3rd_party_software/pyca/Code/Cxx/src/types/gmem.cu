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

#include <gmem.h>
#include <mem.h>
#include <pycaConst.h>
#include <iostream>
#include <Vec2D.h>
#include <stdint.h>

#include <cuda.h>

// TEST make sure boost isn't included in nvcc code
#if defined(BOOST_COMPILER)
int bla[-1];
#endif

namespace PyCA {

template<typename T>
T* createGObj(const T& h_v){
    T* pObj;
    cudaMalloc((void**)&pObj, sizeof(T));
    setGObj(pObj, h_v);
    return pObj;
}

template<typename T>
void setGObj(T* d_obj, const T& h_v){
   cudaMemcpy(d_obj, &h_v, sizeof(T), cudaMemcpyHostToDevice);
}

template<typename T>
T getGObj(T* d_obj){
    T h_v;
    cudaMemcpy(&h_v, d_obj, sizeof(T), cudaMemcpyDeviceToHost);
    return h_v;
}

template Vec2Di* createGObj(const Vec2Di& v);
template void setGObj(Vec2Di* d_obj, const Vec2Di& h_v);
template Vec2Di getGObj(Vec2Di* d_obj);

template Vec2Df* createGObj(const Vec2Df& v);
template void setGObj(Vec2Df* d_obj, const Vec2Df& h_v);
template Vec2Df getGObj(Vec2Df* d_obj);

template float* createGObj(const float& v);
template void setGObj(float* d_obj, const float& h_v);
template float getGObj(float* d_obj);

template int* createGObj(const int& v);
template void setGObj(int* d_obj, const int& h_v);
template int getGObj(int* d_obj);

template uint* createGObj(const uint& v);
template void setGObj(uint* d_obj, const uint& h_v);
template uint getGObj(uint* d_obj);

// template<int mode>
// void IDiv(float* a_o, float a_v){
//     a_o[0] = a_v / a_o[0];
// }

// __global__ void IDiv_kernel(float* d_o, float v) {
//     uint      id = threadIdx.x;
//     if (id == 0)
//         d_o[0] = v / d_o[0];
// }

// template<>
// void IDiv<EXEC_GPU_PARAM>(float* d_o, float v, StreamT stream){
//     dim3 threads(32);
//     dim3 grids(1,0);
//     IDiv_kernel<<<grids, threads, 0, stream>>>(d_o, v);
// }

} // end namespace PyCA


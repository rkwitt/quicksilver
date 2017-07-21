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

#ifndef __PYCA_UTILS_H
#define __PYCA_UTILS_H

#include <pycaConst.h>
#include <conditionMacro.h>

namespace PyCA {

typedef unsigned int uint;

// Check if a point is inside 3D boundary
inline 
__HOSTDEVICE__ 
bool isInside3D(int x, int y, int z, 
		int w, int h, int l)
{
    return ((x >= 0) && (x < w) &&
            (y >= 0) && (y < h) &&
            (z >= 0) && (z < l));
}

inline 
__HOSTDEVICE__ 
int S2p20(float a){
    return int(a* FIX_SCALE_20 + 0.5f);
}

inline 
__HOSTDEVICE__ 
float S2n20(int a){
    return (float)a / FIX_SCALE_20;
}




// divide to the next integer
inline int iDivUp(int a, int b){
    return (a + b - 1) / b;
}

//Align a to nearest higher multiple of b
inline int iAlignUp(int a, int b){
    return (a % b != 0) ?  (a - a % b + b) : a;
}

// bits manipulation
inline int nextPowerOf2(int v){
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    return v + 1;
}

inline unsigned int  log2i(unsigned int v) {
    static const unsigned int MultiplyDeBruijnBitPosition[32] = 
        {
            0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8, 
            31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
        };
    v |= v >> 1; // first round down to power of 2 
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v = (v >> 1) + 1;
    return MultiplyDeBruijnBitPosition[(v * 0x077CB531UL) >> 27];
}

inline int isPowerOf2(int v) {
    return ((v & (v - 1)) == 0);
}

template<class T>
inline 
__HOSTDEVICE__ 
T
PYCAMAX(T a, T b)
{
   return a > b ? a : b;
}

template<class T>
inline 
__HOSTDEVICE__ 
T
PYCAMIN(T a, T b)
{
   return a < b ? a : b;
}

//
// unchecked access to array, need stride&offset to operate on vector fields
//
// non-const version
template< class T, int stride, int off >
inline
__HOSTDEVICE__ 
T& getSOVal(T* array, 
	  int sx, int sy, int sz,
	  int x, int y, int z)
{
  return array[stride*(sx*(z*sy + y) + x) + off];
}

// const version
template< class T, int stride, int off >
inline
__HOSTDEVICE__ 
const T& 
getSOVal(const T* array, 
       int sx, int sy, int sz,
       int x, int y, int z)
{
  return array[stride*(sx*(z*sy + y) + x) + off];
}

//
// versions without stride and offset params
//
// non-const version
template< class T>
inline
__HOSTDEVICE__ 
T& getVal(T* array, 
	  int sx, int sy, int sz,
	  int x, int y, int z)
{
   return getSOVal<T,1,0>(array, sx, sy, sz, x, y, z);
}

// const version
template< class T >
inline
__HOSTDEVICE__ 
const T& 
getVal(const T* array, 
       int sx, int sy, int sz,
       int x, int y, int z)
{
   return getSOVal<T,1,0>(array, sx, sy, sz, x, y, z);
}

template<class T, BackgroundStrategy backgroundStrategy>
static
inline
__HOSTDEVICE__ 
T
getSafeVal(const T *data, 
	   int sx, int sy, int sz,
	   int x, int y, int z,
	   T bgval = static_cast<T>(0))

{
    // MK_CHECK_IMAGE_BACKGROUND(backgroundStrategy)

    if (backgroundStrategy == BACKGROUND_STRATEGY_WRAP){
	// modding x+sx instead of x because modding negative numbers
	// returns a negative number in c
	// jsp 12/2012
	return getVal(data, 
		      sx,sy,sz,
		      ((x+sx)%sx),((y+sy)%sy),((z+sz)%sz));
    }else if (backgroundStrategy == BACKGROUND_STRATEGY_CLAMP){
	// clamp
	return getVal(data,
		      sx,sy,sz,
		      PYCAMIN(PYCAMAX(x, 0),sx-1),
		      PYCAMIN(PYCAMAX(y, 0),sy-1),
		      PYCAMIN(PYCAMAX(z, 0),sz-1));
    }else if (backgroundStrategy == BACKGROUND_STRATEGY_VAL ||
	      backgroundStrategy == BACKGROUND_STRATEGY_ZERO ||
	      backgroundStrategy == BACKGROUND_STRATEGY_PARTIAL_ZERO){

	if(backgroundStrategy == BACKGROUND_STRATEGY_ZERO ||
	   backgroundStrategy == BACKGROUND_STRATEGY_PARTIAL_ZERO){
	    bgval = 0.f;
	}

	if (x >= 0 && x < sx &&
	    y >= 0 && y < sy &&
	    z >= 0 && z < sz){
	    return getVal(data,
			  sx,sy,sz,
			  x,y,z);
	}else{
	    return bgval;
	}
    }else{
	// unknown background strategy, don't allow compilation
	STATIC_ASSERT(backgroundStrategy== BACKGROUND_STRATEGY_WRAP ||
		      backgroundStrategy== BACKGROUND_STRATEGY_CLAMP ||
		      backgroundStrategy== BACKGROUND_STRATEGY_ZERO ||
		      backgroundStrategy== BACKGROUND_STRATEGY_PARTIAL_ZERO ||
		      backgroundStrategy== BACKGROUND_STRATEGY_VAL);
	return 0.f;
    }
}

// compute determinant for 3x3 matrix
inline
__HOSTDEVICE__ 
float det(float a00, float a01, float a02,
	  float a10, float a11, float a12,
	  float a20, float a21, float a22)
{
   return a00 * a11 * a22 + a01 * a12 * a20 + a02 * a10 * a21 -
      a02 * a11 * a20 - a00 * a12 * a21 - a01 * a10 * a22;
}

// compute determinant for 2x2 matrix
inline
__HOSTDEVICE__ 
float det(float a00, float a01,
	  float a10, float a11)
{
   return a00 * a11 - a10 * a01;
}

// ================================================================
// CUDA-only functions
// ================================================================

#ifdef __CUDACC__

inline void checkConfig(dim3& grids){
    int nBlocks = grids.x;
    if (nBlocks >= (1 << 16)){
        int bly = nBlocks;
        int blx = 1;
        while (bly >= (1 << 16)){
            bly >>=1;
            blx <<=1;
        }
        grids.x = blx;
        grids.y = (blx * bly == nBlocks) ? bly : bly + 1;
    }
}

inline dim3 make_grid(uint num_blocks){
    if (num_blocks <= 65535){
        return dim3(num_blocks);
    } else {
        dim3 grids(num_blocks);
        checkConfig(grids);
        return grids;
    }
}

/**
 * Get the block ID from the config
 */
inline __device__ uint get_blockID(){
    return blockIdx.x + blockIdx.y * gridDim.x;
}

/**
 * Get the thread ID from the current config with current block
 */
inline __device__ uint get_threadID(uint blockId){
    return blockId * blockDim.x + threadIdx.x;
}

              
/**
 * @brief Common function to compute the id of the threads in side a block 
 * for 1D/2D setup 
 */
inline __device__ unsigned int get_x_2D(){
    return blockIdx.x * blockDim.x + threadIdx.x;
}

inline __device__ unsigned int get_y_2D(){
    return blockIdx.y * blockDim.y + threadIdx.y;
}

inline __device__ unsigned int get_id_2D(int w){
    return get_y_2D() * w + get_x_2D();
}

inline __device__ __host__ int2 S2p20(float2 a){
    return make_int2((int)(a.x *FIX_SCALE_20 + 0.5f),
                     (int)(a.y *FIX_SCALE_20 + 0.5f));
}

inline __device__ __host__ int3 S2p20(float3 a){
    return make_int3((int)(a.x *FIX_SCALE_20 + 0.5f),
                     (int)(a.y *FIX_SCALE_20 + 0.5f),
                     (int)(a.z *FIX_SCALE_20 + 0.5f));
}

inline __device__ __host__ int4 S2p20(float4 a){
    return make_int4((int)(a.x *FIX_SCALE_20 + 0.5f),
                     (int)(a.y *FIX_SCALE_20 + 0.5f),
                     (int)(a.z *FIX_SCALE_20 + 0.5f),
                     (int)(a.w *FIX_SCALE_20 + 0.5f));
}


inline __device__ __host__ float2 S2n20(int2 a){
    return make_float2((float)a.x / FIX_SCALE_20, (float)a.y / FIX_SCALE_20);
}

inline __device__ __host__ float3 S2n20(int3 a){
    return make_float3((float)a.x / FIX_SCALE_20, (float)a.y / FIX_SCALE_20, (float)a.z / FIX_SCALE_20);
}

inline __device__ __host__ float4 S2n20(int4 a){
    return make_float4((float)a.x / FIX_SCALE_20,
                       (float)a.y / FIX_SCALE_20,
                       (float)a.z / FIX_SCALE_20,
                       (float)a.w / FIX_SCALE_20);
}

#endif


} // end namespace PyCA

#endif


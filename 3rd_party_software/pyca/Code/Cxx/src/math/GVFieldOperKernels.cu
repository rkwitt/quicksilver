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

#include "GVFieldOperKernels.h"
#include "GSplat.h"
#include <pycaUtils.h>

// TEST make sure boost isn't included in nvcc code
#if defined(BOOST_COMPILER)
int bla[-1];
#endif

namespace PyCA {

//////////////////////////////////////////////////////////////////////// 
// set Vfield to identity
// i.e. v(x) = 0
////////////////////////////////////////////////////////////////////////
__global__ void SetToZero_kernel(float* d_vx, float* d_vy, float* d_vz,
                                 int w, int h, int l){
    uint i     = blockIdx.x * blockDim.x + threadIdx.x;
    uint j     = blockIdx.y * blockDim.y + threadIdx.y;
    uint index = j * w + i;
    if (i < w && j < h){
        for (int k=0; k<l; ++k, index+=w*h){
            d_vx[index] = 0;
            d_vy[index] = 0;
            d_vz[index] = 0;
        }
    }
}


void 
SetToZero(float *d_vx, float *d_vy, float *d_vz, 
	  int sz_x, int sz_y, int sz_z,
	  StreamT stream) 
{

    dim3 threads(16,16);
    dim3 grids(iDivUp(sz_x, threads.x), iDivUp(sz_y, threads.y));
    SetToZero_kernel<<<grids, threads, 0, stream>>>
	(d_vx, d_vy, d_vz, sz_x, sz_y, sz_z);
}

//////////////////////////////////////////////////////////////////////// 
// Convert velocity field to hfield
// i.e. h(x) = x + v(x)
////////////////////////////////////////////////////////////////////////
__global__ void velocityToH_kernel(float* d_hx, float* d_hy, float* d_hz,
                                   const float* d_vx, const float* d_vy, const float* d_vz,
                                   float delta, int w, int h, int l,
                                   float ispX, float ispY, float ispZ){
    uint i     = blockIdx.x * blockDim.x + threadIdx.x;
    uint j     = blockIdx.y * blockDim.y + threadIdx.y;
    uint index = j * w + i;
    
    if (i < w && j < h){
        for (int k=0; k<l; ++k, index += w*h){
            d_hx[index] = i + delta * ispX * d_vx[index];
            d_hy[index] = j + delta * ispY * d_vy[index];
            d_hz[index] = k + delta * ispZ * d_vz[index];
        }
    }
}

__device__ __constant__ bool c_delta;

__global__ void velocityToH_const_kernel(
    float* d_hx, float* d_hy, float* d_hz,
    const float* d_vx, const float* d_vy, const float* d_vz,
    int w, int h, int l,
    float ispX, float ispY, float ispZ)
{
    uint i     = blockIdx.x * blockDim.x + threadIdx.x;
    uint j     = blockIdx.y * blockDim.y + threadIdx.y;
    uint index = j * w + i;

    if (i < w && j < h){
        float delta = c_delta;
        for (int k=0; k<l; ++k, index += w*h){
            d_hx[index] = i + delta * ispX * d_vx[index];
            d_hy[index] = j + delta * ispY * d_vy[index];
            d_hz[index] = k + delta * ispZ * d_vz[index];
        }
    }
}


void
toH(float *d_hx, float *d_hy, float *d_hz, 
    float sp_x, float sp_y, float sp_z,
    const float *d_vx, const float *d_vy, const float *d_vz,
    int sz_x, int sz_y, int sz_z,
    const float& delta, 
    StreamT stream, bool onDev)
{
    float isp_x  = 1.f/sp_x;
    float isp_y  = 1.f/sp_y;
    float isp_z  = 1.f/sp_z;

    dim3 threads(16,16);
    dim3 grids(iDivUp(sz_x, threads.x), iDivUp(sz_y, threads.y));
    if (onDev) {
        cudaMemcpyToSymbolAsync(c_delta,&delta,sizeof(float),
                                0,cudaMemcpyDeviceToDevice,stream);
        velocityToH_const_kernel<<<grids, threads, 0, stream>>>
            (d_hx, d_hy, d_hz,
             d_vx, d_vy, d_vz, 
             sz_x, sz_y, sz_z,
             isp_x, isp_y, isp_z);
    } else {
        velocityToH_kernel<<<grids, threads, 0, stream>>>
            (d_hx, d_hy, d_hz,
             d_vx, d_vy, d_vz, delta,
             sz_x, sz_y, sz_z,
             isp_x, isp_y, isp_z);
    }
}

//////////////////////////////////////////////////////////////////////// 
// Convert velocity field to hfield (inplace version)
////////////////////////////////////////////////////////////////////////
__global__ void velocityToH_I_kernel(float* d_hx, float* d_hy, float* d_hz,
                                     float delta,
                                     int w, int h, int l,
                                     float ispX, float ispY, float ispZ){
    uint i     = blockIdx.x * blockDim.x + threadIdx.x;
    uint j     = blockIdx.y * blockDim.y + threadIdx.y;
    uint index = j * w + i;
    
    if (i < w && j < h){
        for (int k=0; k<l; ++k, index+=w*h){
            d_hx[index] = float(i) + delta * ispX * d_hx[index];
            d_hy[index] = float(j) + delta * ispY * d_hy[index];
            d_hz[index] = float(k) + delta * ispZ * d_hz[index];
        }
    }
}

__global__ void velocityToH_I_const_kernel(float* d_hx, float* d_hy, float* d_hz,
                                           int w, int h, int l,
                                           float ispX, float ispY, float ispZ){
    uint i     = blockIdx.x * blockDim.x + threadIdx.x;
    uint j     = blockIdx.y * blockDim.y + threadIdx.y;
    uint index = j * w + i;

    float delta = c_delta;
    if (i < w && j < h){
        for (int k=0; k<l; ++k, index+=w*h){
            d_hx[index] = float(i) + delta * ispX * d_hx[index];
            d_hy[index] = float(j) + delta * ispY * d_hy[index];
            d_hz[index] = float(k) + delta * ispZ * d_hz[index];
        }
    }
}


void
toH_I(float *d_hx, float *d_hy, float *d_hz, 
      float sp_x, float sp_y, float sp_z,
      int sz_x, int sz_y, int sz_z,
      const float& delta, 
      StreamT stream, bool onDev)
{
    float isp_x  = 1.f/sp_x;
    float isp_y  = 1.f/sp_y;
    float isp_z  = 1.f/sp_z;

    dim3 threads(16,16);
    dim3 grids(iDivUp(sz_x, threads.x), iDivUp(sz_y, threads.y));
    if (onDev) {
        cudaMemcpyToSymbolAsync(c_delta,&delta,sizeof(float),
                                0,cudaMemcpyDeviceToDevice,stream);
        velocityToH_I_const_kernel<<<grids, threads, 0, stream>>>
            (d_hx, d_hy, d_hz,
             sz_x, sz_y, sz_z,
             isp_x, isp_y, isp_z);
    } else {
        velocityToH_I_kernel<<<grids, threads, 0, stream>>>
            (d_hx, d_hy, d_hz, delta,
             sz_x, sz_y, sz_z,
             isp_x, isp_y, isp_z);
    }
}

void 
splatUnnormalized
(float* d_iwdx, float* d_iwdy, float* d_iwdz, 
 size_t sizeX, size_t sizeY, size_t sizeZ,
 const float* d_wx, const float* d_wy, const float* d_wz,
 const float* d_px , const float* d_py, const float* d_pz,
 size_t nP, StreamT stream)
{
    int* i_dox = (int*) d_iwdx;
    int* i_doy = (int*) d_iwdy;
    int* i_doz = (int*) d_iwdz;
    
    // Just splat it 
    Splatting::splat3D(i_dox, i_doy, i_doz, 
		       sizeX, sizeY, sizeZ,
		       d_wx, d_wy, d_wz,
		       d_px, d_py, d_pz, nP, stream);
    // and convert directly
    Splatting::FixedToFloating_I(d_iwdx, d_iwdy, d_iwdz, nP, stream);
    
}

void 
splatNormalized
(float* d_iwdx, float* d_iwdy, float* d_iwdz, 
 int *i_dd,
 size_t sizeX, size_t sizeY, size_t sizeZ,
 const float* d_wx, const float* d_wy, const float* d_wz,
 const float* d_px , const float* d_py, const float* d_pz,
 size_t nP, StreamT stream)
{
    int* i_dox = (int*) d_iwdx;
    int* i_doy = (int*) d_iwdy;
    int* i_doz = (int*) d_iwdz;
    
    // Splat to fixed buffer with the distance
    Splatting::splat3D(i_dox, i_doy, i_doz, 
		       i_dd,
		       sizeX, sizeY, sizeZ,
		       d_wx, d_wy, d_wz,
		       d_px, d_py, d_pz, nP, stream);
    //convert to the floating point buffer with weigted distance
    Splatting::convertWeightedDistance_I(d_iwdx, d_iwdy, d_iwdz, 
					 i_dd, nP, stream);
}


} // end namespace PyCA

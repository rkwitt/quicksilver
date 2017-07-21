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

#include "GFieldOperKernels.h"
#include <pycaUtils.h>
// #include <mem.h>
#include "interp.h"
#include "GSplat.h"
#include "GSplat.cuh"
#include <GMemOpers.h>

// #include "GImageFieldOpers.h"
// #include "GImageOpers.h"

// #include "FOpers.h"

#include "FiniteDiff.h"

// TEST make sure boost isn't included in nvcc code
#if defined(BOOST_COMPILER)
int bla[-1];
#endif

namespace PyCA {

////////////////////////////////////////////////////////////////////////////
// compose a velocity and hfield to get an hfield
// h(x) = g(x) + delta * v(g(x))
////////////////////////////////////////////////////////////////////////////
__device__ __constant__ float c_delta;
__device__ __constant__ float c_trans[3];

template<bool fwd, BackgroundStrategy bg>
__global__ void 
ComposeVH_kernel(float* d_hx, 
		 float* d_hy, 
		 float* d_hz,
		 const float* d_vx, 
		 const float* d_vy, 
		 const float* d_vz,
		 const float* d_gx, 
		 const float* d_gy, 
		 const float* d_gz,
		 float delta, 
		 int w, int h, int l,
		 float ispX, float ispY, float ispZ)
{
    uint i = blockIdx.x * blockDim.x + threadIdx.x;
    uint j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i < w && j < h){
        int id = i + w * j;
        for (int k=0; k < l; ++k, id+=w*h){
            
            float gx = d_gx[id];
            float gy = d_gy[id];
            float gz = d_gz[id];
            
            float vgx, vgy, vgz;
            triLerp<bg>(vgx, vgy, vgz,
                        d_vx, d_vy, d_vz,
                        gx, gy, gz, w, h, l);
            if (fwd){
                d_hx[id] = gx + delta * ispX * vgx ;
                d_hy[id] = gy + delta * ispY * vgy;
                d_hz[id] = gz + delta * ispZ * vgz;
            } else {
                d_hx[id] = gx - delta * ispX * vgx;
                d_hy[id] = gy - delta * ispY * vgy;
                d_hz[id] = gz - delta * ispZ * vgz;
            }
        }
    }
}

template<bool fwd, BackgroundStrategy bg>
__global__ void 
ComposeVH_const_kernel(float* d_hx, 
		       float* d_hy, 
		       float* d_hz,
		       const float* d_vx, 
		       const float* d_vy, 
		       const float* d_vz,
		       const float* d_gx, 
		       const float* d_gy, 
		       const float* d_gz,
		       int w, int h, int l,
		       float ispX, float ispY, float ispZ)
{
    uint i = blockIdx.x * blockDim.x + threadIdx.x;
    uint j = blockIdx.y * blockDim.y + threadIdx.y;
    float delta = c_delta;
    if (i < w && j < h){
        int id = i + w * j;
        for (int k=0; k < l; ++k, id+=w*h){
            float gx = d_gx[id];
            float gy = d_gy[id];
            float gz = d_gz[id];
            
            float vgx, vgy, vgz;
            
            triLerp<bg>(vgx, vgy, vgz,
                        d_vx, d_vy, d_vz,
                        gx, gy, gz, w, h, l);
            if (fwd){
                d_hx[id] = gx + delta * ispX * vgx ;
                d_hy[id] = gy + delta * ispY * vgy;
                d_hz[id] = gz + delta * ispZ * vgz;
            } else {
                d_hx[id] = gx - delta * ispX * vgx;
                d_hy[id] = gy - delta * ispY * vgy;
                d_hz[id] = gz - delta * ispZ * vgz;
            }
        }
    }
}

// what is this for?? jsp2014
template<bool fwd, BackgroundStrategy bg>
__global__ void ComposeVH_kernel(float* d_hx,
                                 const float* d_vx, const float* d_gx, int nAlign,
                                 float delta, int w, int h, int l,
                                 float ispX, float ispY, float ispZ)
{
    uint i = blockIdx.x * blockDim.x + threadIdx.x;
    uint j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i < w && j < h){
        int id = i + w * j;
        for (int k=0; k < l; ++k, id+=w*h){
            
            float gx = d_gx[id             ];
            float gy = d_gx[id + nAlign    ];
            float gz = d_gx[id + 2 * nAlign];
            
            float vgx, vgy, vgz;
            triLerp<bg>(vgx, vgy, vgz,
                        d_vx, d_vx + nAlign, d_vx + 2 * nAlign,
                        gx, gy, gz, w, h, l);
            if (fwd){
                d_hx[id            ] = gx + delta * ispX * vgx ;
                d_hx[id + nAlign   ] = gy + delta * ispY * vgy;
                d_hx[id + 2* nAlign] = gz + delta * ispZ * vgz;
            } else {
                d_hx[id            ] = gx - delta * ispX * vgx;
                d_hx[id + nAlign   ] = gy - delta * ispY * vgy;
                d_hx[id + 2* nAlign] = gz - delta * ispZ * vgz;
            }
        }
    }
}

template<bool fwd, BackgroundStrategy bg>
void ComposeVH(float* d_hx, 
	       float* d_hy, 
	       float* d_hz,
               const float* d_vx, 
	       const float* d_vy, 
	       const float* d_vz,
               const float* d_gx, 
	       const float* d_gy, 
	       const float* d_gz,
               const float& delta, 
	       int w, int h, int l,
               float spX, float spY, float spZ, 
	       StreamT stream, bool onDev)
{
    MK_CHECK_VFIELD_BACKGROUND(bg);
    
    dim3 threads(16,16);
    dim3 grids(iDivUp(w, threads.x), iDivUp(h, threads.y));

    if (onDev) {
        cudaMemcpyToSymbolAsync(c_delta, &delta,sizeof(float),
                                0, cudaMemcpyDeviceToDevice,stream);
        ComposeVH_const_kernel<fwd, bg><<<grids, threads, 0, stream>>>
            (d_hx, d_hy, d_hz,
             d_vx, d_vy, d_vz,
             d_gx, d_gy, d_gz,
             w, h, l, 1.f / spX, 1.f / spY, 1.f / spZ);
    } else {
        ComposeVH_kernel<fwd, bg><<<grids, threads, 0, stream>>>
            (d_hx, d_hy, d_hz,
             d_vx, d_vy, d_vz,
             d_gx, d_gy, d_gz,
             delta,  w, h, l, 1.f / spX, 1.f / spY, 1.f / spZ);
    }
}


/**
 * Compose a h field and a velocify field to get an hfield
 * h(x) = g(x+ delta * v(x))
 *
 * davisb 2007
 */
template<bool fwd, BackgroundStrategy bg>
__global__ void ComposeHV_kernel(float* d_hx, float* d_hy, float* d_hz,
                                 const float* d_gx, const float* d_gy, const float* d_gz,
                                 const float* d_vx, const float* d_vy, const float* d_vz,
                                 float delta,
                                 int w, int h, int l,
                                 float ispX, float ispY, float ispZ){
        
    uint i = blockIdx.x * blockDim.x + threadIdx.x;
    uint j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i < w && j < h){
        int id = i + w * j;
        for (int k=0; k < l; ++k, id+=w*h){
            float x,y,z;

            if (fwd){
                x = i + d_vx[id] * ispX * delta;
                y = j + d_vy[id] * ispY * delta;
                z = k + d_vz[id] * ispZ * delta;
            } else {
                x = i - d_vx[id] * ispX * delta;
                y = j - d_vy[id] * ispY * delta;
                z = k - d_vz[id] * ispZ * delta;
            }
            
            float hx, hy, hz;
            triLerp<bg>(hx, hy, hz,
                        d_gx, d_gy, d_gz,
                        x, y, z, w, h, l);
            d_hx[id] = hx;
            d_hy[id] = hy;
            d_hz[id] = hz;
        }
    }
}
template<bool fwd, BackgroundStrategy bg>
__global__ void ComposeHV_const_kernel(
    float* d_hx, float* d_hy, float* d_hz,
    const float* d_gx, const float* d_gy, const float* d_gz,
    const float* d_vx, const float* d_vy, const float* d_vz,
    int w, int h, int l,
    float ispX, float ispY, float ispZ)
{
    uint i = blockIdx.x * blockDim.x + threadIdx.x;
    uint j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i < w && j < h){
        float delta = c_delta;
        int id = i + w * j;
        for (int k=0; k < l; ++k, id+=w*h){
            float x,y,z;
            
            if (fwd){
                x = i + d_vx[id] * ispX * delta;
                y = j + d_vy[id] * ispY * delta;
                z = k + d_vz[id] * ispZ * delta;
            } else {
                x = i - d_vx[id] * ispX * delta;
                y = j - d_vy[id] * ispY * delta;
                z = k - d_vz[id] * ispZ * delta;
            }
            
            float hx, hy, hz;
            triLerp<bg>(hx, hy, hz,
                        d_gx, d_gy, d_gz,
                        x, y, z, w, h, l);
            d_hx[id] = hx;
            d_hy[id] = hy;
            d_hz[id] = hz;
        }
    }
}

template<bool fwd, BackgroundStrategy bg>
void ComposeHV(float* d_hx, 
	       float* d_hy, 
	       float* d_hz,
               const float* d_gx, 
	       const float* d_gy, 
	       const float* d_gz,
               const float* d_vx, 
	       const float* d_vy, 
	       const float* d_vz,
               const float& delta,
               int w, int h, int l,
               float spX, float spY, float spZ, 
	       StreamT stream, bool onDev)
{
    MK_CHECK_HFIELD_BACKGROUND(bg);
    
    dim3 threads(16,16);
    dim3 grids(iDivUp(w, threads.x), iDivUp(h, threads.y));
    if (onDev) {
        cudaMemcpyToSymbolAsync(c_delta, &delta,sizeof(float),
                                0, cudaMemcpyDeviceToDevice,stream);
        ComposeHV_const_kernel<fwd, bg><<<grids, threads,0,stream>>>
            (d_hx, d_hy, d_hz,
             d_gx, d_gy, d_gz,
             d_vx, d_vy, d_vz,
             w, h, l,
             1.f / spX, 1.f / spY, 1.f / spZ);

    } else {
        ComposeHV_kernel<fwd, bg><<<grids, threads,0,stream>>>
            (d_hx, d_hy, d_hz,
             d_gx, d_gy, d_gz,
             d_vx, d_vy, d_vz,
             delta,
             w, h, l,
             1.f / spX, 1.f / spY, 1.f / spZ);
    }
}

template<BackgroundStrategy bg>
__global__ void ComposeTranslation_kernel(
    float* d_ox, float* d_oy, float* d_oz,
    const float* d_ix, const float* d_iy, const float* d_iz,
    const float tx, const float ty, const float tz,
    int sizeX, int sizeY, int sizeZ){
        
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    
    if (i < sizeX && j < sizeY){
        int id = j * sizeX + i;
        for (int k=0; k< sizeZ; ++k, id+= sizeX*sizeY){
            
            float x = i + tx;
            float y = j + ty;
            float z = k + tz;

            float ox, oy, oz;
            triLerp<bg>(ox, oy, oz,
                        d_ix, d_iy, d_iz,
                        x, y, z,
                        sizeX, sizeY, sizeZ);

            d_ox[id] = ox;
            d_oy[id] = oy;
            d_oz[id] = oz;
        }
    }
}

template<BackgroundStrategy bg>
__global__ void ComposeTranslation_const_kernel(
    float* d_ox, float* d_oy, float* d_oz,
    const float* d_ix, const float* d_iy, const float* d_iz,
    int sizeX, int sizeY, int sizeZ){
        
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    
    if (i < sizeX && j < sizeY){
        int id = j * sizeX + i;
        for (int k=0; k< sizeZ; ++k, id+= sizeX*sizeY){
            
            float x = i + c_trans[0];
            float y = j + c_trans[1];
            float z = k + c_trans[2];

            float ox, oy, oz;
            triLerp<bg>(ox, oy, oz,
                        d_ix, d_iy, d_iz,
                        x, y, z,
                        sizeX, sizeY, sizeZ);

            d_ox[id] = ox;
            d_oy[id] = oy;
            d_oz[id] = oz;
        }
    }
}

template<BackgroundStrategy bg>
void
ComposeTranslation(float *d_ox, 
		   float *d_oy, 
		   float *d_oz, 
		   const float *d_ix, 
		   const float *d_iy, 
		   const float *d_iz, 
		   const Vec3Di& sz,
		   const Vec3Df& t, 
		   StreamT stream, 
		   bool onDev) 
{
    dim3 threads(16,16);
    dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));

    if (onDev) {
        cudaMemcpyToSymbolAsync(c_trans, &t.x,sizeof(float) * 3,
                                0, cudaMemcpyDeviceToDevice,stream);
        ComposeTranslation_const_kernel<bg><<<grids, threads, 0, stream>>>
            (d_ox, d_oy, d_oz,
             d_ix, d_iy, d_iz,
             sz.x, sz.y, sz.z);
    } else {
        ComposeTranslation_kernel<bg><<<grids, threads, 0, stream>>>
            (d_ox, d_oy, d_oz,
             d_ix, d_iy, d_iz,
             t.x, t.y, t.z, 
             sz.x, sz.y, sz.z);
    }
}

template<BackgroundStrategy bg>
__global__ void ApplyH_kernel(float* d_hx, float* d_hy, float* d_hz,
                              const float* d_fx, const float* d_fy, const float* d_fz,
                              const float* d_gx, const float* d_gy, const float* d_gz,
                              int w, int h, int l){

    uint i = blockIdx.x * blockDim.x + threadIdx.x;
    uint j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i < w && j < h){
        int id = i + w * j;
        for (int k=0; k < l; ++k, id+=w*h){
            
            float x = d_gx[id];
            float y = d_gy[id];
            float z = d_gz[id];
            
            float hx, hy, hz;
                
            triLerp<bg>(hx, hy, hz,
                        d_fx, d_fy, d_fz,
                        x, y, z,  w, h, l);
            
            d_hx[id] = hx;
            d_hy[id] = hy;
            d_hz[id] = hz;
        }
    }
}


template<BackgroundStrategy bg>
void
ApplyH(float *d_ox, 
       float *d_oy, 
       float *d_oz, 
       const float *d_ix, 
       const float *d_iy, 
       const float *d_iz, 
       const float *d_hx, 
       const float *d_hy, 
       const float *d_hz, 
       const Vec3Di &sz,
       StreamT stream)
{
    dim3 threads(16,16);
    dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));

    ApplyH_kernel<bg><<<grids, threads, 0, stream>>>(
        d_ox, d_oy, d_oz,
        d_ix, d_iy, d_iz,
        d_hx, d_hy, d_hz, 
	sz.x, sz.y, sz.z);
}


template<bool fwd, BackgroundStrategy bg>
__global__ void ApplyV_kernel(float* d_hx, float* d_hy, float* d_hz,
                              const float* d_fx, const float* d_fy, const float* d_fz,
                              const float* d_ux, const float* d_uy, const float* d_uz,
                              float delta, int w, int h, int l,
                              float iSpX, float iSpY, float iSpZ){

    uint i = blockIdx.x * blockDim.x + threadIdx.x;
    uint j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i < w && j < h){
        int id = i + w * j;
        for (int k=0; k < l; ++k, id+=w*h){

            float x, y, z;

            if (fwd) {
                x = i + delta * iSpX * d_ux[id];
                y = j + delta * iSpY * d_uy[id];
                z = k + delta * iSpZ * d_uz[id];
            } else {
                x = i - delta * iSpX * d_ux[id];
                y = j - delta * iSpY * d_uy[id];
                z = k - delta * iSpZ * d_uz[id];
            }

            float hx, hy, hz;
            triLerp<bg>(hx, hy, hz,
                        d_fx, d_fy, d_fz,
                        x, y, z,  w, h, l);

            d_hx[id] = hx;
            d_hy[id] = hy;
            d_hz[id] = hz;
        }
    }
}

template<bool fwd, BackgroundStrategy bg>
__global__ void ApplyV_const_kernel(float* d_hx, float* d_hy, float* d_hz,
                                    const float* d_fx, const float* d_fy, const float* d_fz,
                                    const float* d_ux, const float* d_uy, const float* d_uz,
                                    int w, int h, int l,
                                    float iSpX, float iSpY, float iSpZ){

    uint i = blockIdx.x * blockDim.x + threadIdx.x;
    uint j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i < w && j < h){
        float delta = c_delta;
        int id = i + w * j;
        for (int k=0; k < l; ++k, id+=w*h){

            float x, y, z;
            if (fwd) {
                x = i + delta * iSpX * d_ux[id];
                y = j + delta * iSpY * d_uy[id];
                z = k + delta * iSpZ * d_uz[id];
            } else {
                x = i - delta * iSpX * d_ux[id];
                y = j - delta * iSpY * d_uy[id];
                z = k - delta * iSpZ * d_uz[id];
            }

            float hx, hy, hz;
            triLerp<bg>(hx, hy, hz,
                        d_fx, d_fy, d_fz,
                        x, y, z,  w, h, l);

            d_hx[id] = hx;
            d_hy[id] = hy;
            d_hz[id] = hz;
        }
    }
}

template<bool fwd, BackgroundStrategy bg>
void ApplyV(float *d_ox, 
	    float *d_oy, 
	    float *d_oz, 
	    const float *d_ix, 
	    const float *d_iy, 
	    const float *d_iz, 
	    const float *d_ux, 
	    const float *d_uy, 
	    const float *d_uz, 
	    const Vec3Di &sz,
	    const Vec3Df &sp,
	    const float& delta, 
	    StreamT stream, bool onDev)
{
    Vec3Df iSp  = Vec3Df(1.f / sp.x, 1.f / sp.y, 1.f / sp.z);
    
    dim3 threads(16,16);
    dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));

    if (onDev) {
        cudaMemcpyToSymbolAsync(c_delta, &delta,sizeof(float),
                                0, cudaMemcpyDeviceToDevice,stream);

        ApplyV_const_kernel<fwd, bg><<<grids, threads, 0, stream>>>(
            d_ox, d_oy, d_oz,
            d_ix, d_iy, d_iz,
            d_ux, d_uy, d_uz, 
            sz.x, sz.y, sz.z,
            iSp.x, iSp.y, iSp.z);
    } else {
        ApplyV_kernel<fwd, bg><<<grids, threads, 0, stream>>>(
            d_ox, d_oy, d_oz,
            d_ix, d_iy, d_iz,
            d_ux, d_uy, d_uz, delta,
            sz.x, sz.y, sz.z,
            iSp.x, iSp.y, iSp.z);
    }
}

template<BackgroundStrategy bg, bool rescale>
__global__ void Resampling_kernel(float* d_ox, float* d_oy, float* d_oz,
                                  const float* d_ix, const float* d_iy, const float* d_iz,
                                  int osizeX, int osizeY, int osizeZ,
                                  int isizeX, int isizeY, int isizeZ)
{
    uint x = blockIdx.x * blockDim.x + threadIdx.x;
    uint y = blockIdx.y * blockDim.y + threadIdx.y;
    
    float rX = (float)isizeX / (float)osizeX;
    float rY = (float)isizeY / (float)osizeY;
    float rZ = (float)isizeZ / (float)osizeZ;

    if (x < osizeX && y < osizeY){
        int id = x + osizeX * y;
        
        float i_x =  (rX - 1.f) / 2.f + x * rX;
        float i_y =  (rY - 1.f) / 2.f + y * rY;
        
        for (int z=0; z < osizeZ; ++z, id+=osizeX * osizeY){
            float i_z =  (rZ - 1.f) / 2.f + z * rZ;

            float ox, oy, oz;
            triLerp<bg>(ox, oy, oz,
                        d_ix, d_iy, d_iz,
                        i_x, i_y, i_z,
                        isizeX, isizeY, isizeZ);
                
            if (rescale){
                ox /= rX; oy /= rY; oz /= rZ;
            }
                
            d_ox[id] = ox;
            d_oy[id] = oy;
            d_oz[id] = oz;
        }
    }
}

template<BackgroundStrategy bg,  bool rescaleVector>
void
Resample(float *d_ox, 
	 float *d_oy, 
	 float *d_oz, 
	 const Vec3Di &oSz,
	 const float *d_ix, 
	 const float *d_iy, 
	 const float *d_iz, 
	 const Vec3Di &iSz,
	 StreamT stream)
{
    dim3 threads(16,16);
    dim3 grids(iDivUp(oSz.x, threads.x), iDivUp(oSz.y, threads.y));

    Resampling_kernel<bg, rescaleVector><<<grids, threads, 0, stream>>>
        (d_ox, d_oy, d_oz,
         d_ix, d_iy, d_iz,
         oSz.x, oSz.y, oSz.z,
         iSz.x, iSz.y, iSz.z);

}


__global__ void ReprojectToUnitVec_kernel
(float* d_ox, float* d_oy, float *d_oz,
 int szX, int szY, int szZ)
{
        
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    
    if (i < szX && j < szY){
       int id = j * szX + i;
       
       for (int k=0; k< szZ; ++k, id+= szX*szY){
	  
	  float vx = d_ox[id];
	  float vy = d_oy[id];
	  float vz = d_oz[id];
	  
	  float l = sqrt(vx*vx+vy*vy+vz*vz);
	  
	  if(l>1.0){
	     d_ox[id] = vx/l;
	     d_oy[id] = vy/l;
	     d_oz[id] = vz/l;
	  }
            
        }
    }
}

void
ReprojectToUnitVec(float *d_ox, 
		   float *d_oy, 
		   float *d_oz, 
		   const Vec3Di &sz,
		   StreamT st)
{
   dim3 threads(16,16);
   dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));

   ReprojectToUnitVec_kernel<<<grids, threads, 0, st>>>
      (d_ox, d_oy, d_oz,
       sz.x, sz.y, sz.z);
}

__global__ void NormalizeSafe_kernel
(float* d_ox, float* d_oy, float *d_oz,
 const float* d_ix, const float* d_iy, const float* d_iz,
 float eps,
 int szX, int szY, int szZ)
{
        
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    
    if (i < szX && j < szY){
       int id = j * szX + i;
       
       for (int k=0; k< szZ; ++k, id+= szX*szY){
	  
	  float vx = d_ix[id];
	  float vy = d_iy[id];
	  float vz = d_iz[id];
	  
	  float l = sqrt(vx*vx+vy*vy+vz*vz);
	  
	  if(l>eps){
	     d_ox[id] = vx/l;
	     d_oy[id] = vy/l;
	     d_oz[id] = vz/l;
	  }else{
	     d_ox[id] = 0;
	     d_oy[id] = 0;
	     d_oz[id] = 0;
	  }
            
        }
    }
}

void
NormalizeSafe(float *d_ox, 
	      float *d_oy, 
	      float *d_oz, 
	      const float *d_ix, 
	      const float *d_iy, 
	      const float *d_iz, 
	      const Vec3Di &sz,
	      const float& eps, 
	      StreamT st)
{

   dim3 threads(16,16);
   dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));

   NormalizeSafe_kernel<<<grids, threads, 0, st>>>
      (d_ox, d_oy, d_oz,
       d_ix, d_iy, d_iz,
       eps,
       sz.x, sz.y, sz.z);
}

__global__ void Shrink_kernel
(float* d_ox, float* d_oy, float *d_oz,
 const float* d_ix, const float* d_iy, const float* d_iz,
 float eps,
 int szX, int szY, int szZ)
{
        
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    
    if (i < szX && j < szY){
       int id = j * szX + i;
       
       for (int k=0; k< szZ; ++k, id+= szX*szY){
	  
	  float vx = d_ix[id];
	  float vy = d_iy[id];
	  float vz = d_iz[id];
	  
	  float l = sqrt(vx*vx+vy*vy+vz*vz);
	  float shrink = (l-eps)/l;
	  
	  if(l>eps){
	     d_ox[id] = vx*shrink;
	     d_oy[id] = vy*shrink;
	     d_oz[id] = vz*shrink;
	  }else{
	     d_ox[id] = 0;
	     d_oy[id] = 0;
	     d_oz[id] = 0;
	  }
            
        }
    }
}

void
Shrink(float *d_ox, 
       float *d_oy, 
       float *d_oz, 
       const float *d_ix, 
       const float *d_iy, 
       const float *d_iz, 
       const Vec3Di &sz,
       const float& eps, 
       StreamT st)
{

   dim3 threads(16,16);
   dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));

   Shrink_kernel<<<grids, threads, 0, st>>>
      (d_ox, d_oy, d_oz,
       d_ix, d_iy, d_iz,
       eps,
       sz.x, sz.y, sz.z);
}

template<BackgroundStrategy bg>
__global__ void FixedPointInverse_kernel
(float* d_ginvx, float* d_ginvy, float *d_ginvz,
 const float* d_gx, const float* d_gy, const float *d_gz,
 int szX, int szY, int szZ, unsigned int numIter)
{ // This does the fixed point iteration:
    // g^{-1}_{k+1}(x) = g^{-1}_{k}(x) + (x-g \circ g^{-1}_{k}(x))
    int x = threadIdx.x + blockIdx.x * blockDim.x;
    int y = threadIdx.y + blockIdx.y * blockDim.y;

    if (x < szX && y < szY)
    {
        int id = y * szX + x;

        float ghx, ghy, ghz; // interpolate into g at various places
        for (int z=0; z < szZ; ++z, id+= szX*szY)
        {
            float hx = d_ginvx[id]; // start with given initial estimate ginv(x)
            float hy = d_ginvy[id];
            float hz = d_ginvz[id];
            // this will be the output
            for (unsigned int iter=0;iter < numIter;++iter)
            {
                /*if (hx < 0 || hx > szX-1 ||*/
                    /*hy < 0 || hy > szY-1 ||*/
                    /*hz < 0 || hz > szZ-1)*/
                    /*break;*/
                triLerp<bg>(ghx, ghy, ghz,
                            d_gx, d_gy, d_gz,
                            hx, hy, hz,
                            szX, szY, szZ);

                hx += ((float)x - ghx);
                hy += ((float)y - ghy);
                hz += ((float)z - ghz);
            }
            // set output
            d_ginvx[id] = hx;
            d_ginvy[id] = hy;
            d_ginvz[id] = hz;
        }
    }
}

template<BackgroundStrategy bg> 
void
FixedPointInverse(float *ginvx, 
		  float *ginvy, 
		  float *ginvz, 
		  const float *gx, 
		  const float *gy, 
		  const float *gz, 
		  const Vec3Di &sz,
		  unsigned int numIter, 
		  StreamT stream, bool onDev)
{
    dim3 threads(16,16);
    dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));

    FixedPointInverse_kernel<bg><<<grids, threads, 0, stream>>>
           (ginvx, ginvy, ginvz,
            gx, gy, gz,
            sz.x, sz.y, sz.z,
            numIter);
}

__global__ void updateInverseSubFromIndentity_kernel
(float* d_hx, float* d_hy, float* d_hz,
 int w, int h, int l,
 float ispX, float ispY, float ispZ){
    uint i     = blockIdx.x * blockDim.x + threadIdx.x;
    uint j     = blockIdx.y * blockDim.y + threadIdx.y;
    uint index = j * w + i;
    
    if (i < w && j < h){
        for (int k=0; k<l; ++k, index+=w*h){
            d_hx[index] = float(i) -  ispX * d_hx[index];
            d_hy[index] = float(j) -  ispY * d_hy[index];
            d_hz[index] = float(k) -  ispZ * d_hz[index];
        }
    }
}

void
updateInverseSubFromIndentity(float *d_hx,
			      float *d_hy,
			      float *d_hz,
			      const Vec3Di &sz,
			      const Vec3Df &sp,
			      StreamT stream)
{
    Vec3Df iSp  = Vec3Df(1.f / sp.x, 1.f / sp.y, 1.f / sp.z);
    
    dim3 threads(16,16);
    dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));
    
    updateInverseSubFromIndentity_kernel<<<grids, threads, 0, stream>>>
	(d_hx, d_hy, d_hz,
	 sz.x, sz.y, sz.z,
	 iSp.x, iSp.y, iSp.z);
}

// Lie algebra methods

/*
 * Adjoint action of Diff on its Lie algebra
 * This is just the pushforward
 * Z = Ad_g X = |Dg|\circ g^{-1} X\circ g^{-1}
 */
template<BackgroundStrategy bg>
__global__ void Ad_kernel
(int* d_Zx, int* d_Zy, int *d_Zz,
 const float* d_gx, const float* d_gy, const float *d_gz,
 const float* d_Xx, const float* d_Xy, const float *d_Xz,
 float scalex, float scaley, float scalez,
 int szX, int szY, int szZ)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;

    float Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz;
    float dgvx, dgvy, dgvz;

    if (i < szX && j < szY)
    {
        int id = j * szX + i;

        for (int k=0; k< szZ; ++k, id+= szX*szY)
        {
            // Get Jacobian matrix
            jacobianPoint<float,DIFF_CENTRAL,BC_CLAMP>(d_gx,d_gy,d_gz,i,j,k,Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz,szX,szY,szZ);

            if (szZ == 1)
            { // Special case for flat images
                Jxz = Jyz = Jzx = Jzy = 0;
                Jzz = 1;
            }

            // Compute determinant
            float det = Jxx*(Jyy*Jzz-Jyz*Jzy)
                       -Jxy*(Jyx*Jzz-Jyz*Jzx)
                       +Jxz*(Jyx*Jzy-Jyy*Jzx);

	    // Multiply by det*Dg X
	    dgvx = det*(Jxx*d_Xx[id] + Jxy*d_Xy[id] + Jxz*d_Xz[id]);
	    dgvy = det*(Jyx*d_Xx[id] + Jyy*d_Xy[id] + Jyz*d_Xz[id]);
	    dgvz = det*(Jzx*d_Xx[id] + Jzy*d_Xy[id] + Jzz*d_Xz[id]);
	    
	    // Splat each component (non-normalized)
	    Splatting::atomicSplat(d_Zx, d_Zy, d_Zz, 
				  dgvx, dgvy, dgvz, 
				  d_gx[id], d_gy[id], d_gz[id],
				  szX, szY, szZ);
        }
    }
}

template<BackgroundStrategy bg> 
void
Ad(float *Zx, 
   float *Zy, 
   float *Zz, 
   const float *gx, 
   const float *gy, 
   const float *gz, 
   const float *Xx,
   const float *Xy,
   const float *Xz,
   const Vec3Di &sz,
   const Vec3Df &sp,
   StreamT s,bool onDev)
{
    size_t nVox = sz.x * sz.y * sz.z;

    dim3 threads(16,16);
    dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));

    int* i_Zx =(int*)Zx;
    int* i_Zy =(int*)Zy;
    int* i_Zz =(int*)Zz;
    GMemOpers<int>::SetMem(i_Zx, 0, nVox, s, false);
    GMemOpers<int>::SetMem(i_Zy, 0, nVox, s, false);
    GMemOpers<int>::SetMem(i_Zz, 0, nVox, s, false);

    Ad_kernel<bg><<<grids, threads, 0, s>>>
      (i_Zx, i_Zy, i_Zz,
       gx, gy, gz,
       Xx, Xy, Xz,
       1.f/sp.x, 1.f/sp.y, 1.f/sp.z, // scaling for Jacobian
       sz.x, sz.y, sz.z);

    Splatting::FixedToFloating_I(Zx, Zy, Zz, nVox, s);
}

/*
 * infinitesimal adjoint action
 * Z = ad_X Y = DX Y - DY X
 */
__global__ void ad_kernel
(float* d_Zx, float* d_Zy, float *d_Zz,
 const float* d_Xx, const float* d_Xy, const float *d_Xz,
 const float* d_Yx, const float* d_Yy, const float *d_Yz,
 float scalex, float scaley, float scalez,
 int szX, int szY, int szZ)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;

    float Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz;
    
    if (i < szX && j < szY)
    {
        int id = j * szX + i;

        for (int k=0; k< szZ; ++k, id+= szX*szY)
        {
            // Get DX
	   jacobianPoint<float,DIFF_CENTRAL,BC_CLAMP>(d_Xx,d_Xy,d_Xz,i,j,k,Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz,szX,szY,szZ);

            // Start with DX Y
            d_Zx[id] = Jxx*d_Yx[id] + Jxy*d_Yy[id] + Jxz*d_Yz[id];
            d_Zy[id] = Jyx*d_Yx[id] + Jyy*d_Yy[id] + Jyz*d_Yz[id];
            d_Zz[id] = Jzx*d_Yx[id] + Jzy*d_Yy[id] + Jzz*d_Yz[id];

            // Get DY
            jacobianPoint<float,DIFF_CENTRAL,BC_CLAMP>(d_Yx,d_Yy,d_Yz,i,j,k,Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz,szX,szY,szZ);

            // Subtract DY X
            d_Zx[id] -= Jxx*d_Xx[id] + Jxy*d_Xy[id] + Jxz*d_Xz[id];
            d_Zy[id] -= Jyx*d_Xx[id] + Jyy*d_Xy[id] + Jyz*d_Xz[id];
            d_Zz[id] -= Jzx*d_Xx[id] + Jzy*d_Xy[id] + Jzz*d_Xz[id];
        }
    }
}

void
AdInf(float *Zx, 
      float *Zy, 
      float *Zz, 
      const float *Xx, 
      const float *Xy, 
      const float *Xz, 
      const float *Yx,
      const float *Yy,
      const float *Yz,
      const Vec3Di &sz,
      const Vec3Df &sp,
      StreamT s,bool onDev)
{
   dim3 threads(16,16);
   dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));

   ad_kernel<<<grids, threads, 0, s>>>
      (Zx, Zy, Zz,
       Xx, Xy, Xz,
       Yx, Yy, Yz,
       1.f/sp.x, 1.f/sp.y, 1.f/sp.z, // scaling for Jacobian
       sz.x, sz.y, sz.z);
}

/*
 * Jacobian X times Y
 * Z = DX Y
 */
__global__ void jacobianXY_kernel
(float* d_Zx, float* d_Zy, float *d_Zz,
 const float* d_Xx, const float* d_Xy, const float *d_Xz,
 const float* d_Yx, const float* d_Yy, const float *d_Yz,
 float scalex, float scaley, float scalez,
 int szX, int szY, int szZ)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;

    float Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz;
    
    if (i < szX && j < szY)
    {
        int id = j * szX + i;

        for (int k=0; k< szZ; ++k, id+= szX*szY)
        {
            // Get DX
	   jacobianPoint<float,DIFF_CENTRAL,BC_CLAMP>(d_Xx,d_Xy,d_Xz,i,j,k,Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz,szX,szY,szZ);

            // Compute DX Y
            d_Zx[id] = Jxx*d_Yx[id] + Jxy*d_Yy[id] + Jxz*d_Yz[id];
            d_Zy[id] = Jyx*d_Yx[id] + Jyy*d_Yy[id] + Jyz*d_Yz[id];
            d_Zz[id] = Jzx*d_Yx[id] + Jzy*d_Yy[id] + Jzz*d_Yz[id];
        }
    }
}

void
JacobianXY(float *Zx, 
	   float *Zy, 
	   float *Zz, 
	   const float *Xx, 
	   const float *Xy, 
	   const float *Xz, 
	   const float *Yx,
	   const float *Yy,
	   const float *Yz,
	   const Vec3Di &sz,
	   const Vec3Df &sp,
	   StreamT s, bool onDev)
{

   dim3 threads(16,16);
   dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));

   jacobianXY_kernel<<<grids, threads, 0, s>>>
      (Zx, Zy, Zz,
       Xx, Xy, Xz,
       Yx, Yy, Yz,
       1.f/sp.x, 1.f/sp.y, 1.f/sp.z, // scaling for Jacobian
       sz.x, sz.y, sz.z);
}

/*
 * Jacobian X transpose times Y
 * Z = (DX)' Y
 */
__global__ void jacobianXtY_kernel
(float* d_Zx, float* d_Zy, float *d_Zz,
 const float* d_Xx, const float* d_Xy, const float *d_Xz,
 const float* d_Yx, const float* d_Yy, const float *d_Yz,
 float scalex, float scaley, float scalez,
 int szX, int szY, int szZ)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;

    float Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz;
    
    if (i < szX && j < szY)
    {
        int id = j * szX + i;

        for (int k=0; k< szZ; ++k, id+= szX*szY)
        {
            // Get DX
	   jacobianPoint<float,DIFF_CENTRAL,BC_CLAMP>(d_Xx,d_Xy,d_Xz,i,j,k,Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz,szX,szY,szZ);

            // Compute (DX)' Y
            d_Zx[id] = Jxx*d_Yx[id] + Jyx*d_Yy[id] + Jzx*d_Yz[id];
            d_Zy[id] = Jxy*d_Yx[id] + Jyy*d_Yy[id] + Jzy*d_Yz[id];
            d_Zz[id] = Jxz*d_Yx[id] + Jyz*d_Yy[id] + Jzz*d_Yz[id];
        }
    }
}

void
JacobianXtY(float *Zx, 
	   float *Zy, 
	   float *Zz, 
	   const float *Xx, 
	   const float *Xy, 
	   const float *Xz, 
	   const float *Yx,
	   const float *Yy,
	   const float *Yz,
	   const Vec3Di &sz,
	   const Vec3Df &sp,
	   StreamT s, bool onDev)
{

   dim3 threads(16,16);
   dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));

   jacobianXtY_kernel<<<grids, threads, 0, s>>>
      (Zx, Zy, Zz,
       Xx, Xy, Xz,
       Yx, Yy, Yz,
       1.f/sp.x, 1.f/sp.y, 1.f/sp.z, // scaling for Jacobian
       sz.x, sz.y, sz.z);
}

/*
 * Coadjoint action of Diff on its Lie algebra
 * n = Ad_g^* m = (Dg)^T m\circ g |Dg|
 */
template<BackgroundStrategy bg>
__global__ void CoAd_kernel
(float* d_nx, float* d_ny, float *d_nz,
 const float* d_gx, const float* d_gy, const float *d_gz,
 const float* d_mx, const float* d_my, const float *d_mz,
 float scalex, float scaley, float scalez,
 int szX, int szY, int szZ)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;

    float Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz;
    float mgx, mgy, mgz;
    
    if (i < szX && j < szY)
    {
        int id = j * szX + i;

        for (int k=0; k< szZ; ++k, id+= szX*szY)
        {
            float gx = d_gx[id];
            float gy = d_gy[id];
            float gz = d_gz[id];

            // Get Jacobian matrix
            jacobianPoint<float,DIFF_CENTRAL,BC_CLAMP>(d_gx,d_gy,d_gz,i,j,k,Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz,szX,szY,szZ);

            if (szZ == 1)
            { // Special case for flat images
                Jxz = Jyz = Jzx = Jzy = 0;
                Jzz = 1;
            }

            // Compute determinant
            float det = Jxx*(Jyy*Jzz-Jyz*Jzy)
                       -Jxy*(Jyx*Jzz-Jyz*Jzx)
                       +Jxz*(Jyx*Jzy-Jyy*Jzx);

            // Interpolate m
            triLerp<bg>(mgx, mgy, mgz,
                        d_mx, d_my, d_mz,
                        gx, gy, gz, szX, szY, szZ);

            // Multiply by det*Dg^T
            d_nx[id] = det*(Jxx*mgx + Jyx*mgy + Jzx*mgz);
            d_ny[id] = det*(Jxy*mgx + Jyy*mgy + Jzy*mgz);
            d_nz[id] = det*(Jxz*mgx + Jyz*mgy + Jzz*mgz);
        }
    }
}

template<BackgroundStrategy bg> 
void
CoAd(float *nx, 
     float *ny, 
     float *nz, 
     const float *gx, 
     const float *gy, 
     const float *gz, 
     const float *mx,
     const float *my,
     const float *mz,
     const Vec3Di &sz,
     const Vec3Df &sp,
     StreamT s,bool onDev)
{
    dim3 threads(16,16);
    dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));

    CoAd_kernel<bg><<<grids, threads, 0, s>>>
	(nx, ny, nz,
	 gx, gy, gz,
	 mx, my, mz,
	 1.f/sp.x, 1.f/sp.y, 1.f/sp.z, // scaling for Jacobian
	 sz.x, sz.y, sz.z);
}

/*
 * infinitesimal coadjoint action
 * n = ad_X^* m = (DX)^T m + div(m \otimes X)
 */
__global__ void CoAdInf_kernel
(float* d_nx, float* d_ny, float *d_nz,
 const float* d_Xx, const float* d_Xy, const float *d_Xz,
 const float* d_mx, const float* d_my, const float *d_mz,
 float scalex, float scaley, float scalez,
 int szX, int szY, int szZ)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;

    float Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz;
    
    if (i < szX && j < szY)
    {
        int id = j * szX + i;

        for (int k=0; k< szZ; ++k, id+= szX*szY)
        {
            // Start with the tensor product divergence piece
            divtensorprodPoint<float,DIFF_CENTRAL,BC_CLAMP>(d_mx,d_my,d_mz,d_Xx,d_Xy,d_Xz,i,j,k,d_nx[id],d_ny[id],d_nz[id],szX,szY,szZ);

            // Get DX
            jacobianPoint<float,DIFF_CENTRAL,BC_CLAMP>(d_Xx,d_Xy,d_Xz,i,j,k,Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz,szX,szY,szZ);

            // Add the DX^T m term
            d_nx[id] += Jxx*d_mx[id] + Jyx*d_my[id] + Jzx*d_mz[id];
            d_ny[id] += Jxy*d_mx[id] + Jyy*d_my[id] + Jzy*d_mz[id];
            d_nz[id] += Jxz*d_mx[id] + Jyz*d_my[id] + Jzz*d_mz[id];
        }
    }
}

void
CoAdInf(float *nx, 
	float *ny, 
	float *nz, 
	const float *Xx, 
	const float *Xy, 
	const float *Xz, 
	const float *mx,
	const float *my,
	const float *mz,
	const Vec3Di &sz,
	const Vec3Df &sp,
	StreamT s,bool onDev)
{
   dim3 threads(16,16);
   dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));

   CoAdInf_kernel<<<grids, threads, 0, s>>>
       (nx, ny, nz,
	Xx, Xy, Xz,
	mx, my, mz,
	1.f/sp.x, 1.f/sp.y, 1.f/sp.z, // scaling for Jacobian
	sz.x, sz.y, sz.z);
}


/*
 * computes tensor divergence of outer product of two vector fields 
 * Z = div(X \otimes Y)
 */
__global__ void DivergenceTensor_kernel
(float* d_Zx, float* d_Zy, float *d_Zz,
 const float* d_Xx, const float* d_Xy, const float *d_Xz,
 const float* d_Yx, const float* d_Yy, const float *d_Yz,
 int szX, int szY, int szZ)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
   
    if (i < szX && j < szY)
    {
        int id = j * szX + i;

        for (int k=0; k< szZ; ++k, id+= szX*szY)
        {
            // The tensor product divergence 
            divtensorprodPoint<float,DIFF_CENTRAL,BC_CLAMP>(d_Xx,d_Xy,d_Xz,d_Yx,d_Yy,d_Yz,i,j,k,d_Zx[id],d_Zy[id],d_Zz[id],szX,szY,szZ);
        }
    }
}

void
DivergenceTensor(float *Zx, 
		 float *Zy, 
		 float *Zz, 
		 const float *Xx, 
		 const float *Xy, 
		 const float *Xz, 
		 const float *Yx,
		 const float *Yy,
		 const float *Yz,
		 const Vec3Di &sz,
		 const Vec3Df &sp,
		 StreamT s,bool onDev)
{
   dim3 threads(16,16);
   dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));

   DivergenceTensor_kernel<<<grids, threads, 0, s>>>
       (Zx, Zy, Zz,
	Xx, Xy, Xz,
	Yx, Yy, Yz,
	sz.x, sz.y, sz.z);
}


//Instantiation
#include "GFieldOperKernels_inst.cxx"

} // end namespace PyCA

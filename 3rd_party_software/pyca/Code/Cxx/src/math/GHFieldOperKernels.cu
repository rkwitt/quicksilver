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

#include "GHFieldOperKernels.h"
#include <pycaUtils.h>

// TEST make sure boost isn't included in nvcc code
#if defined(BOOST_COMPILER)
int bla[-1];
#endif

namespace PyCA {

//////////////////////////////////////////////////////////////////////// 
// set hfield to identity
// i.e. h(x) = x
////////////////////////////////////////////////////////////////////////
__global__ void SetToIdentity_kernel(float* d_hx, float* d_hy, float* d_hz,
                                     int w, int h, int l){
    uint i     = blockIdx.x * blockDim.x + threadIdx.x;
    uint j     = blockIdx.y * blockDim.y + threadIdx.y;
    uint index = j * w + i;
    if (i < w && j < h){
        for (int k=0; k<l; ++k, index+=w*h){
            d_hx[index] = i;
            d_hy[index] = j;
            d_hz[index] = k;
        }
    }
}

void
SetToIdentity(float *d_hx, 
	      float *d_hy, 
	      float *d_hz, 
	      const Vec3Di &sz,
	      StreamT stream) 
{
    dim3 threads(16,16);
    dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));
    SetToIdentity_kernel<<<grids, threads, 0, stream>>>
	(d_hx, d_hy, d_hz,
	 sz.x, sz.y, sz.z);
}

//////////////////////////////////////////////////////////////////////// 
// convert hfield to velocity field
// i.e. v(x) = (h(x) - x) * sp
////////////////////////////////////////////////////////////////////////
__global__ void hToVelocity_kernel(float* d_vx, float* d_vy, float* d_vz,
                                   const float* d_hx, const float* d_hy, const float* d_hz,
                                   int w, int h, int l,
                                   float spX, float spY, float spZ)
{
    uint i     = blockIdx.x * blockDim.x + threadIdx.x;
    uint j     = blockIdx.y * blockDim.y + threadIdx.y;
    uint index = j * w + i;
    if (i < w && j < h){
        for (int k=0; k<l; ++k, index+=w*h){
            d_vx[index] = (d_hx[index] - i) * spX;
            d_vy[index] = (d_hy[index] - j) * spY;
            d_vz[index] = (d_hz[index] - k) * spZ;
        }
    }
}


void
toV(float *vx, 
    float *vy, 
    float *vz, 
    const float *hx, 
    const float *hy, 
    const float *hz, 
    const Vec3Di &sz,
    const Vec3Df &sp,
    StreamT stream) 
{
    dim3 threads(16,16);
    dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));
    hToVelocity_kernel<<<grids, threads, 0, stream>>>
	(vx, vy, vz,
	 hx, hy, hz,
	 sz.x, sz.y, sz.z,
	 sp.x, sp.y, sp.z);
}

__global__ void hToVelocity_I_kernel(float* d_vx, float* d_vy, float* d_vz,
                                     int w, int h, int l,
                                     float spX, float spY, float spZ)
{
    uint i     = blockIdx.x * blockDim.x + threadIdx.x;
    uint j     = blockIdx.y * blockDim.y + threadIdx.y;
    uint index = j * w + i;
    if (i < w && j < h){
        for (int k=0; k<l; ++k, index+=w*h){
            d_vx[index] = (d_vx[index] - i) * spX;
            d_vy[index] = (d_vy[index] - j) * spY;
            d_vz[index] = (d_vz[index] - k) * spZ;
        }
    }
}

//////////////////////////////////////////////////////////////////////// 
// convert hfield to velocity field (inplace version)
// i.e. v(x) = (h(x) - x) * sp
////////////////////////////////////////////////////////////////////////
void
toV_I(float *vx, 
      float *vy, 
      float *vz, 
      const Vec3Di &sz,
      const Vec3Df &sp,
      StreamT stream) 
{

    dim3 threads(16,16);
    dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));
    hToVelocity_I_kernel<<<grids, threads, 0, stream>>>
	(vx, vy, vz,
	 sz.x, sz.y, sz.z,
	 sp.x, sp.y, sp.z);
}

//////////////////////////////////////////////////////////////////////// 
// convert displacement field to hfield
// i.e. h(x) = x + u(x)
////////////////////////////////////////////////////////////////////////
__global__ void uToH_I_kernel(float* d_hx, float* d_hy, float* d_hz,
                                     int w, int h, int l)
{
    uint i     = blockIdx.x * blockDim.x + threadIdx.x;
    uint j     = blockIdx.y * blockDim.y + threadIdx.y;
    uint index = j * w + i;
    if (i < w && j < h){
        for (int k=0; k<l; ++k, index+=w*h){
            d_hx[index] += i;
            d_hy[index] += j;
            d_hz[index] += k;
        }
    }
}

void
toH_I(float *vx, 
      float *vy, 
      float *vz, 
      const Vec3Di &sz,
      StreamT stream) 
{
    dim3 threads(16,16);
    dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));
    uToH_I_kernel<<<grids, threads, 0, stream>>>
	(vx, vy, vz,
	 sz.x, sz.y, sz.z);
}


////////////////////////////////////////////////////////////////////////////
// approximate the inverse of an incremental h field using according
// to the following derivation
//
// hInv(x0) = x0 + d
// x0 = h(x0 + d)
// x0 = h(x0) + d // order zero expansion
// d  = x0 - h(x0)
//
// hInv(x0) = x0 + x0 - h(x0)
//
////////////////////////////////////////////////////////////////////////////
__global__ void InverseZerothOrder_kernel(float* hInvx, float* hInvy, float* hInvz,
                                          const float* hx,const float* hy,const float* hz,
                                          int w, int h, int l)
{
    uint i     = blockIdx.x * blockDim.x + threadIdx.x;
    uint j     = blockIdx.y * blockDim.y + threadIdx.y;
    uint id    = j * w + i;
    if (i < w && j < h){
        for (int k=0; k<l; ++k, id+=w*h){
            hInvx[id] = i + i - hx[id];
            hInvy[id] = j + j - hy[id];
            hInvz[id] = k + k - hz[id];
        }
    }
}

void
InverseZerothOrder(float *a_hInvx, 
		   float *a_hInvy, 
		   float *a_hInvz, 
		   const float *a_hx, 
		   const float *a_hy, 
		   const float *a_hz, 
		   const Vec3Di &sz,
		   StreamT stream)
{
    dim3 threads(16,16);
    dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));
    InverseZerothOrder_kernel<<<grids, threads, 0, stream>>>
	(a_hInvx, a_hInvy, a_hInvz,
	 a_hx, a_hy, a_hz,
	 sz.x, sz.y, sz.z);
}

__global__ void InverseZerothOrder_I_kernel(float* hInvx, float* hInvy, float* hInvz,
                                            int w, int h, int l)
{
    uint i     = blockIdx.x * blockDim.x + threadIdx.x;
    uint j     = blockIdx.y * blockDim.y + threadIdx.y;
    uint id    = j * w + i;
    if (i < w && j < h){
        for (int k=0; k<l; ++k, id+=w*h){
            hInvx[id] = i + i - hInvx[id];
            hInvy[id] = j + j - hInvy[id];
            hInvz[id] = k + k - hInvz[id];
        }
    }
}

void
InverseZerothOrder_I(float *a_hInvx, 
		     float *a_hInvy, 
		     float *a_hInvz, 
		     const Vec3Di &sz,
		     StreamT stream)
{
    dim3 threads(16,16);
    dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));
    InverseZerothOrder_I_kernel<<<grids, threads, 0, stream>>>
	(a_hInvx, a_hInvy, a_hInvz,
	 sz.x, sz.y, sz.z);
}

//////////////////////////////////////////////////////////////////////// 
// initialize from affine transformation
// i.e. h(x) = Ax
////////////////////////////////////////////////////////////////////////
__global__ void 
initializeFromAffine_kernel(float* d_hx, float* d_hy, float* d_hz,
			    float a00, float a01, float a02,
			    float a10, float a11, float a12,
			    float a20, float a21, float a22,
			    float t0, float t1, float t2,
			    int w, int h, int l,
			    float spX, float spY, float spZ,
			    float orX, float orY, float orZ)
{
    uint i     = blockIdx.x * blockDim.x + threadIdx.x;
    uint j     = blockIdx.y * blockDim.y + threadIdx.y;
    uint index = j * w + i;
    float x, y, z;
    float xp, yp, zp;
    if (i < w && j < h){
        for (int k=0; k<l; ++k, index+=w*h){
	   x = i*spX + orX;
	   y = j*spY + orY;
	   z = k*spZ + orZ;

	   xp = x*a00 + y*a01 + z*a02 + t0;
	   yp = x*a10 + y*a11 + z*a12 + t1;
	   zp = x*a20 + y*a21 + z*a22 + t2;

	   xp = (xp-orX)/spX;
	   yp = (yp-orY)/spY;
	   zp = (zp-orZ)/spZ;

	   d_hx[index] = xp;
	   d_hy[index] = yp;
	   d_hz[index] = zp;
        }
    }
}

void
initializeFromAffine(float *d_hx, 
		     float *d_hy, 
		     float *d_hz, 
		     const Vec3Di &sz,
		     const Vec3Df &sp,
		     const Vec3Df &org,
		     const Aff3Df &aff_in,
		     bool invertAff,
		     StreamT stream)
{
   dim3 threads(16,16);
   dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));

   Aff3Df aff = aff_in;
   if(invertAff){
      if(!aff.invert()){
	 throw PyCAException(__FILE__,__LINE__,"Error, could not invert affine transform");
      }
   }

   Mat3Df &m = aff.matrix;
   Vec3Df &v = aff.vector;
   initializeFromAffine_kernel<<<grids, threads, 0, stream>>>
      (d_hx, d_hy, d_hz,
       m(0,0), m(0,1), m(0,2),
       m(1,0), m(1,1), m(1,2),
       m(2,0), m(2,1), m(2,2),
       v[0], v[1], v[2],
       sz.x, sz.y, sz.z,
       sp.x, sp.y, sp.z,
       org.x, org.y, org.z);
}

} // end namespace PyCA

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

#include "GImageOperKernels.h"
#include <pycaUtils.h>
#include "interp.h"
#include "GSplat.cuh"

// TEST make sure boost isn't included in nvcc code
#if defined(BOOST_COMPILER)
int bla[-1];
#endif

namespace PyCA {


__global__ void SubVol_kernel(float* d_o, const float* d_i,
			      int osizeX, int osizeY, int osizeZ,
			      int isizeX, int isizeY, int isizeZ,
			      int startX, int startY, int startZ)
{
   uint o_x = blockIdx.x * blockDim.x + threadIdx.x;
   uint o_y = blockIdx.y * blockDim.y + threadIdx.y;

   if (o_x < osizeX && o_y < osizeY){
      int o_id = o_x + osizeX * o_y;
      int o_wh = osizeX*osizeY;

      int i_x = o_x+startX;
      int i_y = o_y+startY;
      if ((i_x >= 0 && i_x < isizeX) && 
	  (i_y >= 0 && i_y < isizeY))
      {
	 int i_wh = isizeX*isizeY;
	 int i_last = i_wh*isizeZ;
	 int i_id = i_x + isizeX * i_y + i_wh*startZ;
	 
	 for (int o_z=0; o_z < osizeZ; 
	      ++o_z, o_id += o_wh, i_id += i_wh)
	 {
	    float v = 0.f;
	    if(i_id >= 0 && i_id < i_last) v = d_i[i_id];
	    d_o[o_id] = v;
	 }
      }else{
	 for (int o_z=0; o_z < osizeZ; 
	      ++o_z, o_id += o_wh)
         {
	    d_o[o_id] = 0.f;
	 }
      }
   }
}

void 
SubVol(float* d_o,const float* d_i, 
       const Vec3Di& oSize, const Vec3Di& iSize, 
       const Vec3Di& start, StreamT st)
{
   dim3 threads(16,16);
   dim3 grids(iDivUp(oSize.x, threads.x), iDivUp(oSize.y, threads.y));
   
   SubVol_kernel<<<grids, threads, 0, st>>>
      (d_o, d_i,
       oSize.x, oSize.y, oSize.z,
       iSize.x, iSize.y, iSize.z,
       start.x, start.y, start.z);
}

__global__ void 
SetSubVol_I_kernel(float* d_o, const float* d_i,
		   int osizeX, int osizeY, int osizeZ,
		   int isizeX, int isizeY, int isizeZ,
		   int startX, int startY, int startZ)
{
   uint i_x = blockIdx.x * blockDim.x + threadIdx.x;
   uint i_y = blockIdx.y * blockDim.y + threadIdx.y;

   if (i_x < isizeX && i_y < isizeY){
      int i_id = i_x + isizeX * i_y;
      int i_wh = isizeX*isizeY;

      int o_x = i_x+startX;
      int o_y = i_y+startY;
      if ((o_x >= 0 && o_x < osizeX) && 
	  (o_y >= 0 && o_y < osizeY))
      {
	 int o_wh = osizeX*osizeY;
	 int istartZ = PYCAMAX(0,-startZ);
	 int iendZ = PYCAMIN(isizeZ,osizeZ-startZ);
	 int o_id = o_x + osizeX * o_y + o_wh*PYCAMAX(startZ,0);
	 
	 for (int i_z=istartZ; i_z < iendZ;
	      ++i_z, i_id += i_wh, o_id += o_wh)
	 {
	    d_o[o_id] = d_i[i_id];
	 }
      }
   }
}

void 
SetSubVol_I(float* d_o,const float* d_i, 
	    const Vec3Di& oSize, const Vec3Di& iSize, 
	    const Vec3Di& start, StreamT st)
{
   dim3 threads(16,16);
   dim3 grids(iDivUp(iSize.x, threads.x), iDivUp(iSize.y, threads.y));
   
   SetSubVol_I_kernel<<<grids, threads, 0, st>>>
      (d_o, d_i,
       oSize.x, oSize.y, oSize.z,
       iSize.x, iSize.y, iSize.z,
       start.x, start.y, start.z);
}

__global__ void Shrink_kernel(float* d_o, const float* d_i,
			      float c,
			      int szX, int szY, int szZ)
{
    uint x = blockIdx.x * blockDim.x + threadIdx.x;
    uint y = blockIdx.y * blockDim.y + threadIdx.y;

    int wh = szX*szY;

    if (x < szX && y < szY){
       int id = x + y*szX;
       for (int z=0; z < szZ; ++z, id += wh){
	  float v = d_i[id];
	  d_o[id] = PYCAMAX(v-c,0.f)+PYCAMIN(v+c,0.f);
       }
    }
}

void 
Shrink(float *d_o, const float *d_i, 
       const Vec3Di &sz, float c, StreamT st)
{
   dim3 threads(16,16);
   dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));
   Shrink_kernel<<<grids, threads, 0, st>>>
      (d_o, d_i, c,
       sz.x, sz.y, sz.z);
}

__global__ void SoftAbs_kernel(float* d_o, const float* d_i,
			       float eps,
			       int szX, int szY, int szZ)
{
   uint x = blockIdx.x * blockDim.x + threadIdx.x;
   uint y = blockIdx.y * blockDim.y + threadIdx.y;
   
   int wh = szX*szY;
   
   if (x < szX && y < szY){
      int id = x + y*szX;
      for (int z=0; z < szZ; ++z, id += wh){
	 
	 float v = d_i[id];
	 if(v < -eps){
	    d_o[id] = -v-eps/2.f;
	 }else if(v > eps){
	    d_o[id] = v-eps/2.f;
	 }else{
	    d_o[id] = (v*v)/(2.f*eps);
	 }
	 
      }
   }
}

void 
SoftAbs(float *d_o, const float *d_i, const Vec3Di &sz,
	float eps, StreamT st)
{
    dim3 threads(16,16);
    dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));
    SoftAbs_kernel<<<grids, threads, 0, st>>>
	(d_o, d_i, eps,
	 sz.x, sz.y, sz.z);
}

__global__ void SoftSgn_kernel(float* d_o, const float* d_i,
			       float eps,
			       int szX, int szY, int szZ)
{
   uint x = blockIdx.x * blockDim.x + threadIdx.x;
   uint y = blockIdx.y * blockDim.y + threadIdx.y;
   
   int wh = szX*szY;
   
   if (x < szX && y < szY){
      int id = x + y*szX;
      for (int z=0; z < szZ; ++z, id += wh){
	 
	 float v = d_i[id];
	 if(v < -eps){
	    d_o[id] = -1.f;
	 }else if(v > eps){
	    d_o[id] = 1.f;
	 }else{
	    d_o[id] = v/eps;
	 }
	 
      }
   }
}

void 
SoftSgn(float *d_o, const float *d_i, const Vec3Di &sz, 
	float eps, StreamT st)
{
    dim3 threads(16,16);
    dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));
    SoftSgn_kernel<<<grids, threads, 0, st>>>
	(d_o, d_i, eps,
	 sz.x, sz.y, sz.z);
}


template<BackgroundStrategy bg, InterpT interp, bool useOriginOffset>
__global__ void Resampling_kernel(float* d_o, const float* d_i,
                                  int osizeX, int osizeY, int osizeZ,
                                  int isizeX, int isizeY, int isizeZ)
{
    uint x = blockIdx.x * blockDim.x + threadIdx.x;
    uint y = blockIdx.y * blockDim.y + threadIdx.y;

    float rX = (float)isizeX / (float)osizeX;
    float rY = (float)isizeY / (float)osizeY;
    float rZ = (float)isizeZ / (float)osizeZ;

    float offX=0.f, offY=0.f, offZ=0.f;

    if(useOriginOffset){
        offX = (rX-1.f)/2.f;
        offY = (rY-1.f)/2.f;
        offZ = (rZ-1.f)/2.f;
    }

    if (x < osizeX && y < osizeY){
        int id = x + osizeX * y;

        float i_x =  x * rX + offX;
        float i_y =  y * rY + offY;

        for (int z=0; z < osizeZ; ++z, id += osizeX * osizeY){
            float i_z = z * rZ + offZ;
            d_o[id] = point_interp<interp, bg>
		(d_i,
		 i_x, i_y, i_z,
		 isizeX, isizeY, isizeZ);
        }
    }
}

template<BackgroundStrategy bg, InterpT interp, bool useOriginOffset>
void 
Resample(float *d_o, const Vec3Di &oSz,
	 const float*d_i, const Vec3Di &iSz,
	 StreamT stream)
{
    MK_CHECK_IMAGE_BACKGROUND(bg);

    dim3 threads(16,16);
    dim3 grids(iDivUp(oSz.x, threads.x), iDivUp(oSz.y, threads.y));
    Resampling_kernel<bg, interp, useOriginOffset>
	<<<grids, threads, 0, stream>>>(d_o, d_i,
					oSz.x ,oSz.y ,oSz.z,
					iSz.x ,iSz.y ,iSz.z);
}

template<BackgroundStrategy bg, InterpT interp>
__global__ void ResampleWorld_kernel(float* d_o, const float* d_i,
				     int oSzX, int oSzY, int oSzZ,
				     float oSpX, float oSpY, float oSpZ,
				     float oOrX, float oOrY, float oOrZ,
				     int iSzX, int iSzY, int iSzZ,
				     float iSpX, float iSpY, float iSpZ,
				     float iOrX, float iOrY, float iOrZ)
{
    uint x = blockIdx.x * blockDim.x + threadIdx.x;
    uint y = blockIdx.y * blockDim.y + threadIdx.y;

    float rX = oSpX/iSpX;
    float rY = oSpY/iSpY;
    float rZ = oSpZ/iSpZ;
    float oX = (oOrX-iOrX)/iSpX;
    float oY = (oOrY-iOrY)/iSpY;
    float oZ = (oOrZ-iOrZ)/iSpZ;

    if (x < oSzX && y < oSzY){
        int id = x + oSzX * y;

        float i_x =  x*rX + oX;
        float i_y =  y*rY + oY;

        for (int z=0; z < oSzZ; ++z, id += oSzX * oSzY){
            float i_z = z*rZ + oZ;
            d_o[id] = point_interp<interp, bg>(d_i, i_x, i_y, i_z, iSzX, iSzY, iSzZ);
        }
    }
}

template<BackgroundStrategy bg, InterpT interp>
void
ResampleWorld(float *d_o, 
	      const Vec3Di &oSz, 
	      const Vec3Df &oSp,
	      const Vec3Df &oOr,
	      const float *d_i, 
	      const Vec3Di &iSz, 
	      const Vec3Df &iSp,
	      const Vec3Df &iOr,
	      StreamT stream)
{
    MK_CHECK_IMAGE_BACKGROUND(bg);
    
    dim3 threads(16,16);
    dim3 grids(iDivUp(oSz.x, threads.x), iDivUp(oSz.y, threads.y));
    ResampleWorld_kernel<bg, interp><<<grids, threads, 0, stream>>>
       (d_o, d_i,
	oSz.x, oSz.y, oSz.z,
	oSp.x, oSp.y, oSp.z,
	oOr.x, oOr.y, oOr.z,
	iSz.x, iSz.y, iSz.z,
	iSp.x, iSp.y, iSp.z,
	iOr.x, iOr.y, iOr.z);
}

template<BackgroundStrategy bg>
__global__ void SplatWorld_kernel(float* d_o, const float* d_i,
				     int oSzX, int oSzY, int oSzZ,
				     float oSpX, float oSpY, float oSpZ,
				     float oOrX, float oOrY, float oOrZ,
				     int iSzX, int iSzY, int iSzZ,
				     float iSpX, float iSpY, float iSpZ,
				     float iOrX, float iOrY, float iOrZ)
{
    uint x = blockIdx.x * blockDim.x + threadIdx.x;
    uint y = blockIdx.y * blockDim.y + threadIdx.y;

    float rX = iSpX/oSpX;
    float rY = iSpY/oSpY;
    float rZ = iSpZ/oSpZ;
    float oX = (iOrX-oOrX)/oSpX;
    float oY = (iOrY-oOrY)/oSpY;
    float oZ = (iOrZ-oOrZ)/oSpZ;

    if (x < iSzX && y < iSzY){
        int id = x + iSzX * y;

        float i_x =  x*rX + oX;
        float i_y =  y*rY + oY;

        for (int z=0; z < iSzZ; ++z, id += iSzX * iSzY){
            float i_z = z*rZ + oZ;
#if __CUDA_ARCH__ >= 200
            // floating point atomics supported only by Fermi and above
            Splatting::atomicSplatFloat(d_o, d_i[id], i_x, i_y, i_z, oSzX, oSzY, oSzZ);
#else
            Splatting::atomicSplat(reinterpret_cast<int*>(d_o), d_i[id], i_x, i_y, i_z, oSzX, oSzY, oSzZ);
#endif
        }
    }
}

template<BackgroundStrategy bg>
__global__ void SplatWorld_kernel(float* d_o, const float* d_i, float* d_w,
				     int oSzX, int oSzY, int oSzZ,
				     float oSpX, float oSpY, float oSpZ,
				     float oOrX, float oOrY, float oOrZ,
				     int iSzX, int iSzY, int iSzZ,
				     float iSpX, float iSpY, float iSpZ,
				     float iOrX, float iOrY, float iOrZ)
{
    uint x = blockIdx.x * blockDim.x + threadIdx.x;
    uint y = blockIdx.y * blockDim.y + threadIdx.y;

    float rX = iSpX/oSpX;
    float rY = iSpY/oSpY;
    float rZ = iSpZ/oSpZ;
    float oX = (iOrX-oOrX)/oSpX;
    float oY = (iOrY-oOrY)/oSpY;
    float oZ = (iOrZ-oOrZ)/oSpZ;

    if (x < iSzX && y < iSzY){
        int id = x + iSzX * y;

        float i_x =  x*rX + oX;
        float i_y =  y*rY + oY;

        for (int z=0; z < iSzZ; ++z, id += iSzX * iSzY){
            float i_z = z*rZ + oZ;
#if __CUDA_ARCH__ >= 200
            // floating point atomics supported only by Fermi and above
            Splatting::atomicSplatFloat(d_o, d_w, d_i[id], i_x, i_y, i_z, oSzX, oSzY, oSzZ);
#else
            Splatting::atomicSplat(reinterpret_cast<int*>(d_o), reinterpret_cast<int*>(d_w), d_i[id], i_x, i_y, i_z, oSzX, oSzY, oSzZ);
#endif
        }
    }
}

template<BackgroundStrategy bg>
void
SplatWorld(float *d_o, 
	   const Vec3Di &oSz,
	   const Vec3Df &oSp,
	   const Vec3Df &oOr,
	   const float *d_i, 
	   const Vec3Di &iSz,
	   const Vec3Df &iSp,
	   const Vec3Df &iOr,
	   StreamT stream)
{
    MK_CHECK_IMAGE_BACKGROUND(bg);
    
    dim3 threads(16,16);
    dim3 grids(iDivUp(oSz.x, threads.x), iDivUp(oSz.y, threads.y));
    SplatWorld_kernel<bg><<<grids, threads, 0, stream>>>
       (d_o, d_i,
	oSz.x, oSz.y, oSz.z,
	oSp.x, oSp.y, oSp.z,
	oOr.x, oOr.y, oOr.z,
	iSz.x, iSz.y, iSz.z,
	iSp.x, iSp.y, iSp.z,
	iOr.x, iOr.y, iOr.z);
}

template<BackgroundStrategy bg>
void
SplatWorld(float *d_o, 
	   const Vec3Di &oSz,
	   const Vec3Df &oSp,
	   const Vec3Df &oOr,
	   const float *d_i, 
	   const Vec3Di &iSz,
	   const Vec3Df &iSp,
	   const Vec3Df &iOr,
	   float *d_w,
	   StreamT stream)
{
    MK_CHECK_IMAGE_BACKGROUND(bg);
    
    dim3 threads(16,16);
    dim3 grids(iDivUp(oSz.x, threads.x), iDivUp(oSz.y, threads.y));
    SplatWorld_kernel<bg><<<grids, threads, 0, stream>>>
       (d_o, d_i, d_w,
	oSz.x, oSz.y, oSz.z,
	oSp.x, oSp.y, oSp.z,
	oOr.x, oOr.y, oOr.z,
	iSz.x, iSz.y, iSz.z,
	iSp.x, iSp.y, iSp.z,
	iOr.x, iOr.y, iOr.z);
}


/**
 * Note: this should be implemented with the kernel copied into shared memory
 * for small kernels, which will be the standard use case
 * jsp2012
 */
__global__ 
void 
Convolve_kernel(float* d_o, 
		const float* d_i,
		const float* d_kernel,
		int iSzX, int iSzY, int iSzZ,
		int kSzX, int kSzY, int kSzZ)
{
   uint cx = blockIdx.x * blockDim.x + threadIdx.x;
   uint cy = blockIdx.y * blockDim.y + threadIdx.y;

   if (cx < iSzX && cy < iSzY){
      
      int halfSzX = kSzX/2;
      int halfSzY = kSzY/2;
      int halfSzZ = kSzZ/2;

      for (int cz=0; cz<iSzZ; ++cz){

	 // loop over offsets in kernel
	 float v = 0.f;
	 for (int oz=-halfSzZ; oz <= halfSzZ; ++oz) {
	    for (int oy=-halfSzY; oy <= halfSzY; ++oy) {
	       for (int ox=-halfSzX; ox <= halfSzX; ++ox) {
		  float kv = getVal<float>
		     (d_kernel, 
		      kSzX,kSzY,kSzZ,
		      ox+halfSzX,oy+halfSzY,oz+halfSzZ);
		  float iv = getSafeVal<float,BACKGROUND_STRATEGY_CLAMP>
		     (d_i, 
		      iSzX,iSzY,iSzZ,
		      cx+ox,cy+oy,cz+oz);
		  v += kv*iv;
	       }
	    }
	 }
	 getVal<float>(d_o,iSzX,iSzY,iSzZ,cx,cy,cz) = v;
	 
      }
   }
}

void 
Convolve(float *d_o, const float *d_i, 
	 const Vec3Di &sz,
	 const float *d_kernel, 
	 const Vec3Di &kSz,
	 StreamT stream)
{
   dim3 threads(16,16);
   dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));
   
   if(kSz.x%2 == 0 || kSz.y%2 == 0 || kSz.z%2 == 0){
      throw PyCAException(__FILE__, __LINE__,
				   "Only odd-sized kernels allowed in convolution");
   }
   
   Convolve_kernel<<<grids, threads, 0, stream>>>
       (d_o, d_i, d_kernel,
	sz.x, sz.y, sz.z,
	kSz.x, kSz.y, kSz.z);
   
}

// template instantiations
#include "GImageOperKernels_inst.cxx"

} // end namespace PyCA

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

#include "GSplat.h"
#include "GSplat.cuh"

#include <pycaUtils.h>
#include <GMemOpers.h>

// TEST make sure boost isn't included in nvcc code
#if defined(BOOST_COMPILER)
int bla[-1];
#endif

namespace PyCA {
  namespace Splatting{


    __device__ void atomicSplatDistance(int* d_d, float x, float y, float z,
					int w, int h, int l)
    {
      int xInt = int(x);
      int yInt = int(y);
      int zInt = int(z);

      if (x < 0 && x != xInt) --xInt;
      if (y < 0 && y != yInt) --yInt;
      if (z < 0 && z != zInt) --zInt;

      float dx = 1.f - (x - xInt);
      float dy = 1.f - (y - yInt);
      float dz = 1.f - (z - zInt);
    
      uint nid = (zInt * h + yInt) * w + xInt;
      float  dist;
      if (isInside3D(xInt, yInt, zInt, w, h, l)){
        dist = dx * dy * dz;
        atomicAdd(&d_d[nid],S2p20(dist));
      }
    
      if (isInside3D(xInt + 1, yInt, zInt, w, h, l)){
        dist = (1.f-dx) * dy * dz;
        atomicAdd(&d_d[nid + 1], S2p20(dist));
      }

      if (isInside3D(xInt, yInt+1, zInt, w, h, l)){
        dist = dx * (1.f - dy) * dz;
        atomicAdd(&d_d[nid + w], S2p20(dist));
      }

      if (isInside3D(xInt+1, yInt+1, zInt, w, h, l)){
        dist = (1.f -dx) * (1.f - dy) * dz;
        atomicAdd(&d_d[nid + w + 1], S2p20(dist));
      } 
            
      nid += w*h;
      if (isInside3D(xInt, yInt, zInt + 1, w, h, l)){
        dist =  dx * dy * (1.f - dz);
        atomicAdd(&d_d[nid],S2p20(dist));
      }
            
      if (isInside3D(xInt + 1, yInt, zInt+1, w, h, l)){
        dist = (1.f-dx) * dy * (1.f -dz);
        atomicAdd(&d_d[nid + 1], S2p20(dist));
      }
    
      if (isInside3D(xInt, yInt+1, zInt+1, w, h, l)){
        dist = dx * (1.f - dy) * (1.f -dz);
        atomicAdd(&d_d[nid + w], S2p20(dist));
      }

      if (isInside3D(xInt+1, yInt+1, zInt+1, w, h, l)){
        dist = (1.f -dx) * (1.f - dy) * (1.f -dz);
        atomicAdd(&d_d[nid + w + 1], S2p20(dist));
      } 
    }

    __global__ void atomicSplatDistance_kernel(
					       int* d_d, int w, int h, int l,
					       const float* d_px, const float* d_py, const float* d_pz, uint nP)
    {
      uint blockId = get_blockID();
      uint id      = get_threadID(blockId);

      if (id >= nP) return;

      float x = d_px[id];
      float y = d_py[id];
      float z = d_pz[id];
      atomicSplatDistance(d_d, x, y, z, w, h, l);
    }


    void splatDistance(int* d_id, size_t sizeX, size_t sizeY, size_t sizeZ,
				  const float* d_px , const float* d_py, const float* d_pz,
				  size_t nP, StreamT stream)
    {
      //1.Init accumulate array 0
      size_t nVox = sizeX * sizeY * sizeZ;
      GMemOpers<int>::SetMem(d_id, 0 , nVox, stream, false);

      //2.Splat value
      dim3 threads(256);
      dim3 grids=make_grid(iDivUp(nP, threads.x));
      atomicSplatDistance_kernel<<<grids, threads, 0, stream>>>
        (d_id, sizeX, sizeY, sizeZ, d_px, d_py, d_pz, nP);
    }

    __global__ void atomicSplatPos_kernel(int* d_wd , int w, int h, int l, 
					  const float* d_w,
					  const float* d_px, const float* d_py, const float* d_pz, int nP)
    {
      uint blockId = get_blockID();
      uint id      = get_threadID(blockId);
    
      if (id < nP){
        float mass = d_w[id];
        float x = d_px[id];
        float y = d_py[id];
        float z = d_pz[id];

	Splatting::atomicSplat(d_wd, mass, x, y, z, w, h, l);
      }
    }


    void splat3D(int* d_iwd, size_t sizeX, size_t sizeY, size_t sizeZ,
			    const float* d_w,
			    const float* d_px , const float* d_py, const float* d_pz,
			    size_t nP, StreamT stream)
    {
      //1.Init accumulate array 0
      size_t nVox = sizeX * sizeY * sizeZ;
      GMemOpers<int>::SetMem(d_iwd, 0, nVox, stream, false);

      //2.Splat value
      dim3 threads(256);
      dim3 grids=make_grid(iDivUp(nP, threads.x));
      atomicSplatPos_kernel<<<grids, threads, 0, stream>>>(d_iwd, sizeX, sizeY, sizeZ,
							   d_w, d_px, d_py, d_pz, nP);
    }


    __global__ void atomicSplatPos_kernel(
					  int* d_wd , int* d_wd1, int* d_wd2, int w, int h, int l, 
					  const float* d_w, const float* d_w1, const float* d_w2,                                     const float* d_px,const float* d_py, const float* d_pz, uint nP)
    {
      uint blockId = get_blockID();
      uint id      = get_threadID(blockId);

      if (id < nP){
        float mass = d_w[id], mass1 = d_w1[id], mass2 = d_w2[id];
                
        float x = d_px[id];
        float y = d_py[id];
        float z = d_pz[id];

	Splatting::atomicSplat(d_wd, d_wd1, d_wd2,
			       mass, mass1, mass2,
			       x, y, z, w, h, l);
      }
    }

    __global__ void atomicSplatPos_kernel(
					  int* d_wd , int* d_wd1, int* d_wd2,
					  int w, int h, int l, 
					  const float* d_w, const float* d_w1, const float* d_w2,                                     const float* d_px,const float* d_py, const float* d_pz)
    {
      uint i = blockIdx.x * blockDim.x + threadIdx.x;
      uint j = blockIdx.y * blockDim.y + threadIdx.y;
    
      if (i < w && j < h){
        int id = i + w * j;
        for (int k=0; k < l; ++k, id+=w*h){
	  float mass = d_w[id], mass1 = d_w1[id], mass2 = d_w2[id];
	  float x = d_px[id];
	  float y = d_py[id];
	  float z = d_pz[id];
	  Splatting::atomicSplat(d_wd, d_wd1, d_wd2,
				 mass, mass1, mass2,
				 x, y, z, w, h, l);
        }
      }
    }


    void splat3D(int* d_iwdx, int* d_iwdy, int* d_iwdz, 
			    size_t sizeX, size_t sizeY, size_t sizeZ,
			    const float* d_wx, const float* d_wy, const float* d_wz,
			    const float* d_px , const float* d_py, const float* d_pz,
			    size_t nP, StreamT stream)
    {
      //1.Init accumulate array 0
      size_t nVox = sizeX * sizeY * sizeZ;
      GMemOpers<int>::SetMem(d_iwdx, 0, nVox, stream, false);
      GMemOpers<int>::SetMem(d_iwdy, 0, nVox, stream, false);
      GMemOpers<int>::SetMem(d_iwdz, 0, nVox, stream, false);
    
      //2.Splat value
      dim3 threads(256);
      dim3 grids=make_grid(iDivUp(nP, threads.x));
      atomicSplatPos_kernel<<<grids, threads, 0, stream>>>(d_iwdx, d_iwdy, d_iwdz,
							   sizeX, sizeY, sizeZ,
							   d_wx, d_wy, d_wz,
							   d_px, d_py, d_pz, nVox);
    }


    __device__ void atomicSplatWeightPos(int* d_wd, int* d_d,
					 float mass, float x, float y, float z,
					 int w, int h, int l)
    {
      int xInt = int(x);
      int yInt = int(y);
      int zInt = int(z);

      if (x < 0 && x != xInt) --xInt;
      if (y < 0 && y != yInt) --yInt;
      if (z < 0 && z != zInt) --zInt;

      float dx = 1.f - (x - xInt);
      float dy = 1.f - (y - yInt);
      float dz = 1.f - (z - zInt);
    
      uint nid = (zInt * h + yInt) * w + xInt;
      float  dist;
      if (isInside3D(xInt, yInt, zInt, w, h, l)){
        dist = dx * dy * dz;
        atomicAdd(&d_wd[nid],S2p20(mass * dist));
        atomicAdd(&d_d[nid],S2p20(dist));
      }
    
      if (isInside3D(xInt + 1, yInt, zInt, w, h, l)){
        dist = (1.f-dx) * dy * dz;
        atomicAdd(&d_wd[nid + 1], S2p20(mass * dist));
        atomicAdd(&d_d[nid + 1], S2p20(dist));
      }

      if (isInside3D(xInt, yInt+1, zInt, w, h, l)){
        dist = dx * (1.f - dy) * dz;
        atomicAdd(&d_wd[nid + w], S2p20(mass * dist));
        atomicAdd(&d_d[nid + w], S2p20(dist));
      }

      if (isInside3D(xInt+1, yInt+1, zInt, w, h, l)){
        dist = (1.f -dx) * (1.f - dy) * dz;
        atomicAdd(&d_wd[nid + w + 1], S2p20(mass * dist));
        atomicAdd(&d_d[nid + w + 1], S2p20(dist));
      } 
            
      nid += w*h;
      if (isInside3D(xInt, yInt, zInt + 1, w, h, l)){
        dist =  dx * dy * (1.f - dz);
        atomicAdd(&d_wd[nid],S2p20(mass * dist));
        atomicAdd(&d_d[nid],S2p20(dist));
      }
            
      if (isInside3D(xInt + 1, yInt, zInt+1, w, h, l)){
        dist = (1.f-dx) * dy * (1.f -dz);
        atomicAdd(&d_wd[nid + 1], S2p20(mass * dist));
        atomicAdd(&d_d[nid + 1], S2p20(dist));
      }
    
      if (isInside3D(xInt, yInt+1, zInt+1, w, h, l)){
        dist = dx * (1.f - dy) * (1.f -dz);
        atomicAdd(&d_wd[nid + w], S2p20(mass * dist));
        atomicAdd(&d_d[nid + w], S2p20(dist));
      }

      if (isInside3D(xInt+1, yInt+1, zInt+1, w, h, l)){
        dist = (1.f -dx) * (1.f - dy) * (1.f -dz);
        atomicAdd(&d_wd[nid + w + 1], S2p20(mass * dist));
        atomicAdd(&d_d[nid + w + 1], S2p20(dist));
      } 
    }

    __global__ void atomicSplatWeightPos_kernel(
						int* d_wd, int* d_d, int w, int h, int l,
						const float* d_w, const float* d_px, const float* d_py, const float* d_pz, uint nP)
    {
      uint blockId = get_blockID();
      uint id      = get_threadID(blockId);

      if (id >= nP) return;

      float mass = d_w[id];
      float x = d_px[id];
      float y = d_py[id];
      float z = d_pz[id];

      atomicSplatWeightPos(d_wd, d_d, mass, x, y, z, w, h, l);
    }

    void splat3D(int* d_iwd, int* d_id,
			    uint sizeX, uint sizeY, uint sizeZ,
			    const float* d_w,
			    const float* d_px , const float* d_py, const float* d_pz,
			    uint nP, StreamT stream)
    {
      //1.Init accumulate array 0
      size_t nVox = sizeX * sizeY * sizeZ;
      GMemOpers<int>::SetMem(d_iwd, 0, nVox, stream, false);
      GMemOpers<int>::SetMem(d_id, 0 , nVox, stream, false);

      //2.Splat value
      dim3 threads(256);
      dim3 grids=make_grid(iDivUp(nP, threads.x));
      atomicSplatWeightPos_kernel<<<grids, threads, 0, stream>>>
        (d_iwd, d_id, sizeX, sizeY, sizeZ,
         d_w, d_px, d_py, d_pz, nP);
    }


    __global__ void atomicSplatV_kernel(
					int* d_wd , int* d_wd1, int* d_wd2,
					int w, int h, int l, 
					const float* d_w, const float* d_w1, const float* d_w2,                                     const float* d_vx,const float* d_vy, const float* d_vz,
					float iSpx, float iSpy, float iSpz)
    {
      uint i = blockIdx.x * blockDim.x + threadIdx.x;
      uint j = blockIdx.y * blockDim.y + threadIdx.y;
    
      if (i < w && j < h){
        int id = i + w * j;
        for (int k=0; k < l; ++k, id+=w*h){
	  float mass = d_w[id], mass1 = d_w1[id], mass2 = d_w2[id];
	  float x = i + d_vx[id] * iSpx;
	  float y = j + d_vy[id] * iSpy;
	  float z = k + d_vz[id] * iSpz;
	  Splatting::atomicSplat(d_wd, d_wd1, d_wd2,
				 mass, mass1, mass2,
				 x, y, z, w, h, l);
        }
      }
    }

    __device__ void atomicSplatWeightPos(int* d_wd, int* d_wd1, int* d_wd2,
					 int* d_d, float mass, float mass1, float mass2,
					 float x, float y, float z,
					 int w, int h, int l)
    {
      int xInt = int(x);
      int yInt = int(y);
      int zInt = int(z);

      if (x < 0 && x != xInt) --xInt;
      if (y < 0 && y != yInt) --yInt;
      if (z < 0 && z != zInt) --zInt;

      float dx = 1.f - (x - xInt);
      float dy = 1.f - (y - yInt);
      float dz = 1.f - (z - zInt);

      uint nid = (zInt * h + yInt) * w + xInt;
      int dist;
      float weight;
    
      if (isInside3D(xInt, yInt, zInt, w, h, l)){
        weight = dx * dy * dz;
        
        atomicAdd(&d_d[nid],S2p20(weight));
        
        dist = S2p20(mass * weight);
        atomicAdd(&d_wd[nid],dist);

        dist = S2p20(mass1 * weight);
        atomicAdd(&d_wd1[nid],dist);

        dist = S2p20(mass2 * weight);
        atomicAdd(&d_wd2[nid],dist);
      }
            
      if (isInside3D(xInt + 1, yInt, zInt, w, h, l)){
        weight = (1.f-dx) * dy * dz;

        atomicAdd(&d_d[nid + 1],S2p20(weight));
        
        dist = S2p20(mass * weight);
        atomicAdd(&d_wd[nid + 1], dist);

        dist = S2p20(mass1 * weight);
        atomicAdd(&d_wd1[nid + 1], dist);
        
        dist = S2p20(mass2 * weight);
        atomicAdd(&d_wd2[nid + 1], dist);
      }

      if (isInside3D(xInt, yInt+1, zInt, w, h, l)){
        weight = dx * (1.f - dy) * dz;

        atomicAdd(&d_d[nid + w],S2p20(weight));
        
        dist = S2p20(mass * weight);
        atomicAdd(&d_wd[nid + w], dist);

        dist = S2p20(mass1 * weight);
        atomicAdd(&d_wd1[nid + w], dist);

        dist = S2p20(mass2 * weight);
        atomicAdd(&d_wd2[nid + w], dist);
      }
    
      if (isInside3D(xInt+1, yInt+1, zInt, w, h, l)){
        weight = (1.f -dx) * (1.f - dy) * dz;

        atomicAdd(&d_d[nid + 1 + w],S2p20(weight));
        
        dist = S2p20(mass * weight);
        atomicAdd(&d_wd[nid + w + 1], dist);

        dist = S2p20(mass1 * weight);
        atomicAdd(&d_wd1[nid + w + 1], dist);

        dist = S2p20(mass2 * weight);
        atomicAdd(&d_wd2[nid + w + 1], dist);
      } 
    
      nid += w*h;

      if (isInside3D(xInt, yInt, zInt + 1, w, h, l)){
        weight = dx * dy * (1.f - dz);

        atomicAdd(&d_d[nid],S2p20(weight));
        
        dist = S2p20(mass * weight);
        atomicAdd(&d_wd[nid],dist);

        dist = S2p20(mass1 * weight);
        atomicAdd(&d_wd1[nid],dist);

        dist = S2p20(mass2 * weight);
        atomicAdd(&d_wd2[nid],dist);
      }
            
      if (isInside3D(xInt + 1, yInt, zInt+1, w, h, l)){
        weight =  (1.f-dx) * dy * (1.f -dz);

        atomicAdd(&d_d[nid + 1],S2p20(weight));
        
        dist = S2p20(mass * weight);
        atomicAdd(&d_wd[nid + 1], dist);

        dist = S2p20(mass1 * weight);
        atomicAdd(&d_wd1[nid + 1], dist);
        
        dist = S2p20(mass2 * weight);
        atomicAdd(&d_wd2[nid + 1], dist);
      }
    
      if (isInside3D(xInt, yInt+1, zInt+1, w, h, l)){
        weight =  dx * (1.f - dy) * (1.f -dz);

        atomicAdd(&d_d[nid + w],S2p20(weight));
        
        dist = S2p20(mass * weight);
        atomicAdd(&d_wd[nid + w], dist);

        dist = S2p20(mass1 * weight);
        atomicAdd(&d_wd1[nid + w], dist);

        dist = S2p20(mass2 * weight);
        atomicAdd(&d_wd2[nid + w], dist);
      }
    
      if (isInside3D(xInt+1, yInt+1, zInt+1, w, h, l)){
        weight = (1.f -dx) * (1.f - dy) * (1.f -dz);

        atomicAdd(&d_d[nid + 1 + w],S2p20(weight));
        
        dist = S2p20(mass * weight);
        atomicAdd(&d_wd[nid + w + 1], dist);

        dist = S2p20(mass1 * weight);
        atomicAdd(&d_wd1[nid + w + 1], dist);

        dist = S2p20(mass2 * weight);
        atomicAdd(&d_wd2[nid + w + 1], dist);
      } 
    }

    __global__ void atomicSplatWeightPos_kernel(
						int* d_wd , int* d_wd1, int* d_wd2, int * d_d,
						int w, int h, int l, 
						const float* d_w, const float* d_w1, const float* d_w2,
						const float* d_px,const float* d_py, const float* d_pz, uint nP)
    {
      uint blockId = get_blockID();
      uint id      = get_threadID(blockId);

      if (id >= nP)
        return;
    
      float mass = d_w[id], mass1 = d_w1[id], mass2 = d_w2[id];
    
      float x = d_px[id];
      float y = d_py[id];
      float z = d_pz[id];

      atomicSplatWeightPos(d_wd, d_wd1, d_wd2, d_d,
			   mass, mass1, mass2,
			   x, y, z, w, h, l);
    }


    void splat3D(int* d_iwdx, int* d_iwdy, int* d_iwdz, int* d_id,
			    size_t sizeX, size_t sizeY, size_t sizeZ,
			    const float* d_wx, const float* d_wy, const float* d_wz,
			    const float* d_px , const float* d_py, const float* d_pz,
			    size_t nP, StreamT stream)
    {
      //1.Init accumulate array 0
      size_t nVox = sizeX * sizeY * sizeZ;
      GMemOpers<int>::SetMem(d_iwdx, 0, nVox, stream, false);
      GMemOpers<int>::SetMem(d_iwdy, 0, nVox, stream, false);
      GMemOpers<int>::SetMem(d_iwdz, 0, nVox, stream, false);
      GMemOpers<int>::SetMem(d_id, 0, nVox, stream, false);
    

      //2.Splat value
      dim3 threads(256);
      dim3 grids=make_grid(iDivUp(nP, threads.x));
      atomicSplatWeightPos_kernel<<<grids, threads, 0, stream>>>(d_iwdx, d_iwdy, d_iwdz, d_id,
								 sizeX, sizeY, sizeZ,
								 d_wx, d_wy, d_wz,
								 d_px, d_py, d_pz, nP);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // The safe version do the normalization on the data first
    // so that the input of the data in the range of [0,1]
    // and scale the data back to the original range
    ////////////////////////////////////////////////////////////////////////////////
    template<bool inverse>
    __global__ void atomicSplatPos_kernel(int* d_wd , int w, int h, int l, 
					  const float* d_w, float max,
					  const float* d_px, const float* d_py, const float* d_pz, int nP)
    {
      uint blockId = get_blockID();
      uint id      = get_threadID(blockId);

      if (id < nP){
        float mass;
        if (inverse){
	  mass = d_w[id] * max; // normalized the mass
        } else {
	  mass = d_w[id] / max; // normalized the mass
        }
        float x = d_px[id];
        float y = d_py[id];
        float z = d_pz[id];

	Splatting::atomicSplat(d_wd, mass, x, y, z, w, h, l);
      }
    }




    ////////////////////////////////////////////////////////////////////////////////
    // Fixed point division
    // 
    ////////////////////////////////////////////////////////////////////////////////
    // d_fwd : floating point weighted distance
    // d_fd  : floating point total linear interpolation distance 
    // d_iwd : fixed point weighted distance
    // d_ifd : fixed point total linear interpolation distance 
    ////////////////////////////////////////////////////////////////////////////////

    __global__ void convertWeightedDistance_kernel(float* d_fwd,
						   const int* d_iwd, const int* d_id, uint n)
    {
      uint blockId = get_blockID();
      uint id      = get_threadID(blockId);

      if (id < n) {
        d_fwd[id] = (d_id[id] == 0) ? 0.f : float(d_iwd[id]) / float(d_id[id]);
      }
    }

    void convertWeightedDistance(float* d_fwd, const int* d_iwd,
					    const int* d_id, size_t n, StreamT stream)
    {
      dim3 threads(256);
      dim3 grids = make_grid(iDivUp(n, threads.x));
      convertWeightedDistance_kernel<<<grids, threads, 0, stream>>>(d_fwd, d_iwd, d_id, n);
    }

    __global__ void convertWeightedDistance_I_kernel(float* d_fwd, const int* d_id, uint n)
    {
      uint blockId = get_blockID();
      uint id      = get_threadID(blockId);

      if (id < n) {
        d_fwd[id] = (d_id[id] == 0) ? 0.f : float(__float_as_int(d_fwd[id])) / float(d_id[id]);
      }
    }

    void convertWeightedDistance_I(float* d_fwd, const int* d_id, size_t n, StreamT stream)
    {
      dim3 threads(256);
      dim3 grids = make_grid(iDivUp(n, threads.x));
      convertWeightedDistance_I_kernel<<<grids, threads, 0, stream>>>(d_fwd, d_id, n);
    }

    __global__ void convertWeightedDistance_I_kernel(float* d_fx, float* d_fy, float* d_fz,
						     const int* d_id, uint n)
    {
      uint blockId = get_blockID();
      uint id      = get_threadID(blockId);

      if (id < n) {
        if (d_id[id] == 0) {
	  d_fx[id] = d_fy[id] = d_fz[id] = 0.f;
        } else {
	  d_fx[id] = float(__float_as_int(d_fx[id])) / float(d_id[id]);
	  d_fy[id] = float(__float_as_int(d_fy[id])) / float(d_id[id]);
	  d_fz[id] = float(__float_as_int(d_fz[id])) / float(d_id[id]);
        }
      }
    }

    void convertWeightedDistance_I(float* d_fx, float* d_fy, float* d_fz,
					      const int* d_id, size_t n, StreamT stream)
    {
      dim3 threads(256);
      dim3 grids = make_grid(iDivUp(n, threads.x));
      convertWeightedDistance_I_kernel<<<grids, threads, 0, stream>>>(d_fx, d_fy, d_fz,
								      d_id, n);
    }


    // __global__ void convertWeightedDistance_kernel(float* d_fwd, float* d_fd,
    //                                                int* d_iwd, int* d_id, uint n)
    // {
    //     uint blockId = get_blockID();
    //     uint id      = get_threadID(blockId);

    //     if (id < n){
    //         if (d_id[id] == 0) {
    //             d_fwd[id] = 0.f;
    //             d_fd[id] = 0.f;
    //         }
    //         else {
    //             d_fwd[id] = float(d_iwd[id]) / float(d_id[id]);
    //             d_fd[id]  = S2n20(d_id[id]);
    //         }
    //     }
    // }

    // void convertWeightedDistance_fixed(float* d_fwd, float* d_fd,
    //                                    int* d_iwd, int* d_id, uint n, StreamT stream)
    // {
    //     dim3 threads(256);
    //     dim3 grids(iDivUp(n, threads.x));
    //     checkConfig(grids);
    //     convertWeightedDistance_kernel<<<grids, threads, 0, stream>>>
    //         (d_fwd, d_fd, d_iwd, d_id, n);
    // }



    __global__ void atomicVelocitySplat_kernel_shared(int* d_wd, const float* d_w,
						      const float* vx, const float* vy, const float* vz,
						      int w, int h, int l)
    {
      __shared__ int s_0[16*16];
      __shared__ int s_1[16*16];
      __shared__ int s_2[16*16];

      const uint wh     = w * h;

      int xc = blockIdx.x * blockDim.x;
      int yc = blockIdx.y * blockDim.y;
    
      int i  = xc + threadIdx.x;
      int j  = yc + threadIdx.y;
    
      if (i < w && j < h){
        uint id       = i + j * w;
        s_0[threadIdx.y * blockDim.x + threadIdx.x] = 0;
        s_1[threadIdx.y * blockDim.x + threadIdx.x] = 0;
        
        int* s_p = s_0, *s_c = s_1, *s_n = s_2;

        for (int k=0; k < l; ++k, id+=wh) {
	  // Initialize the new buffer with zero 
	  s_n[threadIdx.y * blockDim.x + threadIdx.x] = 0;

	  //__syncthreads();

	  float mass = d_w[id];
            
	  float x = i + vx[id];
	  float y = j + vy[id];
	  float z = k + vz[id];

	  int xInt = int(x);
	  int yInt = int(y);
	  int zInt = int(z);

	  if (x < 0 && x != xInt) --xInt;
	  if (y < 0 && y != yInt) --yInt;
	  if (z < 0 && z != zInt) --zInt;
            
	  float dx = 1.f - (x - xInt);
	  float dy = 1.f - (y - yInt);
	  float dz = 1.f - (z - zInt);

	  uint new_id = (zInt * h + yInt) * w + xInt;
	  int dist;
            
	  if (isInside3D(xInt - xc, yInt - yc, zInt + 1 - k,
			 blockDim.x-1,  blockDim.y-1, 2))
            {
	      int* s_l0, *s_l1;
                
	      if (zInt == k){
		s_l0 = s_c;
		s_l1 = s_n;
	      }
	      else {
		s_l0 = s_p;
		s_l1 = s_c;
	      }

	      uint sid = (xInt - xc) + (yInt-yc) * 16;
	      dist = S2p20(mass * dx * dy * dz);
	      atomicAdd(s_l0 + sid, dist);

	      dist = S2p20(mass * (1.f-dx) * dy * dz);
	      atomicAdd(s_l0 + sid + 1, dist);
                
	      dist = S2p20(mass * dx * (1.f - dy) * dz);
	      atomicAdd(s_l0 + sid + 16, dist);

	      dist = S2p20(mass * (1.f -dx) * (1.f - dy) * dz);
	      atomicAdd(s_l0 + sid + 16 +1, dist);
                
	      dist = S2p20(mass * dx * dy * (1-dz));
	      atomicAdd(s_l1 + sid, dist);
                
	      dist = S2p20(mass * (1.f-dx) * dy * (1-dz));
	      atomicAdd(s_l1 + sid + 1, dist);
                
	      dist = S2p20(mass * dx * (1.f - dy) * (1-dz));
	      atomicAdd(s_l1 + sid + 16, dist);
                    
	      dist = S2p20(mass * (1.f -dx) * (1.f - dy) * (1-dz));
	      atomicAdd(s_l1 + sid + 16 +1, dist);
            }else 
#if 1
	    if (isInside3D(xInt, yInt, zInt, w-1, h-1, l-1)){
	      dist = S2p20(mass * dx * dy * dz);
	      atomicAdd(&d_wd[new_id],dist);

	      dist = S2p20(mass * (1.f-dx) * dy * dz);
	      atomicAdd(&d_wd[new_id + 1], dist);

	      dist = S2p20(mass * dx * (1.f - dy) * dz);
	      atomicAdd(&d_wd[new_id + w], dist);

	      dist = S2p20(mass * (1.f -dx) * (1.f - dy) * dz);
	      atomicAdd(&d_wd[new_id + w + 1], dist);

	      new_id += w*h;

	      dist = S2p20(mass * dx * dy * (1.f - dz));
	      atomicAdd(&d_wd[new_id],dist);

	      dist = S2p20(mass * (1.f-dx) * dy * (1.f -dz));
	      atomicAdd(&d_wd[new_id + 1], dist);

	      dist = S2p20(mass * dx * (1.f - dy) * (1.f -dz));
	      atomicAdd(&d_wd[new_id + w], dist);

	      dist = S2p20(mass * (1.f -dx) * (1.f - dy) * (1.f -dz));
	      atomicAdd(&d_wd[new_id + w + 1], dist);
	    }
#else
	  {
	    if (isInside3D(xInt, yInt, zInt, w, h, l)){
	      dist = S2p20(mass * dx * dy * dz);
	      atomicAdd(&d_wd[new_id],dist);
	    }
            
	    if (isInside3D(xInt + 1, yInt, zInt, w, h, l)){
	      dist = S2p20(mass * (1.f-dx) * dy * dz);
	      atomicAdd(&d_wd[new_id + 1], dist);
	    }
	    if (isInside3D(xInt, yInt+1, zInt, w, h, l)){
	      dist = S2p20(mass * dx * (1.f - dy) * dz);
	      atomicAdd(&d_wd[new_id + w], dist);
	    }
	    if (isInside3D(xInt+1, yInt+1, zInt, w, h, l)){
	      dist = S2p20(mass * (1.f -dx) * (1.f - dy) * dz);
	      atomicAdd(&d_wd[new_id + w + 1], dist);
	    } 
	    new_id += w*h;
	    if (isInside3D(xInt, yInt, zInt + 1, w, h, l)){
	      dist = S2p20(mass * dx * dy * (1.f - dz));
	      atomicAdd(&d_wd[new_id],dist);
	    }
	    if (isInside3D(xInt + 1, yInt, zInt+1, w, h, l)){
	      dist = S2p20(mass * (1.f-dx) * dy * (1.f -dz));
	      atomicAdd(&d_wd[new_id + 1], dist);
	    }

	    if (isInside3D(xInt, yInt+1, zInt+1, w, h, l)){
	      dist = S2p20(mass * dx * (1.f - dy) * (1.f -dz));
	      atomicAdd(&d_wd[new_id + w], dist);
	    }
	    if (isInside3D(xInt+1, yInt+1, zInt+1, w, h, l)){
	      dist = S2p20(mass * (1.f -dx) * (1.f - dy) * (1.f -dz));
	      atomicAdd(&d_wd[new_id + w + 1], dist);
	    }
	  }
#endif
	  __syncthreads();
            
	  //write out the previous layer 
	  if( k > 0){
	    atomicAdd(&d_wd[id - wh], s_p[threadIdx.x + threadIdx.y * 16]);
	  }

	  //write out the current layer if it is the last 
	  if ( k == l - 1){
	    atomicAdd(&d_wd[id], s_c[threadIdx.x + threadIdx.y * 16]);
	  }
            
	  int* temp = s_p;
	  s_p = s_c;
	  s_c = s_n;
	  s_n = temp;
        }
      }
    }


    __global__ void FixedToFloating_kernel(float* d_o, const int* d_i, uint n){
      uint blockId = get_blockID();
      uint id      = get_threadID(blockId);
      if (id < n)
        d_o[id] = S2n20(d_i[id]);
    }

    void FixedToFloating(float* d_o, const int* d_i, size_t n, StreamT stream)
    {
      dim3 threads(REG_BLOCK_SIZE);
      dim3 grids=make_grid(iDivUp(n, threads.x));
      FixedToFloating_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, n);
    }

    __global__ void FixedToFloating_kernel(float* d_ox, float* d_oy, float* d_oz,
					   const int* d_ix, const int* d_iy,const int* d_iz,
					   uint n){
      uint blockId = get_blockID();
      uint id      = get_threadID(blockId);
      if (id < n) {
        d_ox[id] = S2n20(d_ix[id]);
        d_oy[id] = S2n20(d_iy[id]);
        d_oz[id] = S2n20(d_iz[id]);
      }
    }

    void FixedToFloating(float* d_ox, float* d_oy, float* d_oz,
				    const int* d_ix, const int* d_iy,const int* d_iz,
				    size_t n, StreamT stream){
      dim3 threads(REG_BLOCK_SIZE);
      dim3 grids=make_grid(iDivUp(n, threads.x));
      FixedToFloating_kernel<<<grids, threads, 0, stream>>>(d_ox, d_oy, d_oz,
							    d_ix, d_iy, d_iz,  n);
    }

    __global__ void FixedToFloating_I_kernel(float* d_o, uint n){
      uint blockId = get_blockID();
      uint id      = get_threadID(blockId);
      if (id < n) {
        int v   = __float_as_int(d_o[id]);
        d_o[id] = S2n20(v);
      }
    }

    void FixedToFloating_I(float* d_o, size_t n, StreamT stream){
      dim3 threads(REG_BLOCK_SIZE);
      dim3 grids=make_grid(iDivUp(n, threads.x));
      FixedToFloating_I_kernel<<<grids, threads, 0, stream>>>(d_o, n);
    }

    __global__ void FixedToFloating_I_kernel(float* d_ox, float* d_oy, float* d_oz, uint n){
      uint blockId = get_blockID();
      uint id      = get_threadID(blockId);
      if (id < n) {
        int vx   = __float_as_int(d_ox[id]);
        int vy   = __float_as_int(d_oy[id]);
        int vz   = __float_as_int(d_oz[id]);
        
        d_ox[id] = S2n20(vx);
        d_oy[id] = S2n20(vy);
        d_oz[id] = S2n20(vz);
      }
    }

    void FixedToFloating_I(float* d_ox, float* d_oy, float* d_oz, size_t n, StreamT stream){
      dim3 threads(REG_BLOCK_SIZE);
      dim3 grids=make_grid(iDivUp(n, threads.x));
      FixedToFloating_I_kernel<<<grids, threads, 0, stream>>>(d_ox, d_oy, d_oz, n);
    }
  } // end namespace Splatting
} // end namespace PyCA

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

#ifndef __GSPLAT_CUH
#define __GSPLAT_CUH

namespace PyCA {
  namespace Splatting{

    inline  __device__ void atomicSplat(int* d_wd, float mass, float x, float y, float z,
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

      int nid = (zInt * h + yInt) * w + xInt;

      int dist;

      if (isInside3D(xInt, yInt, zInt, w, h, l)){
	dist = S2p20(mass * dx * dy * dz);
	atomicAdd(&d_wd[nid],dist);
      }
            
      if (isInside3D(xInt + 1, yInt, zInt, w, h, l)){
	dist = S2p20(mass * (1.f-dx) * dy * dz);
	atomicAdd(&d_wd[nid + 1], dist);
      }

      if (isInside3D(xInt, yInt+1, zInt, w, h, l)){
	dist = S2p20(mass * dx * (1.f - dy) * dz);
	atomicAdd(&d_wd[nid + w], dist);
      }
    
      if (isInside3D(xInt+1, yInt+1, zInt, w, h, l)){
	dist = S2p20(mass * (1.f -dx) * (1.f - dy) * dz);
	atomicAdd(&d_wd[nid + w + 1], dist);
      } 
    
      nid += w*h;

      if (isInside3D(xInt, yInt, zInt + 1, w, h, l)){
	dist = S2p20(mass * dx * dy * (1.f - dz));
	atomicAdd(&d_wd[nid],dist);
      }
            
      if (isInside3D(xInt + 1, yInt, zInt+1, w, h, l)){
	dist = S2p20(mass * (1.f-dx) * dy * (1.f -dz));
	atomicAdd(&d_wd[nid + 1], dist);
      }
    
      if (isInside3D(xInt, yInt+1, zInt+1, w, h, l)){
	dist = S2p20(mass * dx * (1.f - dy) * (1.f -dz));
	atomicAdd(&d_wd[nid + w], dist);
      }
    
      if (isInside3D(xInt+1, yInt+1, zInt+1, w, h, l)){
	dist = S2p20(mass * (1.f -dx) * (1.f - dy) * (1.f -dz));
	atomicAdd(&d_wd[nid + w + 1], dist);
      } 
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Splat mass along with weight vector
    ///  d_wd = sum dist * d_w (weighted distance)
    ////////////////////////////////////////////////////////////////////////////////
    inline  __device__ void atomicSplat(int* d_wd, int* d_ww, float mass, float x, float y, float z,
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

      int nid = (zInt * h + yInt) * w + xInt;

      int dist, ww;

      if (isInside3D(xInt, yInt, zInt, w, h, l)){
	dist = S2p20(mass * dx * dy * dz);
	ww = S2p20(dx * dy * dz);
	atomicAdd(&d_wd[nid],dist);
	atomicAdd(&d_ww[nid],ww);
      }
            
      if (isInside3D(xInt + 1, yInt, zInt, w, h, l)){
	dist = S2p20(mass * (1.f-dx) * dy * dz);
	ww = S2p20((1.f-dx) * dy * dz);
	atomicAdd(&d_wd[nid + 1], dist);
	atomicAdd(&d_ww[nid + 1], ww);
      }

      if (isInside3D(xInt, yInt+1, zInt, w, h, l)){
	dist = S2p20(mass * dx * (1.f - dy) * dz);
	ww = S2p20(dx * (1.f - dy) * dz);
	atomicAdd(&d_wd[nid + w], dist);
	atomicAdd(&d_ww[nid + w], ww);
      }
    
      if (isInside3D(xInt+1, yInt+1, zInt, w, h, l)){
	dist = S2p20(mass * (1.f -dx) * (1.f - dy) * dz);
	ww = S2p20((1.f -dx) * (1.f - dy) * dz);
	atomicAdd(&d_wd[nid + w + 1], dist);
	atomicAdd(&d_ww[nid + w + 1], ww);
      } 
    
      nid += w*h;

      if (isInside3D(xInt, yInt, zInt + 1, w, h, l)){
	dist = S2p20(mass * dx * dy * (1.f - dz));
	ww = S2p20(dx * dy * (1.f - dz));
	atomicAdd(&d_wd[nid], dist);
	atomicAdd(&d_ww[nid], ww);
      }
            
      if (isInside3D(xInt + 1, yInt, zInt+1, w, h, l)){
	dist = S2p20(mass * (1.f-dx) * dy * (1.f -dz));
	ww = S2p20((1.f-dx) * dy * (1.f -dz));
	atomicAdd(&d_wd[nid + 1], dist);
	atomicAdd(&d_ww[nid + 1], ww);
      }
    
      if (isInside3D(xInt, yInt+1, zInt+1, w, h, l)){
	dist = S2p20(mass * dx * (1.f - dy) * (1.f -dz));
	ww = S2p20(dx * (1.f - dy) * (1.f -dz));
	atomicAdd(&d_wd[nid + w], dist);
	atomicAdd(&d_ww[nid + w], ww);
      }
    
      if (isInside3D(xInt+1, yInt+1, zInt+1, w, h, l)){
	dist = S2p20(mass * (1.f -dx) * (1.f - dy) * (1.f -dz));
	ww = S2p20((1.f -dx) * (1.f - dy) * (1.f -dz));
	atomicAdd(&d_wd[nid + w + 1], dist);
	atomicAdd(&d_ww[nid + w + 1], ww);
      } 
    }
    ////////////////////////////////////////////////////////////////////////////////
    /// Vector version: Splat both the function p_u to the neighbor point to the grid
    ///  d_wd = sum dist * d_w (weighted distance)
    ////////////////////////////////////////////////////////////////////////////////
    inline  __device__ void atomicSplat(int* d_wd, int* d_wd1, int* d_wd2,
					float mass, float mass1, float mass2,
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

      int nid = (zInt * h + yInt) * w + xInt;
      int dist;
      float weight;
    
      if (isInside3D(xInt, yInt, zInt, w, h, l)){
	weight = dx * dy * dz;
        
	dist = S2p20(mass * weight);
	atomicAdd(&d_wd[nid],dist);

	dist = S2p20(mass1 * weight);
	atomicAdd(&d_wd1[nid],dist);

	dist = S2p20(mass2 * weight);
	atomicAdd(&d_wd2[nid],dist);
      }
            
      if (isInside3D(xInt + 1, yInt, zInt, w, h, l)){
	weight = (1.f-dx) * dy * dz;
        
	dist = S2p20(mass * weight);
	atomicAdd(&d_wd[nid + 1], dist);

	dist = S2p20(mass1 * weight);
	atomicAdd(&d_wd1[nid + 1], dist);
        
	dist = S2p20(mass2 * weight);
	atomicAdd(&d_wd2[nid + 1], dist);
      }

      if (isInside3D(xInt, yInt+1, zInt, w, h, l)){
	weight = dx * (1.f - dy) * dz;
        
	dist = S2p20(mass * weight);
	atomicAdd(&d_wd[nid + w], dist);

	dist = S2p20(mass1 * weight);
	atomicAdd(&d_wd1[nid + w], dist);

	dist = S2p20(mass2 * weight);
	atomicAdd(&d_wd2[nid + w], dist);
      }
    
      if (isInside3D(xInt+1, yInt+1, zInt, w, h, l)){
	weight = (1.f -dx) * (1.f - dy) * dz;
        
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
        
	dist = S2p20(mass * weight);
	atomicAdd(&d_wd[nid],dist);

	dist = S2p20(mass1 * weight);
	atomicAdd(&d_wd1[nid],dist);

	dist = S2p20(mass2 * weight);
	atomicAdd(&d_wd2[nid],dist);
      }
            
      if (isInside3D(xInt + 1, yInt, zInt+1, w, h, l)){
	weight =  (1.f-dx) * dy * (1.f -dz);
        
	dist = S2p20(mass * weight);
	atomicAdd(&d_wd[nid + 1], dist);

	dist = S2p20(mass1 * weight);
	atomicAdd(&d_wd1[nid + 1], dist);
        
	dist = S2p20(mass2 * weight);
	atomicAdd(&d_wd2[nid + 1], dist);
      }
    
      if (isInside3D(xInt, yInt+1, zInt+1, w, h, l)){
	weight =  dx * (1.f - dy) * (1.f -dz);

	dist = S2p20(mass * weight);
	atomicAdd(&d_wd[nid + w], dist);

	dist = S2p20(mass1 * weight);
	atomicAdd(&d_wd1[nid + w], dist);

	dist = S2p20(mass2 * weight);
	atomicAdd(&d_wd2[nid + w], dist);
      }
    
      if (isInside3D(xInt+1, yInt+1, zInt+1, w, h, l)){
	weight = (1.f -dx) * (1.f - dy) * (1.f -dz);

	dist = S2p20(mass * weight);
	atomicAdd(&d_wd[nid + w + 1], dist);

	dist = S2p20(mass1 * weight);
	atomicAdd(&d_wd1[nid + w + 1], dist);

	dist = S2p20(mass2 * weight);
	atomicAdd(&d_wd2[nid + w + 1], dist);
      } 
    }

#if __CUDA_ARCH__ >= 200
    inline  __device__ void atomicSplatFloat(float* d_wd, float mass, float x, float y, float z,
					int w, int h, int l)
    {
      int xInt = int(x);
      int yInt = int(y);
      int zInt = int(z);

      if (x < 0 && x != xInt) --xInt;
      if (y < 0 && y != yInt) --yInt;
      if (z < 0 && z != zInt) --zInt;
    
      float dx = 1.f - (x - float(xInt));
      float dy = 1.f - (y - float(yInt));
      float dz = 1.f - (z - float(zInt));

      int nid = (zInt * h + yInt) * w + xInt;

      float ww;

      if (isInside3D(xInt, yInt, zInt, w, h, l)){
	ww = dx * dy * dz;
	atomicAdd(&d_wd[nid], mass*ww);
      }
            
      if (isInside3D(xInt + 1, yInt, zInt, w, h, l)){
	ww = (1.f-dx) * dy * dz;
	atomicAdd(&d_wd[nid + 1], mass*ww);
      }

      if (isInside3D(xInt, yInt + 1, zInt, w, h, l)){
	ww = dx * (1.f - dy) * dz;
	atomicAdd(&d_wd[nid + w], mass*ww);
      }
    
      if (isInside3D(xInt+1, yInt+1, zInt, w, h, l)){
	ww = (1.f -dx) * (1.f - dy) * dz;
	atomicAdd(&d_wd[nid + w + 1], mass*ww);
      } 
    
      nid += w*h;

      if (isInside3D(xInt, yInt, zInt + 1, w, h, l)){
	ww = dx * dy * (1.f - dz);
	atomicAdd(&d_wd[nid], mass*ww);
      }
            
      if (isInside3D(xInt + 1, yInt, zInt+1, w, h, l)){
	ww = (1.f-dx) * dy * (1.f -dz);
	atomicAdd(&d_wd[nid + 1], mass*ww);
      }
    
      if (isInside3D(xInt, yInt+1, zInt+1, w, h, l)){
	ww = dx * (1.f - dy) * (1.f -dz);
	atomicAdd(&d_wd[nid + w], mass*ww);
      }
    
      if (isInside3D(xInt+1, yInt+1, zInt+1, w, h, l)){
	ww = (1.f -dx) * (1.f - dy) * (1.f -dz);
	atomicAdd(&d_wd[nid + w + 1], mass*ww);
      } 
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Splat mass along with weight vector
    ///  d_wd = sum dist * d_w (weighted distance)
    ////////////////////////////////////////////////////////////////////////////////
    inline  __device__ void atomicSplatFloat(float* d_wd, float* d_ww, float mass, float x, float y, float z,
					int w, int h, int l)
    {
      int xInt = int(x);
      int yInt = int(y);
      int zInt = int(z);

      if (x < 0 && x != xInt) --xInt;
      if (y < 0 && y != yInt) --yInt;
      if (z < 0 && z != zInt) --zInt;
    
      float dx = 1.f - (x - float(xInt));
      float dy = 1.f - (y - float(yInt));
      float dz = 1.f - (z - float(zInt));

      int nid = (zInt * h + yInt) * w + xInt;

      float ww;

      if (isInside3D(xInt, yInt, zInt, w, h, l)){
	ww = dx * dy * dz;
	atomicAdd(&d_wd[nid], mass*ww);
	atomicAdd(&d_ww[nid], ww);
      }
            
      if (isInside3D(xInt + 1, yInt, zInt, w, h, l)){
	ww = (1.f-dx) * dy * dz;
	atomicAdd(&d_wd[nid + 1], mass*ww);
	atomicAdd(&d_ww[nid + 1], ww);
      }

      if (isInside3D(xInt, yInt + 1, zInt, w, h, l)){
	ww = dx * (1.f - dy) * dz;
	atomicAdd(&d_wd[nid + w], mass*ww);
	atomicAdd(&d_ww[nid + w], ww);
      }
    
      if (isInside3D(xInt+1, yInt+1, zInt, w, h, l)){
	ww = (1.f -dx) * (1.f - dy) * dz;
	atomicAdd(&d_wd[nid + w + 1], mass*ww);
	atomicAdd(&d_ww[nid + w + 1], ww);
      } 
    
      nid += w*h;

      if (isInside3D(xInt, yInt, zInt + 1, w, h, l)){
	ww = dx * dy * (1.f - dz);
	atomicAdd(&d_wd[nid], mass*ww);
	atomicAdd(&d_ww[nid], ww);
      }
            
      if (isInside3D(xInt + 1, yInt, zInt+1, w, h, l)){
	ww = (1.f-dx) * dy * (1.f -dz);
	atomicAdd(&d_wd[nid + 1], mass*ww);
	atomicAdd(&d_ww[nid + 1], ww);
      }
    
      if (isInside3D(xInt, yInt+1, zInt+1, w, h, l)){
	ww = dx * (1.f - dy) * (1.f -dz);
	atomicAdd(&d_wd[nid + w], mass*ww);
	atomicAdd(&d_ww[nid + w], ww);
      }
    
      if (isInside3D(xInt+1, yInt+1, zInt+1, w, h, l)){
	ww = (1.f -dx) * (1.f - dy) * (1.f -dz);
	atomicAdd(&d_wd[nid + w + 1], mass*ww);
	atomicAdd(&d_ww[nid + w + 1], ww);
      } 
    }
#endif

  }// end namespace Splatting
} // end namespace PyCA

#endif // __GSPLAT_CUH

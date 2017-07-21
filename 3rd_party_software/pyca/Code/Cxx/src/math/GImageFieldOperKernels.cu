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

#include "GImageFieldOperKernels.h"
#include <pycaUtils.h>
#include "GSplat.h"
#include "interp.h"
#include <GMemOpers.h>

#include "FiniteDiff.h"

// TEST make sure boost isn't included in nvcc code
#if defined(BOOST_COMPILER)
int bla[-1];
#endif

namespace PyCA {

template<BackgroundStrategy bg, InterpT interp>
__global__ void ApplyH_kernel(float* d_o, const float* d_i,
                              const float* d_hx, const float* d_hy, const float* d_hz,
                              int sizeX, int sizeY, int sizeZ)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    
    if (i < sizeX && j < sizeY){
        int id = j * sizeX + i;
        for (int k=0; k< sizeZ; ++k, id+= sizeX*sizeY){
            float x = d_hx[id];
            float y = d_hy[id];
            float z = d_hz[id];
            
            d_o[id] = point_interp<interp, bg>
		(d_i, x, y, z, sizeX, sizeY, sizeZ);
        }
    }
}

/*
 * apply hField to an image
 *  defImage(x) = image(h(x))
 */
template<BackgroundStrategy bg, InterpT interp>
void ApplyH(float* d_o, const float* d_i,
            const float* d_hx, const float* d_hy, const float* d_hz,
            int sizeX, int sizeY, int sizeZ, StreamT stream)
{
    dim3 threads(16,16);
    dim3 grids(iDivUp(sizeX, threads.x), iDivUp(sizeY, threads.y));

    ApplyH_kernel<bg, interp>
	<<< grids, threads, 0, stream >>>(d_o, d_i,
					  d_hx, d_hy, d_hz,
					  sizeX, sizeY, sizeZ);
 }


/*
 *  apply uField to an image
 *  defImage(x) = image(x + delta * u(x))
 */

__device__ __constant__ float c_delta;
#define MK_FROM_U_TO_H() \
    float x, y, z;                                          \
    if (fwd) {                                              \
        x = i + delta * iSpX * d_uX[id];                    \
        y = j + delta * iSpY * d_uY[id];                    \
        z = k + delta * iSpZ * d_uZ[id];                    \
    } else {                                                \
        x = i - delta * iSpX * d_uX[id];                    \
        y = j - delta * iSpY * d_uY[id];                    \
        z = k - delta * iSpZ * d_uZ[id];                    \
    }                                                       \

#define MK_CORE_APPLY_U()        int id = j * sizeX + i;    \
    for (int k=0; k< sizeZ; ++k, id+= sizeX*sizeY){                 \
        MK_FROM_U_TO_H();                                           \
        d_o[id] = point_interp<interp, bg>                          \
	    (d_i, x, y, z, sizeX, sizeY, sizeZ);		    \
    }

template<bool fwd, BackgroundStrategy bg, InterpT interp>
__global__ void ApplyV_kernel(float* d_o, const float* d_i, 
                              const float* d_uX, const float* d_uY, const float* d_uZ,
                              int sizeX, int sizeY, int sizeZ,
                              float iSpX, float iSpY, float iSpZ){
        
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    if (i >= sizeX || j >= sizeY)
        return;
    float delta = c_delta;
    MK_CORE_APPLY_U();
}

template<bool fwd, BackgroundStrategy bg, InterpT interp>
__global__ void ApplyV_kernel(float* d_o, const float* d_i, 
                              const float* d_uX, const float* d_uY, const float* d_uZ,
                              float delta, 
                              int sizeX, int sizeY, int sizeZ,
                              float iSpX, float iSpY, float iSpZ){
        
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    if (i >= sizeX || j >= sizeY)
        return;
    MK_CORE_APPLY_U();
}

template<bool fwd, BackgroundStrategy bg, InterpT interp>
void ApplyV(float* d_o, const float* d_i,
            const float* d_ux, const float* d_uy, const float* d_uz, const float& delta,
            int sizeX, int sizeY, int sizeZ,
            float spX, float spY, float spZ, StreamT stream, bool onDev)
{
    MK_CHECK_IMAGE_BACKGROUND(bg);
    dim3 threads(16,16);
    dim3 grids(iDivUp(sizeX, threads.x), iDivUp(sizeY, threads.y));
    if (!onDev) {
        ApplyV_kernel<fwd, bg, interp>
	    <<<grids, threads,0,stream>>>
	    (d_o, d_i,
	     d_ux, d_uy, d_uz,
	     delta,
	     sizeX, sizeY, sizeZ,
	     1.f/spX, 1.f/spY, 1.f/spZ);
    } else {
        cudaMemcpyToSymbolAsync(c_delta,&delta,sizeof(float),
                                0,cudaMemcpyDeviceToDevice,stream);
        ApplyV_kernel<fwd, bg, interp>
	    <<<grids, threads,0,stream>>>
	    (d_o, d_i,
	     d_ux, d_uy, d_uz,
	     sizeX, sizeY, sizeZ,
	     1.f/spX, 1.f/spY, 1.f/spZ);
    }
}

__device__ __constant__ float3 c_trans;


template<BackgroundStrategy bg, InterpT interp>
__global__ void ComposeTranslation_kernel(float* d_o, const float* d_i,
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
            
            d_o[id] = point_interp<interp, bg>
		(d_i, x, y, z, sizeX, sizeY, sizeZ);
        }
    }
}

template<BackgroundStrategy bg, InterpT interp>
__global__ void ComposeTranslation_const_kernel(float* d_o, const float* d_i,
                                          int sizeX, int sizeY, int sizeZ){
        
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    
    if (i < sizeX && j < sizeY){
        int id = j * sizeX + i;
        for (int k=0; k< sizeZ; ++k, id+= sizeX*sizeY){
            float x = i + c_trans.x;
            float y = j + c_trans.y;
            float z = k + c_trans.z;

            d_o[id] = point_interp<interp, bg>
		(d_i, x, y, z, sizeX, sizeY, sizeZ);
        }
    }
}

template<BackgroundStrategy bg, InterpT interp>
void ComposeTranslation(float* d_o, const float* d_i,
                        const Vec3Df& t, const Vec3Di& size, StreamT stream, bool onDev)
{
    dim3 threads(16, 16);
    dim3 grids(iDivUp(size.x, threads.x), iDivUp(size.y, threads.y));
    if (onDev) {
        cudaMemcpyToSymbolAsync(c_trans, &t.x,sizeof(float) * 3,
                                0, cudaMemcpyDeviceToDevice, stream);
        ComposeTranslation_const_kernel<bg, interp>
	    <<<grids, threads, 0, stream>>>
            (d_o, d_i, size.x, size.y, size.z);
    } else {
        ComposeTranslation_kernel<bg, interp>
	    <<<grids, threads, 0, stream>>>
            (d_o, d_i, t.x, t.y, t.z, size.x, size.y, size.z);
    }
}

void
Splat(float *d_o, 
      const float *d_hx, 
      const float *d_hy, 
      const float *d_hz, 
      const float *d_i, 
      Vec3Di sz,
      StreamT stream)
{
   size_t nVox = sz.prod();
   int* i_do =(int*)d_o;
   // Splat to fixed buffer
   Splatting::splat3D(i_do, sz.x, sz.y, sz.z,
		      d_i, d_hx, d_hy, d_hz, nVox, stream);
   Splatting::FixedToFloating_I(d_o, nVox, stream);
}

void
SplatAndNormalize(float *d_o, 
		  const float *d_hx, 
		  const float *d_hy, 
		  const float *d_hz, 
		  const float *d_i, 
		  const float *temp,
		  Vec3Di sz,
		  StreamT stream)
{
   size_t nVox = sz.prod();
   int* i_do = (int*)d_o;
   int* i_dd = (int*)temp;
   // Splat to fixed buffer with weighted distance
   Splatting::splat3D(i_do, i_dd, sz.x, sz.y, sz.z,
		      d_i, d_hx, d_hy, d_hz, nVox, stream);
   // convert back to floating point buffer with distance
   Splatting::convertWeightedDistance_I(d_o, i_dd, nVox, stream);
}

template < class T, enum DimT dim, enum DiffT diffType, 
	   enum BoundaryCondT bc,
	   bool accum, enum OpT op >
static
__global__
void
finiteDiff_kernel(T *out, 
		  const T *in, 
		  int sx, int sy, int sz,
		  float scale) // spacing of dim
{
   const int x = blockDim.x * blockIdx.x + threadIdx.x;
   const int y = blockDim.y * blockIdx.y + threadIdx.y;

   if(x < sx && y < sy){

       int index = y*sx + x; 
    
       scale = 1.f/scale;

       const int stridez = sx*sy;
       for(int z = 0; z < sz; ++z){
	   float val = finiteDiff<T, dim, diffType, bc>(in,x,y,z,sx,sy,sz);

	   val *= scale; // take care of voxel scaling

	   if(op == OP_SQR){
	       val = val*val;
	   }
	
	   if(accum)
	       out[index] += val;
	   else
	       out[index] = val;
	
	   index += stridez;
       }
   }
    
}

#define ACCUM_TRUE 1
#define ACCUM_FALSE 0
#define SLICE_TRUE 1
#define SLICE_FALSE 0

template <enum DimT dim, enum DiffT diffType, 
	  enum BoundaryCondT bc, 
	  bool accum, OpT op>
static
inline
void
g_finiteDiff(float* d_o, const float* d_i,
	     int szX, int szY, int szZ,
	     float spX, float spY, float spZ, 
	     StreamT stream)
{
   dim3 threads(16,16);
   dim3 grids(iDivUp(szX, threads.x), iDivUp(szY, threads.y));
   bool slice = (szZ == 1);

   float sp;
   if(dim == DIM_X){
      sp = spX;
   }else if(dim == DIM_Y){
      sp = spY;
   }else if(dim == DIM_Z){ 
      sp = spZ;
   }else{
      throw PyCAException(__FILE__, __LINE__, "unknown DimT");
   }

   if(slice && dim == DIM_Z){ 
      GMemOpers<float>::SetMem(d_o,0.f,szX*szY*szZ,stream,false);
   }else{
      finiteDiff_kernel<float,dim,diffType,bc,accum,op>
	 <<<grids,threads,0,stream>>>
	 (d_o, d_i, szX, szY, szZ, sp);
   }
}

template <enum BoundaryCondT bc, 
	  bool accum, OpT op>
static
inline
void
g_finiteDiff(float* d_o, const float* d_i,
	     int szX, int szY, int szZ,
	     float spX, float spY, float spZ, 
	     DimT dim, DiffT diffType, 
	     StreamT stream)
{
   if(diffType == DIFF_FORWARD){
      if(dim == DIM_X){
	 g_finiteDiff<DIM_X, DIFF_FORWARD, bc, accum, op>
	    (d_o, d_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     stream);
      }else if(dim == DIM_Y){
	 g_finiteDiff<DIM_Y, DIFF_FORWARD, bc, accum, op>
	    (d_o, d_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     stream);
      }else if(dim == DIM_Z){
	 g_finiteDiff<DIM_Z, DIFF_FORWARD, bc, accum, op>
	    (d_o, d_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     stream);
      }else{
	 throw PyCAException(__FILE__, __LINE__, "unknown DimT");
      }
   }else if(diffType == DIFF_BACKWARD){
      if(dim == DIM_X){
	 g_finiteDiff<DIM_X, DIFF_BACKWARD, bc, accum, op>
	    (d_o, d_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     stream);
      }else if(dim == DIM_Y){
	 g_finiteDiff<DIM_Y, DIFF_BACKWARD, bc, accum, op>
	    (d_o, d_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     stream);
      }else if(dim == DIM_Z){
	 g_finiteDiff<DIM_Z, DIFF_BACKWARD, bc, accum, op>
	    (d_o, d_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     stream);
      }else{
	 throw PyCAException(__FILE__, __LINE__, "unknown DimT");
      }
   }else if(diffType == DIFF_CENTRAL){
      if(dim == DIM_X){
	 g_finiteDiff<DIM_X, DIFF_CENTRAL, bc, accum, op>
	    (d_o, d_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     stream);
      }else if(dim == DIM_Y){
	 g_finiteDiff<DIM_Y, DIFF_CENTRAL, bc, accum, op>
	    (d_o, d_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     stream);
      }else if(dim == DIM_Z){
	 g_finiteDiff<DIM_Z, DIFF_CENTRAL, bc, accum, op>
	    (d_o, d_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     stream);
      }else{
	 throw PyCAException(__FILE__, __LINE__, "unknown DimT");
      }
   }else{
      throw PyCAException(__FILE__, __LINE__, "unknown DiffT");
   }
}

template <bool accum, OpT op>
static
inline
void
g_finiteDiff(float* d_o, const float* d_i,
	     int szX, int szY, int szZ,
	     float spX, float spY, float spZ, 
	     DimT dim, DiffT diffType, 
	     enum BoundaryCondT bc, 
	     StreamT stream)
{
   if(bc == BC_APPROX){
      g_finiteDiff<BC_APPROX, accum, op>
	 (d_o, d_i,
	  szX, szY, szZ,
	  spX, spY, spZ, 
	  dim, diffType,
	  stream);
   }else if(bc == BC_WRAP){
      g_finiteDiff<BC_WRAP, accum, op>
	 (d_o, d_i,
	  szX, szY, szZ,
	  spX, spY, spZ, 
	  dim, diffType,
	  stream);
   }else if(bc == BC_CLAMP){
      g_finiteDiff<BC_CLAMP, accum, op>
	 (d_o, d_i,
	  szX, szY, szZ,
	  spX, spY, spZ, 
	  dim, diffType,
	  stream);
   }else{
      throw PyCAException(__FILE__, __LINE__, 
				   "unknown boundary condition (BoundaryCondT)");
   }

}

void
g_finiteDiff(float* d_o, const float* d_i,
	     int szX, int szY, int szZ,
	     float spX, float spY, float spZ, 
	     DimT dim, DiffT diffType, 
	     enum BoundaryCondT bc, 
	     bool accum, OpT op,
	     StreamT stream)
{
   if(accum){
      if(op == OP_VAL){
	 g_finiteDiff<ACCUM_TRUE, OP_VAL>
	    (d_o, d_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     dim, diffType,
	     bc,
	     stream);
      }else if(op == OP_SQR){
	 g_finiteDiff<ACCUM_TRUE, OP_SQR>
	    (d_o, d_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     dim, diffType,
	     bc,
	     stream);
      }else{
	 throw PyCAException(__FILE__, __LINE__, 
				      "unknown OpT");
      }
   }else{
      if(op == OP_VAL){
	 g_finiteDiff<ACCUM_FALSE, OP_VAL>
	    (d_o, d_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     dim, diffType,
	     bc,
	     stream);
      }else if(op == OP_SQR){
	 g_finiteDiff<ACCUM_FALSE, OP_SQR>
	    (d_o, d_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     dim, diffType,
	     bc,
	     stream);
      }else{
	 throw PyCAException(__FILE__, __LINE__, 
				      "unknown OpT");
      }
   }
}

#define BG_CLAMP BACKGROUND_STRATEGY_CLAMP

template<DimT dim, bool slice>
static
__global__
void
upwindDiff_kernel(float *rtn,
		  const float *array,
		  const float *speed,
		  int szX, int szY, int szZ,
		  float spX, float spY, float spZ)
{
   const uint x = blockDim.x * blockIdx.x + threadIdx.x;
   const uint y = blockDim.y * blockIdx.y + threadIdx.y;

   float sp;
   if(dim == DIM_X){
      sp = spX;
   }else if(dim == DIM_Y){
      sp = spY;
   }else if(dim == DIM_Z){
      sp = spZ;
   }

   if(x < szX && y < szY){

      const uint stridez = szX*szY;
  
      uint index = y * szX + x;
      for(uint z = 0; z < szZ; z++){
	 float v = array[index]; // val
	 // previous and next values
	 float vp, vn;
	 if(dim == DIM_X){
	    vp = getSafeVal<float,BG_CLAMP>(array, szX,szY,szZ, x-1,y,z);
	    vn = getSafeVal<float,BG_CLAMP>(array, szX,szY,szZ, x+1,y,z);
	 }else if(dim == DIM_Y){
	    vp = getSafeVal<float,BG_CLAMP>(array, szX,szY,szZ, x,y-1,z);
	    vn = getSafeVal<float,BG_CLAMP>(array, szX,szY,szZ, x,y+1,z);
	 }else if(dim == DIM_Z){
	    if(slice){
	       vp = v;
	       vn = v;
	    }else{
	       vp = getSafeVal<float,BG_CLAMP>(array, szX,szY,szZ, x,y,z-1);
	       vn = getSafeVal<float,BG_CLAMP>(array, szX,szY,szZ, x,y,z+1);
	    }
	 }
     
	 float spd = speed[index];

	 float dx = 0.f;
	 if ( spd < 0.0f){
	    // forward difference
	    dx = (vn - v)/sp;
	 }else{
	    dx = (v - vp)/sp;
	 }
	 rtn[index] = dx;
	 index += stridez;
      }
   }
}

void 
UpwindDiff(float *d_o, const float *d_i, 
	   const float *d_speed,
	   const Vec3Di& sz, 
	   const Vec3Df& sp,
	   DimT dim,
	   StreamT stream)
{
   dim3 threads(16,16);
   dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));
   bool slice = (sz.z == 1);
   
   if(slice){
      if(dim == DIM_X){
	  upwindDiff_kernel<DIM_X, SLICE_TRUE><<<grids,threads,0,stream>>>
	      (d_o, d_i, d_speed,
	       sz.x, sz.y, sz.z, 
	       sp.x, sp.y, sp.z);
      }else if(dim == DIM_Y){
	 upwindDiff_kernel<DIM_Y, SLICE_TRUE><<<grids,threads,0,stream>>>
	     (d_o, d_i, d_speed,
	      sz.x, sz.y, sz.z, 
	      sp.x, sp.y, sp.z);
      }else if(dim == DIM_Z){
	  upwindDiff_kernel<DIM_Z, SLICE_TRUE><<<grids,threads,0,stream>>>
	      (d_o, d_i, d_speed,
	       sz.x, sz.y, sz.z, 
	       sp.x, sp.y, sp.z);
      }else{
	  throw PyCAException(__FILE__, __LINE__, "unknown DimT");
      }
   }else{
       if(dim == DIM_X){
	   upwindDiff_kernel<DIM_X, SLICE_FALSE><<<grids,threads,0,stream>>>
	       (d_o, d_i, d_speed,
		sz.x, sz.y, sz.z, 
		sp.x, sp.y, sp.z);
       }else if(dim == DIM_Y){
	   upwindDiff_kernel<DIM_Y, SLICE_FALSE><<<grids,threads,0,stream>>>
	       (d_o, d_i, d_speed,
		sz.x, sz.y, sz.z, 
		sp.x, sp.y, sp.z);
       }else if(dim == DIM_Z){
	   upwindDiff_kernel<DIM_Z, SLICE_FALSE><<<grids,threads,0,stream>>>
	       (d_o, d_i, d_speed,
		sz.x, sz.y, sz.z, 
		sp.x, sp.y, sp.z);
       }else{
	   throw PyCAException(__FILE__, __LINE__, "unknown DimT");
       }
   }
}

template<bool slice>
static
__global__
void
upwindGradMag_kernel(float *rtn,
		     const float *array,
		     const float *speed,
		     int szX, int szY, int szZ,
		     float spX, float spY, float spZ)
{
   const uint x = blockDim.x * blockIdx.x + threadIdx.x;
   const uint y = blockDim.y * blockIdx.y + threadIdx.y;

   if(x < szX && y < szY){

      const uint stridez = szX*szY;
  
      uint index = y * szX + x;
      for(uint z = 0; z < szZ; z++){
	 float v = array[index]; // val
	 float n = getSafeVal<float,BG_CLAMP>(array, szX,szY,szZ, x,y-1,z); // north
	 float s = getSafeVal<float,BG_CLAMP>(array, szX,szY,szZ, x, y+1, z); // south
	 float e = getSafeVal<float,BG_CLAMP>(array, szX,szY,szZ, x+1, y, z); // east
	 float w = getSafeVal<float,BG_CLAMP>(array, szX,szY,szZ, x-1, y, z); // west
	 float u = 0.f;
	 float d = 0.f;
	 if(!slice){
	    u = getSafeVal<float,BG_CLAMP>(array, szX,szY,szZ, x, y, z+1); // up
	    d = getSafeVal<float,BG_CLAMP>(array, szX,szY,szZ, x, y, z-1); // down
	 }
     
	 float grad_before_x = (x < szX-1 ? e - v : v - w)/spX;
	 float grad_after_x = (x > 0 ? v - w : e - v)/spX;
	 float grad_before_y = (y < szY-1 ? s - v : v - n)/spY;
	 float grad_after_y = (y > 0 ? v - n : s - v)/spY;
	 float grad_before_z = 0.f;
	 float grad_after_z = 0.f;
	 if(!slice){
	    grad_before_z = (index < stridez*(szZ-1) ? u - v : v - d)/spZ;
	    grad_after_z = (index >= stridez ? v - d : u - v)/spZ;
	 }
     
	 float spd = speed[index];
     
	 if ( spd < 0.0f)
	    {
	       grad_before_x = PYCAMIN(grad_before_x, 0.0f);
	       grad_after_x = PYCAMIN(-grad_after_x, 0.0f);
	       grad_before_y = PYCAMIN(grad_before_y, 0.0f);
	       grad_after_y = PYCAMIN(-grad_after_y, 0.0f);
	       if(!slice){
		  grad_before_z = PYCAMIN(grad_before_z, 0.0f);
		  grad_after_z = PYCAMIN(-grad_after_z, 0.0f);
	       }
	    }
	 else
	    {
	       grad_before_x = PYCAMAX(grad_before_x, 0.0f);
	       grad_after_x = PYCAMAX(-grad_after_x, 0.0f);
	       grad_before_y = PYCAMAX(grad_before_y, 0.0f);
	       grad_after_y = PYCAMAX(-grad_after_y, 0.0f);
	       if(!slice){
		  grad_before_z = PYCAMAX(grad_before_z, 0.0f);
		  grad_after_z = PYCAMAX(-grad_after_z, 0.0f);
	       }
	    }
	 float gradmag = 0.f;
	 if(!slice){
	    gradmag = 
	       sqrt(grad_after_x*grad_after_x + grad_before_x*grad_before_x +
		    grad_after_y*grad_after_y + grad_before_y*grad_before_y +
		    grad_after_z*grad_after_z + grad_before_z*grad_before_z);
	 }else{
	    gradmag = 
	       sqrt(grad_after_x*grad_after_x + grad_before_x*grad_before_x +
		    grad_after_y*grad_after_y + grad_before_y*grad_before_y);
	 }
    
	 rtn[index] = gradmag*spd;
	 index += stridez;
      }
   }
}

void 
UpwindGradMag(float *d_o, 
	      const float *d_i, 
	      const float *d_speed, 
	      const Vec3Di& sz,
	      const Vec3Df& sp,
	      StreamT stream)
{
   dim3 threads(16,16);
   dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));
   bool slice = (sz.z == 1);

   if(slice){
      upwindGradMag_kernel<SLICE_TRUE><<<grids,threads,0,stream>>>
	 (d_o, d_i, d_speed,
	  sz.x, sz.y, sz.z, 
	  sp.x, sp.y, sp.z);
   }else{
      upwindGradMag_kernel<SLICE_FALSE><<<grids,threads,0,stream>>>
	 (d_o, d_i, d_speed,
	  sz.x, sz.y, sz.z, 
	  sp.x, sp.y, sp.z);
   }

}

template <enum DiffT diffType, enum BoundaryCondT bc>
inline
void 
g_gradient(float* d_ox, float *d_oy, float* d_oz,
	   const float* d_i,
	   int szX, int szY, int szZ,
	   float spX, float spY, float spZ, StreamT stream)
{
    dim3 threads(16,16);
    dim3 grids(iDivUp(szX, threads.x), iDivUp(szY, threads.y));
    bool slice = (szZ == 1);

    finiteDiff_kernel<float,DIM_X,diffType,bc,
		      ACCUM_FALSE,OP_VAL>
       <<<grids,threads,0,stream>>>
       (d_ox, d_i, szX, szY, szZ, spX);
    
    finiteDiff_kernel<float,DIM_Y,diffType,bc,
		      ACCUM_FALSE,OP_VAL>
       <<<grids,threads,0,stream>>>
       (d_oy, d_i, szX, szY, szZ, spY);
    
    if(slice){
       GMemOpers<float>::SetMem(d_oz,0.f,szX*szY*szZ,stream,false);
    }else{
       finiteDiff_kernel<float,DIM_Z,diffType,bc,
			 ACCUM_FALSE,OP_VAL>
	  <<<grids,threads,0,stream>>>
	  (d_oz, d_i, szX, szY, szZ, spZ);
    }
}

// version of gradient taking boundary condition as a parameter
template <enum DiffT diffType>
inline
void 
g_gradient(float* d_ox, float *d_oy, float* d_oz,
	   const float* d_i,
	   int szX, int szY, int szZ,
	   float spX, float spY, float spZ,
	   BoundaryCondT bc, StreamT stream)
{
   if(bc == BC_APPROX){
      g_gradient<diffType,BC_APPROX>
	 (d_ox,d_oy,d_oz,d_i,szX,szY,szZ,spX,spY,spZ,stream);
   }else if(bc == BC_WRAP){
      g_gradient<diffType,BC_WRAP>
	 (d_ox,d_oy,d_oz,d_i,szX,szY,szZ,spX,spY,spZ,stream);
   }else if(bc == BC_CLAMP){
      g_gradient<diffType,BC_CLAMP>
	 (d_ox,d_oy,d_oz,d_i,szX,szY,szZ,spX,spY,spZ,stream);
   }else{
      throw PyCAException(__FILE__, __LINE__, "unknown BoundaryCondT");
   }
}

template < class T, enum DimT dim, enum DiffT diffType, 
	   enum BoundaryCondT bc,
	   bool accum, enum OpT op >
static
__global__
void
finiteDiffMask_kernel(T *out, 
		      const T *in, const T *mask,
		      int sx, int sy, int sz,
		      float scale) // spacing of dim
{
   const int x = blockDim.x * blockIdx.x + threadIdx.x;
   const int y = blockDim.y * blockIdx.y + threadIdx.y;

   if(x < sx && y < sy){

       int index = y*sx + x; 
    
       scale = 1.f/scale;

       const int stridez = sx*sy;
       for(int z = 0; z < sz; ++z){
	   float val = 
	       finiteDiffMask<T, dim, diffType, bc>
	       (in,mask,x,y,z,sx,sy,sz);

	   val *= scale; // take care of voxel scaling

	   if(op == OP_SQR){
	       val = val*val;
	   }
	
	   if(accum)
	       out[index] += val;
	   else
	       out[index] = val;
	
	   index += stridez;
       }
   }
    
}

template <enum DiffT diffType, enum BoundaryCondT bc>
inline
void 
g_gradientMask(float* d_ox, float *d_oy, float* d_oz,
	       const float* d_i, const float* d_mask,
	       int szX, int szY, int szZ,
	       float spX, float spY, float spZ, 
	       StreamT stream)
{
   dim3 threads(16,16);
   dim3 grids(iDivUp(szX, threads.x), iDivUp(szY, threads.y));
   bool slice = (szZ == 1);
   
   finiteDiffMask_kernel<float,DIM_X,diffType,bc,
			 ACCUM_FALSE,OP_VAL>
      <<<grids,threads,0,stream>>>
      (d_ox, d_i, d_mask, szX, szY, szZ, spX);
   
   finiteDiffMask_kernel<float,DIM_Y,diffType,bc,
			 ACCUM_FALSE,OP_VAL>
      <<<grids,threads,0,stream>>>
      (d_oy, d_i, d_mask, szX, szY, szZ, spY);
   
   if(slice){
      GMemOpers<float>::SetMem(d_oz,0.f,szX*szY*szZ,stream,false);
   }else{
      finiteDiffMask_kernel<float,DIM_Z,diffType,bc,
			    ACCUM_FALSE,OP_VAL>
	 <<<grids,threads,0,stream>>>
	 (d_oz, d_i, d_mask, szX, szY, szZ, spZ);
   }
}

// version of gradient taking boundary condition as a parameter
template <enum DiffT diffType>
inline
void 
g_gradientMask(float* d_ox, float *d_oy, float* d_oz,
	       const float* d_i, const float* d_mask,
	       int szX, int szY, int szZ,
	       float spX, float spY, float spZ,
	       BoundaryCondT bc, StreamT stream)
{
   if(bc == BC_APPROX){
      g_gradientMask<diffType,BC_APPROX>
	 (d_ox,d_oy,d_oz,d_i,d_mask,szX,szY,szZ,spX,spY,spZ,stream);
   }else if(bc == BC_WRAP){
      throw PyCAException(__FILE__, __LINE__, "BC_WRAP boundary condition unimplemented for masked gradient");
   }else if(bc == BC_CLAMP){
      g_gradientMask<diffType,BC_CLAMP>
	 (d_ox,d_oy,d_oz,d_i,d_mask,szX,szY,szZ,spX,spY,spZ,stream);
   }else{
      throw PyCAException(__FILE__, __LINE__, "unknown BoundaryCondT");
   }
}

template<typename T, enum DiffT diffType, enum BoundaryCondT bc>
void g_gradientMag(float* d_o,
                 const T* d_i,
                 int sizeX, int sizeY, int sizeZ,
                 float spX, float spY, float spZ, StreamT stream)
{
   dim3 threads(16,16);
   dim3 grids(iDivUp(sizeX, threads.x), iDivUp(sizeY, threads.y));
   bool slice = (sizeZ == 1);

   finiteDiff_kernel<T,DIM_X,diffType,bc,ACCUM_FALSE,OP_SQR>
      <<<grids,threads,0,stream>>>
      (d_o, d_i, sizeX, sizeY, sizeZ, spX);
   
   finiteDiff_kernel<T,DIM_Y,diffType,bc,ACCUM_TRUE,OP_SQR>
      <<<grids,threads,0,stream>>>
      (d_o, d_i, sizeX, sizeY, sizeZ, spY);
   
   if(!slice){
      finiteDiff_kernel<T,DIM_Z,diffType,bc,ACCUM_TRUE,OP_SQR>
	 <<<grids,threads,0,stream>>>
	 (d_o, d_i, sizeX, sizeY, sizeZ, spZ);
   }

   GMemOpers<float>::Sqrt_I(d_o,sizeX*sizeY*sizeZ,stream);
}

template<typename T, enum DiffT diffType>
void g_gradientMag(float* d_o,
		   const T* d_i,
		   int sizeX, int sizeY, int sizeZ,
		   float spX, float spY, float spZ, 
		   BoundaryCondT bc, StreamT stream)
{
   if(bc == BC_APPROX){
      PyCA::g_gradientMag<T, diffType, BC_APPROX>
	 (d_o, d_i, sizeX, sizeY, sizeZ, spX, spY, spZ, stream);
   }else if(bc == BC_WRAP){
      PyCA::g_gradientMag<T, diffType, BC_WRAP>
	 (d_o, d_i, sizeX, sizeY, sizeZ, spX, spY, spZ, stream);
   }else if(bc == BC_CLAMP){
      PyCA::g_gradientMag<T, diffType, BC_CLAMP>
	 (d_o, d_i, sizeX, sizeY, sizeZ, spX, spY, spZ, stream);
   }else{
      throw PyCAException(__FILE__, __LINE__, "Unknown BoundaryCondT");
   }
}

template <enum DiffT diffType, enum BoundaryCondT bc>
inline
void 
g_divergence(float* d_o,
	   const float* d_ix, const float* d_iy, const float* d_iz,
	   int   sizeX, int sizeY, int sizeZ,
	   float spX, float spY, float spZ, StreamT stream)
{
   dim3 threads(16,16);
   dim3 grids(iDivUp(sizeX, threads.x), iDivUp(sizeY, threads.y));
   bool slice = (sizeZ == 1);

   finiteDiff_kernel<float,DIM_X,diffType,bc,ACCUM_FALSE,OP_VAL>
      <<<grids,threads,0,stream>>>
      (d_o, d_ix, sizeX, sizeY, sizeZ, spX);
   
   finiteDiff_kernel<float,DIM_Y,diffType,bc,ACCUM_TRUE,OP_VAL>
      <<<grids,threads,0,stream>>>
      (d_o, d_iy, sizeX, sizeY, sizeZ, spY);
   
   if(!slice){
      finiteDiff_kernel<float,DIM_Z,diffType,bc,ACCUM_TRUE,OP_VAL>
	 <<<grids,threads,0,stream>>>
	 (d_o, d_iz, sizeX, sizeY, sizeZ, spZ);
   }

}

template <enum DiffT diffType>
inline
void 
g_divergence(float* d_o,
	     const float* d_ix, const float* d_iy, const float* d_iz,
	     int   sizeX, int sizeY, int sizeZ,
	     float spX, float spY, float spZ,
	     BoundaryCondT bc, StreamT stream)
{
   if (bc == BC_APPROX) {
      g_divergence<diffType,BC_APPROX>
	 (d_o,
	  d_ix, d_iy, d_iz,
	  sizeX, sizeY, sizeZ,
	  spX, spY, spZ, stream);
   } else if (bc == BC_WRAP) {
      g_divergence<diffType,BC_WRAP>
	 (d_o,
	  d_ix, d_iy, d_iz,
	  sizeX, sizeY, sizeZ,
	  spX, spY, spZ, stream);
   } else if (bc == BC_CLAMP) {
      g_divergence<diffType,BC_CLAMP>
	 (d_o,
	  d_ix, d_iy, d_iz,
	  sizeX, sizeY, sizeZ,
	  spX, spY, spZ, stream);
   } else {
      throw PyCAException(__FILE__, __LINE__, "Unknownd BoundaryCondT");
   }
   
}

/**
 * Compute the magnitude image
 * d_o[i] = sqrt(d_i[i].x^2 + d_i[i].y^2 + d_i[i].z^2) 
 */

template<bool square_root>
__global__ void Magnitude_kernel(float* d_o, const float* d_ix, const float* d_iy, const float* d_iz, size_t n){
    uint blockId = get_blockID();
    uint id   = get_threadID(blockId);
    
    if (id < n){
        float v = d_ix[id] * d_ix[id] + d_iy[id] * d_iy[id] + d_iz[id] * d_iz[id];
        d_o[id] = (square_root) ? sqrt(v) : v;
    }
}

/**
 * Compute the magnitude array
 * d_o[i] = sqrt(d_i[i].x^2 + d_i[i].y^2 + d_i[i].z^2)
 */
void 
Magnitude(float *d_o, 
	  const float *d_ix, 
	  const float *d_iy, 
	  const float *d_iz,
	  size_t n,
	  StreamT stream)
{
    dim3 threads(256);
    dim3 grids=make_grid(iDivUp(n, threads.x));
    Magnitude_kernel<true><<<grids, threads, 0, stream>>>
	(d_o, d_ix, d_iy, d_iz, n);
}

/**
 * Compute the magnitude array
 * d_o[i] = d_i[i].x^2 + d_i[i].y^2 + d_i[i].z^2 
 */
void 
SqrMagnitude(float *d_o, 
	     const float *d_ix, 
	     const float *d_iy, 
	     const float *d_iz, 
	     size_t n,
	     StreamT stream)
{
    dim3 threads(256);
    dim3 grids=make_grid(iDivUp(n, threads.x));
    Magnitude_kernel<false><<<grids, threads, 0, stream>>>
	(d_o, d_ix, d_iy, d_iz, n);
}

/**
 * Compute dot product array 
 * d_o[i] = d_i[i].x * d_i1[i].x + d_i[i].y * d_i1[i].y + d_i[i].z * d_i1[i].z
 */

__global__ void ComponentDotProd_kernel(float* d_o,
                                        const float* d_i_x, const float* d_i_y, const float* d_i_z,
                                        const float* d_i1_x, const float* d_i1_y, const float* d_i1_z, uint n)
{
    uint blockId = get_blockID();
    uint id      = get_threadID(blockId);
    
    if (id < n){
        d_o[id] = d_i_x[id]*d_i1_x[id] + 
                  d_i_y[id]*d_i1_y[id] +
                  d_i_z[id]*d_i1_z[id];
    }
}

void
ComponentDotProd(float *d_o, 
		 const float *d_ix, 
		 const float *d_iy, 
		 const float *d_iz, 
		 const float *d_i1x, 
		 const float *d_i1y, 
		 const float *d_i1z, 
		 size_t n,
		 StreamT stream)
{
    dim3 threads(256);
    dim3 grids=make_grid(iDivUp(n, threads.x));
    
    ComponentDotProd_kernel<<<grids, threads, 0, stream>>>
	(d_o,
	 d_ix, d_iy, d_iz, 
	 d_i1x, d_i1y, d_i1z, 
	 n);
}

__global__ void Add_kernel(float* d_o_x, float* d_o_y, float* d_o_z,
                           const float* d_i_x, const float* d_i_y, const float* d_i_z,
                           const float* d_i1 , uint n){
    uint blockId = get_blockID();
    uint  id = get_threadID(blockId);
    if (id < n){
        float v            = d_i1[id];
        d_o_x[id          ]= d_i_x[id             ] + v;
        d_o_y[id          ]= d_i_y[id             ] + v;
        d_o_z[id          ]= d_i_z[id             ] + v;
    }
}

/** @brief d_o.x = d_i.x + d_i1, d_o.y = d_i.y + d_i1, d_o.z = d_i.z + d_i1 */
void 
Add(float *d_ox, 
    float *d_oy, 
    float *d_oz, 
    const float *d_ix, 
    const float *d_iy, 
    const float *d_iz, 
    const float *d_i1, 
    size_t n,
    StreamT stream)
{
    dim3 threads(256);
    dim3 grids=make_grid(iDivUp(n, threads.x));
    Add_kernel<<<grids, threads, 0, stream>>>(d_ox, d_oy, d_oz,
                                              d_ix, d_iy, d_iz,
                                              d_i1, n);
}

__global__ void Sub_kernel(float* d_o_x, float* d_o_y, float* d_o_z,
                           const float* d_i_x, const float* d_i_y, const float* d_i_z,
                           const float* d_i1 , uint n){
    uint blockId = get_blockID();
    uint  id = get_threadID(blockId);
    if (id < n){
        float v            = d_i1[id];
        d_o_x[id          ]= d_i_x[id             ] - v;
        d_o_y[id          ]= d_i_y[id             ] - v;
        d_o_z[id          ]= d_i_z[id             ] - v;
    }
}

/** @brief d_o.x = d_i.x - d_i1, d_o.y = d_i.y - d_i1, d_o.z = d_i.z - d_i1 */
void 
Sub(float *d_ox, 
    float *d_oy, 
    float *d_oz, 
    const float *d_ix, 
    const float *d_iy, 
    const float *d_iz, 
    const float *d_i1, 
    size_t n,
    StreamT stream)
{
    dim3 threads(256);
    dim3 grids=make_grid(iDivUp(n, threads.x));
    Sub_kernel<<<grids, threads, 0, stream>>>(d_ox, d_oy, d_oz,
                                              d_ix, d_iy, d_iz,
                                              d_i1, n);
}

__global__ void Mul_kernel(float* d_o_x, float* d_o_y, float* d_o_z,
                           const float* d_i_x, const float* d_i_y, const float* d_i_z,
                           const float* d_i1 , uint n){
    uint blockId = get_blockID();
    uint  id = get_threadID(blockId);
    if (id < n){
        float v            = d_i1[id];
        d_o_x[id          ]= d_i_x[id             ] * v;
        d_o_y[id          ]= d_i_y[id             ] * v;
        d_o_z[id          ]= d_i_z[id             ] * v;
    }
}

/** @brief d_o.x = d_i.x * d_i1, d_o.y = d_i.y * d_i1, d_o.z = d_i.z * d_i1 */
void 
Mul(float *d_ox, 
    float *d_oy, 
    float *d_oz, 
    const float *d_ix, 
    const float *d_iy, 
    const float *d_iz, 
    const float *d_i1, 
    size_t n,
    StreamT stream)
{
    dim3 threads(256);
    dim3 grids=make_grid(iDivUp(n, threads.x));
    Mul_kernel<<<grids, threads, 0, stream>>>(d_ox, d_oy, d_oz,
                                              d_ix, d_iy, d_iz,
                                              d_i1, n);
}

__global__ void Div_kernel(float* d_o_x, float* d_o_y, float* d_o_z,
                           const float* d_i_x, const float* d_i_y, const float* d_i_z,
                           const float* d_i1 , uint n){
    uint blockId = get_blockID();
    uint  id = get_threadID(blockId);
    if (id < n){
        float v            = d_i1[id];
        d_o_x[id          ]= d_i_x[id             ] / v;
        d_o_y[id          ]= d_i_y[id             ] / v;
        d_o_z[id          ]= d_i_z[id             ] / v;
    }
}

/** @brief d_o.x = d_i.x / d_i1, d_o.y = d_i.y / d_i1, d_o.z = d_i.z / d_i1 */
void 
Div(float *d_ox, 
    float *d_oy, 
    float *d_oz, 
    const float *d_ix, 
    const float *d_iy, 
    const float *d_iz, 
    const float *d_i1, 
    size_t n,
    StreamT stream)
{
    dim3 threads(256);
    dim3 grids=make_grid(iDivUp(n, threads.x));
    Div_kernel<<<grids, threads, 0, stream>>>(d_ox, d_oy, d_oz,
                                              d_ix, d_iy, d_iz,
                                              d_i1, n);
}

__global__ void Add_I_kernel(float* d_o_x, float* d_o_y, float* d_o_z,
                             const float* d_i  , uint n){
    uint blockId = get_blockID();
    uint  id = get_threadID(blockId);
    if (id < n){
        float v = d_i[id];
        d_o_x[id] += v;
        d_o_y[id] += v;
        d_o_z[id] += v;
    }
}

/** @brief d_o.x += d_i, d_o.y += d_i, d_o.y += d_i, */
void
Add_I(float *d_ox, 
      float *d_oy, 
      float *d_oz, 
      const float *d_i, 
      size_t n,
      StreamT stream)
{
    dim3 threads(256);
    dim3 grids=make_grid(iDivUp(n, threads.x));
    Add_I_kernel<<<grids, threads, 0, stream>>>(d_ox, d_oy, d_oz,
                                                d_i, n);
}

__global__ void Sub_I_kernel(float* d_o_x, float* d_o_y, float* d_o_z,
                             const float* d_i  , uint n){
    uint blockId = get_blockID();
    uint  id = get_threadID(blockId);
    if (id < n){
        float v = d_i[id];
        d_o_x[id] -= v;
        d_o_y[id] -= v;
        d_o_z[id] -= v;
    }
}

/** @brief d_o.x -= d_i, d_o.y -= d_i, d_o.y -= d_i, */
void
Sub_I(float *d_ox, 
      float *d_oy, 
      float *d_oz, 
      const float *d_i, 
      size_t n,
      StreamT stream)
{
    dim3 threads(256);
    dim3 grids=make_grid(iDivUp(n, threads.x));
    Sub_I_kernel<<<grids, threads, 0, stream>>>(d_ox, d_oy, d_oz,
                                                d_i, n);
}

/** @brief d_o.x *= d_i, d_o.y *= d_i, d_o.y *= d_i, */
__global__ void Mul_I_kernel(float* d_o_x, float* d_o_y, float* d_o_z,
                             const float* d_i, uint n){
    uint blockId = get_blockID();
    uint  id = get_threadID(blockId);
    if (id < n){
        float v = d_i[id];
        d_o_x[id] *= v;
        d_o_y[id] *= v;
        d_o_z[id] *= v;
    }
}

void
Mul_I(float *d_ox, 
      float *d_oy, 
      float *d_oz, 
      const float *d_i, 
      size_t n,
      StreamT stream)
{
    dim3 threads(256);
    dim3 grids=make_grid(iDivUp(n, threads.x));
    Mul_I_kernel<<<grids, threads, 0, stream>>>(d_ox, d_oy, d_oz,
                                                d_i, n);
}

/** @brief d_o.x /= d_i, d_o.y /= d_i, d_o.y /= d_i, */
__global__ void Div_I_kernel(float* d_o_x, float* d_o_y, float* d_o_z,
                           const float* d_i, uint n){
    uint blockId = get_blockID();
    uint  id = get_threadID(blockId);
    if (id < n){
        float v= 1.f / d_i[id];
        d_o_x[id] *= v;
        d_o_y[id] *= v;
        d_o_z[id] *= v;
    }
}

void
Div_I(float *d_ox, 
      float *d_oy, 
      float *d_oz, 
      const float *d_i, 
      size_t n,
      StreamT stream)
{
    dim3 threads(256);
    dim3 grids=make_grid(iDivUp(n, threads.x));
    Div_I_kernel<<<grids, threads, 0, stream>>>(d_ox, d_oy, d_oz,
                                                d_i, n);
}

__global__ void Add_Mul_kernel(float* d_o_x, float* d_o_y, float* d_o_z,
                               const float* d_i_x, const float* d_i_y, const float* d_i_z,
                               const float* d_i1_x, const float* d_i1_y, const float* d_i1_z,
                               const float* d_i2, uint n){
    uint blockId = get_blockID();
    uint  id = get_threadID(blockId);
    if (id < n){
        float v = d_i2[id];
        d_o_x[id] = d_i_x[id]  + d_i1_x[id] * v;
        d_o_y[id] = d_i_y[id]  + d_i1_y[id] * v;
        d_o_z[id] = d_i_z[id]  + d_i1_z[id] * v;
    }
}

/** @brief d_o = d_i + d_i1 * d_i2 */
void
Add_Mul(float *d_ox, 
	float *d_oy, 
	float *d_oz, 
	const float *d_ix, 
	const float *d_iy, 
	const float *d_iz, 
	const float *d_i1x, 
	const float *d_i1y, 
	const float *d_i1z, 
	const float *d_i2, 
	size_t n,
	StreamT stream)
{
    dim3 threads(256);
    dim3 grids=make_grid(iDivUp(n, threads.x));
    Add_Mul_kernel<<<grids, threads, 0, stream>>>
	(d_ox, d_oy, d_oz,
	 d_ix, d_iy, d_iz,
	 d_i1x, d_i1y, d_i1z,
	 d_i2, n);
}

/** @brief d_o = d_o + d_i * d_i1 */
__global__ void Add_Mul_I_kernel(float* d_o_x, float* d_o_y, float* d_o_z,
                               const float* d_i_x, const float* d_i_y, const float* d_i_z,
                               const float* d_i1, uint n){
    uint blockId = get_blockID();
    uint  id = get_threadID(blockId);
    if (id < n){
        float v = d_i1[id];
        d_o_x[id] += d_i_x[id] * v;
        d_o_y[id] += d_i_y[id] * v;
        d_o_z[id] += d_i_z[id] * v;
    }
}

void
Add_Mul_I(float *d_ox, 
	  float *d_oy, 
	  float *d_oz, 
	  const float *d_ix, 
	  const float *d_iy, 
	  const float *d_iz, 
	  const float *d_i1, 
	  size_t n,
	  StreamT stream)
{
    dim3 threads(256);
    dim3 grids=make_grid(iDivUp(n, threads.x));
    Add_Mul_I_kernel<<<grids, threads, 0, stream>>>
	(d_ox, d_oy, d_oz,
	 d_ix, d_iy, d_iz,
	 d_i1, n);
}

__global__ void Sub_Mul_kernel(float* d_o_x, float* d_o_y, float* d_o_z,
                               const float* d_i_x, const float* d_i_y, const float* d_i_z,
                               const float* d_i1_x, const float* d_i1_y, const float* d_i1_z,
                               const float* d_i2, uint n){
    uint blockId = get_blockID();
    uint  id = get_threadID(blockId);
    if (id < n){
        float v = d_i2[id];
        d_o_x[id] = d_i_x[id]  - d_i1_x[id] * v;
        d_o_y[id] = d_i_y[id]  - d_i1_y[id] * v;
        d_o_z[id] = d_i_z[id]  - d_i1_z[id] * v;
    }
}

/** @brief d_o = d_i + d_i1 * d_i2 */
void
Sub_Mul(float *d_ox, 
	float *d_oy, 
	float *d_oz, 
	const float *d_ix, 
	const float *d_iy, 
	const float *d_iz, 
	const float *d_i1x, 
	const float *d_i1y, 
	const float *d_i1z, 
	const float *d_i2, 
	size_t n,
	StreamT stream)
{
    dim3 threads(256);
    dim3 grids=make_grid(iDivUp(n, threads.x));
    Sub_Mul_kernel<<<grids, threads, 0, stream>>>
	(d_ox, d_oy, d_oz,
	 d_ix, d_iy, d_iz,
	 d_i1x, d_i1y, d_i1z,
	 d_i2, n);
}

/** @brief d_o = d_o - d_i * d_i1 */
__global__ void Sub_Mul_I_kernel(float* d_o_x, float* d_o_y, float* d_o_z,
                               const float* d_i_x, const float* d_i_y, const float* d_i_z,
                               const float* d_i1, uint n){
    uint blockId = get_blockID();
    uint  id = get_threadID(blockId);
    if (id < n){
        float v = d_i1[id];
        d_o_x[id] -= d_i_x[id] * v;
        d_o_y[id] -= d_i_y[id] * v;
        d_o_z[id] -= d_i_z[id] * v;
    }
}

void
Sub_Mul_I(float *d_ox, 
	  float *d_oy, 
	  float *d_oz, 
	  const float *d_ix, 
	  const float *d_iy, 
	  const float *d_iz, 
	  const float *d_i1, 
	  size_t n,
	  StreamT stream)
{
    dim3 threads(256);
    dim3 grids=make_grid(iDivUp(n, threads.x));
    Sub_Mul_I_kernel<<<grids, threads, 0, stream>>>
	(d_ox, d_oy, d_oz,
	 d_ix, d_iy, d_iz,
	 d_i1, n);
}

__global__ void MulMulC_kernel(float* d_o_x, float* d_o_y, float* d_o_z,
                               const float* d_i_x, const float* d_i_y, const float* d_i_z,
                               const float* d_i1, float c, uint n){
    uint blockId = get_blockID();
    uint  id = get_threadID(blockId);
    if (id < n){
        float v = d_i1[id] * c;
        d_o_x[id] = d_i_x[id]  * v;
        d_o_y[id] = d_i_y[id]  * v;
        d_o_z[id] = d_i_z[id]  * v;
    }
}

__global__ void MulMulC_const_kernel(float* d_o_x, float* d_o_y, float* d_o_z,
                               const float* d_i_x, const float* d_i_y, const float* d_i_z,
                               const float* d_i1, uint n){
    uint blockId = get_blockID();
    uint  id     = get_threadID(blockId);

    float      c = c_delta;
    if (id < n){
        float v = d_i1[id] * c;
        d_o_x[id] = d_i_x[id]  * v;
        d_o_y[id] = d_i_y[id]  * v;
        d_o_z[id] = d_i_z[id]  * v;
    }
}

/** @brief d_o = d_i * d_i1 * c (d_o.x = d_i.x * d_i1 * c)*/
void
MulMulC(float *d_ox, 
	float *d_oy, 
	float *d_oz, 
	const float *d_ix, 
	const float *d_iy, 
	const float *d_iz, 
	const float *d_i1, 
	const float& c, 
	size_t n,
	StreamT stream, bool onDev)
{
    dim3 threads(256);
    dim3 grids=make_grid(iDivUp(n, threads.x));
    if (onDev) {
        cudaMemcpyToSymbolAsync(c_delta,&c, sizeof(float),
                                0,cudaMemcpyDeviceToDevice,stream);
        MulMulC_const_kernel<<<grids, threads, 0, stream>>>
	    (d_ox, d_oy, d_oz,
	     d_ix, d_iy, d_iz,
	     d_i1, n);
    } else {
        MulMulC_kernel<<<grids, threads, 0, stream>>>
	    (d_ox, d_oy, d_oz,
	     d_ix, d_iy, d_iz,
	     d_i1, c, n);
    }
}

__global__ void MulMulC_I_kernel(float* d_o_x, float* d_o_y, float* d_o_z,
                                 const float* d_i, float c, uint n){
    uint blockId = get_blockID();
    uint  id = get_threadID(blockId);
    if (id < n){
        float v = d_i[id] * c;
        d_o_x[id] *= v;
        d_o_y[id] *= v;
        d_o_z[id] *= v;
    }
}

__global__ void MulMulC_I_const_kernel(float* d_o_x, float* d_o_y, float* d_o_z,
                                       const float* d_i, uint n){
    uint blockId = get_blockID();
    uint  id = get_threadID(blockId);
    
    if (id < n){
        float c = c_delta;
        float v = d_i[id] * c;
        d_o_x[id] *= v;
        d_o_y[id] *= v;
        d_o_z[id] *= v;
    }
}

/** @brief d_o = d_o * d_i * c  (d_o.x = d_o.x * d_i * c)*/
void
MulMulC_I(float *d_ox, 
	  float *d_oy, 
	  float *d_oz, 
	  const float *d_i, 
	  const float& c, 
	  size_t n,
	  StreamT stream, bool onDev)
{
    dim3 threads(256);
    dim3 grids=make_grid(iDivUp(n, threads.x));
    if (onDev) {
        cudaMemcpyToSymbolAsync(c_delta,&c, sizeof(float),
                                0,cudaMemcpyDeviceToDevice,stream);
        MulMulC_I_const_kernel<<<grids, threads, 0, stream>>>
	    (d_ox, d_oy, d_oz,
	     d_i, n);
    } else {
        MulMulC_I_kernel<<<grids, threads, 0, stream>>>
	    (d_ox, d_oy, d_oz,
	     d_i, c, n);
    }
}

__global__ void Add_MulMulC_kernel(float* d_o_x, float* d_o_y, float* d_o_z,
                                   const float* d_i_x, const float* d_i_y, const float* d_i_z,
                                   const float* d_i1_x, const float* d_i1_y, const float* d_i1_z,
                                   const float* d_i2, float c, uint n){
    uint blockId = get_blockID();
    uint  id     = get_threadID(blockId);
    if (id < n){
        float v = d_i2[id] * c;
        d_o_x[id] = d_i_x[id]  + d_i1_x[id] * v;
        d_o_y[id] = d_i_y[id]  + d_i1_y[id] * v;
        d_o_z[id] = d_i_z[id]  + d_i1_z[id] * v;
    }
}

__global__ void Add_MulMulC_const_kernel(float* d_o_x, float* d_o_y, float* d_o_z,
                                         const float* d_i_x, const float* d_i_y, const float* d_i_z,
                                         const float* d_i1_x, const float* d_i1_y, const float* d_i1_z,
                                         const float* d_i2, uint n){
    uint blockId = get_blockID();
    uint  id     = get_threadID(blockId);
    if (id < n){
        float c = c_delta;
        float v = d_i2[id] * c;
        d_o_x[id] = d_i_x[id]  + d_i1_x[id] * v;
        d_o_y[id] = d_i_y[id]  + d_i1_y[id] * v;
        d_o_z[id] = d_i_z[id]  + d_i1_z[id] * v;
    }
}

/** @brief d_o = d_i + d_i1 * d_i2 * c */
void
Add_MulMulC(float *d_ox, 
	    float *d_oy, 
	    float *d_oz, 
	    const float *d_ix, 
	    const float *d_iy, 
	    const float *d_iz, 
	    const float *d_i1x, 
	    const float *d_i1y, 
	    const float *d_i1z, 
	    const float *d_i2, 
	    const float& c,  
	    size_t n,
	    StreamT stream, bool onDev)
{
    dim3 threads(256);
    dim3 grids=make_grid(iDivUp(n, threads.x));
    
    if (onDev) {
        cudaMemcpyToSymbolAsync(c_delta,&c,sizeof(float),
                                0,cudaMemcpyDeviceToDevice,stream);
        Add_MulMulC_const_kernel<<<grids, threads, 0, stream>>>
            (d_ox, d_oy, d_oz, d_ix, d_iy, d_iz, d_i1x, d_i1y, d_i1z,d_i2, n);
    } else {
        Add_MulMulC_kernel<<<grids, threads, 0, stream>>>
            (d_ox, d_oy, d_oz, d_ix, d_iy, d_iz, d_i1x, d_i1y, d_i1z,d_i2, c, n);
    }
}

/** @brief d_o = d_o + d_i * d_i1 * c */
__global__ void Add_MulMulC_I_kernel(float* d_o_x, float* d_o_y, float* d_o_z,
                                     const float* d_i_x, const float* d_i_y, const float* d_i_z,
                                     const float* d_i1, float c, uint n){
    uint blockId = get_blockID();
    uint  id     = get_threadID(blockId);
    if (id < n){
        float v = d_i1[id] * c;
        d_o_x[id] += d_i_x[id] * v;
        d_o_y[id] += d_i_y[id] * v;
        d_o_z[id] += d_i_z[id] * v;
    }
}

__global__ void Add_MulMulC_I_const_kernel(float* d_o_x, float* d_o_y, float* d_o_z,
                                     const float* d_i_x, const float* d_i_y, const float* d_i_z,
                                     const float* d_i1, uint n){
    uint blockId = get_blockID();
    uint  id     = get_threadID(blockId);
    if (id < n){
        float c = c_delta;
        float v = d_i1[id] * c;
        d_o_x[id] += d_i_x[id] * v;
        d_o_y[id] += d_i_y[id] * v;
        d_o_z[id] += d_i_z[id] * v;
    }
}

void 
Add_MulMulC_I(float *d_ox, 
	      float *d_oy, 
	      float *d_oz, 
	      const float *d_ix, 
	      const float *d_iy, 
	      const float *d_iz, 
	      const float *d_i1, 
	      const float& c, 
	      size_t n,
	      StreamT stream, 
	      bool onDev)
{
    dim3 threads(256);
    dim3 grids=make_grid(iDivUp(n, threads.x));
    if (onDev ) {
        cudaMemcpyToSymbolAsync(c_delta,&c,sizeof(float),
                                0,cudaMemcpyDeviceToDevice,stream);
        Add_MulMulC_I_const_kernel<<<grids, threads, 0, stream>>>
            (d_ox, d_oy, d_oz, d_ix, d_iy, d_iz, d_i1, n);
    } else {
        Add_MulMulC_I_kernel<<<grids, threads, 0, stream>>>
            (d_ox, d_oy, d_oz, d_ix, d_iy, d_iz, d_i1, c, n);
    }
}

__global__ void JacDetH2D_kernel(
    float* d_detJ,
    const float* d_Xgx, const float* d_Xgy,
    const float* d_Ygx, const float* d_Ygy,
    uint n)
{
   uint blockId = blockIdx.y * gridDim.x + blockIdx.x;
   uint id      = blockId * blockDim.x + threadIdx.x;
   
   if (id < n){
      d_detJ[id] = d_Xgx[id]*d_Ygy[id] - d_Xgy[id]*d_Ygx[id];
   }
}
                                       
__global__ void JacDetH3D_kernel(
    float* d_detJ,
    const float* d_Xgx, const float* d_Xgy, const float* d_Xgz,
    const float* d_Ygx, const float* d_Ygy, const float* d_Ygz,
    const float* d_Zgx, const float* d_Zgy, const float* d_Zgz, uint n)
{
    uint blockId = blockIdx.y * gridDim.x + blockIdx.x;
    uint id      = blockId * blockDim.x + threadIdx.x;
    
    if (id < n){
        float a00 = d_Xgx[id], a01 = d_Xgy[id], a02 = d_Xgz[id];
        float a10 = d_Ygx[id], a11 = d_Ygy[id], a12 = d_Ygz[id];
        float a20 = d_Zgx[id], a21 = d_Zgy[id], a22 = d_Zgz[id];
        
        d_detJ[id] = det(a00, a01, a02,
                         a10, a11, a12,
                         a20, a21, a22);
    }
}
                                       
void
JacDetH(float *d_detJ,
        const float *d_Xgx, 
        const float *d_Xgy, 
        const float *d_Xgz, 
	const float *d_Ygx, 
	const float *d_Ygy, 
	const float *d_Ygz, 
	const float *d_Zgx, 
	const float *d_Zgy, 
	const float *d_Zgz, 
	size_t n,
	bool slice,
	StreamT stream)
{
    dim3 threads(256);
    dim3 grids=make_grid(iDivUp(n, threads.x));

    if(slice){
       JacDetH2D_kernel<<<grids, threads, 0, stream>>>
	  (d_detJ,
	   d_Xgx, d_Xgy,
	   d_Ygx, d_Ygy,
	   n);
    }else{
       JacDetH3D_kernel<<<grids, threads, 0, stream>>>
	  (d_detJ,
	   d_Xgx, d_Xgy, d_Xgz,
	   d_Ygx, d_Ygy, d_Ygz,
	   d_Zgx, d_Zgy, d_Zgz, n);
    }
}

template<bool fwd>
__global__ void JacDetV_kernel(
    float* d_detJ,
    const float* d_Xgx, const float* d_Xgy, const float* d_Xgz,
    const float* d_Ygx, const float* d_Ygy, const float* d_Ygz,
    const float* d_Zgx, const float* d_Zgy, const float* d_Zgz, uint n)
{
    uint blockId = blockIdx.y * gridDim.x + blockIdx.x;
    uint id      = blockId * blockDim.x + threadIdx.x;

    if (id < n){
        float a00, a01, a02, a10, a11, a12, a20, a21, a22;

        if (fwd) {
            a00 = 1.f + d_Xgx[id], a01 = d_Xgy[id], a02 = d_Xgz[id];
            a10 = d_Ygx[id], a11 = 1.f + d_Ygy[id], a12 = d_Ygz[id];
            a20 = d_Zgx[id], a21 = d_Zgy[id], a22 = 1.f + d_Zgz[id];
        }else {
            a00 = 1.f - d_Xgx[id], a01 = d_Xgy[id], a02 = d_Xgz[id];
            a10 = d_Ygx[id], a11 = 1.f - d_Ygy[id], a12 = d_Ygz[id];
            a20 = d_Zgx[id], a21 = d_Zgy[id], a22 = 1.f - d_Zgz[id];
        }
        
        d_detJ[id] = det(a00, a01, a02,
                         a10, a11, a12,
                         a20, a21, a22);
    }
}

template<bool fwd>
void
JacDetV(float *d_detJ,
        const float *d_Xgx, 
        const float *d_Xgy, 
        const float *d_Xgz, 
	const float *d_Ygx, 
	const float *d_Ygy, 
	const float *d_Ygz, 
	const float *d_Zgx, 
	const float *d_Zgy, 
	const float *d_Zgz, 
	size_t n,
	StreamT stream)
{
    dim3 threads(256);
    dim3 grids=make_grid(iDivUp(n, threads.x));
    JacDetV_kernel<fwd><<<grids, threads, 0, stream>>>
	(d_detJ,
	 d_Xgx, d_Xgy, d_Xgz,
	 d_Ygx, d_Ygy, d_Ygz,
	 d_Zgx, d_Zgy, d_Zgz, n);
}


template<DiffT diffType, BoundaryCondT bc, bool slice>
__global__ void JacDetHPointwise_kernel
(float* d_jdet, 
 const float* d_Hx, const float* d_Hy, const float *d_Hz,
 int szX, int szY, int szZ,
 float ispX, float ispY, float ispZ)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;

    if (i < szX && j < szY)
    {
        int id = j * szX + i;

        for (int k=0; k< szZ; ++k, id+= szX*szY)
        {
	   // from FiniteDiff.h
	   jacDetPoint<float,diffType,bc,slice>
	      (d_jdet[id], 
	       d_Hx,d_Hy,d_Hz,
	       i,j,k,
	       szX,szY,szZ,
	       ispX,ispY,ispZ);
	   
        }
    }
}

template<DiffT diffType, BoundaryCondT bc, bool slice>
void
g_JacDetHPointwise(float *d_jdet,
		   const float *d_hx, 
		   const float *d_hy, 
		   const float *d_hz, 
		   Vec3Di sz,
		   Vec3Df sp,
		   StreamT stream)
{
    Vec3Df isp(1.0/sp.x, 1.0/sp.y, 1.0/sp.z);

    dim3 threads(16,16);
    dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));

    JacDetHPointwise_kernel<diffType, bc, slice>
       <<<grids, threads, 0, stream>>>
	(d_jdet,
	 d_hx, d_hy, d_hz,
	 sz.x, sz.y, sz.z,
	 isp.x, isp.y, isp.z);

}

template<DiffT diffType>
void
g_JacDetHPointwise(float *d_jdet,
		   const float *d_hx, 
		   const float *d_hy, 
		   const float *d_hz, 
		   Vec3Di sz,
		   Vec3Df sp,
		   BoundaryCondT bc,
		   StreamT stream)
{
   bool slice = (sz.z == 1);
   if(bc == BC_APPROX){
      if(slice){
	  PyCA::g_JacDetHPointwise<diffType, BC_APPROX, SLICE_TRUE>
	      (d_jdet, 
	       d_hx, d_hy, d_hz, 
	       sz, sp, stream);
      }else{
	 PyCA::g_JacDetHPointwise<diffType, BC_APPROX, SLICE_FALSE>
	     (d_jdet, 
	      d_hx, d_hy, d_hz, 
	      sz, sp, stream);

      }
   }else if(bc == BC_WRAP){
      if(slice){
	 PyCA::g_JacDetHPointwise<diffType, BC_WRAP, SLICE_TRUE>
	     (d_jdet, 
	      d_hx, d_hy, d_hz, 
	      sz, sp, stream);
      }else{
	 PyCA::g_JacDetHPointwise<diffType, BC_WRAP, SLICE_FALSE>
	     (d_jdet, 
	      d_hx, d_hy, d_hz, 
	      sz, sp, stream);
      }
   }else if(bc == BC_CLAMP){
      if(slice){
	 PyCA::g_JacDetHPointwise<diffType, BC_CLAMP, SLICE_TRUE>
	     (d_jdet, 
	      d_hx, d_hy, d_hz, 
	      sz, sp, stream);
      }else{
	 PyCA::g_JacDetHPointwise<diffType, BC_CLAMP, SLICE_FALSE>
	     (d_jdet, 
	      d_hx, d_hy, d_hz, 
	      sz, sp, stream);
      }
   }else{
      throw PyCAException(__FILE__, __LINE__, "Unknownd BoundaryCondT");
   }

}

void
JacDetH(float *d_jdet, 
	  const float *d_hx, 
	  const float *d_hy, 
	  const float *d_hz, 
	  Vec3Di sz,
	  Vec3Df sp,
	  DiffT diffType, 
	  BoundaryCondT bc, 
	  StreamT stream)
{

   if(diffType == DIFF_FORWARD){
      PyCA::g_JacDetHPointwise<DIFF_FORWARD>
	  (d_jdet, d_hx, d_hy, d_hz, sz, sp, bc, stream);
   }else if(diffType == DIFF_BACKWARD){
      PyCA::g_JacDetHPointwise<DIFF_BACKWARD>
	  (d_jdet, d_hx, d_hy, d_hz, sz, sp, bc, stream);
   }else if(diffType == DIFF_CENTRAL){
      PyCA::g_JacDetHPointwise<DIFF_CENTRAL>
	  (d_jdet, d_hx, d_hy, d_hz, sz, sp, bc, stream);
   }else{
       throw PyCAException(__FILE__, __LINE__, "unknown DiffT");
   }
   
}

// template instantiation
#include "GImageFieldOperKernels_inst.cxx"

    
} // end namespace PyCA

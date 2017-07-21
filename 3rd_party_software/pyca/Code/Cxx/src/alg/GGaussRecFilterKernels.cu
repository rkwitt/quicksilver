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

#include "GGaussRecFilterKernels.h"
#include "pycaUtils.h"

// TEST make sure boost isn't included in nvcc code
#if defined(BOOST_COMPILER)
int bla[-1];
#endif

namespace PyCA {

template<class T, bool clampToEdge>
__global__ void RGFFilter2D_kernel(T* d_o, const T* d_i,
                                   int sizeX, int sizeY, 
                                   float a0, float a1, float a2, float a3, float b1, float b2, float coefp, float coefn)
{
    uint x = blockIdx.x*blockDim.x + threadIdx.x;
    if ((x >= sizeX))
        return;

    d_o += x;
    d_i += x;
    
    T xp = (T)0; // previous input
    T yp = (T)0; // previous output
    T yb = (T)0; // previous output by 2
    
    if (clampToEdge){
        xp = *d_i; yb = coefp*xp; yp = yb;
    }

    for (int y = 0; y < sizeY; y++) {
        float xc = *d_i;
        float yc = a0*xc + a1*xp - b1*yp - b2*yb;
		*d_o = yc;

        //shifting around input output 
        xp = xc; yb = yp; yp = yc;

        // move to the next row
        d_i += sizeX; d_o += sizeX;    // move to next rosizeX
    }

    // reset pointers to point to last element in column
    d_i -= sizeX;
    d_o -= sizeX;

    // reverse pass
    // ensures response is symmetrical
    float xn =0.0f, xa = 0.0f, yn = 0.0f, ya = 0.0f;
    
    if (clampToEdge){
        xn = xa = *d_i; yn = coefn*xn; ya = yn;
    }
    
    for (int y = sizeY-1; y >= 0; y--) {
        float xc = *d_i;
        float yc = a2*xn + a3*xa - b1*yn - b2*ya;
        
		*d_o = *d_o +  yc;

        //shifting around input output 
        xa = xn; xn = xc; ya = yn; yn = yc;
        
        d_o -= sizeX; d_i -= sizeX;  // move to previous row
    }
}

template<class T, bool clampToEdge>
__global__ void RGFFilter3D_kernel(T* d_o, const T* d_i,
                                   int sizeX, int sizeY, int sizeZ,
                                   float a0, float a1, float a2, float a3, float b1, float b2, float coefp, float coefn)
{
    uint x = blockIdx.x*blockDim.x + threadIdx.x;
    uint y = blockIdx.y*blockDim.y + threadIdx.y;

    if ((x >= sizeX) || (y >= sizeY))
        return;

    uint id = x + y * sizeX;
    const uint planeSize = sizeX * sizeY;

    d_o += id;
    d_i += id;
    
    T xp = (T)0; // previous input
    T yp = (T)0; // previous output
    T yb = (T)0; // previous output by 2
    
    if (clampToEdge){
        xp = *d_i; yb = coefp*xp; yp = yb;
    }

    for (int z = 0; z < sizeZ; z++) {
        T xc = *d_i;
        T yc = a0*xc + a1*xp - b1*yp - b2*yb;
		*d_o = yc;

        //shifting around input output 
        xp = xc; yb = yp; yp = yc;

        // move to next plane
        d_i += planeSize;
        d_o += planeSize;    
    }

    // reset pointers to point to last element in column
    d_i -= planeSize;
    d_o -= planeSize;

    // reverse pass
    // ensures response is symmetrical
    T xn = (T)(0.0f);
    T xa = (T)(0.0f);
    T yn = (T)(0.0f);
    T ya = (T)(0.0f);
    
    if (clampToEdge){
        xn = xa = *d_i; yn = coefn*xn; ya = yn;
    }
        
    for (int z = sizeZ-1; z >= 0; z--) {
        T xc = *d_i;
        T yc = a2*xn + a3*xa - b1*yn - b2*ya;
        *d_o = *d_o + yc;

        //shifting around input output 
        xa = xn;
        xn = xc;
        ya = yn;
        yn = yc;

        // move to previous plane
        d_i -= planeSize;
        d_o -= planeSize;  
    }
}

void
ConvolutionX3D(float* d_o, const float* d_i, 
	       const GaussRecFilterBase<EXEC_GPU>::GaussRecParams& p,
               size_t sizeX, size_t sizeY, size_t sizeZ, 
	       StreamT stream)
{
    dim3 threads(16,16);
    dim3 grids(iDivUp(sizeX, threads.x),iDivUp(sizeY, threads.y));

    RGFFilter3D_kernel<float, true><<<grids,threads, 0, stream>>>
        (d_o, d_i,
         sizeX, sizeY, sizeZ,
         p.a0, p.a1, p.a2, p.a3,
         p.b1, p.b2, p.coefp, p.coefn);
}

} // end namespace PyCA

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

#include "GFluidKernelFFTKernels.h"
#include <pycaUtils.h>

// TEST make sure boost isn't included in nvcc code
#if defined(BOOST_COMPILER)
int bla[-1];
#endif

namespace PyCA {

// ================================================================
// Begin Device Code
// ================================================================

__device__ __constant__ float sdCosWX[MAX_FFT_TABLE_LENGTH];
__device__ __constant__ float sdSinWX[MAX_FFT_TABLE_LENGTH];
__device__ __constant__ float sdCosWY[MAX_FFT_TABLE_LENGTH];
__device__ __constant__ float sdSinWY[MAX_FFT_TABLE_LENGTH];
__device__ __constant__ float sdCosWZ[MAX_FFT_TABLE_LENGTH];
__device__ __constant__ float sdSinWZ[MAX_FFT_TABLE_LENGTH];

template<class T>
__device__ 
void 
InverseOperatorMultiply(ComplexT<T>& bX, 
			ComplexT<T>& bY, 
			ComplexT<T>& bZ,
			T L00,
			T L10, T L11,
			T L20, T L21, T L22)
{
    T G00;
    T G10, G11;
    T G20, G21, G22;
    T y0, y1, y2;
    //
    // Given that A is pos-def symetric matrix, solve Ax=b by finding
    // cholesky decomposition GG'=A
    // and then performing 2 back-solves, Gy=b and then G'x=y to get x.
    // 
	   
    // 1. find cholesky decomposition by finding G such that GG'=A.
    //    A must be positive definite symetric (we assume that here)
    //    G is then lower triangular, see algorithm 4.2.1 p142-3
    //    in Golub and VanLoan
    // Note: these are in matlab notation 1:3
    // [ G(1,1)   0      0    ]   [ G(1,1) G(2,1) G(3,1) ]   
    // [ G(2,1) G(2,2)   0    ] * [   0    G(2,2) G(3,2) ] = Amatrix
    // [ G(3,1) G(3,2) G(3,3) ]   [   0      0    G(3,3) ]

    T bRealX = bX.x;
    T bRealY = bY.x;
    T bRealZ = bZ.x;
  
    T bImagX = bX.y;
    T bImagY = bY.y;
    T bImagZ = bZ.y;

    T& vRealX = bX.x;
    T& vRealY = bY.x;
    T& vRealZ = bZ.x;
  
    T& vImagX = bX.y;
    T& vImagY = bY.y;
    T& vImagZ = bZ.y;

    G00 = sqrt(L00);
    G10 = L10 / G00;
    G20 = L20 / G00;

    G11 = L11 - G10 * G10;
    G21 = L21 - G20 * G10;
    G11 = sqrt(G11);
    G21 = G21 / G11;

    G22 = L22 - (G20*G20 + G21*G21);
    G22 = sqrt(G22);

    // back-solve Gy=b to get a temporary vector y
    // back-solve G'x=y to get answer in x
    //
    // Note: these are in matlab notation 1:3
    // [ G(1,1)   0      0    ]   [ y(1) ] = b(1)
    // [ G(2,1) G(2,2)   0    ] * [ y(2) ] = b(2)
    // [ G(3,1) G(3,2) G(3,3) ]   [ y(3) ] = b(3)
    //
    // [ G(1,1) G(2,1) G(3,1) ]   [ x(1) ] = y(1)
    // [   0    G(2,2) G(3,2) ] * [ x(2) ] = y(2)
    // [   0      0    G(3,3) ]   [ x(3) ] = y(3)
    y0 = bRealX / G00;
    y1 = (bRealY - G10*y0) / G11;
    y2 = (bRealZ - G20*y0 - G21*y1) / G22;

    vRealZ = y2 / G22;
    vRealY = (y1 - G21*vRealZ) / G11;
    vRealX = (y0 - G10*vRealY - G20*vRealZ) / G00;

    y0 = bImagX / G00;
    y1 = (bImagY - G10*y0) / G11;
    y2 = (bImagZ - G20*y0 - G21*y1) / G22;

    vImagZ = y2 / G22;
    vImagY = (y1 - G21*vImagZ) / G11;
    vImagX = (y0 - G10*vImagY - G20*vImagZ) / G00;
}

//--------------------------------------------------------------------------------
// General Navier Stoker solver  with beta is different than 0
//
//--------------------------------------------------------------------------------
template<class T>
__device__ 
void 
OperatorMultiply(ComplexT<T>& bX, 
		 ComplexT<T>& bY, 
		 ComplexT<T>& bZ,
		 T L00,
		 T L10, T L11,
		 T L20, T L21, T L22)
{

    T bRealX = bX.x;
    T bRealY = bY.x;
    T bRealZ = bZ.x;
  
    T bImagX = bX.y;
    T bImagY = bY.y;
    T bImagZ = bZ.y;

    T& vRealX = bX.x;
    T& vRealY = bY.x;
    T& vRealZ = bZ.x;
  
    T& vImagX = bX.y;
    T& vImagY = bY.y;
    T& vImagZ = bZ.y;

    vRealX = L00*bRealX + L10*bRealY + L20*bRealZ;
    vRealY = L10*bRealX + L11*bRealY + L21*bRealZ;
    vRealZ = L20*bRealX + L21*bRealY + L22*bRealZ;
  
    vImagX = L00*bImagX + L10*bImagY + L20*bImagZ;
    vImagY = L10*bImagX + L11*bImagY + L21*bImagZ;
    vImagZ = L20*bImagX + L21*bImagY + L22*bImagZ;

}

template<class T>
__device__ inline
void 
ProjectIncomp(ComplexT<T>& bX, ComplexT<T>& bY, ComplexT<T>& bZ,
            T sWXx, T sWYy, T sWZz)
{
    // in Fourier space we project onto (-i*sin(u),-i*sin(v),-i*sin(w)) and remove that component
    // 2008 jdh
  
    T bRealX = bX.x;
    T bRealY = bY.x;
    T bRealZ = bZ.x;
  
    T bImagX = bX.y;
    T bImagY = bY.y;
    T bImagZ = bZ.y;

    T& vRealX = bX.x;
    T& vRealY = bY.x;
    T& vRealZ = bZ.x;
  
    T& vImagX = bX.y;
    T& vImagY = bY.y;
    T& vImagZ = bZ.y;
  
    T nsq = sWXx*sWXx + sWYy*sWYy + sWZz*sWZz; // norm squared of projection vector
  
    // S=(sinwx,sinwy,sinwz)
    // Real part of S dot V in Fourier
    T ReSdotV = ( bRealX*sWXx
    		       +bRealY*sWYy
    		       +bRealZ*sWZz);
    // Imag part of S dot V in Fourier
    T ImSdotV = ( bImagX*sWXx
    		       +bImagY*sWYy
    		       +bImagZ*sWZz);
  
    // Subtract S dot V (normalizing S)
    vRealX = bRealX - ReSdotV*sWXx/nsq;
    vRealY = bRealY - ReSdotV*sWYy/nsq;
    vRealZ = bRealZ - ReSdotV*sWZz/nsq;
  
    vImagX = bImagX - ImSdotV*sWXx/nsq;
    vImagY = bImagY - ImSdotV*sWYy/nsq;
    vImagZ = bImagZ - ImSdotV*sWZz/nsq;
}

template<bool inverseOp, bool incomp, class T>
__global__ 
void 
fullNavierStokesSolver3D_C3_kernel(ComplexT<T>* bX, 
				   ComplexT<T>* bY, 
				   ComplexT<T>* bZ,
				   const T alpha, const T beta, const T gamma,
				   const int sizeX, 
				   const int sizeY, 
				   const int sizeZ)
{
    uint  x = blockIdx.x * blockDim.x + threadIdx.x;
    uint  y = blockIdx.y * blockDim.y + threadIdx.y;

    T lambda;
    T L00;
    T L10, L11;
    T L20, L21, L22;

    uint index     = x + y * sizeX;
    uint planeSize = sizeX * sizeY;
    
    if ( x < sizeX && y < sizeY){
        T wx = sdCosWX[x];
        T wy = sdCosWY[y];
        for (int z=0; z < sizeZ ; ++z, index+=planeSize){
            //
            // compute L (it is symmetric, only need lower triangular part)
            //
            
            lambda = -alpha * (wx + wy + sdCosWZ[z]) + gamma;
            
            L00 = lambda - beta * sdCosWX[x];
            L11 = lambda - beta * sdCosWY[y];
            L22 = lambda - beta * sdCosWZ[z];
            
            L10 = beta * sdSinWX[x] * sdSinWY[y];
            L20 = beta * sdSinWX[x] * sdSinWZ[z];
            L21 = beta * sdSinWY[y] * sdSinWZ[z];
            
            if(inverseOp){
                InverseOperatorMultiply<T>(bX[index], bY[index], bZ[index],
					   L00, L10, L11, L20, L21, L22);
            }else{
                OperatorMultiply<T>(bX[index], bY[index], bZ[index],
				    L00, L10, L11, L20, L21, L22);
            }

            if(incomp){
                if (index > 0) // do nothing with zero component
                    ProjectIncomp<T>(bX[index], bY[index], bZ[index],
                        sdSinWX[x], sdSinWY[y], sdSinWZ[z]);
            }
        }
    }
}

// ================================================================
// End Device Code
// ================================================================

template<class T>
void 
UploadFFTLookupTable3D
(const T* cosWX, const T* cosWY, const T* cosWZ,
 const T* sinWX, const T* sinWY, const T* sinWZ,
 size_t sz_x, size_t sz_y, size_t sz_z)
{
    cudaMemcpyToSymbol(sdCosWX, cosWX, sz_x * sizeof(T));
    cudaMemcpyToSymbol(sdSinWX, sinWX, sz_x * sizeof(T));
    cudaMemcpyToSymbol(sdCosWY, cosWY, sz_y * sizeof(T));
    cudaMemcpyToSymbol(sdSinWY, sinWY, sz_y * sizeof(T));
    cudaMemcpyToSymbol(sdCosWZ, cosWZ, sz_z * sizeof(T));
    cudaMemcpyToSymbol(sdSinWZ, sinWZ, sz_z * sizeof(T));
}

template<class T>
void 
frequencyDomainApply(ComplexT<T>* fftX,
		     ComplexT<T>* fftY,
		     ComplexT<T>* fftZ,
		     float alpha,
		     float beta,
		     float gamma,
		     int cSz_x, int cSz_y, int cSz_z,
		     bool inverseOp, bool divergenceFree,
		     StreamT stream)
{
    dim3 threads(16,16);
    dim3 grids(iDivUp(cSz_x, threads.x),
	       iDivUp(cSz_y, threads.y));


    if (inverseOp == true) {
        if (divergenceFree == true) {
            fullNavierStokesSolver3D_C3_kernel<true, true, T>
            <<<grids, threads, 0, stream>>>
            (fftX, fftY, fftZ,
	     alpha, beta, gamma,
             cSz_x, cSz_y, cSz_z);
         }
        else {
            fullNavierStokesSolver3D_C3_kernel<true, false, T>
            <<<grids, threads, 0, stream>>>
            (fftX, fftY, fftZ,
	     alpha, beta, gamma,
             cSz_x, cSz_y, cSz_z);
         }
    } else {
        if (divergenceFree == true) {
            fullNavierStokesSolver3D_C3_kernel<false, true, T>
            <<<grids, threads, 0, stream>>>
            (fftX, fftY, fftZ,
	     alpha, beta, gamma,
             cSz_x, cSz_y, cSz_z);
         }
        else {
            fullNavierStokesSolver3D_C3_kernel<false, false, T>
            <<<grids, threads, 0, stream>>>
            (fftX, fftY, fftZ,
	     alpha, beta, gamma,
             cSz_x, cSz_y, cSz_z);
         }
    }

}

// template instantiation
template 
void 
UploadFFTLookupTable3D<float>
(const float*, const float*, const float*, 
 const float*, const float*, const float*, 
 size_t, size_t, size_t);

template
void 
frequencyDomainApply<float>
(ComplexT<float>*,
 ComplexT<float>*,
 ComplexT<float>*,
 float alpha,
 float beta,
 float gamma,
 int cSz_x, int cSz_y, int cSz_z,
 bool inverseOp, bool divergenceFree,
 StreamT stream);

} // end namespace PyCA

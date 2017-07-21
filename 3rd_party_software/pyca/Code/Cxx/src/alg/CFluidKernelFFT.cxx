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

#include <CFluidKernelFFT.h>
#include <mem.h>
#include <MemoryManager.h>
#include <MemOpers.h>
#include <FieldOpers.h>

#include "PyCAThread.h"

namespace PyCA {

/* ================================================================
 * FluidKernelFFTWWrapper Float Specialization
 ================================================================ */

FluidKernelFFTWWrapper<float>::
FluidKernelFFTWWrapper() :
    mFFTWForwardPlan(NULL),
    mFFTWBackwardPlan(NULL)
{
}

void
FluidKernelFFTWWrapper<float>::
Initialize(Vec3Di logicalSize, 
	   ComplexT<float> *array,
	   int arrayStride,
	   int nThreads,
	   bool measure)
{
    //
    // create the plans
    //

    // 3D arrays
    int rank = 3;
    // size of real array
    int logicalSizeParam[3];
    logicalSizeParam[0] = logicalSize.z;
    logicalSizeParam[1] = logicalSize.y;
    logicalSizeParam[2] = logicalSize.x;
    
    // number of arrays (one for each dimension)
    int howMany = 3;
    int stride  = 1; // each dimension contiguous in memory
    int dist    = arrayStride; // offset to array for next dimension
    float *dataPtr = (float*)array; // convert from complex* to float*
    
    fftwf_plan_with_nthreads(nThreads);

    unsigned flags = (measure ? FFTW_MEASURE : FFTW_ESTIMATE);
    
    this->mFFTWForwardPlan = 
	fftwf_plan_many_dft_r2c(rank, logicalSizeParam, howMany, 
				dataPtr, 
				NULL, stride, 2*dist, 
				(fftwf_complex*) dataPtr, 
				NULL, stride, dist,
				flags);
  
    this->mFFTWBackwardPlan = 
	fftwf_plan_many_dft_c2r(rank, logicalSizeParam, howMany, 
				(fftwf_complex*) dataPtr,
				NULL, stride, dist, 
				dataPtr,
				NULL, stride, 2*dist,
				flags);
  
    if (!this->mFFTWForwardPlan)
	{
	    throw PyCAException(__FILE__, __LINE__, "FFTW forward plan failed to initialize");
	}
    if (!this->mFFTWBackwardPlan)
	{
	    throw PyCAException(__FILE__, __LINE__, "FFTW backward plan failed to initialize");
	}
  
}

void
FluidKernelFFTWWrapper<float>::
Delete()
{
    if(this->mFFTWForwardPlan){
    	fftwf_destroy_plan(this->mFFTWForwardPlan);
    	this->mFFTWForwardPlan = NULL;
    }
    if(this->mFFTWBackwardPlan){
    	fftwf_destroy_plan(this->mFFTWBackwardPlan);
    	this->mFFTWBackwardPlan = NULL;
    }
}

void
FluidKernelFFTWWrapper<float>::
ExecuteForward(){
    fftwf_execute(this->mFFTWForwardPlan);
}

void
FluidKernelFFTWWrapper<float>::
ExecuteBackward(){
    fftwf_execute(this->mFFTWBackwardPlan);
}

/* ================================================================
 * FluidKernelFFTWWrapper Double Specialization
 ================================================================ */

FluidKernelFFTWWrapper<double>::
FluidKernelFFTWWrapper() :
    mFFTWForwardPlan(NULL),
    mFFTWBackwardPlan(NULL)
{
}

void
FluidKernelFFTWWrapper<double>::
Initialize(Vec3Di logicalSize, 
	   ComplexT<double> *array, 
	   int arrayStride,
	   int nThreads,
	   bool measure)
{
    // create the plans
    int rank = 3;
    int logicalSizeParam[3];
    logicalSizeParam[0] = logicalSize.z;
    logicalSizeParam[1] = logicalSize.y;
    logicalSizeParam[2] = logicalSize.x;
  
    int howMany = 3;
    int stride  = 1; // each dimension contiguous in memory
    int dist    = arrayStride; // offset to array for next dimension
    double *dataPtr = (double*)array;
  
    unsigned flags = (measure ? FFTW_MEASURE : FFTW_ESTIMATE);

    fftw_plan_with_nthreads(nThreads);
    
    this->mFFTWForwardPlan = 
	fftw_plan_many_dft_r2c(rank, logicalSizeParam, howMany, 
			       dataPtr, 
			       NULL, stride, 2*dist, 
			       (fftw_complex*) dataPtr, 
			       NULL, stride, dist,
			       flags);
  
    this->mFFTWBackwardPlan = 
	fftw_plan_many_dft_c2r(rank, logicalSizeParam, howMany, 
			       (fftw_complex*) dataPtr,
			       NULL, stride, dist, 
			       dataPtr,
			       NULL, stride, 2*dist,
			       flags);
  
    if (!this->mFFTWForwardPlan)
	{
	    throw PyCAException(__FILE__, __LINE__, "FFTW forward plan failed to initialize");
	}
    if (!this->mFFTWBackwardPlan)
	{
	    throw PyCAException(__FILE__, __LINE__, "FFTW backward plan failed to initialize");
	}
  
}

void
FluidKernelFFTWWrapper<double>::
Delete()
{
    if(this->mFFTWForwardPlan)
    	fftw_destroy_plan(this->mFFTWForwardPlan);
    if(this->mFFTWBackwardPlan)
    	fftw_destroy_plan(this->mFFTWBackwardPlan);
}

void
FluidKernelFFTWWrapper<double>::
ExecuteForward(){
    fftw_execute(this->mFFTWForwardPlan);
}

void
FluidKernelFFTWWrapper<double>::
ExecuteBackward(){
    fftw_execute(this->mFFTWBackwardPlan);
}

/* ================================================================
 * CFluidKernelFFT_T
 ================================================================ */

// ======== Public Members ======== //

template<class T> 
Mutex
CFluidKernelFFT_T<T>::mFFTWPlanInitializationMutex;

template<class T>
CFluidKernelFFT_T<T>::
CFluidKernelFFT_T()
    : FluidKernelFFTBase<EXEC_CPU, MEM_HOST, T>()
{
}

template<class T>
CFluidKernelFFT_T<T>::
~CFluidKernelFFT_T()
{
    if (this->mHasFFTPlan)
	DestroyFFTPlan();
}

template<class T>
void
CFluidKernelFFT_T<T>::
CreateFFTPlan()
{

    // need to synchronize fftw plan creation
    mutex_lock(&mFFTWPlanInitializationMutex);
    
    mFFTWWrapper.Initialize(this->mSize,
			    this->mFFTArray.get(),
			    this->mAllocateSize);

    // need to synchronize fftw plan creation
    mutex_unlock(&mFFTWPlanInitializationMutex);
    
    this->mHasFFTPlan = true;
  
}

template<class T>
void 
CFluidKernelFFT_T<T>::
DestroyFFTPlan()
{
    mFFTWWrapper.Delete();
    this->mHasFFTPlan = false;
}

template<class T>
void 
CFluidKernelFFT_T<T>::
toFrequencyDomain(const T* hX, const T* hY, const T* hZ, StreamT s)
{

    // number of elements (of type T) in an X-line of real and complex
    // arrays
    int nRealLineVox = this->mSize.x;
    int nCplxLineVox = 2*this->mComplexSize.x;
    
    //
    // copy real data into complex array with appropriate padding for
    // in-place transform
    //
    for(int line=0;line<this->mSize.y*this->mSize.z;++line){
	cpyArrayH2H(&((T*)this->mFFTArrayX)[line*nCplxLineVox], 
			(const T*)&hX[line*nRealLineVox], nRealLineVox);
	cpyArrayH2H(&((T*)this->mFFTArrayY)[line*nCplxLineVox], 
			(const T*)&hY[line*nRealLineVox], nRealLineVox);
	cpyArrayH2H(&((T*)this->mFFTArrayZ)[line*nCplxLineVox], 
			(const T*)&hZ[line*nRealLineVox], nRealLineVox);
    }

    int nRealVox = this->mSize.prod();
    int nCplxVox = 2*this->mComplexSize.prod();

    // Scale to normalize effects of DFT
    MemOpers<EXEC_CPU, T>::MulC_I((T*)this->mFFTArrayX, 1.0f / nRealVox, nCplxVox, NULL);
    MemOpers<EXEC_CPU, T>::MulC_I((T*)this->mFFTArrayY, 1.0f / nRealVox, nCplxVox, NULL);
    MemOpers<EXEC_CPU, T>::MulC_I((T*)this->mFFTArrayZ, 1.0f / nRealVox, nCplxVox, NULL);

    // execute plan
    mFFTWWrapper.ExecuteForward();

    // // TEST
    // boost::shared_ptr<float> tmp = 
    // // END TEST
    
}

template<class T>
void 
CFluidKernelFFT_T<T>::
toSpaceDomain(T* hX, T* hY, T* hZ, StreamT s)
{

    // execute plan
    mFFTWWrapper.ExecuteBackward();

    // number of elements (of type T) in an X-line of real and complex
    // arrays
    int nRealLineVox = this->mSize.x;
    int nCplxLineVox = 2*this->mComplexSize.x;
    
    //
    // copy real data in padded complex array back to real array
    //
    for(int line=0;line<this->mSize.y*this->mSize.z;++line){
	cpyArrayH2H((T*)&hX[line*nRealLineVox], 
			&((T*)this->mFFTArrayX)[line*nCplxLineVox], 
			nRealLineVox);
	cpyArrayH2H((T*)&hY[line*nRealLineVox], 
			&((T*)this->mFFTArrayY)[line*nCplxLineVox], 
			nRealLineVox);
	cpyArrayH2H((T*)&hZ[line*nRealLineVox], 
			&((T*)this->mFFTArrayZ)[line*nCplxLineVox], 
			nRealLineVox);
    }
}

template<class T>
void 
CFluidKernelFFT_T<T>::
frequencyDomainApply(bool inverseOp, StreamT stream)
{
    if(inverseOp){
	this->applyOperatorOnTheFly<true>();
    }else{
	this->applyOperatorOnTheFly<false>();
    }
}

template<class T>
template<bool inverse>
void 
CFluidKernelFFT_T<T>::
applyOperatorOnTheFly()
{
    // apply operator
    double lambda;
    T L00;
    T L10, L11;
    T L20, L21, L22;

#if defined(REPLACE_ZERO_FREQ)
    ComplexT<T> zeroFreqX = this->mFFTArrayX[0];
    ComplexT<T> zeroFreqY = this->mFFTArrayY[0];
    ComplexT<T> zeroFreqZ = this->mFFTArrayZ[0];
#else
    ComplexT<T> zeroFreqX = {static_cast<T>(0),static_cast<T>(0)};
    ComplexT<T> zeroFreqY = {static_cast<T>(0),static_cast<T>(0)};
    ComplexT<T> zeroFreqZ = {static_cast<T>(0),static_cast<T>(0)};
#endif

    for (int z = 0; z < this->mComplexSize.z; ++z)
	{
	    for (int y = 0; y < this->mComplexSize.y; ++y)
		{
		    for (int x = 0; x < this->mComplexSize.x; ++x)
			{
			    //
			    // compute L (it is symmetric, only need lower triangular part)
			    //
	      
			    // maybe lambda should be stored in a lut
			    // it would reduce computation but may cause cache misses
			    lambda = - this->mAlpha
				* (this->mLookupTable.CosWX()[x] + 
				   this->mLookupTable.CosWY()[y] + 
				   this->mLookupTable.CosWZ()[z]) 
				+ this->mGamma;
	      
			    L00 = lambda - 
				this->mBeta * this->mLookupTable.CosWX()[x];
			    L11 = lambda - 
				this->mBeta * this->mLookupTable.CosWY()[y];
			    L22 = lambda - 
				this->mBeta * this->mLookupTable.CosWZ()[z];
			    L10 = this->mBeta * this->mLookupTable.SinWX()[x] * 
				this->mLookupTable.SinWY()[y];
			    L20 = this->mBeta * this->mLookupTable.SinWX()[x] * 
				this->mLookupTable.SinWZ()[z];
			    L21 = this->mBeta * this->mLookupTable.SinWY()[y] * 
				this->mLookupTable.SinWZ()[z];

			    //
			    // compute F = LV (for real and imaginary parts)
			    //
			    size_t idx = 
				this->mComplexSize.x*
				(z*this->mComplexSize.y + y) 
				+ x;

			    if(inverse){
				this->InverseOperatorMultiply(&this->mFFTArrayX[idx],
							      &this->mFFTArrayY[idx],
							      &this->mFFTArrayZ[idx],
							      L00,
							      L10, L11,
							      L20, L21, L22);
			    }else{
				this->OperatorMultiply(&this->mFFTArrayX[idx],
						       &this->mFFTArrayY[idx],
						       &this->mFFTArrayZ[idx],
						       L00,
						       L10, L11,
						       L20, L21, L22);
			    }
	      
			    if (this->mDivergenceFree && (x | y | z))
			    	{
            // Project onto incompressible field
            this->ProjectIncomp(&this->mFFTArrayX[idx],
                                &this->mFFTArrayY[idx],
                                &this->mFFTArrayZ[idx],
                                x,y,z);
			    	}
			}
		}
	}

    if(inverse && this->mGamma == static_cast<T>(0))
      {
      this->mFFTArrayX[0] = zeroFreqX;
      this->mFFTArrayY[0] = zeroFreqY;
      this->mFFTArrayZ[0] = zeroFreqZ;
      }
    

}


template<class T>
inline
void
CFluidKernelFFT_T<T>::
ProjectIncomp(ComplexT<T>* cplx_x, 
	      ComplexT<T>* cplx_y, 
	      ComplexT<T>* cplx_z,
	      unsigned int x, unsigned int y, unsigned int z)
{
    // in Fourier space we project onto (-i*sin(u),-i*sin(v),-i*sin(w)) and remove that component
    // 2008 jdh
  
    T bRealX = cplx_x->x;
    T bRealY = cplx_y->x;
    T bRealZ = cplx_z->x;
	               
    T bImagX = cplx_x->y;
    T bImagY = cplx_y->y;
    T bImagZ = cplx_z->y;

    T& vRealX = cplx_x->x;
    T& vRealY = cplx_y->x;
    T& vRealZ = cplx_z->x;
	                
    T& vImagX = cplx_x->y;
    T& vImagY = cplx_y->y;
    T& vImagZ = cplx_z->y;
  
    FFTLookupTable3D<T> *lut = &this->mLookupTable;
  
    // This is now in LUT
    double nsq = lut->SinWX()[x]*lut->SinWX()[x]
      + lut->SinWY()[y]*lut->SinWY()[y]
      + lut->SinWZ()[z]*lut->SinWZ()[z]; // norm squared of projection vector
  
    // S=(sinwx,sinwy,sinwz)
    // Real part of S dot V in Fourier
    double ReSdotV = ( bRealX*lut->SinWX()[x]
    		       +bRealY*lut->SinWY()[y]
    		       +bRealZ*lut->SinWZ()[z]);
    // Imag part of S dot V in Fourier
    double ImSdotV = ( bImagX*lut->SinWX()[x]
    		       +bImagY*lut->SinWY()[y]
    		       +bImagZ*lut->SinWZ()[z]);
  
    // // Subtract S dot V (normalizing S)
    // vRealX = bRealX - ReSdotV*lut->SinWX()[x]/lut->nsq(x,y,z);
    // vRealY = bRealY - ReSdotV*lut->SinWY()[y]/lut->nsq(x,y,z);
    // vRealZ = bRealZ - ReSdotV*lut->SinWZ()[z]/lut->nsq(x,y,z);
  
    // vImagX = bImagX - ImSdotV*lut->SinWX()[x]/lut->nsq(x,y,z);
    // vImagY = bImagY - ImSdotV*lut->SinWY()[y]/lut->nsq(x,y,z);
    // vImagZ = bImagZ - ImSdotV*lut->SinWZ()[z]/lut->nsq(x,y,z);
    
    vRealX = bRealX - ReSdotV*lut->SinWX()[x]/nsq;
    vRealY = bRealY - ReSdotV*lut->SinWY()[y]/nsq;
    vRealZ = bRealZ - ReSdotV*lut->SinWZ()[z]/nsq;
  
    vImagX = bImagX - ImSdotV*lut->SinWX()[x]/nsq;
    vImagY = bImagY - ImSdotV*lut->SinWY()[y]/nsq;
    vImagZ = bImagZ - ImSdotV*lut->SinWZ()[z]/nsq;

}

template<class T>
inline
void
CFluidKernelFFT_T<T>::
InverseOperatorMultiply(ComplexT<T>* cplx_x, 
			ComplexT<T>* cplx_y, 
			ComplexT<T>* cplx_z,
                        T& L00,
                        T& L10, T& L11,
                        T& L20, T& L21, T& L22)
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

    T bRealX = cplx_x->x;
    T bRealY = cplx_y->x;
    T bRealZ = cplx_z->x;
	               
    T bImagX = cplx_x->y;
    T bImagY = cplx_y->y;
    T bImagZ = cplx_z->y;

    T& vRealX = cplx_x->x;
    T& vRealY = cplx_y->x;
    T& vRealZ = cplx_z->x;
	                
    T& vImagX = cplx_x->y;
    T& vImagY = cplx_y->y;
    T& vImagZ = cplx_z->y;

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

template<class T>
inline
void
CFluidKernelFFT_T<T>::
OperatorMultiply(ComplexT<T>* cplx_x, 
		 ComplexT<T>* cplx_y, 
		 ComplexT<T>* cplx_z,
		 T& L00,
		 T& L10, T& L11,
		 T& L20, T& L21, T& L22)
{
    T bRealX = cplx_x->x;
    T bRealY = cplx_y->x;
    T bRealZ = cplx_z->x;
	               
    T bImagX = cplx_x->y;
    T bImagY = cplx_y->y;
    T bImagZ = cplx_z->y;

    T& vRealX = cplx_x->x;
    T& vRealY = cplx_y->x;
    T& vRealZ = cplx_z->x;
	                
    T& vImagX = cplx_x->y;
    T& vImagY = cplx_y->y;
    T& vImagZ = cplx_z->y;

    vRealX = L00*bRealX + L10*bRealY + L20*bRealZ;
    vRealY = L10*bRealX + L11*bRealY + L21*bRealZ;
    vRealZ = L20*bRealX + L21*bRealY + L22*bRealZ;
  
    vImagX = L00*bImagX + L10*bImagY + L20*bImagZ;
    vImagY = L10*bImagX + L11*bImagY + L21*bImagZ;
    vImagZ = L20*bImagX + L21*bImagY + L22*bImagZ;	
}

template class CFluidKernelFFT_T<float>;

} // end namespace PyCA

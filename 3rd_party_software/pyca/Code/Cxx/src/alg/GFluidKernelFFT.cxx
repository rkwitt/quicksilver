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

#include <GFluidKernelFFT.h>
#include "GFluidKernelFFTKernels.h"
#include <pycaUtils.h>
#include <CudaUtils.h>
#include <GMemOpers.h>

namespace PyCA {

template<class T>
GFluidKernelFFT_T<T>::
GFluidKernelFFT_T()
    : FluidKernelFFTBase<EXEC_GPU, MEM_DEVICE, T>()
{
}


template<class T>
GFluidKernelFFT_T<T>::
~GFluidKernelFFT_T()
{
  if (this->mHasFFTPlan)
    DestroyFFTPlan();
}

template<class T>
void
GFluidKernelFFT_T<T>::
setSize(const GridInfo& g)
{
    // HAS to check befor the information be changed
    Vec3Di newSize   = g.size();
    Vec3Df newSpacing= g.spacing();

    bool changeSize    = (this->mSize    != newSize);
    bool changeSpacing = (this->mSpacing != newSpacing);

    // superclass version
    FluidKernelFFTBase<EXEC_GPU, MEM_DEVICE, T>::setSize(g);

    // recompute the lookup table
    if (changeSize || changeSpacing){
        this->UploadFFTLookupTable3D();
    }
}

template<class T>
void 
GFluidKernelFFT_T<T>::
CreateFFTPlan()
{
    // TEST
    CudaUtils::CheckCUDAError(__FILE__, __LINE__);
    // END TEST
    // ATTENTION : the order is reversed from CUDA FFT documentation
    // this is hard to find bug since it is not failed if 
    // mComplexSize.z = mComplexSize.x
    cufftPlan3d(&mPlanR2C, this->mSize.z, this->mSize.y, this->mSize.x, CUFFT_R2C);
    cufftPlan3d(&mPlanC2R, this->mSize.z, this->mSize.y, this->mSize.x, CUFFT_C2R);
    CudaUtils::CheckCUDAError(__FILE__, __LINE__);
    this->mHasFFTPlan = true;
}

template<class T>
void 
GFluidKernelFFT_T<T>::
DestroyFFTPlan()
{
    cufftDestroy(mPlanR2C);
    cufftDestroy(mPlanC2R);
    CudaUtils::CheckCUDAError(__FILE__, __LINE__);
    this->mHasFFTPlan = false;
}

template<class T>
void 
GFluidKernelFFT_T<T>::
toFrequencyDomain(const T* dX, const T* dY, const T* dZ, StreamT s){
    // TEST
    CudaUtils::CheckCUDAError(__FILE__, __LINE__);
    // END TEST
    cufftExecR2C(mPlanR2C,(cufftReal*)dX,(cufftComplex*)this->mFFTArrayX);
    cufftExecR2C(mPlanR2C,(cufftReal*)dY,(cufftComplex*)this->mFFTArrayY);
    cufftExecR2C(mPlanR2C,(cufftReal*)dZ,(cufftComplex*)this->mFFTArrayZ);
    CudaUtils::CheckCUDAError(__FILE__, __LINE__);
}

template<class T>
void 
GFluidKernelFFT_T<T>::
toSpaceDomain(T* dX, T* dY, T* dZ, StreamT s) {
    cufftExecC2R(mPlanC2R, (cufftComplex*) this->mFFTArrayX, dX);
    cufftExecC2R(mPlanC2R, (cufftComplex*) this->mFFTArrayY, dY);
    cufftExecC2R(mPlanC2R, (cufftComplex*) this->mFFTArrayZ, dZ);
    CudaUtils::CheckCUDAError(__FILE__, __LINE__);

    int nVox = this->mSize.prod();

    // Scale to normalize effects of DFT
    GMemOpers<T>::MulC_I(dX, 1.0f / nVox, nVox, NULL, false);
    GMemOpers<T>::MulC_I(dY, 1.0f / nVox, nVox, NULL, false);
    GMemOpers<T>::MulC_I(dZ, 1.0f / nVox, nVox, NULL, false);

}

template<class T>
void 
GFluidKernelFFT_T<T>::
UploadFFTLookupTable3D()
{
    PyCA::UploadFFTLookupTable3D(this->mLookupTable.CosWX(),
				 this->mLookupTable.CosWY(),
				 this->mLookupTable.CosWZ(),
				 this->mLookupTable.SinWX(),
				 this->mLookupTable.SinWY(),
				 this->mLookupTable.SinWZ(),
				 this->mSize.x,
				 this->mSize.y,
				 this->mSize.z);
    CudaUtils::CheckCUDAError(__FILE__, __LINE__);
}

template<class T>
void 
GFluidKernelFFT_T<T>::
frequencyDomainApply(bool inverseOp, StreamT stream)
{
    // if gamma is zero, save zero frequency component (zero'th
    // element of FFT arrays) and replace it after kernel
    ComplexT<T> zeroFreqX={0.f,0.f}, 
       zeroFreqY={0.f,0.f}, zeroFreqZ={0.f,0.f};
#if defined(REPLACE_ZERO_FREQ)
    if (inverseOp && this->mGamma == static_cast<T>(0)){
       cudaMemcpy(&zeroFreqX,this->mFFTArrayX,
		  sizeof(ComplexT<T>),cudaMemcpyDeviceToHost);
       cudaMemcpy(&zeroFreqY,this->mFFTArrayY,
		  sizeof(ComplexT<T>),cudaMemcpyDeviceToHost);
       cudaMemcpy(&zeroFreqZ,this->mFFTArrayZ,
		  sizeof(ComplexT<T>),cudaMemcpyDeviceToHost);
    }
#endif

    PyCA::frequencyDomainApply(this->mFFTArrayX,
			       this->mFFTArrayY,
			       this->mFFTArrayZ,
			       this->mAlpha,
			       this->mBeta,
			       this->mGamma,
			       this->mComplexSize.x,
			       this->mComplexSize.y,
			       this->mComplexSize.z,
			       inverseOp,
			       this->mDivergenceFree,
			       stream);

    if (inverseOp && this->mGamma == static_cast<T>(0)){
       cudaMemcpy(this->mFFTArrayX,&zeroFreqX,
		  sizeof(ComplexT<T>),cudaMemcpyHostToDevice);
       cudaMemcpy(this->mFFTArrayY,&zeroFreqY,
		  sizeof(ComplexT<T>),cudaMemcpyHostToDevice);
       cudaMemcpy(this->mFFTArrayZ,&zeroFreqZ,
		  sizeof(ComplexT<T>),cudaMemcpyHostToDevice);
    }
}

template class GFluidKernelFFT_T<float>;

} // end namespace PyCA

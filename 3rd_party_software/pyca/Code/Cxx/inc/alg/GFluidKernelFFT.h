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

#ifndef __G_FLUID_KERNEL_FFT_H__
#define __G_FLUID_KERNEL_FFT_H__

#include <pycaConst.h>
#include <FluidKernelFFTBase.h>

#ifdef CUDA_ENABLED
#include <cufft.h>
#else
// it doesn't matter if this is right or not, GFluidKernelFFT.cu is
// never compiled if cuda is not enabled
typedef int cufftHandle;
#endif

namespace PyCA {

// forward declaration
class Field3D;

template<class T>
class GFluidKernelFFT_T : 
    public FluidKernelFFTBase<EXEC_GPU, MEM_DEVICE, T>
{

public:
  
    GFluidKernelFFT_T();
    ~GFluidKernelFFT_T();

protected:

    virtual void setSize(const GridInfo& g);
  
    virtual void UploadFFTLookupTable3D();

    virtual void CreateFFTPlan();
    virtual void DestroyFFTPlan();

    virtual void toFrequencyDomain(const T* dX, const T* dY, const T* dZ, 
				   StreamT s);
    virtual void toSpaceDomain(T* x, T* y, T* z,
			      StreamT s=NULL);
    virtual void frequencyDomainApply(bool inverseOp, StreamT stream);

    // CUDA FFT plans
    cufftHandle mPlanR2C, mPlanC2R;
};

typedef GFluidKernelFFT_T<float> GFluidKernelFFT;

} // end namespace PyCA

#endif //__G_FLUID_KERNEL_FFT_H__


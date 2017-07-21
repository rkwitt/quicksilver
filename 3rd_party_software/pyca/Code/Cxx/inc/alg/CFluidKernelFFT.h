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

#ifndef __C_FLUID_KERNEL_FFT_H__
#define __C_FLUID_KERNEL_FFT_H__

#include <FluidKernelFFTBase.h>

#include <fftw3.h>

#include <pycaConst.h>

#include "PyCAThread.h"

namespace PyCA {

/**
 * DiffOperFFTWrapper is set up to allow templating the FFTW methods
 * of DiffOper over the precision desired (right now just float or
 * double).  The generic template is never implemented, only the
 * specialized versions for float and double
 */
template<class T>
class FluidKernelFFTWWrapper {

protected:
    FluidKernelFFTWWrapper()
    {
    }
    
public:
    void Initialize(Vec3Di logicalSize, 
		    ComplexT<T> *array, 
		    int arrayStride,
		    int nThreads = 1,
		    bool measure = true){}
    
    void ExecuteForward(){}
    
    void ExecuteBackward(){}
    
    void Delete(){}
    
};

/**
 * single-precision specialization of the wrapper, calls fftwf_
 * version of classes
 */  
template<>
class FluidKernelFFTWWrapper<float> {
public:
    FluidKernelFFTWWrapper();
    void Initialize(Vec3Di logicalSize, 
		    ComplexT<float> *array, 
		    int arrayStride,
		    int nThreads = 1,
		    bool measure = true);
    void Delete();
    void ExecuteForward();
    void ExecuteBackward();
  
protected:
    fftwf_plan mFFTWForwardPlan;
    fftwf_plan mFFTWBackwardPlan;
};

/**
 * double-precision specialization of the wrapper, calls fftwf_
 * version of classes
 */  
template<>
class FluidKernelFFTWWrapper<double> {
public:
    FluidKernelFFTWWrapper();
    void Initialize(Vec3Di logicalSize, 
		    ComplexT<double> *array, 
		    int arrayStride,
		    int nThreads = 1,
		    bool measure = true);
    void Delete();
    void ExecuteForward();
    void ExecuteBackward();
    
protected:
    fftw_plan mFFTWForwardPlan;
    fftw_plan mFFTWBackwardPlan;
};

class Field3D;

template<class T>
class CFluidKernelFFT_T : 
    public FluidKernelFFTBase<EXEC_CPU, MEM_HOST, T>
{
    
public:
    
    CFluidKernelFFT_T();
    ~CFluidKernelFFT_T();
    
protected:

    virtual void CreateFFTPlan();
    virtual void DestroyFFTPlan();

    virtual void toFrequencyDomain(const T* hX, const T* hY, const T* hZ, 
				   StreamT s);
    virtual void toSpaceDomain(T* hX, T* hY, T* hZ,
			      StreamT s);
    virtual void frequencyDomainApply(bool inverseOp, StreamT stream);

    /**
     * Apply forward operator to specific frequency location
     */
    void OperatorMultiply(ComplexT<T>* cplx_x, 
			  ComplexT<T>* cplx_y, 
			  ComplexT<T>* cplx_z,
			  T& L00,
			  T& L10, T& L11,
			  T& L20, T& L21, T& L22);
    
    /**
     * Templated internal implementation of frequencyDomainApply
     */
    template<bool inverse>
    void applyOperatorOnTheFly();

    /**
     * Apply inverse operator to specific frequency location
     */
    void InverseOperatorMultiply(ComplexT<T>* cplx_x, 
				 ComplexT<T>* cplx_y, 
				 ComplexT<T>* cplx_z,
				 T& L00,
				 T& L10, T& L11,
				 T& L20, T& L21, T& L22);

    /**
     * Apply divergence-free projection to specific frequency location
     */
    void ProjectIncomp(ComplexT<T>* cplx_x, 
		       ComplexT<T>* cplx_y, 
		       ComplexT<T>* cplx_z,
		       unsigned int x, unsigned int y, unsigned int z);


    // FFT plans
    FluidKernelFFTWWrapper<T> mFFTWWrapper;

    // mutex for plan initialization
    static Mutex mFFTWPlanInitializationMutex;
};

typedef CFluidKernelFFT_T<float> CFluidKernelFFT;

} // end namespace PyCA

#endif // __C_FLUID_KERNEL_FFT_H__

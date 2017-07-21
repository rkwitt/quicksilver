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

%{
#include "Vec3D.h"
#include "GridInfo.h"
#include "Field3D.h"
#include "Image3D.h"

#include "FluidKernelFFT.h" 
#include "FluidKernelFFTInterface.h" 
#ifdef CUDA_ENABLED
#include "GFluidKernelFFT.h" 
#endif // CUDA_ENABLED
#include "CFluidKernelFFT.h" 

#include "IdentityFilter.h"
#include "GaussianFilter.h"
#include "GaussianFilterInterface.h"

#include "MultiscaleManager.h"
#include "MultiscaleResampler.h"
#include "MultiscaleResamplerInterface.h"

%}

// FluidKernelFFT wrapping

// include doxygen docstrings
%include "FluidKernelFFT.i"
 // wrap file
%include "FluidKernelFFT.h"

 // include doxygen docstrings
%include "GaussianFilter.i"
 // wrap file
%include "GaussianFilter.h"

 // include doxygen docstrings
%include "MultiscaleManager.i"
 // wrap file
%include "MultiscaleManager.h"

 // include doxygen docstrings
%include "MultiscaleResampler.i"
 // wrap file
%include "MultiscaleResampler.h"

%template(GaussianFilterCPU) PyCA::GaussianFilter<EXEC_CPU>;

%template(FluidKernelFFTCPU) PyCA::FluidKernelFFT<EXEC_CPU>;

// identity filtered resampler (doesn't blur, just resample)
%template(MultiscaleResamplerCPU) PyCA::MultiscaleResampler<PyCA::IdentityFilter<EXEC_CPU> >;

//%apply const PyCA::Image3D*& OUTPUT { const PyCA::Image3D*& I_ptr };
%template(MultiscaleResamplerGaussCPU) PyCA::MultiscaleResampler<PyCA::GaussianFilter<EXEC_CPU> >;

#ifdef CUDA_ENABLED
%template(GaussianFilterGPU) PyCA::GaussianFilter<EXEC_GPU>;
%template(FluidKernelFFTGPU) PyCA::FluidKernelFFT<EXEC_GPU>;
%template(MultiscaleResamplerGaussGPU) PyCA::MultiscaleResampler<PyCA::GaussianFilter<EXEC_GPU> >;
%template(MultiscaleResamplerGPU) PyCA::MultiscaleResampler<PyCA::IdentityFilter<EXEC_GPU> >;
#endif // CUDA_ENABLED

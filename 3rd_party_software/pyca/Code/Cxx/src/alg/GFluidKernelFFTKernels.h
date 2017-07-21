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

#ifndef __G_FLUID_KERNEL_FFT_KERNELS_H__
#define __G_FLUID_KERNEL_FFT_KERNELS_H__

#include <estream.h>
#include <mem.h>

namespace PyCA {

template<class T>
void 
UploadFFTLookupTable3D
(const T* cosWX, const T* cosWY, const T* cosWZ,
 const T* sinWX, const T* sinWY, const T* sinWZ,
 size_t sz_x, size_t sz_y, size_t sz_z);

template<class T>
void 
frequencyDomainApply(ComplexT<T>* fftX,
		     ComplexT<T>* fftY,
		     ComplexT<T>* fftZ,
		     float alpha,
		     float beta,
		     float gamma,
		     int cSz_x, int cSz_y, int cSz_z,
		     bool inverseOp, 
		     bool divergenceFree,
		     StreamT stream);

} // end namespace PyCA

#endif // __G_FLUID_KERNEL_FFT_KERNELS_H__

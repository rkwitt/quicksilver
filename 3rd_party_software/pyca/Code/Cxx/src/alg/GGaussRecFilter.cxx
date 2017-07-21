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

#include <GGaussRecFilter.h>
#include "GGaussRecFilterKernels.h"

#include <estream.h>
#include <pycaUtils.h>
// #include <MemOpers.h>
#include "SeparableFilter.h"

namespace PyCA {

GGaussRecFilter::GGaussRecFilter()
   : GaussRecFilterBase<EXEC_GPU>()
{
    
}

void
GGaussRecFilter::
ConvolutionX3D(float* d_o, 
	       const float* d_i, 
	       const GaussRecParams& p,
               size_t sizeX, size_t sizeY, size_t sizeZ, 
	       StreamT stream)
{
    PyCA::ConvolutionX3D(d_o, d_i, p,
			 sizeX, sizeY, sizeZ, 
			 stream);
}

void 
GGaussRecFilter::
filter(float *a_o, const float* a_i, float* a_t, StreamT stream){
    SeparableFilter<GGaussRecFilter, EXEC_GPU>(*this, a_o, a_i, a_t, mSize, stream);
}


} // end namespace PyCA

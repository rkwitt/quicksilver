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

#ifndef __G_GAUSS_REC_FILTER_KERNELS_H__
#define __G_GAUSS_REC_FILTER_KERNELS_H__

#include <GaussRecFilterBase.h>
#include <pycaConst.h>
#include <estream.h>

namespace PyCA {

void ConvolutionX3D(float* d_o, 
		    const float* d_i, 
		    const GaussRecFilterBase<EXEC_GPU>::GaussRecParams& p,
		    size_t sizeX, size_t sizeY, size_t sizeZ, 
		    StreamT stream);

} // end namespace PyCA

#endif // __G_GAUSS_REC_FILTER_KERNELS_H__



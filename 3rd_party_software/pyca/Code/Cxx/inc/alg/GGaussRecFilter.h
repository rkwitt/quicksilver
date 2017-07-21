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

#ifndef __G_GAUSS_REC_FILTER_H__
#define __G_GAUSS_REC_FILTER_H__

// #include <cstdlib>
// #include <cstring>
// #include <climits>

#include <estream.h>
#include <Vec3D.h>
#include <vector>
#include <GaussRecFilterBase.h>
#include <boost/shared_ptr.hpp>
#include <pycaConst.h>

// TEST -- make sure filed including boost aren't leaking into
// nvcc-compiled code
#if defined(PYCA_BOOSTTEST)
#if defined(__CUDACC__)
int bla[-1];
#endif
#endif
// END TEST

namespace PyCA {

class GGaussRecFilter : public GaussRecFilterBase<EXEC_GPU> {
public:
   enum { exec_mode = EXEC_GPU };
   GGaussRecFilter();
   ~GGaussRecFilter(){};
   
   void filter(float *a_o, const float* a_i, float* a_t, StreamT stream=NULL);
protected:
   void ConvolutionX3D(float* d_o, const float* d_i, 
		       const GaussRecParams& p,
		       size_t sizeX, size_t sizeY, size_t sizeZ, StreamT stream);

};

} // end namespace PyCA

#endif // __G_GAUSS_REC_FILTER_H__



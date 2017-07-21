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

#ifndef __G_GAUSSIAN_FILTER_H__
#define __G_GAUSSIAN_FILTER_H__

// #include <cstdlib>
// #include <cstring>
// #include <climits>

#include <estream.h>
#include <Vec3D.h>
#include <vector>
#include <GaussianFilterBase.h>
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

class GGaussianFilter : public GaussianFilterBase<EXEC_GPU> {
public:
    enum { exec_mode = EXEC_GPU };
    GGaussianFilter();
    ~GGaussianFilter(){};
    void filter(float *a_o, const float* a_i, float* a_t, 
		StreamT stream=NULL);
    void filter2(float *a_o, const float* a_i, float* a_t, 
		 StreamT stream=NULL);
    void convolutionSingleAxis(float* a_o, const float *a_i,
			       size_t sizeX, size_t sizeY, size_t sizeZ,
			       int axis, StreamT stream);
private:
    void update();
    boost::shared_ptr<float> mdKx, mdKy, mdKz;
    boost::shared_ptr<float> mdSx, mdSy, mdSz;
};

} // end namespace PyCA

#endif // __G_GAUSSIAN_FILTER_H__



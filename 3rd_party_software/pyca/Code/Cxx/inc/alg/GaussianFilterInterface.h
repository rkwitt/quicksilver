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

#ifndef __GAUSSIAN_FILTER_INTERFACE_H__
#define __GAUSSIAN_FILTER_INTERFACE_H__

#include <Vec3D.h>

namespace PyCA {

class GaussianFilterInterface {

public:
    virtual ~GaussianFilterInterface() = 0;
    virtual void updateParams(const Vec3Di& size, 
		      const Vec3Df& sig, 
		      const Vec3Di kRad) = 0;
    
    virtual void filter(float *a_o, const float* a_i, float* a_t, StreamT stream=NULL) = 0;

    virtual void filter2(float *a_o, const float* a_i, float* a_t, StreamT stream=NULL) = 0;

    virtual void convolutionSingleAxis(float* a_o, const float *a_i,
                                       size_t sizeX, size_t sizeY, size_t sizeZ,
                                       int axis, StreamT stream) = 0;

};

 inline GaussianFilterInterface::~GaussianFilterInterface(){}

} // end namespace PyCA

#endif // __GAUSSIAN_FILTER_INTERFACE_H__

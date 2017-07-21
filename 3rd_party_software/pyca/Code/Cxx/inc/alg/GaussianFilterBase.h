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

#ifndef __GAUSSIAN_FILTER_BASE_H__
#define __GAUSSIAN_FILTER_BASE_H__

#include <vector>

#include <Vec3D.h>
#include <estream.h>

#include "GaussianFilterInterface.h"

namespace PyCA {

template<int mode>
class GaussianFilterBase
    : public GaussianFilterInterface
{
public:
    enum { exec_mode = mode };
    
    GaussianFilterBase();
    virtual ~GaussianFilterBase(){};
    void updateParams(const Vec3Di& size, 
		      const Vec3Df& sig, 
		      const Vec3Di kRad);
    
    virtual void filter(float *a_o, const float* a_i, float* a_t, StreamT stream=NULL);
    virtual void filter2(float *a_o, const float* a_i, float* a_t, StreamT stream=NULL);
    virtual void convolutionSingleAxis(float* a_o, const float *a_i,
                                       size_t sizeX, size_t sizeY, size_t sizeZ,
                                       int axis, StreamT stream);

protected:
    virtual void update();

    std::vector<float> mKx, mKy, mKz;
    std::vector<float> mSx, mSy, mSz;

    Vec3Di mSize;
    Vec3Df mSigma;
    Vec3Di mKRadius;
    Vec3Di mKLength;
};

} // end namespace PyCA

#endif // __GAUSSIAN_FILTER_BASE_H__

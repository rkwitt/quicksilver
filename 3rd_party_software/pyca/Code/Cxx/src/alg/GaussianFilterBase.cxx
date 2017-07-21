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

#include <GaussianFilterBase.h>
#include <conditionMacro.h>
#include <GaussUtils.h>
#include "SeparableFilter.h"
#include <GaussParamUtils.h>

namespace PyCA {

template<int mode>
GaussianFilterBase<mode>::
GaussianFilterBase()
    :mSize(0,0,0), mSigma(0.f,0.f,0.f), mKRadius(0,0,0), mKLength(0,0,0)
{
    
}

template<int mode>
void 
GaussianFilterBase<mode>::
updateParams(const Vec3Di& size, 
	     const Vec3Df& sig, 
	     const Vec3Di kRadius)
{
    mSize = size;

    Vec3Di kRad = kRadius;

    // limit radius to size supported by image
    kRad.x = std::min(kRad.x, std::max(size.x/2-1, 0));
    kRad.y = std::min(kRad.y, std::max(size.y/2-1, 0));
    kRad.z = std::min(kRad.z, std::max(size.z/2-1, 0));

    bool changeSigma  = (mSigma   != sig);
    bool changeRadius = (mKRadius != kRad);

    if (changeSigma){ mSigma = sig; }
    if (changeRadius){ mKRadius = kRad; }

    if (changeSigma || changeRadius) {
        update();
    }
}

template<int mode>
void 
GaussianFilterBase<mode>::
update()
{
    mKLength = 2 * mKRadius + 1;
    
    mKx.resize(mKLength.x);
    mKy.resize(mKLength.y);
    mKz.resize(mKLength.z);

    mSx.resize(mKRadius.x + 1);
    mSy.resize(mKRadius.y + 1);
    mSz.resize(mKRadius.z + 1);
    
    GaussUtils::GenerateKernel(&mKx[0],mSigma.x, mKRadius.x);
    GaussUtils::GenerateKernel(&mKy[0],mSigma.y, mKRadius.y);
    GaussUtils::GenerateKernel(&mKz[0],mSigma.z, mKRadius.z);

    GaussUtils::GenerateSupplementKernel(&mSx[0], &mKx[0], mKRadius.x);
    GaussUtils::GenerateSupplementKernel(&mSy[0], &mKy[0], mKRadius.y);
    GaussUtils::GenerateSupplementKernel(&mSz[0], &mKz[0], mKRadius.z);
}

template<int mode>
void 
GaussianFilterBase<mode>::
filter(float *a_o, const float* a_i, float* a_t, StreamT stream){
    GaussUtils::ConvolutionX3D(a_o,a_i,&mKx[0],&mSx[0],mKRadius.x,mSize.x,mSize.y,mSize.z);
    GaussUtils::ConvolutionY3D(a_t,a_o,&mKy[0],&mSy[0],mKRadius.y,mSize.x,mSize.y,mSize.z);
    GaussUtils::ConvolutionZ3D(a_o,a_t,&mKz[0],&mSz[0],mKRadius.z,mSize.x,mSize.y,mSize.z);
}

template<int mode>
void 
GaussianFilterBase<mode>::
convolutionSingleAxis(float* a_o, const float *a_i,
                      size_t sizeX, size_t sizeY, size_t sizeZ,
                      int axis, StreamT stream)
{
   PRECONDITION(sizeX * sizeY * sizeZ == (size_t)mSize.prod(),
		"Incompatible dimensions");

    if (axis == 0) {
	GaussUtils::ConvolutionX3D(a_o,a_i, &mKx[0],&mSx[0], mKRadius.x,
				       sizeX, sizeY, sizeZ);
    }
    else if (axis == 1) {
	GaussUtils::ConvolutionX3D(a_o,a_i, &mKy[0],&mSy[0], mKRadius.y,
				       sizeX, sizeY, sizeZ);
    }else if (axis == 2) {
	GaussUtils::ConvolutionX3D(a_o,a_i, &mKz[0],&mSz[0], mKRadius.z,
				       sizeX, sizeY, sizeZ);
    }
}

template<int mode>
void 
GaussianFilterBase<mode>::
filter2(float *a_o, 
	const float* a_i, 
	float* a_t, 
	StreamT stream)
{
    SeparableFilter<GaussianFilterBase<mode>, mode>(*this, a_o, a_i, a_t, mSize, stream);
}

// template instantiation
template class GaussianFilterBase<EXEC_CPU>;
#ifdef CUDA_ENABLED
template class GaussianFilterBase<EXEC_GPU>;
#endif

/*--------------------------------------------------------------------------------*/

} // end namespace PyCA

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

#include <GGaussianFilter.h>
#include <boost_mem.h>
#include <GGaussUtils.h>
#include <CudaUtils.h>
#include "SeparableFilter.h"

namespace PyCA {

GGaussianFilter::GGaussianFilter(){
}

void GGaussianFilter::update(){

    GaussianFilterBase<exec_mode>::update();

    mdKx = CreateSharedArray<MEM_DEVICE, float>(mKLength.x);
    mdKy = CreateSharedArray<MEM_DEVICE, float>(mKLength.y);
    mdKz = CreateSharedArray<MEM_DEVICE, float>(mKLength.z);

    mdSx = CreateSharedArray<MEM_DEVICE, float>(mKRadius.x + 1);
    mdSy = CreateSharedArray<MEM_DEVICE, float>(mKRadius.y + 1);
    mdSz = CreateSharedArray<MEM_DEVICE, float>(mKRadius.z + 1);

    GGaussUtils::setConvolutionKernelX(&mKx[0], mKLength.x);
    GGaussUtils::setConvolutionKernelY(&mKy[0], mKLength.y);
    GGaussUtils::setConvolutionKernelZ(&mKz[0], mKLength.z);
    
    GGaussUtils::setSupplementKernelX(&mSx[0], mKRadius.x + 1);
    GGaussUtils::setSupplementKernelY(&mSy[0], mKRadius.y + 1);
    GGaussUtils::setSupplementKernelZ(&mSz[0], mKRadius.z + 1);

    cpyArrayH2D(mdKx.get(), &mKx[0], mKLength.x);
    cpyArrayH2D(mdKy.get(), &mKy[0], mKLength.y);
    cpyArrayH2D(mdKz.get(), &mKz[0], mKLength.z);

    cpyArrayH2D(mdSy.get(), &mSy[0], mKRadius.y + 1);
    cpyArrayH2D(mdSx.get(), &mSx[0], mKRadius.x + 1);
    cpyArrayH2D(mdSz.get(), &mSz[0], mKRadius.z + 1);

    CudaUtils::CheckCUDAError(__FILE__, __LINE__);
}

void GGaussianFilter::
filter(float *a_o, const float* a_i, float* a_t, StreamT stream){
    GGaussUtils::ConvolutionX3D(a_o,a_i,mKRadius.x,mSize.x,mSize.y,mSize.z, stream);
    GGaussUtils::ConvolutionY3D(a_t,a_o,mKRadius.y,mSize.x,mSize.y,mSize.z, stream);
    GGaussUtils::ConvolutionZ3D(a_o,a_t,mKRadius.z,mSize.x,mSize.y,mSize.z, stream);

    CudaUtils::CheckCUDAError(__FILE__, __LINE__);
}

void GGaussianFilter::
filter2(float *a_o, const float* a_i, float* a_t, StreamT stream){
    SeparableFilter<GGaussianFilter, EXEC_GPU>(*this, a_o, a_i, a_t, mSize, stream);
    CudaUtils::CheckCUDAError(__FILE__, __LINE__);
}

void GGaussianFilter::
convolutionSingleAxis(float* a_o, const float *a_i,
                      size_t sizeX, size_t sizeY, size_t sizeZ,
                      int axis, StreamT stream)
{
   PRECONDITION(sizeX * sizeY * sizeZ == (size_t)mSize.prod(), 
		"Incompatible dimensions");
    
    if (axis == 0) {
        GGaussUtils::ConvolutionX3D(a_o,a_i,
					mdKx.get(),mdSx.get(),mKRadius.x,
					sizeX, sizeY, sizeZ, stream);
    }
    else if (axis == 1) {
        GGaussUtils::ConvolutionX3D(a_o,a_i,
					mdKy.get(),mdSy.get(),mKRadius.y,
					sizeX, sizeY, sizeZ, stream);
    }else if (axis == 2) {
        GGaussUtils::ConvolutionX3D(a_o,a_i,
					mdKz.get(),mdSz.get(),mKRadius.z,
					sizeX, sizeY, sizeZ, stream);
    }
    
    CudaUtils::CheckCUDAError(__FILE__, __LINE__);
}

} // end namespace PyCA

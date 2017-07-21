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

#ifndef __GAUSS_REC_FILTER_BASE_H__
#define __GAUSS_REC_FILTER_BASE_H__

#include <Vec3D.h>
#include <estream.h>

#include "GaussRecFilterInterface.h"

namespace PyCA {

template<int mode>
class GaussRecFilterBase
    : public GaussRecFilterInterface
{
public:
   enum { exec_mode = mode };
   
   GaussRecFilterBase();
   virtual ~GaussRecFilterBase(){};
   
   void updateParams(const Vec3Di& size, const Vec3Df& sig);
   
   virtual void filter(float *a_o, const float* a_i, float* a_t, StreamT stream=NULL);
   void convolutionSingleAxis(float* a_o, const float *a_i,
			      size_t sizeX, size_t sizeY, size_t sizeZ,
			      int axis, StreamT stream);
   
   typedef struct {
      //a0-a3, b1, b2, coefp, coefn - filter parameters
      float a0, a1, a2, a3;
      float b1, b2;
      float coefp, coefn;
   } GaussRecParams;
   
protected:
   
   virtual void ConvolutionX3D(float* a_o, const float* a_i, 
			       const GaussRecParams& p,
			       size_t sizeX, size_t sizeY, size_t sizeZ, StreamT stream);
   
   Vec3Di mSize;
   Vec3Df mSigma;
   GaussRecParams mParams[3];
};

} // end namespace PyCA

#endif // __GAUSS_REC_FILTER_BASE_H__

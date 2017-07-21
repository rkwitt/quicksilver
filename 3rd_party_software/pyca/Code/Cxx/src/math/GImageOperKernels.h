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

#ifndef __GIMAGE_OPERS_CU_H__
#define __GIMAGE_OPERS_CU_H__

#include <pycaConst.h>
#include <estream.h>
#include <Vec3D.h>

namespace PyCA {

void 
SubVol(float* d_o,const float* d_i, 
       const Vec3Di& oSize, const Vec3Di& iSize, 
       const Vec3Di& start, StreamT st);


void 
SetSubVol_I(float* d_o,const float* d_i, 
	    const Vec3Di& oSize, const Vec3Di& iSize, 
	    const Vec3Di& start, StreamT st);

void 
Shrink(float *d_o, const float *d_i, 
       const Vec3Di &sz, float c, StreamT st);

void 
SoftAbs(float *d_o, const float *d_i, const Vec3Di &sz,
	float eps, StreamT st);

void 
SoftSgn(float *d_o, const float *d_i, const Vec3Di &sz, 
	float eps, StreamT st);

template<BackgroundStrategy bg, InterpT interp, bool useOriginOffset>
void 
Resample(float *d_o, const Vec3Di &oSz,
	 const float*d_i, const Vec3Di &iSz,
	 StreamT stream);

template<BackgroundStrategy bg, InterpT interp>
void
ResampleWorld(float *d_o, 
	      const Vec3Di &oSz, 
	      const Vec3Df &oSp,
	      const Vec3Df &oOr,
	      const float *d_i, 
	      const Vec3Di &iSz, 
	      const Vec3Df &iSp,
	      const Vec3Df &iOr,
	      StreamT stream);

template<BackgroundStrategy bg>
void
SplatWorld(float *d_o, 
	   const Vec3Di &oSz,
	   const Vec3Df &oSp,
	   const Vec3Df &oOr,
	   const float *d_i, 
	   const Vec3Di &iSz,
	   const Vec3Df &iSp,
	   const Vec3Df &iOr,
	   StreamT stream);

template<BackgroundStrategy bg>
void
SplatWorld(float *d_o, 
	   const Vec3Di &oSz,
	   const Vec3Df &oSp,
	   const Vec3Df &oOr,
	   const float *d_i, 
	   const Vec3Di &iSz,
	   const Vec3Df &iSp,
	   const Vec3Df &iOr,
	   float *d_w,
	   StreamT stream);

void 
Convolve(float *d_o, const float *d_i, 
	 const Vec3Di &sz,
	 const float *d_kernel, 
	 const Vec3Di &kSz,
	 StreamT stream);


} // end namespace PyCA

#endif // __GIMAGE_OPERS_CU_H__

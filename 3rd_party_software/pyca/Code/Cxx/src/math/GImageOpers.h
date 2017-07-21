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

#ifndef __GIMAGE_OPERS_H__
#define __GIMAGE_OPERS_H__

#include <pycaConst.h>
#include <estream.h>
#include <Vec3D.h>

namespace PyCA {

// forward declaration
class Image3D;
class Field3D;

class GImageOpers{
public:

   static void SubVol(float* d_o,const float* d_i, 
		      const Vec3Di& oSize, const Vec3Di& iSize, 
		      const Vec3Di& start, StreamT st = NULL);

   static void SubVol(Image3D& d_o,const Image3D& d_i,
		      const Vec3Di& start, StreamT st=NULL);

   static void SetSubVol_I(Image3D& d_o,const Image3D& d_i,
			   const Vec3Di& start, StreamT st);
   
   static void SetSubVol_I(float* d_o,const float* d_i, 
			   const Vec3Di& oSize, const Vec3Di& iSize, 
			   const Vec3Di& start, StreamT st);

   static void Shrink_I(Image3D& d_o, float c, StreamT st = NULL);

   static void Shrink(Image3D& d_o, const Image3D& d_i, 
		      float c, StreamT st = NULL);

   static void SoftAbs_I(Image3D& d_o, 
			 float eps, StreamT st);
   static void SoftAbs(Image3D& d_o, const Image3D& d_i, 
		       float eps, StreamT st);

   static void SoftSgn_I(Image3D& d_o, 
			 float eps, StreamT st);
   static void SoftSgn(Image3D& d_o, const Image3D& d_i, 
		       float eps, StreamT st);

    template<BackgroundStrategy bg, InterpT interp, bool useOriginOffset>
    static void Resample(Image3D& d_o, const Image3D& d_i, StreamT s=NULL);

   template<BackgroundStrategy bg, InterpT interp> 
   static void ResampleWorld(Image3D& d_o, const Image3D& d_i, StreamT stream=NULL);

   template<BackgroundStrategy bg> 
   static void SplatWorld(Image3D& d_o, const Image3D& d_i, StreamT stream=NULL);

   template<BackgroundStrategy bg> 
   static void SplatWorld(Image3D& d_o, const Image3D& d_i, Image3D& d_w, StreamT stream=NULL);

   static void Convolve(Image3D& d_o, const Image3D& d_i, 
			const Image3D& d_kernel, 
			StreamT stream=NULL);
};

} // end namespace PyCA

#endif // __GIMAGE_OPERS_H__

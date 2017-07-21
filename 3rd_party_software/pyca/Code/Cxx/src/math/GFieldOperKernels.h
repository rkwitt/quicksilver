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

#ifndef __GFIELD_OPER_KERNELS_H
#define __GFIELD_OPER_KERNELS_H

#include <pycaConst.h>
#include <estream.h>
#include <Vec3D.h>

namespace PyCA {

template<bool fwd, BackgroundStrategy bg>
void ComposeVH(float* d_hx, 
	       float* d_hy, 
	       float* d_hz,
               const float* d_vx, 
	       const float* d_vy, 
	       const float* d_vz,
               const float* d_gx, 
	       const float* d_gy, 
	       const float* d_gz,
               const float& delta, 
	       int w, int h, int l,
               float spX, float spY, float spZ, 
	       StreamT stream, bool onDev);

template<bool fwd, BackgroundStrategy bg>
void ComposeHV(float* d_hx, 
	       float* d_hy, 
	       float* d_hz,
               const float* d_gx, 
	       const float* d_gy, 
	       const float* d_gz,
               const float* d_vx, 
	       const float* d_vy, 
	       const float* d_vz,
               const float& delta,
               int w, int h, int l,
               float spX, float spY, float spZ, 
	       StreamT stream, bool onDev);

template<BackgroundStrategy bg>
void
ComposeTranslation(float *d_ox, 
		   float *d_oy, 
		   float *d_oz, 
		   const float *d_ix, 
		   const float *d_iy, 
		   const float *d_iz, 
		   const Vec3Di& sz,
		   const Vec3Df& t, 
		   StreamT stream, 
		   bool onDev) ;

template<BackgroundStrategy bg>
void
ApplyH(float *d_ox, 
       float *d_oy, 
       float *d_oz, 
       const float *d_ix, 
       const float *d_iy, 
       const float *d_iz, 
       const float *d_hx, 
       const float *d_hy, 
       const float *d_hz, 
       const Vec3Di &sz,
       StreamT stream);

template<bool fwd, BackgroundStrategy bg>
void 
ApplyV(float *d_ox, 
       float *d_oy, 
       float *d_oz, 
       const float *d_ix, 
       const float *d_iy, 
       const float *d_iz, 
       const float *d_ux, 
       const float *d_uy, 
       const float *d_uz, 
       const Vec3Di &sz,
       const Vec3Df &sp,
       const float& delta, 
       StreamT stream, bool onDev);

template<BackgroundStrategy bg,  bool rescaleVector>
void
Resample(float *d_ox, 
	 float *d_oy, 
	 float *d_oz, 
	 const Vec3Di &oSz,
	 const float *d_ix, 
	 const float *d_iy, 
	 const float *d_iz, 
	 const Vec3Di &iSz,
	 StreamT stream);

void
ReprojectToUnitVec(float *d_ox, 
		   float *d_oy, 
		   float *d_oz, 
		   const Vec3Di &sz,
		   StreamT st);

void
NormalizeSafe(float *d_ox, 
	      float *d_oy, 
	      float *d_oz, 
	      const float *d_ix, 
	      const float *d_iy, 
	      const float *d_iz, 
	      const Vec3Di &sz,
	      const float& eps, 
	      StreamT st);

void
Shrink(float *d_ox, 
       float *d_oy, 
       float *d_oz, 
       const float *d_ix, 
       const float *d_iy, 
       const float *d_iz, 
       const Vec3Di &sz,
       const float& eps, 
       StreamT st);

template<BackgroundStrategy bg> 
void
FixedPointInverse(float *ginvx, 
		  float *ginvy, 
		  float *ginvz, 
		  const float *gx, 
		  const float *gy, 
		  const float *gz, 
		  const Vec3Di &sz,
		  unsigned int numIter, 
		  StreamT stream, bool onDev);

void
updateInverseSubFromIndentity(float *d_hx,
			      float *d_hy,
			      float *d_hz,
			      const Vec3Di &sz,
			      const Vec3Df &sp,
			      StreamT stream);

template<BackgroundStrategy bg> 
void
Ad(float *Zx, 
   float *Zy, 
   float *Zz, 
   const float *gx, 
   const float *gy, 
   const float *gz, 
   const float *Xx,
   const float *Xy,
   const float *Xz,
   const Vec3Di &sz,
   const Vec3Df &sp,
   StreamT s,bool onDev);

void
AdInf(float *Zx, 
      float *Zy, 
      float *Zz, 
      const float *Xx, 
      const float *Xy, 
      const float *Xz, 
      const float *Yx,
      const float *Yy,
      const float *Yz,
      const Vec3Di &sz,
      const Vec3Df &sp,
      StreamT s,bool onDev);

void
JacobianXY(float *Zx, 
	   float *Zy, 
	   float *Zz, 
	   const float *Xx, 
	   const float *Xy, 
	   const float *Xz, 
	   const float *Yx,
	   const float *Yy,
	   const float *Yz,
	   const Vec3Di &sz,
	   const Vec3Df &sp,
	   StreamT s, bool onDev);

void
JacobianXtY(float *Zx, 
	   float *Zy, 
	   float *Zz, 
	   const float *Xx, 
	   const float *Xy, 
	   const float *Xz, 
	   const float *Yx,
	   const float *Yy,
	   const float *Yz,
	   const Vec3Di &sz,
	   const Vec3Df &sp,
	   StreamT s, bool onDev);

template<BackgroundStrategy bg> 
void
CoAd(float *nx, 
     float *ny, 
     float *nz, 
     const float *gx, 
     const float *gy, 
     const float *gz, 
     const float *mx,
     const float *my,
     const float *mz,
     const Vec3Di &sz,
     const Vec3Df &sp,
     StreamT s,bool onDev);

void
CoAdInf(float *nx, 
	float *ny, 
	float *nz, 
	const float *Xx, 
	const float *Xy, 
	const float *Xz, 
	const float *mx,
	const float *my,
	const float *mz,
	const Vec3Di &sz,
	const Vec3Df &sp,
	StreamT s,bool onDev);

void
DivergenceTensor(float *Zx, 
		 float *Zy, 
		 float *Zz, 
		 const float *Xx, 
		 const float *Xy, 
		 const float *Xz, 
		 const float *Yx,
		 const float *Yy,
		 const float *Yz,
		 const Vec3Di &sz,
		 const Vec3Df &sp,
		 StreamT s,bool onDev);
} // end namespace PyCA

#endif // __GFIELD_OPER_KERNELS_H

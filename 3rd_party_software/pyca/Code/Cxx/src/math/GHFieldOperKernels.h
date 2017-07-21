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

#ifndef __GHFIELD_OPER_KERNELS_H
#define __GHFIELD_OPER_KERNELS_H

#include <estream.h>
#include <pycaConst.h>
#include <Vec3D.h>
#include <Aff3D.h>

namespace PyCA {

void
SetToIdentity(float *d_hx, 
	      float *d_hy, 
	      float *d_hz, 
	      const Vec3Di &sz,
	      StreamT stream);

void
toV(float *vx, 
    float *vy, 
    float *vz, 
    const float *hx, 
    const float *hy, 
    const float *hz, 
    const Vec3Di &sz,
    const Vec3Df &sp,
    StreamT stream);

void
toV_I(float *vx, 
      float *vy, 
      float *vz, 
      const Vec3Di &sz,
      const Vec3Df &sp,
      StreamT stream);

void
toH_I(float *vx, 
      float *vy, 
      float *vz, 
      const Vec3Di &sz,
      StreamT stream);

void
InverseZerothOrder(float *a_hInvx, 
		   float *a_hInvy, 
		   float *a_hInvz, 
		   const float *a_hx, 
		   const float *a_hy, 
		   const float *a_hz, 
		   const Vec3Di &sz,
		   StreamT stream);

void
InverseZerothOrder_I(float *a_hInvx, 
		     float *a_hInvy, 
		     float *a_hInvz, 
		     const Vec3Di &sz,
		     StreamT stream);

void
initializeFromAffine(float *d_hx, 
		     float *d_hy, 
		     float *d_hz, 
		     const Vec3Di &sz,
		     const Vec3Df &sp,
		     const Vec3Df &org,
		     const Aff3Df &aff_in,
		     bool invertAff,
		     StreamT stream);

} // end namespace PyCA
    
#endif // __GHFIELD_OPER_KERNELS_H

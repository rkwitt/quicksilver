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

#ifndef __GVFIELD_OPER_KERNELS_H
#define __GVFIELD_OPER_KERNELS_H

#include <estream.h>
#include <pycaConst.h>

namespace PyCA {

void 
SetToZero(float *d_vx, float *d_vy, float *d_vz, 
	  int sz_x, int sz_y, int sz_z,
	  StreamT stream);

void
toH(float *d_hx, float *d_hy, float *d_hz, 
    float sp_x, float sp_y, float sp_z,
    const float *d_vx, const float *d_vy, const float *d_vz,
    int sz_x, int sz_y, int sz_z,
    const float& delta, 
    StreamT stream, bool onDev);

void
toH_I(float *d_hx, float *d_hy, float *d_hz, 
      float sp_x, float sp_y, float sp_z,
      int sz_x, int sz_y, int sz_z,
      const float& delta, 
      StreamT stream, bool onDev);

    // these are called internally by Splat depending on 'normalize'
    // parameter, should not be called directly by users
void 
splatUnnormalized(float* d_iwdx, float* d_iwdy, float* d_iwdz, 
		  size_t sizeX, size_t sizeY, size_t sizeZ,
		  const float* d_wx, const float* d_wy, const float* d_wz,
		  const float* d_px, const float* d_py, const float* d_pz,
		  size_t nP, StreamT stream);
    
void 
splatNormalized(float* d_iwdx, float* d_iwdy, float* d_iwdz,
		int *i_dd,
		size_t sizeX, size_t sizeY, size_t sizeZ,
		const float* d_wx, const float* d_wy, const float* d_wz,
		const float* d_px, const float* d_py, const float* d_pz,
		size_t nP, StreamT stream);
    
} // end namespace PyCA

#endif // __GVFIELD_OPER_KERNELS_H


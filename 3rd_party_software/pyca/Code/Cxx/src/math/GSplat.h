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

#ifndef __GPU_SPLAT_H
#define __GPU_SPLAT_H

#include <estream.h>
#include <pycaUtils.h>

namespace PyCA {
  namespace Splatting{

    void splatDistance(int* d_id, size_t sizeX, size_t sizeY, size_t sizeZ,
		       const float* d_px , const float* d_py, const float* d_pz,
		       size_t nP, StreamT stream);

    void splat3D(int* d_iwd,
		 size_t sizeX, size_t sizeY, size_t sizeZ,
		 const float* d_w,
		 const float* d_px , const float* d_py, const float* d_pz,
		 size_t nP, StreamT stream);

    void splat3D(int* d_iwdx, int* d_iwdy, int* d_iwdz, 
		 size_t sizeX, size_t sizeY, size_t sizeZ,
		 const float* d_wx, const float* d_wy, const float* d_wz,
		 const float* d_px , const float* d_py, const float* d_pz,
		 size_t nP, StreamT stream);
    
    void splat3D(int* d_iwd, int* d_id,
		 uint sizeX, uint sizeY, uint sizeZ,
		 const float* d_w,
		 const float* d_px , const float* d_py, const float* d_pz,
		 uint nP, StreamT stream);

    void splat3D(int* d_iwdx, int* d_iwdy, int* d_iwdz, int* d_id,
		 size_t sizeX, size_t sizeY, size_t sizeZ,
		 const float* d_wx, const float* d_wy, const float* d_wz,
		 const float* d_px , const float* d_py, const float* d_pz,
		 size_t nP, StreamT stream);

    void convertWeightedDistance_I(float* d_fwd, const int* d_id, size_t n, StreamT stream);
    void convertWeightedDistance_I(float* d_fx, float* d_fy, float* d_fz, const int* d_id, size_t n, StreamT stream);
    
    void convertWeightedDistance(float* d_fwd, const int* d_iwd,
				 const int* d_id, size_t n, StreamT stream);

    void FixedToFloating(float* d_o, const int* d_i, size_t n, StreamT stream);
    void FixedToFloating(float* d_ox, float* d_oy, float* d_oz,
			 const int* d_ix, const int* d_iy,const int* d_iz, size_t n, StreamT stream);

    void FixedToFloating_I(float* d_o, size_t n, StreamT stream);
    void FixedToFloating_I(float* d_ox, float* d_oy, float* d_oz, size_t n, StreamT stream);

  }// end namespace Splatting
} // end namespace PyCA

#endif

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

#ifndef __GIMAGE_FIELD_OPER_KERNELS_H
#define __GIMAGE_FIELD_OPER_KERNELS_H

#include<pycaConst.h>
#include<estream.h>
#include<Vec3D.h>

namespace PyCA {

template<BackgroundStrategy bg, InterpT interp>
void ApplyH(float* d_o, const float* d_i,
	    const float* d_hx, const float* d_hy, const float* d_hz,
            int sizeX, int sizeY, int sizeZ, StreamT stream);

template<bool fwd, BackgroundStrategy bg, InterpT interp>
void ApplyV(float* d_o, const float* d_i,
            const float* d_ux, const float* d_uy, const float* d_uz, const float& delta,
            int sizeX, int sizeY, int sizeZ,
            float spX, float spY, float spZ, StreamT stream, bool onDev);

    
template<BackgroundStrategy bg, InterpT interp>
void ComposeTranslation(float* d_o, const float* d_i,
                        const Vec3Df& t, const Vec3Di& size, StreamT stream, bool onDev);

void
Splat(float *d_o, 
      const float *d_hx, 
      const float *d_hy, 
      const float *d_hz, 
      const float *d_i, 
      Vec3Di sz,
      StreamT stream);

void
SplatAndNormalize(float *d_o, 
		  const float *d_hx, 
		  const float *d_hy, 
		  const float *d_hz, 
		  const float *d_i, 
		  const float *temp,
		  Vec3Di sz,
		  StreamT stream);

void
g_finiteDiff(float* d_o, const float* d_i,
	     int szX, int szY, int szZ,
	     float spX, float spY, float spZ, 
	     DimT dim, DiffT diffType, 
	     enum BoundaryCondT bc, 
	     bool accum, OpT op,
	     StreamT stream);

void 
UpwindDiff(float *d_o, const float *d_i, 
	   const float *d_speed,
	   const Vec3Di& sz, 
	   const Vec3Df& sp,
	   DimT dim,
	   StreamT stream);

void 
UpwindGradMag(float *d_o, 
	      const float *d_i, 
	      const float *d_speed, 
	      const Vec3Di& sz,
	      const Vec3Df& sp,
	      StreamT stream);

template <enum DiffT diffType>
void 
g_gradient(float* d_ox, float *d_oy, float* d_oz,
	   const float* d_i,
	   int szX, int szY, int szZ,
	   float spX, float spY, float spZ,
	   BoundaryCondT bc, StreamT stream);

template <enum DiffT diffType>
void 
g_gradientMask(float* d_ox, float *d_oy, float* d_oz,
	       const float* d_i, const float* d_mask,
	       int szX, int szY, int szZ,
	       float spX, float spY, float spZ,
	       BoundaryCondT bc, StreamT stream);

template<typename T, enum DiffT diffType>
void g_gradientMag(float* d_o,
		   const T* d_i,
		   int sizeX, int sizeY, int sizeZ,
		   float spX, float spY, float spZ, 
		   BoundaryCondT bc, StreamT stream);

template <enum DiffT diffType>
void 
g_divergence(float* d_o,
	     const float* d_ix, const float* d_iy, const float* d_iz,
	     int   sizeX, int sizeY, int sizeZ,
	     float spX, float spY, float spZ,
	     BoundaryCondT bc, StreamT stream);

void 
Magnitude(float *d_o, 
	  const float *d_ix, 
	  const float *d_iy, 
	  const float *d_iz,
	  size_t n,
	  StreamT stream);

/**
 * Compute the magnitude array
 * d_o[i] = d_i[i].x^2 + d_i[i].y^2 + d_i[i].z^2 
 */
void 
SqrMagnitude(float *d_o, 
	     const float *d_ix, 
	     const float *d_iy, 
	     const float *d_iz, 
	     size_t n,
	     StreamT stream);

void
ComponentDotProd(float *d_o, 
		 const float *d_ix, 
		 const float *d_iy, 
		 const float *d_iz, 
		 const float *d_i1x, 
		 const float *d_i1y, 
		 const float *d_i1z, 
		 size_t n,
		 StreamT stream);

void 
Add(float *d_ox, 
    float *d_oy, 
    float *d_oz, 
    const float *d_ix, 
    const float *d_iy, 
    const float *d_iz, 
    const float *d_i1, 
    size_t n,
    StreamT stream);

void 
Sub(float *d_ox, 
    float *d_oy, 
    float *d_oz, 
    const float *d_ix, 
    const float *d_iy, 
    const float *d_iz, 
    const float *d_i1, 
    size_t n,
    StreamT stream);

void 
Mul(float *d_ox, 
    float *d_oy, 
    float *d_oz, 
    const float *d_ix, 
    const float *d_iy, 
    const float *d_iz, 
    const float *d_i1, 
    size_t n,
    StreamT stream);

void 
Div(float *d_ox, 
    float *d_oy, 
    float *d_oz, 
    const float *d_ix, 
    const float *d_iy, 
    const float *d_iz, 
    const float *d_i1, 
    size_t n,
    StreamT stream);

void
Add_I(float *d_ox, 
      float *d_oy, 
      float *d_oz, 
      const float *d_i, 
      size_t n,
      StreamT stream);

void
Sub_I(float *d_ox, 
      float *d_oy, 
      float *d_oz, 
      const float *d_i, 
      size_t n,
      StreamT stream);

void
Mul_I(float *d_ox, 
      float *d_oy, 
      float *d_oz, 
      const float *d_i, 
      size_t n,
      StreamT stream);
void
Div_I(float *d_ox, 
      float *d_oy, 
      float *d_oz, 
      const float *d_i, 
      size_t n,
      StreamT stream);

/** @brief d_o = d_i + d_i1 * d_i2 */
void
Add_Mul(float *d_ox, 
	float *d_oy, 
	float *d_oz, 
	const float *d_ix, 
	const float *d_iy, 
	const float *d_iz, 
	const float *d_i1x, 
	const float *d_i1y, 
	const float *d_i1z, 
	const float *d_i2, 
	size_t n,
	StreamT stream);

/** @brief d_o = d_o + d_i * d_i1 */
void
Add_Mul_I(float *d_ox, 
	  float *d_oy, 
	  float *d_oz, 
	  const float *d_ix, 
	  const float *d_iy, 
	  const float *d_iz, 
	  const float *d_i1, 
	  size_t n,
	  StreamT stream);

/** @brief d_o = d_i + d_i1 * d_i2 */
void
Sub_Mul(float *d_ox, 
	float *d_oy, 
	float *d_oz, 
	const float *d_ix, 
	const float *d_iy, 
	const float *d_iz, 
	const float *d_i1x, 
	const float *d_i1y, 
	const float *d_i1z, 
	const float *d_i2, 
	size_t n,
	StreamT stream);

/** @brief d_o = d_o - d_i * d_i1 */
void
Sub_Mul_I(float *d_ox, 
	  float *d_oy, 
	  float *d_oz, 
	  const float *d_ix, 
	  const float *d_iy, 
	  const float *d_iz, 
	  const float *d_i1, 
	  size_t n,
	  StreamT stream);

/** @brief d_o = d_i * d_i1 * c (d_o.x = d_i.x * d_i1 * c)*/
void
MulMulC(float *d_ox, 
	float *d_oy, 
	float *d_oz, 
	const float *d_ix, 
	const float *d_iy, 
	const float *d_iz, 
	const float *d_i1, 
	const float& c, 
	size_t n,
	StreamT stream, bool onDev);

/** @brief d_o = d_o * d_i * c  (d_o.x = d_o.x * d_i * c)*/
void
MulMulC_I(float *d_ox, 
	  float *d_oy, 
	  float *d_oz, 
	  const float *d_i, 
	  const float& c, 
	  size_t n,
	  StreamT stream, bool onDev);


/** @brief d_o = d_i + d_i1 * d_i2 * c */
void
Add_MulMulC(float *d_ox, 
	    float *d_oy, 
	    float *d_oz, 
	    const float *d_ix, 
	    const float *d_iy, 
	    const float *d_iz, 
	    const float *d_i1x, 
	    const float *d_i1y, 
	    const float *d_i1z, 
	    const float *d_i2, 
	    const float& c,  
	    size_t n,
	    StreamT stream, bool onDev);

/** @brief d_o = d_o + d_i * d_i1 * c */
void 
Add_MulMulC_I(float *d_ox, 
	      float *d_oy, 
	      float *d_oz, 
	      const float *d_ix, 
	      const float *d_iy, 
	      const float *d_iz, 
	      const float *d_i1, 
	      const float& c, 
	      size_t n,
	      StreamT stream, 
	      bool onDev);


void
JacDetH(float *d_detJ,
        const float *d_Xgx, 
        const float *d_Xgy, 
        const float *d_Xgz, 
	const float *d_Ygx, 
	const float *d_Ygy, 
	const float *d_Ygz, 
	const float *d_Zgx, 
	const float *d_Zgy, 
	const float *d_Zgz, 
	size_t n,
	bool slice,
	StreamT stream);

template<bool fwd>
void
JacDetV(float *d_detJ,
        const float *d_Xgx, 
        const float *d_Xgy, 
        const float *d_Xgz, 
	const float *d_Ygx, 
	const float *d_Ygy, 
	const float *d_Ygz, 
	const float *d_Zgx, 
	const float *d_Zgy, 
	const float *d_Zgz, 
	size_t n,
	StreamT stream);

void
JacDetH(float *d_jdet, 
	  const float *d_hx, 
	  const float *d_hy, 
	  const float *d_hz, 
	  Vec3Di sz,
	  Vec3Df sp,
	  DiffT diffType, 
	  BoundaryCondT bc, 
	  StreamT stream);

} // end namespace PyCA

#endif // __GIMAGE_FIELD_OPER_KERNELS_H


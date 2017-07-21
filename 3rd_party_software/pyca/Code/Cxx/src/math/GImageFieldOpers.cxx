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

#include "GImageFieldOpers.h"
#include "GImageFieldOperKernels.h"
#include <pycaUtils.h>
#include "interp.h"
#include <GMemOpers.h>
#include <Image3D.h>
#include <Field3D.h>

#include "FiniteDiff.h"

#include "CudaUtils.h"

namespace PyCA {

template<BackgroundStrategy bg, InterpT interp>
void GImageFieldOpers::
ApplyH(Image3D& d_o, const Image3D& d_i, const Field3D& d_h, StreamT stream){
    MK_CHECK_IMAGE_BACKGROUND(bg);
    MK_CHECK3_SIZE(d_o, d_i, d_h);

    Vec3Di size = d_o.grid().size();
    PyCA::ApplyH<bg, interp>
	(d_o.get(), d_i.get(), d_h.x, d_h.y, d_h.z,
	 size.x, size.y, size.z,  stream);
}

/////////////////////////////////////////////////////////////////////////////
// apply uField to an image
// defImage(x) = image(x + delta * u(x))
/////////////////////////////////////////////////////////////////////////////

template<BackgroundStrategy bg, InterpT interp>
void GImageFieldOpers::
ApplyV(Image3D& d_o, const Image3D& d_i, const Field3D& u, const float& delta, StreamT s,
       bool onDev)
{
    MK_CHECK3_SIZE(d_o, d_i, u);

    Vec3Di size = d_o.grid().size();
    Vec3Df sp   = d_o.grid().spacing();
    PyCA::ApplyV<true, bg, interp>
	(d_o.get(), d_i.get(),
	 u.x, u.y, u.z, delta,
	 size.x, size.y, size.z, 
	 sp.x, sp.y, sp.z, s, onDev);
        
}

/////////////////////////////////////////////////////////////////////////////
// apply uField to an image
// defImage(x) = image(x - delta * u(x))
/////////////////////////////////////////////////////////////////////////////

template<BackgroundStrategy bg, InterpT interp>
void GImageFieldOpers::
ApplyVInv(Image3D& d_o, const Image3D& d_i, const Field3D& u, const float& delta, StreamT s,
          bool onDev)
{
    MK_CHECK3_SIZE(d_o, d_i, u);
    
    Vec3Di size = d_o.grid().size();
    Vec3Df sp   = d_o.grid().spacing();
    PyCA::ApplyV<false, bg, interp>
	(d_o.get(), d_i.get(),
	 u.x, u.y, u.z, delta,
	 size.x, size.y, size.z, 
	 sp.x, sp.y, sp.z, s, onDev);
};

template<BackgroundStrategy bg, InterpT interp>
void GImageFieldOpers::
ComposeTranslation(Image3D& d_o, const Image3D& d_i, const Vec3Df& t, StreamT s, bool onDev){
    MK_CHECK_IMAGE_BACKGROUND(bg);
    MK_CHECK2_SIZE(d_o, d_i);
    PyCA::ComposeTranslation<bg, interp>(d_o.get(), d_i.get(), t, d_o.size(), s, onDev);
}


void GImageFieldOpers::
Splat(Image3D& d_o, const Field3D& d_h, const Image3D& d_i, bool normalize, StreamT stream)
{
   GImageFieldOpers::Splat(d_o.get(), d_h, d_i.get(), normalize, stream);
}

void GImageFieldOpers::
Splat(float *d_o, const Field3D& d_h, const float *d_i, bool normalize, StreamT stream)
{
   Vec3Di sz = d_h.size();
   if (normalize) {
       Image3D temp(d_h.grid(), MEM_DEVICE);
       PyCA::SplatAndNormalize(d_o, 
			       d_h.x, d_h.y, d_h.z,
			       d_i, temp.get(),
			       sz, stream);
       CudaUtils::CheckCUDAError(__FILE__,__LINE__,"Error after SplatAndNormalize");
   } else {
       PyCA::Splat(d_o, 
		   d_h.x, d_h.y, d_h.z,
		   d_i,
		   sz, stream);
       CudaUtils::CheckCUDAError(__FILE__,__LINE__,"Error after Splat");
   }
}

#define SLICE_TRUE 1
#define SLICE_FALSE 0

void GImageFieldOpers
::FiniteDiff(Image3D& d_o, const Image3D& d_i, 
	     DimT dim, DiffT diffType, 
	     enum BoundaryCondT bc, 
	     bool accum, OpT op,
	     StreamT stream)
{
   MK_CHECK2_SIZE(d_o, d_i);
   Vec3Di sz = d_o.size();
   Vec3Df sp = d_o.spacing();

   PyCA::g_finiteDiff(d_o.get(), d_i.get(), 
		    sz.x, sz.y, sz.z, 
		    sp.x, sp.y, sp.z,
		    dim, diffType,
		    bc, accum, op, stream);

}

void GImageFieldOpers
::UpwindDiff(Image3D& d_o, const Image3D& d_i, 
	     const Image3D& d_speed,
	     DimT dim,
	     StreamT stream)
{
   MK_CHECK3_SIZE(d_o, d_i, d_speed);
   Vec3Di sz = d_i.size();
   Vec3Df sp = d_i.spacing();

   PyCA::UpwindDiff(d_o.get(), d_i.get(), d_speed.get(), 
		    sz, sp, dim, stream);

}

void GImageFieldOpers
::UpwindGradMag(Image3D& d_o, 
		const Image3D& d_i, 
		const Image3D& d_speed, 
		StreamT stream)
{
   MK_CHECK3_SIZE(d_o, d_i, d_speed);
   Vec3Di sz = d_i.size();
   Vec3Df sp = d_i.spacing();

   PyCA::UpwindGradMag(d_o.get(), d_i.get(), d_speed.get(),
		       sz, sp, stream);
}

void GImageFieldOpers
::Gradient(Field3D& d_o, const float* d_i, DiffT diffType,
	   BoundaryCondT bc, StreamT stream)
{
    Vec3Di size = d_o.size();
    Vec3Df sp   = d_o.spacing();
    if(diffType == DIFF_FORWARD){
       PyCA::g_gradient<DIFF_FORWARD>
	  (d_o.x, d_o.y, d_o.z, d_i,
	   size.x, size.y, size.z,
	   sp.x, sp.y, sp.z, bc, stream);
    }else if(diffType == DIFF_BACKWARD){
       PyCA::g_gradient<DIFF_BACKWARD>
	  (d_o.x, d_o.y, d_o.z, d_i,
	   size.x, size.y, size.z,
	   sp.x, sp.y, sp.z, bc, stream);
    }else if(diffType == DIFF_CENTRAL){
       PyCA::g_gradient<DIFF_CENTRAL>
	  (d_o.x, d_o.y, d_o.z, d_i,
	   size.x, size.y, size.z,
	   sp.x, sp.y, sp.z, bc, stream);
    }else{
      throw PyCAException(__FILE__, __LINE__, "unknown DiffT");
    }
}


void GImageFieldOpers::
Gradient(Field3D& d_o, const Image3D& d_i, DiffT diffType, 
	 BoundaryCondT bc, StreamT stream)
{
    MK_CHECK2_SIZE(d_o, d_i);
    GImageFieldOpers::Gradient(d_o, d_i.get(), diffType, bc, stream);
}

void GImageFieldOpers::
Gradient2(Field3D& d_o, const Image3D& d_i, DiffT diffType, BoundaryCondT bc, StreamT stream){
   throw PyCAException(__FILE__, __LINE__, "GPU Gradient2 unimplemented");
}

void GImageFieldOpers::
GradFor(Field3D& d_o, const Image3D& d_i, BoundaryCondT bc, StreamT stream){
    MK_CHECK2_SIZE(d_o, d_i);
    Vec3Di size = d_o.size();
    Vec3Df sp   = d_o.spacing();
    PyCA::g_gradient<DIFF_FORWARD>
       (d_o.x, d_o.y, d_o.z, d_i.get(), 
	size.x, size.y, size.z,
	sp.x, sp.y, sp.z,
	bc, stream);
}

void GImageFieldOpers
::GradientMask(Field3D& d_o, 
	       const Image3D& d_i, const Image3D& d_mask,
	       DiffT diffType, BoundaryCondT bc, 
	       StreamT stream)
{
    Vec3Di size = d_o.size();
    Vec3Df sp   = d_o.spacing();
    if(diffType == DIFF_FORWARD){
       PyCA::g_gradientMask<DIFF_FORWARD>
	  (d_o.x, d_o.y, d_o.z, 
	   d_i.get(), d_mask.get(),
	   size.x, size.y, size.z,
	   sp.x, sp.y, sp.z, bc, stream);
    }else if(diffType == DIFF_BACKWARD){
       PyCA::g_gradientMask<DIFF_BACKWARD>
	  (d_o.x, d_o.y, d_o.z, 
	   d_i.get(), d_mask.get(),
	   size.x, size.y, size.z,
	   sp.x, sp.y, sp.z, bc, stream);
    }else if(diffType == DIFF_CENTRAL){
       PyCA::g_gradientMask<DIFF_CENTRAL>
	  (d_o.x, d_o.y, d_o.z, 
	   d_i.get(), d_mask.get(),
	   size.x, size.y, size.z,
	   sp.x, sp.y, sp.z, bc, stream);
    }else{
      throw PyCAException(__FILE__, __LINE__, "unknown DiffT");
    }
}

void GImageFieldOpers::
GradientMag(Image3D& d_o, const Image3D& d_i, DiffT diffType, 
	    BoundaryCondT bc, StreamT stream)
{

    MK_CHECK2_SIZE(d_o, d_i);

    Vec3Di size = d_o.size();
    Vec3Df sp   = d_o.spacing();
    
    if(diffType == DIFF_FORWARD){
       PyCA::g_gradientMag<float,DIFF_FORWARD>
	  (d_o.get(), d_i.get(),
	   size.x,size.y,size.z,
	   sp.x,sp.y,sp.z, 
	   bc, stream);
    }else if(diffType == DIFF_BACKWARD){
       PyCA::g_gradientMag<float,DIFF_BACKWARD>
	  (d_o.get(), d_i.get(),
	   size.x,size.y,size.z,
	   sp.x,sp.y,sp.z, 
	   bc, stream);
    }else if(diffType == DIFF_CENTRAL){
       PyCA::g_gradientMag<float,DIFF_CENTRAL>
	  (d_o.get(), d_i.get(),
	   size.x,size.y,size.z,
	   sp.x,sp.y,sp.z, 
	   bc, stream);
    }else{
      throw PyCAException(__FILE__, __LINE__, "Unknown DiffT");
    }
}

void GImageFieldOpers::
GradForMag(Image3D& d_o, const Image3D& d_i, StreamT stream){
    MK_CHECK2_SIZE(d_o, d_i);

    Vec3Di size = d_o.size();
    Vec3Df sp   = d_o.spacing();
    
    PyCA::g_gradientMag<float,DIFF_FORWARD>
       (d_o.get(), d_i.get(),
	size.x,size.y,size.z,
	sp.x,sp.y,sp.z,
	BC_APPROX,stream);
}

void GImageFieldOpers::
Divergence(Image3D& d_o, const Field3D& d_i, DiffT diffType, 
	   BoundaryCondT bc, StreamT stream)
{
   MK_CHECK2_SIZE(d_o, d_i);
   
   Vec3Di size = d_o.size();
   Vec3Df sp   = d_o.spacing();

   if(diffType == DIFF_FORWARD){
      PyCA::g_divergence<DIFF_FORWARD>
	 (d_o.get(),
	  d_i.x, d_i.y, d_i.z,
	  size.x, size.y, size.z,
	  sp.x, sp.y, sp.z, bc, stream);
   }else if(diffType == DIFF_BACKWARD){
      PyCA::g_divergence<DIFF_BACKWARD>
	 (d_o.get(),
	  d_i.x, d_i.y, d_i.z,
	  size.x, size.y, size.z,
	  sp.x, sp.y, sp.z, bc, stream);
   }else if(diffType == DIFF_CENTRAL){
      PyCA::g_divergence<DIFF_CENTRAL>
	 (d_o.get(),
	  d_i.x, d_i.y, d_i.z,
	  size.x, size.y, size.z,
	  sp.x, sp.y, sp.z, bc, stream);
   }else{
      throw PyCAException(__FILE__, __LINE__, "Unknown DiffT");
   }
   
}

void GImageFieldOpers::
DivBack(Image3D& d_o, const Field3D& d_i, BoundaryCondT bc, StreamT stream){
   MK_CHECK2_SIZE(d_o, d_i);
   
   Vec3Di size = d_o.size();
   Vec3Df sp   = d_o.spacing();
   
   PyCA::g_divergence<DIFF_BACKWARD>
      (d_o.get(),
       d_i.x, d_i.y, d_i.z,
       size.x, size.y, size.z,
       sp.x, sp.y, sp.z, bc, stream);
}

/**
 * Compute the magnitude image
 * d_o[i] = sqrt(d_i[i].x^2 + d_i[i].y^2 + d_i[i].z^2) 
 */
void GImageFieldOpers::
Magnitude(Image3D& d_o, const Field3D& d_i, StreamT stream){
    MK_CHECK2_SIZE(d_o, d_i);
    size_t n = d_o.nVox();

    PyCA::Magnitude
	(d_o.get(), d_i.x, d_i.y, d_i.z, n, stream);
}

/**
 * Compute the magnitude array
 * d_o[i] = d_i[i].x^2 + d_i[i].y^2 + d_i[i].z^2 
 */
void GImageFieldOpers::
SqrMagnitude(Image3D& d_o, const Field3D& d_i, StreamT stream){
    MK_CHECK2_SIZE(d_o, d_i);
    size_t n = d_o.nVox();
    PyCA::SqrMagnitude
	(d_o.get(), d_i.x, d_i.y, d_i.z, n, stream);
}

void GImageFieldOpers::
ComponentDotProd(Image3D& d_o, const Field3D& d_i, const Field3D& d_i1, StreamT stream){
    MK_CHECK3_SIZE(d_o, d_i, d_i1);
    size_t n = d_o.nVox();

    PyCA::ComponentDotProd(d_o.get(),
			   d_i.getX(), d_i.getY(), d_i.getZ(), 
			   d_i1.getX(), d_i1.getY(), d_i1.getZ(), 
			   n, stream);
}

/** @brief d_o.x = d_i.x + d_i1, d_o.y = d_i.y + d_i1, d_o.z = d_i.z + d_i1 */
void GImageFieldOpers::
Add(Field3D& d_o, const Field3D& d_i, const Image3D& d_i1, StreamT stream){
    MK_CHECK3_SIZE(d_o, d_i, d_i1);
    size_t n = d_o.nVox();
    PyCA::Add(d_o.x, d_o.y, d_o.z,
	      d_i.x, d_i.y, d_i.z,
	      d_i1.get(),
	      n, stream);
}

/** @brief d_o.x = d_i.x - d_i1, d_o.y = d_i.y - d_i1, d_o.z = d_i.z - d_i1 */
void GImageFieldOpers::
Sub(Field3D& d_o, const Field3D& d_i, const Image3D& d_i1, StreamT stream){
    MK_CHECK3_SIZE(d_o, d_i, d_i1);
    size_t n = d_o.nVox();
    PyCA::Sub(d_o.x, d_o.y, d_o.z,
	      d_i.x, d_i.y, d_i.z,
	      d_i1.get(), 
	      n, stream);
}


/** @brief d_o.x = d_i.x * d_i1, d_o.y = d_i.y * d_i1, d_o.z = d_i.z * d_i1 */
void GImageFieldOpers::
Mul(Field3D& d_o, const Field3D& d_i, const Image3D& d_i1, StreamT stream){
    MK_CHECK3_SIZE(d_o, d_i, d_i1);
    size_t n = d_o.nVox();
    PyCA::Mul(d_o.x, d_o.y, d_o.z,
	      d_i.x, d_i.y, d_i.z,
	      d_i1.get(), 
	      n, stream);
}

/** @brief d_o.x = d_i.x / d_i1, d_o.y = d_i.y / d_i1, d_o.z = d_i.z / d_i1 */
void GImageFieldOpers::
Div(Field3D& d_o, const Field3D& d_i, const Image3D& d_i1, StreamT stream){
    MK_CHECK3_SIZE(d_o, d_i, d_i1);
    size_t n = d_o.nVox();
    PyCA::Div(d_o.x, d_o.y, d_o.z,
	      d_i.x, d_i.y, d_i.z,
	      d_i1.get(), 
	      n, stream);
}

/** @brief d_o.x += d_i, d_o.y += d_i, d_o.y += d_i, */
void GImageFieldOpers::
Add_I(Field3D& d_o, const Image3D& d_i, StreamT stream){
    MK_CHECK2_SIZE(d_o, d_i);
    size_t n = d_o.nVox();
    PyCA::Add_I(d_o.x, d_o.y, d_o.z,
		d_i.get(), 
		n, stream);
}

/** @brief d_o.x -= d_i, d_o.y -= d_i, d_o.y -= d_i, */
void GImageFieldOpers::
Sub_I(Field3D& d_o, const Image3D& d_i, StreamT stream){
    MK_CHECK2_SIZE(d_o, d_i);
    size_t n = d_o.nVox();
    PyCA::Sub_I(d_o.x, d_o.y, d_o.z,
		d_i.get(), 
		n, stream);
    
}

void GImageFieldOpers::
Mul_I(Field3D& d_o, const Image3D& d_i, StreamT stream){
    MK_CHECK2_SIZE(d_o, d_i);
    size_t n = d_o.nVox();
    PyCA::Mul_I(d_o.x, d_o.y, d_o.z, 
		d_i.get(), 
		n, stream);
}

void GImageFieldOpers::
Div_I(Field3D& d_o, const Image3D& d_i, StreamT stream){
    MK_CHECK2_SIZE(d_o, d_i);
    size_t n = d_o.nVox();
    PyCA::Div_I(d_o.x, d_o.y, d_o.z,
		d_i.get(), 
		n, stream);
}

/** @brief d_o = d_i + d_i1 * d_i2 */
void GImageFieldOpers::
Add_Mul(Field3D& d_o, const Field3D& d_i, const Field3D& d_i1, const Image3D& d_i2, StreamT stream){
    MK_CHECK4_SIZE(d_o, d_i, d_i1, d_i2);
    size_t n = d_o.nVox();
    PyCA::Add_Mul(d_o.x, d_o.y, d_o.z,
		  d_i.x, d_i.y, d_i.z,
		  d_i1.x, d_i1.y, d_i1.z,
		  d_i2.get(), 
		  n, stream);
}

void GImageFieldOpers::
Add_Mul_I(Field3D& d_o, const Field3D& d_i, const Image3D& d_i1, StreamT stream){
    MK_CHECK3_SIZE(d_o, d_i, d_i1);
    size_t n = d_o.nVox();
    PyCA::Add_Mul_I(d_o.x, d_o.y, d_o.z,
		    d_i.x, d_i.y, d_i.z,
		    d_i1.get(), 
		    n, stream);
}

/** @brief d_o = d_i + d_i1 * d_i2 */
void GImageFieldOpers::
Sub_Mul(Field3D& d_o, const Field3D& d_i, const Field3D& d_i1, const Image3D& d_i2, StreamT stream){
    MK_CHECK4_SIZE(d_o, d_i, d_i1, d_i2);
    size_t n = d_o.nVox();
    PyCA::Sub_Mul(d_o.x, d_o.y, d_o.z,
		  d_i.x, d_i.y, d_i.z,
		  d_i1.x, d_i1.y, d_i1.z,
		  d_i2.get(), 
		  n, stream);
}

void GImageFieldOpers::
Sub_Mul_I(Field3D& d_o, const Field3D& d_i, const Image3D& d_i1, StreamT stream){
    MK_CHECK3_SIZE(d_o, d_i, d_i1);
    size_t n = d_o.nVox();
    PyCA::Sub_Mul_I(d_o.x, d_o.y, d_o.z,
		    d_i.x, d_i.y, d_i.z,
		    d_i1.get(), n, stream);
}

/** @brief d_o = d_i * d_i1 * c (d_o.x = d_i.x * d_i1 * c)*/
void GImageFieldOpers::
MulMulC(Field3D& d_o, const Field3D& d_i, const Image3D& d_i1, const float& c, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(d_o, d_i, d_i1);
    size_t n = d_o.nVox();
    PyCA::MulMulC(d_o.x, d_o.y, d_o.z,
		  d_i.x, d_i.y, d_i.z,
		  d_i1.get(), c, 
		  n, stream, onDev);
}

/** @brief d_o = d_o * d_i * c  (d_o.x = d_o.x * d_i * c)*/
void GImageFieldOpers::
MulMulC_I(Field3D& d_o, const Image3D& d_i, const float& c, StreamT stream, bool onDev){
    MK_CHECK2_SIZE(d_o, d_i);
    size_t n = d_o.nVox();
    PyCA::MulMulC_I(d_o.x, d_o.y, d_o.z, 
		    d_i.get(), c, 
		    n, stream, onDev);
}

/** @brief d_o = d_i + d_i1 * d_i2 * c */
void GImageFieldOpers::
Add_MulMulC(Field3D& d_o, const Field3D& d_i, const Field3D& d_i1, const Image3D& d_i2, const float& c,  StreamT stream, bool onDev){
    MK_CHECK4_SIZE(d_o, d_i, d_i1, d_i2);
    size_t n = d_o.nVox();
    PyCA::Add_MulMulC
	(d_o.x, d_o.y, d_o.z, 
	 d_i.x, d_i.y, d_i.z, 
	 d_i1.x, d_i1.y, d_i1.z,
	 d_i2.get(), c, 
	 n, stream, onDev);
}

void GImageFieldOpers::
Add_MulMulC_I(const Field3D& d_o, const Field3D& d_i, const Image3D& d_i1, const float& c, StreamT stream, bool onDev){

    MK_CHECK3_SIZE(d_o, d_i, d_i1);
    size_t n = d_o.nVox();

    PyCA::Add_MulMulC_I
	(d_o.x, d_o.y, d_o.z, 
	 d_i.x, d_i.y, d_i.z, 
	 d_i1.get(), c, 
	 n, stream, onDev);
}

void GImageFieldOpers::
JacDetH(Image3D& d_detJ,
        const Field3D& d_Xg, const Field3D& d_Yg, const Field3D& d_Zg, StreamT stream)
{
    size_t n = d_detJ.nVox();
    bool slice = (d_Xg.size().z == 1);
    PyCA::JacDetH
	  (d_detJ.get(),
	   d_Xg.x, d_Xg.y, d_Xg.z,
	   d_Yg.x, d_Yg.y, d_Yg.z,
	   d_Zg.x, d_Zg.y, d_Zg.z, 
	   n, slice, stream);
}

void GImageFieldOpers::
JacDetV(Image3D& d_detJ,
        const Field3D& d_Xg, const Field3D& d_Yg, const Field3D& d_Zg, StreamT stream){

    MK_CHECK2_SIZE(d_detJ, d_Xg);
    size_t  n= d_detJ.nVox();
        
    PyCA::JacDetV<true>(d_detJ.get(),
			 d_Xg.x, d_Xg.y, d_Xg.z,
			 d_Yg.x, d_Yg.y, d_Yg.z,
			 d_Zg.x, d_Zg.y, d_Zg.z, 
			 n, stream);
}


void GImageFieldOpers::
JacDetVInv(Image3D& d_detJ,
           const Field3D& d_Xg, const Field3D& d_Yg, const Field3D& d_Zg, StreamT stream){

    MK_CHECK2_SIZE(d_detJ, d_Xg);
    size_t  n= d_detJ.nVox();
        
    PyCA::JacDetV<false>(d_detJ.get(),
			 d_Xg.x, d_Xg.y, d_Xg.z,
			 d_Yg.x, d_Yg.y, d_Yg.z,
			 d_Zg.x, d_Zg.y, d_Zg.z, 
			 n, stream);
}

void GImageFieldOpers
::JacDetH(Image3D& d_jdet, const Field3D& d_h, 
	  DiffT diffType, BoundaryCondT bc, 
	  StreamT stream)
{

    MK_CHECK2_SIZE(d_jdet, d_h);
    Vec3Di sz = d_h.size();
    Vec3Df sp = d_h.spacing();

    PyCA::JacDetH(d_jdet.get(),
		  d_h.x, d_h.y, d_h.z,
		  sz, sp,
		  diffType, bc,
		  stream);

    CudaUtils::CheckCUDAError(__FILE__,__LINE__,"Error after JacDetHPointwise");
}

// template instantiation
#include "GImageFieldOpers_inst.cxx"

    
} // end namespace PyCA

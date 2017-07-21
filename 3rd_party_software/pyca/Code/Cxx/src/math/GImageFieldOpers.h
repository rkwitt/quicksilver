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

#ifndef __GIMAGE_FIELD_OPER_H
#define __GIMAGE_FIELD_OPER_H

#include<pycaConst.h>
#include<estream.h>
#include<Vec3D.h>

namespace PyCA {

class Image3D;
class Field3D;

class GImageFieldOpers{
public:
/*
 * apply hField to an image
 *  defImage(x) = image(h(x))
 */
    template<BackgroundStrategy bg, InterpT interp>
    static void ApplyH(Image3D& d_o, const Image3D& d_i, const Field3D& d_h, StreamT s);

/*
 *  apply uField to an image
 *  defImage(x) = image(x + delta * u(x))
 */
    template<BackgroundStrategy bg, InterpT interp>
    static void ApplyV(Image3D& d_o, const Image3D& d_i, const Field3D& d_u, const float& delta, StreamT s, bool onDev);
    
    /*
     *  apply uField to an image
     *  defImage(x) = image(x - delta * u(x))
     */
    template<BackgroundStrategy bg, InterpT interp>
    static void ApplyVInv(Image3D& d_o, const Image3D& d_i, const Field3D& d_u, const float& delta, StreamT s, bool onDev);

    /*
     * d_o = d_i ( x + t)
     */
    template<BackgroundStrategy bg, InterpT interp>
    static void ComposeTranslation(Image3D& d_o, const Image3D& d_i, const Vec3Df& t, StreamT s, bool onDev);
    
    static void Splat(float *d_o, const Field3D& d_h, const float* d_i, bool normalize, StreamT s);

    static void Splat(Image3D& d_o, const Field3D& d_h, const Image3D& d_i, bool normalize, StreamT s);
    
   /**
    * Compute finite difference in the dimension indicated
    *  d_o = diff(d_i)
    */
   static void
   FiniteDiff(Image3D& d_o, const Image3D& d_i, 
	      DimT dim,  DiffT diffType=DIFF_CENTRAL, 
	      enum BoundaryCondT bc=BC_CLAMP, 
	      bool accum=false, OpT op=OP_VAL,
	      StreamT stream=NULL);
   
/**
 * Compute the gradient of an image
 *  d_o = grad (d_i)
 */
    static void Gradient(Field3D& d_o, const Image3D& d_i, 
			 DiffT diffType=DIFF_CENTRAL, 
			 BoundaryCondT bc=BC_CLAMP, 
			 StreamT s=NULL);

    static void Gradient2(Field3D& d_o, const Image3D& d_i, 
			  DiffT diffType=DIFF_CENTRAL, 
			  BoundaryCondT bc=BC_CLAMP, 
			  StreamT s=NULL);

   /** for taking gradient of vector field component */
   static void Gradient(Field3D& d_o, const float* d_i, 
			DiffT diffType=DIFF_CENTRAL, 
			BoundaryCondT bc=BC_CLAMP, 
			StreamT s=NULL);
			


/**
 * Masked gradient computation
 */
   static 
   void 
   GradientMask(Field3D& d_o, 
		const Image3D& d_i, 
		const Image3D& d_mask,
		DiffT diffType=DIFF_CENTRAL, 
		BoundaryCondT bc=BC_CLAMP, 
		StreamT s=NULL);

   /**
    * Gradient via forward differences
    *  d_o = grad (d_i)
    */
   static void GradFor(Field3D& d_o, const Image3D& d_i, BoundaryCondT bc, StreamT s);
   
   static void GradForMag(Image3D& d_o, const Image3D& d_i, StreamT stream);

    static void GradientMag(Image3D& d_o, const Image3D& d_i, 
			    DiffT diffType=DIFF_CENTRAL, 
			    BoundaryCondT bc=BC_CLAMP, 
			    StreamT s=NULL);


   /**
    * returns deriv(h_i), where deriv is upwind derivative based on
    * d_speed in dimension 'dim'
    */
   static void 
   UpwindDiff(Image3D& d_o, const Image3D& d_i, 
	      const Image3D& d_speed,
	      DimT dim,
	      StreamT stream=NULL);
   
   /**
    * returns |\nabla d_i|*d_speed, where \nabla taken with upwind
    * differences based on d_speed
    */
   static void UpwindGradMag(Image3D& d_o, const Image3D& d_i,
			     const Image3D& d_speed, StreamT stream=NULL);
/**
 * Compute the diverence of a field
 *  d_o = div (d_i)
 */
   static void Divergence(Image3D& d_o, const Field3D& d_i, 
			  DiffT diffType=DIFF_CENTRAL, 
			  BoundaryCondT bc=BC_CLAMP, 
			  StreamT s=NULL);

/**
 * Divergence via backward differences
 *  d_o = div (d_i)
 */
   static void DivBack(Image3D& d_o, const Field3D& d_i, BoundaryCondT bc, StreamT s);

/**
 * Compute the magnitude image
 * d_o[i] = sqrt(d_i[i].x^2 + d_i[i].y^2 + d_i[i].z^2) 
 */
    static void Magnitude(Image3D& d_o, const Field3D& d_i, StreamT s);

/**
 * Compute the magnitude array
 * d_o[i] = d_i[i].x^2 + d_i[i].y^2 + d_i[i].z^2 
 */
    static void SqrMagnitude(Image3D& d_o, const Field3D& d_i, StreamT s);

    /**
     * Compute dot product array 
     * d_o[i] = d_i[i].x * d_i1[i].x + d_i[i].y * d_i1[i].y + d_i[i].z * d_i1[i].z
     */
    static void ComponentDotProd(Image3D& d_o, const Field3D& d_i, const Field3D& d_i1, StreamT s);
    
/** @brief d_o.x = d_i.x + d_i1, d_o.y = d_i.y + d_i1, d_o.z = d_i.z + d_i1 */
    static void Add(Field3D& d_o, const Field3D& d_i, const Image3D& d_i1, StreamT s);

/** @brief d_o.x = d_i.x - d_i1, d_o.y = d_i.y - d_i1, d_o.z = d_i.z - d_i1 */
    static void Sub(Field3D& d_o, const Field3D& d_i, const Image3D& d_i1, StreamT s);

/** @brief d_o.x = d_i.x * d_i1, d_o.y = d_i.y * d_i1, d_o.z = d_i.z * d_i1 */
    static void Mul(Field3D& d_o, const Field3D& d_i, const Image3D& d_i1, StreamT s);

/** @brief d_o.x = d_i.x / d_i1, d_o.y = d_i.y / d_i1, d_o.z = d_i.z / d_i1 */
    static void Div(Field3D& d_o, const Field3D& d_i, const Image3D& d_i1, StreamT s);

/** @brief d_o.x += d_i, d_o.y += d_i, d_o.y += d_i, */
    static void Add_I(Field3D& d_o, const Image3D& d_i, StreamT s);

/** @brief d_o.x -= d_i, d_o.y -= d_i, d_o.y -= d_i, */
    static void Sub_I(Field3D& d_o, const Image3D& d_i, StreamT s);

/** @brief d_o.x *= d_i, d_o.y *= d_i, d_o.y *= d_i, */
    static void Mul_I(Field3D& d_o, const Image3D& d_i, StreamT s);

/** @brief d_o.x /= d_i, d_o.y = d_i.y / d_i, d_o.z = d_i.z / d_i */
    static void Div_I(Field3D& d_o, const Image3D& d_i, StreamT s);

/** @brief d_o = d_i + d_i1 * d_i2 */
    static void Add_Mul(Field3D& d_o, const Field3D& d_i, const Field3D& d_i1, const Image3D& d_i2, StreamT s);

/** @brief d_o = d_o + d_i * d_i1 */
    static void Add_Mul_I(Field3D& d_o, const Field3D& d_i, const Image3D& d_i1, StreamT s);

/** @brief d_o = d_i - d_i1 * d_i2 */
    static void Sub_Mul(Field3D& d_o, const Field3D& d_i, const Field3D& d_i1, const Image3D& d_i2, StreamT s);

/** @brief d_o = d_o - d_i * d_i1 */
    static void Sub_Mul_I(Field3D& d_o, const Field3D& d_i, const Image3D& d_i1, StreamT s);

/** @brief d_o = d_i * d_i1 * c (d_o.x = d_i.x * d_i1 * c)*/
    static void MulMulC(Field3D& d_o, const Field3D& d_i, const Image3D& d_i1, const float& c, StreamT s, bool onDev);

/** @brief d_o = d_o * d_i * c  (d_o.x = d_o.x * d_i * c)*/
    static void MulMulC_I(Field3D& d_o, const Image3D& d_i, const float& c, StreamT s, bool onDev);

/** @brief d_o = d_i + d_i1 * d_i2 * c */
    static void Add_MulMulC(Field3D& d_o, const Field3D& d_i, const Field3D& d_i1, const Image3D& d_i2, const float& c,  StreamT s, bool onDev);

/** @brief d_o = d_o + d_i * d_i1 * c */
    static void Add_MulMulC_I(const Field3D& d_o, const Field3D& d_i, const Image3D& d_i1, const float& c,  StreamT s, bool onDev);


    static void JacDetH(Image3D& d_detJ,
                        const Field3D& d_Xg, const Field3D& d_Yg, const Field3D& d_Zg, StreamT stream);

    static void JacDetV(Image3D& d_detJ,
                        const Field3D& d_Xg, const Field3D& d_Yg, const Field3D& d_Zg, StreamT stream);

    static void JacDetVInv(Image3D& d_detJ,
                           const Field3D& d_Xg, const Field3D& d_Yg, const Field3D& d_Zg, StreamT stream);

   static void JacDetH(Image3D& d_jdet, const Field3D& d_h, 
		       DiffT diffType, BoundaryCondT bc, 
		       StreamT stream);
    
};

} // end namespace PyCA

#endif

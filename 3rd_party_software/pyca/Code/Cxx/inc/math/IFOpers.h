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

#ifndef __IF_OPERS_H
#define __IF_OPERS_H

#ifndef SWIG
#include<pycaConst.h>
#include<estream.h>
#include<Vec3D.h>
#endif // SWIG


namespace PyCA {

class Image3D;
class Field3D;

namespace Opers{

/**
 * Copy a component of a vector field to an image(where dim = 0, 1, or 2)
 */
void Copy(Image3D& a_o, const Field3D& a_i, int dim, StreamT stream=NULL);
/**
 * Copy an image to a component of a vector field (where dim = 0, 1, or 2)
 */
void Copy(Field3D& a_o, const Image3D& a_i, int dim, StreamT stream=NULL);

/**
@brief
    apply hField to an image

    a_o(x) = a_i(h(x))
@param a_o
    calculated deformed image
@param a_i
    input image
@param a_h
    h-field that applies the deformation
*/
template<BackgroundStrategy bg, InterpT interp>
void ApplyH(Image3D& a_o, const Image3D& a_i, const Field3D& a_h, StreamT s=NULL);

// Non-template version
void ApplyH(Image3D& a_o, const Image3D& a_i, const Field3D& a_h, 
            BackgroundStrategy bg = BACKGROUND_STRATEGY_CLAMP,
	    InterpT interp = DEFAULT_INTERP_METHOD,
            StreamT s=NULL);

/**
@brief
    apply uField to an image

    a_o(x) = a_i(x + delta*a_u(x))
@param a_o
    calculated deformed image
@param a_i
    input image
@param a_u
    u-field that applies the deformation
@param delta
    float that scales the deformation (typically want 1.0)
*/
template<BackgroundStrategy bg, InterpT interp>
void ApplyV(Image3D& a_o, const Image3D& a_i,
            const Field3D& a_u, const float& delta=1.0,
            StreamT s=NULL, bool onDev=false);

// Non-template version
void ApplyV(Image3D& a_o, const Image3D& a_i, const Field3D& a_u, const float& delta=1.0, 
            BackgroundStrategy bg = BACKGROUND_STRATEGY_CLAMP,
	    InterpT interp = DEFAULT_INTERP_METHOD,
            StreamT s=NULL, bool onDev=false);
    
template<BackgroundStrategy bg, InterpT interp>
void ApplyVInv(Image3D& a_o, const Image3D& a_i,
               const Field3D& a_u, const float& delta=1.0,
               StreamT s=NULL, bool onDev=false);

// Non-template version
void ApplyVInv(Image3D& a_o, const Image3D& a_i, const Field3D& a_u, const float& delta=1.0, 
               BackgroundStrategy bg = BACKGROUND_STRATEGY_CLAMP,
	       InterpT interp = DEFAULT_INTERP_METHOD,
               StreamT s=NULL, bool onDev=false);

/**
@brief
   Applies translation to an image

   a_o = a_i ( x + t)
@param a_o
   output (transformed) image
@param a_i
   input image
@param t
   translation (Vec3Df)
*/
template<BackgroundStrategy bg, InterpT interp>
void ComposeTranslation(Image3D& a_o, const Image3D& a_i, const Vec3Df& t,
                        StreamT s = NULL, bool onDev = false);
   
// Non-template version
void ComposeTranslation(Image3D& a_o, const Image3D& a_i, const Vec3Df& t, 
                        BackgroundStrategy bg = BACKGROUND_STRATEGY_CLAMP,
                        InterpT interp = DEFAULT_INTERP_METHOD,
                        StreamT s = NULL, bool onDev = false);
    
void Splat(Image3D& a_o, const Field3D& a_h, const Image3D& a_i, bool normalize=false, StreamT s = NULL);
    
/**
* Compute finite difference in the dimension indicated
*  a_o = diff(a_i)
*/
void FiniteDiff(Image3D& a_o, const Image3D& a_i, 
                DimT dim,  DiffT diffType=DIFF_CENTRAL, 
                enum BoundaryCondT bc=BC_CLAMP, 
                bool accum=false, OpT op=OP_VAL,
                StreamT stream=NULL);

/**
* Compute the gradient of an image
*  d_o = grad (a_i)
*/
void Gradient(Field3D& a_o, const float* a_i, 
              DiffT diffType = DIFF_CENTRAL, 
              BoundaryCondT bc=BC_CLAMP, 
              StreamT s = NULL);
void Gradient(Field3D& a_o, const Image3D& a_i, 
              DiffT diffType = DIFF_CENTRAL, 
              BoundaryCondT bc=BC_CLAMP, 
              StreamT s = NULL);
void Gradient2(Field3D& a_o, const Image3D& a_i, 
               DiffT diffType = DIFF_CENTRAL, 
               BoundaryCondT bc=BC_CLAMP, 
               StreamT s = NULL);

/**
 * Masked gradient computation
 */
   
void GradientMask(Field3D& a_o, 
                  const Image3D& a_i, 
                  const Image3D& a_mask,
                  DiffT diffType=DIFF_CENTRAL, 
                  BoundaryCondT bc=BC_CLAMP, 
                  StreamT s=NULL);
/**
 * Compute the gradient of an image using forward differences
 *  d_o = grad (a_i)
 */
void GradFor(Field3D& a_o, const Image3D& a_i, 
             BoundaryCondT bc=BC_CLAMP, StreamT s = NULL);
void GradForMag(Image3D& a_o, const Image3D& a_i, StreamT s = NULL);
void GradientMag(Image3D& a_o, const Image3D& a_i, 
                 DiffT diffType = DIFF_CENTRAL, 
                 BoundaryCondT bc=BC_CLAMP, 
                 StreamT s = NULL);

/**
* returns deriv(h_i), where deriv is upwind derivative based on
* a_speed in dimension 'dim'
*/
void UpwindDiff(Image3D& a_o, const Image3D& a_i, 
                const Image3D& a_speed,
                DimT dim,
                StreamT s = NULL);

/**
  returns |\nabla a_i|*a_speed, where \nabla taken with upwind
  differences based on a_speed
*/
void UpwindGradMag(Image3D& a_o, const Image3D& a_i,
                   const Image3D& a_speed, StreamT s=NULL);
    
/**
 * Compute the diverence of a field
 *  d_o = div (a_i)
 */
void Divergence(Image3D& a_o, const Field3D& a_i, 
                DiffT diffType = DIFF_CENTRAL, 
                BoundaryCondT bc=BC_CLAMP, 
                StreamT s = NULL);

/**
 * Compute the diverence of a field using backward differences
 *  d_o = div (a_i)
 */
void DivBack(Image3D& a_o, const Field3D& a_i,
             BoundaryCondT bc=BC_CLAMP, StreamT s = NULL);

/**
 * Compute the magnitude image
 * a_o[i] = sqrt(a_i[i].x^2 + a_i[i].y^2 + a_i[i].z^2) 
 */
void Magnitude(Image3D& a_o, const Field3D& a_i, StreamT s = NULL);

/**
@brief
   Compute the squared magnitude image from a vector field

   a_o[i] = a_i[i].x^2 + a_i[i].y^2 + a_i[i].z^2 
*/
void SqrMagnitude(Image3D& a_o, const Field3D& a_i, StreamT s = NULL);

/**
@brief
   Compute the dot product between two fields 

   d_o[i] = d_i[i].x * d_i1[i].x + d_i[i].y * d_i1[i].y + d_i[i].z * d_i1[i].z
*/
void ComponentDotProd(Image3D& d_o, const Field3D& d_i, const Field3D& d_i1, StreamT s = NULL);
    
/** @brief a_o.x = a_i.x + a_i1, a_o.y = a_i.y + a_i1, a_o.z = a_i.z + a_i1 */
void Add(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, StreamT s = NULL);

/** @brief a_o.x = a_i.x - a_i1, a_o.y = a_i.y - a_i1, a_o.z = a_i.z - a_i1 */
void Sub(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, StreamT s = NULL);

/** @brief a_o.x = a_i.x * a_i1, a_o.y = a_i.y * a_i1, a_o.z = a_i.z * a_i1 */
void Mul(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, StreamT s = NULL);

/** @brief a_o.x = a_i.x / a_i1, a_o.y = a_i.y / a_i1, a_o.z = a_i.z / a_i1 */
void Div(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, StreamT s = NULL);

/** @brief a_o.x += a_i, a_o.y += a_i, a_o.y += a_i, */
void Add_I(Field3D& a_o, const Image3D& a_i, StreamT s = NULL);

/** @brief a_o.x -= a_i, a_o.y -= a_i, a_o.y -= a_i, */
void Sub_I(Field3D& a_o, const Image3D& a_i, StreamT s = NULL);

/** @brief a_o.x *= a_i, a_o.y *= a_i, a_o.y *= a_i, */
void Mul_I(Field3D& a_o, const Image3D& a_i, StreamT s = NULL);

/** @brief a_o.x /= a_i, a_o.y = a_i.y / a_i, a_o.z = a_i.z / a_i */
void Div_I(Field3D& a_o, const Image3D& a_i, StreamT s = NULL);

/** @brief a_o = a_i + a_i1 * a_i2 */
void Add_Mul(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const Image3D& a_i2, StreamT s = NULL);

/** @brief a_o = a_o + a_i * a_i1 */
void Add_Mul_I(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, StreamT s = NULL);

/** @brief a_o = a_i - a_i1 * a_i2 */
void Sub_Mul(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const Image3D& a_i2, StreamT s = NULL);

/** @brief a_o = a_o - a_i * a_i1 */
void Sub_Mul_I(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, StreamT s = NULL);

/** @brief a_o = a_i * a_i1 * c (a_o.x = a_i.x * a_i1 * c)*/
void MulMulC(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, const float& c, StreamT s = NULL, bool onDev = false);

/** @brief a_o = a_o * a_i * c  (a_o.x = a_o.x * a_i * c)*/
void MulMulC_I(Field3D& a_o, const Image3D& a_i, const float& c, StreamT s = NULL, bool onDev = false);

/** @brief a_o = a_i + a_i1 * a_i2 * c */
void Add_MulMulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const Image3D& a_i2, const float& c,  StreamT s = NULL, bool onDev = false);

// // unimplemented
// void Mul_Add_Mul(Field3D& a_fo, const Field3D& a_fi, const Image3D& a_ii,
//                         const Field3D& a_fi1, const Image3D& a_ii1, StreamT s = NULLtream);

/** @brief a_o = a_o + a_i * a_i1 * c */
void Add_MulMulC_I(const Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, const float& c,  StreamT s = NULL, bool onDev = false);
    
void JacDetH(Image3D& a_jdet, const Field3D& a_h, 
             DiffT diffType=DIFF_CENTRAL, 
             BoundaryCondT bc=BC_APPROX, 
             StreamT s=NULL);

void JacDetH(Image3D& a_detJ,const Field3D& a_Xg, const Field3D& a_Yg, const Field3D& a_Zg, StreamT s = NULL);

void JacDetV(Image3D& a_detJ,const Field3D& a_Xg, const Field3D& a_Yg, const Field3D& a_Zg, StreamT s = NULL);

void JacDetVInv(Image3D& a_detJ,const Field3D& a_Xg,const Field3D& a_Yg,const Field3D& a_Zg,StreamT s = NULL);

} // end namespace Opers

} // end namespace PyCA

#endif

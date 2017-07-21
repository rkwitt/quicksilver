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

#ifndef __F_OPERS_H
#define __F_OPERS_H

#ifndef SWIG
#include <pycaConst.h>
#include <estream.h>
#include <Vec3D.h>
#include <MemOpers.h>
#endif // SWIG

#include <FieldOpers.h>

namespace PyCA {

class Image3D;
class Field3D;

bool Is1D(const Field3D& a_i, size_t n);
bool Is1D(const Field3D& a_i, const Field3D& a_i1, size_t n);
bool Is1D(const Field3D& a_i, const Field3D& a_i1, const Field3D& a_i2, size_t n);
bool Is1D(const Field3D& a_i, const Field3D& a_i1, const Field3D& a_i2, const Field3D& a_i3, size_t n);


/*In H-V computation 
 *  - h is always defined in the normal coordinate 
 *  - v is in the regular space coordinate 
 *
 *  H-V computation is always performed using following conversion
 *
 *  - img(h = I) = img;
 *  - v = (h - I) * sp
 *  - h = v / sp = v * iSp + I
 *
 * The diplacement (u) is compute by the velocity in the time unit
 * so we also have to take into account space for u computation
 * u = v * delta
 */

/*
  Admissible trilinear interpolation based on the bondary condition of the field
  - HField: BACKGROUND_STRATEGY_ID, BACKGROUND_STRATEGY_PARTIAL_ID
  - VField: BACKGROUND_STRATEGY_ZERO, BACKGROUND_STRATEGY_PARTIAL_ZERO
  - Image : BACKGROUND_STRATEGY_ZERO, BACKGROUND_STRATEGY_PARTIAL_ZERO,
  BACKGROUND_STRATEGY_CLAMP, BACKGROUND_STRATEGY_WARP
*/
namespace Opers{

/**
 * Copy a_i to a_o
 */
    void Copy(Field3D& a_o, const Field3D& a_i, StreamT s = NULL);
/**
 * Set all the value of the Vector field with single value Vector3f(v)
 */
    void SetMem(Field3D& a_o, const Vec3Df& v, StreamT s = NULL,bool onDev=false);
    void SetMem(Field3D& a_o, const Vec3Df& v, Image3D& mask, 
		StreamT s = NULL,bool onDev=false);
/**  
 * Set all the value of the Vector field with single float value
 * (normaly used to zero out the memory)
 */
    void SetMem(Field3D& a_o, const float& c, StreamT s = NULL,bool onDev=false);
    void SetMem(Field3D& a_o, const float& c, Image3D& mask, 
		StreamT s = NULL,bool onDev=false);
    
/** @defgroup operator
 * Basic operator on the Vector field
 * @{
 */


//
// Add
//
/** @brief a_o = a_i + c */
FIELD_BINARY_ARRAY_OPC_VEC_DEC(Add)
/** @brief a_o = a_i + c */
FIELD_BINARY_ARRAY_OPC_FLOAT_DEC(Add)
/** @brief a_o = a_i + c */
FIELD_BINARY_ARRAY_OPC_I_VEC_DEC(Add)
/** @brief a_o = a_i + c */
FIELD_BINARY_ARRAY_OPC_I_FLOAT_DEC(Add)
/** @brief a_o = a_i + a_i1 */
FIELD_BINARY_ARRAY_OP_FIELD_DEC(Add)
/** @brief a_o += a_i1 */
FIELD_BINARY_ARRAY_OP_I_FIELD_DEC(Add)

//
// Sub
//
/** @brief a_o = a_i - c */
FIELD_BINARY_ARRAY_OPC_VEC_DEC(Sub)
/** @brief a_o = a_i - c */
FIELD_BINARY_ARRAY_OPC_FLOAT_DEC(Sub)
/** @brief a_o = a_i - c */
FIELD_BINARY_ARRAY_OPC_I_VEC_DEC(Sub)
/** @brief a_o = a_i - c */
FIELD_BINARY_ARRAY_OPC_I_FLOAT_DEC(Sub)
/** @brief a_o = a_i - a_i1 */
FIELD_BINARY_ARRAY_OP_FIELD_DEC(Sub)
/** @brief a_o -= a_i1 */
FIELD_BINARY_ARRAY_OP_I_FIELD_DEC(Sub)

//
// Mul
//
/** @brief a_o = a_i * c */
FIELD_BINARY_ARRAY_OPC_VEC_DEC(Mul)
/** @brief a_o = a_i * c */
FIELD_BINARY_ARRAY_OPC_FLOAT_DEC(Mul)
/** @brief a_o = a_i * c */
FIELD_BINARY_ARRAY_OPC_I_VEC_DEC(Mul)
/** @brief a_o = a_i * c */
FIELD_BINARY_ARRAY_OPC_I_FLOAT_DEC(Mul)
/** @brief a_o = a_i * a_i1 */
FIELD_BINARY_ARRAY_OP_FIELD_DEC(Mul)
/** @brief a_o *= a_i1 */
FIELD_BINARY_ARRAY_OP_I_FIELD_DEC(Mul)

//
// Div
//
/** @brief a_o = a_i / c */
FIELD_BINARY_ARRAY_OPC_VEC_DEC(Div)
/** @brief a_o = a_i / c */
FIELD_BINARY_ARRAY_OPC_FLOAT_DEC(Div)
/** @brief a_o = a_i / c */
FIELD_BINARY_ARRAY_OPC_I_VEC_DEC(Div)
/** @brief a_o = a_i / c */
FIELD_BINARY_ARRAY_OPC_I_FLOAT_DEC(Div)
/** @brief a_o = a_i / a_i1 */
FIELD_BINARY_ARRAY_OP_FIELD_DEC(Div)
/** @brief a_o /= a_i1 */
FIELD_BINARY_ARRAY_OP_I_FIELD_DEC(Div)

//
// Max
//
/** @brief a_o = [max(a_i.x, c.x),max(a_i.y, c.y),max(a_i.z, c.z)] */
FIELD_BINARY_ARRAY_OPC_VEC_DEC(Max)
/** @brief a_o = [max(a_i.x, c),max(a_i.y, c),max(a_i.z, c)] */
FIELD_BINARY_ARRAY_OPC_FLOAT_DEC(Max)
/** @brief a_o = [max(a_o.x, c.x),max(a_o.y, c.y),max(a_o.z, c.z)] */
FIELD_BINARY_ARRAY_OPC_I_VEC_DEC(Max)
/** @brief a_o = [max(a_o.x, c),max(a_o.y, c),max(a_o.z, c)] */
FIELD_BINARY_ARRAY_OPC_I_FLOAT_DEC(Max)
/** @brief a_o = [max(a_i.x, a_i1.x),max(a_i.y, a_i1.y),max(a_i.z, a_i1.z)] */
FIELD_BINARY_ARRAY_OP_FIELD_DEC(Max)
/** @brief a_o = [max(a_o.x, a_i.x),max(a_o.y, a_i.y),max(a_o.z, a_i.z)] */
FIELD_BINARY_ARRAY_OP_I_FIELD_DEC(Max)

//
// Min
//
/** @brief a_o = [min(a_i.x, c.x),min(a_i.y, c.y),min(a_i.z, c.z)] */
FIELD_BINARY_ARRAY_OPC_VEC_DEC(Min)
/** @brief a_o = [min(a_i.x, c),min(a_i.y, c),min(a_i.z, c)] */
FIELD_BINARY_ARRAY_OPC_FLOAT_DEC(Min)
/** @brief a_o = [min(a_o.x, c.x),min(a_o.y, c.y),min(a_o.z, c.z)] */
FIELD_BINARY_ARRAY_OPC_I_VEC_DEC(Min)
/** @brief a_o = [min(a_o.x, c),min(a_o.y, c),min(a_o.z, c)] */
FIELD_BINARY_ARRAY_OPC_I_FLOAT_DEC(Min)
/** @brief a_o = [min(a_i.x, a_i1.x),min(a_i.y, a_i1.y),min(a_i.z, a_i1.z)] */
FIELD_BINARY_ARRAY_OP_FIELD_DEC(Min)
/** @brief a_o = [min(a_o.x, a_i.x),min(a_o.y, a_i.y),min(a_o.z, a_i.z)] */
FIELD_BINARY_ARRAY_OP_I_FIELD_DEC(Min)

/** @brief a_o = (a_i + a_i1) * c */
    void AddMulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const float& c, StreamT s = NULL,bool onDev=false);
/** @brief a_o = (a_i - a_i1) * c */
    void SubMulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const float& c, StreamT s = NULL,bool onDev=false);
/** @brief a_o = (a_i * a_i1) * c */
    void MulMulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const float& c, StreamT s = NULL,bool onDev=false);

    
/** @brief a_o = (a_o + a_i) * c */
    void AddMulC_I(Field3D& a_o, const Field3D& a_i, const float& c, StreamT s = NULL,bool onDev=false);

/** @brief a_o = (a_o - a_i) * c */
    void SubMulC_I(Field3D& a_o, const Field3D& a_i, const float& c, StreamT s = NULL,bool onDev=false);

/** @brief a_o = (a_o * a_i) * c */
    void MulMulC_I(Field3D& a_o, const Field3D& a_i, const float& c, StreamT s = NULL,bool onDev=false);


/** @brief a_o = a_i + a_i1 * c */
    void Add_MulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const float& c, StreamT s = NULL,bool onDev=false);
/** @brief a_o = a_o+ a_i1 * c */
    void Add_MulC_I(Field3D& a_o, const Field3D& a_i1, const float& c, StreamT s = NULL,bool onDev=false);

/** @brief a_o = a_i * c + a_i1 */
    void MulCAdd(Field3D& a_o, const Field3D& a_i, const float& c, const Field3D& a_i1, StreamT s = NULL,bool onDev=false);
                 
/** @brief a_o = a_i * c - a_i1 */
    void MulCSub(Field3D& a_o, const Field3D& a_i, const float& c, const Field3D& a_i1, StreamT s = NULL,bool onDev=false);

/** @brief a_o = a_o * c + a_i */
    void MulCAdd_I(Field3D& a_o, const float& c, const Field3D& a_i, StreamT s = NULL,bool onDev=false);
/** @brief a_o = a_o * c - a_i */
    void MulCSub_I(Field3D& a_o, const float& c, const Field3D& a_i, StreamT s = NULL,bool onDev=false);
/** @brief a_o = (a_i + a) * b */
    void AddCMulC(Field3D& a_o, const Field3D& a_i, const Vec3Df& a, const Vec3Df& b, StreamT s = NULL,bool onDev=false);
/** @brief a_o = (a_i * a) + b */
    void MulCAddC(Field3D& a_o, const Field3D& a_i, const Vec3Df& a, const Vec3Df& b, StreamT s = NULL,bool onDev=false);
/** @brief a_o = (a_o + a) * b */
    void AddCMulC_I(Field3D& a_o, const Vec3Df& a, const Vec3Df& b, StreamT s = NULL,bool onDev=false);
/** @brief a_o = (a_o * a) + b */
    void MulCAddC_I(Field3D& a_o, const Vec3Df& a, const Vec3Df& b, StreamT s = NULL,bool onDev=false);

/** @brief a_o = a_i * ca + a_i1 * cb */
    void MulC_Add_MulC(Field3D& a_o,
                              const Field3D& a_i, const float& ca,
                              const Field3D& a_i1, const float& cb, StreamT s = NULL,bool onDev=false);
/** @brief a_o = a_o * co + a_i1 * cb */
    void MulC_Add_MulC_I(Field3D& a_o, const float& co,
                         const Field3D& a_i1, const float& cb, StreamT s = NULL,bool onDev=false);
/** @brief a_o = (a_i + a) * b + c */
    void AddCMulCAddC(Field3D& a_o, const Field3D& a_i,
                      const float& a, const float& b, const float& c, StreamT s = NULL,bool onDev=false);
/** @brief a_o = (a_o + a) * b + c */
    void AddCMulCAddC_I(Field3D& a_o, const float& a, const float& b, const float& c, StreamT s = NULL,bool onDev=false);

/** @brief a_o = (a_i + a) * b + c */
    void AddCMulCAddC(Field3D& a_o, const Field3D& a_i, const Vec3Df& a, const Vec3Df& b, const Vec3Df& c, StreamT s = NULL,bool onDev=false);
    
/** @brief a_o = (a_o + a) * b + c */
    void AddCMulCAddC_I(Field3D& a_o, const Vec3Df& a, const Vec3Df& b, const Vec3Df& c, StreamT s = NULL,bool onDev=false);
    
/** @brief a_o = a_i + a_i1 * a_i2 * c */
    void Add_MulMulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1,
                     const Field3D& a_i2, const float& c, StreamT s = NULL,bool onDev=false);
/** @brief a_o = a_o + a_i * a_i1 * c */
    void Add_MulMulC_I(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const float& c, StreamT s = NULL,bool onDev=false);

/** @brief a_o = a_i - a_i1 * a_i2 * c */
    void Sub_MulMulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1,
                     const Field3D& a_i2, const float& c, StreamT s = NULL,bool onDev=false);
/** @brief a_o = a_o - a_i * a_i1 * c */
    void Sub_MulMulC_I(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const float& c, StreamT s = NULL,bool onDev=false);

/*
 * Advance function
 */ 

/**
compose two hfields to get an hfield

f(x) = g(h(x))
*/
    template<BackgroundStrategy bg>
    void ComposeHH(Field3D& a_f, const Field3D& a_g,
		   const Field3D& a_h, StreamT s = NULL);

    // non-template version
    void ComposeHH(Field3D& a_f, const Field3D& a_g,
	       const Field3D& a_h,
	       BackgroundStrategy bg = BACKGROUND_STRATEGY_PARTIAL_ID,
	       StreamT s = NULL);
/**
@brief
   compose a velocity and h field to get an hfield

   h(x) = g(x) + delta * v(g(x))
*/
    template<BackgroundStrategy bg>
    void ComposeVH(Field3D& a_h, const Field3D& a_v,
                          const Field3D& a_g, const float& delta=1.f, StreamT s = NULL,bool onDev=false);

   // non-template version
    void ComposeVH(Field3D& a_h, const Field3D& a_v,
                          const Field3D& a_g, const float& delta=1.f, 
			  BackgroundStrategy bg = BACKGROUND_STRATEGY_PARTIAL_ZERO,
			  StreamT s = NULL,bool onDev=false);

/**
@brief
   compose a inverse velocity and h field to get an hfield

   h(x) = g(x) - delta * v(g(x))
*/
template<BackgroundStrategy bg>
void ComposeVInvH(Field3D& a_h, const Field3D& a_v, const Field3D& a_g,
                  const float& delta=1.f, StreamT s = NULL,bool onDev=false);

// non-template version
void ComposeVInvH(Field3D& a_h, const Field3D& a_v, const Field3D& a_g, const float& delta=1.f, 
                  BackgroundStrategy bg = BACKGROUND_STRATEGY_PARTIAL_ZERO,
                  StreamT s = NULL,bool onDev=false);

/**
@brief
   compose a h field and a velocify field to get an hfield

   h(x) = g(x + delta * v(x))
*/
template<BackgroundStrategy bg>
void ComposeHV(Field3D& a_h, const Field3D& a_g, const Field3D& a_v,
               const float& delta=1.f, StreamT s = NULL,bool onDev=false);

// non-template version
void ComposeHV(Field3D& a_h, const Field3D& a_g, 
               const Field3D& a_v, const float& delta=1.f, 
               BackgroundStrategy bg = BACKGROUND_STRATEGY_PARTIAL_ID,
               StreamT s = NULL,bool onDev=false);

/**
@brief
   compose a h field and an inverse velocity field to get an hfield

   h(x) = g(x - delta * v(x))
*/
template<BackgroundStrategy bg>
void ComposeHVInv(Field3D& a_h, const Field3D& a_g, const Field3D& a_v,
                  const float& delta=1.f, StreamT s = NULL,bool onDev=false);

// non-template version
void ComposeHVInv(Field3D& a_h, const Field3D& a_g, const Field3D& a_v,const float& delta=1.f, 
                  BackgroundStrategy bg = BACKGROUND_STRATEGY_PARTIAL_ID,
                  StreamT s = NULL,bool onDev=false);

/**
@brief
   precompose h field with translation

   h(x) = f(x + t)
*/
template<BackgroundStrategy bg>
void ComposeTranslation(Field3D& a_h, const Field3D& a_f, const Vec3Df& t, StreamT s = NULL,bool onDev=false);

// non-template version
void ComposeTranslation(Field3D& a_h, const Field3D& a_f, const Vec3Df& t, 
                        BackgroundStrategy bg = BACKGROUND_STRATEGY_CLAMP,
                        StreamT s = NULL,bool onDev=false);

void Jacobian(Field3D& a_Xg, Field3D& a_Yg, Field3D& a_Zg, 
              const Field3D& a_h, 
              DiffT diffType = DIFF_CENTRAL,
              BoundaryCondT bc = BC_APPROX, 
              StreamT s = NULL);
    
/**
deform a field using hfield

a_o(x) = a_i(h(x))
*/
template<BackgroundStrategy bg>
void ApplyH(Field3D& a_o, const Field3D& a_i, const Field3D& h, StreamT s = NULL);

// non-template version
void ApplyH(Field3D& a_o, const Field3D& a_i, const Field3D& h, 	
            BackgroundStrategy bg, StreamT s = NULL);

/**
deform a field using displacement field

a_o(x) = a_i(x + delta * u)
*/
template<BackgroundStrategy bg>
void ApplyV(Field3D& a_o, const Field3D& a_i, const Field3D& u,
            const float& delta=1.f, StreamT s = NULL,bool onDev=false);

// non-template version
void ApplyV(Field3D& a_o, const Field3D& a_i, 
            const Field3D& u, 
            BackgroundStrategy bg, 
            const float& delta=1.f, 
            StreamT s = NULL,bool onDev=false);

/**
deform a field using displacement field

a_o(x) = a_i(x - delta * u)
*/
template<BackgroundStrategy bg>
void ApplyVInv(Field3D& a_o, const Field3D& a_i, const Field3D& u,
               const float& delta=1.f, StreamT s = NULL,bool onDev=false);

// non-template version
void ApplyVInv(Field3D& a_o, const Field3D& a_i, const 
               Field3D& u, 
               BackgroundStrategy bg,
               const float& delta=1.f, 
               StreamT s = NULL,bool onDev=false);

void SplatField(Field3D& a_o, const Field3D& a_i, 
                const Field3D& h, StreamT stream=NULL);

void SubVol(Field3D& a_o,const Field3D& a_i,
            const Vec3Di& start, StreamT st=NULL);

void SetSubVol_I(Field3D& a_o,const Field3D& a_i,
                 const Vec3Di& start, StreamT st=NULL);

template<BackgroundStrategy bg,  bool rescaleVector>
void Resample(Field3D& a_o, const Field3D& a_i, StreamT s = NULL);

// non-template version
void Resample(Field3D& a_o, const Field3D& a_i, 
              BackgroundStrategy bg = BACKGROUND_STRATEGY_CLAMP,
              bool rescaleVector = false, StreamT s = NULL);

void ReprojectToUnitVec(Field3D& a_o, StreamT st=NULL);

// normalize all vectors to have unit length, except those less
// than length eps which will be set to zero
void NormalizeSafe(Field3D& a_o, const Field3D& a_i, const float& eps, StreamT st=NULL);
void NormalizeSafe_I(Field3D& a_o, const float& eps, StreamT st=NULL);

void Shrink(Field3D& a_o, const Field3D& a_i, const float& eps, StreamT st=NULL);
void Shrink_I(Field3D& a_o, const float& eps, StreamT st=NULL);

// invert a deformation field g, given displacement and previous estimate of inverse
template<BackgroundStrategy bg>
void FixedPointInverse(Field3D &ginv, const Field3D& g, unsigned int numIter=2,
                       StreamT stream = NULL, bool onDev=false);

// non-template version
void FixedPointInverse(Field3D &ginv, const Field3D& g, unsigned int numIter=2,
                       StreamT stream = NULL, bool onDev=false);

// update inverse according to time stepping
void UpdateInverse(Field3D &ginv0t1,Field3D &scratchV, const Field3D& ginv0t, const Field3D& w, unsigned int numIter = 3, StreamT stream = NULL, bool onDev = false);

/**
Adjoint action of diffeomorphism g on momentum vector field X, which is just
the pushforward of the vector field X:

Z = Ad_g X = |Dg|\circ g^{-1} X\circ g^{-1}
*/
template<BackgroundStrategy bg>
void Ad(Field3D& Z, const Field3D& g, const Field3D& X,
        StreamT stream = NULL, bool onDev=false);
// non-template version
void Ad(Field3D& Z, const Field3D& g, const Field3D& X,
        StreamT stream = NULL, BackgroundStrategy bg = BACKGROUND_STRATEGY_PARTIAL_ZERO,
        bool onDev=false);
/**
 * Infinitesimal adjoint action of vector field X on vector field Y,
 * which is the negative Jacobi-Lie bracket
 * Z = ad_X Y = DX Y - DY X
 */
void AdInf(Field3D& Z, const Field3D& X, const Field3D& Y, StreamT stream = NULL, bool onDev=false);

/**
Coadjoint action of diffeomorphism g on momentum vector field m
This is the dual of the adjoint action Ad_g, given by the formula:

n(x) = Ad_g^* m(x) = (Dg(x))^T m(g(x)) |Dg(x)|
 */
template<BackgroundStrategy bg>
void CoAd(Field3D& n, const Field3D& g, const Field3D& m, StreamT stream = NULL, bool onDev=false);
// non-template version
void CoAd(Field3D& n, const Field3D& g, const Field3D& m, StreamT stream = NULL,
          BackgroundStrategy bg = BACKGROUND_STRATEGY_PARTIAL_ZERO, bool onDev=false);
/**
 * Infinitesimal coadjoint action of vector field X on momentum vector field m,
 which is the negative Jacobi-Lie bracket

 n = ad_X^* m = (DX)^T m + div(m \otimes X)
*/
    void CoAdInf(Field3D& n, const Field3D& X, const Field3D& m, StreamT stream = NULL, bool onDev=false);

/**
 * Computes tensor divergence of outer product of two vector fields,
 which is also the second term in coadjoint action of vector fields as
 in CoAdInf.

 Z = div(X \otimes Y)
*/
    void DivergenceTensor(Field3D& Z, const Field3D& X, const Field3D& Y, StreamT stream = NULL, bool onDev=false);
    void JacobianXY(Field3D& Z, const Field3D& X, const Field3D& Y, StreamT st=NULL, bool onDev=false);
    void JacobianXtY(Field3D& Z, const Field3D& X, const Field3D& Y, StreamT st=NULL, bool onDev=false);

} // end namespace Opers

} // end namespace PyCA

#endif

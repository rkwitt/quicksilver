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

#ifndef __FIELD3D_OPERS_H
#define __FIELD3D_OPERS_H

#ifndef SWIG
#include <pycaConst.h>
#include <estream.h>
#include <Vec3D.h>
#include <MemOpers.h>
#endif // SWIG

#define FIELD_BINARY_ARRAY_OPC_VEC_DEC(OP)     	       	       	       	   \
void OP##C(Field3D& a_o, const Field3D& a_i, const Vec3Df& c,		   \
	   StreamT st=NULL,bool onDev=false);

#define FIELD_BINARY_ARRAY_OPC_FLOAT_DEC(OP)				   \
void OP##C(Field3D& a_o, const Field3D& a_i, const float& c,		   \
	   StreamT st=NULL,bool onDev=false);

#define FIELD_BINARY_ARRAY_OPC_I_VEC_DEC(OP)				   \
void OP##C_I(Field3D& a_o, const Vec3Df& c,				   \
	     StreamT st=NULL, bool onDev=false);

#define FIELD_BINARY_ARRAY_OPC_I_FLOAT_DEC(OP)				   \
void OP##C_I(Field3D& a_o, const float& c,				   \
	     StreamT st=NULL,bool onDev=false);

#define FIELD_BINARY_ARRAY_OP_FIELD_DEC(OP)				   \
void OP(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1,		   \
	StreamT st = NULL);

#define FIELD_BINARY_ARRAY_OP_I_FIELD_DEC(OP)				   \
void OP##_I(Field3D& a_o, const Field3D& a_i1,				   \
            StreamT st = NULL);

namespace PyCA {

class Image3D;
class Field3D;

bool Is1D(const Field3D& d_i, size_t n);
bool Is1D(const Field3D& d_i, const Field3D& d_i1, size_t n);
bool Is1D(const Field3D& d_i, const Field3D& d_i1, const Field3D& d_i2, size_t n);
bool Is1D(const Field3D& d_i, const Field3D& d_i1, const Field3D& d_i2, const Field3D& d_i3, size_t n);


/* In H-V computation 
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
 *  Admissible trilinear interpolation based on the bondary condition of the field
 *  - HField: BACKGROUND_STRATEGY_ID, BACKGROUND_STRATEGY_PARTIAL_ID
 *  - VField: BACKGROUND_STRATEGY_ZERO, BACKGROUND_STRATEGY_PARTIAL_ZERO
 *  - Image : BACKGROUND_STRATEGY_ZERO, BACKGROUND_STRATEGY_PARTIAL_ZERO,
 *  BACKGROUND_STRATEGY_CLAMP, BACKGROUND_STRATEGY_WRAP
 */
class CFieldOpers;
class GFieldOpers;

template<int mode>
class FieldOpers{
public:
  typedef typename BinSelector<mode, CFieldOpers, GFieldOpers>::Result Executer;
  enum { exec_mode = mode };

/** Set all the value of the Vector field with single Vector3f(v) */
  static void SetMem(Field3D& a_o, const Vec3Df& v, StreamT st,bool onDev=false);
  static void SetMem(Field3D& a_o, const Vec3Df& v, const Image3D &mask,
		     StreamT st,bool onDev=false);
/**  
 * Set all the value of the Vector field with single float value
 * (normally used to zero out the memory)
 */
  static void SetMem(Field3D& a_o, const float& c, StreamT st,bool onDev=false);
  static void SetMem(Field3D& a_o, const float& c, const Image3D &mask,
		     StreamT st,bool onDev=false);
  
/** @defgroup operator
 * Basic operator on the Vector field
 * @{
 */
   //
   // Add
   //
   /** @brief a_o = a_i + c */
   static FIELD_BINARY_ARRAY_OPC_VEC_DEC(Add)
   /** @brief a_o = a_i + c */
   static FIELD_BINARY_ARRAY_OPC_FLOAT_DEC(Add)
   /** @brief a_o = a_o + c */
   static FIELD_BINARY_ARRAY_OPC_I_VEC_DEC(Add)
   /** @brief a_o = a_o + c */
   static FIELD_BINARY_ARRAY_OPC_I_FLOAT_DEC(Add)
   /** @brief a_o = a_o + a_i1 */
   static FIELD_BINARY_ARRAY_OP_FIELD_DEC(Add)
   /** @brief a_o += a_i1 */
   static FIELD_BINARY_ARRAY_OP_I_FIELD_DEC(Add)

   //
   // Sub
   //
   /** @brief a_o = a_i - c */
   static FIELD_BINARY_ARRAY_OPC_VEC_DEC(Sub)
   /** @brief a_o = a_i - c */
   static FIELD_BINARY_ARRAY_OPC_FLOAT_DEC(Sub)
   /** @brief a_o = a_o - c */
   static FIELD_BINARY_ARRAY_OPC_I_VEC_DEC(Sub)
   /** @brief a_o = a_o - c */
   static FIELD_BINARY_ARRAY_OPC_I_FLOAT_DEC(Sub)
   /** @brief a_o = a_i - a_i1 */
   static FIELD_BINARY_ARRAY_OP_FIELD_DEC(Sub)
   /** @brief a_o -= a_i1 */
   static FIELD_BINARY_ARRAY_OP_I_FIELD_DEC(Sub)

   //
   // Mul
   //
   /** @brief a_o = a_i * c */
   static FIELD_BINARY_ARRAY_OPC_VEC_DEC(Mul)
   /** @brief a_o = a_i * c */
   static FIELD_BINARY_ARRAY_OPC_FLOAT_DEC(Mul)
   /** @brief a_o = a_o * c */
   static FIELD_BINARY_ARRAY_OPC_I_VEC_DEC(Mul)
   /** @brief a_o = a_o * c */
   static FIELD_BINARY_ARRAY_OPC_I_FLOAT_DEC(Mul)
   /** @brief a_o = a_i * a_i1 */
   static FIELD_BINARY_ARRAY_OP_FIELD_DEC(Mul)
   /** @brief a_o *= a_i1 */
   static FIELD_BINARY_ARRAY_OP_I_FIELD_DEC(Mul)

   //
   // Div
   //
   /** @brief a_o = a_i / c */
   static FIELD_BINARY_ARRAY_OPC_VEC_DEC(Div)
   /** @brief a_o = a_i / c */
   static FIELD_BINARY_ARRAY_OPC_FLOAT_DEC(Div)
   /** @brief a_o = a_o / c */
   static FIELD_BINARY_ARRAY_OPC_I_VEC_DEC(Div)
   /** @brief a_o = a_o / c */
   static FIELD_BINARY_ARRAY_OPC_I_FLOAT_DEC(Div)
   /** @brief a_o = a_i / a_i1 */
   static FIELD_BINARY_ARRAY_OP_FIELD_DEC(Div)
   /** @brief a_o /= a_i1 */
   static FIELD_BINARY_ARRAY_OP_I_FIELD_DEC(Div)

   //
   // Max
   //
   /** @brief a_o = [max(a_i.x, c.x),max(a_i.y, c.y),max(a_i.z, c.z)] */
   static FIELD_BINARY_ARRAY_OPC_VEC_DEC(Max)
   /** @brief a_o = [max(a_i.x, c),max(a_i.y, c),max(a_i.z, c)] */
   static FIELD_BINARY_ARRAY_OPC_FLOAT_DEC(Max)
   /** @brief a_o = [max(a_o.x, c.x),max(a_o.y, c.y),max(a_o.z, c.z)] */
   static FIELD_BINARY_ARRAY_OPC_I_VEC_DEC(Max)
   /** @brief a_o = [max(a_o.x, c),max(a_o.y, c),max(a_o.z, c)] */
   static FIELD_BINARY_ARRAY_OPC_I_FLOAT_DEC(Max)
   /** @brief a_o = [max(a_i.x, a_i1.x),max(a_i.y, a_i1.y),max(a_i.z, a_i1.z)] */
   static FIELD_BINARY_ARRAY_OP_FIELD_DEC(Max)
   /** @brief a_o = [max(a_o.x, a_i.x),max(a_o.y, a_i.y),max(a_o.z, a_i.z)] */
   static FIELD_BINARY_ARRAY_OP_I_FIELD_DEC(Max)

   //
   // Min
   //
   /** @brief a_o = [min(a_i.x, c.x),min(a_i.y, c.y),min(a_i.z, c.z)] */
   static FIELD_BINARY_ARRAY_OPC_VEC_DEC(Min)
   /** @brief a_o = [min(a_i.x, c),min(a_i.y, c),min(a_i.z, c)] */
   static FIELD_BINARY_ARRAY_OPC_FLOAT_DEC(Min)
   /** @brief a_o = [min(a_o.x, c.x),min(a_o.y, c.y),min(a_o.z, c.z)] */
   static FIELD_BINARY_ARRAY_OPC_I_VEC_DEC(Min)
   /** @brief a_o = [min(a_o.x, c),min(a_o.y, c),min(a_o.z, c)] */
   static FIELD_BINARY_ARRAY_OPC_I_FLOAT_DEC(Min)
   /** @brief a_o = [min(a_i.x, a_i1.x),min(a_i.y, a_i1.y),min(a_i.z, a_i1.z)] */
   static FIELD_BINARY_ARRAY_OP_FIELD_DEC(Min)
   /** @brief a_o = [min(a_o.x, a_i.x),min(a_o.y, a_i.y),min(a_o.z, a_i.z)] */
   static FIELD_BINARY_ARRAY_OP_I_FIELD_DEC(Min)

/** @brief a_o = (a_i + a_i1) * c */
  static void AddMulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const float& c, StreamT st,bool onDev=false);
/** @brief a_o = (a_i - a_i1) * c */
  static void SubMulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const float& c, StreamT st,bool onDev=false);
/** @brief a_o = (a_i * a_i1) * c */
  static void MulMulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const float& c, StreamT st,bool onDev=false);
    
/** @brief a_o = (a_o + a_i) * c */
  static void AddMulC_I(Field3D& a_o, const Field3D& a_i, const float& c, StreamT st,bool onDev=false);
/** @brief a_o = (a_o - a_i) * c */
  static void SubMulC_I(Field3D& a_o, const Field3D& a_i, const float& c, StreamT st,bool onDev=false);
/** @brief a_o = (a_o * a_i) * c */
  static void MulMulC_I(Field3D& a_o, const Field3D& a_i, const float& c, StreamT st,bool onDev=false);

/** @brief a_o = a_i + a_i1 * c */
  static void Add_MulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const float& c, StreamT st,bool onDev=false);
/** @brief a_o = a_o+ a_i1 * c */
  static void Add_MulC_I(Field3D& a_o, const Field3D& a_i1, const float& c, StreamT st,bool onDev=false);

/** @brief a_o = a_i * c + a_i1 */
  static void MulCAdd(Field3D& a_o, const Field3D& a_i, const float& c, const Field3D& a_i1, StreamT st,bool onDev=false);
/** @brief a_o = a_i * c - a_i1 */
  static void MulCSub(Field3D& a_o, const Field3D& a_i, const float& c, const Field3D& a_i1, StreamT st,bool onDev=false);

/** @brief a_o = a_o * c + a_i */
  static void MulCAdd_I(Field3D& a_o, const float& c, const Field3D& a_i, StreamT st,bool onDev=false);
/** @brief a_o = a_o * c - a_i */
  static void MulCSub_I(Field3D& a_o, const float& c, const Field3D& a_i, StreamT st,bool onDev=false);
/** @brief a_o = (a_i + a) * b */
  static void AddCMulC(Field3D& a_o, const Field3D& a_i, const Vec3Df& a, const Vec3Df& b, StreamT st,bool onDev=false);
/** @brief a_o = (a_i * a) + b */
  static void MulCAddC(Field3D& a_o, const Field3D& a_i, const Vec3Df& a, const Vec3Df& b, StreamT st,bool onDev=false);
/** @brief a_o = (a_o + a) * b */
  static void AddCMulC_I(Field3D& a_o, const Vec3Df& a, const Vec3Df& b, StreamT st,bool onDev=false);
/** @brief a_o = (a_o * a) + b */
  static void MulCAddC_I(Field3D& a_o, const Vec3Df& a, const Vec3Df& b, StreamT st,bool onDev=false);

/** @brief a_o = a_i * ca + a_i1 * cb */
  static void MulC_Add_MulC(Field3D& a_o, const Field3D& a_i, const float& ca,
                            const Field3D& a_i1, const float& cb, StreamT st,bool onDev=false);
/** @brief a_o = a_o * co + a_i1 * cb */
  static void MulC_Add_MulC_I(Field3D& a_o, const float& co,
                              const Field3D& a_i1, const float& cb, StreamT st,bool onDev=false);
/** @brief a_o = (a_i + a) * b + c */
  static void AddCMulCAddC(Field3D& a_o, const Field3D& a_i,
                           const float& a, const float& b, const float& c, StreamT st,bool onDev=false);
/** @brief a_o = (a_o + a) * b + c */
  static void AddCMulCAddC_I(Field3D& a_o, const float& a, const float& b, const float& c, StreamT st,bool onDev=false);
/** @brief a_o = (a_i + a) * b + c */
  static void AddCMulCAddC(Field3D& a_o, const Field3D& a_i, const Vec3Df& a,
                           const Vec3Df& b, const Vec3Df& c, StreamT st,bool onDev=false);
/** @brief a_o = (a_o + a) * b + c */
    static void AddCMulCAddC_I(Field3D& a_o, const Vec3Df& a, const Vec3Df& b, const Vec3Df& c, StreamT st,bool onDev=false);
/** @brief a_o = a_i + a_i1 * a_i2 * c */
    static void Add_MulMulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1,
                     const Field3D& a_i2, const float& c, StreamT st,bool onDev=false);
/** @brief a_o = a_o + a_i * a_i1 * c */
    static void Add_MulMulC_I(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const float& c, StreamT st,bool onDev=false);
/** @brief a_o = a_i - a_i1 * a_i2 * c */
    static void Sub_MulMulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1,
                     const Field3D& a_i2, const float& c, StreamT st,bool onDev=false);
/** @brief a_o = a_o - a_i * a_i1 * c */
  static void Sub_MulMulC_I(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1,
                            const float& c, StreamT st,bool onDev=false);

/*
 * Advance function
 */

/**
 * @brief
 * compose two hfields to get an hfield
 * 
 * f(x) = g(h(x))
 */
    template<BackgroundStrategy bg>
    static void ComposeHH(Field3D& a_f, const Field3D& a_g,
                          const Field3D& a_h, StreamT s);

/**
 * compose a velocity and h field to get an hfield
 * 
 * h(x) = g(x) + delta * v(g(x))
 */
    template<BackgroundStrategy bg>
    static void ComposeVH(Field3D& a_h, const Field3D& a_v,
                          const Field3D& a_g, const float& delta, StreamT s,bool onDev=false);

/**
 * compose a inverse velocity and h field to get an hfield
 * 
 * h(x) = g(x) - delta * v(g(x))
 */

    template<BackgroundStrategy bg>
    static void ComposeVInvH(Field3D& a_h, const Field3D& a_v, const Field3D& a_g, const float& delta, StreamT s,bool onDev=false);

/**
 * compose a h field and a velocify field to get an hfield
 * 
 * h(x) = g(x + delta * v(x))
 */
    template<BackgroundStrategy bg>
    static void ComposeHV(Field3D& a_h, const Field3D& a_g, const Field3D& a_v,
                          const float& delta, StreamT s,bool onDev=false);

/**
 * @brief compose a h field and an inverse velocity field to get an hfield
 * 
 * h(x) = g(x - delta * v(x))
 */
    template<BackgroundStrategy bg>
    static void ComposeHVInv(Field3D& a_h, const Field3D& a_g, const Field3D& a_v,
                             const float& delta, StreamT s,bool onDev=false);

/**
 * @brief precompose h field with translation
 *
 * creating h(x) = f(x + t)
 */
    template<BackgroundStrategy bg>
    static void ComposeTranslation(Field3D& a_h, const Field3D& a_f, const Vec3Df& t, StreamT s,bool onDev=false);

    static void Jacobian(Field3D& d_Xg, Field3D& d_Yg, Field3D& d_Zg,
                         const Field3D& d_h,
			 DiffT diffType = DIFF_CENTRAL,
			 BoundaryCondT bc = BC_APPROX,
			 StreamT s=NULL);
/**
 * @brief deform a field using hfield
 *
 * a_o(x) = a_i(h(x))
 */
    template<BackgroundStrategy bg>
    static void ApplyH(Field3D& a_o, const Field3D& a_i, const Field3D& h, StreamT s);
/**
 * @brief deform a field using displacement field
 * 
 * a_o(x) = a_i(x + delta * u)
 */
    template<BackgroundStrategy bg>
    static void ApplyV(Field3D& a_o, const Field3D& a_i, const Field3D& u,
                       const float& delta, StreamT s,bool onDev=false);
/**
 *
 * @brief deform a field using displacement field
 *
 * a_o(x) = a_i(x - delta * u)
 */
    template<BackgroundStrategy bg>
    static void ApplyVInv(Field3D& a_o, const Field3D& a_i, const Field3D& u,
                          const float& delta, StreamT s,bool onDev=false);

   static void SplatField(Field3D& a_o, const Field3D& a_i,
				 const Field3D& h, StreamT stream=NULL);

   static void SubVol(Field3D& a_o,const Field3D& a_i,
		      const Vec3Di& start, StreamT st=NULL);

   static void SetSubVol_I(Field3D& a_o,const Field3D& a_i,
			   const Vec3Di& start, StreamT st=NULL);


    template<BackgroundStrategy bg,  bool rescaleVector>
    static void Resample(Field3D& a_o, const Field3D& a_i, StreamT stream);

   static void ReprojectToUnitVec(Field3D& a_o, StreamT st=NULL);

   // normalize all vectors to have unit length, except those less
   // than length eps which will be set to zero
   static void NormalizeSafe(Field3D& a_o, const Field3D& a_i, const float& eps, StreamT st=NULL);
   static void NormalizeSafe_I(Field3D& a_o, const float& eps, StreamT st=NULL);

   static void Shrink(Field3D& a_o, const Field3D& a_i, const float& eps, StreamT st=NULL);
   static void Shrink_I(Field3D& a_o, const float& eps, StreamT st=NULL);

    // Given a deformation g, a displacement vector field v, and an inverse ginv,
    // of v, invert the deformation gt(x) = g(x) + dt v(g(x))
    template<BackgroundStrategy bg>
    static void FixedPointInverse(Field3D &ginv, const Field3D& g, unsigned int numIter=2, StreamT stream = NULL, bool onDev=false);

    static void UpdateInverse(Field3D &ginv0t1,Field3D &scratchV, const Field3D& ginv0t, const Field3D& w, unsigned int numIter = 3, StreamT stream = NULL, bool onDev = false);
////////////////////////////////////////////////////////////////////////////
// Lie algebra operations
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
// adjoint action of Diff on its Lie algebra
// Z = Ad_g X = |Dg|\circ g^{-1} X\circ g^{-1}
////////////////////////////////////////////////////////////////////////////
    template<BackgroundStrategy bg>
      static void Ad(Field3D& Z, const Field3D& g, const Field3D& X,
          StreamT stream, bool onDev=false);
////////////////////////////////////////////////////////////////////////////
// infinitesimal adjoint action, equal to negative Jacobi-Lie bracket
// Z = ad_X Y = DX Y - DY X
////////////////////////////////////////////////////////////////////////////
    static void AdInf(Field3D& Z, const Field3D& X, const Field3D& Y,
              StreamT stream, bool onDev=false);
////////////////////////////////////////////////////////////////////////////
// coadjoint action of Diff on its Lie coalgebra
// n = Ad_g^* m = (Dg)^T m\circ g |Dg|
////////////////////////////////////////////////////////////////////////////
    template<BackgroundStrategy bg>
    static void CoAd(Field3D& n, const Field3D& g, const Field3D& m,
              StreamT stream, bool onDev=false);
////////////////////////////////////////////////////////////////////////////
// infinitesimal coadjoint action
// n = ad_X^* m = (DX)^T m + div(m \otimes X)
////////////////////////////////////////////////////////////////////////////
    static void CoAdInf(Field3D& n, const Field3D& X, const Field3D& m,
              StreamT stream, bool onDev=false);
////////////////////////////////////////////////////////////////////////////
// computes tensor divergence of outer product of two vector fields 
// Z = div(X \otimes Y)
////////////////////////////////////////////////////////////////////////////
    static void DivergenceTensor(Field3D& Z, const Field3D& X, const Field3D& Y,
              StreamT stream, bool onDev=false);
////////////////////////////////////////////////////////////////////////////
// computes jacobian of X times Y
// Z = DX Y
////////////////////////////////////////////////////////////////////////////
    static void JacobianXY(Field3D& Z, const Field3D& X, const Field3D& Y, StreamT stream=NULL, bool onDev=false) ;
////////////////////////////////////////////////////////////////////////////
// computes jacobian of X transpose times Y
// Z = (DX)' Y
////////////////////////////////////////////////////////////////////////////
    static void JacobianXtY(Field3D& Z, const Field3D& X, const Field3D& Y, StreamT stream=NULL, bool onDev=false) ;

};

} // end namespace PyCA

#endif

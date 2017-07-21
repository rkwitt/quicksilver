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

#ifndef __CFIELD3D_OPERS_H
#define __CFIELD3D_OPERS_H

#include <pycaConst.h>
#include <estream.h>
#include <Vec3D.h>

namespace PyCA {

class Image3D;
class Field3D;

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
  BACKGROUND_STRATEGY_CLAMP, BACKGROUND_STRATEGY_WRAP
*/
class CFieldOpers{
public:
/*
 * Set all the value of the Vector field with single value Vector3f(v)
 */
    static void SetMem(Field3D& a_o, const Vec3Df& v, StreamT st,bool onDev=false);
/**  
 * Set all the value of the Vector field with single float value
 * (normaly used to zero out the memory)
 */
    static void SetMem(Field3D& a_o, const float& c, StreamT st,bool onDev=false);
    
/** @defgroup operator
 * Basic operator on the Vector field
 * @{
 */
/** @brief a_o = a_i + c */
    static void AddC(Field3D& a_o, const Field3D& a_i, const Vec3Df& c, StreamT st,bool onDev=false);
/** @brief a_o = a_i - c */
    static void SubC(Field3D& a_o, const Field3D& a_i, const Vec3Df& c, StreamT st,bool onDev=false);
/** @brief a_o = a_i * c */
    static void MulC(Field3D& a_o, const Field3D& a_i, const Vec3Df& c, StreamT st,bool onDev=false);
/** @brief a_o = a_i * c */
    static void MulC(Field3D& a_o, const Field3D& a_i, const float& c, StreamT st,bool onDev=false);

/** @brief a_o = a_i + c */
    static void AddC_I(Field3D& a_o, const Vec3Df& c, StreamT st,bool onDev=false);
/** @brief a_o -= c */
    static void SubC_I(Field3D& a_o, const Vec3Df& c, StreamT st,bool onDev=false);
/** @brief a_o *= c */
    static void MulC_I(Field3D& a_o, const Vec3Df& c, StreamT st,bool onDev=false);
/** @brief a_o *= c */
    static void MulC_I(Field3D& a_o, const float& c, StreamT st,bool onDev=false);

/** @brief a_o = a_i + a_i1 */
    static void Add(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, StreamT st);
/** @brief a_o = a_i - a_i1 */
    static void Sub(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, StreamT st);
/** @brief a_o = a_i * a_i1 */
    static void Mul(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, StreamT st);

/** @brief a_o += a_i1 */
    static void Add_I(Field3D& a_o, const Field3D& a_i1, StreamT st);
/** @brief a_o -= a_i1 */
    static void Sub_I(Field3D& a_o, const Field3D& a_i1, StreamT st);
/** @brief a_o *= a_i1 */
    static void Mul_I(Field3D& a_o, const Field3D& a_i1, StreamT st);

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
    static void MulC_Add_MulC(Field3D& a_o,
                              const Field3D& a_i, const float& ca,
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
    static void AddCMulCAddC(Field3D& a_o, const Field3D& a_i, const Vec3Df& a, const Vec3Df& b, const Vec3Df& c, StreamT st,bool onDev=false);
    
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
    static void Sub_MulMulC_I(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const float& c, StreamT st,bool onDev=false);

/*
 * Advance function
 */ 

/*
 * compose two hfields to get an hfield
 *  f(x) = g(h(x))
 */
    template<BackgroundStrategy bg>
    static void ComposeHH(Field3D& a_f, const Field3D& a_g,
                          const Field3D& a_h, StreamT s);

/*
 * compose a velocity and h field to get an hfield
 *  h(x) = g(x) + delta * v(g(x))
 */
    template<BackgroundStrategy bg>
    static void ComposeVH(Field3D& a_h, const Field3D& a_v,
                          const Field3D& a_g, const float& delta, StreamT s,bool onDev=false);

/*
 * compose a inverse velocity and h field to get an hfield
 *  h(x) = g(x) - delta * v(g(x))
 */

    template<BackgroundStrategy bg>
    static void ComposeVInvH(Field3D& a_h, const Field3D& a_v, const Field3D& a_g, const float& delta, StreamT s,bool onDev=false);

/*
 * compose a h field and a velocify field to get an hfield
 * h(x) = g(x + delta * v(x))
 */
    template<BackgroundStrategy bg>
    static void ComposeHV(Field3D& a_h, const Field3D& a_g, const Field3D& a_v,
                          const float& delta, StreamT s,bool onDev=false);

/*
 * compose a h field and an inverse velocity field to get an hfield
 * h(x) = g(x - delta * v(x))
 */
    template<BackgroundStrategy bg>
    static void ComposeHVInv(Field3D& a_h, const Field3D& a_g, const Field3D& a_v,
                             const float& delta, StreamT s,bool onDev=false);

/*
 * precompose h field with translation
 * creating h(x) = f(x + t)
 */
    template<BackgroundStrategy bg>
    static void ComposeTranslation(Field3D& a_h, const Field3D& a_f, const Vec3Df& t, StreamT s,bool onDev=false);

/*
 * deform a field using hfield
 * a_o(x) = a_i(h(x))
 */
    template<BackgroundStrategy bg>
    static void ApplyH(Field3D& a_o, const Field3D& a_i, const Field3D& h, StreamT s);
/*
 * deform a field using displacement field
 * a_o(x) = a_i(x + delta * u)
 */
    template<BackgroundStrategy bg>
    static void ApplyV(Field3D& a_o, const Field3D& a_i, const Field3D& u,
                       const float& delta, StreamT s,bool onDev=false);
/*
 * deform a field using displacement field
 * a_o(x) = a_i(x - delta * u)
 */
    template<BackgroundStrategy bg>
    static void ApplyVInv(Field3D& a_o, const Field3D& a_i, const Field3D& u,
                          const float& delta, StreamT s,bool onDev=false);

   static void SplatField(Field3D& a_o, const Field3D& a_i, 
				 const Field3D& h, StreamT st);

   static void SubVol(Field3D& a_o, const Field3D& a_i, 
		      const Vec3Di& start, StreamT st=NULL);

   static void SetSubVol_I(Field3D& a_o, const Field3D& a_i, 
			   const Vec3Di& start, StreamT st=NULL);
   
    template<BackgroundStrategy bg,  bool rescaleVector>
    static void Resample(Field3D& a_o, const Field3D& a_i, StreamT stream);

   static void ReprojectToUnitVec(Field3D& a_o, StreamT st=NULL);

   static void NormalizeSafe(Field3D& a_o, const Field3D& a_i, const float& eps, StreamT st=NULL);
   static void NormalizeSafe_I(Field3D& a_o, const float& eps, StreamT st=NULL);

   static void Shrink(Field3D& a_o, const Field3D& a_i, const float& eps, StreamT st=NULL);
   static void Shrink_I(Field3D& a_o, const float& eps, StreamT st=NULL);

    template<BackgroundStrategy bg>
    static void FixedPointInverse(Field3D &ginv, const Field3D& g, unsigned int numIter=2, StreamT stream = NULL, bool onDev=false);
    static void UpdateInverse(Field3D &ginv0t1,Field3D &scratchV, const Field3D& ginv0t, const Field3D& w, unsigned int numIter = 3, StreamT stream = NULL, bool onDev=false);

   /*
    * Lie algebra methods
    */

/*
 * Adjoint action of Diff on its Lie algebra
 * Z = Ad_g X = |Dg|\circ g^{-1} X\circ g^{-1}
 */
    template<BackgroundStrategy bg>
      static void Ad(Field3D& Z, const Field3D& g, const Field3D& X,
                          StreamT s,bool onDev=false);
/*
 * infinitesimal adjoint action
 * Z = ad_X Y = DX Y - DY X
 */
    static void AdInf(Field3D& Z, const Field3D& X, const Field3D& Y,
                          StreamT s,bool onDev=false);
/*
 * Coadjoint action of Diff on its Lie algebra
 * n = Ad_g^* m = (Dg)^T m\circ g |Dg|
 */
    template<BackgroundStrategy bg>
    static void CoAd(Field3D& n, const Field3D& g, const Field3D& m,
                          StreamT s,bool onDev=false);
/*
 * infinitesimal coadjoint action
 * n = ad_X^* m = (DX)^T m + div(m \otimes X)
 */
    static void CoAdInf(Field3D& n, const Field3D& X, const Field3D& m,
                          StreamT s,bool onDev=false);

/*
 * computes tensor divergence of outer product of two vector fields 
 * Z = div(X \otimes Y)
 */
    static void DivergenceTensor(Field3D& Z, const Field3D& X, const Field3D& Y,
				 StreamT s,bool onDev=false);


/*
 * Jacobian X times Y
 * Z = DX Y
 */
    static void JacobianXY(Field3D& Z, const Field3D& X, const Field3D& Y, StreamT s=NULL,bool onDev=false);

/*
 * Jacobian X transpose times Y
 * Z = (DX)' Y
 */
    static void JacobianXtY(Field3D& Z, const Field3D& X, const Field3D& Y, StreamT s=NULL,bool onDev=false);

};

} // end namespace PyCA

#endif

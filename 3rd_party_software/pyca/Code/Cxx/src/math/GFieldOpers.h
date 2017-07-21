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

#ifndef __GFIELD3D_OPERS_H
#define __GFIELD3D_OPERS_H

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
  BACKGROUND_STRATEGY_CLAMP, BACKGROUND_STRATEGY_WARP
*/
class GFieldOpers{
public:
/*
 * Copy a_i to a_o
 */
    static void Copy(Field3D& a_o, const Field3D& a_i, StreamT stream);
/*
 * Set all the value of the Vector field with single value Vector3f(v)
 */
    static void SetMem(Field3D& a_o, const Vec3Df& v, StreamT st,bool onDev);
/**  
 * Set all the value of the Vector field with single float value
 * (normaly used to zero out the memory)
 */
    static void SetMem(Field3D& a_o, const float& c, StreamT st,bool onDev);
    
/*
 * Advance function
 */ 

/*
 * Compose two hfields to get an hfield
 *  f(x) = g(h(x))
 */
    template<BackgroundStrategy bg>
    static void ComposeHH(Field3D& d_f, const Field3D& d_g,
                          const Field3D& d_h, StreamT s);

/*
 * Compose a velocity and h field to get an hfield
 *  h(x) = g(x) + delta * v(g(x))
 */
    template<BackgroundStrategy bg>
    static void ComposeVH(Field3D& a_h, const Field3D& a_v,
                          const Field3D& a_g, const float& delta, StreamT s,bool onDev);

/*
 * Compose a inverse velocity and h field to get an hfield
 *  h(x) = g(x) - delta * v(g(x))
 */

    template<BackgroundStrategy bg>
    static void ComposeVInvH(Field3D& a_h, const Field3D& a_v, const Field3D& a_g, const float& delta, StreamT s,bool onDev);

/*
 * Compose a h field and a velocify field to get an hfield
 * h(x) = g(x + delta * v(x))
 */
    template<BackgroundStrategy bg>
    static void ComposeHV(Field3D& a_h, const Field3D& a_g, const Field3D& a_v,
                          const float& delta, StreamT s,bool onDev);

/*
 * Compose a h field and an inverse velocity field to get an hfield
 * h(x) = g(x - delta * v(x))
 */
    template<BackgroundStrategy bg>
    static void ComposeHVInv(Field3D& a_h, const Field3D& a_g, const Field3D& a_v,
                             const float& delta, StreamT s,bool onDev);

/*
 * preCompose h field with translation
 * creating h(x) = f(x + t)
 */
    template<BackgroundStrategy bg>
    static void ComposeTranslation(Field3D& a_h, const Field3D& a_f, const Vec3Df& t, StreamT s,bool onDev);

    static void Jacobian(Field3D& d_Xg, Field3D& d_Yg, Field3D& d_Zg, 
                         const Field3D& d_h, StreamT s);
    
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
                       const float& delta, StreamT s,bool onDev);
/*
 * deform a field using displacement field
 * a_o(x) = a_i(x - delta * u)
 */
    template<BackgroundStrategy bg>
    static void ApplyVInv(Field3D& a_o, const Field3D& a_i, const Field3D& u,
                          const float& delta, StreamT s,bool onDev);

   static void SplatField(Field3D& d_o, const Field3D& d_i, 
				 const Field3D& d_h, StreamT stream);

   static void SubVol(Field3D& d_o, const Field3D& d_i, 
		      const Vec3Di& start, StreamT st=NULL);

   static void SetSubVol_I(Field3D& d_o, const Field3D& d_i, 
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

    static void UpdateInverse(Field3D &ginv0t1,Field3D &scratchV, const Field3D& ginv0t, const Field3D& w, unsigned int numIter = 3, StreamT stream = NULL, bool onDev = false);

// Lie algebra methods

/*
 * Adjoint action of Diff on its Lie algebra
 * This is just the pushforward
 * Z = Ad_g X = |Dg|\circ g^{-1} X\circ g^{-1}
 */
    template<BackgroundStrategy bg> 
      static void Ad(Field3D& Z, const Field3D& g, const Field3D& X, StreamT s,bool onDev=true);
/*
 * infinitesimal adjoint action
 * Z = ad_X Y = DX Y - DY X
 */
    static void AdInf(Field3D& Z, const Field3D& X, const Field3D& Y, StreamT s,bool onDev=true);
/*
 * Coadjoint action of Diff on its Lie algebra
 * n = Ad_g^* m = (Dg)^T m\circ g |Dg|
 */
    template<BackgroundStrategy bg> 
    static void CoAd(Field3D& n, const Field3D& g, const Field3D& m, StreamT s,bool onDev=true);

/*
 * infinitesimal coadjoint action
 * n = ad_X^* m = (DX)^T m + div(m \otimes X)
 */
    static void CoAdInf(Field3D& n, const Field3D& X, const Field3D& m, StreamT s,bool onDev=true);

/*
 * computes tensor divergence of outer product of two vector fields 
 * Z = div(X \otimes Y)
 */
    static void DivergenceTensor(Field3D& Z, const Field3D& X, const Field3D& Y, StreamT s,bool onDev=true);


/*
 * Jacobian X times Y
 * Z = DX Y
 */
    static void JacobianXY(Field3D& Z, const Field3D& X, const Field3D& Y, StreamT s=NULL,bool onDev=true);

/*
 * Jacobian X transpose times Y
 * Z = (DX)' Y
 */
    static void JacobianXtY(Field3D& Z, const Field3D& X, const Field3D& Y, StreamT s=NULL,bool onDev=true);

};

} // end namespace PyCA

#endif

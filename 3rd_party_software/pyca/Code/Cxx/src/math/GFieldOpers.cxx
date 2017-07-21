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

#include "GFieldOpers.h"
#include "GFieldOperKernels.h"
#include <pycaUtils.h>
#include <mem.h>
#include "interp.h"
#include "GSplat.h"
#include <MemOpers.h>
#include <Field3D.h>

#include "GImageFieldOpers.h"
#include "GImageOpers.h"

#include "FOpers.h"

#include "FiniteDiff.h"

namespace PyCA {

////////////////////////////////////////////////////////////////////////////
// compose a velocity and hfield to get an hfield
// h(x) = g(x) + delta * v(g(x))
////////////////////////////////////////////////////////////////////////////
__device__ __constant__ float c_delta;
__device__ __constant__ float c_trans[3];

////////////////////////////////////////////////////////////////////////////
// compose a two hfields to get an hfield
// f(x) = g(h(x))
////////////////////////////////////////////////////////////////////////////
template<BackgroundStrategy bg>
void GFieldOpers::
ComposeHH(Field3D& d_f, const Field3D& d_g,
	       const Field3D& d_h, StreamT s)
{
   MK_CHECK_HFIELD_BACKGROUND(bg);
   ApplyH<bg>(d_f, d_g, d_h, s);
}

template<BackgroundStrategy bg>
void GFieldOpers::
ComposeVH(Field3D& d_h, const Field3D& d_v, const Field3D& d_g, const float& delta, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(d_h, d_v, d_g);
    Vec3Di size = d_h.size();
    Vec3Df sp   = d_h.spacing();
    PyCA::ComposeVH<true, bg>(d_h.x, d_h.y, d_h.z,
                          d_v.x, d_v.y, d_v.z,
                          d_g.x, d_g.y, d_g.z, delta,
                          size.x, size.y, size.z,
                          sp.x, sp.y, sp.z, stream, onDev);
}


template<BackgroundStrategy bg>
void GFieldOpers::
ComposeVInvH(Field3D& d_h, const Field3D& d_v, const Field3D& d_g, const float& delta, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(d_h, d_v, d_g);
    Vec3Di size = d_h.size();
    Vec3Df sp   = d_h.spacing();
    PyCA::ComposeVH<false, bg>(d_h.x, d_h.y, d_h.z,
                           d_v.x, d_v.y, d_v.z,
                           d_g.x, d_g.y, d_g.z, delta,
                           size.x, size.y, size.z,
                           sp.x, sp.y, sp.z, stream, onDev);
}

template<BackgroundStrategy bg>
void GFieldOpers::
ComposeHV(Field3D& d_h, const Field3D& d_g, const Field3D& d_v, const float& delta, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(d_h, d_g, d_v);
    Vec3Di size = d_h.size();
    Vec3Df sp   = d_h.spacing();

    PyCA::ComposeHV<true, bg>
	(d_h.x, d_h.y, d_h.z,
	 d_g.x, d_g.y, d_g.z,
	 d_v.x, d_v.y, d_v.z, delta,
	 size.x, size.y, size.z,
	 sp.x, sp.y, sp.z, stream, onDev);
}


template<BackgroundStrategy bg>
void GFieldOpers::
ComposeHVInv(Field3D& d_h, const Field3D& d_g, const Field3D& d_v, const float& delta, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(d_h, d_g, d_v);
    Vec3Di size = d_h.size();
    Vec3Df sp   = d_h.spacing();

    PyCA::ComposeHV<false, bg>
	(d_h.x, d_h.y, d_h.z,
	 d_g.x, d_g.y, d_g.z,
	 d_v.x, d_v.y, d_v.z, delta,
	 size.x, size.y, size.z,
	 sp.x, sp.y, sp.z, 
	 stream, onDev);
}

template<BackgroundStrategy bg>
void GFieldOpers::
ComposeTranslation(Field3D& d_o, const Field3D& d_i, const Vec3Df& t, StreamT stream, bool onDev) {
    MK_CHECK2_SIZE(d_o, d_i);
    Vec3Di sz = d_o.size();

    PyCA::ComposeTranslation<bg>
	(d_o.x, d_o.y, d_o.z,
	 d_i.x, d_i.y, d_i.z,
	 sz, t, 
	 stream, onDev);
}

template<BackgroundStrategy bg>
void GFieldOpers::
ApplyH(Field3D& d_o, const Field3D& d_i, const Field3D& d_h, StreamT stream) {
    MK_CHECK3_SIZE(d_o, d_i, d_h);
    Vec3Di sz = d_o.size();

    PyCA::ApplyH<bg>
	(d_o.x, d_o.y, d_o.z,
	 d_i.x, d_i.y, d_i.z,
	 d_h.x, d_h.y, d_h.z, 
	 sz, 
	 stream);
}


template<BackgroundStrategy bg>
void GFieldOpers::
ApplyV(Field3D& d_o, const Field3D& d_i, const Field3D& d_u, const float& delta, StreamT stream, bool onDev) {

    MK_CHECK3_SIZE(d_o, d_i, d_u);
    Vec3Di sz = d_o.size();
    Vec3Df sp   = d_o.spacing();

    PyCA::ApplyV<true, bg>(d_o.x, d_o.y, d_o.z, 
			   d_i.x, d_i.y, d_i.z, 
			   d_u.x, d_u.y, d_u.z, 
			   sz, sp, delta, 
			   stream, onDev);
}


template<BackgroundStrategy bg>
void GFieldOpers::
ApplyVInv(Field3D& d_o, const Field3D& d_i, const Field3D& d_u, const float& delta, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(d_o, d_i, d_u);
    Vec3Di sz = d_o.size();
    Vec3Df sp   = d_o.spacing();

    PyCA::ApplyV<false, bg>(d_o.x, d_o.y, d_o.z, 
			   d_i.x, d_i.y, d_i.z, 
			   d_u.x, d_u.y, d_u.z, 
			   sz, sp, delta, 
			   stream, onDev);
}

void GFieldOpers::
SplatField(Field3D& d_o, const Field3D& d_i, 
		  const Field3D& d_h, StreamT stream) 
{
   GImageFieldOpers::Splat(d_o.x, d_h, d_i.x, false, stream);
   GImageFieldOpers::Splat(d_o.y, d_h, d_i.y, false, stream);
   GImageFieldOpers::Splat(d_o.z, d_h, d_i.z, false, stream);
}

void GFieldOpers::
SubVol(Field3D& d_o, const Field3D& d_i, 
       const Vec3Di& start, StreamT st)
{
   Vec3Di oSize = d_o.size();
   Vec3Di iSize = d_i.size();
   GImageOpers::SubVol(d_o.x, d_i.x, oSize, iSize, start, st);
   GImageOpers::SubVol(d_o.y, d_i.y, oSize, iSize, start, st);
   GImageOpers::SubVol(d_o.z, d_i.z, oSize, iSize, start, st);
}

void GFieldOpers::
SetSubVol_I(Field3D& d_o, const Field3D& d_i, 
	    const Vec3Di& start, StreamT st)
{
   Vec3Di oSize = d_o.size();
   Vec3Di iSize = d_i.size();
   GImageOpers::SetSubVol_I(d_o.x, d_i.x, oSize, iSize, start, st);
   GImageOpers::SetSubVol_I(d_o.y, d_i.y, oSize, iSize, start, st);
   GImageOpers::SetSubVol_I(d_o.z, d_i.z, oSize, iSize, start, st);
}

template<BackgroundStrategy bg,  bool rescaleVector>
void GFieldOpers::
Resample(Field3D& d_o, const Field3D& d_i, StreamT stream)
{
    if (d_o.size() == d_i.size()) {
        Opers::Copy(d_o, d_i, stream);
        return;
    }
    Vec3Di oSz = d_o.size();
    Vec3Di iSz = d_i.size();
    
    PyCA::Resample<bg, rescaleVector>
        (d_o.x, d_o.y, d_o.z, oSz, 
         d_i.x, d_i.y, d_i.z, iSz,
	 stream);

}

void GFieldOpers::
ReprojectToUnitVec(Field3D& d_o, StreamT st)
{

   Vec3Di sz = d_o.size();

   PyCA::ReprojectToUnitVec
       (d_o.x, d_o.y, d_o.z, sz, st);
}

void GFieldOpers::
NormalizeSafe(Field3D& d_o, const Field3D& d_i, const float& eps, 
	      StreamT st)
{

   Vec3Di sz = d_o.size();

   PyCA::NormalizeSafe
       (d_o.x, d_o.y, d_o.z,
	d_i.x, d_i.y, d_i.z,
	sz, eps, st);
}

void GFieldOpers::
NormalizeSafe_I(Field3D& d_o, const float& eps, 
		StreamT st)
{
   NormalizeSafe(d_o, d_o, eps, st);
}

void GFieldOpers::
Shrink(Field3D& d_o, const Field3D& d_i, const float& eps, 
	      StreamT st)
{

   Vec3Di sz = d_o.size();

   PyCA::Shrink
       (d_o.x, d_o.y, d_o.z,
	d_i.x, d_i.y, d_i.z,
	sz, eps, st);
}

void GFieldOpers::
Shrink_I(Field3D& d_o, const float& eps, 
		StreamT st)
{
   Shrink(d_o, d_o, eps, st);
}

template<BackgroundStrategy bg> 
void GFieldOpers::
FixedPointInverse(Field3D &ginv, const Field3D& g, unsigned int numIter, StreamT stream, bool onDev)
{
    MK_CHECK_HFIELD_BACKGROUND(bg);
    MK_CHECK2_SIZE(ginv,g);

    Vec3Di sz = g.size();

    PyCA::FixedPointInverse<bg>
	(ginv.x, ginv.y, ginv.z,
	 g.x, g.y, g.z,
	 sz, numIter,
	 stream, onDev);
    CudaUtils::CheckCUDAError(__FILE__,__LINE__, "Fixed point kernel");
}

void GFieldOpers::
UpdateInverse(Field3D &ginv0t1,Field3D &scratchV, const Field3D& ginv0t, const Field3D& w, unsigned int numIter, StreamT stream, bool onDev)
{
  MK_CHECK4_SIZE(ginv0t1,scratchV,ginv0t,w);

  Vec3Di sz = ginv0t.grid().size();
  Vec3Df sp   = ginv0t.grid().spacing();

  dim3 threads(16,16);
  dim3 grids(iDivUp(sz.x, threads.x), iDivUp(sz.y, threads.y));

  Opers::Copy(ginv0t1, ginv0t, stream);
  for(unsigned int iter=0;iter<numIter;iter++){    
    // compute w\circ g_{t+1,0}
    ApplyH<BACKGROUND_STRATEGY_PARTIAL_ZERO>(scratchV, w, ginv0t1, stream);

    // compute Id - w\circ g_{t+1,0}    
    PyCA::updateInverseSubFromIndentity
	(scratchV.x, scratchV.y, scratchV.z,
	 sz, sp, stream);
    CudaUtils::CheckCUDAError(__FILE__,__LINE__, "Update inverse: SubFromIdentity kernel");
    // compute g_{t+1,0} = g_{t,0}\circ(Id - w\circ g_{t+1,0})
    ComposeHH<BACKGROUND_STRATEGY_PARTIAL_ID>(ginv0t1,ginv0t,scratchV, stream);
  }
}
// Lie algebra methods

/*
 * Adjoint action of Diff on its Lie algebra
 * This is just the pushforward
 * Z = Ad_g X = |Dg|\circ g^{-1} X\circ g^{-1}
 */
template<BackgroundStrategy bg> 
void GFieldOpers::
Ad(Field3D& Z, const Field3D& g, const Field3D& X,
                          StreamT s,bool onDev)
{
    MK_CHECK3_SIZE(Z,g,X);
    MK_CHECK_VFIELD_BACKGROUND(bg);

    Vec3Di sz = g.size();
    Vec3Df sp = g.spacing();
    PyCA::Ad<bg>
      (Z.x, Z.y, Z.z,
       g.x, g.y, g.z,
       X.x, X.y, X.z,
       sz, sp,
       s, onDev);
    CudaUtils::CheckCUDAError(__FILE__,__LINE__, "Adjoint kernel");
}

/*
 * infinitesimal adjoint action
 * Z = ad_X Y = DX Y - DY X
 */
void GFieldOpers::
AdInf(Field3D& Z, const Field3D& X, const Field3D& Y,
                          StreamT s,bool onDev)
{
    MK_CHECK3_SIZE(Z,X,Y);

    Vec3Di sz = Z.size();
    Vec3Df sp = Z.spacing();

    PyCA::AdInf
	(Z.x, Z.y, Z.z,
	 X.x, X.y, X.z,
	 Y.x, Y.y, Y.z,
	 sz, sp, 
	 s, onDev);

    CudaUtils::CheckCUDAError(__FILE__,__LINE__, "Little adjoint kernel");
}

/*
 * Jacobian X times Y
 * Z = DX Y
 */
void GFieldOpers::
JacobianXY(Field3D& Z, const Field3D& X, const Field3D& Y,
                          StreamT s,bool onDev)
{
    MK_CHECK3_SIZE(Z,X,Y);

   Vec3Di sz = Z.size();
   Vec3Df sp = Z.spacing();

   PyCA::JacobianXY
       (Z.x, Z.y, Z.z,
	X.x, X.y, X.z,
	Y.x, Y.y, Y.z,
	sz, sp,
	s, onDev);

    CudaUtils::CheckCUDAError(__FILE__,__LINE__, "JacobianXY kernel");

}


/*
 * Jacobian X transpose times Y
 * Z = (DX)' Y
 */
void GFieldOpers::
JacobianXtY(Field3D& Z, const Field3D& X, const Field3D& Y,
                          StreamT s,bool onDev)
{
    MK_CHECK3_SIZE(Z,X,Y);

   Vec3Di sz = Z.size();
   Vec3Df sp = Z.spacing();

   PyCA::JacobianXtY
       (Z.x, Z.y, Z.z,
	X.x, X.y, X.z,
	Y.x, Y.y, Y.z,
	sz, sp,
	s, onDev);

    CudaUtils::CheckCUDAError(__FILE__,__LINE__, "JacobianXtY kernel");

}

/*
 * Coadjoint action of Diff on its Lie algebra
 * n = Ad_g^* m = (Dg)^T m\circ g |Dg|
 */

template<BackgroundStrategy bg> 
void GFieldOpers::
CoAd(Field3D& n, const Field3D& g, const Field3D& m,
                          StreamT s,bool onDev)
{
    MK_CHECK3_SIZE(n,g,m);
    MK_CHECK_VFIELD_BACKGROUND(bg);

   Vec3Di sz = g.size();
   Vec3Df sp = g.spacing();

   PyCA::CoAd<bg>
       (n.x, n.y, n.z,
	g.x, g.y, g.z,
	m.x, m.y, m.z,
	sz, sp,
	s, onDev);

    CudaUtils::CheckCUDAError(__FILE__,__LINE__, "Coadjoint kernel");
}

/*
 * infinitesimal coadjoint action
 * n = ad_X^* m = (DX)^T m + div(m \otimes X)
 */
void GFieldOpers::
CoAdInf(Field3D& n, const Field3D& X, const Field3D& m,
                          StreamT s,bool onDev)
{
    MK_CHECK3_SIZE(n,X,m);
    
    Vec3Di sz = n.size();
    Vec3Df sp = n.spacing();
    
    PyCA::CoAdInf
	(n.x, n.y, n.z,
	 X.x, X.y, X.z,
	 m.x, m.y, m.z,
	 sz, sp,
	 s, onDev);

    CudaUtils::CheckCUDAError(__FILE__,__LINE__, "Little coadjoint kernel");

}


/*
 * computes tensor divergence of outer product of two vector fields 
 * Z = div(X \otimes Y)
 */
void GFieldOpers::
DivergenceTensor(Field3D& Z, const Field3D& X, const Field3D& Y,
                          StreamT s,bool onDev)
{
    MK_CHECK3_SIZE(Z,X,Y);
    
    Vec3Di sz = Z.size();
    Vec3Df sp = Z.spacing();
    
    PyCA::DivergenceTensor
	(Z.x, Z.y, Z.z,
	 X.x, X.y, X.z,
	 Y.x, Y.y, Y.z,
	 sz, sp,
	 s, onDev);

    CudaUtils::CheckCUDAError(__FILE__,__LINE__, "Divergence of outer product of two vectors kernel");

}

//Instantiation
#include "GFieldOpers_inst.cxx"

} // end namespace PyCA

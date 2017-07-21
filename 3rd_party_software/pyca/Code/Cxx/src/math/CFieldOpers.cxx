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

#include "CFieldOpers.h"
#include <conditionMacro.h>
#include <Field3D.h>
#include "interp.h"
#include "CSplat.h"
#include <MemOpers.h>
#include "CImageFieldOpers.h"
#include "CImageOpers.h"
#include "CFiniteDiff.h"
#include "FiniteDiff.h"
#include "FOpers.h"

namespace PyCA {

////////////////////////////////////////////////////////////////////////////
// compose a two hfields to get an hfield
// f(x) = g(h(x))
////////////////////////////////////////////////////////////////////////////
template<BackgroundStrategy bg>
void CFieldOpers::
ComposeHH(Field3D& a_f, const Field3D& a_g,
	  const Field3D& a_h, StreamT s)
{
   MK_CHECK_HFIELD_BACKGROUND(bg);
   ApplyH<bg>(a_f, a_g, a_h, s);
}

////////////////////////////////////////////////////////////////////////////
// compose a velocity and hfield to get an hfield
// h(x) = g(x) + delta * v(g(x))
////////////////////////////////////////////////////////////////////////////
template<bool fwd, BackgroundStrategy bg>
void ComposeVH(Field3D& h, const Field3D& v, const Field3D& g, float delta){
    MK_CHECK3_SIZE(h, v, g);
    MK_CHECK_VFIELD_BACKGROUND(bg);
    
    Vec3Di size = v.grid().size();
    Vec3Df sp   = v.grid().spacing();
    Vec3Df iSp  = Vec3Df(1.f/sp.x, 1.f/sp.y, 1.f/sp.z);
    
    size_t id = 0;
    for (int k=0; k < size.z; ++k)
        for (int j=0; j < size.y; ++j)
            for (int i=0; i < size.x; ++i, ++id) {
                float gx = g.x[id], gy = g.y[id], gz = g.z[id];
                float vgx, vgy, vgz;
                triLerp<bg>(vgx, vgy, vgz,
                            v.x, v.y, v.z,
                            gx, gy, gz,
                            size.x, size.y, size.z);
                if (fwd) {
                    h.x[id] = gx + vgx * iSp.x * delta;
                    h.y[id] = gy + vgy * iSp.y * delta;
                    h.z[id] = gz + vgz * iSp.z * delta;
                } else {
                    h.x[id] = gx - vgx * iSp.x * delta;
                    h.y[id] = gy - vgy * iSp.y * delta;
                    h.z[id] = gz - vgz * iSp.z * delta;
                }
            }
}

template<BackgroundStrategy bg>
void CFieldOpers::
ComposeVH(Field3D& h, const Field3D& v, const Field3D& g, const float& delta, StreamT st, bool onDev) {
    PyCA::ComposeVH<true, bg>(h, v, g, delta);
}



template<BackgroundStrategy bg>
void CFieldOpers::
ComposeVInvH(Field3D& h, const Field3D& v, const Field3D& g, const float& delta, StreamT st, bool onDev) {
    PyCA::ComposeVH<false, bg>(h, v, g, delta);
}

/**
 * compose a h field and a velocify field to get an hfield
 * h(x) = g(x+ delta * v(x))
 *
 * davisb 2007
 */
template<bool fwd, BackgroundStrategy bg>
void ComposeHV(Field3D& h, const Field3D& g, const Field3D& v, const float& delta){
    MK_CHECK3_SIZE(h, g, v);
    MK_CHECK_HFIELD_BACKGROUND(bg);
    
    Vec3Di size = v.grid().size();
    Vec3Df sp   = v.grid().spacing();
    Vec3Df iSp  = Vec3Df(1.f/sp.x, 1.f/sp.y, 1.f/sp.z);
        
    size_t id = 0;
    float x,y,z;
    
    for (int k=0; k < size.z; ++k)
        for (int j=0; j < size.y; ++j)
            for (int i=0; i < size.x; ++i, ++id) {
		if(fwd){
		    x = i + v.x[id] * iSp.x * delta;
		    y = j + v.y[id] * iSp.y * delta;
		    z = k + v.z[id] * iSp.z * delta;
		}else{
		    x = i - v.x[id] * iSp.x * delta;
		    y = j - v.y[id] * iSp.y * delta;
		    z = k - v.z[id] * iSp.z * delta;
		}

                float hx, hy, hz;
                triLerp<bg>(hx, hy, hz,
                            g.x, g.y, g.z,
                            x, y, z,
                            size.x, size.y, size.z);
                
                h.x[id] = hx;
                h.y[id] = hy;
                h.z[id] = hz;
            }
}

template<BackgroundStrategy bg>
void CFieldOpers::
ComposeHV(Field3D& h, const Field3D& g, const Field3D& v, const float& delta, StreamT st, bool onDev){
    PyCA::ComposeHV<true, bg>(h, g, v, delta);
}

template<BackgroundStrategy bg>
void CFieldOpers::
ComposeHVInv(Field3D& h, const Field3D& g, const Field3D& v, const float& delta, StreamT st, bool onDev){
    PyCA::ComposeHV<false, bg>(h, g, v, delta);
}

////////////////////////////////////////////////////////////////////////////
// compose field with translation
// creating a_o(x) = a_i(x + t)
////////////////////////////////////////////////////////////////////////////
template<BackgroundStrategy bg>
void CFieldOpers::
ComposeTranslation(Field3D& a_o, const Field3D& a_i, const Vec3Df& t, StreamT s, bool onDev)
{
    MK_CHECK2_SIZE(a_o, a_i);
    Vec3Di size = a_o.grid().size();
    
    size_t id = 0;
    for (int k=0; k < size.z; ++k)
        for (int j=0; j < size.y; ++j)
            for (int i=0; i < size.x; ++i, ++id) {
                float x = i + t.x;
                float y = j + t.y;
                float z = k + t.z;

                float ax, ay, az;
                triLerp<bg>(ax, ay, az,
                            a_i.x, a_i.y, a_i.z,
                            x, y, z,
                            size.x, size.y, size.z);

                a_o.x[id] = ax;
                a_o.y[id] = ay;
                a_o.z[id] = az;
            }
}

/*
 * Private field :: field function to call from other HField or VField to have the
 * syntax check
 */
 
template<BackgroundStrategy bg>
void CFieldOpers::
ApplyH(Field3D& a_o, const Field3D& a_i, const Field3D& h, StreamT st) {
    Vec3Di size = h.grid().size();
    size_t id = 0;
    for (int k=0; k < size.z; ++k)
        for (int j=0; j < size.y; ++j)
            for (int i=0; i < size.x; ++i, ++id) {
                float x = h.x[id], y = h.y[id], z = h.z[id];
                
                float ax, ay, az;
                triLerp<bg>(ax, ay, az,
                            a_i.x, a_i.y, a_i.z,
                            x, y, z,
                            size.x, size.y, size.z);
                
                a_o.x[id] = ax;
                a_o.y[id] = ay;
                a_o.z[id] = az;
            }
}

template<bool fwd, BackgroundStrategy bg>
void ApplyV(Field3D& a_o, const Field3D& a_i, const Field3D& u,
            float delta, StreamT s){

    MK_CHECK3_SIZE(a_o, a_i, u);
        
    Vec3Di size = a_o.grid().size();
    Vec3Df sp   = a_o.grid().spacing();
    Vec3Df iSp  = Vec3Df(1.f/sp.x, 1.f/sp.y, 1.f/sp.z);

    size_t id = 0;
    for (int k=0; k < size.z; ++k)
        for (int j=0; j < size.y; ++j)
            for (int i=0; i < size.x; ++i, ++id) {

                float x,y,z;
                if (fwd) {
                    x = i + delta * iSp.x * u.x[id];
                    y = j + delta * iSp.y * u.y[id];
                    z = k + delta * iSp.z * u.z[id];
                } else {
                    x = i - delta * iSp.x * u.x[id];
                    y = j - delta * iSp.y * u.y[id];
                    z = k - delta * iSp.z * u.z[id];
                }

                float ax, ay, az;
                triLerp<bg>(ax, ay, az,
                            a_i.x, a_i.y, a_i.z,
                            x, y, z,
                            size.x, size.y, size.z);
                
                a_o.x[id] = ax;
                a_o.y[id] = ay;
                a_o.z[id] = az;
            }
}

template<BackgroundStrategy bg> 
void CFieldOpers
::ApplyV(Field3D& a_o, const Field3D& a_i, const Field3D& u, const float& delta, StreamT s, bool onDev)
{
    PyCA::ApplyV<true, bg>(a_o, a_i, u, delta, s);
}

template<BackgroundStrategy bg> 
void CFieldOpers::
ApplyVInv(Field3D& a_o, const Field3D& a_i, const Field3D& u, const float& delta, StreamT s, bool onDev)
{
    PyCA::ApplyV<false, bg>(a_o, a_i, u, delta, s);
}


void CFieldOpers::
SplatField(Field3D& a_o, const Field3D& a_i, 
		  const Field3D& h, StreamT st)
{
   Vec3Di size = h.grid().size();
   CImageFieldOpers::Splat(a_o.x, h, a_i.x, false, st);
   CImageFieldOpers::Splat(a_o.y, h, a_i.y, false, st);
   CImageFieldOpers::Splat(a_o.z, h, a_i.z, false, st);
}

void CFieldOpers::
SubVol(Field3D& a_o, const Field3D& a_i, 
       const Vec3Di& start, StreamT st)
{
   Vec3Di oSize = a_o.size();
   Vec3Di iSize = a_i.size();
   CImageOpers::SubVol(a_o.x, a_i.x, oSize, iSize, start, st);
   CImageOpers::SubVol(a_o.y, a_i.y, oSize, iSize, start, st);
   CImageOpers::SubVol(a_o.z, a_i.z, oSize, iSize, start, st);
}

void CFieldOpers::
SetSubVol_I(Field3D& a_o, const Field3D& a_i, 
	    const Vec3Di& start, StreamT st)
{
   Vec3Di oSize = a_o.size();
   Vec3Di iSize = a_i.size();
   CImageOpers::SetSubVol_I(a_o.x, a_i.x, oSize, iSize, start, st);
   CImageOpers::SetSubVol_I(a_o.y, a_i.y, oSize, iSize, start, st);
   CImageOpers::SetSubVol_I(a_o.z, a_i.z, oSize, iSize, start, st);
}

template<BackgroundStrategy bg,  bool rescaleVector>
void CFieldOpers::
Resample(Field3D& a_o, const Field3D& a_i, StreamT stream)
{
    if (a_o.size() == a_i.size()){
       Opers::Copy(a_o, a_i, stream);
       return;
    }

    Vec3Di oSize = a_o.size();
    Vec3Di iSize = a_i.size();
    
    size_t id = 0;

    float rX = (float)iSize.x/ (float)oSize.x;
    float rY = (float)iSize.y/ (float)oSize.y;
    float rZ = (float)iSize.z/ (float)oSize.z;

    for (int k=0; k < oSize.z; ++k) {
        float i_z =  (rZ - 1.f) / 2.f + k * rZ;
        for (int j=0; j < oSize.y; ++j) {
            float i_y =  (rY - 1.f) / 2.f + j * rY;
            for (int i=0; i < oSize.x; ++i, ++id) {
                float i_x =  (rX - 1.f) / 2.f + i * rX;
                
                float ox, oy, oz;
                triLerp<bg>(ox, oy, oz,
                            a_i.x, a_i.y, a_i.z,
                            i_x, i_y, i_z,
                            iSize.x, iSize.y, iSize.z);
                
                if (rescaleVector){
                    ox /= rX; oy /= rY; oz /= rZ;
                }

                a_o.x[id] = ox;
                a_o.y[id] = oy;
                a_o.z[id] = oz;
            }
        }
    }
}

void CFieldOpers::
ReprojectToUnitVec(Field3D& a_o, StreamT st)
{
   Vec3Di sz = a_o.size();
   unsigned int numElements = sz.prod();
   for (unsigned int i = 0; i < numElements; ++i){
      Vec3Df v = a_o.get(i);
      double l = v.length();
      if(l > 1.0) v /= l;
      a_o.set(i,v);
   }
   
}

void CFieldOpers::
NormalizeSafe(Field3D& a_o, const Field3D& a_i, const float& eps, 
	      StreamT st)
{
   Vec3Di sz = a_o.size();
   unsigned int numElements = sz.prod();
   for (unsigned int i = 0; i < numElements; ++i){
      Vec3Df v = a_i.get(i);
      double l = v.length();
      if(l > eps) 
	 v /= l;
      else
	 v.x = v.y = v.z = 0.f;
      a_o.set(i,v);
   }
   
}

void CFieldOpers::
NormalizeSafe_I(Field3D& a_o, const float& eps, 
	      StreamT st)
{
   NormalizeSafe(a_o, a_o, eps, st);
}

void CFieldOpers::
Shrink(Field3D& a_o, const Field3D& a_i, const float& eps, 
       StreamT st)
{
   Vec3Di sz = a_o.size();
   unsigned int numElements = sz.prod();
   for (unsigned int i = 0; i < numElements; ++i){
      Vec3Df v = a_i.get(i);
      double l = v.length();
      if(l > eps){
	  v *= (l-eps)/l;
      }else{
	  v.x = v.y = v.z = 0.f;
      }
      a_o.set(i,v);
   }
   
}

void CFieldOpers::
Shrink_I(Field3D& a_o, const float& eps, 
	 StreamT st)
{
   Shrink(a_o, a_o, eps, st);
}


template<BackgroundStrategy bg> 
void CFieldOpers::
FixedPointInverse(Field3D &ginv, const Field3D& g, unsigned int numIter, StreamT stream, bool onDev)
{ // iterates, computing ginv -> ginv + x - g(ginv(x))
    Vec3Di sz = g.size();
    float ghx, ghy, ghz;
    int id = 0;
    for (int k=0; k < sz.z; ++k)
        for (int j=0; j < sz.y; ++j)
            for (int i=0; i < sz.x; ++i, ++id)
                for (unsigned int iter=0;iter < numIter;++iter)
                {
                    triLerp<bg>(ghx, ghy, ghz,
                                g.x, g.y, g.z,
				ginv.x[id], ginv.y[id], ginv.z[id],	
                                sz.x, sz.y, sz.z);
                    ginv.x[id] += 0.1*(float(i) - ghx);
                    ginv.y[id] += 0.1*(float(j) - ghy);
                    ginv.z[id] += 0.1*(float(k) - ghz);
                }
}
/*TODO: check why this version fails
void CFieldOpers::UpdateInverse(Field3D &ginv0t1,Field3D &scratchV, const Field3D& ginv0t, const Field3D& w, unsigned int numIter, StreamT stream, bool onDev)
{ // iterates, computing g_{t+1,0} given g_{t,0} as g_{t+1,0} = g_{t,0}\circ(Id - w\circ g_{t+1,0})

  //initialize
  Vec3Di sz = ginv0t.size();
  size_t nEl = ginv0t.nVox();
  MemOpers<EXEC_CPU, float>::Copy(ginv0t1.get(), ginv0t.get(), nEl, stream);
  float trlrpx, trlrpy, trlrpz;
  float trlrpx1, trlrpy1, trlrpz1;

  for(unsigned int iter=0;iter<numIter;iter++){
    int id = 0;
    for (int k=0; k < sz.z; ++k){
      for (int j=0; j < sz.y; ++j){
	for (int i=0; i < sz.x; ++i, ++id){
	  // compute w\circ g_{t+1,0}
	  triLerp<BACKGROUND_STRATEGY_PARTIAL_ZERO>(trlrpx, trlrpy, trlrpz,
						    w.x, w.y, w.z,
						    ginv0t1.x[id], ginv0t1.y[id], ginv0t1.z[id],	
						    sz.x, sz.y, sz.z);    
	  // compute Id - w\circ g_{t+1,0}
	  trlrpx = i - trlrpx;
	  trlrpy = j - trlrpy;
	  trlrpz = k - trlrpz;

	  // compute g_{t+1,0} = g_{t,0}\circ(Id - w\circ g_{t+1,0})
	  triLerp<BACKGROUND_STRATEGY_PARTIAL_ID>(trlrpx1, trlrpy1,trlrpz1,	
						  ginv0t.x, ginv0t.y, ginv0t.z,
						  trlrpx, trlrpy,trlrpz,	
						  sz.x, sz.y, sz.z);    

	  scratchV.x[id] = trlrpx1;	  
	  scratchV.y[id] = trlrpy1;
	  scratchV.z[id] = trlrpz1;
	}
      }
    }
    
    MemOpers<EXEC_CPU, float>::Copy(ginv0t1.get(), scratchV.get(), nEl, stream);
  }  
    
}
*/

void CFieldOpers::UpdateInverse(Field3D &ginv0t1,Field3D &scratchV, const Field3D& ginv0t, const Field3D& w, unsigned int numIter, StreamT stream, bool onDev)
{ // iterates, computing g_{t+1,0} given g_{t,0} as g_{t+1,0} = g_{t,0}\circ(Id - w\circ g_{t+1,0})

  //initialize
  Vec3Di size = ginv0t.size();
  Opers::Copy(ginv0t1, ginv0t, stream);

  for(unsigned int iter=0;iter<numIter;iter++){    
    // compute w\circ g_{t+1,0}
    ApplyH<BACKGROUND_STRATEGY_PARTIAL_ZERO>(scratchV, w, ginv0t1, stream);

    // compute Id - w\circ g_{t+1,0}    
    size_t id = 0;
    for (int k=0; k < size.z; ++k)
      for (int j=0; j < size.y; ++j)
	for (int i=0; i < size.x; ++i, ++id) {
	  scratchV.x[id] = i - scratchV.x[id];
	  scratchV.y[id] = j - scratchV.y[id];
	  scratchV.z[id] = k - scratchV.z[id];
	}

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
void CFieldOpers::
Ad(Field3D& Z, const Field3D& g, const Field3D& X,
                          StreamT s,bool onDev)
{
    Vec3Di size = X.size();

    float Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz;
    float dgvx, dgvy, dgvz;
    int id = 0;
    Opers::SetMem(Z, 0.f, s, false);

    for (int k=0; k < size.z; ++k)
        for (int j=0; j < size.y; ++j)
            for (int i=0; i < size.x; ++i, ++id)
            {
                // Get Jacobian matrix
                jacobianPoint<float,DIFF_CENTRAL,BC_CLAMP>(g.x,g.y,g.z,i,j,k,Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz,size.x,size.y,size.z);

                if (size.z == 1)
                { // Special case for flat images
                    Jxz = Jyz = Jzx = Jzy = 0;
                    Jzz = 1;
                }

                // Compute determinant
                float det = Jxx*(Jyy*Jzz-Jyz*Jzy)
                           -Jxy*(Jyx*Jzz-Jyz*Jzx)
                           +Jxz*(Jyx*Jzy-Jyy*Jzx);

                // Multiply by det*Dg X
                dgvx = det*(Jxx*X.x[id] + Jxy*X.y[id] + Jxz*X.z[id]);
                dgvy = det*(Jyx*X.x[id] + Jyy*X.y[id] + Jyz*X.z[id]);
                dgvz = det*(Jzx*X.x[id] + Jzy*X.y[id] + Jzz*X.z[id]);
		
		// Splat each component (non-normalized)
		Splatting::splatPoint(g.x[id], g.y[id], g.z[id], dgvx, Z.x, size.x, size.y, size.z);
		Splatting::splatPoint(g.x[id], g.y[id], g.z[id], dgvy, Z.y, size.x, size.y, size.z);
		Splatting::splatPoint(g.x[id], g.y[id], g.z[id], dgvz, Z.z, size.x, size.y, size.z);    
		
            }
}
/*
 * infinitesimal adjoint action
 * Z = ad_X Y = DX Y - DY X
 */
void CFieldOpers::
AdInf(Field3D& Z, const Field3D& X, const Field3D& Y,
                          StreamT s,bool onDev)
{
    Vec3Di sz = Y.size();
    float Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz;
    int id = 0;
    for (int k=0; k < sz.z; ++k)
        for (int j=0; j < sz.y; ++j)
            for (int i=0; i < sz.x; ++i, ++id)
            {
                // Get DX
	       jacobianPoint<float,DIFF_CENTRAL,BC_CLAMP>(X.x,X.y,X.z,i,j,k,Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz,sz.x,sz.y,sz.z);

                // Start with DX Y
                Z.x[id] = Jxx*Y.x[id] + Jxy*Y.y[id] + Jxz*Y.z[id];
                Z.y[id] = Jyx*Y.x[id] + Jyy*Y.y[id] + Jyz*Y.z[id];
                Z.z[id] = Jzx*Y.x[id] + Jzy*Y.y[id] + Jzz*Y.z[id];

                // Get DY
                jacobianPoint<float,DIFF_CENTRAL,BC_CLAMP>(Y.x,Y.y,Y.z,i,j,k,Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz,sz.x,sz.y,sz.z);

                // Subtract DY X
                Z.x[id] -= Jxx*X.x[id] + Jxy*X.y[id] + Jxz*X.z[id];
                Z.y[id] -= Jyx*X.x[id] + Jyy*X.y[id] + Jyz*X.z[id];
                Z.z[id] -= Jzx*X.x[id] + Jzy*X.y[id] + Jzz*X.z[id];
            }
}
/*
 * Coadjoint action of Diff on its Lie algebra
 * n = Ad_g^* m = (Dg)^T m\circ g |Dg|
 */
template<BackgroundStrategy bg> 
void CFieldOpers::
CoAd(Field3D& n, const Field3D& g, const Field3D& m,
                          StreamT s,bool onDev)
{
    Vec3Di sz = m.size();
    float Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz;
    float mgx,mgy,mgz;
    int id = 0;
    for (int k=0; k < sz.z; ++k)
        for (int j=0; j < sz.y; ++j)
            for (int i=0; i < sz.x; ++i, ++id)
            {
                // Get Jacobian matrix
                jacobianPoint<float,DIFF_CENTRAL,BC_CLAMP>(g.x,g.y,g.z,i,j,k,Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz,sz.x,sz.y,sz.z);

                if (sz.z == 1)
                { // Special case for flat images
                    Jxz = Jyz = Jzx = Jzy = 0;
                    Jzz = 1;
                }

                // Compute determinant
                float det = Jxx*(Jyy*Jzz-Jyz*Jzy)
                           -Jxy*(Jyx*Jzz-Jyz*Jzx)
                           +Jxz*(Jyx*Jzy-Jyy*Jzx);

                // Interpolate m
                triLerp<bg>(mgx, mgy, mgz,
                            m.x, m.y, m.z,
                            g.x[id], g.y[id], g.z[id],
                            sz.x, sz.y, sz.z);

                // Multiply by det*Dg^T
                n.x[id] = det*(Jxx*mgx + Jyx*mgy + Jzx*mgz);
                n.y[id] = det*(Jxy*mgx + Jyy*mgy + Jzy*mgz);
                n.z[id] = det*(Jxz*mgx + Jyz*mgy + Jzz*mgz);
            }
}
/*
 * infinitesimal coadjoint action
 * n = ad_X^* m = (DX)^T m + div(m \otimes X)
 */
void CFieldOpers::
CoAdInf(Field3D& n, const Field3D& X, const Field3D& m,
                          StreamT s,bool onDev)
{
    Vec3Di sz = m.size();
    float Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz;
    int id = 0;
    for (int k=0; k < sz.z; ++k)
        for (int j=0; j < sz.y; ++j)
            for (int i=0; i < sz.x; ++i, ++id)
            {
                // Start with the tensor product divergence piece
                divtensorprodPoint<float,DIFF_CENTRAL,BC_CLAMP>(m.x,m.y,m.z,X.x,X.y,X.z,i,j,k,n.x[id],n.y[id],n.z[id],sz.x,sz.y,sz.z);

                // Get DX
                jacobianPoint<float,DIFF_CENTRAL,BC_CLAMP>(X.x,X.y,X.z,i,j,k,Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz,sz.x,sz.y,sz.z);

                // Add the DX^T m term
                n.x[id] += Jxx*m.x[id] + Jyx*m.y[id] + Jzx*m.z[id];
                n.y[id] += Jxy*m.x[id] + Jyy*m.y[id] + Jzy*m.z[id];
                n.z[id] += Jxz*m.x[id] + Jyz*m.y[id] + Jzz*m.z[id];
            }
}

/*
 * computes tensor divergence of outer product of two vector fields 
 * Z = div(X \otimes Y)
 */
void CFieldOpers::
DivergenceTensor(Field3D& Z, const Field3D& X, const Field3D& Y,
                          StreamT s, bool onDev)
{
    Vec3Di sz = Y.size();
    int id = 0;
    for (int k=0; k < sz.z; ++k)
        for (int j=0; j < sz.y; ++j)
            for (int i=0; i < sz.x; ++i, ++id)
            {
                // The tensor product divergence 
                divtensorprodPoint<float,DIFF_CENTRAL,BC_CLAMP>(X.x,X.y,X.z,Y.x,Y.y,Y.z,i,j,k,Z.x[id],Z.y[id],Z.z[id],sz.x,sz.y,sz.z);
            }
}


/*
 * Jacobian X times Y
 * Z = DX Y
 */
void CFieldOpers::
JacobianXY(Field3D& Z, const Field3D& X, const Field3D& Y,
                          StreamT s,bool onDev)
{
    Vec3Di sz = Y.size();
    float Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz;
    int id = 0;
    for (int k=0; k < sz.z; ++k)
        for (int j=0; j < sz.y; ++j)
            for (int i=0; i < sz.x; ++i, ++id)
            {
                // Get DX
	       jacobianPoint<float,DIFF_CENTRAL,BC_CLAMP>(X.x,X.y,X.z,i,j,k,Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz,sz.x,sz.y,sz.z);

                // Start with DX Y
                Z.x[id] = Jxx*Y.x[id] + Jxy*Y.y[id] + Jxz*Y.z[id];
                Z.y[id] = Jyx*Y.x[id] + Jyy*Y.y[id] + Jyz*Y.z[id];
                Z.z[id] = Jzx*Y.x[id] + Jzy*Y.y[id] + Jzz*Y.z[id];
            }
}


/*
 * Jacobian X times Y
 * Z = (DX)' Y
 */
void CFieldOpers::
JacobianXtY(Field3D& Z, const Field3D& X, const Field3D& Y,
                          StreamT s,bool onDev)
{
    Vec3Di sz = Y.size();
    float Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz;
    int id = 0;
    for (int k=0; k < sz.z; ++k)
        for (int j=0; j < sz.y; ++j)
            for (int i=0; i < sz.x; ++i, ++id)
            {
                // Get DX
	       jacobianPoint<float,DIFF_CENTRAL,BC_CLAMP>(X.x,X.y,X.z,i,j,k,Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz,sz.x,sz.y,sz.z);

                // Compute (DX)' Y
                Z.x[id] = Jxx*Y.x[id] + Jyx*Y.y[id] + Jzx*Y.z[id];
                Z.y[id] = Jxy*Y.x[id] + Jyy*Y.y[id] + Jzy*Y.z[id];
                Z.z[id] = Jxz*Y.x[id] + Jyz*Y.y[id] + Jzz*Y.z[id];
            }
}

// Instantiations
template void CFieldOpers::
ComposeHH<BACKGROUND_STRATEGY_ID>(Field3D&, const Field3D&, const Field3D&, StreamT);
template void CFieldOpers::
ComposeHH<BACKGROUND_STRATEGY_PARTIAL_ID>(Field3D&, const Field3D&, const Field3D&, StreamT);
template void CFieldOpers::
ComposeHH<BACKGROUND_STRATEGY_CLAMP>(Field3D&, const Field3D&, const Field3D&, StreamT);

template void CFieldOpers::
ComposeVH<BACKGROUND_STRATEGY_ZERO>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);
template void CFieldOpers::
ComposeVH<BACKGROUND_STRATEGY_PARTIAL_ZERO>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);
template void CFieldOpers::
ComposeVH<BACKGROUND_STRATEGY_WRAP>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);
template void CFieldOpers::
ComposeVH<BACKGROUND_STRATEGY_CLAMP>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);

template void CFieldOpers::
ComposeVInvH<BACKGROUND_STRATEGY_ZERO>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);
template void CFieldOpers::
ComposeVInvH<BACKGROUND_STRATEGY_PARTIAL_ZERO>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);
template void CFieldOpers::
ComposeVInvH<BACKGROUND_STRATEGY_WRAP>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);
template void CFieldOpers::
ComposeVInvH<BACKGROUND_STRATEGY_CLAMP>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);

template void CFieldOpers::
ComposeHV<BACKGROUND_STRATEGY_ID>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);
template void CFieldOpers::
ComposeHV<BACKGROUND_STRATEGY_PARTIAL_ID>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);
template void CFieldOpers::
ComposeHV<BACKGROUND_STRATEGY_CLAMP>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);

template void CFieldOpers::
ComposeHVInv<BACKGROUND_STRATEGY_ID>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);
template void CFieldOpers::
ComposeHVInv<BACKGROUND_STRATEGY_PARTIAL_ID>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);
template void CFieldOpers::
ComposeHVInv<BACKGROUND_STRATEGY_CLAMP>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);

template void CFieldOpers::
ComposeTranslation<BACKGROUND_STRATEGY_ZERO>(Field3D&, const Field3D&, const Vec3Df&, StreamT,bool);
template void CFieldOpers::
ComposeTranslation<BACKGROUND_STRATEGY_PARTIAL_ZERO>(Field3D&, const Field3D&, const Vec3Df&, StreamT,bool);
template void CFieldOpers::
ComposeTranslation<BACKGROUND_STRATEGY_ID>(Field3D&, const Field3D&, const Vec3Df&, StreamT,bool);
template void CFieldOpers::
ComposeTranslation<BACKGROUND_STRATEGY_PARTIAL_ID>(Field3D&, const Field3D&, const Vec3Df&, StreamT,bool);
template void CFieldOpers::
ComposeTranslation<BACKGROUND_STRATEGY_WRAP>(Field3D&, const Field3D&, const Vec3Df&, StreamT,bool);
template void CFieldOpers::
ComposeTranslation<BACKGROUND_STRATEGY_CLAMP>(Field3D&, const Field3D&, const Vec3Df&, StreamT,bool);

template void CFieldOpers::
ApplyH<BACKGROUND_STRATEGY_ZERO>(Field3D&, const Field3D&, const Field3D&, StreamT st);
template void CFieldOpers::
ApplyH<BACKGROUND_STRATEGY_PARTIAL_ZERO>(Field3D&, const Field3D&, const Field3D&, StreamT st);
template void CFieldOpers::
ApplyH<BACKGROUND_STRATEGY_ID>(Field3D&, const Field3D&, const Field3D&, StreamT st);
template void CFieldOpers::
ApplyH<BACKGROUND_STRATEGY_PARTIAL_ID>(Field3D&, const Field3D&, const Field3D&, StreamT st);
template void CFieldOpers::
ApplyH<BACKGROUND_STRATEGY_WRAP>(Field3D&, const Field3D&, const Field3D&, StreamT st);
template void CFieldOpers::
ApplyH<BACKGROUND_STRATEGY_CLAMP>(Field3D&, const Field3D&, const Field3D&, StreamT st);


template void CFieldOpers::
ApplyV<BACKGROUND_STRATEGY_ZERO>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT s,bool);
template void CFieldOpers::
ApplyV<BACKGROUND_STRATEGY_PARTIAL_ZERO>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT s,bool);
template void CFieldOpers::
ApplyV<BACKGROUND_STRATEGY_ID>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT s,bool);
template void CFieldOpers::
ApplyV<BACKGROUND_STRATEGY_PARTIAL_ID>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT s,bool);
template void CFieldOpers::
ApplyV<BACKGROUND_STRATEGY_WRAP>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT s,bool);
template void CFieldOpers::
ApplyV<BACKGROUND_STRATEGY_CLAMP>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT s,bool);

template void CFieldOpers::
ApplyVInv<BACKGROUND_STRATEGY_ZERO>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT s,bool);
template void CFieldOpers::
ApplyVInv<BACKGROUND_STRATEGY_PARTIAL_ZERO>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT s,bool);
template void CFieldOpers::
ApplyVInv<BACKGROUND_STRATEGY_ID>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT s,bool);
template void CFieldOpers::
ApplyVInv<BACKGROUND_STRATEGY_PARTIAL_ID>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT s,bool);
template void CFieldOpers::
ApplyVInv<BACKGROUND_STRATEGY_WRAP>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT s,bool);
template void CFieldOpers::
ApplyVInv<BACKGROUND_STRATEGY_CLAMP>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT s,bool);

template void CFieldOpers::
Resample<BACKGROUND_STRATEGY_ZERO, true>(Field3D&, const Field3D&, StreamT);
template void CFieldOpers::
Resample<BACKGROUND_STRATEGY_PARTIAL_ZERO,true>(Field3D&, const Field3D&, StreamT);
template void CFieldOpers::
Resample<BACKGROUND_STRATEGY_ID, true>(Field3D&, const Field3D&, StreamT);
template void CFieldOpers::
Resample<BACKGROUND_STRATEGY_PARTIAL_ID,true>(Field3D&, const Field3D&, StreamT);
template void CFieldOpers::
Resample<BACKGROUND_STRATEGY_CLAMP, true>(Field3D&, const Field3D&, StreamT);
template void CFieldOpers::
Resample<BACKGROUND_STRATEGY_WRAP, true>(Field3D&, const Field3D&, StreamT);
template void CFieldOpers::
Resample<BACKGROUND_STRATEGY_ZERO, false>(Field3D&, const Field3D&, StreamT);
template void CFieldOpers::
Resample<BACKGROUND_STRATEGY_PARTIAL_ZERO,false>(Field3D&, const Field3D&, StreamT);
template void CFieldOpers::
Resample<BACKGROUND_STRATEGY_ID, false>(Field3D&, const Field3D&, StreamT);
template void CFieldOpers::
Resample<BACKGROUND_STRATEGY_PARTIAL_ID,false>(Field3D&, const Field3D&, StreamT);
template void CFieldOpers::
Resample<BACKGROUND_STRATEGY_CLAMP, false>(Field3D&, const Field3D&, StreamT);
template void CFieldOpers::
Resample<BACKGROUND_STRATEGY_WRAP, false>(Field3D&, const Field3D&, StreamT);

template void CFieldOpers::
FixedPointInverse<BACKGROUND_STRATEGY_ID>(Field3D &ginv, const Field3D& g, unsigned int numIter, StreamT stream, bool onDev);
template void CFieldOpers::
FixedPointInverse<BACKGROUND_STRATEGY_PARTIAL_ID>(Field3D &ginv, const Field3D& g, unsigned int numIter, StreamT stream, bool onDev);
template void CFieldOpers::
FixedPointInverse<BACKGROUND_STRATEGY_WRAP>(Field3D &ginv, const Field3D& g, unsigned int numIter, StreamT stream, bool onDev);

// #ifndef _MSC_VER
// // Fix MSVC complaining about a redeclaration here. Not sure why it is -JDH Nov 2013
// void CFieldOpers::UpdateInverse(Field3D &ginv0t1,Field3D &scratchV, const Field3D& ginv0t, const Field3D& w, unsigned int numIter, StreamT stream, bool onDev);
// #endif

template void CFieldOpers::
Ad<BACKGROUND_STRATEGY_ZERO>(Field3D& Z, const Field3D& g, const Field3D& X, StreamT s,bool onDev);
template void CFieldOpers::
Ad<BACKGROUND_STRATEGY_PARTIAL_ZERO>(Field3D& Z, const Field3D& g, const Field3D& X, StreamT s,bool onDev);
template void CFieldOpers::
Ad<BACKGROUND_STRATEGY_WRAP>(Field3D& Z, const Field3D& g, const Field3D& X, StreamT s,bool onDev);
template void CFieldOpers::
Ad<BACKGROUND_STRATEGY_CLAMP>(Field3D& Z, const Field3D& g, const Field3D& X, StreamT s,bool onDev);
template void CFieldOpers::
CoAd<BACKGROUND_STRATEGY_ZERO>(Field3D& n, const Field3D& g, const Field3D& m, StreamT s,bool onDev);
template void CFieldOpers::
CoAd<BACKGROUND_STRATEGY_PARTIAL_ZERO>(Field3D& n, const Field3D& g, const Field3D& m, StreamT s,bool onDev);
template void CFieldOpers::
CoAd<BACKGROUND_STRATEGY_WRAP>(Field3D& n, const Field3D& g, const Field3D& m, StreamT s,bool onDev);
template void CFieldOpers::
CoAd<BACKGROUND_STRATEGY_CLAMP>(Field3D& n, const Field3D& g, const Field3D& m, StreamT s,bool onDev);
} // end namespace PyCA

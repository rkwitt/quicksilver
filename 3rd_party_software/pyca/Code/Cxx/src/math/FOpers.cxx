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

#include "FOpers.h"

#include "FieldOpers.h"
#include <Field3D.h>
#include <MemOpers.h>


#define FIELD_BINARY_ARRAY_OPC_VEC_DEF(OP)     	       	       	       	   \
void Opers::								   \
OP##C(Field3D& a_o, const Field3D& a_i, const Vec3Df& c,		   \
     StreamT st, bool onDev)					           \
{									   \
   MK_CHECK2_ALL(a_o, a_i);						   \
   AUTO_EXEC(a_o.memType(), FieldOpers,					   \
	     OP##C(a_o, a_i, c, st, onDev));				   \
}

#define FIELD_BINARY_ARRAY_OPC_FLOAT_DEF(OP)				   \
void Opers::								   \
OP##C(Field3D& a_o, const Field3D& a_i, const float& c,			   \
     StreamT st, bool onDev)					           \
{									   \
   MK_CHECK2_ALL(a_o, a_i);						   \
   AUTO_EXEC(a_o.memType(), FieldOpers,					   \
	     OP##C(a_o, a_i, c, st, onDev));				   \
}

#define FIELD_BINARY_ARRAY_OPC_I_VEC_DEF(OP)				   \
void Opers::								   \
OP##C_I(Field3D& a_o, const Vec3Df& c,					   \
       StreamT st, bool onDev)					           \
{									   \
   AUTO_EXEC(a_o.memType(), FieldOpers,					   \
	     OP##C_I(a_o, c, st, onDev));				   \
}

#define FIELD_BINARY_ARRAY_OPC_I_FLOAT_DEF(OP)				   \
void Opers::								   \
OP##C_I(Field3D& a_o, const float& c,					   \
       StreamT st, bool onDev)					           \
{									   \
   AUTO_EXEC(a_o.memType(), FieldOpers,					   \
	     OP##C_I(a_o, c, st, onDev));				   \
}

#define FIELD_BINARY_ARRAY_OP_FIELD_DEF(OP)				   \
void Opers::								   \
OP(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1,		   \
    StreamT st)							           \
{									   \
   MK_CHECK3_ALL(a_o, a_i, a_i1);					   \
   AUTO_EXEC(a_o.memType(), FieldOpers,					   \
	     OP(a_o, a_i, a_i1, st));					   \
}

#define FIELD_BINARY_ARRAY_OP_I_FIELD_DEF(OP)				   \
void Opers::								   \
OP##_I(Field3D& a_o, const Field3D& a_i,				   \
      StreamT st)							   \
{									   \
   MK_CHECK2_ALL(a_o, a_i);						   \
   AUTO_EXEC(a_o.memType(), FieldOpers,					   \
	     OP##_I(a_o, a_i, st));					   \
}

#define FIELD_BINARY_ARRAY_OP_DEFS(OP)					   \
FIELD_BINARY_ARRAY_OPC_VEC_DEF(OP)					   \
FIELD_BINARY_ARRAY_OPC_FLOAT_DEF(OP)					   \
FIELD_BINARY_ARRAY_OPC_I_VEC_DEF(OP)					   \
FIELD_BINARY_ARRAY_OPC_I_FLOAT_DEF(OP)					   \
FIELD_BINARY_ARRAY_OP_FIELD_DEF(OP)					   \
FIELD_BINARY_ARRAY_OP_I_FIELD_DEF(OP)

namespace PyCA {

void Opers::
Copy(Field3D& a_o, const Field3D& a_i, StreamT stream){
    MK_CHECK2_SIZE(a_o, a_i);
    Copy(a_o.getMemPoolX(), a_i.getMemPoolX(), a_o.nVox(), stream);
    Copy(a_o.getMemPoolY(), a_i.getMemPoolY(), a_o.nVox(), stream);
    Copy(a_o.getMemPoolZ(), a_i.getMemPoolZ(), a_o.nVox(), stream);
}

////////////////////////////////////////////////////////////////////////////
// compose two hfields to get an hfield
// f(x) = g(h(x))
////////////////////////////////////////////////////////////////////////////

template<BackgroundStrategy bg>
void Opers::
ComposeHH(Field3D& f, const Field3D& g, const Field3D& h, StreamT stream)
{
   MK_CHECK3_ALL(f, g, h);
   AUTO_EXEC(h.memType(), FieldOpers, template ComposeHH<bg>(f, g, h, stream));
}

// non-template version
void Opers::
ComposeHH(Field3D& f, const Field3D& g,
	  const Field3D& h,
	  BackgroundStrategy bg,
	  StreamT s)
{
   if(bg == BACKGROUND_STRATEGY_ID){
      ComposeHH<BACKGROUND_STRATEGY_ID>(f, g, h, s);
   }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ID){
      ComposeHH<BACKGROUND_STRATEGY_PARTIAL_ID>(f, g, h, s);
   }else if(bg == BACKGROUND_STRATEGY_CLAMP){
      ComposeHH<BACKGROUND_STRATEGY_CLAMP>(f, g, h, s);
   }else{
      throw PyCAException(__FILE__, __LINE__, "Unknown/unsupported background strategy");
   }
}

////////////////////////////////////////////////////////////////////////////
// compose a velocity and hfield to get an hfield
// h(x) = g(x) + delta * v(g(x))
////////////////////////////////////////////////////////////////////////////

template<BackgroundStrategy bg>
void Opers::
ComposeVH(Field3D& h, const Field3D& v, const Field3D& g, const float& delta,
          StreamT stream, bool onDev) {
   MK_CHECK3_ALL(h, v, g);
   AUTO_EXEC(h.memType(), FieldOpers,
	     template ComposeVH<bg>(h, v, g, delta, stream, onDev));
}

// non-template version
void Opers::
ComposeVH(Field3D& h, const Field3D& v,
	  const Field3D& g, const float& delta,
	  BackgroundStrategy bg,
	  StreamT s,bool onDev)
{
   if(bg == BACKGROUND_STRATEGY_ZERO){
      ComposeVH<BACKGROUND_STRATEGY_ZERO>(h, v, g, delta, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ZERO){
      ComposeVH<BACKGROUND_STRATEGY_PARTIAL_ZERO>(h, v, g, delta, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_CLAMP){
      ComposeVH<BACKGROUND_STRATEGY_CLAMP>(h, v, g, delta, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_WRAP){
      ComposeVH<BACKGROUND_STRATEGY_WRAP>(h, v, g, delta, s, onDev);
   }else{
      throw PyCAException(__FILE__, __LINE__, "Unknown/unsupported background strategy");
   }
}

////////////////////////////////////////////////////////////////////////////
// compose an inverse velocity and hfield to get an hfield
// h(x) = g(x) - delta * v(g(x))
////////////////////////////////////////////////////////////////////////////

template<BackgroundStrategy bg>
void Opers::
ComposeVInvH(Field3D& h, const Field3D& v, const Field3D& g, const float& delta, StreamT stream, bool onDev) {
   MK_CHECK3_ALL(h, v, g);
   AUTO_EXEC(h.memType(), FieldOpers,
	     template ComposeVInvH<bg>(h, v, g, delta, stream, onDev));
}

// non-template version
void Opers::
ComposeVInvH(Field3D& h, const Field3D& v,
	  const Field3D& g, const float& delta,
	  BackgroundStrategy bg,
	  StreamT s,bool onDev)
{
   if(bg == BACKGROUND_STRATEGY_ZERO){
      ComposeVInvH<BACKGROUND_STRATEGY_ZERO>(h, v, g, delta, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ZERO){
      ComposeVInvH<BACKGROUND_STRATEGY_PARTIAL_ZERO>(h, v, g, delta, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_CLAMP){
      ComposeVInvH<BACKGROUND_STRATEGY_CLAMP>(h, v, g, delta, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_WRAP){
      ComposeVInvH<BACKGROUND_STRATEGY_WRAP>(h, v, g, delta, s, onDev);
   }else{
      throw PyCAException(__FILE__, __LINE__, "Unknown/unsupported background strategy");
   }
}

/**
 * compose hfield and vfield to get hfield
 * h(x) = g(x+ delta * v(x))
 */

template<BackgroundStrategy bg>
void Opers::
ComposeHV(Field3D& h, const Field3D& g, const Field3D& v, const float& delta, StreamT stream, bool onDev){
   MK_CHECK3_ALL(h, g, v);
   AUTO_EXEC(h.memType(), FieldOpers,
	     template ComposeHV<bg>(h, g, v, delta, stream, onDev));
}

// non-template version
void Opers::
ComposeHV(Field3D& h, const Field3D& g,
	  const Field3D& v, const float& delta,
	  BackgroundStrategy bg,
	  StreamT s,bool onDev)
{
   if(bg == BACKGROUND_STRATEGY_ID){
      ComposeHV<BACKGROUND_STRATEGY_ID>(h, g, v, delta, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ID){
      ComposeHV<BACKGROUND_STRATEGY_PARTIAL_ID>(h, g, v, delta, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_CLAMP){
      ComposeHV<BACKGROUND_STRATEGY_CLAMP>(h, g, v, delta, s, onDev);
   }else{
      throw PyCAException(__FILE__, __LINE__, "Unknown/unsupported background strategy");
   }
}

template<BackgroundStrategy bg>
void Opers::
ComposeHVInv(Field3D& h, const Field3D& g, const Field3D& v, const float& delta, StreamT stream, bool onDev){
   MK_CHECK3_ALL(h, g, v);
   AUTO_EXEC(h.memType(), FieldOpers,
	     template ComposeHVInv<bg>(h, g, v, delta, stream, onDev));
}

// non-template version
void Opers::
ComposeHVInv(Field3D& h, const Field3D& g,
	  const Field3D& v, const float& delta,
	  BackgroundStrategy bg,
	  StreamT s,bool onDev)
{
   if(bg == BACKGROUND_STRATEGY_ID){
      ComposeHVInv<BACKGROUND_STRATEGY_ID>(h, g, v, delta, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ID){
      ComposeHVInv<BACKGROUND_STRATEGY_PARTIAL_ID>(h, g, v, delta, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_CLAMP){
      ComposeHVInv<BACKGROUND_STRATEGY_CLAMP>(h, g, v, delta, s, onDev);
   }else{
      throw PyCAException(__FILE__, __LINE__, "Unknown/unsupported background strategy");
   }
}


////////////////////////////////////////////////////////////////////////////
// Compose field with translation
// creating a_o(x) = a_i(x + t)
////////////////////////////////////////////////////////////////////////////

template<BackgroundStrategy bg>
void Opers::
ComposeTranslation(Field3D& a_o, const Field3D& a_i, const Vec3Df& t, StreamT stream, bool onDev){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     template ComposeTranslation<bg>(a_o, a_i, t, stream, onDev));
}

// non-template version
void Opers::
ComposeTranslation(Field3D& a_o, const Field3D& a_i, const Vec3Df& t,
		   BackgroundStrategy bg, StreamT s, bool onDev)
{
   if(bg == BACKGROUND_STRATEGY_ID){
      ComposeTranslation<BACKGROUND_STRATEGY_ID>(a_o, a_i, t, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ID){
      ComposeTranslation<BACKGROUND_STRATEGY_PARTIAL_ID>(a_o, a_i, t, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_ZERO){
      ComposeTranslation<BACKGROUND_STRATEGY_ZERO>(a_o, a_i, t, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ZERO){
      ComposeTranslation<BACKGROUND_STRATEGY_PARTIAL_ZERO>(a_o, a_i, t, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_CLAMP){
      ComposeTranslation<BACKGROUND_STRATEGY_CLAMP>(a_o, a_i, t, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_WRAP){
      ComposeTranslation<BACKGROUND_STRATEGY_WRAP>(a_o, a_i, t, s, onDev);
   }else{
      throw PyCAException(__FILE__, __LINE__, "Unknown/unsupported background strategy");
   }
}


template<BackgroundStrategy bg>
void Opers::
ApplyH(Field3D& a_o, const Field3D& a_i, const Field3D& h, StreamT stream) {
   MK_CHECK3_ALL(a_o, a_i, h);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     template ApplyH<bg>(a_o, a_i, h, stream));
}

// non-template version
void Opers::
ApplyH(Field3D& a_o, const Field3D& a_i, const Field3D& h,
       BackgroundStrategy bg, StreamT s)
{
   if(bg == BACKGROUND_STRATEGY_ID){
      ApplyH<BACKGROUND_STRATEGY_ID>(a_o, a_i, h, s);
   }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ID){
      ApplyH<BACKGROUND_STRATEGY_PARTIAL_ID>(a_o, a_i, h, s);
   }else if(bg == BACKGROUND_STRATEGY_ZERO){
      ApplyH<BACKGROUND_STRATEGY_ZERO>(a_o, a_i, h, s);
   }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ZERO){
      ApplyH<BACKGROUND_STRATEGY_PARTIAL_ZERO>(a_o, a_i, h, s);
   }else if(bg == BACKGROUND_STRATEGY_CLAMP){
      ApplyH<BACKGROUND_STRATEGY_CLAMP>(a_o, a_i, h, s);
   }else if(bg == BACKGROUND_STRATEGY_WRAP){
      ApplyH<BACKGROUND_STRATEGY_WRAP>(a_o, a_i, h, s);
   }else{
      throw PyCAException(__FILE__, __LINE__, "Unknown/unsupported background strategy");
   }
}


template<BackgroundStrategy bg>
void Opers::
ApplyV(Field3D& a_o, const Field3D& a_i, const Field3D& u, const float& delta, StreamT stream, bool onDev) {
   MK_CHECK3_ALL(a_o, a_i, u);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     template ApplyV<bg>(a_o, a_i, u, delta, stream, onDev));
}

// nont-template version
void Opers::
ApplyV(Field3D& a_o, const Field3D& a_i, const Field3D& u,
       BackgroundStrategy bg,
       const float& delta, StreamT s, bool onDev)
{
   if(bg == BACKGROUND_STRATEGY_ID){
      ApplyV<BACKGROUND_STRATEGY_ID>(a_o, a_i, u, delta, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ID){
      ApplyV<BACKGROUND_STRATEGY_PARTIAL_ID>(a_o, a_i, u, delta, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_ZERO){
      ApplyV<BACKGROUND_STRATEGY_ZERO>(a_o, a_i, u, delta, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ZERO){
      ApplyV<BACKGROUND_STRATEGY_PARTIAL_ZERO>(a_o, a_i, u, delta, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_CLAMP){
      ApplyV<BACKGROUND_STRATEGY_CLAMP>(a_o, a_i, u, delta, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_WRAP){
      ApplyV<BACKGROUND_STRATEGY_WRAP>(a_o, a_i, u, delta, s, onDev);
   }else{
      throw PyCAException(__FILE__, __LINE__, "Unknown/unsupported background strategy");
   }
}

template<BackgroundStrategy bg>
void Opers::
ApplyVInv(Field3D& a_o, const Field3D& a_i, const Field3D& u, const float& delta, StreamT stream, bool onDev) {
   MK_CHECK3_ALL(a_o, a_i, u);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     template ApplyVInv<bg>(a_o, a_i, u, delta, stream, onDev));
}

// non-template version
void Opers::
ApplyVInv(Field3D& a_o, const Field3D& a_i, const Field3D& u,
	  BackgroundStrategy bg,
	  const float& delta, StreamT s, bool onDev)
{
   if(bg == BACKGROUND_STRATEGY_ID){
      ApplyVInv<BACKGROUND_STRATEGY_ID>(a_o, a_i, u, delta, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ID){
      ApplyVInv<BACKGROUND_STRATEGY_PARTIAL_ID>(a_o, a_i, u, delta, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_ZERO){
      ApplyVInv<BACKGROUND_STRATEGY_ZERO>(a_o, a_i, u, delta, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ZERO){
      ApplyVInv<BACKGROUND_STRATEGY_PARTIAL_ZERO>(a_o, a_i, u, delta, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_CLAMP){
      ApplyVInv<BACKGROUND_STRATEGY_CLAMP>(a_o, a_i, u, delta, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_WRAP){
      ApplyVInv<BACKGROUND_STRATEGY_WRAP>(a_o, a_i, u, delta, s, onDev);
   }else{
      throw PyCAException(__FILE__, __LINE__, "Unknown/unsupported background strategy");
   }
}


void Opers::
Jacobian(Field3D& a_Xg, Field3D& a_Yg, Field3D& a_Zg,
         const Field3D& a_h,
	 DiffT diffType, BoundaryCondT bc, StreamT stream)
{
   MK_CHECK4_ALL(a_Xg, a_Yg, a_Zg, a_h);
   AUTO_EXEC(a_Xg.memType(), FieldOpers,
	     Jacobian(a_Xg, a_Yg, a_Zg, a_h, diffType, bc, stream));
}


template<BackgroundStrategy bg, bool rescaleVector>
void Opers::
Resample(Field3D& a_o, const Field3D& a_i, StreamT stream){
   MK_CHECK2_MEM(a_o, a_i);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     template Resample<bg COMMA rescaleVector>(a_o, a_i, stream));
}

// non-template version
void Opers::
Resample(Field3D& a_o, const Field3D& a_i,
	 BackgroundStrategy bg, bool rescaleVector,
	 StreamT stream)
{
   if(rescaleVector){
      if(bg == BACKGROUND_STRATEGY_ID){
	 Resample<BACKGROUND_STRATEGY_ID, true>(a_o, a_i, stream);
      }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ID){
	 Resample<BACKGROUND_STRATEGY_PARTIAL_ID, true>(a_o, a_i, stream);
      }else if(bg == BACKGROUND_STRATEGY_ZERO){
	 Resample<BACKGROUND_STRATEGY_ZERO, true>(a_o, a_i, stream);
      }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ZERO){
	 Resample<BACKGROUND_STRATEGY_PARTIAL_ZERO, true>(a_o, a_i, stream);
      }else if(bg == BACKGROUND_STRATEGY_CLAMP){
	 Resample<BACKGROUND_STRATEGY_CLAMP, true>(a_o, a_i, stream);
      }else if(bg == BACKGROUND_STRATEGY_WRAP){
	 Resample<BACKGROUND_STRATEGY_WRAP, true>(a_o, a_i, stream);
      }else{
	 throw PyCAException(__FILE__, __LINE__, "Unknown/unsupported background strategy");
      }
   }else{

      if(bg == BACKGROUND_STRATEGY_ID){
	 Resample<BACKGROUND_STRATEGY_ID, false>(a_o, a_i, stream);
      }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ID){
	 Resample<BACKGROUND_STRATEGY_PARTIAL_ID, false>(a_o, a_i, stream);
      }else if(bg == BACKGROUND_STRATEGY_ZERO){
	 Resample<BACKGROUND_STRATEGY_ZERO, false>(a_o, a_i, stream);
      }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ZERO){
	 Resample<BACKGROUND_STRATEGY_PARTIAL_ZERO, false>(a_o, a_i, stream);
      }else if(bg == BACKGROUND_STRATEGY_CLAMP){
	 Resample<BACKGROUND_STRATEGY_CLAMP, false>(a_o, a_i, stream);
      }else if(bg == BACKGROUND_STRATEGY_WRAP){
	 Resample<BACKGROUND_STRATEGY_WRAP, false>(a_o, a_i, stream);
      }else{
	 throw PyCAException(__FILE__, __LINE__, "Unknown/unsupported background strategy");
      }
   }
}

void Opers::
SplatField(Field3D& a_o, const Field3D& a_i, const Field3D& h, StreamT stream) {
   MK_CHECK3_ALL(a_o, a_i, h);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     SplatField(a_o, a_i, h, stream));
}

void Opers::
SubVol(Field3D& a_o, const Field3D& a_i,
       const Vec3Di& start, StreamT stream)
{
   MK_CHECK2_MEM(a_o, a_i);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     SubVol(a_o, a_i, start, stream));
}

void Opers::
SetSubVol_I(Field3D& a_o, const Field3D& a_i,
	    const Vec3Di& start, StreamT stream)
{
   MK_CHECK2_MEM(a_o, a_i);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     SetSubVol_I(a_o, a_i, start, stream));
}

void Opers::
ReprojectToUnitVec(Field3D& a_o, StreamT st)
{
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     ReprojectToUnitVec(a_o, st));
}

void Opers::
NormalizeSafe(Field3D& a_o, const Field3D& a_i, const float& eps,
	      StreamT st)
{
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     NormalizeSafe(a_o, a_i, eps, st));
}

void Opers::
NormalizeSafe_I(Field3D& a_o, const float& eps, StreamT st)
{
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     NormalizeSafe_I(a_o, eps, st));
}

void Opers::
Shrink(Field3D& a_o, const Field3D& a_i, const float& eps,
	      StreamT st)
{
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     Shrink(a_o, a_i, eps, st));
}

void Opers::
Shrink_I(Field3D& a_o, const float& eps, StreamT st)
{
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     Shrink_I(a_o, eps, st));
}

template<BackgroundStrategy bg>
void Opers::
FixedPointInverse(Field3D &ginv, const Field3D& g, unsigned int numIter, StreamT stream, bool onDev)
{
   MK_CHECK2_ALL(ginv, g);
   AUTO_EXEC(ginv.memType(), FieldOpers,
	     template FixedPointInverse<bg>(ginv,g,numIter,stream,onDev));
}

void Opers::
FixedPointInverse(Field3D &ginv, const Field3D& g, unsigned int numIter, StreamT stream, bool onDev)
{
    FixedPointInverse<BACKGROUND_STRATEGY_PARTIAL_ID>(ginv,g,numIter,stream,onDev);
}

void Opers::
UpdateInverse(Field3D &ginv0t1,Field3D &scratchV, const Field3D& ginv0t, const Field3D& w, unsigned int numIter, StreamT stream, bool onDev){
  MK_CHECK4_ALL(ginv0t1,scratchV, ginv0t, w);
  AUTO_EXEC(ginv0t.memType(), FieldOpers, 
	    UpdateInverse(ginv0t1,scratchV,ginv0t,w,numIter,stream,onDev));
}

template<BackgroundStrategy bg>
void Opers::
Ad(Field3D& Z, const Field3D& g, const Field3D& X, StreamT st, bool onDev)
{
   MK_CHECK3_ALL(Z, g, X);
   AUTO_EXEC(Z.memType(), FieldOpers,
	     template Ad<bg>(Z,g,X,st,onDev));
}
// non-template version
void Opers::
Ad(Field3D& Z, const Field3D& g, const Field3D& X, StreamT st, BackgroundStrategy bg, bool onDev)
{
  if (bg==BACKGROUND_STRATEGY_ZERO) Ad<BACKGROUND_STRATEGY_ZERO>(Z,g,X,st,onDev);
    else if (bg==BACKGROUND_STRATEGY_PARTIAL_ZERO) Ad<BACKGROUND_STRATEGY_PARTIAL_ZERO>(Z,g,X,st,onDev);
    else if (bg==BACKGROUND_STRATEGY_WRAP) Ad<BACKGROUND_STRATEGY_WRAP>(Z,g,X,st,onDev);
    else if (bg==BACKGROUND_STRATEGY_CLAMP) Ad<BACKGROUND_STRATEGY_CLAMP>(Z,g,X,st,onDev);
    else throw PyCAException(__FILE__, __LINE__, "Unknown/unsupported background strategy");
}
void Opers::
AdInf(Field3D& Z, const Field3D& X, const Field3D& Y, StreamT st, bool onDev)
{
   MK_CHECK3_ALL(Z, X, Y);
   AUTO_EXEC(Z.memType(), FieldOpers,
	     AdInf(Z,X,Y,st,onDev));
}
template<BackgroundStrategy bg>
void Opers::
CoAd(Field3D& n, const Field3D& g, const Field3D& m, StreamT st, bool onDev)
{
   MK_CHECK3_ALL(n, g, m);
   AUTO_EXEC(n.memType(), FieldOpers,
	     template CoAd<bg>(n,g,m,st,onDev));
}
// non-template version
void Opers::
CoAd(Field3D& n, const Field3D& g, const Field3D& m, StreamT st, BackgroundStrategy bg, bool onDev)
{
    if (bg==BACKGROUND_STRATEGY_ZERO) CoAd<BACKGROUND_STRATEGY_ZERO>(n,g,m,st,onDev);
    else if (bg==BACKGROUND_STRATEGY_PARTIAL_ZERO) CoAd<BACKGROUND_STRATEGY_PARTIAL_ZERO>(n,g,m,st,onDev);
    else if (bg==BACKGROUND_STRATEGY_WRAP) CoAd<BACKGROUND_STRATEGY_WRAP>(n,g,m,st,onDev);
    else if (bg==BACKGROUND_STRATEGY_CLAMP) CoAd<BACKGROUND_STRATEGY_CLAMP>(n,g,m,st,onDev);
    else throw PyCAException(__FILE__, __LINE__, "Unknown/unsupported background strategy");
}
void Opers::
CoAdInf(Field3D& n, const Field3D& X, const Field3D& m, StreamT st, bool onDev)
{
   MK_CHECK3_ALL(n, X, m);
   AUTO_EXEC(n.memType(), FieldOpers,
	     CoAdInf(n,X,m,st,onDev));
}
void Opers::
DivergenceTensor(Field3D& Z, const Field3D& X, const Field3D& Y, StreamT st, bool onDev)
{
   MK_CHECK3_ALL(Z, X, Y);
   AUTO_EXEC(Z.memType(), FieldOpers,
	     DivergenceTensor(Z,X,Y,st,onDev));
}

void Opers::
JacobianXY(Field3D& Z, const Field3D& X, const Field3D& Y, StreamT st, bool onDev)
{
   MK_CHECK3_ALL(Z, X, Y);
   AUTO_EXEC(Z.memType(), FieldOpers, 
	     JacobianXY(Z,X,Y,st,onDev));
}

void Opers::
JacobianXtY(Field3D& Z, const Field3D& X, const Field3D& Y, StreamT st, bool onDev)
{
   MK_CHECK3_ALL(Z, X, Y);
   AUTO_EXEC(Z.memType(), FieldOpers, 
	     JacobianXtY(Z,X,Y,st,onDev));
}

/*
 * Set all the value of the Vector field with single value Vector3f(v)
 */
void Opers::
SetMem(Field3D& a_o, const Vec3Df& v, StreamT st,bool onDev){
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     SetMem(a_o, v, st, onDev));
}
void Opers::
SetMem(Field3D& a_o, const Vec3Df& v, Image3D& mask, StreamT st,bool onDev){
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     SetMem(a_o, v, mask, st, onDev));
}
/**
 * Set all the value of the Vector field with single float value
 * (normaly used to zero out the memory)
 */
void Opers::
SetMem(Field3D& a_o, const float& c, StreamT st,bool onDev){
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     SetMem(a_o, c, st, onDev));
}
void Opers::
SetMem(Field3D& a_o, const float& c, Image3D& mask, StreamT st,bool onDev){
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     SetMem(a_o, c, mask, st, onDev));
}

/** @defgroup operator
 * Basic operator on the Vector field
 * @{
 */

FIELD_BINARY_ARRAY_OP_DEFS(Add)

FIELD_BINARY_ARRAY_OP_DEFS(Sub)

FIELD_BINARY_ARRAY_OP_DEFS(Mul)

FIELD_BINARY_ARRAY_OP_DEFS(Div)

FIELD_BINARY_ARRAY_OP_DEFS(Max)

FIELD_BINARY_ARRAY_OP_DEFS(Min)

/** @brief a_o = (a_i + a_i1) * c */
void Opers::
AddMulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const float& c, StreamT st,bool onDev){
   MK_CHECK3_ALL(a_o, a_i, a_i1);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     AddMulC(a_o, a_i, a_i1, c, st, onDev));
}
/** @brief a_o = (a_i - a_i1) * c */
void Opers::
SubMulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const float& c, StreamT st,bool onDev){
   MK_CHECK3_ALL(a_o, a_i, a_i1);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     SubMulC(a_o, a_i, a_i1, c, st, onDev));
}
/** @brief a_o = (a_i * a_i1) * c */
void Opers::
MulMulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const float& c, StreamT st,bool onDev){
   MK_CHECK3_ALL(a_o, a_i, a_i1);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     MulMulC(a_o, a_i, a_i1, c, st, onDev));
}


/** @brief a_o = (a_o + a_i) * c */
void Opers::
AddMulC_I(Field3D& a_o, const Field3D& a_i, const float& c, StreamT st,bool onDev){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     AddMulC_I(a_o, a_i, c, st, onDev));
}

/** @brief a_o = (a_o - a_i) * c */
void Opers::
SubMulC_I(Field3D& a_o, const Field3D& a_i, const float& c, StreamT st,bool onDev){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     SubMulC_I(a_o, a_i, c, st, onDev));
}

/** @brief a_o = (a_o * a_i) * c */
void Opers::
MulMulC_I(Field3D& a_o, const Field3D& a_i, const float& c, StreamT st,bool onDev){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     MulMulC_I(a_o, a_i, c, st, onDev));
}


/** @brief a_o = a_i + a_i1 * c */
void Opers::
Add_MulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const float& c, StreamT st,bool onDev){
   MK_CHECK3_ALL(a_o, a_i, a_i1);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     Add_MulC(a_o, a_i, a_i1, c, st, onDev));
}
/** @brief a_o = a_o+ a_i * c */
void Opers::
Add_MulC_I(Field3D& a_o, const Field3D& a_i, const float& c, StreamT st,bool onDev){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     Add_MulC_I(a_o, a_i, c, st, onDev));
}

/** @brief a_o = a_i * c + a_i1 */
void Opers::
MulCAdd(Field3D& a_o, const Field3D& a_i, const float& c, const Field3D& a_i1, StreamT st,bool onDev){
   MK_CHECK3_ALL(a_o, a_i, a_i1);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     MulCAdd(a_o, a_i, c, a_i1, st, onDev));
}

/** @brief a_o = a_i * c - a_i1 */
void Opers::
MulCSub(Field3D& a_o, const Field3D& a_i, const float& c, const Field3D& a_i1, StreamT st,bool onDev){
   MK_CHECK3_ALL(a_o, a_i, a_i1);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     MulCSub(a_o, a_i, c, a_i1, st, onDev));
}

/** @brief a_o = a_o * c + a_i */
void Opers::
MulCAdd_I(Field3D& a_o, const float& c, const Field3D& a_i, StreamT st,bool onDev){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     MulCAdd_I(a_o, c, a_i, st, onDev));
}
/** @brief a_o = a_o * c - a_i */
void Opers::
MulCSub_I(Field3D& a_o, const float& c, const Field3D& a_i, StreamT st,bool onDev){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     MulCSub_I(a_o, c, a_i, st, onDev));
}
/** @brief a_o = (a_i + a) * b */
void Opers::
AddCMulC(Field3D& a_o, const Field3D& a_i, const Vec3Df& a, const Vec3Df& b, StreamT st,bool onDev){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     AddCMulC(a_o, a_i, a, b, st, onDev));
}
/** @brief a_o = (a_i * a) + b */
void Opers::
MulCAddC(Field3D& a_o, const Field3D& a_i, const Vec3Df& a, const Vec3Df& b, StreamT st,bool onDev){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     MulCAddC(a_o, a_i, a, b, st, onDev));
}
/** @brief a_o = (a_o + a) * b */
void Opers::
AddCMulC_I(Field3D& a_o, const Vec3Df& a, const Vec3Df& b, StreamT st,bool onDev){
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     AddCMulC_I(a_o, a, b, st, onDev));
}
/** @brief a_o = (a_o * a) + b */
void Opers::
MulCAddC_I(Field3D& a_o, const Vec3Df& a, const Vec3Df& b, StreamT st,bool onDev){
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     MulCAddC_I(a_o, a, b, st, onDev));
}

/** @brief a_o = a_i * ca + a_i1 * cb */
void Opers::
MulC_Add_MulC(Field3D& a_o,
              const Field3D& a_i, const float& ca,
              const Field3D& a_i1, const float& cb, StreamT st,bool onDev){
   MK_CHECK3_ALL(a_o, a_i, a_i1);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     MulC_Add_MulC(a_o, a_i, ca, a_i1, cb, st, onDev));
}
/** @brief a_o = a_o * co + a_i * cb */
void Opers::
MulC_Add_MulC_I(Field3D& a_o, const float& co,
                const Field3D& a_i, const float& cb, StreamT st,bool onDev){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     MulC_Add_MulC_I(a_o, co, a_i, cb, st, onDev));
}
/** @brief a_o = (a_i + a) * b + c */
void Opers::
AddCMulCAddC(Field3D& a_o, const Field3D& a_i,
             const float& a, const float& b, const float& c, StreamT st,bool onDev){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     AddCMulCAddC(a_o, a_i, a, b, c, st, onDev));
}
/** @brief a_o = (a_o + a) * b + c */
void Opers::
AddCMulCAddC_I(Field3D& a_o, const float& a, const float& b, const float& c, StreamT st,bool onDev){
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     AddCMulCAddC_I(a_o, a, b, c, st, onDev));
}

/** @brief a_o = (a_i + a) * b + c */
void Opers::
AddCMulCAddC(Field3D& a_o, const Field3D& a_i, const Vec3Df& a, const Vec3Df& b, const Vec3Df& c, StreamT st,bool onDev){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     AddCMulCAddC(a_o, a_i, a, b, c, st, onDev));
}

/** @brief a_o = (a_o + a) * b + c */
void Opers::
AddCMulCAddC_I(Field3D& a_o, const Vec3Df& a, const Vec3Df& b, const Vec3Df& c, StreamT st,bool onDev){
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     AddCMulCAddC_I(a_o, a, b, c, st, onDev));
}

/** @brief a_o = a_i + a_i1 * a_i2 * c */
void Opers::
Add_MulMulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1,
            const Field3D& a_i2, const float& c, StreamT st,bool onDev){
   MK_CHECK4_ALL(a_o, a_i, a_i1, a_i2);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     Add_MulMulC(a_o, a_i, a_i1, a_i2, c, st, onDev));
}
/** @brief a_o = a_o + a_i * a_i1 * c */
void Opers::
Add_MulMulC_I(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const float& c, StreamT st,bool onDev){
   MK_CHECK3_ALL(a_o, a_i, a_i1);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     Add_MulMulC_I(a_o, a_i, a_i1, c, st, onDev));
}

/** @brief a_o = a_i - a_i1 * a_i2 * c */
void Opers::
Sub_MulMulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1,
            const Field3D& a_i2, const float& c, StreamT st,bool onDev){
   MK_CHECK4_ALL(a_o, a_i, a_i1, a_i2);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     Sub_MulMulC(a_o, a_i, a_i1, a_i2, c, st, onDev));
}
/** @brief a_o = a_o - a_i * a_i1 * c */
void Opers::
Sub_MulMulC_I(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const float& c, StreamT st,bool onDev){
   MK_CHECK3_ALL(a_o, a_i, a_i1);
   AUTO_EXEC(a_o.memType(), FieldOpers,
	     Sub_MulMulC_I(a_o, a_i, a_i1, c, st, onDev));
}

// template instantiations
#include "FOpers_inst.cxx"

} // end namespace PyCA

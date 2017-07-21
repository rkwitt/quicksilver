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

#include <HFOpers.h>

#include <FOpers.h>
#include <FieldOpers.h>
#include <HFieldOpers.h>

#include <Field3D.h>

namespace PyCA {

//////////////////////////////////////////////////////////////////////// 
// set hfield to identity
// i.e. h(x) = x
////////////////////////////////////////////////////////////////////////
void Opers::
SetToIdentity(Field3D& h, StreamT st){
   AUTO_EXEC(h.memType(), HFieldOpers, 
	     SetToIdentity(h, st));
}

/////////////////////////////////////////////////////////////////////// 
// Convert hfield to velocity field
// i.e. v(x) = (h(x) - I)
////////////////////////////////////////////////////////////////////////

void Opers::
HtoV(Field3D& v, const Field3D& h, StreamT st) {
   MK_CHECK2_ALL(v, h);
   AUTO_EXEC(v.memType(), HFieldOpers, 
	     toV(v, h, st));
}

////////////////////////////////////////////////////////////////////////////
// Hfield to displacement (inplace version)
// i.e. v(x) = v(x) - x
////////////////////////////////////////////////////////////////////////////
void Opers::
HtoV_I(Field3D& v, StreamT st) {
   AUTO_EXEC(v.memType(), HFieldOpers, 
	     toV_I(v, st));
}

template<BackgroundStrategy bg>
void Opers::
ComposeHTranslation(Field3D& h, const Field3D& f, const Vec3Df& t, StreamT s, bool onDev) {
   MK_CHECK_HFIELD_BACKGROUND(bg);
   Opers::ComposeTranslation<bg>(h, f, t, s, onDev);
}

// non-templated version
void Opers::
ComposeHTranslation(Field3D& h, const Field3D& f, const Vec3Df& t, 
			       BackgroundStrategy bg, StreamT s, bool onDev)
{
   MK_CHECK_HFIELD_BACKGROUND(bg);
   Opers::ComposeTranslation(h, f, t, bg, s, onDev);
}

template<BackgroundStrategy bg,bool rescaleVector>
void Opers::ResampleH(Field3D& a_o, const Field3D& a_i, StreamT stream) {
    MK_CHECK_HFIELD_BACKGROUND(bg);
    Opers::Resample<bg, rescaleVector>(a_o, a_i, stream);
}

// non-templated version
void Opers::ResampleH(Field3D& a_o, const Field3D& a_i, 
		       BackgroundStrategy bg, bool rescaleVector, 
		       StreamT stream) 
{
    MK_CHECK_HFIELD_BACKGROUND(bg);
    Opers::Resample(a_o, a_i, bg, rescaleVector, stream);
}

////////////////////////////////////////////////////////////////////////////
// approximate the inverse of an incremental h field using according
// to the following derivation
//
// hInv(x0) = x0 + d
// x0 = h(x0 + d)
// x0 = h(x0) + d // order zero expansion
// d  = x0 - h(x0)
//
// hInv(x0) = x0 + x0 - h(x0)
//
////////////////////////////////////////////////////////////////////////////
void Opers::InverseZerothOrder(Field3D& hInv, const Field3D& h, StreamT s){
   MK_CHECK2_ALL(hInv, h);
   AUTO_EXEC(hInv.memType(), HFieldOpers, 
	     InverseZerothOrder(hInv, h, s));
}

void Opers::InverseZerothOrder_I(Field3D& hInv, StreamT s){
   AUTO_EXEC(hInv.memType(), HFieldOpers, 
	     InverseZerothOrder_I(hInv, s));
}

void Opers::
SplatH(Field3D& a_o, const Field3D& a_h, const Field3D& a_i, bool normalize, StreamT stream)
{
   MK_CHECK3_ALL(a_o, a_h, a_i);
   AUTO_EXEC(a_o.memType(), HFieldOpers, 
	     Splat(a_o, a_h, a_i, normalize, stream));
}

void Opers::
InitializeHFromAffine(Field3D& a_h, const Aff3Df &aff_in,
		     bool invertAff, StreamT s)
{
   AUTO_EXEC(a_h.memType(), HFieldOpers, 
	     initializeFromAffine(a_h, aff_in, invertAff, s));
}

// template instantiations

#include "HFOpers_inst.cxx"

} // end namespace PyCA

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

#include <VFOpers.h>

#include <VFieldOpers.h>
#include <conditionMacro.h>
#include <Field3D.h>

#include "FOpers.h"

namespace PyCA {

//////////////////////////////////////////////////////////////////////// 
// set vfield to zero
// i.e. v(x) = 0
////////////////////////////////////////////////////////////////////////

void Opers::SetToZero(Field3D& v, StreamT s) {
   AUTO_EXEC(v.memType(), VFieldOpers, 
	     SetToZero(v, s));
}

//////////////////////////////////////////////////////////////////////// 
// convert velocity field to HField
// i.e. h(x) = x + v(x) * delta
////////////////////////////////////////////////////////////////////////

void Opers::VtoH(Field3D& h, const Field3D& v, const float& delta ,StreamT s, bool onDev) {
   MK_CHECK2_ALL(h, v)
   AUTO_EXEC(v.memType(), VFieldOpers, 
	     toH(h,v, delta, s, onDev));
}

//////////////////////////////////////////////////////////////////////// 
// convert velocity field to HField (inplace version)
// i.e. h(x) = x + v(x) * delta
////////////////////////////////////////////////////////////////////////

void Opers::VtoH_I(Field3D& h, const float& delta ,StreamT s, bool onDev) {
   AUTO_EXEC(h.memType(), VFieldOpers, 
	     toH_I(h, delta, s, onDev));
}

////////////////////////////////////////////////////////////////////////////
// compose Velocity field with translation
// v(x) = g(x + t)
////////////////////////////////////////////////////////////////////////////

template<BackgroundStrategy bg>
void Opers::ComposeVTranslation(Field3D& v, const Field3D& g, const Vec3Df& t, StreamT s, bool onDev) {
    MK_CHECK_VFIELD_BACKGROUND(bg);
    Opers::ComposeTranslation<bg>(v, g, t, s, onDev);
}

// non-templated version
void Opers::ComposeVTranslation(Field3D& v, const Field3D& g, const Vec3Df& t, 
				 BackgroundStrategy bg, StreamT s, bool onDev) {
    MK_CHECK_VFIELD_BACKGROUND(bg);
    Opers::ComposeTranslation(v, g, t, bg, s, onDev);
}

void Opers::SplatV(Field3D& a_o, const Field3D& a_h, const Field3D& a_i,
                                        bool normalize, StreamT stream) {
   MK_CHECK3_ALL(a_o, a_h, a_i)
   AUTO_EXEC(a_o.memType(), VFieldOpers, 
	     Splat(a_o, a_h, a_i, normalize, stream));
}

template<BackgroundStrategy bg,bool rescaleVector>
void Opers::ResampleV(Field3D& a_o, const Field3D& a_i, StreamT stream) {
    MK_CHECK_VFIELD_BACKGROUND(bg);
    Opers::Resample<bg, rescaleVector>(a_o, a_i, stream);
}

// non-templated version
void Opers::ResampleV(Field3D& a_o, const Field3D& a_i, 
		       BackgroundStrategy bg, bool rescaleVector, 
		       StreamT stream) 
{
    MK_CHECK_VFIELD_BACKGROUND(bg);
    Opers::Resample(a_o, a_i, bg, rescaleVector, stream);
}

// template instantiations
#include "VFOpers_inst.cxx"

} // end namespace PyCA

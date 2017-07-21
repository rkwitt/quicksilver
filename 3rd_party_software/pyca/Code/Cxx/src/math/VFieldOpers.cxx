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

#include <FieldOpers.h>
#include <VFieldOpers.h>
#include <Image3D.h>
#include <Field3D.h>

// Local
#include "CVFieldOpers.h"
#include "GVFieldOpers.h"

namespace PyCA {

//////////////////////////////////////////////////////////////////////// 
// set vfield to zero
// i.e. v(x) = 0
////////////////////////////////////////////////////////////////////////
template<int mode>
void VFieldOpers<mode>::SetToZero(Field3D& v, StreamT s) {
    Executer::SetToZero(v, s);
}

//////////////////////////////////////////////////////////////////////// 
// convert velocity field to HField
// i.e. h(x) = x + v(x) * delta
////////////////////////////////////////////////////////////////////////
template<int mode>
void VFieldOpers<mode>::toH(Field3D& h, const Field3D& v, const float& delta ,StreamT s, bool onDev) {
    Executer::toH(h,v, delta, s, onDev);
}

//////////////////////////////////////////////////////////////////////// 
// convert velocity field to HField (inplace version)
// i.e. h(x) = x + v(x) * delta
////////////////////////////////////////////////////////////////////////
template<int mode>
void VFieldOpers<mode>::toH_I(Field3D& h, const float& delta ,StreamT s, bool onDev) {
    Executer::toH_I(h, delta, s, onDev);
}

template<int mode>
void VFieldOpers<mode>:: Splat(Field3D& a_o, const Field3D& a_h, const Field3D& a_i,
                                        bool normalize, StreamT stream) {
    Executer::Splat(a_o, a_h, a_i, normalize, stream);
}

// template instantiations
#include "VFieldOpers_inst.cxx"

} // end namespace PyCA

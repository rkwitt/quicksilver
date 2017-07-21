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

#include <HFieldOpers.h>

#include <FieldOpers.h>
#include <VFieldOpers.h>

#include <conditionMacro.h>
#include <Field3D.h>
#include <MemoryManager.h>

#include <ThreadSingletonC.h>

// local
#include "CHFieldOpers.h"
#include "GHFieldOpers.h"

namespace PyCA {

//////////////////////////////////////////////////////////////////////// 
// set hfield to identity
// i.e. h(x) = x
////////////////////////////////////////////////////////////////////////
template<int mode>
void HFieldOpers<mode>::
SetToIdentity(Field3D& h, StreamT s) {
    Executer::SetToIdentity(h, s);
}

////////////////////////////////////////////////////////////////////////////
// Hfield to displacement 
// i.e. v(x) = h(x) - x
////////////////////////////////////////////////////////////////////////////

template<int mode>
void HFieldOpers<mode>::
toV(Field3D& a_v, const Field3D& a_h, StreamT s) {
    Executer::toV(a_v, a_h, s);
}
    
////////////////////////////////////////////////////////////////////////////
// Hfield to displacement (inplace version)
// i.e. v(x) = v(x) - x
////////////////////////////////////////////////////////////////////////////
template<int mode>
void HFieldOpers<mode>::
toV_I(Field3D& a_v, StreamT s) {
    Executer::toV_I(a_v, s);
}
    
////////////////////////////////////////////////////////////////////////////
// Hfield to displacement 
// i.e. u(x) = h(x) - x
////////////////////////////////////////////////////////////////////////////
template<int mode>
void HFieldOpers<mode>::
toU(Field3D& a_u, const Field3D& a_h, StreamT s) {
    Executer::toV(a_u, a_h, s);
}
    
////////////////////////////////////////////////////////////////////////////
// Hfield to displacement (place version)
// i.e. u(x) = u(x) - x
////////////////////////////////////////////////////////////////////////////
template<int mode>
void HFieldOpers<mode>::
toU_I(Field3D& a_u, StreamT s) {
    Executer::toV_I(a_u, s);
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
	
template<int mode>
void HFieldOpers<mode>::
InverseZerothOrder(Field3D& a_hInv, const Field3D& a_h, StreamT s) {
    Executer::InverseZerothOrder(a_hInv, a_h, s);
}

    
template<int mode>
void HFieldOpers<mode>::
InverseZerothOrder_I(Field3D& a_hInv, StreamT s) {
    Executer::InverseZerothOrder_I(a_hInv, s);
}
    
template<int mode>
void HFieldOpers<mode>::
Splat(Field3D& a_o, const Field3D& a_h, const Field3D& a_i, bool normalize, StreamT stream)
{
    MK_CHECK3_SIZE(a_o, a_h, a_i);

    ManagedField3D a_vi(a_o.grid(), a_o.memType());

    // Convert from HField to VField
    toV(a_vi, a_i, stream);

    VFieldOpers<mode>::Splat(a_o, a_h, a_vi, normalize, stream);

    // Back from VField
    VFieldOpers<mode>::toH_I(a_o, 1.f, stream);
}

template<int mode>
void HFieldOpers<mode>::
initializeFromAffine(Field3D& a_h, const Aff3Df &aff_in,
		     bool invertAff, StreamT s){
   Executer::initializeFromAffine(a_h, aff_in, invertAff, s);
}

// template instantiations
#include "HFieldOpers_inst.cxx"

} // end namespace PyCA

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

#include "GHFieldOpers.h"
#include "GHFieldOperKernels.h"
#include <Field3D.h>
#include <pycaUtils.h>

namespace PyCA {

//////////////////////////////////////////////////////////////////////// 
// set hfield to identity
// i.e. h(x) = x
////////////////////////////////////////////////////////////////////////
void GHFieldOpers::
SetToIdentity(Field3D& d_h, StreamT stream) {
    Vec3Di sz = d_h.grid().size();
    PyCA::SetToIdentity(d_h.x, d_h.y, d_h.z, sz,
			stream);
}

//////////////////////////////////////////////////////////////////////// 
// convert hfield to velocity field
// i.e. v(x) = (h(x) - x) * sp
////////////////////////////////////////////////////////////////////////
void GHFieldOpers::
toV(Field3D& v, const Field3D& h, StreamT stream) {
    MK_CHECK2_SIZE(v, h);
    Vec3Di sz = v.grid().size();
    Vec3Df sp   = v.grid().spacing();
    PyCA::toV(v.x, v.y, v.z,
	      h.x, h.y, h.z,
	      sz, sp, stream);
}

//////////////////////////////////////////////////////////////////////// 
// convert hfield to velocity field (inplace version)
// i.e. v(x) = (h(x) - x) * sp
////////////////////////////////////////////////////////////////////////
void GHFieldOpers::
toV_I(Field3D& v, StreamT stream) {
    Vec3Di sz = v.grid().size();
    Vec3Df sp   = v.grid().spacing();
    
    PyCA::toV_I(v.x, v.y, v.z,
		sz, sp, 
		stream);

}

//////////////////////////////////////////////////////////////////////// 
// convert displacement field to hfield
// i.e. h(x) = x + u(x)
////////////////////////////////////////////////////////////////////////
void GHFieldOpers::
toH_I(Field3D& v, StreamT stream) {

    Vec3Di sz = v.grid().size();

    PyCA::toH_I(v.x, v.y, v.z, sz, stream);

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
void GHFieldOpers::
InverseZerothOrder(Field3D& a_hInv, const Field3D& a_h, StreamT stream){
    MK_CHECK2_SIZE(a_hInv, a_h);
    Vec3Di sz = a_hInv.grid().size();
    PyCA::InverseZerothOrder
	(a_hInv.x, a_hInv.y, a_hInv.z,
	 a_h.x, a_h.y, a_h.z,
	 sz, stream);
}

void GHFieldOpers::
InverseZerothOrder_I(Field3D& a_hInv, StreamT stream){
    Vec3Di sz = a_hInv.grid().size();

    PyCA::InverseZerothOrder_I
	(a_hInv.x, a_hInv.y, a_hInv.z, sz, stream);
}

//////////////////////////////////////////////////////////////////////// 
// initialize from affine transformation
// i.e. h(x) = Ax
////////////////////////////////////////////////////////////////////////
void GHFieldOpers::
initializeFromAffine(Field3D& d_h, 
		     const Aff3Df &aff_in,
		     bool invertAff,
		     StreamT stream)
{
   Vec3Di sz = d_h.size();
   Vec3Df sp = d_h.spacing();
   Vec3Df org = d_h.origin();
   PyCA::initializeFromAffine
       (d_h.x, d_h.y, d_h.z,
	sz, sp, org,
	aff_in, invertAff,
	stream);
}

} // end namespace PyCA

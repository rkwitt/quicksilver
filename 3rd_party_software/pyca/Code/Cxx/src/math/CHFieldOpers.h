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

#ifndef __CHFIELD3D_OPERS_H
#define __CHFIELD3D_OPERS_H

#include <estream.h>
#include <pycaConst.h>
#include <Vec3D.h>
#include <Aff3D.h>

namespace PyCA {

class Image3D;
class Field3D;

class CHFieldOpers {
public:
    //////////////////////////////////////////////////////////////////////// 
    // set hfield to identity
    // i.e. h(x) = x
    ////////////////////////////////////////////////////////////////////////
    static void SetToIdentity(Field3D& h, StreamT s);

    ////////////////////////////////////////////////////////////////////////////
    // Hfield to displacement 
    // i.e. v(x) = h(x) - x
    ////////////////////////////////////////////////////////////////////////////
    static void toV(Field3D& a_v, const Field3D& a_h, StreamT s);
    
    ////////////////////////////////////////////////////////////////////////////
    // Hfield to displacement (inplace version)
    // i.e. v(x) = v(x) - x
    ////////////////////////////////////////////////////////////////////////////
    static void toV_I(Field3D& a_v, StreamT s);
    
    ////////////////////////////////////////////////////////////////////////////
    // displacement to Hfield (inplace version)
    // i.e. u(x) = u(x) + x
    ////////////////////////////////////////////////////////////////////////////
    static void toH_I(Field3D& a_u, StreamT s=NULL);
    
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
	static void InverseZerothOrder(Field3D& a_hInv, const Field3D& a_h, StreamT s);
    static void InverseZerothOrder_I(Field3D& a_hInv, StreamT s);

   static void 
   initializeFromAffine(Field3D& h_h, 
			const Aff3Df &aff_in,
			bool invertAff,
			StreamT s = NULL);

};

} // end namespace PyCA
    
#endif

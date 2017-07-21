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

#ifndef __HFIELD3D_OPERS_H
#define __HFIELD3D_OPERS_H

#include <estream.h>
#include <pycaConst.h>
#include <Vec3D.h>
#include <Aff3D.h>
#include <Selector.h>

namespace PyCA {

class Image3D;
class Field3D;
class CHFieldOpers;
class GHFieldOpers;

template<int mode>
class HFieldOpers {
public:
    enum { exec_mode = mode };
    typedef typename BinSelector<mode, CHFieldOpers, GHFieldOpers>::Result Executer;

    //////////////////////////////////////////////////////////////////////// 
// set hfield to identity
// i.e. h(x) = x
////////////////////////////////////////////////////////////////////////
    static void SetToIdentity(Field3D& h, StreamT s=NULL);
    
////////////////////////////////////////////////////////////////////////////
// Hfield to displacement 
// i.e. v(x) = h(x) - x
////////////////////////////////////////////////////////////////////////////
    static void toV(Field3D& a_v, const Field3D& a_h, StreamT s=NULL);
    
////////////////////////////////////////////////////////////////////////////
// Hfield to displacement (inplace version)
// i.e. v(x) = v(x) - x
////////////////////////////////////////////////////////////////////////////
    static void toV_I(Field3D& a_v, StreamT s=NULL);
    
////////////////////////////////////////////////////////////////////////////
// Hfield to displacement 
// i.e. u(x) = h(x) - x
////////////////////////////////////////////////////////////////////////////
    static void toU(Field3D& a_u, const Field3D& a_h, StreamT s=NULL);
    
////////////////////////////////////////////////////////////////////////////
// Hfield to displacement (place version)
// i.e. u(x) = u(x) - x
////////////////////////////////////////////////////////////////////////////
    static void toU_I(Field3D& a_u, StreamT s=NULL);
    
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
    static void InverseZerothOrder(Field3D& a_hInv, const Field3D& a_h, StreamT s=NULL);
    
    static void InverseZerothOrder_I(Field3D& a_hInv, StreamT s=NULL);
    
    static void Splat(Field3D& d_vo, const Field3D& d_h, const Field3D& d_vi,
                             bool normalize, StreamT s=NULL);


   static void initializeFromAffine(Field3D& a_h, const Aff3Df &aff_in,
				    bool invertAff, StreamT s = NULL);
};

} // end namespace PyCA
    
#endif

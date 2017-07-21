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

#ifndef __HF_OPERS_H
#define __HF_OPERS_H

#ifndef SWIG
#include<pycaConst.h>
#include <estream.h>
#include <Vec3D.h>
#include <Aff3D.h>
#endif // SWIG

namespace PyCA {

class Image3D;
class Field3D;

namespace Opers {

    //////////////////////////////////////////////////////////////////////// 
// set hfield to identity
// i.e. h(x) = x
////////////////////////////////////////////////////////////////////////
    void SetToIdentity(Field3D& h, StreamT s=NULL);

////////////////////////////////////////////////////////////////////////////
// Hfield to displacement 
// i.e. v(x) = h(x) - x
////////////////////////////////////////////////////////////////////////////

    void HtoV(Field3D& a_v, const Field3D& a_h, StreamT s=NULL);
    
////////////////////////////////////////////////////////////////////////////
// Hfield to displacement (inplace version)
// i.e. v(x) = v(x) - x
////////////////////////////////////////////////////////////////////////////
    
    void HtoV_I(Field3D& a_v, StreamT s=NULL);

    template<BackgroundStrategy bg>
    void ComposeHTranslation(Field3D& h, const Field3D& f, const Vec3Df& t, StreamT s=NULL, bool onDev=false);

   // non-templated version
    void ComposeHTranslation(Field3D& h, const Field3D& f, const Vec3Df& t, 
				   BackgroundStrategy bg = BACKGROUND_STRATEGY_PARTIAL_ID, 
				   StreamT s=NULL, bool onDev=false);

    template<BackgroundStrategy bg,  bool rescaleVector>
    void ResampleH(Field3D& a_o, const Field3D& a_i, 
			 StreamT s = NULL);

   // non-tempated version
    void ResampleH(Field3D& a_o, const Field3D& a_i, 
			 BackgroundStrategy bg = BACKGROUND_STRATEGY_PARTIAL_ID,
			 bool rescaleVector = false,
			 StreamT s = NULL);

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
    void InverseZerothOrder(Field3D& a_hInv, const Field3D& a_h, StreamT s=NULL);

    void InverseZerothOrder_I(Field3D& a_hInv, StreamT s=NULL);
    
    void SplatH(Field3D& d_vo, const Field3D& d_h, const Field3D& d_vi,
                             bool normalize, StreamT s=NULL);


   void InitializeHFromAffine(Field3D& a_h, const Aff3Df &aff_in,
				    bool invertAff, StreamT s = NULL);
} // end namespace Opers

} // end namespace PyCA
    
#endif

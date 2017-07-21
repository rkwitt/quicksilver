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

#ifndef __VF_OPERS_H
#define __VF_OPERS_H

#ifndef SWIG
#include <estream.h>
#include <pycaConst.h>
#include <Vec3D.h>
#endif // SWIG

namespace PyCA {

class Image3D;
class Field3D;
/*
 * This is the operation on the VField or Displacement Field
 *
 */
namespace Opers{

//////////////////////////////////////////////////////////////////////// 
// set vfield to zero
// i.e. v(x) = 0
////////////////////////////////////////////////////////////////////////
    void SetToZero(Field3D& v, StreamT s = NULL);

//////////////////////////////////////////////////////////////////////// 
// convert velocity field to HField
// i.e. h(x) = x + v(x) * delta
////////////////////////////////////////////////////////////////////////
    void VtoH(Field3D& h, const Field3D& v, const float& delta=1.f ,StreamT s = NULL, bool onDev=false);

//////////////////////////////////////////////////////////////////////// 
// convert velocity field to HField (inplace version)
// i.e. h(x) = x + v(x) * delta
////////////////////////////////////////////////////////////////////////
    void VtoH_I(Field3D& h, const float& delta=1.f ,StreamT s = NULL, bool onDev=false);

////////////////////////////////////////////////////////////////////////////
// compose Velocity field with translation
// v(x) = g(x + t)
////////////////////////////////////////////////////////////////////////////
    template<BackgroundStrategy bg>
    void ComposeVTranslation(Field3D& v, const Field3D& g, const Vec3Df& t, 
				   StreamT s = NULL, bool onDev=false);

   // non-templated version
    void ComposeVTranslation(Field3D& v, const Field3D& g, const Vec3Df& t, 
				   BackgroundStrategy bg = BACKGROUND_STRATEGY_PARTIAL_ZERO,
				   StreamT s = NULL, bool onDev=false);

    void SplatV(Field3D& a_o, const Field3D& a_h, const Field3D& a_i,
                             bool normalize, StreamT s = NULL);
    
    template<BackgroundStrategy bg,  bool rescaleVector>
    void ResampleV(Field3D& a_o, const Field3D& a_i, 
			 StreamT s = NULL);

   // non-tempated version
    void ResampleV(Field3D& a_o, const Field3D& a_i, 
			 BackgroundStrategy bg = BACKGROUND_STRATEGY_PARTIAL_ZERO,
			 bool rescaleVector = false,
			 StreamT s = NULL);

} // end namespace Opers

} // end namespace PyCA

#endif

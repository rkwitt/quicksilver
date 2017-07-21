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

#include "CHFieldOpers.h"
#include <Field3D.h>

namespace PyCA {

//////////////////////////////////////////////////////////////////////// 
// set hfield to identity
// i.e. h(x) = x
////////////////////////////////////////////////////////////////////////
void CHFieldOpers::
SetToIdentity(Field3D& h, StreamT st){
    Vec3Di size = h.grid().size();
    size_t id = 0;
    for (int k=0; k < size.z; ++k)
        for (int j=0; j < size.y; ++j)
            for (int i=0; i < size.x; ++i, ++id) {
                h.x[id] = i;
                h.y[id] = j;
                h.z[id] = k;
            }
}

/////////////////////////////////////////////////////////////////////// 
// Convert hfield to velocity field
// i.e. v(x) = (h(x) - I)
////////////////////////////////////////////////////////////////////////
void CHFieldOpers::
toV(Field3D& v, const Field3D& h, StreamT st) {
    MK_CHECK2_SIZE(v, h);

    Vec3Di size = v.grid().size();
    Vec3Df sp   = v.grid().spacing();

    size_t id = 0;
    for (int k=0; k < size.z; ++k)
        for (int j=0; j < size.y; ++j)
            for (int i=0; i < size.x; ++i, ++id) {
                v.x[id] = (h.x[id] - i) * sp.x;
                v.y[id] = (h.y[id] - j) * sp.y;
                v.z[id] = (h.z[id] - k) * sp.z;
            }
}


void CHFieldOpers::
toV_I(Field3D& v, StreamT st) {
    Vec3Di size = v.grid().size();
    Vec3Df sp   = v.grid().spacing();
    
    size_t id = 0;
    for (int k=0; k < size.z; ++k)
        for (int j=0; j < size.y; ++j)
            for (int i=0; i < size.x; ++i, ++id) {
                v.x[id] = (v.x[id] - i) * sp.x;
                v.y[id] = (v.y[id] - j) * sp.y;
                v.z[id] = (v.z[id] - k) * sp.z;
            }
}

void CHFieldOpers::
toH_I(Field3D& v, StreamT st) {
    Vec3Di size = v.grid().size();
    Vec3Df sp   = v.grid().spacing();
    
    size_t id = 0;
    for (int k=0; k < size.z; ++k)
        for (int j=0; j < size.y; ++j)
            for (int i=0; i < size.x; ++i, ++id) {
                v.x[id] += i;
                v.y[id] += j;
                v.z[id] += k;
            }
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
void CHFieldOpers::
InverseZerothOrder(Field3D& hInv, const Field3D& h, StreamT s){
    MK_CHECK2_SIZE(hInv, h);
    Vec3Di size = hInv.grid().size();
    size_t id = 0;
    for (int k=0; k < size.z; ++k)
        for (int j=0; j < size.y; ++j)
            for (int i=0; i < size.x; ++i, ++id) {
                hInv.x[id] = i + i - h.x[id];
                hInv.y[id] = j + j - h.y[id];
                hInv.z[id] = k + k - h.z[id];
            }
}


void CHFieldOpers::
InverseZerothOrder_I(Field3D& hInv, StreamT s){
    Vec3Di size = hInv.grid().size();

    size_t id = 0;
    for (int k=0; k < size.z; ++k)
        for (int j=0; j < size.y; ++j)
            for (int i=0; i < size.x; ++i, ++id) {
                hInv.x[id] = i + i - hInv.x[id];
                hInv.y[id] = j + j - hInv.y[id];
                hInv.z[id] = k + k - hInv.z[id];
            }
}

void CHFieldOpers::
initializeFromAffine(Field3D& h_h, 
		     const Aff3Df &aff_in,
		     bool invertAff,
		     StreamT s)
{
   Vec3Di size = h_h.size();
   Vec3Df sp = h_h.spacing();
   Vec3Df org = h_h.origin();
    
   Aff3Df aff = aff_in;
   if(invertAff){
      if(!aff.invert()){
	 throw PyCAException(__FILE__,__LINE__,"Error, could not invert affine transform");
      }
   }

  for (int z = 0; z < size.z; ++z) {
      for (int y = 0; y < size.y; ++y) {
	  for (int x = 0; x < size.x; ++x) {
	      Vec3Df p,tp;
	      // convert to world space
	      p[0] = x * sp[0] + org[0];
	      p[1] = y * sp[1] + org[1];
	      p[2] = z * sp[2] + org[2];
	      // transform via affine transform
	      aff.transformVector(p,tp);
	      // convert back to index coords
	      tp -= org;
	      tp /= sp;
	      // assign result
	      h_h.set(x,y,z,tp);
	    }
	}
    }
}
  
} // end namespace PyCA

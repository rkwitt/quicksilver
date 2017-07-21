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

#include "GVFieldOpers.h"
#include "GVFieldOperKernels.h"
#include <MemoryManager.h>
#include <pycaUtils.h>
#include <boost_mem.h>

namespace PyCA {

void GVFieldOpers::
SetToZero(Field3D& d_v, StreamT stream) {
    
    Vec3Di sz = d_v.grid().size();
    PyCA::SetToZero(d_v.x, d_v.y, d_v.z,
		    sz.x, sz.y, sz.z,
		    stream);

}

void GVFieldOpers::
toH(Field3D& d_h, const Field3D& d_v, const float& delta, StreamT stream, bool onDev){

    Vec3Di sz = d_h.grid().size();
    Vec3Df sp = d_h.grid().spacing();
    PyCA::toH(d_h.x, d_h.y, d_h.z,
	      sp.x, sp.y, sp.z,
	      d_v.x, d_v.y, d_v.z,
	      sz.x, sz.y, sz.z,
	      delta, 
	      stream, onDev);
}

void GVFieldOpers::
toH_I(Field3D& d_h, const float& delta, StreamT stream, bool onDev){

    Vec3Di sz = d_h.grid().size();
    Vec3Df sp = d_h.grid().spacing();

    PyCA::toH_I(d_h.x, d_h.y, d_h.z,
		sp.x, sp.y, sp.z,
		sz.x, sz.y, sz.z,
		delta,
		stream, onDev);

}


// this needs to be in .cxx file to sequester reference to boost
// shared_ptr
void GVFieldOpers::
Splat(Field3D& d_o, const Field3D& d_h, const Field3D& d_i,
             bool normalize, StreamT stream)
{
    MK_CHECK3_SIZE(d_o, d_h, d_i);

    Vec3Di size = d_o.grid().size();
    size_t nVox = size.prod();

    STATIC_ASSERT(sizeof(int) == sizeof(float));
    
    if (normalize) {
        boost::shared_ptr<int> distBuffer = CreateSharedArray<MEM_DEVICE, int>(nVox);
        int* i_dd = (int*) distBuffer.get();

        splatNormalized(d_o.x, d_o.y, d_o.z, i_dd,
			size.x, size.y, size.z,
			d_i.x, d_i.y, d_i.z,
			d_h.x, d_h.y, d_h.z, nVox, stream);
    } else {
        splatUnnormalized(d_o.x, d_o.y, d_o.z,
			  size.x, size.y, size.z,
			  d_i.x, d_i.y, d_i.z,
			  d_h.x, d_h.y, d_h.z, nVox, stream);

    }
}

} // end namespace PyCA

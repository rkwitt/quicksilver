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

#include "CVFieldOpers.h"
#include <vector>
#include <Field3D.h>
#include <MemOpers.h>

namespace PyCA {

//////////////////////////////////////////////////////////////////////// 
// set Vfield to zero
// i.e. v(x) = 0
////////////////////////////////////////////////////////////////////////
void CVFieldOpers
::SetToZero(Field3D& v, StreamT st){
    size_t n = v.nVox();
    for (size_t id = 0; id < n; ++id) {
        v.x[id] = 0;
        v.y[id] = 0;
        v.z[id] = 0;
    }
}

//////////////////////////////////////////////////////////////////////// 
// Convert velocity field to hfield
// i.e. h(x) = x + v(x) * delta
// delta = 1.f by default
////////////////////////////////////////////////////////////////////////

void CVFieldOpers::
toH(Field3D& h, const Field3D& v, const float& delta, StreamT st, bool onDev){
    MK_CHECK2_SIZE(h, v);
    Vec3Di size = h.grid().size();
    Vec3Df sp   = h.grid().spacing();
    Vec3Df iSp  = Vec3Df(1.f/sp.x, 1.f/sp.y, 1.f/sp.z);

    size_t id = 0;
    for (int k=0; k < size.z; ++k)
        for (int j=0; j < size.y; ++j)
            for (int i=0; i < size.x; ++i, ++id) {
                h.x[id] = i + v.x[id] * (delta * iSp.x);
                h.y[id] = j + v.y[id] * (delta * iSp.y);
                h.z[id] = k + v.z[id] * (delta * iSp.z);
            }
}

//////////////////////////////////////////////////////////////////////// 
// Convert velocity field to hfield (inplace version)
////////////////////////////////////////////////////////////////////////

void CVFieldOpers::
toH_I(Field3D& h, const float& delta, StreamT st, bool onDev){
    Vec3Di size = h.grid().size();
    Vec3Df sp   = h.grid().spacing();
    Vec3Df iSp  = Vec3Df(1.f/sp.x, 1.f/sp.y, 1.f/sp.z);
        
    size_t id = 0;
    for (int k=0; k < size.z; ++k)
        for (int j=0; j < size.y; ++j)
            for (int i=0; i < size.x; ++i, ++id) {
                h.x[id] = i + h.x[id] * (delta * iSp.x);
                h.y[id] = j + h.y[id] * (delta * iSp.y);
                h.z[id] = k + h.z[id] * (delta * iSp.z);
            }
}

void CVFieldOpers::
Splat(Field3D& a_o, const Field3D& a_h, const Field3D& a_i, bool normalize,
               StreamT stream)
{
    Vec3Di size = a_o.grid().size();
    size_t nVox = size.prod();

    size_t id = 0;
    size_t planeSize = size.x * size.y;

    std::vector<float> a_w(nVox);

    size_t nEl = a_o.nVox();
    MemOpers<EXEC_CPU, float>::SetMem(a_o.x, 0.f, nEl, stream);
    MemOpers<EXEC_CPU, float>::SetMem(a_o.y, 0.f, nEl, stream);
    MemOpers<EXEC_CPU, float>::SetMem(a_o.z, 0.f, nEl, stream);
    std::fill(a_w.begin(), a_w.end(), 0.f);
        
    for (int z = 0; z < size.z; ++z) {
        for (int y = 0; y < size.y; ++y) {
            for (int x = 0; x < size.x; ++x, ++id) {
                float ix = a_i.x[id], iy = a_i.y[id], iz = a_i.x[id];
                
                // get h field value---where this intensity should go
                float hx = a_h.x[id];
                float hy = a_h.y[id];
                float hz = a_h.z[id];
                
                // this is a fast version of the floor function
                int fx = static_cast<int>(hx);
                int fy = static_cast<int>(hy);
                int fz = static_cast<int>(hz);
                
                if (hx < 0 && hx != static_cast<int>(hx)) --fx;
                if (hy < 0 && hy != static_cast<int>(hy)) --fy;
                if (hz < 0 && hz != static_cast<int>(hz)) --fz;

                if (hx > -1 && hx < (int) size.x &&  // inside vol with 1 pix border
                    hy > -1 && hy < (int) size.y &&
                    hz > -1 && hz < (int) size.z)
                {
                    // compute trilinear weights
                    float rx = hx - fx, ry = hy - fy, rz = hz - fz;
                    float w000 = (1.0 - rx) * (1.0 - ry) * (1.0 - rz);
                    float w001 = (1.0 - rx) * (1.0 - ry) * (rz);
                    float w010 = (1.0 - rx) * (ry)       * (1.0 - rz);
                    float w011 = (1.0 - rx) * (ry)       * (rz);
                    float w100 = (rx)       * (1.0 - ry) * (1.0 - rz);
                    float w101 = (rx)       * (1.0 - ry) * (rz);
                    float w110 = (rx)       * (ry)       * (1.0 - rz);
                    float w111 = (rx)       * (ry)       * (rz);

                    // see which corners of cube are valid
                    bool floorXIn = (fx >= 0), ceilXIn = (fx < (int) size.x-1),
                         floorYIn = (fy >= 0), ceilYIn = (fy < (int) size.y-1),
                         floorZIn = (fz >= 0), ceilZIn = (fz < (int) size.z-1);

                    size_t nid = (fz * size.y + fy) * size.x + fx;
                    if (floorXIn && floorYIn && floorZIn) {
                        a_o.x[nid]     += w000 * ix;
                        a_o.y[nid]     += w000 * iy;
                        a_o.z[nid]     += w000 * iz;
                        a_w[nid]       += w000;
                    }
                    if (floorXIn && floorYIn && ceilZIn){
                        a_o.x[nid + planeSize]  += w001 * ix;
                        a_o.y[nid + planeSize]  += w001 * iy;
                        a_o.z[nid + planeSize]  += w001 * iz;
                        a_w[nid + planeSize]  += w001;
                    }
                    if (floorXIn && ceilYIn && floorZIn){
                        a_o.x[nid + size.x]  += w010 * ix;
                        a_o.y[nid + size.x]  += w010 * iy;
                        a_o.z[nid + size.x]  += w010 * iz;
                        a_w[nid + size.x]  += w010;
                    }
                    if (floorXIn && ceilYIn && ceilZIn) {
                        a_o.x[nid + size.x + planeSize] += w011 * ix;
                        a_o.y[nid + size.x + planeSize] += w011 * iy;
                        a_o.z[nid + size.x + planeSize] += w011 * iz;
                        a_w[nid + size.x + planeSize] += w011;
                    }
                    if (ceilXIn && floorYIn && floorZIn) {
                        a_o.x[nid + 1] += w100 * ix;
                        a_o.y[nid + 1] += w100 * iy;
                        a_o.z[nid + 1] += w100 * iz;
                        a_w[nid + 1] += w100;
                    }
                    if (ceilXIn && floorYIn && ceilZIn) {
                        a_o.x[nid + 1 + planeSize] += w101 * ix;
                        a_o.y[nid + 1 + planeSize] += w101 * iy;
                        a_o.z[nid + 1 + planeSize] += w101 * iz;
                        a_w[nid + 1 + planeSize]   += w101;
                    }
                    if (ceilXIn && ceilYIn && floorZIn)
                    {
                        a_o.x[nid + 1 + size.x] += w110 * ix;
                        a_o.y[nid + 1 + size.x] += w110 * iy;
                        a_o.z[nid + 1 + size.x] += w110 * iz;
                        a_w[nid + 1 + size.x] += w110;
                    }
                    if (ceilXIn && ceilYIn && ceilZIn)
                    {
                        a_o.x[nid + 1 + size.x + planeSize] += w111 * ix;
                        a_o.y[nid + 1 + size.x + planeSize] += w111 * iy;
                        a_o.z[nid + 1 + size.x + planeSize] += w111 * iz;
                        
                        a_w[nid + 1 + size.x + planeSize] += w111;
                    }
                }
            }
        }
    }
    // normalize counts (NOTE: no rounding for integer types)
    if(normalize){
        size_t numElements = a_o.nVox();
        for (size_t e = 0; e < numElements; ++e) {
            if (a_w[e] > 0) {
                a_o.x[e] = a_o.x[e]/ a_w[e];
                a_o.y[e] = a_o.y[e]/ a_w[e];
                a_o.z[e] = a_o.z[e]/ a_w[e];
            }
            else {
                a_o.x[e] = a_o.y[e] = a_o.z[e] = 0.f;
            }
            
        }
    }
}
} // end namespace PyCA

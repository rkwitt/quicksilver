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

#include "GImageOpers.h"
#include "GImageOperKernels.h"
#include <pycaUtils.h>
#include "interp.h"
#include <Image3D.h>
#include <Field3D.h>
#include <MemOpers.h>

namespace PyCA {

void 
GImageOpers::
SubVol(Image3D& d_o,const Image3D& d_i,
       const Vec3Di& start, StreamT st)
{
   SubVol(d_o.get(), d_i.get(), d_o.size(), d_i.size(), start, st);
}

void 
GImageOpers::
SubVol(float* d_o,const float* d_i, 
       const Vec3Di& oSize, const Vec3Di& iSize, 
       const Vec3Di& start, StreamT st)
{

    PyCA::SubVol(d_o, d_i, oSize, iSize, start, st);
}

void 
GImageOpers::
SetSubVol_I(Image3D& d_o,const Image3D& d_i,
	    const Vec3Di& start, StreamT st)
{
    SetSubVol_I(d_o.get(), d_i.get(), d_o.size(), d_i.size(), start, st);
}

void 
GImageOpers::
SetSubVol_I(float* d_o,const float* d_i, 
	    const Vec3Di& oSize, const Vec3Di& iSize, 
	    const Vec3Di& start, StreamT st)
{
    PyCA::SetSubVol_I(d_o, d_i, oSize, iSize, start, st);
}

void GImageOpers::
Shrink_I(Image3D& d_o, float c, StreamT st)
{
   Shrink(d_o, d_o, c, st);
}

void GImageOpers::
Shrink(Image3D& d_o, const Image3D& d_i, float c, StreamT st)
{
   Vec3Di sz = d_i.size();
   PyCA::Shrink(d_o.get(), d_i.get(), sz, c, st);
}

void GImageOpers::
SoftAbs_I(Image3D& d_o, float eps, StreamT st)
{
   SoftAbs(d_o, d_o, eps, st);
}

void GImageOpers::
SoftAbs(Image3D& d_o, const Image3D& d_i, float eps, StreamT st)
{
    Vec3Di sz = d_i.size();
    PyCA::SoftAbs(d_o.get(), d_i.get(), sz, eps, st);
}

void GImageOpers::
SoftSgn_I(Image3D& d_o, float eps, StreamT st)
{
   SoftSgn(d_o, d_o, eps, st);
}

void GImageOpers::
SoftSgn(Image3D& d_o, const Image3D& d_i, float eps, StreamT st)
{
   Vec3Di sz = d_i.size();
   PyCA::SoftSgn(d_o.get(), d_i.get(), sz, eps, st);
}

template<BackgroundStrategy bg, InterpT interp, bool useOriginOffset>
void GImageOpers::
Resample(Image3D& d_o, const Image3D& d_i, StreamT stream)
{
    Vec3Di oSz = d_o.size();
    Vec3Di iSz = d_i.size();

    if (oSz == iSz){
	size_t nEl = oSz.prod();
	MemOpers<EXEC_GPU, float>::Copy(d_o.get(), 
					d_i.get(), 
					nEl, stream);
	return;
    }
    
    PyCA::Resample<bg, interp, useOriginOffset>
	(d_o.get(), oSz,
	 d_i.get(), iSz, 
	 stream);
}

template<BackgroundStrategy bg, InterpT interp>
void GImageOpers::
ResampleWorld(Image3D& d_o, const Image3D& d_i, StreamT stream)
{

    if (d_o.grid() == d_i.grid()){
	size_t nEl = d_o.nVox();
	MemOpers<EXEC_GPU, float>::Copy(d_o.get(), 
					d_i.get(), 
					nEl, stream);
	return;
    }

    Vec3Di oSz = d_o.size();
    Vec3Di iSz = d_i.size();
    Vec3Df oSp = d_o.spacing();
    Vec3Df iSp = d_i.spacing();
    Vec3Df oOr = d_o.origin();
    Vec3Df iOr = d_i.origin();
    PyCA::ResampleWorld<bg, interp>
	(d_o.get(), oSz, oSp, oOr,
	 d_i.get(), iSz, iSp, iOr,
	 stream);
}

template<BackgroundStrategy bg>
void GImageOpers::
SplatWorld(Image3D& d_o, const Image3D& d_i, StreamT stream)
{

    Vec3Di oSz = d_o.size();
    Vec3Di iSz = d_i.size();
    Vec3Df oSp = d_o.spacing();
    Vec3Df iSp = d_i.spacing();
    Vec3Df oOr = d_o.origin();
    Vec3Df iOr = d_i.origin();

    PyCA::SplatWorld<bg>
	(d_o.get(), oSz, oSp, oOr,
	 d_i.get(), iSz, iSp, iOr,
	 stream);
}

template<BackgroundStrategy bg>
void GImageOpers::
SplatWorld(Image3D& d_o, const Image3D& d_i, Image3D& d_w, StreamT stream)
{

    Vec3Di oSz = d_o.size();
    Vec3Di iSz = d_i.size();
    Vec3Df oSp = d_o.spacing();
    Vec3Df iSp = d_i.spacing();
    Vec3Df oOr = d_o.origin();
    Vec3Df iOr = d_i.origin();
    dim3 threads(16,16);
    dim3 grids(iDivUp(oSz.x, threads.x), iDivUp(oSz.y, threads.y));
    PyCA::SplatWorld<bg>
	(d_o.get(), oSz, oSp, oOr,
	 d_i.get(), iSz, iSp, iOr,
	 d_w.get(),
	 stream);
}

void GImageOpers::
Convolve(Image3D& d_o, const Image3D& d_i, 
	 const Image3D& d_kernel, 
	 StreamT stream)
{
   MK_CHECK2_SIZE(d_o, d_i);
   Vec3Di sz = d_i.size();
   Vec3Di kSz = d_kernel.size();
   PyCA::Convolve(d_o.get(),
		  d_i.get(), sz,
		  d_kernel.get(), kSz,
		  stream);
   
}






// template instantiations
#include "GImageOpers_inst.cxx"

} // end namespace PyCA

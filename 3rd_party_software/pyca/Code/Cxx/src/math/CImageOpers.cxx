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

#include "CImageOpers.h"
#include <pycaUtils.h>
#include "interp.h"
#include <Image3D.h>
#include <conditionMacro.h>
#include <MemOpers.h>
#include <algorithm>

#include "CSplat.h"

namespace PyCA {

void CImageOpers::
SubVol(Image3D& a_o,const Image3D& a_i, 
       const Vec3Di& start, StreamT st)
{
   SubVol(a_o.get(), a_i.get(), a_o.size(), a_i.size(), start, st);
}

void CImageOpers::
SubVol(float* a_o,const float* a_i, 
       const Vec3Di& oSize, const Vec3Di& iSize, 
       const Vec3Di& start, StreamT st)
{
    // make sure extracted reguion is within a_i bounds
   Vec3Di oStart = std::max(Vec3Di(0,0,0),-start);
   Vec3Di oEnd = std::min(oSize,iSize-start);

    MemOpers<EXEC_CPU,float>::SetMem(a_o, 0.f, oSize.prod());
    for (int k=oStart.z; k < oEnd.z; ++k) {
       for (int j=oStart.y; j < oEnd.y; ++j) {
	  for (int i=oStart.x; i < oEnd.x; ++i) {
	     int o_id = oSize.x*(k*oSize.y + j) + i;
	     int i_id = 
		iSize.x*((start.z+k)*iSize.y + (start.y+j)) 
		+ start.x+i;
	     a_o[o_id] = a_i[i_id];
	  }
       }
    }
}

void CImageOpers::
SetSubVol_I(Image3D& a_o,const Image3D& a_i, 
       const Vec3Di& start, StreamT st)
{
   SetSubVol_I(a_o.get(), a_i.get(), a_o.size(), a_i.size(), start, st);
}

void CImageOpers::
SetSubVol_I(float* a_o,const float* a_i, 
	    const Vec3Di& oSize, const Vec3Di& iSize, 
	    const Vec3Di& start, StreamT st)
{
    // make sure extracted reguion is within a_i bounds
   Vec3Di iStart = std::max(Vec3Di(0,0,0),-start);
   Vec3Di iEnd = std::min(iSize,oSize-start);

   for (int k=iStart.z; k < iEnd.z; ++k) {
      for (int j=iStart.y; j < iEnd.y; ++j) {
	 for (int i=iStart.x; i < iEnd.x; ++i) {
	    int i_id = iSize.x*(k*iSize.y + j) + i;
	    int o_id = 
	       oSize.x*((start.z+k)*oSize.y + (start.y+j)) 
	       + start.x+i;
	    a_o[o_id] = a_i[i_id];
	 }
       }
   }
}

void CImageOpers::
Shrink_I(Image3D& a_o, float c, StreamT st)
{
   Shrink(a_o, a_o, c, st);
}

// safe for a_o == a_i
void CImageOpers::
Shrink(Image3D& a_o, const Image3D& a_i, float c, StreamT st)
{
   MK_CHECK2_SIZE(a_o, a_i);
   Vec3Di sz = a_i.size();

   for (int k=0; k < sz.z; ++k) {
      for (int j=0; j < sz.y; ++j) {
	 for (int i=0; i < sz.x; ++i) {
	    int id = sz.x*(k*sz.y + j) + i;
	    float v = a_i[id];
	    a_o[id] = std::max(v-c,0.f)+std::min(v+c,0.f);
	 }
      }
   }
}

void CImageOpers::
SoftAbs_I(Image3D& a_o, float eps, StreamT st)
{
   SoftAbs(a_o, a_o, eps, st);
}

// safe for a_o == a_i
void CImageOpers::
SoftAbs(Image3D& a_o, const Image3D& a_i, float eps, StreamT st)
{
   MK_CHECK2_SIZE(a_o, a_i);

   size_t nVox = a_o.nVox();

   for(size_t i=0;i<nVox;++i){
      float i_v = a_i[i];
      if(i_v < -eps){
	 a_o[i] = -i_v-eps/2.f;
      }else if(i_v > eps){
	 a_o[i] = i_v-eps/2.f;
      }else{
	 a_o[i] = (i_v*i_v)/(2.f*eps);
      }
   }
}

void CImageOpers::
SoftSgn_I(Image3D& a_o, float eps, StreamT st)
{
   SoftSgn(a_o, a_o, eps, st);
}

// safe for a_o == a_i
void CImageOpers::
SoftSgn(Image3D& a_o, const Image3D& a_i, float eps, StreamT st)
{
   MK_CHECK2_SIZE(a_o, a_i);

   size_t nVox = a_o.nVox();

   for(size_t i=0;i<nVox;++i){
      float i_v = a_i[i];
      if(i_v < -eps){
	 a_o[i] = -1.f;
      }else if(i_v > eps){
	 a_o[i] = 1.f;
      }else{
	 a_o[i] = i_v/eps;
      }
   }
}

template<BackgroundStrategy bg, InterpT interp, bool useOriginOffset>
void CImageOpers::
Resample(Image3D& a_o, const Image3D& a_i, StreamT stream){
    MK_CHECK_IMAGE_BACKGROUND(bg);
    if (a_o.size() == a_i.size()){
       size_t nEl = a_o.nVox();
       MemOpers<EXEC_CPU, float>::Copy(a_o.get(), a_i.get(), nEl, stream);
       return;
    }
    Vec3Di oSize = a_o.size();
    Vec3Di iSize = a_i.size();

    size_t id = 0;
    
    float rX = (float)iSize.x/ (float)oSize.x;
    float rY = (float)iSize.y/ (float)oSize.y;
    float rZ = (float)iSize.z/ (float)oSize.z;

    float offX=0.f, offY=0.f, offZ=0.f;
    
    if(useOriginOffset){
        offX = (rX-1.f)/2.f;
        offY = (rY-1.f)/2.f;
        offZ = (rZ-1.f)/2.f;
    }

    for (int k=0; k < oSize.z; ++k) {
        float i_z =  offZ + k * rZ;
        for (int j=0; j < oSize.y; ++j) {
            float i_y =  offY + j * rY;
            for (int i=0; i < oSize.x; ++i, ++id) {
                float i_x =  offX + i * rX;
                
                a_o[id] = point_interp<interp, bg>
		    (a_i.get(),
		     i_x, i_y, i_z,
		     iSize.x, iSize.y, iSize.z);
            }
        }
    }
}

template<BackgroundStrategy bg, InterpT interp>
void CImageOpers::
ResampleWorld(Image3D& a_o, const Image3D& a_i, StreamT stream)
{
   Vec3Di i_size = a_i.size();
   Vec3Di o_size = a_o.size();
   Vec3Df i_sp = a_i.spacing();
   Vec3Df o_sp = a_o.spacing();
   Vec3Df i_org = a_i.origin();
   Vec3Df o_org = a_o.origin();
   Vec3Df world;
   Vec3Df voxel;
   //Vec3Df samplingOrigin = (o_sp - i_sp)/2.0;
   size_t id = 0;
   for (int z = 0; z < o_size.z; ++z) {
      for (int y = 0; y < o_size.y; ++y) {
	 for (int x = 0; x < o_size.x; ++x, ++id) {
	    world.x = x*o_sp.x + o_org.x;
	    world.y = y*o_sp.y + o_org.y;
	    world.z = z*o_sp.z + o_org.z;
	    
	    // offset for proper centering
		//world = world + samplingOrigin;
	    
	    voxel.x = (world.x-i_org.x)/i_sp.x;
	    voxel.y = (world.y-i_org.y)/i_sp.y;
	    voxel.z = (world.z-i_org.z)/i_sp.z;
	    
	    a_o[id] = point_interp<interp, bg>
		(a_i.get(), 
		 voxel.x, voxel.y, voxel.z, 
		 i_size.x, i_size.y, i_size.z);
	 }
      }
   }
}

template<BackgroundStrategy bg>
void CImageOpers::
SplatWorld(Image3D& a_o, const Image3D& a_i, StreamT stream)
{
    Vec3Di i_size = a_i.size();
    Vec3Di o_size = a_o.size();
    Vec3Df i_sp = a_i.spacing();
    Vec3Df o_sp = a_o.spacing();
    Vec3Df i_org = a_i.origin();
    Vec3Df o_org = a_o.origin();
    Vec3Df world;
    Vec3Df voxel;
    size_t id = 0;
    for (int z = 0; z < i_size.z; ++z) {
        for (int y = 0; y < i_size.y; ++y) {
            for (int x = 0; x < i_size.x; ++x, ++id) {
                world.x = x*i_sp.x + i_org.x;
                world.y = y*i_sp.y + i_org.y;
                world.z = z*i_sp.z + i_org.z;

                voxel.x = (world.x-o_org.x)/o_sp.x;
                voxel.y = (world.y-o_org.y)/o_sp.y;
                voxel.z = (world.z-o_org.z)/o_sp.z;

                Splatting::splatPoint(voxel.x, voxel.y, voxel.z, a_i[id], 
                        a_o.get(),
                        o_size.x, o_size.y, o_size.z);
            }
        }
    }
}

template<BackgroundStrategy bg>
void CImageOpers::
SplatWorld(Image3D& a_o, const Image3D& a_i, Image3D& a_w, StreamT stream)
{
    Vec3Di i_size = a_i.size();
    Vec3Di o_size = a_o.size();
    Vec3Df i_sp = a_i.spacing();
    Vec3Df o_sp = a_o.spacing();
    Vec3Df i_org = a_i.origin();
    Vec3Df o_org = a_o.origin();
    Vec3Df world;
    Vec3Df voxel;
    size_t id = 0;
    for (int z = 0; z < i_size.z; ++z) {
        for (int y = 0; y < i_size.y; ++y) {
            for (int x = 0; x < i_size.x; ++x, ++id) {
                world.x = x*i_sp.x + i_org.x;
                world.y = y*i_sp.y + i_org.y;
                world.z = z*i_sp.z + i_org.z;

                voxel.x = (world.x-o_org.x)/o_sp.x;
                voxel.y = (world.y-o_org.y)/o_sp.y;
                voxel.z = (world.z-o_org.z)/o_sp.z;

                Splatting::splatPoint(voxel.x, voxel.y, voxel.z, a_i[id], 
                        a_o.get(), a_w.get(),
                        o_size.x, o_size.y, o_size.z);
            }
        }
    }
}

void CImageOpers::
Convolve(Image3D& h_o, const Image3D& h_i, 
	 const Image3D& h_kernel, 
	 StreamT stream)
{
   MK_CHECK2_SIZE(h_o, h_i);
   Vec3Di iSz = h_i.size();
   Vec3Di kSz = h_kernel.size();

   if(kSz.x%2 == 0 || kSz.y%2 == 0 || kSz.z%2 == 0){
      throw PyCAException(__FILE__, __LINE__,
				   "Only odd-sized kernels allowed in convolution");
   }

   Vec3Di halfSz = kSz/2;

   // loop over pixels in image
   for (int cz=0; cz < iSz.z; ++cz) {
      for (int cy=0; cy < iSz.y; ++cy) {
	 for (int cx=0; cx < iSz.x; ++cx) {
	    
	    // loop over offsets in kernel
	    float v = 0.f;
	    for (int oz=-halfSz.z; oz <= halfSz.z; ++oz) {
	       for (int oy=-halfSz.y; oy <= halfSz.y; ++oy) {
		  for (int ox=-halfSz.x; ox <= halfSz.x; ++ox) {
		     float kv = getVal<float>
			(h_kernel.get(), 
			 kSz.x,kSz.y,kSz.z,
			 ox+halfSz.x,oy+halfSz.y,oz+halfSz.z);
		     float iv = getSafeVal<float,BACKGROUND_STRATEGY_CLAMP>
			(h_i.get(), 
			 iSz.x,iSz.y,iSz.z,
			 cx+ox,cy+oy,cz+oz);
		     v += kv*iv;
		  }
	       }
	    }
	    getVal<float>(h_o.get(),iSz.x,iSz.y,iSz.z,cx,cy,cz) = v;
	 }
      }
   }
}

// template instantiations
#include "CImageOpers_inst.cxx"

} // end namespace PyCA

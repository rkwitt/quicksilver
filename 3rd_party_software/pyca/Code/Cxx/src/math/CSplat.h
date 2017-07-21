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

#ifndef __CPU_SPLAT_H
#define __CPU_SPLAT_H

#include <pycaConst.h>
#include <vector>

namespace PyCA {
  namespace Splatting{
    //////////////////////////////////////////////////////////////////////////////////////////
    // splatPoint function for single array input with no weighted normalization
    //////////////////////////////////////////////////////////////////////////////////////////

    inline void splatPoint(const float hx, const float hy, const float hz, const float ix ,float* outIm, int sizeX, int sizeY, int sizeZ){
      // this is a fast version of the floor function
    
      int fx = static_cast<int>(hx);
      int fy = static_cast<int>(hy);
      int fz = static_cast<int>(hz);
      size_t planeSize = sizeX * sizeY;
                
      if (hx < 0 && hx != static_cast<int>(hx)) --fx;
      if (hy < 0 && hy != static_cast<int>(hy)) --fy;
      if (hz < 0 && hz != static_cast<int>(hz)) --fz;

      if (hx > -1 && hx < (int) sizeX &&  // inside vol with 1 pix border
	  hy > -1 && hy < (int) sizeY &&
	  hz > -1 && hz < (int) sizeZ)
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
	  bool floorXIn = (fx >= 0), ceilXIn = (fx < (int) sizeX-1),
	    floorYIn = (fy >= 0), ceilYIn = (fy < (int) sizeY-1),
	    floorZIn = (fz >= 0), ceilZIn = (fz < (int) sizeZ-1);

	  size_t nid = (fz * sizeY + fy) * sizeX + fx;
	  if (floorXIn && floorYIn && floorZIn) {
	    outIm[nid]       += w000 * ix;
	  }
	  if (floorXIn && floorYIn && ceilZIn){
	    outIm[nid + planeSize]  += w001 * ix;
	  }
	  if (floorXIn && ceilYIn && floorZIn){
	    outIm[nid + sizeX]  += w010 * ix;
	  }
	  if (floorXIn && ceilYIn && ceilZIn) {
	    outIm[nid + sizeX + planeSize] += w011 * ix;
	  }
	  if (ceilXIn && floorYIn && floorZIn) {
	    outIm[nid + 1] += w100 * ix;
	  }
	  if (ceilXIn && floorYIn && ceilZIn) {
	    outIm[nid + 1 + planeSize] += w101 * ix;
	  }
	  if (ceilXIn && ceilYIn && floorZIn)
	    {
	      outIm[nid + 1 + sizeX] += w110 * ix;
	    }
	  if (ceilXIn && ceilYIn && ceilZIn)
	    {
	      outIm[nid + 1 + sizeX + planeSize] += w111 * ix;
	    }
	
	}
    }

    //////////////////////////////////////////////////////////////////////////////////////////
    // splatPoint function for single array input with weighted normalization
    //////////////////////////////////////////////////////////////////////////////////////////
  
    inline void splatPoint(const float hx, const float hy, const float hz, const float ix ,float* outIm, float* outWt, int sizeX, int sizeY, int sizeZ){
      // this is a fast version of the floor function
    
      int fx = static_cast<int>(hx);
      int fy = static_cast<int>(hy);
      int fz = static_cast<int>(hz);
      size_t planeSize = sizeX * sizeY;
                
      if (hx < 0 && hx != static_cast<int>(hx)) --fx;
      if (hy < 0 && hy != static_cast<int>(hy)) --fy;
      if (hz < 0 && hz != static_cast<int>(hz)) --fz;

      if (hx > -1 && hx < (int) sizeX &&  // inside vol with 1 pix border
	  hy > -1 && hy < (int) sizeY &&
	  hz > -1 && hz < (int) sizeZ)
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
	  bool floorXIn = (fx >= 0), ceilXIn = (fx < (int) sizeX-1),
	    floorYIn = (fy >= 0), ceilYIn = (fy < (int) sizeY-1),
	    floorZIn = (fz >= 0), ceilZIn = (fz < (int) sizeZ-1);

	  size_t nid = (fz * sizeY + fy) * sizeX + fx;
	  if (floorXIn && floorYIn && floorZIn) {
	    outIm[nid]       += w000 * ix;
	    outWt[nid]       += w000;
	  }
	  if (floorXIn && floorYIn && ceilZIn){
	    outIm[nid + planeSize]  += w001 * ix;
	    outWt[nid + planeSize]  += w001;
	  }
	  if (floorXIn && ceilYIn && floorZIn){
	    outIm[nid + sizeX]  += w010 * ix;
	    outWt[nid + sizeX]  += w010;
	  }
	  if (floorXIn && ceilYIn && ceilZIn) {
	    outIm[nid + sizeX + planeSize] += w011 * ix;
	    outWt[nid + sizeX + planeSize] += w011;
	  }
	  if (ceilXIn && floorYIn && floorZIn) {
	    outIm[nid + 1] += w100 * ix;
	    outWt[nid + 1] += w100;
	  }
	  if (ceilXIn && floorYIn && ceilZIn) {
	    outIm[nid + 1 + planeSize] += w101 * ix;
	    outWt[nid + 1 + planeSize] += w101;
	  }
	  if (ceilXIn && ceilYIn && floorZIn)
	    {
	      outIm[nid + 1 + sizeX] += w110 * ix;
	      outWt[nid + 1 + sizeX] += w110;
	    }
	  if (ceilXIn && ceilYIn && ceilZIn)
	    {
	      outIm[nid + 1 + sizeX + planeSize] += w111 * ix;
	      outWt[nid + 1 + sizeX + planeSize] += w111;
	    }	
	}
    }
  }// end namespace Splatting
} // end namespace PyCA

#endif

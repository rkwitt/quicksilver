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

#ifndef __TRI_LERP_H
#define __TRI_LERP_H

#include <pycaConst.h>
#include <pycaUtils.h>

namespace PyCA {

//////////////////////////////////////////////////////////////////////////////////////////
// Get the pixel from 3D array 
//////////////////////////////////////////////////////////////////////////////////////////
// template<class T>
// inline __HOSTDEVICE__ T GetPixel(int x, const T* d_i, int sizeX)
// {
//     return d_i[x];
// }

// template<class T>
// inline __HOSTDEVICE__ T GetPixel(int x, int y, const T* d_i, int sizeX, int sizeY)
// {
//     int id = x + y * sizeX;
//     return d_i[id];
// }

// template<class T>
// inline __HOSTDEVICE__ T GetPixel(int x, int y, int z, const T* d_i, int sizeX, int sizeY, int sizeZ)
// {
//     int id = x + (y + z * sizeY) * sizeX;
//     return d_i[id];
// }

inline __HOSTDEVICE__ float getPixel3D(int x, int y, int z,
                                            const float* d_i,
                                            int sizeX, int sizeY, int sizeZ)
{
    int index = (z * sizeY + y) * sizeX + x;
    return d_i[index];
}

inline __HOSTDEVICE__ void getPixel3D(float& vx, float& vy, float& vz,
                                           int x, int y, int z,
                                           const float* d_iX, const float* d_iY, const float* d_iZ,
                                           int sizeX, int sizeY, int sizeZ)
{
    int index = (z * sizeY + y) * sizeX + x;

    vx = d_iX[index];
    vy = d_iY[index];
    vz = d_iZ[index];
}

//////////////////////////////////////////////////////////////////////////////////////////
// Wrap around strategy 
//////////////////////////////////////////////////////////////////////////////////////////
// Old version
// __HOSTDEVICE__ void wrap(int& r, int b){
//     if (r < 0) r += b;
//     else if (r > b) r %= b;
// }

inline __HOSTDEVICE__ int safe_mod(int r, int b){
    int m = r % b;
    return (m < 0) ? m + b : m;
}

inline __HOSTDEVICE__ void wrap(int& r, int b){
    r %= b;
    if (r < 0) {
        r += b;
    }
}

inline __HOSTDEVICE__ void wrapBackground(int& floorX, int& ceilX, int sizeX){
    wrap(floorX, sizeX);
    wrap(ceilX, sizeX);
}

inline __HOSTDEVICE__ void wrapBackground(int& floorX,int& floorY,int& floorZ,
                                               int& ceilX, int& ceilY, int& ceilZ,
                                               int  sizeX, int  sizeY, int  sizeZ){

    wrapBackground(floorX, ceilX, sizeX);
    wrapBackground(floorY, ceilY, sizeY);
    wrapBackground(floorZ, ceilZ, sizeZ);
}

//////////////////////////////////////////////////////////////////////////////////////////
// Clamp strategy
//////////////////////////////////////////////////////////////////////////////////////////
inline __HOSTDEVICE__ void clamp(int& r, int b){
    if (r < 0) r = 0;
    else if (r >= b) r = b - 1;
}
    
inline __HOSTDEVICE__ void clampBackground(int& floorX, int& ceilX, int sizeX){
    if(floorX < 0){
      floorX = 0;
      if(ceilX < 0) ceilX = 0;
    }

    if(ceilX >= sizeX){
      ceilX = sizeX-1;
      if(floorX >= sizeX) floorX = sizeX-1;
    }
}

inline __HOSTDEVICE__ void clampBackground(int& floorX,int& floorY,int& floorZ,
                                int& ceilX, int& ceilY, int& ceilZ,
                                int  sizeX, int  sizeY, int  sizeZ){

    clampBackground(floorX, ceilX, sizeX);
    clampBackground(floorY, ceilY, sizeY);
    clampBackground(floorZ, ceilZ, sizeZ);
}


//////////////////////////////////////////////////////////////////////////////////////////
// Check if the point is completely inside the boundary
//////////////////////////////////////////////////////////////////////////////////////////

inline __HOSTDEVICE__ bool isInside(int floorX,int floorY,int floorZ,
                                  int ceilX, int ceilY, int ceilZ,
                                  int  sizeX, int  sizeY, int  sizeZ){
    
    return (floorX >= 0 && ceilX < sizeX &&
            floorY >= 0 && ceilY < sizeY &&
            floorZ >= 0 && ceilZ < sizeZ);
}

//////////////////////////////////////////////////////////////////////////////////////////
// Trilerp function for Field input 
//////////////////////////////////////////////////////////////////////////////////////////
template<BackgroundStrategy backgroundStrategy>
inline  __HOSTDEVICE__ void triLerp(float& hx, float& hy, float& hz,
				    const float* imgX, const float* imgY, const float* imgZ,
				    float x, float y, float z,
				    int sizeX, int sizeY, int sizeZ){
    
    int floorX = (int)(x);
    int floorY = (int)(y);
    int floorZ = (int)(z);

    if (x < 0 && x != (int)(x)) --floorX;
    if (y < 0 && y != (int)(y)) --floorY;
    if (z < 0 && z != (int)(z)) --floorZ;

    // this is not truly ceiling, but floor + 1, which is usually ceiling    
    int ceilX = floorX + 1;
    int ceilY = floorY + 1;
    int ceilZ = floorZ + 1;

    float t = x - floorX;
	float u = y - floorY;
	float v = z - floorZ;

    float oneMinusT = 1.f - t;
	float oneMinusU = 1.f - u;
    float oneMinusV = 1.f - v;

    float v0X=0.f, v0Y=0.f, v0Z=0.f;
    float v1X=0.f, v1Y=0.f, v1Z=0.f;
    float v2X=0.f, v2Y=0.f, v2Z=0.f;
    float v3X=0.f, v3Y=0.f, v3Z=0.f;
    float v4X=0.f, v4Y=0.f, v4Z=0.f;
    float v5X=0.f, v5Y=0.f, v5Z=0.f;
    float v6X=0.f, v6Y=0.f, v6Z=0.f;
    float v7X=0.f, v7Y=0.f, v7Z=0.f;

    int inside = 1;

    // adjust the position of the sample point if required
    if (backgroundStrategy == BACKGROUND_STRATEGY_WRAP){
        wrapBackground(floorX, floorY, floorZ,
                       ceilX, ceilY, ceilZ,
                       sizeX, sizeY, sizeZ);
    }
    else if (backgroundStrategy == BACKGROUND_STRATEGY_CLAMP){
        clampBackground(floorX, floorY, floorZ,
                        ceilX, ceilY, ceilZ,
                        sizeX, sizeY, sizeZ);
    }
    else {
        inside = isInside(floorX, floorY, floorZ,
                          ceilX, ceilY, ceilZ,
                          sizeX, sizeY, sizeZ);
    }

    
    if (inside){
        getPixel3D(v0X, v0Y, v0Z, floorX, floorY, floorZ, imgX, imgY, imgZ, sizeX, sizeY, sizeZ);
        getPixel3D(v1X, v1Y, v1Z, ceilX, floorY, floorZ,  imgX, imgY, imgZ, sizeX, sizeY, sizeZ);
        getPixel3D(v2X, v2Y, v2Z, ceilX, ceilY, floorZ,  imgX, imgY, imgZ, sizeX, sizeY, sizeZ);
        getPixel3D(v3X, v3Y, v3Z, floorX, ceilY, floorZ,  imgX, imgY, imgZ, sizeX, sizeY, sizeZ);

        getPixel3D(v4X, v4Y, v4Z, floorX, ceilY, ceilZ,  imgX, imgY, imgZ, sizeX, sizeY, sizeZ);
        getPixel3D(v5X, v5Y, v5Z, ceilX, ceilY, ceilZ,  imgX, imgY, imgZ, sizeX, sizeY, sizeZ);
        getPixel3D(v6X, v6Y, v6Z, ceilX, floorY, ceilZ,  imgX, imgY, imgZ, sizeX, sizeY, sizeZ);
        getPixel3D(v7X, v7Y, v7Z, floorX, floorY, ceilZ,  imgX, imgY, imgZ, sizeX, sizeY, sizeZ);
    }else if (backgroundStrategy == BACKGROUND_STRATEGY_ID){
            //
            // coordinate is not inside volume, return identity
            //
            hx = x; hy = y; hz = z;
            return;
    } else if (backgroundStrategy == BACKGROUND_STRATEGY_ZERO) {
            hx = 0; hy = 0; hz = 0;
            return;
    } else if (backgroundStrategy == BACKGROUND_STRATEGY_PARTIAL_ID ||
               backgroundStrategy == BACKGROUND_STRATEGY_PARTIAL_ZERO)
        //
        // coordinate is not inside volume; initialize cube
        // corners to identity then set any corners of cube that
        // fall on volume boundary
    {
        if (backgroundStrategy == BACKGROUND_STRATEGY_PARTIAL_ID) {
            v0X = floorX; v0Y = floorY; v0Z = floorZ;
            v1X = ceilX;  v1Y = floorY; v1Z = floorZ;
            v2X = ceilX;  v2Y = ceilY;  v2Z = floorZ;
            v3X = floorX; v3Y = ceilY;  v3Z = floorZ;
            v4X = floorX; v4Y = ceilY;  v4Z = ceilZ;
            v5X = ceilX;  v5Y = ceilY;  v5Z = ceilZ;
            v6X = ceilX;  v6Y = floorY; v6Z = ceilZ;
            v7X = floorX; v7Y = floorY; v7Z = ceilZ;
        } else {
             v0X = 0; v0Y = 0; v0Z = 0;
             v1X = 0; v1Y = 0; v1Z = 0;
             v2X = 0; v2Y = 0; v2Z = 0;
             v3X = 0; v3Y = 0; v3Z = 0;
             v4X = 0; v4Y = 0; v4Z = 0;
             v5X = 0; v5Y = 0; v5Z = 0;
             v6X = 0; v6Y = 0; v6Z = 0;
             v7X = 0; v7Y = 0; v7Z = 0;
        }
        
        bool floorXIn = (floorX >= 0) && (floorX < sizeX);
        bool floorYIn = (floorY >= 0) && (floorY < sizeY);
        bool floorZIn = (floorZ >= 0) && (floorZ < sizeZ);

        bool ceilXIn = (ceilX >= 0) && (ceilX < sizeX);
        bool ceilYIn = (ceilY >= 0) && (ceilY < sizeY);
        bool ceilZIn = (ceilZ >= 0) && (ceilZ < sizeZ);

        if (floorXIn && floorYIn && floorZIn)
            getPixel3D(v0X, v0Y, v0Z, floorX, floorY, floorZ,
                       imgX, imgY, imgZ, sizeX, sizeY, sizeZ);
        if (ceilXIn && floorYIn && floorZIn)
            getPixel3D(v1X, v1Y, v1Z, ceilX, floorY, floorZ,
                       imgX, imgY, imgZ, sizeX, sizeY, sizeZ);
        if (ceilXIn && ceilYIn && floorZIn)
            getPixel3D(v2X, v2Y, v2Z,ceilX, ceilY, floorZ,
                       imgX, imgY, imgZ, sizeX, sizeY, sizeZ);

        if (floorXIn && ceilYIn && floorZIn)
            getPixel3D(v3X, v3Y, v3Z, floorX, ceilY, floorZ,
                       imgX, imgY, imgZ, sizeX, sizeY, sizeZ);

        if (floorXIn && ceilYIn && ceilZIn)
            getPixel3D(v4X, v4Y, v4Z, floorX, ceilY, ceilZ,
                       imgX, imgY, imgZ, sizeX, sizeY, sizeZ);

        if (ceilXIn && ceilYIn && ceilZIn)
            getPixel3D(v5X, v5Y, v5Z, ceilX, ceilY, ceilZ,
                       imgX, imgY, imgZ, sizeX, sizeY, sizeZ);

        if (ceilXIn && floorYIn && ceilZIn)
            getPixel3D(v6X, v6Y, v6Z, ceilX, floorY, ceilZ,
                       imgX, imgY, imgZ, sizeX, sizeY, sizeZ);

        if (floorXIn && floorYIn && ceilZIn)
            getPixel3D(v7X, v7Y, v7Z, floorX, floorY, ceilZ,
                       imgX, imgY, imgZ, sizeX, sizeY, sizeZ);
    }
    
    //
    // do trilinear interpolation
    //
    
    //
    // this is the basic trilerp function...
    //
    //     h = 
    //       v0 * (1 - t) * (1 - u) * (1 - v) +
    //       v1 * t       * (1 - u) * (1 - v) +
    //       v2 * t       * u       * (1 - v) +
    //       v3 * (1 - t) * u       * (1 - v) +
    //       v4 * (1 - t) * u       * v       +
    //       v5 * t       * u       * v       +
    //       v6 * t       * (1 - u) * v       +
    //       v7 * (1 - t) * (1 - u) * v;
    //
    // the following nested version saves 30 multiplies.
    //
    hx = 
        oneMinusT * (oneMinusU * (v0X * oneMinusV + v7X * v)  +
                     u         * (v3X * oneMinusV + v4X * v)) +
        t         * (oneMinusU * (v1X * oneMinusV + v6X * v)  +
                     u         * (v2X * oneMinusV + v5X * v));
    
    hy = 
        oneMinusT * (oneMinusU * (v0Y * oneMinusV + v7Y * v)  +
                     u         * (v3Y * oneMinusV + v4Y * v)) +
        t         * (oneMinusU * (v1Y * oneMinusV + v6Y * v)  +
                     u         * (v2Y * oneMinusV + v5Y * v));
    
    hz = 
        oneMinusT * (oneMinusU * (v0Z * oneMinusV + v7Z * v)  +
                     u         * (v3Z * oneMinusV + v4Z * v)) +
        t         * (oneMinusU * (v1Z * oneMinusV + v6Z * v)  +
                     u         * (v2Z * oneMinusV + v5Z * v));

}

//////////////////////////////////////////////////////////////////////////////////////////
// Trilerp function for single array input 
//////////////////////////////////////////////////////////////////////////////////////////

template<BackgroundStrategy backgroundStrategy>
inline __HOSTDEVICE__ 
float 
triLerp(const float* img,
	float x, float y, float z,
	int sizeX, int sizeY, int sizeZ,
	float background = 0.f)
{
    int floorX = (int)(x);
    int floorY = (int)(y);
    int floorZ = (int)(z);

    if (x < 0 && x != (int)(x)) --floorX;
    if (y < 0 && y != (int)(y)) --floorY;
    if (z < 0 && z != (int)(z)) --floorZ;

    // this is not truly ceiling, but floor + 1, which is usually ceiling    
    int ceilX = floorX + 1;
    int ceilY = floorY + 1;
    int ceilZ = floorZ + 1;

    float t = x - floorX;
	float u = y - floorY;
	float v = z - floorZ;

    float oneMinusT = 1.f- t;
	float oneMinusU = 1.f- u;
    float oneMinusV = 1.f- v;

    float v0, v1, v2, v3, v4, v5, v6, v7;
    int inside = 1;

    // adjust the position of the sample point if required
    if (backgroundStrategy == BACKGROUND_STRATEGY_WRAP){
        wrapBackground(floorX, floorY, floorZ,
                       ceilX, ceilY, ceilZ,
                       sizeX, sizeY, sizeZ);
    }
    else if (backgroundStrategy == BACKGROUND_STRATEGY_CLAMP){
        clampBackground(floorX, floorY, floorZ,
                        ceilX, ceilY, ceilZ,
                        sizeX, sizeY, sizeZ);
    }
    else if (backgroundStrategy == BACKGROUND_STRATEGY_VAL || 
	     backgroundStrategy == BACKGROUND_STRATEGY_ZERO || 
	     backgroundStrategy == BACKGROUND_STRATEGY_PARTIAL_ZERO){

	if(backgroundStrategy == BACKGROUND_STRATEGY_ZERO || 
	   backgroundStrategy == BACKGROUND_STRATEGY_PARTIAL_ZERO){
	    background = 0.f;
	}

        inside = isInside(floorX, floorY, floorZ,
                          ceilX, ceilY, ceilZ,
                          sizeX, sizeY, sizeZ);
    }else{
	// unknown background strategy, don't allow compilation
	STATIC_ASSERT(backgroundStrategy== BACKGROUND_STRATEGY_WRAP ||
		      backgroundStrategy== BACKGROUND_STRATEGY_CLAMP ||
		      backgroundStrategy== BACKGROUND_STRATEGY_ZERO ||
		      backgroundStrategy== BACKGROUND_STRATEGY_PARTIAL_ZERO ||
		      backgroundStrategy== BACKGROUND_STRATEGY_VAL);
	return 0.f;
    }

    
    if (inside){
        v0 = getPixel3D(floorX, floorY, floorZ, img, sizeX, sizeY, sizeZ);
        v1 = getPixel3D(ceilX, floorY, floorZ,  img, sizeX, sizeY, sizeZ);
        v2 = getPixel3D(ceilX, ceilY, floorZ,   img, sizeX, sizeY, sizeZ);
        v3 = getPixel3D(floorX, ceilY, floorZ,  img, sizeX, sizeY, sizeZ);

        v4 = getPixel3D(floorX, ceilY, ceilZ,  img, sizeX, sizeY, sizeZ);
        v5 = getPixel3D(ceilX, ceilY, ceilZ,   img, sizeX, sizeY, sizeZ);
        v6 = getPixel3D(ceilX, floorY, ceilZ,  img, sizeX, sizeY, sizeZ);
        v7 = getPixel3D(floorX, floorY, ceilZ, img, sizeX, sizeY, sizeZ);
    }else {
        bool floorXIn = floorX >= 0 && floorX < sizeX;
        bool floorYIn = floorY >= 0 && floorY < sizeY;
        bool floorZIn = floorZ >= 0 && floorZ < sizeZ;
      
        bool ceilXIn = ceilX >= 0 && ceilX < sizeX;
        bool ceilYIn = ceilY >= 0 && ceilY < sizeY;
        bool ceilZIn = ceilZ >= 0 && ceilZ < sizeZ;

        v0 = (floorXIn && floorYIn && floorZIn) ? getPixel3D(floorX, floorY, floorZ, img, sizeX, sizeY, sizeZ): background;
        v1 = (ceilXIn && floorYIn && floorZIn)  ? getPixel3D(ceilX, floorY, floorZ,  img, sizeX, sizeY, sizeZ): background;
        v2 = (ceilXIn && ceilYIn && floorZIn)   ? getPixel3D(ceilX, ceilY, floorZ,   img, sizeX, sizeY, sizeZ): background;
        v3 = (floorXIn && ceilYIn && floorZIn)  ? getPixel3D(floorX, ceilY, floorZ,  img, sizeX, sizeY, sizeZ): background;
        
        v4 = (floorXIn && ceilYIn && ceilZIn)   ? getPixel3D(floorX, ceilY, ceilZ,  img, sizeX, sizeY, sizeZ): background;
        v5 = (ceilXIn && ceilYIn && ceilZIn)    ? getPixel3D(ceilX, ceilY, ceilZ,   img, sizeX, sizeY, sizeZ): background;
        v6 = (ceilXIn && floorYIn && ceilZIn)   ? getPixel3D(ceilX, floorY, ceilZ,  img, sizeX, sizeY, sizeZ): background;
        v7 = (floorXIn && floorYIn && ceilZIn)  ? getPixel3D(floorX, floorY, ceilZ, img, sizeX, sizeY, sizeZ): background;
    }

    //
    // do trilinear interpolation
    //
    
    //
    // this is the basic trilerp function...
    //
    //     h = 
    //       v0 * (1 - t) * (1 - u) * (1 - v) +
    //       v1 * t       * (1 - u) * (1 - v) +
    //       v2 * t       * u       * (1 - v) +
    //       v3 * (1 - t) * u       * (1 - v) +
    //       v4 * (1 - t) * u       * v       +
    //       v5 * t       * u       * v       +
    //       v6 * t       * (1 - u) * v       +
    //       v7 * (1 - t) * (1 - u) * v;
    //
    // the following nested version saves 30 multiplies.
    //
    return  oneMinusT * (oneMinusU * (v0 * oneMinusV + v7 * v)  +
                         u         * (v3 * oneMinusV + v4 * v)) +
        t         * (oneMinusU * (v1 * oneMinusV + v6 * v)  +
                     u         * (v2 * oneMinusV + v5 * v));
}

//
// TRI-CUBIC INTERPOLATION
//

//
// Cubic Spline interpolation w.r.to Catmull-Rom
//
inline __HOSTDEVICE__
float 
cubicsplineinterp(float *cubestart, float x,float y, float z)
{
	float xcube = x*x*x;
	float xsquare = x*x;
	float ycube = y*y*y;
	float ysquare = y*y;
	float zcube = z*z*z;
	float zsquare = z*z;
// along X-direction	
	float v[4][4];	
	for (int r=0;r<4;r++){
		for (int c=0;c<4;c++){
			float x0 = *(cubestart+16*r+4*c);
			float x1 = *(cubestart+16*r+4*c+1);
			float x2 = *(cubestart+16*r+4*c+2);
               float x3 = *(cubestart+16*r+4*c+3);
			float a0 = x3 - x2 - x0 + x1;
			float a1 = x0 - x1 - a0;
			float a2 = x2 - x0;
			float a3 = x1;
			v[r][c] = a0 * xcube + a1 * xsquare + a2 * x + a3;  
		}
	}
// along Y-direction
	float u[4];	
	for (int r=0;r<4;r++){
		float y0 = v[r][0];
		float y1 = v[r][1];
		float y2 = v[r][2];
		float y3 = v[r][3];
		float a0 = y3 - y2 - y0 + y1;
		float a1 = y0 - y1 - a0;
		float a2 = y2 - y0;
		float a3 = y1;
 		u[r] = a0 * ycube + a1 * ysquare + a2 * y + a3;
	}
// along Z-direction
		float z0 = u[0];
		float z1 = u[1];
		float z2 = u[2];
		float z3 = u[3];
		float a0 = z3 - z2 - z0 + z1;
		float a1 = z0 - z1 - a0;
		float a2 = z2 - z0;
		float a3 = z1;
 		float retval = a0 * zcube + a1 * zsquare + a2 * z + a3;
		//if (retval < 0) retval = 0; // making negative values = 0
	return(retval);
	
}

template<BackgroundStrategy backgroundStrategy>
inline __HOSTDEVICE__ 
float 
triCubic(const float* img,
	 float x, float y, float z,
	 int sizeX, int sizeY, int sizeZ,
	 float background=0.f)
{
  int floorX = static_cast<int>(x);
  int floorY = static_cast<int>(y);
  int floorZ = static_cast<int>(z);
  if (x < 0 && x != static_cast<int>(x)) --floorX;
  if (y < 0 && y != static_cast<int>(y)) --floorY;
  if (z < 0 && z != static_cast<int>(z)) --floorZ;
		
  const float t = x - floorX;
  const float u = y - floorY;
  const float v = z - floorZ; 

  // array of coefficients
  float sm_array[4][4][4];

  for(int sm_a=0;sm_a<4;sm_a++){
    for(int sm_b=0;sm_b<4;sm_b++){
      for(int sm_c=0;sm_c<4;sm_c++){
	// current image co-ordinates of an image
	int arrayX=floorX+sm_a-1;
	int arrayY=floorY+sm_b-1;
	int arrayZ=floorZ+sm_c-1;
	// checking boundary condition
	sm_array[sm_a][sm_b][sm_c] = 
	    getSafeVal<float, backgroundStrategy>(img, 
						  sizeX, sizeY, sizeZ,
						  arrayX, arrayY, arrayZ,
						  background);
      }
    }
  }

  float retval =  cubicsplineinterp(&sm_array[0][0][0],v,u,t);

  // can't throw an exception from device function
  // if (retval != retval)
  //   {
  // 	throw std::runtime_error("Found a NAN in cubicSplineinterp.");
  //   } 	

  return retval;
}

// NEAREST NEIGHBOR INTERP
template<BackgroundStrategy backgroundStrategy>
inline __HOSTDEVICE__
float 
nearestNeighbor(const float *img,
		float x, float y, float z,
		int sizeX, int sizeY, int sizeZ,
		float background=0.0)
{
    int xIndex = x > 0 ? static_cast<int>(x + 0.5) : static_cast<int>(x - 0.5);
    int yIndex = y > 0 ? static_cast<int>(y + 0.5) : static_cast<int>(y - 0.5);
    int zIndex = z > 0 ? static_cast<int>(z + 0.5) : static_cast<int>(z - 0.5);

    return getSafeVal<float, backgroundStrategy>(img, 
						 sizeX, sizeY, sizeZ,
						 xIndex, yIndex, zIndex,
						 background);
}


template <InterpT interp_method,
	  BackgroundStrategy backgroundStrategy>
inline __HOSTDEVICE__
float 
point_interp(const float *img,
	     float x, float y, float z,
	     int sizeX, int sizeY, int sizeZ,
	     float background=0.0)
{
    if(interp_method == INTERP_NN){
	return 
	    nearestNeighbor<backgroundStrategy>
	    (img, x, y, z, sizeX, sizeY, sizeZ, background);
    }else if(interp_method == INTERP_LINEAR){
	return 
	    triLerp<backgroundStrategy>
	    (img, x, y, z, sizeX, sizeY, sizeZ, background);
    }else if(interp_method == INTERP_CUBIC){
	return 
	    triCubic<backgroundStrategy>
	    (img, x, y, z, sizeX, sizeY, sizeZ, background);
    }else{
	// unknown interp method, don't allow to compile
	STATIC_ASSERT(interp_method == INTERP_NN ||
		      interp_method == INTERP_LINEAR ||
		      interp_method == INTERP_CUBIC);
	return 0.f;
    }
}

template <BackgroundStrategy backgroundStrategy>
inline __HOSTDEVICE__
float 
point_interp(const float *img,
       float x, float y, float z,
       int sizeX, int sizeY, int sizeZ,
       float background=0.0)
{
	return 
	    point_interp<DEFAULT_INTERP_METHOD, backgroundStrategy>
	    (img, x, y, z, sizeX, sizeY, sizeZ, background);
}



} // end namespace PyCA

#endif

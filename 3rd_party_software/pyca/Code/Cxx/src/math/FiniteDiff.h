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

#ifndef __FINITEDIFF_H
#define __FINITEDIFF_H

#include "pycaUtils.h"

namespace PyCA {

// NOTE: no scaling is done here
template < class T, enum DimT dim, enum DiffT diffType, BoundaryCondT bc>
__HOSTDEVICE__ inline T finiteDiff(const T *in, int x, int y, int z, int sx, int sy, int sz)
{
    if(x >= sx || y >= sy || z >= sz)
       return 0;

    int xp, yp, zp, xn, yn, zn;
    xp = xn = x; 
    yp = yn = y; 
    zp = zn = z;

    if(diffType == DIFF_FORWARD || diffType == DIFF_CENTRAL)
    {
        if(dim == DIM_X)
           xn++;
        else if(dim == DIM_Y)
           yn++;
        else if(dim == DIM_Z)
           zn++;
    }
    if(diffType == DIFF_BACKWARD || diffType == DIFF_CENTRAL)
    {
        if(dim == DIM_X)
           xp--;
        else if(dim == DIM_Y)
           yp--;
        else if(dim == DIM_Z)
           zp--;
    }

    // this is all compiled out for non-BC_APPROX versions
    bool bcApproxEdge = 
       (diffType == DIFF_CENTRAL && bc == BC_APPROX
	&& ((dim == DIM_X && (x == 0 || x == sx-1)) ||
	    (dim == DIM_Y && (y == 0 || y == sy-1)) ||
	    (dim == DIM_Z && (z == 0 || z == sz-1))));

    // avaiable at compile time
    const BackgroundStrategy bs = (bc == BC_WRAP ? 
				   BACKGROUND_STRATEGY_WRAP : 
				   BACKGROUND_STRATEGY_CLAMP);
       
    if (diffType == DIFF_CENTRAL && !bcApproxEdge){
       // scale by two
       return 0.5*(getSafeVal<T,bs>(in,sx,sy,sz,xn,yn,zn)-
                    getSafeVal<T,bs>(in,sx,sy,sz,xp,yp,zp));
    }else{
       return getSafeVal<T,bs>(in,sx,sy,sz,xn,yn,zn)-
	  getSafeVal<T,bs>(in,sx,sy,sz,xp,yp,zp);
    }
}

// Finite difference of a product of two arrays
// This is a bit more precise than expanding using the product rule
template < class T, enum DimT dim, enum DiffT diffType, BoundaryCondT bc>
__HOSTDEVICE__ inline T finiteDiffProd(const T *X, const T *Y, int x, int y, int z, int sx, int sy, int sz)
{
    if(x >= sx && y >= sy && z >= sz)
       return 0;

    int xp, yp, zp, xn, yn, zn;
    xp = xn = x; 
    yp = yn = y; 
    zp = zn = z;

    if(diffType == DIFF_FORWARD || diffType == DIFF_CENTRAL)
    {
        if(dim == DIM_X)
           xn++;
        else if(dim == DIM_Y)
           yn++;
        else if(dim == DIM_Z)
           zn++;
    }
    if(diffType == DIFF_BACKWARD || diffType == DIFF_CENTRAL)
    {
        if(dim == DIM_X)
           xp--;
        else if(dim == DIM_Y)
           yp--;
        else if(dim == DIM_Z)
           zp--;
    }

    // this is all compiled out for non-BC_APPROX versions
    // Note: I think this is right, but have't tested -- jacob
    // you might want to have a look (jsp 12/2012)
    bool bcApproxEdge = 
       (diffType == DIFF_CENTRAL && bc == BC_APPROX
	&& ((dim == DIM_X && (x == 0 || x == sx-1)) ||
	    (dim == DIM_Y && (y == 0 || y == sy-1)) ||
	    (dim == DIM_Z && (z == 0 || z == sz-1))));

    // avaiable at compile time
    const BackgroundStrategy bs = (bc == BC_WRAP ? 
				   BACKGROUND_STRATEGY_WRAP : 
				   BACKGROUND_STRATEGY_CLAMP);

    if (diffType == DIFF_CENTRAL && !bcApproxEdge) // scale by two
        return 0.5*(getSafeVal<T,bs>(X,sx,sy,sz,xn,yn,zn)*
                    getSafeVal<T,bs>(Y,sx,sy,sz,xn,yn,zn)
                    -
                    getSafeVal<T,bs>(X,sx,sy,sz,xp,yp,zp)*
                    getSafeVal<T,bs>(Y,sx,sy,sz,xp,yp,zp));
    else
        return getSafeVal<T,bs>(X,sx,sy,sz,xn,yn,zn)*
	    getSafeVal<T,bs>(Y,sx,sy,sz,xn,yn,zn)
	    -
	    getSafeVal<T,bs>(X,sx,sy,sz,xp,yp,zp)*
	    getSafeVal<T,bs>(Y,sx,sy,sz,xp,yp,zp);
}

// Masked finite difference -- clamp/approx at boundary of mask.  Wrapping not supported.
template < class T, enum DimT dim, enum DiffT diffType, BoundaryCondT bc>
__HOSTDEVICE__ inline T finiteDiffMask(const T *in, const T*mask, int x, int y, int z, int sx, int sy, int sz)
{
    if(x >= sx && y >= sy && z >= sz)
       return 0;

    int xp, yp, zp, xn, yn, zn;
    xp = xn = x; 
    yp = yn = y; 
    zp = zn = z;

    // zero diff if this pixel is masked
    if(getSafeVal<T,BACKGROUND_STRATEGY_CLAMP>(mask,sx,sy,sz,x,y,z) == (T)0)
       return (T)0;

    if(diffType == DIFF_FORWARD || diffType == DIFF_CENTRAL)
    {
        if(dim == DIM_X)
           xn++;
        else if(dim == DIM_Y)
           yn++;
        else if(dim == DIM_Z)
           zn++;
    }
    if(diffType == DIFF_BACKWARD || diffType == DIFF_CENTRAL)
    {
        if(dim == DIM_X)
           xp--;
        else if(dim == DIM_Y)
           yp--;
        else if(dim == DIM_Z)
           zp--;
    }

    // mask values at prev/next points
    T mp = getSafeVal<T, BACKGROUND_STRATEGY_CLAMP>(mask,sx,sy,sz,xp,yp,zp);
    T mn = getSafeVal<T, BACKGROUND_STRATEGY_CLAMP>(mask,sx,sy,sz,xn,yn,zn);
    
    T vp = getSafeVal<T,BACKGROUND_STRATEGY_CLAMP>(in,sx,sy,sz,xp,yp,zp);
    T v = getSafeVal<T,BACKGROUND_STRATEGY_CLAMP>(in,sx,sy,sz,x,y,z);
    T vn = getSafeVal<T,BACKGROUND_STRATEGY_CLAMP>(in,sx,sy,sz,xn,yn,zn);

    if(mn == (T)0) vn = v;
    if(mp == (T)0) vp = v;
    
    // this is all compiled out for non-BC_APPROX versions.
    // Not sure this is implemented correctly.
    bool bcApproxEdge = 
       (diffType == DIFF_CENTRAL && bc == BC_APPROX
	&& 
	( // border pixel
	 ((dim == DIM_X && (x == 0 || x == sx-1)) ||
	  (dim == DIM_Y && (y == 0 || y == sy-1)) ||
	  (dim == DIM_Z && (z == 0 || z == sz-1)))
	 ||
	 // or borders a masked pixel
	 (mn == (T)0 || mp == (T)0)
	 )
	);
       
    if (diffType == DIFF_CENTRAL && !bcApproxEdge){
       return 0.5*(vn-vp);
    }else{
       return vn-vp;
    }
}

// NOTE: no scaling is done here
template < class T, enum DiffT diffType, BoundaryCondT bc>
__HOSTDEVICE__ inline void gradientPoint(const T *in, int x, int y, int z, float& gx, float& gy, float& gz, int sx, int sy, int sz)
{
   gx = finiteDiff<T,DIM_X,diffType,bc>(in,x,y,z,sx,sy,sz);
   gy = finiteDiff<T,DIM_Y,diffType,bc>(in,x,y,z,sx,sy,sz);
   gz = finiteDiff<T,DIM_Z,diffType,bc>(in,x,y,z,sx,sy,sz);
}

template < class T, enum DiffT diffType, BoundaryCondT bc>
__HOSTDEVICE__ inline void jacobianPoint(const T *vx, const T *vy, const T *vz, int x, int y, int z, float& Jxx, float& Jxy, float& Jxz,
          float& Jyx, float& Jyy, float& Jyz,
          float& Jzx, float& Jzy, float& Jzz,
          int sx, int sy, int sz)
{ // Full jacobianPoint matrix computation
   gradientPoint<T,diffType,bc>(vx,x,y,z,Jxx,Jxy,Jxz,sx,sy,sz);
   gradientPoint<T,diffType,bc>(vy,x,y,z,Jyx,Jyy,Jyz,sx,sy,sz);
   gradientPoint<T,diffType,bc>(vz,x,y,z,Jzx,Jzy,Jzz,sx,sy,sz);
}

template < class T, enum DiffT diffType, BoundaryCondT bc>
__HOSTDEVICE__ inline void divtensorprodPoint(const T *Xx, const T *Xy, const T *Xz, const T *Yx, const T *Yy, const T *Yz, int x, int y, int z, float& ox, float& oy, float& oz, int sx, int sy, int sz)
{ // take the (row-wise) divergence of a tensor product X\otimes Y
   ox = finiteDiffProd<T,DIM_X,diffType,bc>(Xx,Yx,x,y,z,sx,sy,sz)
      +finiteDiffProd<T,DIM_Y,diffType,bc>(Xx,Yy,x,y,z,sx,sy,sz)
      +finiteDiffProd<T,DIM_Z,diffType,bc>(Xx,Yz,x,y,z,sx,sy,sz);
   oy = finiteDiffProd<T,DIM_X,diffType,bc>(Xy,Yx,x,y,z,sx,sy,sz)
      +finiteDiffProd<T,DIM_Y,diffType,bc>(Xy,Yy,x,y,z,sx,sy,sz)
      +finiteDiffProd<T,DIM_Z,diffType,bc>(Xy,Yz,x,y,z,sx,sy,sz);
   oz = finiteDiffProd<T,DIM_X,diffType,bc>(Xz,Yx,x,y,z,sx,sy,sz)
      +finiteDiffProd<T,DIM_Y,diffType,bc>(Xz,Yy,x,y,z,sx,sy,sz)
      +finiteDiffProd<T,DIM_Z,diffType,bc>(Xz,Yz,x,y,z,sx,sy,sz);
}

// scaled versions
template < class T, enum DiffT diffType, BoundaryCondT bc>
__HOSTDEVICE__ inline void gradientPoint(const T *in, int x, int y, int z, float& gx, float& gy, float& gz, int sx, int sy, int sz,float scalex, float scaley, float scalez)
{
   gradientPoint<T,diffType,bc>(in,x,y,z,gx,gy,gz,sx,sy,sz);
   gx *= scalex;
   gy *= scaley;
   gz *= scalez;
}

template < class T, enum DiffT diffType, BoundaryCondT bc>
__HOSTDEVICE__ inline void jacobianPoint(const T *vx, const T *vy, const T *vz, int x, int y, int z, float& Jxx, float& Jxy, float& Jxz,
          float& Jyx, float& Jyy, float& Jyz,
          float& Jzx, float& Jzy, float& Jzz,
          int sx, int sy, int sz,
          float scalex, float scaley, float scalez)
{
   gradientPoint<T,diffType,bc>(vx,x,y,z,Jxx,Jxy,Jxz,sx,sy,sz,scalex,scaley,scalez);
   gradientPoint<T,diffType,bc>(vy,x,y,z,Jyx,Jyy,Jyz,sx,sy,sz,scalex,scaley,scalez);
   gradientPoint<T,diffType,bc>(vz,x,y,z,Jzx,Jzy,Jzz,sx,sy,sz,scalex,scaley,scalez);
}

template < class T, enum DiffT diffType, BoundaryCondT bc, bool slice>
__HOSTDEVICE__ inline void jacDetPoint(T &jdet, 
				       const T *vx, const T *vy, const T *vz, 
				       int x, int y, int z,
				       int sx, int sy, int sz,
				       float scalex=1.f, float scaley=1.f, float scalez=1.f)
{ 
   T Jxx, Jxy, Jxz, Jyx, Jyy, Jyz, Jzx, Jzy, Jzz;
   gradientPoint<T,diffType,bc>(vx,x,y,z,Jxx,Jxy,Jxz,sx,sy,sz,scalex,scaley,scalez);
   gradientPoint<T,diffType,bc>(vy,x,y,z,Jyx,Jyy,Jyz,sx,sy,sz,scalex,scaley,scalez);
   if(slice){
      jdet = det(Jxx, Jxy,
		 Jyx, Jyy);
   }else{
      gradientPoint<T,diffType,bc>(vz,x,y,z,Jzx,Jzy,Jzz,sx,sy,sz,scalex,scaley,scalez);
      jdet = det(Jxx, Jxy, Jxz,
		 Jyx, Jyy, Jyz,
		 Jzx, Jzy, Jzz);
   }
   
}

template < class T, enum DiffT diffType, BoundaryCondT bc>
__HOSTDEVICE__ 
inline 
void 
gradientPointMask(const T *in, const T *mask, 
		  int x, int y, int z, 
		  float& gx, float& gy, float& gz, 
		  int sx, int sy, int sz,
		  float scalex, float scaley, float scalez)
{
   gx = finiteDiffMask<T,DIM_X,diffType,bc>(in,mask,x,y,z,sx,sy,sz);
   gy = finiteDiffMask<T,DIM_Y,diffType,bc>(in,mask,x,y,z,sx,sy,sz);
   gz = finiteDiffMask<T,DIM_Z,diffType,bc>(in,mask,x,y,z,sx,sy,sz);
   gx *= scalex;
   gy *= scaley;
   gz *= scalez;
}


}

#endif // __FINITEDIFF_H

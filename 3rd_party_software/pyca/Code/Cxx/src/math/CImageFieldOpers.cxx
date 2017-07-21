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

#include "CImageFieldOpers.h"
#include <vector>
#include "interp.h"
#include "CSplat.h"
#include <algorithm>
#include "pycaUtils.h"
#include "CFiniteDiff.h"
#include "FiniteDiff.h"
#include <Image3D.h>
#include <Field3D.h>
#include <MemOpers.h>

namespace PyCA {

/////////////////////////////////////////////////////////////////////////////
// apply hField to an image
// defImage(x) = image(h(x))
/////////////////////////////////////////////////////////////////////////////
template<BackgroundStrategy bg, InterpT interp>
void ApplyH(float* a_o, const float* a_i, const Field3D& h, StreamT s){
    Vec3Di size = h.size();
    
    size_t id = 0;
    for (int k=0; k < size.z; ++k)
        for (int j=0; j < size.y; ++j)
            for (int i=0; i < size.x; ++i, ++id) {
                float x = h.x[id], y = h.y[id], z = h.z[id];
                
                a_o[id] = 
		    point_interp<interp, bg>
		    (a_i,
		     x, y, z,
		     size.x, size.y, size.z);
            }
}

template<BackgroundStrategy bg, InterpT interp>
void CImageFieldOpers::
ApplyH(Image3D& a_o, const Image3D& a_i, const Field3D& h, StreamT s){
    MK_CHECK_IMAGE_BACKGROUND(bg);
    MK_CHECK3_SIZE(a_o, a_i, h);

    PyCA::ApplyH<bg, interp>(a_o.get(), a_i.get(), h, s);
}

template<bool fwd, BackgroundStrategy bg, InterpT interp>
void ApplyV(float *a_o, const float* a_i, const Field3D& u,
            float delta, StreamT s){

    MK_CHECK_IMAGE_BACKGROUND(bg);
    Vec3Di size = u.size();
    Vec3Df sp   = u.spacing();
    Vec3Df iSp  = Vec3Df(1.f/sp.x, 1.f/sp.y, 1.f/sp.z);

    size_t id = 0;
    for (int k=0; k < size.z; ++k)
        for (int j=0; j < size.y; ++j)
            for (int i=0; i < size.x; ++i, ++id) {
                float x,y,z;
                if (fwd) {
                    x = i + delta * iSp.x * u.x[id];
                    y = j + delta * iSp.y * u.y[id];
                    z = k + delta * iSp.z * u.z[id];
                } else {
                    x = i - delta * iSp.x * u.x[id];
                    y = j - delta * iSp.y * u.y[id];
                    z = k - delta * iSp.z * u.z[id];
                }

                a_o[id] = 
		    point_interp<interp, bg>
		    (a_i, x, y, z,
		     size.x, size.y, size.z);
            }
}

/////////////////////////////////////////////////////////////////////////////
// apply uField to an image
// defImage(x) = image(x + delta * u(x))
/////////////////////////////////////////////////////////////////////////////

template<BackgroundStrategy bg, InterpT interp>
void CImageFieldOpers::
ApplyV(Image3D& a_o, const Image3D& a_i, const Field3D& u, const float& delta, StreamT s, bool onDev)
{
    MK_CHECK3_SIZE(a_o, a_i, u);
    
    PyCA::ApplyV<true, bg, interp>(a_o.get(), a_i.get(), u, delta, s);
}

/////////////////////////////////////////////////////////////////////////////
// apply uField to an image
// defImage(x) = image(x - delta * u(x))
/////////////////////////////////////////////////////////////////////////////

template<BackgroundStrategy bg, InterpT interp>
void CImageFieldOpers::
ApplyVInv(Image3D& a_o, const Image3D& a_i, const Field3D& u, const float& delta, StreamT s, bool onDev)
{
    MK_CHECK3_SIZE(a_o, a_i, u);
    PyCA::ApplyV<false, bg, interp>(a_o.get(), a_i.get(), u, delta, s);
};

template<BackgroundStrategy bg, InterpT interp>
void ComposeTranslation(float* a_o, const float* a_i, const Vec3Df& t,
                        const Vec3Di& size, StreamT s)
{
    size_t id = 0;
    for (int k=0; k < size.z; ++k)
        for (int j=0; j < size.y; ++j)
            for (int i=0; i < size.x; ++i, ++id) {
                
                float x = i + t.x;
                float y = j + t.y;
                float z = k + t.z;

                a_o[id] = 
		    point_interp<interp, bg>
		    (a_i, x, y, z,
		     size.x, size.y, size.z);
            }
}


template<BackgroundStrategy bg, InterpT interp>
void CImageFieldOpers::
ComposeTranslation(Image3D& a_o, const Image3D& a_i, const Vec3Df& t, StreamT s, bool onDev){
    MK_CHECK_IMAGE_BACKGROUND(bg);
    MK_CHECK2_SIZE(a_o, a_i);
    PyCA::ComposeTranslation<bg, interp>
	(a_o.get(), a_i.get(), t, a_o.size(), s);
}

void CImageFieldOpers
::Splat(Image3D& a_o, const Field3D& a_h, const Image3D& a_i,
               bool normalize, StreamT stream)
{
   Splat(a_o.get(), a_h, a_i.get(), normalize, stream);
}

void CImageFieldOpers
::Splat(float* a_o, const Field3D& a_h, const float* a_i,
               bool normalize, StreamT stream)
{
    if (normalize){
      Vec3Di size = a_h.size();
      size_t nVox = size.prod();
      float* a_w = new float[nVox];
      Splat(a_o, a_w, a_h, a_i, stream);
    }
    else{
      Splat(a_o, a_h, a_i, stream);
    }
}

void CImageFieldOpers::
Splat(float* a_o, float* a_w, const Field3D& a_h, const float* a_i, StreamT stream){
    
    Vec3Di size = a_h.size();
    size_t nVox = size.prod();

    size_t id = 0;

    CMemOpers<float>::SetMem(a_o, 0.f, nVox, stream, false);
    CMemOpers<float>::SetMem(a_w, 0.f, nVox, stream, false);
        
    for (int z = 0; z < size.z; ++z) {
        for (int y = 0; y < size.y; ++y) {
            for (int x = 0; x < size.x; ++x, ++id) {		
	      Splatting::splatPoint(a_h.x[id], a_h.y[id], a_h.z[id], a_i[id] , a_o, a_w, size.x, size.y, size.z);		
            }
        }
    }    
    for (size_t e = 0; e < nVox; ++e) {
      a_o[e] = a_w[e] > 0 ? a_o[e] / a_w[e] : 0;
    } 
    
} 

void CImageFieldOpers::
Splat(float* a_o, const Field3D& a_h, const float* a_i, StreamT stream){
    
    Vec3Di size = a_h.size();
    size_t nVox = size.prod();

    size_t id = 0;

    CMemOpers<float>::SetMem(a_o, 0.f, nVox, stream, false);
        
    for (int z = 0; z < size.z; ++z) {
        for (int y = 0; y < size.y; ++y) {
            for (int x = 0; x < size.x; ++x, ++id) {		
	      Splatting::splatPoint(a_h.x[id], a_h.y[id], a_h.z[id], a_i[id] , a_o, size.x, size.y, size.z);
            }
        }
    } 
      
}
   
// Implementation matching GPU kernel, good for testing
// void CImageFieldOpers
// ::Splat(float* a_o, const Field3D& a_h, const float* a_i,
//                bool normalize, StreamT stream)
// {
//     Vec3Di size = a_h.size();
//     int w = size.x;
//     int h = size.y;
//     int l = size.z;
//     size_t nVox = size.prod();

//     size_t id = 0;
//     size_t planeSize = size.x * size.y;

//     std::vector<float> a_w(nVox);

//     CMemOpers<float>::SetMem(a_o, 0.f, nVox, stream, false);
//     std::fill(a_w.begin(), a_w.end(), 0.f);

//     for (size_t z = 0; z < size.z; ++z) {
//         for (size_t y = 0; y < size.y; ++y) {
//             for (size_t x = 0; x < size.x; ++x, ++id) {

// 	       PYCA_ASSERT(id < nVox);
// 	       float mass = a_i[id];
// 	       float x = a_h.x[id];
// 	       float y = a_h.y[id];
// 	       float z = a_h.z[id];

// 	       int xInt = int(x);
// 	       int yInt = int(y);
// 	       int zInt = int(z);

// 	       if (x < 0 && x != xInt) --xInt;
// 	       if (y < 0 && y != yInt) --yInt;
// 	       if (z < 0 && z != zInt) --zInt;

// 	       float dx = 1.f - (x - xInt);
// 	       float dy = 1.f - (y - yInt);
// 	       float dz = 1.f - (z - zInt);

// 	       int nid = (zInt * h + yInt) * w + xInt;
// 	       float dist;
	       
// 	       if (isInside3D(xInt, yInt, zInt, w, h, l)){
// 		  dist = mass * dx * dy * dz;
// 		  PYCA_ASSERT(nid < nVox);
// 		  a_o[nid] += dist;
// 	       }
            
// 	       if (isInside3D(xInt + 1, yInt, zInt, w, h, l)){
// 		  dist = mass * (1.f-dx) * dy * dz;
// 		  PYCA_ASSERT(nid + 1 < nVox);
// 		  a_o[nid + 1] += dist;
// 	       }

// 	       if (isInside3D(xInt, yInt+1, zInt, w, h, l)){
// 		  dist = mass * dx * (1.f - dy) * dz;
// 		  PYCA_ASSERT(nid + w < nVox);
// 		  a_o[nid + w] += dist;
// 	       }
    
// 	       if (isInside3D(xInt+1, yInt+1, zInt, w, h, l)){
// 		  dist = mass * (1.f -dx) * (1.f - dy) * dz;
// 		  PYCA_ASSERT(nid + w + 1 < nVox);
// 		  a_o[nid + w + 1] += dist;
// 	       } 
    
// 	       nid += w*h;

// 	       if (isInside3D(xInt, yInt, zInt + 1, w, h, l)){
// 		  dist = mass * dx * dy * (1.f - dz);
// 		  PYCA_ASSERT(nid < nVox);
// 		  a_o[nid] +=dist;
// 	       }
            
// 	       if (isInside3D(xInt + 1, yInt, zInt+1, w, h, l)){
// 		  dist = mass * (1.f-dx) * dy * (1.f -dz);
// 		  PYCA_ASSERT(nid + 1 < nVox);
// 		  a_o[nid + 1] += dist;
// 	       }
    
// 	       if (isInside3D(xInt, yInt+1, zInt+1, w, h, l)){
// 		  dist = mass * dx * (1.f - dy) * (1.f -dz);
// 		  PYCA_ASSERT(nid + w < nVox);
// 		  a_o[nid + w] += dist;
// 	       }
    
// 	       if (isInside3D(xInt+1, yInt+1, zInt+1, w, h, l)){
// 		  dist = mass * (1.f -dx) * (1.f - dy) * (1.f -dz);
// 		  PYCA_ASSERT(nid + w + 1 < nVox);
// 		  a_o[nid + w + 1] += dist;
// 	       }
// 	    }
// 	}
//     }
// }

inline
void 
c_gradient(float* a_ox, float *a_oy, float* a_oz,
	   const float* a_i,
	   int szX, int szY, int szZ,
	   float spX, float spY, float spZ, 
	   DiffT diffType,
	   enum BoundaryCondT bc, 
	   StreamT stream)
{
   bool slice = (szZ==1);

   CFiniteDiff::
      FiniteDiff(a_ox, a_i, szX, szY, szZ, spX, spY, spZ,
		 DIM_X, diffType, bc, ACCUM_FALSE, OP_VAL, stream);
   
   CFiniteDiff::
      FiniteDiff(a_oy, a_i, szX, szY, szZ, spX, spY, spZ,
		 DIM_Y, diffType, bc, ACCUM_FALSE, OP_VAL, stream);
   
   if(slice){
      unsigned int nEl = szX*szY*szZ;
      for(unsigned int idx=0;idx<nEl;++idx){
	 a_oz[idx] = 0;
      }
   }else{
      CFiniteDiff::
	 FiniteDiff(a_oz, a_i, szX, szY, szZ, spX, spY, spZ,
		    DIM_Z, diffType, bc, ACCUM_FALSE, OP_VAL, stream);
   }
}

template <enum DiffT diffType, BoundaryCondT bc>
inline
void 
c_gradient2(float* a_ox, float *a_oy, float* a_oz,
	    const float* a_i,
	    int szX, int szY, int szZ,
	    float spX, float spY, float spZ, 
	    StreamT stream)
{
   int id=0;
   for(int z=0;z<szZ;++z){
      for(int y=0;y<szY;++y){
	 for(int x=0;x<szX;++x,++id){
	    gradientPoint<float, diffType, bc>
	       (a_i, 
		x,y,z,
		a_ox[id], a_oy[id], a_oz[id],
		szX, szY, szZ,
		1.f/spX, 1.f/spY, 1.f/spZ);
	 }
      }
   }
}

template <enum DiffT diffType>
inline
void 
c_gradient2(float* a_ox, float *a_oy, float* a_oz,
	    const float* a_i,
	    int szX, int szY, int szZ,
	    float spX, float spY, float spZ, 
	    enum BoundaryCondT bc, 
	    StreamT stream)
{
   if(bc == BC_APPROX){
      PyCA::c_gradient2<diffType, BC_APPROX>
	 (a_ox, a_oy, a_oz, a_i, szX, szY, szZ, spX, spY, spZ, stream);
   }else if(bc == BC_WRAP){
      PyCA::c_gradient2<diffType, BC_WRAP>
	 (a_ox, a_oy, a_oz, a_i, szX, szY, szZ, spX, spY, spZ, stream);
   }else if(bc == BC_CLAMP){
      PyCA::c_gradient2<diffType, BC_CLAMP>
	 (a_ox, a_oy, a_oz, a_i, szX, szY, szZ, spX, spY, spZ, stream);
   }else{
      throw PyCAException(__FILE__, __LINE__, "unknown boundary condition");
   }
}

inline
void 
c_gradient2(float* a_ox, float *a_oy, float* a_oz,
	    const float* a_i,
	    int szX, int szY, int szZ,
	    float spX, float spY, float spZ, 
	    DiffT diffType,
	    BoundaryCondT bc, 
	    StreamT stream)
{
   if(diffType == DIFF_FORWARD){
      PyCA::c_gradient2<DIFF_FORWARD>
	 (a_ox, a_oy, a_oz, a_i, szX, szY, szZ, spX, spY, spZ, bc, stream);
   }else if(diffType == DIFF_BACKWARD){
      PyCA::c_gradient2<DIFF_BACKWARD>
	 (a_ox, a_oy, a_oz, a_i, szX, szY, szZ, spX, spY, spZ, bc, stream);
      
   }else if(diffType == DIFF_CENTRAL){
      PyCA::c_gradient2<DIFF_CENTRAL>
	 (a_ox, a_oy, a_oz, a_i, szX, szY, szZ, spX, spY, spZ, bc, stream);
   }else{
      throw PyCAException(__FILE__, __LINE__, "unknown DiffT");
   }
}

void 
CImageFieldOpers::
FiniteDiff(Image3D& a_o, const Image3D& a_i, 
	   DimT dim,  DiffT diffType, 
	   enum BoundaryCondT bc, 
	   bool accum, OpT op,
	   StreamT stream)
{
   CFiniteDiff::
      FiniteDiff(a_o, a_i, 
		 dim, diffType, bc, accum, op, 
		 stream);
}

void 
CImageFieldOpers::
UpwindDiff(Image3D& h_o, const Image3D& h_i,
	   const Image3D& h_speed, 
	   DimT dim, StreamT stream)
{
   CFiniteDiff::
      UpwindDiff(h_o, h_i, h_speed, dim, stream);
}

void 
CImageFieldOpers::
UpwindGradMag(Image3D& h_o, const Image3D& h_i,
	      const Image3D& h_speed, StreamT stream)
{
   CFiniteDiff::
      UpwindGradMag(h_o, h_i, 
		    h_speed, stream);
}

void CImageFieldOpers
::Gradient(Field3D& a_o, const float* a_i, DiffT diffType, 
	   BoundaryCondT bc, StreamT stream)
{
   Vec3Di size = a_o.size();
   Vec3Df sp   = a_o.spacing();
   PyCA::c_gradient
      (a_o.x, a_o.y, a_o.z, a_i,
       size.x, size.y, size.z,
       sp.x, sp.y, sp.z, 
       diffType, bc,
       stream);
}


void CImageFieldOpers::
Gradient(Field3D& a_o, const Image3D& a_i, DiffT diffType, 
	 BoundaryCondT bc, StreamT stream)
{
    MK_CHECK2_SIZE(a_o, a_i);
    Vec3Di size = a_o.size();
    Vec3Df sp   = a_o.spacing();
    PyCA::c_gradient
       (a_o.x, a_o.y, a_o.z, a_i.get(),
	size.x, size.y, size.z,
	sp.x, sp.y, sp.z, 
	diffType, bc,
	stream);
}

void CImageFieldOpers::
Gradient2(Field3D& a_o, const Image3D& a_i, DiffT diffType, 
	  BoundaryCondT bc, StreamT stream)
{
    MK_CHECK2_SIZE(a_o, a_i);
    Vec3Di size = a_o.size();
    Vec3Df sp   = a_o.spacing();
    PyCA::c_gradient2
       (a_o.x, a_o.y, a_o.z, a_i.get(),
	size.x, size.y, size.z,
	sp.x, sp.y, sp.z, 
	diffType, bc,
	stream);
}


void CImageFieldOpers::
GradFor(Field3D& a_o, const Image3D& a_i, BoundaryCondT bc, StreamT stream){
    MK_CHECK2_SIZE(a_o, a_i);
    Vec3Di size = a_o.size();
    Vec3Df sp   = a_o.spacing();
    PyCA::c_gradient
       (a_o.x, a_o.y, a_o.z, a_i.get(),
	size.x, size.y, size.z,
	sp.x, sp.y, sp.z, 
	DIFF_FORWARD, bc,
	stream);
}

template <enum DiffT diffType, BoundaryCondT bc>
inline
void 
c_gradientMask(float* a_ox, float *a_oy, float* a_oz,
	       const float* a_i, const float* a_mask,
	       int szX, int szY, int szZ,
	       float spX, float spY, float spZ, 
	       StreamT stream)
{
   int id=0;
   for(int z=0;z<szZ;++z){
      for(int y=0;y<szY;++y){
	 for(int x=0;x<szX;++x,++id){
	    gradientPointMask<float, diffType, bc>
	       (a_i, a_mask,
		x,y,z,
		a_ox[id], a_oy[id], a_oz[id],
		szX, szY, szZ,
		1.f/spX, 1.f/spY, 1.f/spZ);
	 }
      }
   }
}

template <enum DiffT diffType>
inline
void 
c_gradientMask(float* a_ox, float *a_oy, float* a_oz,
	       const float* a_i, const float* a_mask,
	       int szX, int szY, int szZ,
	       float spX, float spY, float spZ, 
	       enum BoundaryCondT bc, 
	       StreamT stream)
{
   if(bc == BC_APPROX){
      PyCA::c_gradientMask<diffType, BC_APPROX>
	 (a_ox, a_oy, a_oz, a_i, a_mask, szX, szY, szZ, spX, spY, spZ, stream);
   }else if(bc == BC_WRAP){
      throw PyCAException(__FILE__, __LINE__, "BC_WRAP boundary condition unimplemented for masked gradient");
   }else if(bc == BC_CLAMP){
      PyCA::c_gradientMask<diffType, BC_CLAMP>
	 (a_ox, a_oy, a_oz, a_i, a_mask, szX, szY, szZ, spX, spY, spZ, stream);
   }else{
      throw PyCAException(__FILE__, __LINE__, "unknown boundary condition");
   }
}

inline
void 
c_gradientMask(float* a_ox, float *a_oy, float* a_oz,
	       const float* a_i, const float* a_mask,
	       int szX, int szY, int szZ,
	       float spX, float spY, float spZ, 
	       DiffT diffType,
	       BoundaryCondT bc, 
	       StreamT stream)
{
   if(diffType == DIFF_FORWARD){
      PyCA::c_gradientMask<DIFF_FORWARD>
	 (a_ox, a_oy, a_oz, a_i, a_mask, szX, szY, szZ, spX, spY, spZ, bc, stream);
   }else if(diffType == DIFF_BACKWARD){
      PyCA::c_gradientMask<DIFF_BACKWARD>
	 (a_ox, a_oy, a_oz, a_i, a_mask, szX, szY, szZ, spX, spY, spZ, bc, stream);
      
   }else if(diffType == DIFF_CENTRAL){
      PyCA::c_gradientMask<DIFF_CENTRAL>
	 (a_ox, a_oy, a_oz, a_i, a_mask, szX, szY, szZ, spX, spY, spZ, bc, stream);
   }else{
      throw PyCAException(__FILE__, __LINE__, "unknown DiffT");
   }
}

void CImageFieldOpers::
GradientMask(Field3D& a_o, 
	     const Image3D& a_i, const Image3D& a_mask, 
	     DiffT diffType, BoundaryCondT bc, 
	     StreamT stream)
{
    MK_CHECK2_SIZE(a_o, a_i);
    Vec3Di size = a_o.size();
    Vec3Df sp   = a_o.spacing();
    PyCA::c_gradientMask
       (a_o.x, a_o.y, a_o.z, 
	a_i.get(), a_mask.get(),
	size.x, size.y, size.z,
	sp.x, sp.y, sp.z, 
	diffType, bc,
	stream);
}

void c_gradientMag(float* a_o,
		   const float* a_i,
		   int szX, int szY, int szZ,
		   float spX, float spY, float spZ, 
		   DiffT diffType, BoundaryCondT bc, 
		   StreamT stream)
{
   bool slice = (szZ == 1);

   CFiniteDiff::
      FiniteDiff(a_o, a_i, szX, szY, szZ, spX, spY, spZ,
		 DIM_X, diffType, bc, ACCUM_FALSE, OP_SQR, stream);
   
   CFiniteDiff::
      FiniteDiff(a_o, a_i, szX, szY, szZ, spX, spY, spZ,
		 DIM_Y, diffType, bc, ACCUM_TRUE, OP_SQR, stream);
   
   if(!slice){
      CFiniteDiff::
	 FiniteDiff(a_o, a_i, szX, szY, szZ, spX, spY, spZ,
		    DIM_Z, diffType, bc, ACCUM_TRUE, OP_SQR, stream);
   }

   unsigned int nEl = szX*szY*szZ;
   for(unsigned int i=0;i<nEl;++i){
      a_o[i] = sqrt(a_o[i]);
   }
}

void CImageFieldOpers::
GradForMag(Image3D& a_o, const Image3D& a_i, StreamT stream){
    MK_CHECK2_SIZE(a_o, a_i);

    Vec3Di size = a_o.size();
    Vec3Df sp   = a_o.spacing();
    
    PyCA::c_gradientMag
       (a_o.get(), a_i.get(),
	size.x,size.y,size.z,
	sp.x,sp.y,sp.z, 
	DIFF_FORWARD,BC_APPROX,
	stream);
}

void CImageFieldOpers::
GradientMag(Image3D& a_o, const Image3D& a_i, 
	    DiffT diffType, BoundaryCondT bc, 
	    StreamT stream)
{
    MK_CHECK2_SIZE(a_o, a_i);

    Vec3Di size = a_o.size();
    Vec3Df sp   = a_o.spacing();
    
    PyCA::c_gradientMag
       (a_o.get(), a_i.get(),
	size.x,size.y,size.z,
	sp.x,sp.y,sp.z, 
	diffType, bc,
	stream);
}

inline
void 
c_divergence(float* a_o,
	     const float* a_ix, const float* a_iy, const float* a_iz,
	     int   szX, int szY, int szZ,
	     float spX, float spY, float spZ,
	     DiffT diffType, BoundaryCondT bc, 
	     StreamT stream)
{
   bool slice = (szZ == 1);
   CFiniteDiff::
      FiniteDiff(a_o, a_ix, szX, szY, szZ, spX, spY, spZ,
		 DIM_X, diffType, bc, ACCUM_FALSE, OP_VAL, stream);
   
   CFiniteDiff::
      FiniteDiff(a_o, a_iy, szX, szY, szZ, spX, spY, spZ,
		 DIM_Y, diffType, bc, ACCUM_TRUE, OP_VAL, stream);
   
   if(!slice){
      CFiniteDiff::
	 FiniteDiff(a_o, a_iz, szX, szY, szZ, spX, spY, spZ,
		    DIM_Z, diffType, bc, ACCUM_TRUE, OP_VAL, stream);
   }

}

void CImageFieldOpers::
Divergence(Image3D& a_o, const Field3D& a_i, DiffT diffType, 
	   BoundaryCondT bc, StreamT stream)
{
    MK_CHECK2_SIZE(a_o, a_i);
    
    Vec3Di size = a_o.size();
    Vec3Df sp   = a_o.spacing();

    PyCA::c_divergence
       (a_o.get(),
	a_i.x, a_i.y, a_i.z,
	size.x, size.y, size.z,
	sp.x, sp.y, sp.z, 
	diffType, bc, 
	stream);
}

void CImageFieldOpers::
DivBack(Image3D& a_o, const Field3D& a_i, BoundaryCondT bc, StreamT stream){
   MK_CHECK2_SIZE(a_o, a_i);
   
   Vec3Di size = a_o.size();
   Vec3Df sp   = a_o.spacing();
   
   PyCA::c_divergence
      (a_o.get(),
       a_i.x, a_i.y, a_i.z,
       size.x, size.y, size.z,
       sp.x, sp.y, sp.z, 
       DIFF_BACKWARD, bc, 
       stream);
}

/**
 * Compute the magnitude image
 * a_o[i] = sqrt(a_i[i].x^2 + a_i[i].y^2 + a_i[i].z^2) 
 */

void CImageFieldOpers::
Magnitude(Image3D& a_o, const Field3D& a_i, StreamT stream){
    MK_CHECK2_SIZE(a_o, a_i);
    size_t n = a_o.nVox();
    for (size_t id=0; id< n; ++id) {
        a_o[id] = sqrt(a_i.x[id] * a_i.x[id] + a_i.y[id] * a_i.y[id] + a_i.z[id] * a_i.z[id]);
    }
}

/**
 * Compute the magnitude array
 * a_o[i] = a_i[i].x^2 + a_i[i].y^2 + a_i[i].z^2 
 */

void CImageFieldOpers::
SqrMagnitude(Image3D& a_o, const Field3D& a_i, StreamT stream){
    MK_CHECK2_SIZE(a_o, a_i);
    size_t n = a_o.nVox();
    for (size_t id=0; id< n; ++id) {
        a_o[id] = a_i.x[id] * a_i.x[id] + a_i.y[id] * a_i.y[id] + a_i.z[id] * a_i.z[id];
    }
}

 /**
 * Compute dot product image
 * a_o[i] = a_i[i].x * a_i1[i].x + a_i[i].y * a_i1[i].y + a_i[i].z * a_i1[i].z
 */

void CImageFieldOpers::
ComponentDotProd(Image3D& a_o, const Field3D& a_i, const Field3D& a_i1, StreamT stream) {
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t n = a_o.nVox();
    for (size_t id=0; id< n; ++id) {
        a_o[id] = a_i.x[id] * a_i1.x[id] + a_i.y[id] * a_i1.y[id] + a_i.z[id] * a_i1.z[id];
    }
}

/** @brief a_o.x = a_i.x + a_i1, a_o.y = a_i.y + a_i1, a_o.z = a_i.z + a_i1 */

void CImageFieldOpers::
Add(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t n = a_o.nVox();

    for (size_t id=0; id< n; ++id) {
        float v = a_i1[id];
        a_o.x[id] = a_i.x[id] + v;
        a_o.y[id] = a_i.y[id] + v;
        a_o.z[id] = a_i.z[id] + v;
    }
}

/** @brief a_o.x = a_i.x - a_i1, a_o.y = a_i.y - a_i1, a_o.z = a_i.z - a_i1 */

void CImageFieldOpers::
Sub(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t n = a_o.nVox();
    for (size_t id=0; id< n; ++id) {
        float v = a_i1[id];
        a_o.x[id] = a_i.x[id] - v;
        a_o.y[id] = a_i.y[id] - v;
        a_o.z[id] = a_i.z[id] - v;
    }

}
/** @brief a_o.x = a_i.x * a_i1, a_o.y = a_i.y * a_i1, a_o.z = a_i.z * a_i1 */

void CImageFieldOpers::
Mul(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t n = a_o.nVox();
    for (size_t id=0; id< n; ++id) {
        float v = a_i1[id];
        a_o.x[id] = a_i.x[id] * v;
        a_o.y[id] = a_i.y[id] * v;
        a_o.z[id] = a_i.z[id] * v;
    }
}

/** @brief a_o.x = a_i.x / a_i1, a_o.y = a_i.y / a_i1, a_o.z = a_i.z / a_i1 */

void CImageFieldOpers::
Div(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t n = a_o.nVox();
    for (size_t id=0; id< n; ++id) {
        float v = 1.f / a_i1[id];
        a_o.x[id] = a_i.x[id] * v;
        a_o.y[id] = a_i.y[id] * v;
        a_o.z[id] = a_i.z[id] * v;
    }
}

/** @brief a_o.x += a_i, a_o.y += a_i, a_o.y += a_i, */

void CImageFieldOpers::
Add_I(Field3D& a_o, const Image3D& a_i, StreamT stream){
    MK_CHECK2_SIZE(a_o, a_i);
    size_t n = a_o.nVox();
    for (size_t id=0; id< n; ++id) {
        float v = a_i[id];
        a_o.x[id] += v;
        a_o.y[id] += v;
        a_o.z[id] += v;
    }    
}
/** @brief a_o.x -= a_i, a_o.y -= a_i, a_o.y -= a_i, */

void CImageFieldOpers::
Sub_I(Field3D& a_o, const Image3D& a_i, StreamT stream){
    MK_CHECK2_SIZE(a_o, a_i);
    size_t n = a_o.nVox();
    for (size_t id=0; id< n; ++id) {
        float v = a_i[id];
        a_o.x[id] -= v;
        a_o.y[id] -= v;
        a_o.z[id] -= v;
    }    
}
/** @brief a_o.x *= a_i, a_o.y *= a_i, a_o.y *= a_i, */

void CImageFieldOpers::
Mul_I(Field3D& a_o, const Image3D& a_i, StreamT stream){
    MK_CHECK2_SIZE(a_o, a_i);
    size_t n = a_o.nVox();
    for (size_t id=0; id< n; ++id) {
        float v = a_i[id];
        a_o.x[id] *= v;
        a_o.y[id] *= v;
        a_o.z[id] *= v;
    }    
}

/** @brief a_o.x /= a_i, a_o.y /= a_i, a_o.y /= a_i*/

void CImageFieldOpers::
Div_I(Field3D& a_o, const Image3D& a_i, StreamT stream){
    MK_CHECK2_SIZE(a_o, a_i);
    size_t n = a_o.nVox();
    for (size_t id=0; id< n; ++id) {
        float v = 1.f / a_i[id];
        a_o.x[id] *= v;
        a_o.y[id] *= v;
        a_o.z[id] *= v;
    }    
}

/** @brief a_o = a_i + a_i1 * a_i2 */

void CImageFieldOpers::
Add_Mul(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const Image3D& a_i2, StreamT stream){
    MK_CHECK4_SIZE(a_o, a_i, a_i1, a_i2);
    size_t n = a_o.nVox();
    for (size_t id=0; id< n; ++id) {
        float v = a_i2[id];
        a_o.x[id] = a_i.x[id] + a_i1.x[id] * v;
        a_o.y[id] = a_i.y[id] + a_i1.y[id] * v;
        a_o.z[id] = a_i.z[id] + a_i1.z[id] * v;
    }
}
/** @brief a_o = a_o + a_i * a_i1 */

void CImageFieldOpers::
Add_Mul_I(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t n = a_o.nVox();
    for (size_t id=0; id< n; ++id) {
        float v = a_i1[id];
        a_o.x[id] += a_i.x[id] * v;
        a_o.y[id] += a_i.y[id] * v;
        a_o.z[id] += a_i.z[id] * v;
    }
}

/** @brief a_o = a_i - a_i1 * a_i2 */

void CImageFieldOpers::
Sub_Mul(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const Image3D& a_i2, StreamT stream){
    MK_CHECK4_SIZE(a_o, a_i, a_i1, a_i2);
    size_t n = a_o.nVox();
    for (size_t id=0; id< n; ++id) {
        float v = a_i2[id];
        a_o.x[id] = a_i.x[id] - a_i1.x[id] * v;
        a_o.y[id] = a_i.y[id] - a_i1.y[id] * v;
        a_o.z[id] = a_i.z[id] - a_i1.z[id] * v;
    }
}
/** @brief a_o = a_o - a_i * a_i1 */

void CImageFieldOpers::
Sub_Mul_I(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t n = a_o.nVox();
    for (size_t id=0; id< n; ++id) {
        float v = a_i1[id];
        a_o.x[id] -= a_i.x[id] * v;
        a_o.y[id] -= a_i.y[id] * v;
        a_o.z[id] -= a_i.z[id] * v;
    }
}

/** @brief a_o = a_i * a_i1 * c (a_o.x = a_i.x * a_i1 * c)*/

void CImageFieldOpers::
MulMulC(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, const float& c, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t n = a_o.nVox();

    for (size_t id=0; id< n; ++id) {
        float v = a_i1[id] * c;
        a_o.x[id] = a_i.x[id] * v;
        a_o.y[id] = a_i.y[id] * v;
        a_o.z[id] = a_i.z[id] * v;
    }
}
/** @brief a_o = a_o * a_i * c  (a_o.x = a_o.x * a_i * c)*/

void CImageFieldOpers::
MulMulC_I(Field3D& a_o, const Image3D& a_i, const float& c, StreamT stream, bool onDev){
    MK_CHECK2_SIZE(a_o, a_i);
    size_t n = a_o.nVox();

    for (size_t id=0; id< n; ++id) {
        float v = a_i[id] * c;
        a_o.x[id] *= v;
        a_o.y[id] *= v;
        a_o.z[id] *= v;
    }
}

/** @brief a_o = a_i + a_i1 * a_i2 * c */

void CImageFieldOpers::
Add_MulMulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const Image3D& a_i2, const float& c,  StreamT stream, bool onDev){
    MK_CHECK4_SIZE(a_o, a_i, a_i1, a_i2);
    size_t n = a_o.nVox();
    for (size_t id=0; id< n; ++id) {
        float v = a_i2[id] * c;
        a_o.x[id] = a_i.x[id] + a_i1.x[id] * v;
        a_o.y[id] = a_i.y[id] + a_i1.y[id] * v;
        a_o.z[id] = a_i.z[id] + a_i1.z[id] * v;
    }    
}

/** @brief a_o = a_o + a_i * a_i1 * c */

void CImageFieldOpers::
Add_MulMulC_I(const Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, const float& c,  StreamT stream, bool onDev){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t n = a_o.nVox();
    for (size_t id=0; id< n; ++id) {
        float v = a_i1[id] * c;
        a_o.x[id] += a_i.x[id] * v;
        a_o.y[id] += a_i.y[id] * v;
        a_o.z[id] += a_i.z[id] * v;
    }
}

void CImageFieldOpers
::JacDetH(Image3D& a_detJ,
          const Field3D& a_Xg, const Field3D& a_Yg, const Field3D& a_Zg, StreamT stream)
{
    MK_CHECK2_SIZE(a_detJ, a_Xg);
    size_t n = a_detJ.nVox();
    bool slice = (a_Xg.size().z == 1);
    
    if(slice){
       // 2x2 matrix determinant
       for (size_t id=0; id< n; ++id) {
	  a_detJ[id] = a_Xg.x[id]*a_Yg.y[id]-a_Yg.x[id]*a_Xg.y[id];
       }
    }else{
       // 3x3 matrix determinant
       for (size_t id=0; id< n; ++id) {
	  float a00 = a_Xg.x[id], a01 = a_Xg.y[id], a02 = a_Xg.z[id];
	  float a10 = a_Yg.x[id], a11 = a_Yg.y[id], a12 = a_Yg.z[id];
	  float a20 = a_Zg.x[id], a21 = a_Zg.y[id], a22 = a_Zg.z[id];
	  
	  a_detJ[id] = det(a00, a01, a02,
			   a10, a11, a12,
			   a20, a21, a22);
	  
       }
    }
}

void CImageFieldOpers
::JacDetV(Image3D& a_detJ,
          const Field3D& a_Xg, const Field3D& a_Yg, const Field3D& a_Zg, StreamT stream)
{
    MK_CHECK2_SIZE(a_detJ, a_Xg);
    size_t n = a_detJ.nVox();
    for (size_t id=0; id< n; ++id) {
        float a00 = 1.f + a_Xg.x[id], a01 = a_Xg.y[id], a02 = a_Xg.z[id];
        float a10 = a_Yg.x[id], a11 = 1.f + a_Yg.y[id], a12 = a_Yg.z[id];
        float a20 = a_Zg.x[id], a21 = a_Zg.y[id], a22 = 1 + a_Zg.z[id];
        
        a_detJ[id] = det(a00, a01, a02,
                         a10, a11, a12,
                         a20, a21, a22);

    }
}

void CImageFieldOpers
::JacDetVInv(Image3D& a_detJ,
             const Field3D& a_Xg, const Field3D& a_Yg, const Field3D& a_Zg, StreamT stream)
{
    MK_CHECK2_SIZE(a_detJ, a_Xg);
    size_t n = a_detJ.nVox();
    for (size_t id=0; id< n; ++id) {
        float a00 = 1.f - a_Xg.x[id], a01 = a_Xg.y[id], a02 = a_Xg.z[id];
        float a10 = a_Yg.x[id], a11 = 1.f - a_Yg.y[id], a12 = a_Yg.z[id];
        float a20 = a_Zg.x[id], a21 = a_Zg.y[id], a22 = 1.f - a_Zg.z[id];
        
        a_detJ[id] = det(a00, a01, a02,
                         a10, a11, a12,
                         a20, a21, a22);

    }
}


template<enum DiffT diffType, BoundaryCondT bc, bool slice>
inline
static
void
c_JacDetH(Image3D& a_jdet, const Field3D& a_h)
{
    MK_CHECK2_SIZE(a_jdet, a_h);
    Vec3Di sz = a_h.size();
    Vec3Df sp = a_h.spacing();
    Vec3Df isp(1.0/sp.x, 1.0/sp.y, 1.0/sp.z);
    int id = 0;
    for (int k=0; k < sz.z; ++k){
       for (int j=0; j < sz.y; ++j){
            for (int i=0; i < sz.x; ++i, ++id){
	       jacDetPoint<float,diffType,bc,slice>
		  (a_jdet[id],
		   a_h.x,a_h.y,a_h.z,
		   i,j,k,
		   sz.x,sz.y,sz.z,
		   isp.x,isp.y,isp.z);
            }
       }
    }
}

template<enum DiffT diffType>
inline
static
void
c_JacDetH(Image3D& a_jdet, const Field3D& a_h, BoundaryCondT bc)
{
   bool slice = (a_jdet.size().z == 1);
   if(bc == BC_APPROX){
      if(slice){
	 PyCA::c_JacDetH<diffType, BC_APPROX, SLICE_TRUE>(a_jdet, a_h);
      }else{
	 PyCA::c_JacDetH<diffType, BC_APPROX, SLICE_FALSE>(a_jdet, a_h);
      }
   }else if(bc == BC_WRAP){
      if(slice){
	 PyCA::c_JacDetH<diffType, BC_WRAP, SLICE_TRUE>(a_jdet, a_h);
      }else{
	 PyCA::c_JacDetH<diffType, BC_WRAP, SLICE_FALSE>(a_jdet, a_h);
      }
   }else if(bc == BC_CLAMP){
      if(slice){
	 PyCA::c_JacDetH<diffType, BC_CLAMP, SLICE_TRUE>(a_jdet, a_h);
      }else{
	 PyCA::c_JacDetH<diffType, BC_CLAMP, SLICE_FALSE>(a_jdet, a_h);
      }
   }else{
      throw PyCAException(__FILE__, __LINE__, "unknown BoundaryCondT");
   }
}

void CImageFieldOpers
::JacDetH(Image3D& a_jdet, const Field3D& a_h, 
	  DiffT diffType, BoundaryCondT bc, 
	  StreamT stream)
{
   if(diffType == DIFF_FORWARD){
      PyCA::c_JacDetH<DIFF_FORWARD>(a_jdet, a_h, bc);
   }else if(diffType == DIFF_BACKWARD){
      PyCA::c_JacDetH<DIFF_BACKWARD>(a_jdet, a_h, bc);
   }else if(diffType == DIFF_CENTRAL){
      PyCA::c_JacDetH<DIFF_CENTRAL>(a_jdet, a_h, bc);
   }else{
      throw PyCAException(__FILE__, __LINE__, "unknown DiffT");
   }
}

#include "CImageFieldOpers_inst.cxx"
    
} // end namespace PyCA

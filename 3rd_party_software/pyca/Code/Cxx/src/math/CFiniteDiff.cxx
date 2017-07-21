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


#include "CFiniteDiff.h"
#include <Image3D.h>
#include <CMemOpers.h>

namespace PyCA {

// template params for slice,diffType allow compiler to optimize
// out any conditional statements based on these vars 
// jsp 2012

template< class T, bool accum, OpT op, int stride, int off >
inline
void setVal(T* array, 
	    int w, int h, int d,
	    int x, int y, int z,
	    const T& val)
{
   if(accum){
      if(op == OP_SQR){
	 getSOVal<T,stride,off>(array,w,h,d,x,y,z) += val*val;
      }else{
	 getSOVal<T,stride,off>(array,w,h,d,x,y,z) += val;
      }
   }else{
      if(op == OP_SQR){
	 getSOVal<T,stride,off>(array,w,h,d,x,y,z) = val*val;
      }else{
	 getSOVal<T,stride,off>(array,w,h,d,x,y,z) = val;
      }
   }
}

template < class U, class T, 
	   DimT dim, DiffT diffType, 
	   BoundaryCondT bc,
	   bool accum, OpT op,
	   int inStride, int inOff,
	   int outStride, int outOff >
void 
c_finiteDiff(U* diff, 
	     const T* array, 
	     int w, int h, int d,
	     T spacing)
{
  int i, j, k;
  // next and prev coords for each dimension
  int xn, xp, yn, yp, zn, zp; 
  // min and max points to process
  int xinit, xfin, yinit, yfin, zinit, zfin; 
  xinit = yinit = zinit = 0;
  xfin = w-1;
  yfin = h-1;
  zfin = d-1;

  U val;
   
  // invert spacing
  spacing = 1.f/spacing;
  T edgeSpacing = spacing;
  if(diffType == DIFF_CENTRAL){
    spacing *= 0.5f;
    if(bc != BC_APPROX){
      edgeSpacing *= 0.5f;
    }
  }
   
  if(diffType == DIFF_FORWARD || diffType == DIFF_CENTRAL){
    if(dim == DIM_X) xfin--;
    if(dim == DIM_Y) yfin--;
    if(dim == DIM_Z) zfin--;
  }
  if(diffType == DIFF_BACKWARD || diffType == DIFF_CENTRAL){
    if(dim == DIM_X) xinit++;
    if(dim == DIM_Y) yinit++;
    if(dim == DIM_Z) zinit++;
  }
   
  // process interior pixels
  for (k = zinit; k <= zfin; k++){
    for (j = yinit; j <= yfin; j++){
      for (i = xinit; i <= xfin; i++){
	xp = xn = i; 
	yp = yn = j; 
	zp = zn = k;
	    
	if(diffType == DIFF_FORWARD || diffType == DIFF_CENTRAL){
	  if(dim == DIM_X){
	    xn++;
	  }else if(dim == DIM_Y){
	    yn++;
	  }else if(dim == DIM_Z){
	    zn++;
	  }
	}
	if(diffType == DIFF_BACKWARD || diffType == DIFF_CENTRAL){
	  if(dim == DIM_X){
	    xp--;
	  }else if(dim == DIM_Y){
	    yp--;
	  }else if(dim == DIM_Z){
	    zp--;
	  }
	}
	// calculate and set value
	val = spacing*
	   (getSOVal<T,inStride,inOff>(array,w,h,d,xn,yn,zn) - 
	    getSOVal<T,inStride,inOff>(array,w,h,d,xp,yp,zp));
	setVal<U,accum,op,outStride,outOff>(diff,w,h,d,i,j,k,val);
      }
    }
  }
   
  // deal with borders
  if(dim == DIM_X){
    for (k = zinit; k <= zfin; k++){
      for (j = yinit; j <= yfin; j++){
	if(diffType == DIFF_BACKWARD){
	  if (bc == BC_WRAP){
	     val = spacing*
		(getSOVal<T,inStride,inOff>(array,w,h,d,0,j,k)-
		 getSOVal<T,inStride,inOff>(array,w,h,d,w-1,j,k));
	     setVal<U,accum,op,outStride,outOff>(diff,w,h,d,0,j,k,val);
	  }else{
	     setVal<U,accum,op,outStride,outOff>(diff,w,h,d,0,j,k,(U)0);
	  }
	}else if(diffType == DIFF_CENTRAL){
	  if (bc == BC_WRAP){
	     val = spacing*
		(getSOVal<T,inStride,inOff>(array,w,h,d,1,j,k)-
		 getSOVal<T,inStride,inOff>(array,w,h,d,w-1,j,k));
	     setVal<U,accum,op,outStride,outOff>(diff,w,h,d,0,j,k,val);
	     val = spacing*
		(getSOVal<T,inStride,inOff>(array,w,h,d,0,j,k)-
		 getSOVal<T,inStride,inOff>(array,w,h,d,w-2,j,k));
	     setVal<U,accum,op,outStride,outOff>(diff,w,h,d,w-1,j,k,val);
	  }else{
	     val = edgeSpacing*
		(getSOVal<T,inStride,inOff>(array,w,h,d,1,j,k)-
		 getSOVal<T,inStride,inOff>(array,w,h,d,0,j,k));
	     setVal<U,accum,op,outStride,outOff>(diff,w,h,d,0,j,k,val);
	     val = edgeSpacing*
		(getSOVal<T,inStride,inOff>(array,w,h,d,w-1,j,k)-
		 getSOVal<T,inStride,inOff>(array,w,h,d,w-2,j,k));
	     setVal<U,accum,op,outStride,outOff>(diff,w,h,d,w-1,j,k,val);
	  }
	}else if(diffType == DIFF_FORWARD){
	   if (bc == BC_WRAP){
	      val = spacing*
		 (getSOVal<T,inStride,inOff>(array,w,h,d,0,j,k)-
		  getSOVal<T,inStride,inOff>(array,w,h,d,w-1,j,k));
	      setVal<U,accum,op,outStride,outOff>(diff,w,h,d,w-1,j,k,val);
	   }else{
	      setVal<U,accum,op,outStride,outOff>(diff,w,h,d,w-1,j,k,(U)0);
	   }
	}
      }
    }
  }else if(dim == DIM_Y){
    for (k = zinit; k <= zfin; k++){
      for (i = xinit; i <= xfin; i++){
	if(diffType == DIFF_BACKWARD){
	  if (bc == BC_WRAP){
	     val = spacing*
		(getSOVal<T,inStride,inOff>(array,w,h,d,i,0,k)-
		 getSOVal<T,inStride,inOff>(array,w,h,d,i,h-1,k));
	     setVal<U,accum,op,outStride,outOff>(diff,w,h,d,i,0,k,val);
	  }else{
	     setVal<U,accum,op,outStride,outOff>(diff,w,h,d,i,0,k,(U)0);
	  }
	}else if(diffType == DIFF_CENTRAL){
	  if (bc == BC_WRAP){
	     val = spacing*
		(getSOVal<T,inStride,inOff>(array,w,h,d,i,1,k)-
		 getSOVal<T,inStride,inOff>(array,w,h,d,i,h-1,k));
	     setVal<U,accum,op,outStride,outOff>(diff,w,h,d,i,0,k,val);
	     val = spacing*
		(getSOVal<T,inStride,inOff>(array,w,h,d,i,0,k)-
		 getSOVal<T,inStride,inOff>(array,w,h,d,i,h-2,k));
	     setVal<U,accum,op,outStride,outOff>(diff,w,h,d,i,h-1,k,val);
	  }else{
	     val = edgeSpacing*
		(getSOVal<T,inStride,inOff>(array,w,h,d,i,1,k)-
		 getSOVal<T,inStride,inOff>(array,w,h,d,i,0,k));
	     setVal<U,accum,op,outStride,outOff>(diff,w,h,d,i,0,k,val);
	     val = edgeSpacing*
		(getSOVal<T,inStride,inOff>(array,w,h,d,i,h-1,k)-
		 getSOVal<T,inStride,inOff>(array,w,h,d,i,h-2,k));
	     setVal<U,accum,op,outStride,outOff>(diff,w,h,d,i,h-1,k,val);
	  }
	}else if(diffType == DIFF_FORWARD){
	  if (bc == BC_WRAP){
	     val = spacing*
		(getSOVal<T,inStride,inOff>(array,w,h,d,i,0,k)-
		 getSOVal<T,inStride,inOff>(array,w,h,d,i,h-1,k));
	     setVal<U,accum,op,outStride,outOff>(diff,w,h,d,i,h-1,k,val);
	  }else{
	     setVal<U,accum,op,outStride,outOff>(diff,w,h,d,i,h-1,k,(U)0);
	  }
	}
      }
    }
  }else if(dim == DIM_Z){
    for (j = yinit; j <= yfin; j++){
      for (i = xinit; i <= xfin; i++){
	if(diffType == DIFF_BACKWARD){
	  if (bc == BC_WRAP){
	     val = spacing*
		(getSOVal<T,inStride,inOff>(array,w,h,d,i,j,0)-
		 getSOVal<T,inStride,inOff>(array,w,h,d,i,j,d-1));
	     setVal<U,accum,op,outStride,outOff>(diff,w,h,d,i,j,0,val);
	  }else{
	     setVal<U,accum,op,outStride,outOff>(diff,w,h,d,i,j,0,(U)0);
	  }
	}else if(diffType == DIFF_CENTRAL){
	  if (bc == BC_WRAP){
	     val = spacing*
		(getSOVal<T,inStride,inOff>(array,w,h,d,i,j,1)-
		 getSOVal<T,inStride,inOff>(array,w,h,d,i,j,d-1));
	     setVal<U,accum,op,outStride,outOff>(diff,w,h,d,i,j,0,val);
	     val = spacing*
		(getSOVal<T,inStride,inOff>(array,w,h,d,i,j,0)-
		 getSOVal<T,inStride,inOff>(array,w,h,d,i,j,d-2));
	     setVal<U,accum,op,outStride,outOff>(diff,w,h,d,i,j,d-1,val);
	  }else{
	     val = edgeSpacing*
		(getSOVal<T,inStride,inOff>(array,w,h,d,i,j,1)-
		 getSOVal<T,inStride,inOff>(array,w,h,d,i,j,0));
	     setVal<U,accum,op,outStride,outOff>(diff,w,h,d,i,j,0,val);
	     val = edgeSpacing*
		(getSOVal<T,inStride,inOff>(array,w,h,d,i,j,d-1)-
		 getSOVal<T,inStride,inOff>(array,w,h,d,i,j,d-2));
	     setVal<U,accum,op,outStride,outOff>(diff,w,h,d,i,j,d-1,val);
	  }
	}else if(diffType == DIFF_FORWARD){
	  if (bc == BC_WRAP){
	     val = spacing*
		(getSOVal<T,inStride,inOff>(array,w,h,d,i,j,0)-
		 getSOVal<T,inStride,inOff>(array,w,h,d,i,j,d-1));
	     setVal<U,accum,op,outStride,outOff>(diff,w,h,d,i,j,d-1,val);
	  }else{
	     setVal<U,accum,op,outStride,outOff>(diff,w,h,d,i,j,d-1,(U)0);
	  }
	}
      }
    }
  }
}

template <enum DimT dim, enum DiffT diffType, 
	  enum BoundaryCondT bc, 
	  bool accum, OpT op>
static
inline
void
c_finiteDiff(float* h_o, const float* h_i,
	     int szX, int szY, int szZ,
	     float spX, float spY, float spZ, 
	     StreamT stream)
{
   bool slice = (szZ == 1);

   float sp;
   if(dim == DIM_X){
      sp = spX;
   }else if(dim == DIM_Y){
      sp = spY;
   }else if(dim == DIM_Z){ 
      sp = spZ;
   }else{
      throw PyCAException(__FILE__, __LINE__, "unknown DimT");
   }

   if(slice && dim == DIM_Z){ 
      CMemOpers<float>::SetMem(h_o,0.f,szX*szY*szZ,stream,false);
   }else{
      c_finiteDiff<float,float,dim,diffType,bc,accum,op,1,0,1,0>
	 (h_o, h_i, szX, szY, szZ, sp);
   }
}

template <enum BoundaryCondT bc, 
	  bool accum, OpT op>
static
inline
void
c_finiteDiff(float* h_o, const float* h_i,
	     int szX, int szY, int szZ,
	     float spX, float spY, float spZ, 
	     DimT dim, DiffT diffType, 
	     StreamT stream)
{
   if(diffType == DIFF_FORWARD){
      if(dim == DIM_X){
	 c_finiteDiff<DIM_X, DIFF_FORWARD, bc, accum, op>
	    (h_o, h_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     stream);
      }else if(dim == DIM_Y){
	 c_finiteDiff<DIM_Y, DIFF_FORWARD, bc, accum, op>
	    (h_o, h_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     stream);
      }else if(dim == DIM_Z){
	 c_finiteDiff<DIM_Z, DIFF_FORWARD, bc, accum, op>
	    (h_o, h_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     stream);
      }else{
	 throw PyCAException(__FILE__, __LINE__, "unknown DimT");
      }
   }else if(diffType == DIFF_BACKWARD){
      if(dim == DIM_X){
	 c_finiteDiff<DIM_X, DIFF_BACKWARD, bc, accum, op>
	    (h_o, h_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     stream);
      }else if(dim == DIM_Y){
	 c_finiteDiff<DIM_Y, DIFF_BACKWARD, bc, accum, op>
	    (h_o, h_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     stream);
      }else if(dim == DIM_Z){
	 c_finiteDiff<DIM_Z, DIFF_BACKWARD, bc, accum, op>
	    (h_o, h_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     stream);
      }else{
	 throw PyCAException(__FILE__, __LINE__, "unknown DimT");
      }
   }else if(diffType == DIFF_CENTRAL){
      if(dim == DIM_X){
	 c_finiteDiff<DIM_X, DIFF_CENTRAL, bc, accum, op>
	    (h_o, h_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     stream);
      }else if(dim == DIM_Y){
	 c_finiteDiff<DIM_Y, DIFF_CENTRAL, bc, accum, op>
	    (h_o, h_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     stream);
      }else if(dim == DIM_Z){
	 c_finiteDiff<DIM_Z, DIFF_CENTRAL, bc, accum, op>
	    (h_o, h_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     stream);
      }else{
	 throw PyCAException(__FILE__, __LINE__, "unknown DimT");
      }
   }else{
      throw PyCAException(__FILE__, __LINE__, "unknown DiffT");
   }
}

template <bool accum, OpT op>
static
inline
void
c_finiteDiff(float* h_o, const float* h_i,
	     int szX, int szY, int szZ,
	     float spX, float spY, float spZ, 
	     DimT dim, DiffT diffType, 
	     enum BoundaryCondT bc, 
	     StreamT stream)
{
   if(bc == BC_APPROX){
      c_finiteDiff<BC_APPROX, accum, op>
	 (h_o, h_i,
	  szX, szY, szZ,
	  spX, spY, spZ, 
	  dim, diffType,
	  stream);
   }else if(bc == BC_WRAP){
	 c_finiteDiff<BC_WRAP, accum, op>
	    (h_o, h_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     dim, diffType,
	     stream);
   }else if(bc == BC_CLAMP){
	 c_finiteDiff<BC_CLAMP, accum, op>
	    (h_o, h_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     dim, diffType,
	     stream);
   }else{
      throw PyCAException(__FILE__, __LINE__, 
				   "unknown boundary condition (BoundaryCondT)");
   }
}

void CFiniteDiff
::FiniteDiff(float* h_o, const float* h_i,
	     int szX, int szY, int szZ,
	     float spX, float spY, float spZ, 
	     DimT dim, DiffT diffType, 
	     enum BoundaryCondT bc, 
	     bool accum, OpT op,
	     StreamT stream)
{
   if(accum){
      if(op == OP_VAL){
	 c_finiteDiff<ACCUM_TRUE, OP_VAL>
	    (h_o, h_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     dim, diffType,
	     bc,
	     stream);
      }else if(op == OP_SQR){
	 c_finiteDiff<ACCUM_TRUE, OP_SQR>
	    (h_o, h_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     dim, diffType,
	     bc,
	     stream);
      }else{
	 throw PyCAException(__FILE__, __LINE__, 
				      "unknown OpT");
      }
   }else{
      if(op == OP_VAL){
	 c_finiteDiff<ACCUM_FALSE, OP_VAL>
	    (h_o, h_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     dim, diffType,
	     bc,
	     stream);
      }else if(op == OP_SQR){
	 c_finiteDiff<ACCUM_FALSE, OP_SQR>
	    (h_o, h_i,
	     szX, szY, szZ,
	     spX, spY, spZ, 
	     dim, diffType,
	     bc,
	     stream);
      }else{
	 throw PyCAException(__FILE__, __LINE__, 
				      "unknown OpT");
      }
   }
}

void CFiniteDiff
::FiniteDiff(Image3D& h_o, const Image3D& h_i, 
	     DimT dim, DiffT diffType, 
	     enum BoundaryCondT bc, 
	     bool accum, OpT op,
	     StreamT stream)
{
   MK_CHECK2_SIZE(h_o, h_i);
   Vec3Di sz = h_o.size();
   Vec3Df sp = h_o.spacing();
   
   CFiniteDiff::FiniteDiff(h_o.get(), h_i.get(), 
			   sz.x, sz.y, sz.z, 
			   sp.x, sp.y, sp.z,
			   dim, diffType,
			   bc, accum, op, stream);
   
}

#define BG_CLAMP BACKGROUND_STRATEGY_CLAMP

template < DimT dim, bool slice>
void 
c_upwindDiff(float* rtn, 
	     const float* array, 
	     const float* speed, 
	     int szX, int szY, int szZ,
	     float spX, float spY, float spZ)
{

   float sp;
   if(dim == DIM_X){
      sp = spX;
   }else if(dim == DIM_Y){
      sp = spY;
   }else if(dim == DIM_Z){
      sp = spZ;
   }else{
      throw PyCAException(__FILE__, __LINE__, "unknown DimT");
   }

   for (int z=0; z<szZ; z++){
      for (int y=0; y<szY; y++){
	 for (int x=0; x<szX; x++){

	    float v = getVal<float>(array, szX, szY, szZ, x, y, z); // val
	    // previous and next values
	    float vp, vn;
	    if(dim == DIM_X){
	       vp = getSafeVal<float,BG_CLAMP>
		  (array, szX, szY, szZ, x-1, y, z);
	       vn = getSafeVal<float,BG_CLAMP>
		  (array, szX, szY, szZ, x+1, y, z);
	    }else if(dim == DIM_Y){
	       vp = getSafeVal<float,BG_CLAMP>
		  (array, szX, szY, szZ, x, y-1, z);
	       vn = getSafeVal<float,BG_CLAMP>
		  (array, szX, szY, szZ, x, y+1, z);
	    }else if(dim == DIM_Z){
	       if(slice){
		  vp = v;
		  vn = v;
	       }else{
		  vp = getSafeVal<float,BG_CLAMP>
		     (array, szX, szY, szZ, x, y, z-1);
		  vn = getSafeVal<float,BG_CLAMP>
		     (array, szX, szY, szZ, x, y, z+1);
	       }
	    }
	    
	    
	    float spd = getVal<float>(speed, szX, szY, szZ, x, y, z); // speed
	    
	    float dx = 0.f;
	    if (spd < 0.0f){
	       // forward difference
	       dx = (vn - v)/sp;
	    }else{
	       // backward difference
	       dx = (v - vp)/sp;
	    }
	    
	    getVal<float>(rtn, szX, szY, szZ, x, y, z) =  dx;
	    
	 }
      }
   }
   
}

void CFiniteDiff
::UpwindDiff(Image3D& h_o, const Image3D& h_i, 
	     const Image3D& h_speed,
	     DimT dim,
	     StreamT stream)
{
   MK_CHECK3_SIZE(h_o, h_i, h_speed);
   Vec3Di sz = h_i.size();
   Vec3Df sp = h_i.spacing();
   bool slice = (sz.z == 1);
   
   if(slice){
      if(dim == DIM_X){
	 PyCA::c_upwindDiff<DIM_X, SLICE_TRUE>
	    (h_o.get(), h_i.get(), h_speed.get(),
	     sz.x, sz.y, sz.z, 
	     sp.x, sp.y, sp.z);
      }else if(dim == DIM_Y){
	 PyCA::c_upwindDiff<DIM_Y, SLICE_TRUE>
	    (h_o.get(), h_i.get(), h_speed.get(),
	     sz.x, sz.y, sz.z, 
	     sp.x, sp.y, sp.z);
      }else if(dim == DIM_Z){
	 PyCA::c_upwindDiff<DIM_Z, SLICE_TRUE>
	    (h_o.get(), h_i.get(), h_speed.get(),
	     sz.x, sz.y, sz.z, 
	     sp.x, sp.y, sp.z);
      }else{
	 throw PyCAException(__FILE__, __LINE__, "unknown DimT");
      }
   }else{
      if(dim == DIM_X){
	 PyCA::c_upwindDiff<DIM_X, SLICE_FALSE>
	    (h_o.get(), h_i.get(), h_speed.get(),
	     sz.x, sz.y, sz.z, 
	     sp.x, sp.y, sp.z);
      }else if(dim == DIM_Y){
	 PyCA::c_upwindDiff<DIM_Y, SLICE_FALSE>
	    (h_o.get(), h_i.get(), h_speed.get(),
	     sz.x, sz.y, sz.z, 
	     sp.x, sp.y, sp.z);
      }else if(dim == DIM_Z){
	 PyCA::c_upwindDiff<DIM_Z, SLICE_FALSE>
	    (h_o.get(), h_i.get(), h_speed.get(),
	     sz.x, sz.y, sz.z, 
	     sp.x, sp.y, sp.z);
      }else{
	 throw PyCAException(__FILE__, __LINE__, "unknown DimT");
      }
   }
}

template<bool slice>
void 
c_upwindGradMag(float* rtn, 
		const float* array, 
		const float* speed, 
		int szX, int szY, int szZ,
		float spX, float spY, float spZ)
{

   for (int z=0; z<szZ; z++){
      for (int y=0; y<szY; y++){
	 for (int x=0; x<szX; x++){

	    float v = getVal<float>(array, szX, szY, szZ, x, y, z); // val
	    float n = getSafeVal<float,BG_CLAMP>(array, szX, szY, szZ, 
						       x, y-1, z); // north
	    float s = getSafeVal<float,BG_CLAMP>(array, szX, szY, szZ, 
						       x, y+1, z); // south
	    float e = getSafeVal<float,BG_CLAMP>(array, szX, szY, szZ, 
						       x+1, y, z); // east
	    float w = getSafeVal<float,BG_CLAMP>(array, szX, szY, szZ, 
						       x-1, y, z); // west
	    float u = 0.f;
	    float d = 0.f;
	    if(!slice){
	       u = getSafeVal<float,BG_CLAMP>(array, szX, szY, szZ, 
						    x, y, z+1); // up
	       d = getSafeVal<float,BG_CLAMP>(array, szX, szY, szZ, 
						    x, y, z-1); // down
	    }

	    float grad_before_x = (x < szX-1 ? e - v : v - w)/spX;
	    float grad_after_x = (x > 0 ? v - w : e - v)/spX;
	    float grad_before_y = (y < szY-1 ? s - v : v - n)/spY;
	    float grad_after_y = (y > 0 ? v - n : s - v)/spY;
	    float grad_before_z = 0.f;
	    float grad_after_z = 0.f;
	    if(!slice){
	       grad_before_z = (z < szZ-1 ? u - v : v - d)/spZ;
	       grad_after_z = (z > 0 ? v - d : u - v)/spZ;
	    }
	    
	    float spd = getVal<float>(speed, szX, szY, szZ, x, y, z); // speed
	    
	    if (spd < 0.0f)
	       {
		  grad_before_x = PYCAMIN(grad_before_x, 0.0f);
		  grad_after_x = PYCAMIN(-grad_after_x, 0.0f);
		  grad_before_y = PYCAMIN(grad_before_y, 0.0f);
		  grad_after_y = PYCAMIN(-grad_after_y, 0.0f);
		  if(!slice){
		     grad_before_z = PYCAMIN(grad_before_z, 0.0f);
		     grad_after_z = PYCAMIN(-grad_after_z, 0.0f);
		  }
	       }
	    else
	       {
		  grad_before_x = PYCAMAX(grad_before_x, 0.0f);
		  grad_after_x = PYCAMAX(-grad_after_x, 0.0f);
		  grad_before_y = PYCAMAX(grad_before_y, 0.0f);
		  grad_after_y = PYCAMAX(-grad_after_y, 0.0f);
		  if(!slice){
		     grad_before_z = PYCAMAX(grad_before_z, 0.0f);
		     grad_after_z = PYCAMAX(-grad_after_z, 0.0f);
		  }
	       }

	    float gradmag=0.f;
	    if(!slice){
	       gradmag = 
		  sqrt(grad_after_x*grad_after_x + grad_before_x*grad_before_x +
		       grad_after_y*grad_after_y + grad_before_y*grad_before_y +
		       grad_after_z*grad_after_z + grad_before_z*grad_before_z);
	    }else{
	       gradmag = 
		  sqrt(grad_after_x*grad_after_x + grad_before_x*grad_before_x +
		       grad_after_y*grad_after_y + grad_before_y*grad_before_y);
	    }
	    
	    getVal<float>(rtn, szX, szY, szZ, x, y, z) =  gradmag*spd;
	    
	 }
      }
   }
   
}

void CFiniteDiff
::UpwindGradMag(Image3D& h_o, const Image3D& h_i, 
		const Image3D& h_speed,
		StreamT stream)
{
   MK_CHECK3_SIZE(h_o, h_i, h_speed);
   Vec3Di sz = h_i.size();
   Vec3Df sp = h_i.spacing();
   bool slice = (sz.z == 1);
   
   if(slice){
      PyCA::c_upwindGradMag<SLICE_TRUE>(h_o.get(), h_i.get(), h_speed.get(),
				      sz.x, sz.y, sz.z, 
				      sp.x, sp.y, sp.z);
   }else{
      PyCA::c_upwindGradMag<SLICE_FALSE>(h_o.get(), h_i.get(), h_speed.get(),
				       sz.x, sz.y, sz.z, 
				       sp.x, sp.y, sp.z);
   }
   
}

} // end namespace PyCA

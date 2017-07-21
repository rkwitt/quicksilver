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

%{
#include "mem.h"
#include "pycaConst.h"
#include "IOpers.h"
#include "IFOpers.h"
#include "FOpers.h"
#include "HFOpers.h"
#include "VFOpers.h"
#include "Reduction.h"
#include "Vec2D.h"
#include "Vec3D.h"
%}

// include generated docstrings
%include "Opers.i"

%include "IOpers.h"
%include "IFOpers.h"
%include "FOpers.h"
%include "HFOpers.h"
%include "VFOpers.h"

 // for Reduction
%apply float& OUTPUT { float& FLOAT_OUT };

%ignore MaxMin;

%include "Reduction.h"

 // fix maxmin/sum functions to return multiple values
%apply float *OUTPUT{ float *MIN_OUT, float *MAX_OUT};

//
// Add MinMax function for mapped for python
//
namespace PyCA {
namespace Opers{

   void MinMax(float *MIN_OUT,float *MAX_OUT,const Image3D& a_i,
	       bool update=false,StreamT stream=NULL,bool onDev=false);

   void MinMax(float *MIN_OUT,float *MAX_OUT,const Field3D& a_i,
	       bool update=false,StreamT stream=NULL,bool onDev=false);

}
} 
  
%{
namespace PyCA {
namespace Opers{

   void MinMax(float *MIN_OUT,float *MAX_OUT,const Image3D& a_i,
	       bool update=false,StreamT stream=NULL,bool onDev=false)
   {
      Vec2D<float> maxmin;
      MaxMin(maxmin, a_i, update, stream, onDev);
      *MIN_OUT = maxmin.y;
      *MAX_OUT = maxmin.x;
   }

   void MinMax(float *MIN_OUT,float *MAX_OUT,const Field3D& a_i,
	       bool update=false,StreamT stream=NULL,bool onDev=false)
   {
      Vec2D<float> maxmin;
      MaxMin(maxmin, a_i, update, stream, onDev);
      *MIN_OUT = maxmin.y;
      *MAX_OUT = maxmin.x;
   }

}
}   
%}

// add version of Sum for Field3Ds that returns a Vec3Df

namespace PyCA {
namespace Opers{

    Vec3D<float> SumComp(const Field3D& a_i, StreamT stream=NULL);
}
}

%{
namespace PyCA {
namespace Opers{

    Vec3D<float> SumComp(const Field3D& a_i, StreamT stream=NULL)
   {
      Vec3D<float> sum;
      Sum(sum, a_i, stream);
      return sum;
   }

}
}   
%}


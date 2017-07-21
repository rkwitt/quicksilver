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

#ifndef __CFINITE_DIFF_H
#define __CFINITE_DIFF_H

#include<pycaConst.h>
#include<estream.h>
#include<Vec3D.h>
#include<pycaUtils.h>

namespace PyCA {

class Image3D;
class Field3D;

class CFiniteDiff{
public:
   
   /**
    * Compute finite difference in the dimension indicated
    *  a_o = diff(a_i)
    */
   static void
   FiniteDiff(float* h_o, const float* h_i,
	      int szX, int szY, int szZ,
	      float spX, float spY, float spZ, 
	      DimT dim,  DiffT diffType=DIFF_CENTRAL, 
	      enum BoundaryCondT bc=BC_CLAMP, 
	      bool accum=false, OpT op=OP_VAL,
	      StreamT stream=NULL);
   
   static void 
   FiniteDiff(Image3D& a_o, const Image3D& a_i, 
	      DimT dim,  DiffT diffType=DIFF_CENTRAL, 
	      enum BoundaryCondT bc=BC_CLAMP, 
	      bool accum=false, OpT op=OP_VAL,
	      StreamT stream=NULL);


   /**
    * returns deriv(h_i)*a_speed, where deriv is upwind
    * derivative based on a_speed in dimension 'dim'
    */
   static void 
   UpwindDiff(Image3D& h_o, const Image3D& h_i, 
	      const Image3D& h_speed,
	      DimT dim,
	      StreamT stream=NULL);

   /**
    * returns |\nabla a_i|*a_speed, where \nabla taken with upwind
    * differences based on a_speed
    */
   static void 
   UpwindGradMag(Image3D& h_o, const Image3D& h_i,
		 const Image3D& h_speed, StreamT stream=NULL);

};

} // end namespace PyCA

#endif // __CFINITE_DIFF_H

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

#ifndef __REDUCEOB_H__
#define __REDUCEOB_H__

#ifndef SWIG
#include <estream.h>
#endif // !SWIG
namespace PyCA {

class Image3D;
class Field3D;
template<typename T> class Vec2D;
template<typename T> class Vec3D;

template<int mode>
class ReduceOb {
public:
   // Operation on the image
   static void Max(float& FLOAT_OUT,const Image3D& a_i,
		   bool update,StreamT stream,bool onDev);
   static void Min(float& FLOAT_OUT,const Image3D& a_i,
		   bool update,StreamT stream,bool onDev);
   static void Sum(float& FLOAT_OUT,const Image3D& a_i,
		   bool update,StreamT stream,bool onDev);
   static void LInf(float& FLOAT_OUT,const Image3D& a_i,
		    bool update,StreamT stream,bool onDev);
   static void L1(float& FLOAT_OUT,const Image3D& a_i,
		  bool update,StreamT stream,bool onDev);
   static void Sum2(float& FLOAT_OUT,const Image3D& a_i,
		    bool update,StreamT stream,bool onDev);
   static void MaxMin(Vec2D<float>& a_o,const Image3D& a_i,
		      bool update,StreamT stream,bool onDev);
   static void Dot(float& FLOAT_OUT,const Image3D& a_i,const Image3D& a_i1,
		   bool update,StreamT stream,bool onDev);
   
   // Operation on the field
   static void Max(float& FLOAT_OUT,const Field3D& a_i,
		   bool update,StreamT stream,bool onDev);
   static void Min(float& FLOAT_OUT,const Field3D& a_i,
		   bool update,StreamT stream,bool onDev);
   static void Sum(float& FLOAT_OUT,const Field3D& a_i,
		   bool update,StreamT stream,bool onDev);
   static void Sum(Vec3D<float>& a_o,const Field3D& a_i, 
                   StreamT stream);
   static void LInf(float& FLOAT_OUT,const Field3D& a_i,
		    bool update,StreamT stream,bool onDev);
   static void L1(float& FLOAT_OUT,const Field3D& a_i,
		  bool update,StreamT stream,bool onDev);
   static void Sum2(float& FLOAT_OUT,const Field3D& a_i,
		    bool update,StreamT stream,bool onDev);
   static void Dot(float& FLOAT_OUT,const Field3D& a_i,const Field3D& a_i1,
		   bool update,StreamT stream,bool onDev);
   static void MaxMin(Vec2D<float>& a_o,const Field3D& a_i,
		      bool update,StreamT stream,bool onDev);
   
}; // end class ReduceOb

} // end namespace PyCA

#endif // __REDUCEOB_H__

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

#ifndef __REDUCTION_H
#define __REDUCTION_H

#ifndef SWIG
#include <estream.h>
#endif // !SWIG
namespace PyCA {

class Image3D;
class Field3D;
template<typename T> class Vec2D;
template<typename T> class Vec3D;

namespace Opers {

// Operation on the image
void Max(float& FLOAT_OUT,const Image3D& a_i,
	 bool update=false,StreamT stream=NULL,bool onDev=false);
void Min(float& FLOAT_OUT,const Image3D& a_i,
	 bool update=false,StreamT stream=NULL,bool onDev=false);
void Sum(float& FLOAT_OUT,const Image3D& a_i,
	 bool update=false,StreamT stream=NULL,bool onDev=false);
void LInf(float& FLOAT_OUT,const Image3D& a_i,
	  bool update=false,StreamT stream=NULL,bool onDev=false);
void L1(float& FLOAT_OUT,const Image3D& a_i,
	bool update=false,StreamT stream=NULL,bool onDev=false);
void Sum2(float& FLOAT_OUT,const Image3D& a_i,
	  bool update=false,StreamT stream=NULL,bool onDev=false);
void MaxMin(Vec2D<float>& a_o,const Image3D& a_i,
	    bool update=false,StreamT stream=NULL,bool onDev=false);
void Dot(float& FLOAT_OUT,const Image3D& a_i,const Image3D& a_i1,
	 bool update=false,StreamT stream=NULL,bool onDev=false);

// Operation on the field
void Max(float& FLOAT_OUT,const Field3D& a_i,
	 bool update=false,StreamT stream=NULL,bool onDev=false);
void Min(float& FLOAT_OUT,const Field3D& a_i,
	 bool update=false,StreamT stream=NULL,bool onDev=false);
void Sum(float& FLOAT_OUT,const Field3D& a_i,
	 bool update=false,StreamT stream=NULL,bool onDev=false);
void Sum(Vec3D<float>& a_o,const Field3D& a_i, 
	 StreamT stream=NULL);
void LInf(float& FLOAT_OUT,const Field3D& a_i,
	  bool update=false,StreamT stream=NULL,bool onDev=false);
void L1(float& FLOAT_OUT,const Field3D& a_i,
	bool update=false,StreamT stream=NULL,bool onDev=false);
void Sum2(float& FLOAT_OUT,const Field3D& a_i,
	  bool update=false,StreamT stream=NULL,bool onDev=false);
void Dot(float& FLOAT_OUT,const Field3D& a_i,const Field3D& a_i1,
	 bool update=false,StreamT stream=NULL,bool onDev=false);
void MaxMin(Vec2D<float>& a_o,const Field3D& a_i,
	    bool update=false,StreamT stream=NULL,bool onDev=false);

}; // end namespace Opers

} // end namespace PyCA

#endif

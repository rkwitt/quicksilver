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

#ifndef __FLUID_KERNEL_FFT_INTERFACE_H__
#define __FLUID_KERNEL_FFT_INTERFACE_H__

#ifndef SWIG

#include <estream.h>
#include <GridInfo.h>
#include <Field3D.h>

#endif // !SWIG

namespace PyCA {

// forward declaration

class FluidKernelFFTInterface {
public:
   // parameters:
   virtual ~FluidKernelFFTInterface() = 0;
   virtual void setAlpha(float alpha) = 0;
   virtual float getAlpha() = 0;
   virtual void setBeta(float beta) = 0;
   virtual float getBeta() = 0;
   virtual void setGamma(float gamma) = 0;
   virtual float getGamma() = 0;
   virtual void setLPow(int lpow) = 0;
   virtual int getLPow() = 0;
   virtual void setDivergenceFree(bool divfree) = 0;
   virtual bool getDivergenceFree() = 0;

   // execution:
   virtual void setGrid(const GridInfo& grid) = 0;
   
   virtual void applyOperator(Field3D& f, const Field3D& v, 
			      StreamT stream) = 0;
   
   virtual void applyInverseOperator(Field3D& v, const Field3D& f, 
				     StreamT stream) = 0;
   
   virtual void applyOperator(Field3D& f, StreamT stream) = 0;
   
   virtual void applyInverseOperator(Field3D& v, StreamT stream) = 0;

};

 inline FluidKernelFFTInterface::~FluidKernelFFTInterface(){}

} // end namespace PyCA

#endif // __FLUID_KERNEL_FFT_INTERFACE_H__

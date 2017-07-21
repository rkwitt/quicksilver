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

#ifndef __MULTISCALE_RESAMPLER_INTERFACE_H__
#define __MULTISCALE_RESAMPLER_INTERFACE_H__


#ifndef SWIG

#include <boost/shared_ptr.hpp>

// TEST -- make sure filed including boost aren't leaking into
// nvcc-compiled code
#if defined(PYCA_BOOSTTEST)
#if defined(__CUDACC__)
int bla[-1];
#endif
#endif
// END TEST

#include <MultiscaleManager.h>

#endif // !SWIG

namespace PyCA {

/**
 * Just a non-templated pure virtual interface for MultiscaleResampler
 */
class MultiscaleResamplerInterface
{
public:
    virtual ~MultiscaleResamplerInterface() = 0;
    virtual const MultiscaleManager &getScaleManager() const = 0;

   /**
    * Downsample image from orig to current scale
    */
   virtual void downsampleImage(Image3D& sI, 
				const Image3D& orgI, 
				StreamT stream=NULL) const = 0;
   
   /**
    * Upsample image to current scale
    */
   virtual void upsampleImage(Image3D& sI, 
			      StreamT stream=NULL) const = 0;
      
   /**
    * Downsample field from orig to current scale
    */
   virtual void downsampleField(Field3D& sF, 
				const Field3D& orgF, 
				StreamT stream=NULL) const = 0;

    virtual void updateImage(const Image3D*& I_ptr, 
			     const Image3D& orgI, 
			     Image3D& sI, 
			     StreamT stream=NULL) const = 0;
    virtual void updateHField(Field3D& h, 
			      StreamT stream=NULL) const = 0;
    virtual void updateVField(Field3D& v, 
			      StreamT stream=NULL) const = 0;
};

 inline MultiscaleResamplerInterface::~MultiscaleResamplerInterface(){}

} // end namespace PyCA

#endif // __MULTISCALE_RESAMPLER_INTERFACE_H__

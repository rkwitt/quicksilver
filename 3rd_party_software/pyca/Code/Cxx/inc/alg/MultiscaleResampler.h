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

#ifndef __MULTISCALE_RESAMPLER_H__
#define __MULTISCALE_RESAMPLER_H__


#include <MultiscaleResamplerInterface.h>

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

template<class FilterT>
class MultiscaleResampler :
    public MultiscaleResamplerInterface
{
public:
    typedef FilterT  SmoothFilterType;
    
    enum { exec_mode = FilterT::exec_mode};
    
    // Running parameter information
    MultiscaleResampler(const GridInfo& origGrid);

    const MultiscaleManager &getScaleManager() const 
    { return *mScaleManager; }

    void setScaleLevel(const MultiscaleManager &scaleManager,
		       StreamT stream=NULL);

   /**
    * Downsample image from orig to current scale
    */
   void downsampleImage(Image3D& sI, 
			const Image3D& orgI, 
			StreamT stream=NULL) const;

   /**
    * Downsample field from orig to current scale
    */
   void downsampleField(Field3D& sF, 
			const Field3D& orgF, 
			StreamT stream=NULL) const;
   
   /**
    * Upsample image to current scale
    */
   void upsampleImage(Image3D& sI, 
		      StreamT stream=NULL) const;
   
    void updateImage(const Image3D*& I_ptr, 
		     const Image3D& orgI, 
		     Image3D& sI, 
		     StreamT stream=NULL) const;
    void updateHField(Field3D& h, StreamT stream=NULL) const;
    void updateVField(Field3D& v, StreamT stream=NULL) const;
protected:    
    boost::shared_ptr<SmoothFilterType> mSmoothFilter;

    // just reference to scale manager
    const MultiscaleManager *mScaleManager;
    // original size/spacing info
    GridInfo mOrigGrid;
    // cur scale size/spacing info
    GridInfo mCurGrid;
};

} // end namespace PyCA

#endif // __MULTISCALE_RESAMPLER_H__

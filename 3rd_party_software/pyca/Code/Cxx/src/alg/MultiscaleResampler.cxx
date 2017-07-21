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


#include "MultiscaleResampler.h"

#include <MemoryManager.h>
#include <GaussianFilter.h>
#include <IdentityFilter.h>
#include "DownSampleParamUtils.h"

#include <pycaConst.h>
#include <IOpers.h>
#include <FOpers.h>
#include <ImageOpers.h>
#include <FieldOpers.h>
#include <VFieldOpers.h>
#include <HFieldOpers.h>

namespace PyCA {

template<class FilterT>
MultiscaleResampler<FilterT>::
MultiscaleResampler(const GridInfo& grid)
   :  mSmoothFilter(new FilterT()),
      mScaleManager(NULL),
      mOrigGrid(grid),
      mCurGrid(grid)
{
}

template<class FilterT>
void 
MultiscaleResampler<FilterT>::
setScaleLevel(const MultiscaleManager &scaleManager,
	      StreamT stream)
{
    // keep internal reference
    mScaleManager = &scaleManager;
    mCurGrid = scaleManager.getCurGrid();

    // Init smooth filter for downsampler from the orginal image
    if (!mScaleManager->isIdentityScale()) {
	Vec3Df sig;
	Vec3Di kRad;
	DownSampleParamUtils::
	    SigmaRadiusFromScale(sig, kRad, 
				 mScaleManager->getCurFactor());
        mSmoothFilter->updateParams(mOrigGrid.size(), sig, kRad);
    }
}

template<class FilterT>
void MultiscaleResampler<FilterT>::
downsampleImage(Image3D& sI, 
		const Image3D& orgI, 
		StreamT stream) const
{
   PYCA_ASSERT(sI.get());
   PYCA_ASSERT(mScaleManager);
   
   if (mScaleManager->isIdentityScale()){
      sI.setGrid(mCurGrid);
      Opers::Copy(sI, orgI, stream);
   } else {
      
      ManagedImage3D tempI0(mOrigGrid, sI.memType());
      ManagedImage3D tempI1(mOrigGrid, sI.memType());
      
      PYCA_ASSERT(orgI.grid() == tempI0.grid());
      sI.setGrid(mCurGrid);
      
      mSmoothFilter->filter(tempI0.get(), orgI.get(), tempI1.get(), stream);
      Opers::Resample<BACKGROUND_STRATEGY_CLAMP, 
		      DEFAULT_INTERP_METHOD, 
		      true>(sI, tempI0, stream);
      
   }
}

template<class FilterT>
void MultiscaleResampler<FilterT>::
downsampleField(Field3D& sF, 
		const Field3D& orgF, 
		StreamT stream) const
{
   PYCA_ASSERT(sF.getX());
   PYCA_ASSERT(sF.getY());
   PYCA_ASSERT(sF.getZ());
   PYCA_ASSERT(mScaleManager);
   
   if (mScaleManager->isIdentityScale()){
      sF.setGrid(mCurGrid);
      Opers::Copy(sF, orgF, stream);
   } else {
      
      ManagedImage3D tempI(mOrigGrid, sF.memType());
      ManagedField3D tempF(mOrigGrid, sF.memType());
      
      PYCA_ASSERT(orgF.grid() == tempI.grid());
      sF.setGrid(mCurGrid);
      
      mSmoothFilter->filter(tempF.x, orgF.x, tempI.get(), stream);
      mSmoothFilter->filter(tempF.y, orgF.y, tempI.get(), stream);
      mSmoothFilter->filter(tempF.z, orgF.z, tempI.get(), stream);
      Opers::Resample<BACKGROUND_STRATEGY_CLAMP, 
		      true>(sF, tempF, stream);
      
   }
}

template<class FilterT>
void MultiscaleResampler<FilterT>::
upsampleImage(Image3D& sI, 
	      StreamT stream) const
{
   //2. Upsamping the image from the previous
   // start with zero
   PYCA_ASSERT(sI.get());
   PYCA_ASSERT(mScaleManager);
   if (mScaleManager->isFirstScale()){
      sI.setGrid(mCurGrid);
      ImageOpers<exec_mode>::SetMem(sI, 0.f, stream);
   }else {
      
      PYCA_ASSERT((sI.grid() == mScaleManager->getPreGrid()));
      
      ManagedImage3D im(mCurGrid, sI.memType());
      Opers::Resample<BACKGROUND_STRATEGY_CLAMP, 
		      DEFAULT_INTERP_METHOD,
		      true>
	  (im, sI, stream);

      // worried about this with new ManagedImage3D
      //sI.swap(im);
      sI.setGrid(mCurGrid);
      Opers::Copy(sI, im);
      
   }

}


template<class FilterT>
void MultiscaleResampler<FilterT>::
updateImage(const Image3D*& I_ptr, 
	    const Image3D& orgI, 
	    Image3D& sI, 
	    StreamT stream) const
{
    PYCA_ASSERT(orgI.get());
    PYCA_ASSERT(mScaleManager);
    
    if (mScaleManager->isIdentityScale()){
        I_ptr = const_cast<Image3D*>(&orgI);
    } else {
       this->downsampleImage(sI, orgI, stream);
       I_ptr = &sI;
    }
}

template<class FilterT>
void MultiscaleResampler<FilterT>::
updateHField(Field3D& h, StreamT stream) const
{
    PYCA_ASSERT(h.getX());
    PYCA_ASSERT(h.getY());
    PYCA_ASSERT(h.getZ());
    PYCA_ASSERT(mScaleManager);

    
    //2. Upsamping the Hfield from the previous
    // start with identity
    if (mScaleManager->isFirstScale()) {
        h.setGrid(mCurGrid);
	HFieldOpers<exec_mode>::SetToIdentity(h, stream);
    }
    else {
        PYCA_ASSERT((h.grid() == mScaleManager->getPreGrid()));
    ManagedField3D tempH(mCurGrid, h.memType());

	HFieldOpers<exec_mode>::toV_I(h, stream);

        FieldOpers<exec_mode>::
	    template Resample<BACKGROUND_STRATEGY_CLAMP, 
			      false>(tempH, h, stream);

	// worried about this with new ManagedField3D
	// h.swap(tempH);
	h.setGrid(mCurGrid);
	Opers::Copy(h, tempH);

	VFieldOpers<exec_mode>::toH_I(h, 1.f, stream);
    }
}

template<class FilterT>
void MultiscaleResampler<FilterT>::
updateVField(Field3D& v, StreamT stream) const
{
    //2. Upsamping the VField from the previous
    // start with zero
    PYCA_ASSERT(v.getX());
    PYCA_ASSERT(v.getY());
    PYCA_ASSERT(v.getZ());
    PYCA_ASSERT(mScaleManager);
    if (mScaleManager->isFirstScale()){
        v.setGrid(mCurGrid);
	VFieldOpers<exec_mode>::SetToZero(v, stream);
    }else {

        PYCA_ASSERT((v.grid() == mScaleManager->getPreGrid()));
        
    ManagedField3D tempV(mCurGrid, v.memType());

	FieldOpers<exec_mode>::
            template Resample<BACKGROUND_STRATEGY_CLAMP, 
			      false>(tempV, v, stream);
	// worried about this with new ManagedField3D
        // v.swap(tempV);
	v.setGrid(mCurGrid);
	Opers::Copy(v, tempV);

    }
}

// template instantiations
template class MultiscaleResampler<GaussianFilter<EXEC_CPU> >;
template class MultiscaleResampler<IdentityFilter<EXEC_CPU> >;
#ifdef CUDA_ENABLED
template class MultiscaleResampler<GaussianFilter<EXEC_GPU> >;
template class MultiscaleResampler<IdentityFilter<EXEC_GPU> >;
#endif
} // end namespace PyCA

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


#include "MultiscaleManager.h"
#include "PyCAException.h"

#include <MemoryManager.h>
#include "Vec3D.h"

#include <pycaConst.h>

namespace PyCA {

// static methods for downsampling grid

PyCA::GridInfo 
Multiscale::
downScale(const PyCA::GridInfo &origGrid, int f)
{
   PyCA::Vec3Di newSize = componentMax(origGrid.size() / f, PyCA::Vec3Di(1,1,1));
   PyCA::Vec3Df factor = origGrid.size();
   factor /= newSize;
   PyCA::Vec3Df newSpacing = origGrid.spacing()*factor;
   PyCA::GridInfo newGrid = PyCA::GridInfo(newSize, newSpacing, origGrid.origin());
   return newGrid;
}


MultiscaleManager::
MultiscaleManager(const PyCA::GridInfo& grid)
    : mFactor(1.0,1.0,1.0)
{
    mPreGrid = mCurGrid = mOrgGrid = grid;
    mCurLevel = 0;
}

void MultiscaleManager::updateLevel(){
    int f = this->getScaleFactor();
    
    // Save the previous grid
    mPreGrid = mCurGrid;
    
    mCurGrid = Multiscale::downScale(mOrgGrid, f);
    
    mFactor = PyCA::Vec3Df(((float)mPreGrid.size().x)/mCurGrid.size().x,
			 ((float)mPreGrid.size().y)/mCurGrid.size().y,
			 ((float)mPreGrid.size().z)/mCurGrid.size().z);

    std::cout << "level is " << mCurLevel << ", grid is " << std::cout << mCurGrid << std::endl;

}

void 
MultiscaleManager::
addScaleLevel(unsigned int downsampleFactor)
{
   // insert the new scale
   std::vector<unsigned int>::iterator it = mScaleFactors.begin();
   for(;it!=mScaleFactors.end();++it){
      
      if(*it < downsampleFactor) break;
      
      if(*it == downsampleFactor){
	 throw PyCAException(__FILE__, __LINE__, "Error, cannot add multiple "
				      "scale levels with the same downsample factor");
      }
      
   }
   mScaleFactors.insert(it,downsampleFactor);

   it = mScaleFactors.begin();
   std::cout << "Scale factor list: ";
   for(;it!=mScaleFactors.end();++it){
      std::cout << *it << ", ";
   }
   std::cout << std::endl;

}

} // end namespace PyCA

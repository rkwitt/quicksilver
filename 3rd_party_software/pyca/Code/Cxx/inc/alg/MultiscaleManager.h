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

#ifndef __MULTISCALE_MANAGER_BASE_H__
#define __MULTISCALE_MANAGER_BASE_H__

#ifndef SWIG

#include <vector>

#include <Vec3D.h>
#include <GridInfo.h>
#include <estream.h>

#endif // !SWIG

namespace PyCA {

// forward declarations
class Image3D;
class Field3D;

// static method for downsampling grid
class Multiscale{
public:
    static GridInfo downScale(const GridInfo &orig, int factor);
};

class MultiscaleManager
{
public:
    
    // Running parameter information
    MultiscaleManager(const GridInfo& origGrid);

   void addScaleLevel(unsigned int downsampleFactor);

    void start() { 
       if(nScale() == 0) this->addScaleLevel(1);
       mCurLevel = 0; 
       updateLevel(); 
    }

    void next()  { if (!isLastScale()) { ++mCurLevel; updateLevel();}}

    void set(size_t level) 
    { 
	mCurLevel = level >= nScale() ? nScale()-1 : level; 
	updateLevel(); 
    }

    size_t getCurLevel()     const { return mCurLevel; };
    size_t nScale()         const { return mScaleFactors.size(); };
    size_t getScaleFactor()  const { return mScaleFactors[mCurLevel];};
    
    bool isLastScale()       const { return (mCurLevel == nScale() - 1);}
    bool isFirstScale()      const { return (mCurLevel == 0); }
    bool isIdentityScale()   const { return getScaleFactor() == 1; };

    const Vec3Df& getCurFactor() const { return mFactor; }

    // Grid information
    const GridInfo& getOrgGrid() const { return mOrgGrid; }
    const GridInfo& getCurGrid() const { return mCurGrid; }
    const GridInfo& getPreGrid() const { return mPreGrid; }

    const Vec3Di& getOrgSize()  const { return mOrgGrid.size(); }
    const Vec3Di& getCurSize()  const { return getCurGrid().size();};
    
    const Vec3Df& getOrgSpace() const { return mOrgGrid.spacing(); }
    const Vec3Df& getCurSpace() const { return getCurGrid().spacing(); }

    size_t origNVox()    const { return mOrgGrid.nVox();  }
    size_t curNVox()    const { return getCurGrid().nVox();  }
    float  origVoxelVol() const { return mOrgGrid.voxelVol(); }
    float  curVoxelVol() const { return getCurGrid().voxelVol(); }
    
   //void print(std::ostream& oss) const;

protected:    

    std::vector<unsigned int> mScaleFactors;
    size_t             mCurLevel; // current level
    GridInfo           mOrgGrid;
    GridInfo           mPreGrid;
    GridInfo           mCurGrid;
    Vec3Df             mFactor;
    virtual void updateLevel();
};

} // end namespace PyCA

#endif // __MULTISCALE_MANAGER_BASE_H__

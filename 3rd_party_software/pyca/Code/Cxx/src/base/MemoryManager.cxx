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

#include <MemPool.h>
#include <Image3D.h>
#include <Field3D.h>
#include <MemoryManager.h>
#include "PyCAException.h"

#include <algorithm>

#undef max
#undef min

namespace PyCA {

__thread MemoryManager* ThreadMemoryManager::mMemMan;

ManagedImage3D::ManagedImage3D(const GridInfo& g, MemoryType mt)
    : Image3D(mt) {
    if (ThreadMemoryManager::isInitialized()) {
        MemoryManager& mm = ThreadMemoryManager::instance();
        if (mm.memType() != mt) {  // can't use mm so don't manage
            mManager = NULL;
            mId = MM_ID_NONE;
            mData = MemPool<float>(g.nVox(), mt);
        } else {  // if types are compatible go ahead and use the mm
            boost::shared_ptr<MemPool<float> > p;
            mManager = &mm;
            mId = mm.getPool(g, p);
            mData = MemPool<float>(p->getSharedPtr(), g.nVox(), mt);
        }
    } else {  // default constructor
        mManager = NULL;
        mId = MM_ID_NONE;
        mData = MemPool<float>(g.nVox(), mt);
    }
    setGrid(g);  // have to use the public accessor here
}

ManagedImage3D::ManagedImage3D(const GridInfo& g, MemoryType mt, MemoryManager &mm)
    : Image3D(mt) {
    if (mm.memType() != mt) {  // can't use mm so don't manage
        throw PyCAException(__FILE__, __LINE__, "MemoryManager memType does not match");
    }

    // if types are compatible go ahead and use the mm
    boost::shared_ptr<MemPool<float> > p;
    mManager = &mm;
    mId = mm.getPool(g, p);
    mData = MemPool<float>(p->getSharedPtr(), g.nVox(), mt);

    setGrid(g);  // have to use the public accessor here
}

ManagedImage3D::~ManagedImage3D() {
    if (mManager) {  // This image is being managed, release it
        mManager->releasePool(mId);
    }
    mType = MEM_UNINITIALIZED;
}

ManagedField3D::ManagedField3D(const GridInfo& g, MemoryType mt)
    : Field3D(mt) {
    if (ThreadMemoryManager::isInitialized()) {
        MemoryManager& mm = ThreadMemoryManager::instance();
        if (mm.memType() != mt) {  // can't use mm so don't manage
            mManager = NULL;
            mIdX = mIdY = mIdZ = MM_ID_NONE;
            mDataX = MemPool<float>(g.nVox(), mt);
            mDataY = MemPool<float>(g.nVox(), mt);
            mDataZ = MemPool<float>(g.nVox(), mt);
        } else {  // if types are compatible go ahead and use the mm
            boost::shared_ptr<MemPool<float> > p;
            mManager = &mm;
            mIdX = mm.getPool(g, p);
            mDataX = MemPool<float>(p->getSharedPtr(), g.nVox(), mt);
            mIdY = mm.getPool(g, p);
            mDataY = MemPool<float>(p->getSharedPtr(), g.nVox(), mt);
            mIdZ = mm.getPool(g, p);
            mDataZ = MemPool<float>(p->getSharedPtr(), g.nVox(), mt);
        }
    } else {  // default constructor
        mManager = NULL;
        mIdX = mIdY = mIdZ = MM_ID_NONE;
        mDataX = MemPool<float>(g.nVox(), mt);
        mDataY = MemPool<float>(g.nVox(), mt);
        mDataZ = MemPool<float>(g.nVox(), mt);
    }
    setGrid(g);  // have to use the public accessor here
}

ManagedField3D::ManagedField3D(const GridInfo& g, MemoryType mt, MemoryManager& mm)
    : Field3D(mt) {
    if (mm.memType() != mt) {  // can't use mm so don't manage
        throw PyCAException(__FILE__, __LINE__, "MemoryManager memType does not match");
    }
    boost::shared_ptr<MemPool<float> > p;
    mManager = &mm;
    mIdX = mm.getPool(g, p);
    mDataX = MemPool<float>(p->getSharedPtr(), g.nVox(), mt);
    mIdY = mm.getPool(g, p);
    mDataY = MemPool<float>(p->getSharedPtr(), g.nVox(), mt);
    mIdZ = mm.getPool(g, p);
    mDataZ = MemPool<float>(p->getSharedPtr(), g.nVox(), mt);

    setGrid(g);  // have to use the public accessor here
}

ManagedField3D::~ManagedField3D() {
    if (mManager) {  // This image is being managed, release it
        mManager->releasePool(mIdX);
        mManager->releasePool(mIdY);
        mManager->releasePool(mIdZ);
    }
    mType = MEM_UNINITIALIZED;
}

MemoryManager::MemoryManager(const GridInfo& g, MemoryType type, size_t nP):
    mType(type), mGrid(g), mMaxUsed(0)
{
    for (size_t i=0; i < nP; ++i)
	this->addPool();
}

size_t MemoryManager::getPool(const GridInfo& g, boost::shared_ptr<MemPool<float> >& out){
    PRECONDITION(g.nVox()<= mGrid.nVox(), "Grid too large");
    PRECONDITION(mPool.size() < mPool.max_size(), 
		 "Maximum number of MemPools reached");

    size_t id = mPool.size();
    
    // see if there's an unused mempool already allocated
    for(size_t i=0;i<mLocked.size();++i){
	if(!mLocked[i]){
	    id = i;
	}
    }
    
    // if no unused mempools, allocate a new one
    if(id == mPool.size()){
	this->addPool();
    }
    
    out = mPool[id];
    mLocked[id] = true;
    
    mMaxUsed = std::max(mMaxUsed, id);
    
    return id;
}

void MemoryManager::acquirePool(size_t id){
    PYCA_ASSERT(id < this->getNumPools());
    mLocked[id] = true;
}

void MemoryManager::releasePool(size_t id){
    PYCA_ASSERT(id < this->getNumPools());
    mLocked[id] = false;
}
    
MemoryManager::~MemoryManager(){
    std::cerr << "Maximum temporary Pools " << mMaxUsed + 1 << std::endl;
}
  
void MemoryManager::addPool(){
    mPool.push_back(boost::shared_ptr<MemPool<float> >
		    (new MemPool<float>(mGrid.nVox(), mType)));
    mLocked.push_back(false);
}

// static 
void 
ThreadMemoryManager::
init(const GridInfo& g, MemoryType type, size_t nP)
{
   mMemMan = new MemoryManager(g, type, nP);
}

// static 
MemoryManager& 
ThreadMemoryManager::
instance()
{
   // if this assert is failing through SWIG wrapping, see note in
   // alg.i
   PYCA_ASSERT(mMemMan);
   return *mMemMan;
}

// static 
void 
ThreadMemoryManager::
destroy()
{
   if (mMemMan){
      delete mMemMan;
      mMemMan = NULL;
   }
}

} // end namespace PyCA

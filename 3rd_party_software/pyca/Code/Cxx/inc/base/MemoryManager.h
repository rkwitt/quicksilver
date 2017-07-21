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

#ifndef __MEMORY_MANAGER
#define __MEMORY_MANAGER


#ifndef SWIG

#include <boost/shared_ptr.hpp>
#include <mem.h>
#include <vector>
#include <GridInfo.h>
#include <Image3D.h>
#include <Field3D.h>

// TEST -- make sure filed including boost aren't leaking into
// nvcc-compiled code
#if defined(PYCA_BOOSTTEST)
#if defined(__CUDACC__)
int bla[-1];
#endif
#endif
// END TEST

#endif // SWIG

#include <PyCAThread.h>

#define MM_ID_NONE ((size_t)-1)

namespace PyCA {

class MemoryManager;
  
/**
 * Th type can be used exactly how an Image3D currently is.  If a memory
 * manager is provided, or if a ThreadMemoryManager has been initialized, it
 * will use managed memory.  Otherwise it will allocate and deallocate its own
 * memory.
 */
class ManagedImage3D : public Image3D {
  private:
    ManagedImage3D();  // hide default constructor
    ManagedImage3D(const ManagedImage3D& other);  // no copys or assignment either
    ManagedImage3D& operator=(const ManagedImage3D& other);
  public:
    // Constructor that either calls ThreadMemoryManager or creates an Image3D
    ManagedImage3D(const GridInfo& grid, MemoryType mType);
    ManagedImage3D(const GridInfo& grid, MemoryType mType, MemoryManager& mm);

    ~ManagedImage3D();  // Just a hook to release from the mm

  private:
    size_t               mId;
    MemoryManager*    mManager;
};


/**
 * Th type can be used exactly how a Field3D currently is.  If a memory
 * manager is provided, or if a ThreadMemoryManager has been initialized, it
 * will use managed memory.  Otherwise it will allocate and deallocate its own
 * memory.
 */
class ManagedField3D : public Field3D {
  private:
    ManagedField3D();  // hide default constructor
    ManagedField3D(const ManagedField3D& other);  // no copys or assignment either
    ManagedField3D& operator=(const ManagedField3D& other);
  public:
    // Constructor that either calls ThreadMemoryManager or creates an Field3D
    ManagedField3D(const GridInfo& grid, MemoryType mType);
    ManagedField3D(const GridInfo& grid, MemoryType mType, MemoryManager& mm);

    ~ManagedField3D();  // Just a hook to release from the mm

  private:
    size_t            mIdX, mIdY, mIdZ;
    MemoryManager*    mManager;
};


/**
 * MemoryManager holds 'scratch' Image3Ds and Field3Ds.  Algorithms
 * that need scratch space call CreateImage(...) or CreateField(...)
 * to acquire these resources.  Image/Field3Ds are created as needed.
 * Destroying the returned object (Image3DBuffer/Field3DBuffer)
 * returns control of the underlying Image/Field3D to the
 * MemoryManager without destroying the underlying data structure.
 */
class MemoryManager {
public:
    
    /**
     * Construct a MemoryManager with the given initial number of
     * pools, with the default grid size g and memory type.
     */
    MemoryManager(const GridInfo& g, MemoryType type, size_t nP);
    ~MemoryManager();
    
    size_t getPool(const GridInfo& g, boost::shared_ptr<MemPool<float> >& out);

    // accessors
    size_t getNumPools() const { return mPool.size();}
    MemoryType memType() const { return mType;}
    const GridInfo& grid() const { return mGrid;}
    
private:
    MemoryManager& operator =(const MemoryManager& rhs);
    MemoryManager(MemoryManager& rhs);
    
    void acquirePool(size_t id);
    void releasePool(size_t id);
    
    void addPool();
    
    MemoryType mType;
    GridInfo   mGrid;
    
    std::vector<boost::shared_ptr<MemPool<float> > > mPool;
    std::vector<bool> mLocked;
    
    size_t mMaxUsed;
    
    friend class ManagedImage3D;
    friend class ManagedField3D;
};

/**
 * Creates a per-thread MemoryManager
 */
class ThreadMemoryManager {
public:
   static void init(const GridInfo& g, MemoryType type, size_t nP);
   
   static bool isInitialized() { return (mMemMan != NULL); };
   static MemoryManager& instance();
   
   static void destroy();

private:
    static __thread MemoryManager* mMemMan;
};

} // end namespace PyCA

#endif

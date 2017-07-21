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

#ifndef __MEM_POOL_H
#define __MEM_POOL_H

#include <mem.h>

#ifndef SWIG
#include <boost/shared_ptr.hpp>
#endif // !SWIG

// TEST -- make sure filed including boost aren't leaking into
// nvcc-compiled code
#if defined(PYCA_BOOSTTEST)
#if defined(__CUDACC__)
int bla[-1];
#endif
#endif
// END TEST

namespace PyCA {

/*
 * memory object with the information of capacity and memory type
 * use shared_ptr for resource manager for the data
 */
template<typename T>
class MemPool{
public:
    typedef T element_type;

    MemPool(const MemPool& rhs): mType(rhs.memType()),
                                 mCap(rhs.capacity()),
                                 mData(rhs.mData){};

    explicit MemPool(MemoryType type);
    MemPool(size_t n, MemoryType type);
    MemPool(element_type* p, size_t n, MemoryType type);
    MemPool(boost::shared_ptr<element_type> const r, size_t n, MemoryType type);
public:
    void print(std::ostream& oss) const;

    // get memory type
    MemoryType memType() const { return mType; };
    
    // get the basic properties
    size_t capacity() const { return mCap; };

    // get the raw pointer to data
    const element_type* get() const { return mData.get(); };
    element_type* get()             { return mData.get(); };


    // get element
    const T& get(size_t idx) const { 
       PYCA_ASSERT(mType != MEM_DEVICE); // Should we allow MEM_HOST_PINNED?
       return mData.get()[idx]; 
    };
    T& get(size_t idx) { 
      PYCA_ASSERT(mType != MEM_DEVICE); // Should we allow MEM_HOST_PINNED?
      return mData.get()[idx]; 
    };

#ifndef SWIG
    const T& operator[] (size_t idx) const {
       return this->get(idx);
    }

    T& operator[] (size_t idx) {
       return this->get(idx);
    }
#endif // !SWIG

    // get the smart pointer to data
    boost::shared_ptr<element_type>& getSharedPtr()             { return mData; };
    const boost::shared_ptr<element_type>& getSharedPtr() const { return mData; };

    void swap(MemPool<T>& rhs);

    // convert to different memory type
    void toType(MemoryType memType, StreamT stream = NULL);
    // resize the capacity, destroying current contents
   void resize(size_t cap, bool preserveData = false, StreamT stream = NULL);
   void clone(const MemPool<T> &other, StreamT stream=NULL);

private:
    MemoryType     mType;
    size_t               mCap;
    boost::shared_ptr<T> mData;

   boost::shared_ptr<T> allocMem(size_t nEl, MemoryType memType);
};

} // end namespace PyCA

#endif

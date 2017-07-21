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
#include <pycaConst.h>
#include <pycaUtils.h>
#include <iostream>
#include <mem.h>
#include <CudaUtils.h>
#include <PyCAException.h>

namespace PyCA {

// private function for memory allocation
template<typename T>
boost::shared_ptr<T> 
MemPool<T>::allocMem(size_t nEl, MemoryType memType)
{
   boost::shared_ptr<T> p;
    if (memType == MEM_HOST) {
        p = boost::shared_ptr<T>(ArrayAllocator<MEM_HOST, T>(nEl),
                                     ArrayDeleter<MEM_HOST, T>);
    } else if (memType == MEM_HOST_PINNED) {
        p = boost::shared_ptr<T>(ArrayAllocator<MEM_HOST_PINNED, T>(nEl),
                                     ArrayDeleter<MEM_HOST_PINNED, T>);
    } else if (memType == MEM_DEVICE) {
        p = boost::shared_ptr<T>(ArrayAllocator<MEM_DEVICE, T>(nEl),
                                     ArrayDeleter<MEM_DEVICE, T>);
    } else {
       throw PyCAException(__FILE__,__LINE__,
				    "Error, Unknown Memory Type");
    }

    if(p.get() == NULL){
       throw PyCAException(__FILE__,__LINE__,
				    "Error, memory allocation failed");
    }

    return p;
}

// Round up the memory to the aligned size
template<typename T>
MemPool<T>::MemPool(MemoryType type)
    :mType(type), mCap(0), mData((T*)NULL)
{
   
}

template<typename T>
MemPool<T>::MemPool(size_t n, MemoryType type)
    :mType(type), mCap(iAlignUp(n, BLOCK_ALIGN))
{
   mData = this->allocMem(mCap, mType);
}

template<typename T>
void MemPool<T>::print(std::ostream& oss) const {
    oss << "Address " << mData.get() << " capacity "; this->capacity();
    oss << " (in bytes "  << this->capacity() * sizeof(T) << ")" << std::endl;
}

template<typename T>
MemPool<T>::MemPool(T* p, size_t n, MemoryType type)
    :mType(type), mCap(n)
{
    if (mType == MEM_HOST) {
        mData = boost::shared_ptr<T>(p, ArrayDeleter<MEM_HOST, T>);
    } else if (mType == MEM_HOST_PINNED) {
        mData = boost::shared_ptr<T>(p, ArrayDeleter<MEM_HOST_PINNED, T>);
    } else {
        mData = boost::shared_ptr<T>(p, ArrayDeleter<MEM_DEVICE, T>);
    }
}

template<typename T>
MemPool<T>::MemPool(boost::shared_ptr<T> const r, size_t n, MemoryType type)
    :mType(type), mCap(n), mData(r)
{
    
}

template<typename T>
void MemPool<T>::swap(MemPool<T>& rhs) {
    PYCA_ASSERT(mType == rhs.mType);
    std::swap(mCap, rhs.mCap);
    mData.swap(rhs.mData);
}

template<typename T>
void MemPool<T>::resize(size_t cap, bool preserveData, StreamT stream)
{
    if(cap == mCap) return;

    // TEST
    CudaUtils::CheckCUDAError(__FILE__, __LINE__);
    // END TEST

    boost::shared_ptr<T> tmp = this->allocMem(cap, mType);
    
    if(preserveData){
       memCopy(tmp.get(), mData.get(), std::min(cap,mCap), mType, mType, stream);
    }
    
    mData = tmp;
    mCap = cap;
}

template<typename T>
void MemPool<T>::toType(MemoryType memType, StreamT stream)
{
    if(memType == mType) return;

    // TEST
    CudaUtils::CheckCUDAError(__FILE__, __LINE__);
    // END TEST

    boost::shared_ptr<T> tmp = this->allocMem(mCap, memType);
    
    memCopy(tmp.get(), mData.get(), mCap, memType, mType, stream);
    
    mData = tmp;
    mType = memType;
}

template<typename T>
void MemPool<T>::clone(const MemPool<T> &other, StreamT stream)
{
   if(&other == this) return;

   MemoryType otherType = other.memType();
   size_t otherCap = other.capacity();
   if(mType != otherType || mCap != otherCap){

      boost::shared_ptr<T> tmp = this->allocMem(otherCap, otherType);

      mData = tmp;
      mType = otherType;
      mCap = otherCap;
   }

   memCopy(this->get(), other.get(), mCap, mType, mType, stream);
   
}

template class MemPool<float>;
template class MemPool<short>;
template class MemPool<int>;
template class MemPool<unsigned int>;

} // end namespace PyCA

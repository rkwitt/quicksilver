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

#ifndef __BOOST_MEM_INC
#define __BOOST_MEM_INC

namespace PyCA {

template<MemoryType mType, class T>
inline boost::shared_ptr<T> CreateSharedArray(size_t n) {
    return boost::shared_ptr<T>(ArrayAllocator<mType, T>(n), ArrayDeleter<mType, T>);
}

// Compatible function between regular pointer and shared_pointer
// So that we can exchange the implementation between
// shared pointer and regular pointer
template<typename T>
inline void memAlloc(boost::shared_ptr<T>& p, size_t n, MemoryType mType){
    if (mType == MEM_HOST)
        p = CreateSharedArray<MEM_HOST, T>(n);
    else if (mType == MEM_HOST_PINNED)
        p = CreateSharedArray<MEM_HOST_PINNED, T>(n);
    else
        p = CreateSharedArray<MEM_DEVICE, T>(n);
}

template<typename T>
inline void memFree(boost::shared_ptr<T>& p, MemoryType mType){
    // Just to keep to make it compatible with raw pointer
}

} // end namespace PyCA

#endif // __MEM_INC

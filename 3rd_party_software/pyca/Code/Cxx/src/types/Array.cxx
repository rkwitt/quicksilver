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

#include <Array.h>
#include <pycaConst.h>
#include <pycaUtils.h>
#include <boost_mem.h>

namespace PyCA {

template<typename T>
Array<T>::Array(size_t n, MemoryType type):mType(type),
                                           mN(n),
                                           mCap(iAlignUp(n, BLOCK_ALIGN)){
    if (mType == MEM_HOST)
        mData = CreateSharedArray<MEM_HOST, T>(mCap);
    else if (mType == MEM_HOST_PINNED)
        mData = CreateSharedArray<MEM_HOST_PINNED, T>(mCap);
    else
        mData = CreateSharedArray<MEM_DEVICE, T>(mCap);
}

template class Array<float>;
template class Array<int>;
template class Array<unsigned int>;
} // end namespace PyCA

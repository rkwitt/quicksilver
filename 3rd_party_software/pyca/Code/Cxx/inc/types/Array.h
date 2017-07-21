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

#ifndef __ARRAY_H
#define __ARRAY_H

#include <iosfwd>
#include <boost/shared_ptr.hpp>
#include <mem.h>

namespace PyCA {

template<typename T>
class Array{
private:
public:
    typedef T element_type;
    Array(MemoryType type):mType(type), mN(0), mCap(0){
    }
    Array(size_t n, MemoryType type);
public:
    void print(std::ostream& oss) const;

    // get the basic properties
    size_t size() const { return mN; };
    size_t capacity() const { return mCap; };

    // get the raw pointer to data
    const T* get() const { return mData.get(); };
    T* get()             { return mData.get(); };
private:
    MemoryType mType;
    size_t     mN;
    size_t     mCap;
    boost::shared_ptr<T> mData;
};

} // end namespace PyCA

#endif

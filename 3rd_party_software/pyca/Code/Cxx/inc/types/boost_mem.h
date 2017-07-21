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

#ifndef __BOOST_MEM_H
#define __BOOST_MEM_H

#ifndef SWIG
#include <boost/shared_ptr.hpp>
#endif // SWIG

#include "mem.h"

namespace PyCA {

template<MemoryType mType, class T>
boost::shared_ptr<T> CreateSharedArray(size_t n);

template<typename T>
void memAlloc(boost::shared_ptr<T>& p, size_t n, MemoryType mType);

template<typename T>
void memFree(boost::shared_ptr<T>& p, MemoryType mType);

} // end namespace PyCA

#ifndef SWIG
#include "boost_mem.txx"
#endif // SWIG

#endif // __BOOST_MEM_H

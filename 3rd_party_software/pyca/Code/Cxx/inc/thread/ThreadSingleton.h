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

#ifndef __THREAD_SINGLETON_H
#define __THREAD_SINGLETON_H

#define MAX_NUMBER_DEVICES 256
#include <ThreadIdentity.h>
#include <boost/shared_ptr.hpp>

// TEST -- make sure filed including boost aren't leaking into
// nvcc-compiled code
#if defined(PYCA_BOOSTTEST)
#if defined(__CUDACC__)
int bla[-1];
#endif
#endif
// END TEST

namespace PyCA {

template <typename T>
class ThreadSingleton
{
public:
    static T& Instance() {
        size_t id = GetThreadId();
        if (insA[id].get()==NULL) {
            insA[id] = boost::shared_ptr<T>(new T);
        }
        return *insA[id];
    }
private:
    ThreadSingleton();
    ~ThreadSingleton();
    ThreadSingleton(ThreadSingleton const&);    // copy ctor hidden
    ThreadSingleton& operator=(ThreadSingleton const&);  // assign op hidden
    
    static boost::shared_ptr<T> insA[MAX_NUMBER_DEVICES];
};


template<typename T>
boost::shared_ptr<T> ThreadSingleton<T>::insA[MAX_NUMBER_DEVICES];

} // end namespace PyCA

#endif

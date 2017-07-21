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

#ifndef __THREAD_SINGLETON_COMPLEX_PARAMETER_ONE_H
#define __THREAD_SINGLETON_COMPLEX_PARAMETER_ONE_H

#define MAX_NUMBER_DEVICES 256

#include <ThreadIdentity.h>
#include <boost/shared_ptr.hpp>
#include <typeinfo>

// TEST -- make sure filed including boost aren't leaking into
// nvcc-compiled code
#if defined(PYCA_BOOSTTEST)
#if defined(__CUDACC__)
int bla[-1];
#endif
#endif
// END TEST

namespace PyCA {

 // struct ObjectGenerator {
 //    ObjectGenerator(complexParameter); // seting the parametter to create object
 //    Object* operator()() {
 //       return new Object(complexParameter);
 //    }
 // }
 // ..........
 // MemoryManager& objIns = ThreadSingletonC<Object>::Instance();
 // try {
 //        ThreadSingletonC<Object>::Instance();
 //    }
 //    catch (const char* e){
 //        std::cerr<< "Error " << e << std::endl;
 //        ThreadSingletonC<MemoryManager>::Create(ObjectGenerator(complexParameter));
 //    }
 // MemoryManager& objIns = ThreadSingletonC<Object>::Instance();
    

template<typename T>
class ThreadSingletonC
{
public:
    static T& Instance() {
        size_t id = GetThreadId();
        if (insA[id].get()==NULL) {
            throw "Singleton has not been create Call Instance(Generator) instead";
        }
        return *insA[id];
    }

    template<class Gen>
    static T& Create(Gen g) {
        size_t id = GetThreadId();
        std::cerr << "Thread id " << id << std::endl;
        if (insA[id].get()!=NULL)
            throw "Singleton has been create. Call Instance() instead";

        insA[id] = boost::shared_ptr<T>(g());
        return *insA[id];
    }

private:
    ThreadSingletonC();
    ~ThreadSingletonC();
    ThreadSingletonC(ThreadSingletonC const&);    // copy ctor hidden
    ThreadSingletonC& operator=(ThreadSingletonC const&);  // assign op hidden
    
    static boost::shared_ptr<T> insA[MAX_NUMBER_DEVICES];
};


template<typename T>
boost::shared_ptr<T> ThreadSingletonC<T>::insA[MAX_NUMBER_DEVICES];

} // end namespace PyCA

#endif

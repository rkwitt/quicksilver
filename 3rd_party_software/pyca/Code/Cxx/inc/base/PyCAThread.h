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


#ifndef __PYCA_THREAD_H__
#define __PYCA_THREAD_H__

// swig doesn't like __thread
#ifdef SWIG
#define __thread
#endif
// MSVC prefers __declspec(thread) instead of __thread
#ifdef _MSC_VER
#define __thread __declspec(thread)
#endif

// the rest of this isn't needed by SWIG
#ifndef SWIG

#include <string>
#include <exception>

#if defined(_MSC_VER)
#include <windows.h>
#elif defined(__GNUC__) || defined(__clang__)
#include <pthread.h>
#endif

namespace PyCA {

#if defined _MSC_VER
typedef DWORD ThreadID;
typedef HANDLE Mutex;
#define PYCA_MUTEX_INITIALIZER CreateMutex(NULL,FALSE,NULL)
#elif defined __GNUC__ || defined __clang__
typedef pthread_t ThreadID;
typedef pthread_mutex_t Mutex;
#define PYCA_MUTEX_INITIALIZER PTHREAD_MUTEX_INITIALIZER
#endif

ThreadID current_thread_id();
bool compare_threads(ThreadID t1, ThreadID t2);
void mutex_lock(Mutex *m);
void mutex_unlock(Mutex *m);

} // end namespace PyCA

#endif // SWIG

#endif // __PYCA_THREAD_H__

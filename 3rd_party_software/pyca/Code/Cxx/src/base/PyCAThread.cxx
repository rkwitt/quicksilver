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

#if defined(_MSC_VER)
#include <windows.h>
#elif defined(__GNUC__)
#include <pthread.h>
#endif

#include "PyCAThread.h"

namespace PyCA {

#if defined(_MSC_VER)
ThreadID current_thread_id(){ return GetCurrentThreadId(); }
bool compare_threads(ThreadID t1, ThreadID t2){ return t1 == t2; }

void mutex_lock(Mutex *m){ WaitForSingleObject(*m, INFINITE); }
void mutex_unlock(Mutex *m){ ReleaseMutex(*m); }
#elif defined(__GNUC__) || defined(__clang__)
ThreadID current_thread_id(){ return pthread_self(); }
bool compare_threads(ThreadID t1, ThreadID t2){ return pthread_equal(t1, t2); }

void mutex_lock(Mutex *m){ pthread_mutex_lock(m); }
void mutex_unlock(Mutex *m){ pthread_mutex_unlock(m); }
#endif

}

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

#include <ThreadIdentity.h>
#include <PyCAThread.h>

#if defined _MSC_VER

#include <windows.h>
#define ATOMICINC(x) InterlockedIncrement(&x);

#elif defined __clang__

#include <libkern/OSAtomic.h>
#define ATOMICINC(x) OSAtomicIncrement32Barrier(reinterpret_cast<int*>(&x));

#elif defined __GNUC__

#if (__GNUC__ > 4) || ((__GNUC__ == 4) && (__GNUC_MINOR__ >= 2))
#include <ext/atomicity.h>
#else
#include <bits/atomicity.h>
#endif
#define ATOMICINC(x) __sync_fetch_and_add(&x, 1);

#endif

namespace PyCA {

size_t& GetThreadIdRef(){
    static __thread size_t perThreadId=0;
    return perThreadId;
}

void InitThreadId() {
    static size_t threadCnt = 0 ;
    size_t& id = GetThreadIdRef();
    if (id == 0) {
        id = ATOMICINC(threadCnt);
    }
}

size_t GetThreadId() {
    return GetThreadIdRef();
}
} // end namespace PyCA

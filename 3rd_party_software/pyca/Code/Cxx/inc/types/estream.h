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

#ifndef __EXECUTION_STREAM_H
#define __EXECUTION_STREAM_H


#ifndef SWIG

#include <cstdlib>

#ifdef CUDA_ENABLED
#include <cuda_runtime.h>
#endif

#endif // SWIG

namespace PyCA {

#ifdef CUDA_ENABLED
typedef cudaStream_t StreamT;
#else // CUDA_ENABLED
typedef void* StreamT;
#endif // CUDA_ENABLED


bool CheckStreamCondition(size_t nstreams);
void StreamCreate(int nStream);
void StreamDestroy();
StreamT GetStream(int id);

#define STM_NULL NULL
#define STM_H2D  (GetStream(0))
#define STM_D2D  (GetStream(1))
#define STM_D2H  (GetStream(2))
#define STM_Di2H (GetStream(3))
#define STM_H2Di (GetStream(4))

#define STREAM_WARNING

void StreamIncompatibleWarning(const char* what, const char* file, int line);
void StreamIncompatibleWarning(const char* what, const char* file, int line, const char* sol);

} // end namespace PyCA

#endif

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

#include <estream.h>

#include <iostream>

#include <conditionMacro.h>

#include <PyCAThread.h>

#ifdef CUDA_ENABLED
#include <cuda_runtime_api.h>
#endif

namespace PyCA {

// This variable is global per-thread
__thread unsigned int streamCnt = 0;
__thread StreamT* streams = NULL;

bool CheckStreamCondition(size_t nstreams) {
    return (nstreams <= streamCnt);
}

StreamT GetStream(int id)
{
   if (id == 0 && streamCnt == 0)
      return StreamT (NULL);
   else {
      PYCA_ASSERT(id < (int)streamCnt);
      return streams[id];
   }
}

#ifdef CUDA_ENABLED

void StreamCreate(int n)
{
    streamCnt = n;
    streams   = new StreamT [n];
    for(int i = 0; i < n; i++)
        cudaStreamCreate(&(streams[i]));
}

void StreamDestroy()
{
   for (int i=0; i< (int)streamCnt; ++i)
      cudaStreamDestroy(streams[i]);
   delete []streams;
   streamCnt = 0;
   streams = NULL;
}

#else // CUDA_ENABLED

void StreamCreate(int n)
{
   throw PyCAException(__FILE__, __LINE__, "Cuda support not compiled");
}

void StreamDestroy()
{
   throw PyCAException(__FILE__, __LINE__, "Cuda support not compiled");
}

#endif // CUDA_ENABLED

void StreamIncompatibleWarning(const char* what){
    std::cerr << " WARNING: This function " << what << " does not support stream mode " << std::endl;
}

void StreamIncompatibleWarning(const char* what, const char* file, int line) {
    std::cerr << " WARNING: at " << file << ":" << line << " the function " << what << "does not support stream mode " << std::endl;
}

void StreamIncompatibleWarning(const char* what, const char* file, int line, const char* sol) {
    StreamIncompatibleWarning(what, file, line);
    std::cerr << " Call function " << sol << " instead " << std::endl;
}
} // end namespace PyCA

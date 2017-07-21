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


#include <string>

#include "CudaUtils.h"
#include "StringUtils.h"
#include "PyCAException.h"

#ifdef CUDA_ENABLED
#include <cuda.h>
#include <cuda_runtime.h>
#endif // CUDA_ENABLED

#include <log.h>

namespace PyCA {

#ifdef CUDA_ENABLED

unsigned int 
CudaUtils::
GetNumberOfCUDADevices()
{
    int deviceCount =0;
    cudaGetDeviceCount(&deviceCount);
    return deviceCount;
}

void 
CudaUtils::
SetCUDADevice(unsigned int deviceNum)
{
    // determine the number of GPUs
    int systemGPUs =0;
    cudaGetDeviceCount(&systemGPUs);
    if(systemGPUs < 1){
	throw PyCAException(__FILE__, __LINE__, "Error, no CUDA capable devices found.");
    }
    if((int)deviceNum >= systemGPUs){
	std::string err = StringUtils::strPrintf("Error, cannot select device %d, only %d devices available.",deviceNum, systemGPUs);
	throw PyCAException(__FILE__, __LINE__, err);
    }
    
    cudaSetDevice(deviceNum);
}

void 
CudaUtils::
CheckCUDAError(const char *file, int line, const char *msg)
{
   static const char *defaultMsg = "Error ";
   if(msg == NULL) msg = defaultMsg;
   std::string errMsg = std::string(msg);
   cudaError_t err;
   if (ERRLOG_MAX_LEVEL < logDEBUG){
      err = cudaGetLastError();
   }else{
      err = cudaThreadSynchronize();
   }
   if(cudaSuccess != err){
      if (ERRLOG_MAX_LEVEL < logDEBUG){
	 errMsg = std::string(" -- Cuda last error: ");
      }else{
	 errMsg = std::string("-- Cuda synch error: ");
      }
      errMsg = errMsg + cudaGetErrorString(err);
      throw PyCAException(file, line, errMsg.c_str());
   }
}

std::string
CudaUtils::
GetCUDAMemUsage()
{
    size_t free, total;
    cuMemGetInfo(&free, &total);
    std::string rtn = StringUtils::strPrintf("CUDA Memory Usage: used %u of %u", (total - free) >> 20, total >> 20);
    return rtn;
}

void
CudaUtils::
AssertMinCUDACapabilityVersion(int devNum, int major, int minor)
{
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, devNum);
  if(deviceProp.major == 9999 && deviceProp.minor == 9999){
    std::string err = StringUtils::strPrintf("error, cuda-capable device %d not found", devNum);
    throw PyCAException(__FILE__, __LINE__, err);
  }
  if(deviceProp.major < major || (deviceProp.major == major && deviceProp.minor < minor)){
    std::string err = StringUtils::strPrintf("error, cuda device %d supports cuda capability version %d.%d, "
					     "but %d.%d required", 
					     devNum, deviceProp.major, deviceProp.minor, major, minor);
    throw PyCAException(__FILE__, __LINE__, err);
  }
}

void
CudaUtils::
AssertMinCUDACapabilityVersion(int major, int minor)
{
  int devNum;
  cudaGetDevice(&devNum);
  AssertMinCUDACapabilityVersion(devNum, major, minor);
}

#else // CUDA_ENABLED

// implementations that should relate the fact that there is no cuda support

unsigned int 
CudaUtils::
GetNumberOfCUDADevices()
{
    return 0;
}

void 
CudaUtils::
SetCUDADevice(unsigned int deviceNum)
{
   throw PyCAException(__FILE__, __LINE__, "Error, no CUDA capable devices found (GPU support not compiled)");
}

void 
CudaUtils::
CheckCUDAError(const char *file, int line, const char *msg)
{
}

std::string
CudaUtils::
GetCUDAMemUsage()
{
   std::string rtn("CUDA support not compiled");
   return rtn;
}

void
CudaUtils::
AssertMinCUDACapabilityVersion(int devNum, int major, int minor)
{
    throw PyCAException(__FILE__, __LINE__, "Cuda support not compiled");
}

void
CudaUtils::
AssertMinCUDACapabilityVersion(int major, int minor)
{
    throw PyCAException(__FILE__, __LINE__, "Cuda support not compiled");
}

#endif // CUDA_ENABLED

} // end namespace PyCA


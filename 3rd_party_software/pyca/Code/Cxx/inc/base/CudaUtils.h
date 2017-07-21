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


#ifndef __CUDA_UTILS_H__
#define __CUDA_UTILS_H__

namespace PyCA {

namespace CudaUtils {

    /**
     * Get the number of cuda-capable devices on the current machine
     */
    unsigned int GetNumberOfCUDADevices();
    
    /**
     * Set the CUDA device to be used, throw exception if it is not a
     * valid device.
     */
    void SetCUDADevice(unsigned int deviceNum);

    /**
     * Check for a CUDA error, and throw an exception
     * (PyCAException) if one is found.  If compiled with
     * debug-level logging enabled, this will call
     * cudaThreadSynchronize, which may significantly affect
     * performance.
     */
   void CheckCUDAError(const char *file, int line, const char * = NULL);

    /**
     * Return CUDA memory usage as a string
     */
    std::string GetCUDAMemUsage();

    /**
     * Test for minimum cuda capability version, throw an exception if
     * it is not met
     */
    void AssertMinCUDACapabilityVersion(int devNum, int major, int minor);

    /**
     * Test current device for minimum cuda capability version, throw an
     * exception if it is not met
     */
    void AssertMinCUDACapabilityVersion(int major, int minor);

} // end namespace CudaUtils

} // end namespace PyCA

#endif //  __CUDA_UTILS_H__

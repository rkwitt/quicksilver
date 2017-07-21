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

#ifndef __GGAUSS_UTILS_H
#define __GGAUSS_UTILS_H

#include <estream.h>

namespace PyCA {

namespace GGaussUtils {
    void setConvolutionKernelX(const float *h_Kernel, size_t kLength);
    void setConvolutionKernelY(const float *h_Kernel, size_t kLength);
    void setConvolutionKernelZ(const float *h_Kernel, size_t kLength);
    
    void setSupplementKernelX(const float *h_Kernel, size_t kRadius);
    void setSupplementKernelY(const float *h_Kernel, size_t kRadius);
    void setSupplementKernelZ(const float *h_Kernel, size_t kRadius);

    void ConvolutionX3D(float* d_o, const float* d_i, int kRadius,
                        int sizeX, int sizeY, int sizeZ, StreamT stream);
    void ConvolutionY3D(float* d_o, const float* d_i, int kRadius,
                        int sizeX, int sizeY, int sizeZ, StreamT stream);
    void ConvolutionZ3D(float* d_o, const float* d_i, int kRadius,
                        int sizeX, int sizeY, int sizeZ, StreamT stream);

    void ConvolutionX3D(float* d_o, const float* d_i,
                        const float* d_Kernel, const float* d_sKernel, int kRadius,
                        int sizeX, int sizeY, int sizeZ, StreamT stream);
}

} // end namespace PyCA

#endif

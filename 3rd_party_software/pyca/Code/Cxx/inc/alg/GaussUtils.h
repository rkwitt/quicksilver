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

#ifndef __GAUSS_UTILS_H
#define __GAUSS_UTILS_H

#include <stddef.h>

namespace PyCA {

namespace GaussUtils {
    void GenerateKernel(float* h_filter, float sigma, size_t kRadius);
    void GenerateSupplementKernel(float* h_s, const float* h_filter, size_t kRadius);

    void ConvolutionX3D(float* h_o, const float* h_i, const float* h_k, const float* h_sK,
                        size_t kRadius, size_t sizeX, size_t sizeY, size_t sizeZ);
    void ConvolutionY3D(float* h_o, const float* h_i, const float* h_k, const float* h_sK,
                        size_t kRadius, size_t sizeX, size_t sizeY, size_t sizeZ);
    void ConvolutionZ3D(float* h_o, const float* h_i, const float* h_k, const float* h_sK,
                        size_t kRadius, size_t sizeX, size_t sizeY, size_t sizeZ);
};

} // end namespace PyCA

#endif

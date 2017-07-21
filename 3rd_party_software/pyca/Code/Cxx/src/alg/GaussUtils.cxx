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

#include <GaussUtils.h>
#include <math.h>
#include <iostream>

namespace PyCA {

namespace GaussUtils {
    void GenerateKernel(float* h_filter, float sigma, size_t kRadius){
        size_t kLength = kRadius * 2 + 1;
        float scale = 1.f / (2.f * sigma * sigma);

        // Compute the Gaussian coefficient 
        float sum  = h_filter[kRadius] = 1.f;
        for (size_t i=1; i <= kRadius; ++i){
            float fval = exp( -(float)(i * i) * scale);
            h_filter[kRadius - i] = h_filter[kRadius + i] = fval;
            sum += 2 * fval;
        }

        // Normalize the Gaussian coefficient
        sum = 1.f / sum;
        for (size_t i=0; i < kLength; ++i)
            h_filter[i] *= sum;
    }

    void GenerateSupplementKernel(float* h_s, const float* h_filter, size_t kRadius){
        float sum    = 0.f;
        h_s[0] = 1.f;
        // Calculate the coefficient factor for the edge
        for (size_t i=1; i<= kRadius; ++i){
            sum   += h_filter[i - 1];
            h_s[i] = 1.f / (1.f - sum);
        }
    }

    void ConvolutionX3D(float* h_o, const float* h_i,
                        const float* h_kernel, const float* h_sKernel,
                        size_t kRadius, size_t sizeX, size_t sizeY, size_t sizeZ){

        size_t planeSize = sizeX * sizeY;
        size_t id = 0;
        for (size_t iy=0; iy < sizeY; ++iy)
            for (size_t ix=0; ix < sizeX; ++ix, ++id) {
                size_t nId = id;
                for (size_t iz=0; iz < sizeZ; ++iz, nId+= planeSize ) {
                    float f = 1.f;
                    if (ix < kRadius)
                        f = h_sKernel[kRadius - ix];
                    else if (ix + kRadius > (sizeX-1))
                        f = h_sKernel[(ix + kRadius)-(sizeX - 1)];

                    float sum = 0;
                    for (int k=-kRadius; k <= (int)kRadius; ++k){
                        int x = ix + k;
                        if (x >=0 && x < (int)sizeX){
                            sum += h_i[nId + k] * h_kernel[kRadius - k];
                        }
                    }
                    h_o[nId] = sum * f;
                }
            }
    }

    void ConvolutionY3D(float* h_o, const float* h_i,
                        const float* h_kernel, const float* h_sKernel,
                        size_t kRadius, size_t sizeX, size_t sizeY, size_t sizeZ){

        size_t planeSize = sizeX * sizeY;
        size_t id = 0;
        for (size_t iy=0; iy < sizeY; ++iy)
            for (size_t ix=0; ix < sizeX; ++ix, ++id) {
                size_t nId = id;
                for (size_t iz=0; iz < sizeZ; ++iz, nId+= planeSize ) {
                    // compute the supplement elements
                    float f = 1.f;
                    if (iy < kRadius)
                        f = h_sKernel[kRadius - iy];
                    else if (iy + kRadius > (sizeY-1))
                        f = h_sKernel[(iy + kRadius) - (sizeY - 1)];
                    
                    float sum = 0;
                    for (int k=-kRadius; k <= (int)kRadius; ++k){
                        int y = iy + k;
                        if (y >=0 && y < (int)sizeY){
                            sum += h_i[nId + k * sizeX ] * h_kernel[kRadius - k];
                        }
                    }
                    h_o[nId] = sum * f;
                }
        
            }
    }

    void ConvolutionZ3D(float* h_o, const float* h_i,
                        const float* h_kernel, const float* h_sKernel,
                        size_t kRadius, size_t sizeX, size_t sizeY, size_t sizeZ){
        size_t planeSize = sizeX * sizeY;
        size_t id = 0;
        for (size_t iy=0; iy < sizeY; ++iy)
            for (size_t ix=0; ix < sizeX; ++ix, ++id) {
                size_t nId = id;
                for (size_t iz=0; iz < sizeZ; ++iz, nId+= planeSize ) {
                    float f = 1.f;
                    if (iz < kRadius)
                        f = h_sKernel[kRadius - iz];
                    else if (iz + kRadius > (sizeZ - 1))
                        f = h_sKernel[(iz + kRadius) - (sizeZ - 1)];

                    // 
                    float sum = 0;
                    for (int k=-kRadius; k <= (int)kRadius; ++k){
                        int z = iz + k;
                        if (z >=0 && z < (int)sizeZ){
                            sum += h_i[nId + k * sizeX * sizeY ] * h_kernel[kRadius - k];
                        }
                    }
                    h_o[nId] = sum * f;
                }
            }
    }
}
} // end namespace PyCA

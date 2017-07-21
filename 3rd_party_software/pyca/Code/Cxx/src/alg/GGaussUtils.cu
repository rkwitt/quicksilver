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

#include <GGaussUtils.h>
#include <pycaUtils.h>
#include <gcache.h>

// TEST make sure boost isn't included in nvcc code
#if defined(BOOST_COMPILER)
int bla[-1];
#endif

namespace PyCA {

namespace GGaussUtils {
#define MAX_KERNEL_LENGTH 1025
#define MAX_KERNEL_RADIUS (MAX_KERNEL_LENGTH / 2)
    // kernel coefficient function
    __device__ __constant__ float c_x_Kernel[MAX_KERNEL_LENGTH];
    __device__ __constant__ float c_y_Kernel[MAX_KERNEL_LENGTH];
    __device__ __constant__ float c_z_Kernel[MAX_KERNEL_LENGTH];

// supplement kernel for the boundary
    __device__ __constant__ float c_sx_Kernel[MAX_KERNEL_RADIUS+1];
    __device__ __constant__ float c_sy_Kernel[MAX_KERNEL_RADIUS+1];
    __device__ __constant__ float c_sz_Kernel[MAX_KERNEL_RADIUS+1];

    void setConvolutionKernelX(const float *h_Kernel, size_t kLength){
        cudaMemcpyToSymbol(c_x_Kernel, h_Kernel, kLength * sizeof(float));
    }

    void setConvolutionKernelY(const float *h_Kernel, size_t kLength){
        cudaMemcpyToSymbol(c_y_Kernel, h_Kernel, kLength * sizeof(float));
    }

    void setConvolutionKernelZ(const float *h_Kernel, size_t kLength){
        cudaMemcpyToSymbol(c_z_Kernel, h_Kernel, kLength * sizeof(float));
    }

    void setSupplementKernelX(const float *h_Kernel, size_t kRadiusPlusOne){
        cudaMemcpyToSymbol(c_sx_Kernel, h_Kernel, kRadiusPlusOne * sizeof(float));
    }

    void setSupplementKernelY(const float *h_Kernel, size_t kRadiusPlusOne){
        cudaMemcpyToSymbol(c_sy_Kernel, h_Kernel, kRadiusPlusOne * sizeof(float));
    }

    void setSupplementKernelZ(const float *h_Kernel, size_t kRadiusPlusOne){
        cudaMemcpyToSymbol(c_sz_Kernel, h_Kernel, kRadiusPlusOne * sizeof(float));
    }

    __global__ void ConvolutionX3D_kernel(float* d_o, const float* d_i,
                                          uint kRadius,
                                          uint sizeX, uint sizeY, uint sizeZ){

        const uint ix = blockDim.x * blockIdx.x + threadIdx.x;
        const uint iy = blockDim.y * blockIdx.y + threadIdx.y;
    
        if (ix >= sizeX || iy >= sizeY)
            return;
        
        uint id = iy * sizeX + ix;
        for (uint iz = 0; iz < sizeZ; ++iz, id+= sizeX * sizeY){
            // compute the supplement elements
            float f = 1.f;
            if (ix < kRadius)
                f = c_sx_Kernel[kRadius - ix];
            else if (ix + kRadius > (sizeX-1))
                f = c_sx_Kernel[(ix + kRadius)-(sizeX - 1)];

            float sum = 0;
            for (int k=-kRadius; k <= (int)kRadius; ++k){
                int x = ix + k;
                if (x >=0 && x < sizeX){
                    sum += d_i[id + k] * c_x_Kernel[kRadius - k];
                }
            }
            d_o[id] = sum * f;
        }
    }

    void ConvolutionX3D(float* d_o, const float* d_i,
                        int kRadius,
                        int sizeX, int sizeY, int sizeZ, StreamT stream){

        dim3 threads(16,16);
        dim3 grids(iDivUp(sizeX, threads.x),iDivUp(sizeY, threads.y));
        ConvolutionX3D_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, kRadius, sizeX, sizeY, sizeZ);
    }

//////////////////////////////////////////////////////////////////////////////////
// Convolution column
//////////////////////////////////////////////////////////////////////////////////
    __global__ void ConvolutionY3D_kernel(float* d_o, const float* d_i,
                                          int kRadius,
                                          int sizeX, int sizeY, int sizeZ){

        const int ix = blockDim.x * blockIdx.x + threadIdx.x;
        const int iy = blockDim.y * blockIdx.y + threadIdx.y;
    
        if (ix >= sizeX || iy >= sizeY)
            return;

        int id = iy * sizeX + ix;
        for (int iz = 0; iz < sizeZ; ++iz){
            // compute the supplement elements

            float f = 1.f;
            if (iy - kRadius < 0)
                f = c_sy_Kernel[-(iy-kRadius)];
            else if (iy + kRadius > (sizeY-1))
                f = c_sy_Kernel[(iy + kRadius) - (sizeY - 1)];

            // 
            float sum = 0;
            for (int k=-kRadius; k <= (int)kRadius; ++k){
                int y = iy + k;
                if (y >=0 && y < sizeY){
                    sum += d_i[id + k * sizeX ] * c_y_Kernel[kRadius - k];
                }
            }
            d_o[id] = sum * f;
            id     += sizeX * sizeY;
        }
    }

    void ConvolutionY3D(float* d_o, const float* d_i,
                        int kRadius,
                        int sizeX, int sizeY, int sizeZ, StreamT stream){

        dim3 threads(16,16);
        dim3 grids(iDivUp(sizeX, threads.x),iDivUp(sizeY, threads.y));
        ConvolutionY3D_kernel<<<grids, threads, 0, stream>>>
            (d_o, d_i, kRadius, sizeX, sizeY, sizeZ);
    }



    __global__ void ConvolutionZ3D_kernel(float* d_o, const float* d_i,
                                          int kRadius,
                                          int sizeX, int sizeY, int sizeZ){

        const int ix = blockDim.x * blockIdx.x + threadIdx.x;
        const int iy = blockDim.y * blockIdx.y + threadIdx.y;
    
        if (ix >= sizeX || iy >= sizeY)
            return;

        int id = iy * sizeX + ix;
    
        for (int iz = 0; iz < sizeZ; ++iz){
            // compute the supplement elements

            float f = 1.f;
            if (iz - kRadius < 0)
                f = c_sz_Kernel[-(iz-kRadius)];
            else if (iz + kRadius > (sizeZ - 1))
                f = c_sz_Kernel[(iz + kRadius) - (sizeZ - 1)];

            // 
            float sum = 0;
            for (int k=-kRadius; k <= (int)kRadius; ++k){
                int z = iz + k;
                if (z >=0 && z < sizeZ){
                    sum += d_i[id + k * sizeX * sizeY ] * c_z_Kernel[kRadius - k];
                }
            }
            d_o[id] = sum * f;
            id     += sizeX * sizeY;
        }
    }

    void ConvolutionZ3D(float* d_o, const float* d_i,
                        int kRadius,
                        int sizeX, int sizeY, int sizeZ, StreamT stream){
        dim3 threads(16,16);
        dim3 grids(iDivUp(sizeX, threads.x),iDivUp(sizeY, threads.y));
        ConvolutionZ3D_kernel<<<grids, threads, 0, stream>>>
            (d_o, d_i, kRadius, sizeX, sizeY, sizeZ);
    }


    __global__ void ConvolutionX3D_kernel2(float* d_o, const float* d_i,
                                           uint kRadius,
                                           uint sizeX, uint sizeY, uint sizeZ){

        const uint ix = blockDim.x * blockIdx.x + threadIdx.x;
        const uint iy = blockDim.y * blockIdx.y + threadIdx.y;
    
        if (ix >= sizeX || iy >= sizeY)
            return;
        
        uint id = iy * sizeX + ix;
        for (uint iz = 0; iz < sizeZ; ++iz, id+= sizeX * sizeY){
            // compute the supplement elements
            float f = 1.f;
            if (ix < kRadius)
                f = fetch_x(kRadius - ix, (float*)NULL);
            else if (ix + kRadius > (sizeX-1))
                f = fetch_x(ix + kRadius-(sizeX - 1), (float*)NULL);

            float sum = 0;
            for (int k=-kRadius; k <= (int)kRadius; ++k){
                int x = ix + k;
                if (x >=0 && x < sizeX){
                    sum += d_i[id + k] * fetch_y(kRadius - k, (float*)NULL);
                }
            }
            d_o[id] = sum * f;
        }
    }

    void ConvolutionX3D(float* d_o, const float* d_i,
                        const float* d_Kernel, const float* d_sKernel,
                        int kRadius,
                        int sizeX, int sizeY, int sizeZ, StreamT stream){

        dim3 threads(16,16);
        dim3 grids(iDivUp(sizeX, threads.x),iDivUp(sizeY, threads.y));
        cache_bind(d_sKernel);
        cache_bind_y(d_Kernel);
        ConvolutionX3D_kernel2<<<grids, threads, 0, stream>>>(
            d_o, d_i, kRadius, sizeX, sizeY, sizeZ);
    }
};
} // end namespace PyCA

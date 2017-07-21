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

#include <cstdio>

#include "GReduce.h"
#include <pycaUtils.h>
#include <mem.h>
#include "ReduceKernel.cu"
#include <gmem.h>
#include <Vec2D.h>
#include "CudaUtils.h"

// TEST make sure boost isn't included in nvcc code
#if defined(BOOST_COMPILER)
int bla[-1];
#endif

namespace PyCA {

template<typename T, typename opers>
T accumulate(T* h_temp, size_t N){
    T sum = opers::identity();
    for (size_t i=0; i< N; ++i)
        opers::iop(sum, h_temp[i]);
    return sum;
}

GReduce::GReduce() {
    uint size = MAX_NUMBER_REDUCE_BLOCKS * MAX_NUMBER_REDUCE_STREAMS;
    dmemAlloc(d_temp, size);
    phmemAlloc(h_temp, size);

    CudaUtils::CheckCUDAError(__FILE__,__LINE__,
			      "GReduce::GReduce");
}

GReduce::~GReduce() {
    if (d_temp)
        dmemFree(d_temp);
    
    if (h_temp)
        phmemFree(h_temp);
}

void getNumBlocksAndThreads(int n, int maxBlocks, int maxThreads, int &blocks, int &threads)
{
    threads = nextPowerOf2(iDivUp(n, 2));
    threads = (threads < maxThreads) ? threads : maxThreads;
    
    blocks = iDivUp(n, threads * 2);
    blocks = blocks < maxBlocks ? blocks : maxBlocks;
}


template <typename T, typename oper>
void reduce(uint n, int threads, int blocks, const T *d_i, T *d_o, StreamT stream){
    dim3 dimBlock(threads);
    dim3 dimGrid(blocks);
    switch (threads)  {
        case 512: reduce_kernel<T, oper, 512><<< dimGrid, dimBlock, 0, stream>>>(d_i, d_o, n);
            break;
        case 256: reduce_kernel<T, oper, 256><<< dimGrid, dimBlock, 0, stream>>>(d_i, d_o, n);
            break;
        case 128: reduce_kernel<T, oper, 128><<< dimGrid, dimBlock, 0, stream>>>(d_i, d_o, n);
            break;
        case 64:  reduce_kernel<T, oper, 64><<< dimGrid, dimBlock, 0, stream>>>(d_i, d_o, n);
            break;
        case 32:  reduce_kernel<T, oper, 32><<< dimGrid, dimBlock, 0, stream>>>(d_i, d_o, n);
            break;
        case 16:  reduce_kernel<T, oper, 16><<< dimGrid, dimBlock, 0, stream>>>(d_i, d_o, n);
            break;
        case  8:  reduce_kernel<T, oper, 8><<< dimGrid, dimBlock, 0, stream>>>(d_i, d_o, n);
            break;
        case  4:  reduce_kernel<T, oper, 4><<< dimGrid, dimBlock, 0, stream>>>(d_i, d_o, n);
            break;
        case  2:  reduce_kernel<T, oper, 2><<< dimGrid, dimBlock, 0, stream>>>(d_i, d_o, n);
            break;
        case  1:  reduce_kernel<T, oper, 1><<< dimGrid, dimBlock, 0, stream>>>(d_i, d_o, n);
            break;
    }
}


template <typename T, MATH_OPS op>
void GReduce::reduce(T* h_o, const T* d_i, size_t n, bool update, StreamT stream){
    int blocks, threads;

    typedef class MOpers<T, op> oper;
    
    getNumBlocksAndThreads(n, MAX_NUMBER_REDUCE_BLOCKS, MAX_NUMBER_REDUCE_THREADS, blocks, threads);
    PyCA::reduce<T, oper> (n, threads, blocks, d_i, (T*)d_temp, stream);

    cpyArrayD2H((T*)h_temp, (T*) d_temp, blocks);
    T s = accumulate<T, oper>((T*)h_temp, blocks);
    h_o[0] = update ? oper::op(h_o[0], s) : s;
}

template <typename T, typename oper, typename oper1>
inline void reduce(size_t n, int threads, int blocks, const T *d_i, T *d_o, StreamT stream){
    dim3 dimBlock(threads, 1, 1);
    dim3 dimGrid(blocks, 1, 1);
    int smem = threads * sizeof(T);
    switch (threads)
    {
        case 512:  reduce_kernel<T, oper, oper1, 512><<< dimGrid, dimBlock, smem, stream>>>(d_i, d_o, n);
            break;
        case 256:  reduce_kernel<T, oper, oper1, 256><<< dimGrid, dimBlock, smem, stream>>>(d_i, d_o, n);
            break;
        case 128:  reduce_kernel<T, oper, oper1, 128><<< dimGrid, dimBlock, smem, stream>>>(d_i, d_o, n);
            break;
        case 64:   reduce_kernel<T, oper, oper1, 64><<< dimGrid, dimBlock, smem, stream>>>(d_i, d_o, n);
            break;
        case 32:   reduce_kernel<T, oper, oper1, 32><<< dimGrid, dimBlock, smem, stream>>>(d_i, d_o, n);
            break;
        case 16:   reduce_kernel<T, oper, oper1, 16><<< dimGrid, dimBlock, smem, stream>>>(d_i, d_o, n);
            break;
        case  8:   reduce_kernel<T, oper, oper1, 8><<< dimGrid, dimBlock, smem, stream>>>(d_i, d_o, n);
            break;
        case  4:   reduce_kernel<T, oper, oper1, 4><<< dimGrid, dimBlock, smem, stream>>>(d_i, d_o, n);
            break;
        case  2:   reduce_kernel<T, oper, oper1, 2><<< dimGrid, dimBlock, smem, stream>>>(d_i, d_o, n);
            break;
        case  1:   reduce_kernel<T, oper, oper1, 1><<< dimGrid, dimBlock, smem, stream>>>(d_i, d_o, n);
            break;
    }
}

template <typename T, MATH_OPS op, MATH_OPS op1>
void GReduce::compReduce(T* h_o, const T* d_i, size_t n, bool update, StreamT stream){
    typedef class MOpers<T, op> oper;
    typedef class MOpers<T, op1> oper1;

    int blocks, threads;
    getNumBlocksAndThreads(n, MAX_NUMBER_REDUCE_BLOCKS, MAX_NUMBER_REDUCE_THREADS, blocks, threads);
    PyCA::reduce<T, oper, oper1>(n, threads, blocks, d_i, (T*) d_temp, stream);

    cudaMemcpy(h_temp, d_temp, sizeof(T) * blocks, cudaMemcpyDeviceToHost);
    T s = accumulate<T, oper>((T*)h_temp, blocks);
    h_o[0] = update ? oper::op(h_o[0], s) : s;
}

/**
 * @brief Perform cuda kernel parallel reductin 
 *        s = a1 + a2 + a3 + .... + an
 * @param[in]  T         Input data type (currently int, float)
 *             oper    Binary operation (+, max, min)
 *             oper1    Self data function (square, cude, sqrt, abs) 
 *             blockSize Size of block (related to optimize problem)
 *             d_i   Input data
 *             n         Size of the input
 * @param[out] array of output redution perform for each block
 *
*/


template <typename T, MATH_OPS op, MATH_OPS op1>
void GReduce::bireduce(T* h_o, const T* d_i, size_t n, bool update, StreamT stream)
{
    typedef class MOpers<T, op> oper;
    typedef class MOpers<T, op1> oper1;

    const uint blockSize = MAX_NUMBER_REDUCE_THREADS;

    dim3 threads(blockSize);
    uint nBlocks = fminf(iDivUp(n, 2 * blockSize),  MAX_NUMBER_REDUCE_BLOCKS);
    dim3 grids(nBlocks);

    T* d_rd0 = (T*) d_temp;
    T* d_rd1 = d_rd0 + MAX_NUMBER_REDUCE_BLOCKS;

    T* h_rd0 = (T*) h_temp;
    T* h_rd1 = h_rd0 + MAX_NUMBER_REDUCE_BLOCKS;
    biReduce_kernel<T, oper, oper1, blockSize><<<grids, threads, 0, stream>>>(d_rd0, d_rd1, d_i, n);
    
    cudaMemcpy(h_rd0, d_rd0, sizeof(T) * (nBlocks + MAX_NUMBER_REDUCE_BLOCKS), cudaMemcpyDeviceToHost);
    
    T rd0 = oper::identity();
    T rd1 = oper1::identity();

    for (uint i=0; i< nBlocks; ++i){
        oper::iop(rd0, h_rd0[i]);
        oper1::iop(rd1, h_rd1[i]);
    }

    h_o[0] = (update) ? oper::op(h_o[0], rd0) : rd0;
    h_o[1] = (update) ? oper1::op(h_o[1], rd1) : rd1;
}

template <typename T, MATH_OPS op, MATH_OPS op1>
void GReduce::product(T* h_o, const T*d_i, const T*d_i1, size_t n, bool update,StreamT stream)
{
    typedef class MOpers<T, op> oper;
    typedef class MOpers<T, op1> oper1;

    const uint blockSize = MAX_NUMBER_REDUCE_THREADS;
    dim3 threads(blockSize);
    uint nBlocks = fminf(iDivUp(n, 2 * blockSize),  MAX_NUMBER_REDUCE_BLOCKS);

    dim3 grids(nBlocks);
    reduceProduct_kernel<T, oper, oper1, blockSize><<<grids, threads, 0, stream>>>((T*)d_temp, d_i, d_i1, n);

    cudaMemcpy(h_temp, d_temp, sizeof(T) * nBlocks, cudaMemcpyDeviceToHost);
    
    T s = accumulate<T, oper>((T*)h_temp, nBlocks);
    h_o[0] = (update) ? oper::op(h_o[0], s) : s;
}

////////////////////////////////////////////////////////////////////////////////
// Instantiate for implementation
////////////////////////////////////////////////////////////////////////////////
template<typename T>
void GReduce::Max(T& h_o, const T* d_i, size_t n, bool update, StreamT stream){
    reduce<T, MATH_Max>((T*)&h_o, d_i, n, update, stream);
}

template void GReduce::Max(float& h_o, const float* d_i, size_t n, bool update, StreamT stream);
template void GReduce::Max(int& h_o, const int* d_i, size_t n, bool update, StreamT stream);
template void GReduce::Max(uint& h_o, const uint* d_i, size_t n, bool update, StreamT stream);

////////////////////////////////////////////////////////////////////////////////
template<typename T>
void GReduce::Min(T& h_o, const T* d_i, size_t n, bool update, StreamT stream){
    reduce<T, MATH_Min>((T*)&h_o, d_i, n, update, stream);
}

template void GReduce::Min(float& h_o, const float* d_i, size_t n, bool update, StreamT stream);
template void GReduce::Min(int& h_o, const int* d_i, size_t n, bool update, StreamT stream);
template void GReduce::Min(uint& h_o, const uint* d_i, size_t n, bool update, StreamT stream);

////////////////////////////////////////////////////////////////////////////////
template<typename T>
void GReduce::Sum(T& h_o, const T* d_i, size_t n, bool update, StreamT stream){
    reduce<T, MATH_Add>((T*)&h_o, d_i, n, update, stream);
}

template void GReduce::Sum(float& h_o, const float* d_i, size_t n, bool update, StreamT stream);
template void GReduce::Sum(int& h_o, const int* d_i, size_t n, bool update, StreamT stream);
template void GReduce::Sum(uint& h_o, const uint* d_i, size_t n, bool update, StreamT stream);


template<typename T>
void GReduce::LInf(T& h_o, const T* d_i, size_t n, bool update, StreamT stream){
    compReduce<T, MATH_Max, MATH_Abs>((T*)&h_o, d_i, n, update, stream);
};

template void GReduce::LInf(float& h_o, const float* d_i, size_t n, bool update, StreamT stream);
template void GReduce::LInf(int& h_o, const int* d_i, size_t n, bool update, StreamT stream);
template void GReduce::LInf(uint& h_o, const uint* d_i, size_t n, bool update, StreamT stream);

template<typename T>
void GReduce::L1(T& h_o, const T* d_i, size_t n, bool update, StreamT stream){
    compReduce<T, MATH_Add, MATH_Abs>((T*)&h_o, d_i, n, update, stream);
};

template void GReduce::L1(float& h_o, const float* d_i, size_t n, bool update, StreamT stream);
template void GReduce::L1(int& h_o, const int* d_i, size_t n, bool update, StreamT stream);
template void GReduce::L1(uint& h_o, const uint* d_i, size_t n, bool update, StreamT stream);

template<typename T>
void GReduce::Sum2(T& h_o, const T* d_i, size_t n, bool update, StreamT stream){
    compReduce<T, MATH_Add, MATH_Sqr>((T*)&h_o, d_i, n, update, stream);
};

template void GReduce::Sum2(float& h_o, const float* d_i, size_t n, bool update, StreamT stream);
template void GReduce::Sum2(int& h_o, const int* d_i, size_t n, bool update, StreamT stream);
template void GReduce::Sum2(uint& h_o, const uint* d_i, size_t n, bool update, StreamT stream);

template<typename T>
void GReduce::MaxMin(Vec2D<T>&h_o, const T* d_i, size_t n, bool update, StreamT stream){
    bireduce<T, MATH_Max, MATH_Min>((T*)&h_o.x, d_i, n, update, stream);
}

template void GReduce::MaxMin(Vec2D<float>& h_o, const float* d_i, size_t n, bool update, StreamT stream);
template void GReduce::MaxMin(Vec2D<int>& h_o, const int* d_i, size_t n, bool update, StreamT stream);
template void GReduce::MaxMin(Vec2D<uint>& h_o, const uint* d_i, size_t n, bool update, StreamT stream);

template<typename T>
void GReduce::Dot(T& h_o, const T* d_i, const T* d_i1, size_t n, bool update, StreamT stream){
    product<T, MATH_Add, MATH_Mul>((T*)&h_o, d_i, d_i1, n, update, stream);
}

template
void GReduce::Dot(float& h_o, const float* d_i, const float* d_i1, size_t n, bool update, StreamT stream);
template
void GReduce::Dot(int& h_o, const int* d_i, const int* d_i1, size_t n, bool update, StreamT stream);
template
void GReduce::Dot(uint& h_o, const uint* d_i, const uint* d_i1, size_t n, bool update, StreamT stream);

bool GReduce::selfTest(size_t n){
    int test = true;
    int* h_i  = new int [n];
    int* h_i1 = new int [n];

    for (size_t j=0; j< n; ++j) h_i[j] = (rand() & 0xff);
    for (size_t j=0; j< n; ++j) h_i1[j] = (rand() & 0xff);
    
    int *d_i;
    dmemAlloc(d_i, n);
    cpyArrayH2D(d_i, h_i, n);
    
    int *d_i1;
    dmemAlloc(d_i1,n);
    cpyArrayH2D(d_i1,h_i1,n);

    int *d_o;
    dmemAlloc(d_o, 256);
    
    int h_max    = -INT_MAX, h_min = INT_MAX;
    int h_LInf   = 0;
    int h_sum2   = 0;
    int h_sum    = 0;
    int h_L1     = 0;
    int h_dot    = 0;
    
    for (size_t i=0; i< n; ++i)
    {
        h_max  = fmaxf(h_max, h_i[i]);
        h_min  = fminf(h_min, h_i[i]);
        h_sum += h_i[i];

        h_LInf = fmaxf(h_LInf, h_i[i]);
        h_L1  += fabsf(h_i[i]);
        h_sum2+= h_i[i]*h_i[i];
        
        h_dot += h_i1[i] * h_i[i];
    }

    int d_max = -INT_MAX, d_min = INT_MAX;
    int d_LInf= 0;
    int d_L1= 0;
    int d_sum2= 0;
    int d_sum = 0;
    int d_dot = 0;
    Vec2D<int> d_pair;
    
    this->Sum(d_sum, d_i, n);
    this->Max(d_max, d_i, n);
    this->Min(d_min, d_i, n);
    this->LInf(d_LInf,d_i, n);
    this->L1(d_L1, d_i, n);
    this->Sum2(d_sum2,d_i, n);
    this->Dot(d_dot,d_i,d_i1, n);
    this->MaxMin(d_pair, d_i, n);

    fprintf(stderr, "Maximum value from CPU %d from GPU %d\n",h_max, d_max);
    fprintf(stderr, "Minumum value from CPU %d from GPU %d\n",h_min, d_min);
    fprintf(stderr, "Total value from CPU %d from GPU %d\n",h_sum, d_sum);
    
    fprintf(stderr, "Maximum abosulte value from CPU %d from GPU %d\n",h_LInf, d_LInf);
    fprintf(stderr, "Total square value from CPU %d from GPU %d\n",h_sum2, d_sum2);
    fprintf(stderr, "Dot product value from CPU %d from GPU %d\n",h_dot, d_dot);
    fprintf(stderr, "Max Min value from CPU %d %d from GPU %d %d\n",h_max, h_min, d_pair.x, d_pair.y);
    
    //Extensive test
    h_max = -INT_MAX, h_min = INT_MAX;
    h_LInf = 0;
    h_sum2 = 0;
    h_sum = 0;
    h_dot = 0;
    h_L1 = 0;
    
    for (int l=1; l < 10001;++l){
        h_max = fmaxf(h_max, h_i[l-1]);
        h_LInf = fmaxf(h_LInf, h_i[l-1]);
        h_min = fminf(h_min, h_i[l-1]);
        h_sum += h_i[l-1];
        h_sum2 += h_i[l-1]*h_i[l-1];
        h_dot += h_i1[l-1] * h_i[l-1];
        h_L1  += fabsf(h_i[l-1]);

        this->Sum(d_sum, d_i, l);
        this->Max(d_max, d_i, l);
        this->Min(d_min, d_i, l);
        this->LInf(d_LInf,d_i, l);
        this->L1(d_L1, d_i, l);
        this->Sum2(d_sum2,d_i, l);
        this->Dot(d_dot,d_i,d_i1, l);
        this->MaxMin(d_pair, d_i, l);
                
        if (d_max != h_max){
            fprintf(stderr, "Max Test FAILED at %d GPU %d CPU %d\n", l, d_max, h_max );
            test = false;
        }

        if (d_min != h_min){
            fprintf(stderr, "Min Test FAILED at %d GPU %d CPU %d\n", l, d_min, h_min );
            test = false;
        }

        if (d_LInf!= h_LInf){
            fprintf(stderr, "MaxAbs Test FAILED at %d GPU %d CPU %d\n", l, d_LInf, h_LInf );
            test = false;
        }

        if (d_sum!= h_sum){
            fprintf(stderr, "Sum Test FAILED at %d GPU %d CPU %d\n", l, d_sum, h_sum );
            test = false;
        }

        if (d_sum2!= h_sum2){
            fprintf(stderr, "Sum SQR Test FAILED at %d GPU %d CPU %d\n", l, d_sum2, h_sum2 );
            test = false;
        }

        if (d_dot!= h_dot){
            fprintf(stderr, "Dot Test FAILED at %d GPU %d CPU %d\n", l, d_dot, h_dot );
            test = false;
        }
        
        if ( d_pair.x != h_max || d_pair.y != h_min){
            fprintf(stderr, "Max Min Test FAILED at %d GPU %d %d CPU %d %d\n",
                    l, d_pair.x, d_pair.y, h_max, h_min);
            test = false;
        }
        if (test == false)
            break;
    }

    if (test)
        fprintf(stderr, "Test PASSED  \n");
    
    delete []h_i1;
    delete []h_i;
    cudaFree(d_i);
    cudaFree(d_i1);
    return test;
}

}

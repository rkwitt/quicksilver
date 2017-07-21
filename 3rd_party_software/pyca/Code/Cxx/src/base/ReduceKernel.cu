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

#define MAX_NUMBER_REDUCE_STREAMS     4
#define MAX_NUMBER_REDUCE_THREADS     128
#define MAX_NUMBER_REDUCE_BLOCKS      128

namespace PyCA {

template<typename T, typename traits, uint blockSize>
__inline__  __device__ void reduce_shared(T& mySum, volatile T* sdata, uint tid)
{
    // do reduction in shared mem
    if (blockSize >= 512) {
        if (tid < 256) { sdata[tid] = mySum = traits::op(mySum, sdata[tid + 256]); } __syncthreads();
    }
    if (blockSize >= 256) {
        if (tid < 128) { sdata[tid] = mySum = traits::op(mySum, sdata[tid + 128]); } __syncthreads();
    }
    if (blockSize >= 128) {
        if (tid <  64) { sdata[tid] = mySum = traits::op(mySum, sdata[tid +  64]); } __syncthreads();
    }
                            
    if (tid < 32)
    {
        if (blockSize >=  64) sdata[tid] = mySum = traits::op(mySum, sdata[tid + 32]);
        if (blockSize >=  32) sdata[tid] = mySum = traits::op(mySum, sdata[tid + 16]);
        if (blockSize >=  16) sdata[tid] = mySum = traits::op(mySum, sdata[tid +  8]);
        if (blockSize >=   8) sdata[tid] = mySum = traits::op(mySum, sdata[tid +  4]);
        if (blockSize >=   4) sdata[tid] = mySum = traits::op(mySum, sdata[tid +  2]);
        if (blockSize >=   2) sdata[tid] = mySum = traits::op(mySum, sdata[tid +  1]);
    }
}

/**
 * @brief Perform cuda kernel parallel reductin 
 *        s = a1 + a2 + a3 + .... + an
 * @param[in]  T         Input data type (currently int, float)
 *             traits    Binary operation (+, max, min)
 *             blockSize Size of block (related to optimize problem)
 *             d_i   Input data
 *             n         Size of the input
 * @param[out] array of output redution perform for each block
 *
*/
template <typename T, typename traits, uint blockSize>
__global__ void reduce_kernel(const T *d_i, T *d_o, uint n)
{
    volatile __shared__ T sdata[MAX_NUMBER_REDUCE_THREADS];
    // reading from global memory, writing to shared memory
    uint tid      = threadIdx.x;
    uint i        = blockIdx.x*(blockSize*2) + tid;
    uint gridSize = blockSize * 2 * gridDim.x;
    
    T mySum = traits::identity();;
    while (i + blockSize < n )
    {
        traits::iop(mySum,traits::op(d_i[i],d_i[i+blockSize]));
        i += gridSize;
    }
    if ( i < n) traits::iop(mySum, d_i[i]);
    sdata[tid] = mySum;

    __syncthreads();
    reduce_shared<T, traits, blockSize>(mySum, sdata, tid);

    // write result for this block to global mem 
    if (tid == 0) {
        d_o[blockIdx.x] = mySum;
    }
}


/**
 * @brief Perform cuda kernel parallel reductin 
 *        s = a1 + a2 + a3 + .... + an
 * @param[in]  T         Input data type (currently int, float)
 *             traits    Binary operation (+, max, min)
 *             traits1    Self data function (square, cude, sqrt, abs) 
 *             blockSize Size of block (related to optimize problem)
 *             d_i   Input data
 *             n         Size of the input
 * @param[out] array of output redution perform for each block
 *
*/
// Correct version of reduce function
template <typename T, typename traits, typename traits1,  uint blockSize>
__global__ void reduce_kernel(const T *d_i, T *d_o, uint n)
{
    volatile __shared__ T sdata[MAX_NUMBER_REDUCE_THREADS];
    
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    uint tid = threadIdx.x;
    uint i   = blockIdx.x*(blockSize*2) + tid;
    uint gridSize = blockSize*2*gridDim.x;

    T mySum = traits::identity();
    // we reduce multiple elements per thread.  The number is determined by the 
    // number of active thread blocks (via gridSize).  More blocks will result
    // in a larger gridSize and therefore fewer elements per thread
    while (i + blockSize < n ) {
        traits::iop(mySum,traits::op(traits1::op(d_i[i]),traits1::op(d_i[i+blockSize])));
        i += gridSize;
    }

    if ( i < n)
        traits::iop(mySum, traits1::op(d_i[i]));
    sdata[tid] = mySum;
    
    __syncthreads();

    reduce_shared<T, traits, blockSize>(mySum, sdata, tid);
    // write result for this block to global mem 
    if (tid == 0)
        d_o[blockIdx.x] = mySum;
}


////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename oper, typename oper1, uint blockSize>
__global__ void biReduce_kernel(T *d_o, T *d_o1, const T *d_i, uint n)
{
    volatile __shared__ T s0[MAX_NUMBER_REDUCE_THREADS];
    volatile __shared__ T s1[MAX_NUMBER_REDUCE_THREADS];
        
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    uint tid      = threadIdx.x;
    uint i        = blockIdx.x*(blockSize*2) + tid;
    uint gridSize = blockSize * 2 * gridDim.x;

    T  mySum  = oper::identity();
    T  mySum1 = oper1::identity();

    while (i + blockSize < n )
    {
        oper::iop(mySum, oper::op(d_i[i],d_i[i+blockSize]));
        oper1::iop(mySum1 ,oper1::op(d_i[i],d_i[i+blockSize]));
        i += gridSize;
    }

    if ( i < n){
        oper::iop(mySum, d_i[i]);
        oper1::iop(mySum1, d_i[i]);
    }

    s0[tid] = mySum;
    s1[tid] = mySum1;
    
    __syncthreads();

    // do reduction in shared mem
    if (blockSize >= 512) {
        if (tid < 256) {
            s0[tid] = mySum = oper::op(mySum, s0[tid + 256]);
            s1[tid] = mySum1 = oper1::op(mySum1, s1[tid + 256]);
        }
        __syncthreads();
    }

    if (blockSize >= 256) {
        if (tid < 128) {
            s0[tid] = mySum = oper::op(mySum, s0[tid + 128]);
            s1[tid] = mySum1 = oper1::op(mySum1, s1[tid + 128]);
        }
        __syncthreads();
    }

    if (blockSize >= 128) {
        if (tid < 64) {
            s0[tid] = mySum = oper::op(mySum, s0[tid + 64]);
            s1[tid] = mySum1 = oper1::op(mySum1, s1[tid + 64]);
        }
        __syncthreads();
    }

    if (tid < 32)
    {
        if (blockSize >=  64) {
            s0[tid] = mySum = oper::op(mySum, s0[tid + 32]);
            s1[tid] = mySum1 = oper1::op(mySum1, s1[tid + 32]);
        }

        if (blockSize >=  32) {
            s0[tid] = mySum = oper::op(mySum, s0[tid + 16]);
            s1[tid] = mySum1 = oper1::op(mySum1, s1[tid + 16]);
        }

        if (blockSize >=  16) {
            s0[tid] = mySum = oper::op(mySum, s0[tid + 8]);
            s1[tid] = mySum1 = oper1::op(mySum1, s1[tid + 8]);
        }

        if (blockSize >=  8) {
            s0[tid] = mySum = oper::op(mySum, s0[tid + 4]);
            s1[tid] = mySum1 = oper1::op(mySum1, s1[tid + 4]);
        }

        if (blockSize >=  4) {
            s0[tid] = mySum = oper::op(mySum, s0[tid + 2]);
            s1[tid] = mySum1 = oper1::op(mySum1, s1[tid + 2]);
        }

        if (blockSize >=  2) {
            s0[tid] = mySum = oper::op(mySum, s0[tid + 1]);
            s1[tid] = mySum1 = oper1::op(mySum1, s1[tid + 1]);
        }
    }
    // write result for this block to global mem 
    if (tid == 0){
        d_o[blockIdx.x] = mySum;
        d_o1[blockIdx.x] = mySum1;
    }
}

////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename opers, typename opers1, uint blockSize>
__global__ void reduceProduct_kernel(T *d_o, const T*d_i, const T*d_i1, size_t n)
{
    volatile __shared__ T sdata[MAX_NUMBER_REDUCE_THREADS];
        
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    uint blockId  = blockIdx.x;
    uint tid      = threadIdx.x;
    uint i        = blockId * (blockSize * 2) + tid;
    uint gridSize = (blockSize * 2) * gridDim.x;

    T mySum = opers::identity();
    
    // we reduce multiple elements per thread.  The number is determined by the 
    // number of active thread blocks (via gridSize).  More blocks will result
    // in a larger gridSize and therefore fewer elements per thread
    while (i + blockSize < n ) {
        T t1 = opers1::op(d_i[i], d_i1[i]);
        T t2 = opers1::op(d_i[i + blockSize], d_i1[i + blockSize]);
        opers::iop(mySum,opers::op(t1, t2));
        i += gridSize;
    }

    if ( i < n) {
        T t1 = opers1::op(d_i[i], d_i1[i]);
        opers::iop(mySum,t1);
    }

    sdata[tid] = mySum;
    
    __syncthreads();

    // do reduction in shared mem
    reduce_shared<T, opers, blockSize>(mySum, sdata, tid);
    
    // write result for this block to global mem 
    if (tid == 0)
        d_o[blockIdx.x] = sdata[0];
}
} // end namespace PyCA

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

#ifndef __REDUCE_STREAM_KERNEL_CU
#define __REDUCE_STREAM_KERNEL_CU

#if __CUDA_ARCH__ != 200
#define NUMBER_OF_CORES               (240 * 16)
#else
#define NUMBER_OF_CORES               (240)
#endif

#define OPTIMAL_THREAD_SIZE           64
#define MAX_NUMBER_OF_REDUCE_STREAM   8
#define M    NUMBER_OF_CORES          // Number of core 
#define N    OPTIMAL_THREAD_SIZE      // Best number of thread per block 
#define REDUCE_STREAM_SIZE            (M * N)

namespace PyCA {

template<typename T, class op>
__global__ void reduce_L1_kernel(T *res, const T* d_i, int size) {
    __shared__ T shm[N];
    int idx=blockIdx.x*blockDim.x+threadIdx.x;

    T s=op::identity();
    int j=idx;
    for(; j + M * N <size; j+=M*N*2)
        op::iop(s, op::op(d_i[j], d_i[j + M*N]));
    shm[threadIdx.x]= (j < size) ? op::op(s, d_i[j]) : s;
    __syncthreads();
    if(threadIdx.x==0) {
        T s=op::identity();
        for(int i=0; i<N; i++)
            op::iop(s, shm[i]);
        res[blockIdx.x]=s;
    }
}

template<typename T, class op, class op1>
__global__ void compReduce_L1_kernel(T *res, const T* d_i, int size) {
    __shared__ T shm[N];
    int idx=blockIdx.x*blockDim.x+threadIdx.x;

    T s=op::identity();
    int j=idx;
    for(; j + M * N <size; j+=M*N*2)
        op::iop(s, op::op(op1::op(d_i[j]), op1::op(d_i[j + M*N])));
    shm[threadIdx.x]= (j < size) ? op::op(s, op1::op(d_i[j])) : s;
    __syncthreads();
    if(threadIdx.x==0) {
        T s=op::identity();
        for(int i=0; i<N; i++)
            op::iop(s, shm[i]);
        res[blockIdx.x]=s;
    }
}

template<typename T, class op, class op1>
__global__ void product_L1_kernel(T *res, const T* d_i0, const T* d_i1, int size) {
    __shared__ T shm[N];
    int idx=blockIdx.x*blockDim.x+threadIdx.x;
    T s=op::identity();
    int j=idx;
    for(; j + M * N <size; j+=M*N*2)
        op::iop(s, op::op(op1::op(d_i0[j],d_i1[j]), op1::op(d_i0[j + M*N], d_i1[j + M*N])));
    shm[threadIdx.x]= (j < size) ? op::op(s, op1::op(d_i0[j],d_i1[j])) : s;
    __syncthreads();
    if(threadIdx.x==0) {
        T s=op::identity();
        for(int i=0; i<N; i++)
            op::iop(s, shm[i]);
        res[blockIdx.x]=s;
    }
}


////////////////////////////////////////////////////////////////////////////////
template<typename T, class op, bool accumulate>
__device__ void reduce_L2_dev(T& d_o, const T* d_i) {
    __shared__ T shm[N];
    uint idx=threadIdx.x;
    T s=op::identity();
#if M%(2*N)==0
    for(uint j=idx; j<M; j+=N*2)
        op::iop(s, op::op(d_i[j], d_i[j + N]));
#else
    for(uint j=idx; j<M; j+=N)
        op::iop(s, d_i[j]);
#endif
    shm[threadIdx.x]=s;
    __syncthreads();
    if(idx==0) {
        T s=op::identity();
        for(uint i=0; i<N; i++)
            op::iop(s, shm[i]);

        d_o=(accumulate) ? op::op(d_o, s) : s;
    }
}

template<typename T, class op>
__global__ void reduce_L2_kernel(T *d_i) {
    reduce_L2_dev<T, op, false>(d_i[0], d_i);
}

template<typename T, class op, bool accumulate>
__global__ void reduce_L2_kernel(T* d_o, const T* d_i) {
    reduce_L2_dev<T, op, accumulate>(d_o[0], d_i);
}

////////////////////////////////////////////////////////////////////////////////
template<typename T, class op, class op1>
__device__ void bireduce_L1_dev(T& res, T& res1, const T* d_i, int size) {
    __shared__ T shm[N];
    __shared__ T shm1[N];
    
    int idx=blockIdx.x*blockDim.x+threadIdx.x;

    T s =op::identity();
    T s1=op1::identity();
    
    int j=idx;
    for(; j + M * N <size; j+=M*N*2) {
        op::iop(s, op::op(d_i[j], d_i[j + M*N]));
        op1::iop(s1, op1::op(d_i[j], d_i[j + M*N]));
    }

    shm[threadIdx.x]  = (j < size) ? op::op(s, d_i[j]) : s;
    shm1[threadIdx.x] = (j < size) ? op1::op(s1, d_i[j]) : s1;
    
    __syncthreads();
#if 1
    if(threadIdx.x==0) {
        T s =op::identity();
        T s1=op1::identity();
        
        for(int i=0; i<N; i++){
            op::iop(s, shm[i]);
            op1::iop(s1, shm1[i]);
        }
        res  = s;
        res1 = s1;
    }
#else
    if(threadIdx.x==0) {
        T s =op::identity();
        for(int i=0; i<N; i++)
            op::iop(s, shm[i]);
        res = s;
    }
    if(threadIdx.x==1) {
        T s1 =op1::identity();
        for(int i=0; i<N; i++)
            op1::iop(s1, shm1[i]);
        res1 = s1;
    }
#endif
}

template<typename T, class op, class op1>
__global__ void bireduce_L1_kernel(T *res, T* res1, const T* d_i,  int size) {
    int blockId = blockIdx.x;
    bireduce_L1_dev<T, op, op1>
        (res[blockId], res1[blockId], d_i, size);
}

template<typename T, class  op, class op1, bool accumulate>
__device__ void reduce_ip2_op2_L2_dev(T& d_o, T& d_o1, const T* d_i, const T* d_i1) {
    __shared__ T shm[N];
    __shared__ T shm1[N];
    uint idx=threadIdx.x;
    T s =op::identity();
    T s1=op1::identity();
#if M%(2*N)==0
    for(uint j=idx; j<M; j+=N*2){
        op::iop(s, op::op(d_i[j], d_i[j + N]));
        op1::iop(s1, op1::op(d_i1[j], d_i1[j + N]));
    }
#else
    for(uint j=idx; j<M; j+=N){
        op::iop(s, d_i[j]);
        op1::iop(s1, d_i1[j]);
    }
#endif
    shm[threadIdx.x]=s;
    shm1[threadIdx.x]=s1;
    __syncthreads();
    if(idx==0) {
        T s=op::identity();
        T s1=op1::identity();
        for(uint i=0; i<N; i++){
            op::iop(s, shm[i]);
            op1::iop(s1, shm1[i]);
        }
        d_o = accumulate ? op::op(d_o, s) : s;
        d_o1 = accumulate ? op1::op(d_o1, s1) : s1;
    }
}

template<typename T, class  op, class op1, bool accumulate>
__global__ void reduce_ip2_op2_L2_kernel(T* d_o, T* d_o1, const T* d_i, const T* d_i1) {
    reduce_ip2_op2_L2_dev<T, op, op1, accumulate>(d_o[0], d_o1[0], d_i, d_i1);
}

template<typename T, class  op, class op1, bool accumulate>
__global__ void reduce_ip2_op2_L2_kernel(T* d_o, T *d_i, const T* d_i1) {
    reduce_ip2_op2_L2_dev<T, op, op1, accumulate>(d_o[0], d_o[1], d_i, d_i1);
}

////////////////////////////////////////////////////////////////////////////////
template<typename T, class op, class op1, class op2>
__device__ void trireduce_L1_dev(T& res, T& res1, T& res2, const T* d_i, int size) {
    __shared__ T shm[N];
    __shared__ T shm1[N];
    __shared__ T shm2[N];
    
    int idx=blockIdx.x*blockDim.x+threadIdx.x;

    T s =op::identity();
    T s1=op1::identity();
    T s2=op2::identity();
    
    int j=idx;
    for(; j + M * N <size; j+=M*N*2) {
        op::iop(s, op::op(d_i[j], d_i[j + M*N]));
        op1::iop(s1, op1::op(d_i[j], d_i[j + M*N]));
        op2::iop(s2, op2::op(d_i[j], d_i[j + M*N]));
    }

    shm[threadIdx.x]  = (j < size) ? op::op(s, d_i[j]) : s;
    shm1[threadIdx.x] = (j < size) ? op1::op(s1, d_i[j]) : s1;
    shm2[threadIdx.x] = (j < size) ? op2::op(s2, d_i[j]) : s2;
    
    __syncthreads();
#if 1
    if(threadIdx.x==0) {
        T s =op::identity();
        T s1=op1::identity();
        T s2=op2::identity();
        
        for(int i=0; i<N; i++){
            op::iop(s, shm[i]);
            op1::iop(s1, shm1[i]);
            op2::iop(s2, shm2[i]);
        }
        res  = s;
        res1 = s1;
        res2 = s2;
    }
#else
    if(threadIdx.x==0) {
        T s =op::identity();
        for(int i=0; i<N; i++)
            op::iop(s, shm[i]);
        res = s;
    }
    if(threadIdx.x==1) {
        T s1 =op1::identity();
        for(int i=0; i<N; i++)
            op1::iop(s1, shm1[i]);
        res1 = s1;
    }
    if(threadIdx.x==2) {
        T s2 =op2::identity();
        for(int i=0; i<N; i++)
            op2::iop(s2, shm2[i]);
        res2 = s2;
    }
#endif
}

template<typename T, class op, class op1, class op2, bool accumulate>
__global__ void trireduce_L1_kernel(T *res, T* res1, T* res2,
                                    const T* d_i,  int size) {
    int blockId = blockIdx.x;
    trireduce_L1_dev<T, op, op1, op2, accumulate>
        (res[blockId], res1[blockId], res2[blockId], d_i, size);
}

template<typename T, class  op, class op1, class op2, bool accumulate>
__device__ void reduce_ip3_op3_L2_dev(T& d_o, T& d_o1, T& d_o2,
                                      const T* d_i, const T* d_i1, const T* d_i2) {
    __shared__ T shm[N];
    __shared__ T shm1[N];
    __shared__ T shm2[N];
    uint idx=threadIdx.x;
    
    T s =op::identity();
    T s1=op1::identity();
    T s2=op2::identity();
#if M%(2*N)==0
    for(uint j=idx; j<M; j+=N*2){
        op::iop(s, op::op(d_i[j], d_i[j + N]));
        op1::iop(s1, op1::op(d_i1[j], d_i1[j + N]));
        op2::iop(s2, op2::op(d_i2[j], d_i2[j + N]));
    }
#else
    for(uint j=idx; j<M; j+=N){
        op::iop(s, d_i[j]);
        op1::iop(s1, d_i1[j]);
        op2::iop(s2, d_i2[j]);
    }
#endif
    shm[threadIdx.x]=s;
    shm1[threadIdx.x]=s1;
    shm2[threadIdx.x]=s2;
    __syncthreads();
    if(idx==0) {
        T s=op::identity();
        T s1=op1::identity();
        T s2=op2::identity();
        
        for(uint i=0; i<N; i++){
            op::iop(s, shm[i]);
            op1::iop(s1, shm1[i]);
            op2::iop(s2, shm2[i]);
        }
        d_o = accumulate ? op::op(d_o, s) : s;
        d_o1 = accumulate ? op1::op(d_o1, s1) : s1;
        d_o2 = accumulate ? op2::op(d_o2, s2) : s2;
    }
}

template<typename T, class  op, class op1, bool accumulate>
__global__ void reduce_ip3_op3_L3_kernel(T* d_o, T* d_o1, T* d_o2,
                                         const T* d_i, const T* d_i1, const T* d_i2) {
    reduce_ip3_op3_L2_dev<T, op, op1, accumulate>(d_o[0], d_o1[0], d_o2[0],
                                                  d_i, d_i1, d_i2);
}

template<typename T, class  op, class op1, bool accumulate>
__global__ void reduce_ip3_op3_L3_kernel(T* d_o, const T* d_i, const T* d_i1, const T* d_i2) {
    reduce_ip3_op3_L2_dev<T, op, op1, accumulate>(d_o[0], d_o[1], d_o[2],
                                                  d_i, d_i1, d_i2);
}

#endif
} // end namespace PyCA

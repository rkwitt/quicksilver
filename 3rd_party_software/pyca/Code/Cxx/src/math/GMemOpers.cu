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

#include <GMemOpers.h>

#include <cuda_runtime.h>
#include <pycaUtils.h>
#include <gcache.h>
#include "PyCAException.h"
#include "mem.h"

#include "MemOperDefs.h"

// TEST make sure boost isn't included in nvcc code
#if defined(BOOST_COMPILER)
int bla[-1];
#endif

#define GMEM_UNARY_OPERS(OP) MEM_UNARY_OPERS(GMemOpers, OP)
#define GMEM_BINARY_OPERS(OP) MEM_BINARY_OPERS(GMemOpers, OP)

#define GMEM_BINARY_OPERC(OP) MEM_BINARY_OPERC(GMemOpers, OP)
#define GMEM_BINARY_OPERC_I(OP) MEM_BINARY_OPERC_I(GMemOpers, OP)
#define GMEM_BINARY_OPER_NOC(OP) MEM_BINARY_OPER_NOC(GMemOpers, OP)
#define GMEM_BINARY_OPER_NOC_I(OP) MEM_BINARY_OPER_NOC_I(GMemOpers, OP)
#define GMEM_BINARY_OPER(OP) MEM_BINARY_OPER(GMemOpers, OP)
#define GMEM_BINARY_OPER_I(OP) MEM_BINARY_OPER_I(GMemOpers, OP)

namespace PyCA {

template<typename T>
void GMemOpers<T>::Copy(T* d_o, const T* d_i, size_t n, StreamT stream){
   acpyArrayD2D(d_o, d_i, n, stream);
}

// masked version requires kernel
template<typename T>
__global__ void Copy_masked_kernel(T* d_o, const T* d_i, const T* d_mask, 
				   uint n)
{
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
	if(d_mask[id])
	    d_o[id] = d_i[id];
}

template<typename T>
void GMemOpers<T>::Copy(T* d_o, const T* d_i, const T* d_mask, 
			size_t n, StreamT stream)
{
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids = make_grid(iDivUp(n, threads.x));
    Copy_masked_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_mask, n);
}

/////////////////////////////////////////////////////////////////////////////////
//  the function provided by CUDA is slow and is not flexible enough to set
//  value with different type
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void SetMem_kernel(T* d_o, T c, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = c;
}

template<class T>
__global__ void SetMem_kernel(T* d_o, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = fetch(0, (T*)NULL);
}

template<typename T>
void GMemOpers<T>::SetMem(T* d_o, const T& c, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids = make_grid(iDivUp(n, threads.x));

    if (onDev) {
        cache_bind(&c);
        SetMem_kernel<<<grids, threads, 0, stream>>>(d_o, n);
    } else {
        SetMem_kernel<<<grids, threads, 0, stream>>>(d_o, c, n);
    }


}

template<typename T>
__global__ void SetMem_masked_kernel(T* d_o, T c, const T* d_mask, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
	if(d_mask[id])
	    d_o[id] = c;
}

template<class T>
__global__ void SetMem_masked_kernel(T* d_o, const T* d_mask, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
	if(d_mask[id])
	    d_o[id] = fetch(0, (T*)NULL);
}

template<typename T>
void GMemOpers<T>::SetMem(T* d_o, const T& c, const T* d_mask, size_t n, 
			  StreamT stream, bool onDev)
{
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids = make_grid(iDivUp(n, threads.x));

    if (onDev) {
        cache_bind(&c);
        SetMem_masked_kernel<<<grids, threads, 0, stream>>>(d_o, d_mask, n);
    } else {
        SetMem_masked_kernel<<<grids, threads, 0, stream>>>(d_o, c, d_mask, n);
    }


}


// /////////////////////////////////////////////////////////////////////////////////
// //  SetLinear
// //  Initiate the value of an unsigned array with a linear ram
// //////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void SetLinear_kernel(T* d_o, uint n){
    uint blockId = get_blockID();
    uint id      = get_threadID(blockId);
    if (id < n)
        d_o[id] = (T)id;
}

template<typename T>
void GMemOpers<T>::SetLinear(T* d_o, size_t n, StreamT stream) {
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids = make_grid(iDivUp(n,threads.x));
    SetLinear_kernel<<<grids, threads, 0, stream>>>(d_o, n);
}


/////////////////////////////////////////////////////////////////////////////////
// SetLinearDown
//  Initiate the value of an unsigned array with a linear ram down
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void SetLinearDown_kernel(T* d_o, uint n){
    uint blockId = get_blockID();
    uint id      = get_threadID(blockId);
    if (id < n)
        d_o[id] = (T)(n - id - 1);
}

template<typename T>
void GMemOpers<T>::SetLinearDown(T* d_o,  size_t n, StreamT stream) {
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));
    SetLinearDown_kernel<<<grids, threads, 0, stream>>>(d_o, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Comp_unary
//  Return the result on a single operation on the input : abs, negative, sqrt ...
//////////////////////////////////////////////////////////////////////////////////
template<typename T, class trait>
__global__ void Comp_unary_kernel(T* d_o, const T* d_i, int n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] = trait::op(d_i[id]);
    }
}

template<typename T>
template<MATH_OPS op>
void GMemOpers<T>::Comp_unary(T* d_o, const T* d_i, size_t n, StreamT stream){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));
    if(!MOpers<T, op>::valid()){
       throw PyCAException(__FILE__, __LINE__, "Unimplemented math operation for specified type");
    }
    Comp_unary_kernel<T, MOpers<T, op> ><<<grids, threads, 0, stream>>>(d_o, d_i, n);
}

// masked version

template<typename T, class trait>
__global__ 
void 
Comp_unary_masked_kernel(T* d_o, const T* d_i, const T* d_mask, int n)
{
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
	if(d_mask[id])
	    d_o[id] = trait::op(d_i[id]);
    }
}

template<typename T>
template<MATH_OPS op>
void GMemOpers<T>::
Comp_unary(T* d_o, const T* d_i, const T* d_mask,
	   size_t n, StreamT stream)
{
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));
    if(!MOpers<T, op>::valid()){
       throw PyCAException(__FILE__, __LINE__, 
			   "Unimplemented math operation for specified type");
    }
    Comp_unary_masked_kernel<T, MOpers<T, op> >
	<<<grids, threads, 0, stream>>>
	(d_o, d_i, d_mask, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Comp_unary Inplace version
//////////////////////////////////////////////////////////////////////////////////
template<typename T, class trait>
__global__ void Comp_unary_kernel_I(T* d_o, int n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] = trait::op(d_o[id]);
    }
}

template<typename T>
template<MATH_OPS op>
void GMemOpers<T>::Comp_unary_I(T* d_o, size_t n, StreamT stream){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));
    checkConfig(grids);
    if(!MOpers<T, op>::valid()){
       throw PyCAException(__FILE__, __LINE__, "Unimplemented math operation for specified type");
    }
    Comp_unary_kernel_I<T, MOpers<T, op> ><<<grids, threads, 0, stream>>>(d_o, n);
}

// masked version

template<typename T, class trait>
__global__ void 
Comp_unary_masked_kernel_I(T* d_o, const T* d_mask, int n)
{
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
	if(d_mask[id])
	    d_o[id] = trait::op(d_o[id]);
    }
}

template<typename T>
template<MATH_OPS op>
void GMemOpers<T>::
Comp_unary_I(T* d_o, const T* d_mask,
	     size_t n, StreamT stream)
{
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));
    checkConfig(grids);
    if(!MOpers<T, op>::valid()){
       throw PyCAException(__FILE__, __LINE__, 
			   "Unimplemented math operation for specified type");
    }
    Comp_unary_masked_kernel_I<T, MOpers<T, op> >
	<<<grids, threads, 0, stream>>>
	(d_o, d_mask, n);
}

///////////////////////////////////////////////////////////////////////////////
// Generated Unary Function Definitions
///////////////////////////////////////////////////////////////////////////////

GMEM_UNARY_OPERS(Abs)
GMEM_UNARY_OPERS(Cube)
GMEM_UNARY_OPERS(Exp)
GMEM_UNARY_OPERS(Log)
GMEM_UNARY_OPERS(Sqr)
GMEM_UNARY_OPERS(Neg)
GMEM_UNARY_OPERS(Ramp)
GMEM_UNARY_OPERS(Sgn)
GMEM_UNARY_OPERS(Step)
GMEM_UNARY_OPERS(Sqrt)
GMEM_UNARY_OPERS(Inv)
GMEM_UNARY_OPERS(Sin)
GMEM_UNARY_OPERS(Asin)
GMEM_UNARY_OPERS(Cos)
GMEM_UNARY_OPERS(Acos)
GMEM_UNARY_OPERS(Tan)
GMEM_UNARY_OPERS(Atan)
GMEM_UNARY_OPERS(Csc)
GMEM_UNARY_OPERS(Sec)
GMEM_UNARY_OPERS(Cot)
GMEM_UNARY_OPERS(Ceil)
GMEM_UNARY_OPERS(Floor)
GMEM_UNARY_OPERS(Round)

/////////////////////////////////////////////////////////////////////////////////
// Binary function with constant
//////////////////////////////////////////////////////////////////////////////////
template<typename T, class op>
__global__ void binaryC_kernel(T* d_o, const T* d_i, T c, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] = op::op(d_i[id], c);
    }
}

template<class T, class op>
__global__ void binaryC_kernel(T* d_o, const T* d_i, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] = op::op(d_i[id], fetch(0,(T*)NULL));
    }
}


template<typename T>
template<MATH_OPS op>
void GMemOpers<T>::binaryC(T* d_o, const T* d_i, const T& c, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));
    if(!MOpers<T, op>::valid()){
       throw PyCAException(__FILE__, __LINE__, "Unimplemented math operation for specified type");
    }
    if (onDev) {
        cache_bind(&c);
        binaryC_kernel<T, MOpers<T, op> ><<<grids, threads, 0, stream>>>(d_o, d_i, n);
    } else {
        binaryC_kernel<T, MOpers<T, op> ><<<grids, threads, 0, stream>>>(d_o, d_i, c, n);
    }
}

// masked versions

template<typename T, class op>
__global__ void 
binaryC_masked_kernel(T* d_o, const T* d_i, T c, const T* d_mask,
		      uint n)
{
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
	if(d_mask[id])
	    d_o[id] = op::op(d_i[id], c);
    }
}

template<class T, class op>
__global__ void 
binaryC_masked_kernel(T* d_o, const T* d_i, const T* d_mask,
		      uint n)
{
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
	if(d_mask[id])
	    d_o[id] = op::op(d_i[id], fetch(0,(T*)NULL));
    }
}


template<typename T>
template<MATH_OPS op>
void GMemOpers<T>::
binaryC(T* d_o, const T* d_i, const T& c, const T* d_mask,
	size_t n, StreamT stream, bool onDev)
{
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));
    if(!MOpers<T, op>::valid()){
       throw PyCAException(__FILE__, __LINE__, 
			   "Unimplemented math operation for specified type");
    }
    if (onDev) {
        cache_bind(&c);
        binaryC_masked_kernel<T, MOpers<T, op> >
	    <<<grids, threads, 0, stream>>>(d_o, d_i, d_mask, n);
    } else {
        binaryC_masked_kernel<T, MOpers<T, op> >
	    <<<grids, threads, 0, stream>>>(d_o, d_i, c, d_mask, n);
    }
}


/////////////////////////////////////////////////////////////////////////////////
// Binary in-place function with constant
//////////////////////////////////////////////////////////////////////////////////
template<typename T, class op>
__global__ void binaryC_I_kernel(T* d_o, T c, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        op::iop(d_o[id], c);
    }
}

template<class T, class op>
__global__ void binaryC_I_kernel(T* d_o, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        op::iop(d_o[id], fetch(0,(T*)NULL));
    }
}


template<typename T>
template<MATH_OPS op>
void GMemOpers<T>::binaryC_I(T* d_o, const T& c, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));
    if(!MOpers<T, op>::valid()){
       throw PyCAException(__FILE__, __LINE__, "Unimplemented math operation for specified type");
    }
    if (onDev)  {
        cache_bind(&c);
        binaryC_I_kernel<T, MOpers<T, op> ><<<grids, threads, 0, stream>>>(d_o, n);
    } else {
        binaryC_I_kernel<T, MOpers<T, op> ><<<grids, threads, 0, stream>>>(d_o, c, n);
    }

}

// masked versions

template<typename T, class op>
__global__ void 
binaryC_I_masked_kernel(T* d_o, T c, const T* d_mask, 
		 uint n)
{
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
	if(d_mask[id])
	    op::iop(d_o[id], c);
    }
}

template<class T, class op>
__global__ void 
binaryC_I_masked_kernel(T* d_o, const T* d_mask, 
		 uint n)
{
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
	if(d_mask[id])
	    op::iop(d_o[id], fetch(0,(T*)NULL));
    }
}


template<typename T>
template<MATH_OPS op>
void GMemOpers<T>::
binaryC_I(T* d_o, const T& c, const T* d_mask,
	  size_t n, StreamT stream, bool onDev)
{
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));
    if(!MOpers<T, op>::valid()){
       throw PyCAException(__FILE__, __LINE__, 
			   "Unimplemented math operation for specified type");
    }
    if (onDev)  {
        cache_bind(&c);
        binaryC_I_masked_kernel<T, MOpers<T, op> >
	    <<<grids, threads, 0, stream>>>(d_o, d_mask, n);
    } else {
        binaryC_I_masked_kernel<T, MOpers<T, op> >
	    <<<grids, threads, 0, stream>>>(d_o, c, d_mask, n);
    }

}


/////////////////////////////////////////////////////////////////////////////////
// Binary function with image
//////////////////////////////////////////////////////////////////////////////////
template<typename T, class op>
__global__ void binary_kernel(T* d_o, const T* d_i, const T* d_i1, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] = op::op(d_i[id], d_i1[id]);
    }
}

template<typename T>
template<MATH_OPS op>
void GMemOpers<T>::binary(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT stream){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));
    if(!MOpers<T, op>::valid()){
       throw PyCAException(__FILE__, __LINE__, "Unimplemented math operation for specified type");
    }
    binary_kernel<T, MOpers<T, op> ><<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);
}

// masked version

template<typename T, class op>
__global__ void 
binary_masked_kernel(T* d_o, const T* d_i, const T* d_i1, const T* d_mask,
	      uint n)
{
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
	if(d_mask[id])
	    d_o[id] = op::op(d_i[id], d_i1[id]);
    }
}

template<typename T>
template<MATH_OPS op>
void GMemOpers<T>::
binary(T* d_o, const T* d_i, const T* d_i1, const T* d_mask,
       size_t n, StreamT stream)
{
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));
    if(!MOpers<T, op>::valid()){
       throw PyCAException(__FILE__, __LINE__, 
			   "Unimplemented math operation for specified type");
    }
    binary_masked_kernel<T, MOpers<T, op> >
	<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, d_mask, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Binary in-place function with image
//////////////////////////////////////////////////////////////////////////////////
template<typename T, class op>
__global__ void binary_I_kernel(T* d_o, const T* d_i, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        op::iop(d_o[id], d_i[id]);
    }
}

template<typename T>
template<MATH_OPS op>
void GMemOpers<T>::binary_I(T* d_o, const T* d_i, size_t n, StreamT stream){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));
    if(!MOpers<T, op>::valid()){
       throw PyCAException(__FILE__, __LINE__, "Unimplemented math operation for specified type");
    }
    binary_I_kernel<T, MOpers<T, op> ><<<grids, threads, 0, stream>>>(d_o, d_i, n);
}

// masked versions

template<typename T, class op>
__global__ void 
binary_I_masked_kernel(T* d_o, const T* d_i, const T* d_mask,
		       uint n)
{
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
	if(d_mask[id])
	    op::iop(d_o[id], d_i[id]);
    }
}

template<typename T>
template<MATH_OPS op>
void GMemOpers<T>::binary_I(T* d_o, const T* d_i, const T* d_mask,
			    size_t n, StreamT stream)
{
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));
    if(!MOpers<T, op>::valid()){
       throw PyCAException(__FILE__, __LINE__, "Unimplemented math operation for specified type");
    }
    binary_I_masked_kernel<T, MOpers<T, op> >
	<<<grids, threads, 0, stream>>>(d_o, d_i, d_mask, n);
}

///////////////////////////////////////////////////////////////////////////////
// Generated Binary Function Definitions
///////////////////////////////////////////////////////////////////////////////

GMEM_BINARY_OPERS(Add)
GMEM_BINARY_OPERS(Sub)
GMEM_BINARY_OPERS(Mul)
GMEM_BINARY_OPERS(Div)
GMEM_BINARY_OPERS(Max)
GMEM_BINARY_OPERS(Min)
GMEM_BINARY_OPERS(LT)
GMEM_BINARY_OPERS(LTE)
GMEM_BINARY_OPERS(EQ)
GMEM_BINARY_OPERS(NEQ)
GMEM_BINARY_OPERS(GT)
GMEM_BINARY_OPERS(GTE)
GMEM_BINARY_OPERS(Atan2)

// Pow doesn't have image version
GMEM_BINARY_OPERC(Pow)
GMEM_BINARY_OPERC_I(Pow)
GMEM_BINARY_OPER_NOC(Pow)
GMEM_BINARY_OPER_NOC_I(Pow)


/////////////////////////////////////////////////////////////////////////////////
// Absolute value of the difference
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void AbsDiff_kernel(T* d_o , const T* d_i, const T* d_i1, uint n){
    uint blockId = get_blockID();
    uint id      = get_threadID(blockId);
    if (id < n){
        d_o[id] = (d_i[id] >= d_i1[id]) ? (d_i[id] - d_i1[id]) : (d_i1[id] - d_i[id]);
    }
}
template<typename T>
void GMemOpers<T>::AbsDiff(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT stream){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));
    AbsDiff_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Absolute value of the difference
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void AbsDiff_kernel_I(T* d_o, const T* d_i, uint n){
    uint blockId = get_blockID();
    uint id      = get_threadID(blockId);
    if (id < n){
        d_o[id] = (d_o[id] >= d_i[id]) ? d_o[id] - d_i[id] : d_o[id] = d_i[id] - d_o[id];
    }
}

template<typename T>
void GMemOpers<T>::AbsDiff_I(T* d_o, const T* d_i, size_t n, StreamT stream){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));
    AbsDiff_kernel_I<<<grids, threads, 0, stream>>>(d_o, d_i, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Square value of the difference
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void SqrDiff_kernel(T* d_o , const T* d_i, const T* d_i1, uint n){
    uint blockId = get_blockID();
    uint id      = get_threadID(blockId);
    if (id < n){
        d_o[id] = (d_i[id] - d_i1[id]) * (d_i[id] - d_i1[id]);
    }
}

template<typename T>
void GMemOpers<T>::SqrDiff(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT stream){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));
    SqrDiff_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Square value of the difference
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void SqrDiff_kernel_I(T* d_o, const T* d_i, uint n){
    uint blockId = get_blockID();
    uint id      = get_threadID(blockId);
    if (id < n){
        d_o[id] = (d_o[id] - d_i[id]) * (d_o[id] - d_i[id]);
    }
}

template<typename T>
void GMemOpers<T>::SqrDiff_I(T* d_o, const T* d_i, size_t n, StreamT stream){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));
    SqrDiff_kernel_I<<<grids, threads, 0, stream>>>(d_o, d_i, n);
}

/////////////////////////////////////////////////////////////////////////////////
// d_o = d_i * (d_i1 * c)
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void MulMulC_kernel(T* d_o, const T* d_i, const T* d_i1, T c, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id <n){
        d_o[id] = d_i[id] * (d_i1[id] * c);
    }
}

template<class T>
__global__ void MulMulC_kernel(T* d_o, const T* d_i, const T* d_i1, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id <n){
        d_o[id] = d_i[id] * (d_i1[id] * fetch(0,(T*)NULL));
    }
}

template<typename T>
void GMemOpers<T>::MulMulC(T* d_o, const T* d_i, const T* d_i1, const T& c, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n, threads.x));

    if (onDev) {
        cache_bind(&c);
        MulMulC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);
    } else
        MulMulC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, c, n);
}

template<class T>
__global__ void MulMulC_I_kernel(T* d_o, const T* d_i, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id <n){
        d_o[id] *= (d_i[id] * fetch(0,(T*)NULL));
    }
}

template<typename T>
__global__ void MulMulC_I_kernel(T* d_o, const T* d_i, T c, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id <n){
        d_o[id] *= (d_i[id] * c);
    }
}

template<typename T>
void GMemOpers<T>::MulMulC_I(T* d_o, const T* d_i, const T& c, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n, threads.x));

    if (onDev)
    {
        cache_bind(&c);
        MulMulC_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, n);
    }else
        MulMulC_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, c, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Multiply three arrays together
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void MulMul_kernel(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id <n){
        d_o[id] = d_i[id] * (d_i1[id] * d_i2[id]);
    }
}

template<typename T>
void GMemOpers<T>::MulMul(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, size_t n, StreamT stream){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));

    MulMul_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, d_i2, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Multiply three arrays together inplace version
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void MulMul_I_kernel(T* d_o, const T* d_i, const T* d_i1, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id <n){
        d_o[id] *=  (d_i[id] * d_i1[id]);
    }
}

template<typename T>
void GMemOpers<T>::MulMul_I(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT stream){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));

    MulMul_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);
}


/////////////////////////////////////////////////////////////////////////////////
// Add two array together and divide by the third array
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void AddDiv_kernel(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = (d_i[id] + d_i1[id]) / d_i2[id];
}

template<typename T>
void GMemOpers<T>::AddDiv(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, size_t n, StreamT stream){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));

    AddDiv_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, d_i2, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Add two array together and divide by the third array inplace
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void AddDiv_I_kernel(T* d_o, const T* d_i, const T* d_i1, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = (d_o[id] + d_i[id]) / d_i1[id];
}

template<typename T>
void GMemOpers<T>::AddDiv_I(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT stream){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));
    AddDiv_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Sub two array together and divide by the third array
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void SubDiv_kernel(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = (d_i[id] - d_i1[id]) / d_i2[id];
}

template<typename T>
void GMemOpers<T>::SubDiv(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, size_t n, StreamT stream){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));

    SubDiv_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, d_i2, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Sub two array together and divide by the third array inplace
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void SubDiv_I_kernel(T* d_o, const T* d_i, const T* d_i1, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = (d_o[id] - d_i[id]) / d_i1[id];
}

template<typename T>
void GMemOpers<T>::SubDiv_I(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT stream){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));
    SubDiv_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);
}


/////////////////////////////////////////////////////////////////////////////////
// Multiply two array together and add the third one
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void MulAdd_kernel(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = d_i[id] * d_i1[id] + d_i2[id];
}

template<typename T>
void GMemOpers<T>::MulAdd(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, size_t n, StreamT stream){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));

    MulAdd_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, d_i2, n);
}


/////////////////////////////////////////////////////////////////////////////////
// Multiply two array together and add the third one inplace version
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void MulAdd_I_kernel(T* d_o, const T* d_i, const T* d_i1, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = d_o[id] * d_i[id] + d_i1[id];
}

template<typename T>
void GMemOpers<T>::MulAdd_I(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT stream){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));

    MulAdd_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Multiply two array together and sub the third one
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void MulSub_kernel(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = d_i[id] * d_i1[id] - d_i2[id];
}

template<typename T>
void GMemOpers<T>::MulSub(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, size_t n, StreamT stream){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));

    MulSub_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, d_i2, n);
}


/////////////////////////////////////////////////////////////////////////////////
// Multiply two array together and sub the third one inplace version
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void MulSub_I_kernel(T* d_o, const T* d_i, const T* d_i1, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = d_o[id] * d_i[id] - d_i1[id];
}

template<typename T>
void GMemOpers<T>::MulSub_I(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT stream){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));

    MulSub_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Add two arrays and multiply by the third one
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void AddMul_kernel(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = (d_i[id] + d_i1[id]) * d_i2[id];
}

template<typename T>
void GMemOpers<T>::AddMul(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, size_t n, StreamT stream){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));

    AddMul_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, d_i2, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Add two array together and multiply by the third one inplace version
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void AddMul_I_kernel(T* d_o, const T* d_i, const T* d_i1, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = (d_o[id] + d_i[id]) * d_i1[id];
}

template<typename T>
void GMemOpers<T>::AddMul_I(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT stream){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));

    AddMul_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Sub two arrays and multiply by the third one
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void SubMul_kernel(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = (d_i[id] - d_i1[id]) * d_i2[id];
}

template<typename T>
void GMemOpers<T>::SubMul(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, size_t n, StreamT stream){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));

    SubMul_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, d_i2, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Sub two array together and multiply by the third one inplace version
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void SubMul_I_kernel(T* d_o, const T* d_i, const T* d_i1, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = (d_o[id] - d_i[id]) * d_i1[id];
}

template<typename T>
void GMemOpers<T>::SubMul_I(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT stream){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));

    SubMul_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Add an array by a constant then multiply by other constant
// Used with normalized function
//////////////////////////////////////////////////////////////////////////////////
template<class T>
__global__ void AddCMulC_kernel(T* d_o, const T* d_i, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = (d_i[id] + fetch_y(0, (T*)NULL)) * fetch_z(0, (T*)NULL);
}

template<typename T>
__global__ void AddCMulC_kernel(T* d_o, const T* d_i, T a, T b, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = (d_i[id] + a) * b;
}

template<typename T>
void GMemOpers<T>::AddCMulC(T* d_o, const T* d_i, const T& a, const T& b, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));

    if (onDev) {
        cache_bind_y(&a);
        cache_bind_z(&b);
        AddCMulC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, n);
    } else {
        AddCMulC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, a, b, n);
    }
}

template<typename T>
__global__ void AddCMulC_I_kernel(T* d_o, T a, T b, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = (d_o[id] + a) * b;
}

template<class T>
__global__ void AddCMulC_I_kernel(T* d_o, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = (d_o[id] + fetch_y(0, (T*)NULL)) * fetch_z(0, (T*)NULL);
}

template<typename T>
void GMemOpers<T>::AddCMulC_I(T* d_o, const T& a, const T& b, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));

    if (onDev) {
        cache_bind_y(&a);
        cache_bind_z(&b);
        AddCMulC_I_kernel<<<grids, threads, 0, stream>>>(d_o, n);
    } else
        AddCMulC_I_kernel<<<grids, threads, 0, stream>>>(d_o, a, b, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Multiply an array by a constant and add the second one
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void MulCAdd_kernel(T* d_o, const T* d_i, T c, const T* d_i1, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = d_i[id] * c + d_i1[id];
}

template<class T>
__global__ void MulCAdd_kernel(T* d_o, const T* d_i, const T* d_i1, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = d_i[id] * fetch(0, (T*)NULL) + d_i1[id];
}

template<typename T>
void GMemOpers<T>::MulCAdd(T* d_o, const T* d_i, const T& c, const T* d_i1, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));
    if (onDev) {
        cache_bind(&c);
        MulCAdd_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);
    } else
        MulCAdd_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, c, d_i1, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Multiply an array by a constant and add the second one inplace version
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void MulCAdd_I_kernel(T* d_o, T c, const T* d_i, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = d_o[id] * c + d_i[id];
}
template<class T>
__global__ void MulCAdd_I_kernel(T* d_o, const T* d_i, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = d_o[id] * fetch(0, (T*)NULL) + d_i[id];
}

template<typename T>
void GMemOpers<T>::MulCAdd_I(T* d_o, const T& c, const T* d_i, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));

    if (onDev) {
        cache_bind(&c);
        MulCAdd_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, n);
    } else
        MulCAdd_I_kernel<<<grids, threads, 0, stream>>>(d_o, c, d_i, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Multiply an array by a constant and sub the second one
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void MulCSub_kernel(T* d_o, const T* d_i, T c, const T* d_i1, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = d_i[id] * c - d_i1[id];
}

template<class T>
__global__ void MulCSub_kernel(T* d_o, const T* d_i, const T* d_i1, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = d_i[id] * fetch(0, (T*)NULL) - d_i1[id];
}

template<typename T>
void GMemOpers<T>::MulCSub(T* d_o, const T* d_i, const T& c, const T* d_i1, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));
    if (onDev) {
        cache_bind(&c);
        MulCSub_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);
    } else
        MulCSub_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, c, d_i1, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Multiply an array by a constant and sub the second one inplace version
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void MulCSub_I_kernel(T* d_o, T c, const T* d_i, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = d_o[id] * c - d_i[id];
}

template<class T>
__global__ void MulCSub_I_kernel(T* d_o, const T* d_i, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = d_o[id] * fetch(0, (T*)NULL) - d_i[id];
}

template<typename T>
void GMemOpers<T>::MulCSub_I(T* d_o, const T& c, const T* d_i, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));
    if (onDev) {
        cache_bind(&c);
        MulCSub_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, n);
    } else
        MulCSub_I_kernel<<<grids, threads, 0, stream>>>(d_o, c, d_i, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Multiply an array by a constant and add the second constant
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void MulCAddC_kernel(T* d_o, const T* d_i, T a, T b, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = d_i[id] * a + b;
}

template<class T>
__global__ void MulCAddC_kernel(T* d_o, const T* d_i, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = d_i[id] * fetch_y(0, (T*)NULL) + fetch_z(0, (T*)NULL);
}

template<typename T>
void GMemOpers<T>::MulCAddC(T* d_o, const T* d_i, const T& a, const T& b, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));

    if (onDev){
        cache_bind_y(&a);
        cache_bind_z(&b);
        MulCAddC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, n);
    } else
        MulCAddC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, a, b, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Multiply an array by a constant and add the second constant inplace version
//////////////////////////////////////////////////////////////////////////////////

template<typename T>
__global__ void MulCAddC_I_kernel(T* d_o, T a, T b, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = d_o[id] * a + b;
}

template<class T>
__global__ void MulCAddC_I_kernel(T* d_o, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = d_o[id]  * fetch_y(0, (T*)NULL) + fetch_z(0, (T*)NULL);
}

template<typename T>
void GMemOpers<T>::MulCAddC_I(T* d_o, const T& a, const T& b, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));

    if (onDev)
    {
        cache_bind_y(&a);
        cache_bind_z(&b);
        MulCAddC_I_kernel<<<grids, threads, 0, stream>>>(d_o, n);
    }else
        MulCAddC_I_kernel<<<grids, threads, 0, stream>>>(d_o, a, b, n);
}


/////////////////////////////////////////////////////////////////////////////////
// Add two array then multiply by a constant
//////////////////////////////////////////////////////////////////////////////////
template<class T>
__global__ void AddMulC_kernel(T* d_o, const T* d_i, const T* d_i1, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = (d_i[id] + d_i1[id]) * fetch(0,(T*)NULL);
}

template<typename T>
__global__ void AddMulC_kernel(T* d_o, const T* d_i, const T* d_i1, T c, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = (d_i[id] + d_i1[id]) * c;
}

template<typename T>
void GMemOpers<T>::AddMulC(T* d_o, const T* d_i, const T* d_i1, const T& c, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));

    if (onDev) {
        cache_bind(&c);
        AddMulC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);
    } else
        AddMulC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, c, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Add two array then multiply by a constant inplace version
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void AddMulC_I_kernel(T* d_o, const T* d_i, T c, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = (d_o[id] + d_i[id]) * c;
}

template<class T>
__global__ void AddMulC_I_kernel(T* d_o, const T* d_i, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = (d_o[id] + d_i[id]) * fetch(0,(T*)NULL);
}

template<typename T>
void GMemOpers<T>::AddMulC_I(T* d_o, const T* d_i, const T& c, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));

    if (onDev) {
        cache_bind(&c);
        AddMulC_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, n);
    } else
        AddMulC_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, c, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Sub two array then multiply by a constant
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void SubMulC_kernel(T* d_o, const T* d_i, const T* d_i1, T c, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = (d_i[id] - d_i1[id]) * c;
}

template<class T>
__global__ void SubMulC_kernel(T* d_o, const T* d_i, const T* d_i1, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = (d_i[id] - d_i1[id]) * fetch(0,(T*)NULL);
}

template<typename T>
void GMemOpers<T>::SubMulC(T* d_o, const T* d_i, const T* d_i1, const T& c, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));

    if (onDev) {
        cache_bind(&c);
        SubMulC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);
    }else
        SubMulC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, c, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Sub two array then multiply by a constant inplace version
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void SubMulC_I_kernel(T* d_o, const T* d_i, T c, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = (d_o[id] - d_i[id]) * c;
}
template<class T>
__global__ void SubMulC_I_kernel(T* d_o, const T* d_i, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = (d_o[id] - d_i[id]) * fetch(0,(T*)NULL);
}

template<typename T>
void GMemOpers<T>::SubMulC_I(T* d_o, const T* d_i, const T& c, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));

    if (onDev) {
        cache_bind(&c);
        SubMulC_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, n);
    } else
        SubMulC_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, c, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Add an arrays with another array that multiply by a constant
//  d_o = d_i + d_i1 * c
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void Add_MulC_kernel(T* d_o, const T* d_i, const T* d_i1, T c, uint n) {
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = d_i[id] + d_i1[id] * c;
}

template<class T>
__global__ void Add_MulC_kernel(T* d_o, const T* d_i, const T* d_i1, uint n) {
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = d_i[id] + d_i1[id] * fetch(0,(T*)NULL);
}

template<typename T>
void GMemOpers<T>::Add_MulC(T* d_o, const T* d_i, const T* d_i1, const T& c, size_t n, StreamT stream, bool onDev) {
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n, threads.x));

    if (onDev) {
        cache_bind(&c);
        Add_MulC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);
    } else
        Add_MulC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, c, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Add an arrays with another array that multiply by a constant
//  d_o = d_i + d_i1 * c
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void Add_MulC_kernel_I(T* d_o, const T* d_i, T c, uint n) {
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] += d_i[id] * c;
}

template<class T>
__global__ void Add_MulC_kernel_I(T* d_o, const T* d_i, uint n) {
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] += d_i[id] * fetch(0,(T*)NULL);
}

template<typename T>
void GMemOpers<T>::Add_MulC_I(T* d_o, const T* d_i, const T& c, size_t n, StreamT stream, bool onDev) {
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n, threads.x));

    if (onDev) {
        cache_bind(&c);
        Add_MulC_kernel_I<<<grids, threads, 0, stream>>>(d_o, d_i, n);
    } else
        Add_MulC_kernel_I<<<grids, threads, 0, stream>>>(d_o, d_i, c, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Add an arrays with another array that multiply by an array
//  d_o = d_i + d_i1 * d_i2
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void Add_Mul_kernel(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, uint n) {
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = d_i[id] + d_i1[id] * d_i2[id];
}

template<typename T>
void GMemOpers<T>::Add_Mul(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, size_t n, StreamT stream) {
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n, threads.x));

    Add_Mul_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, d_i2, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Add an arrays with another array that multiply by an array
//  d_o = d_o + d_i * d_i1
/////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void Add_Mul_kernel_I(T* d_o, const T* d_i, const T* d_i1, uint n) {
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] += d_i[id] * d_i1[id];
}

template<typename T>
void GMemOpers<T>::Add_Mul_I(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT stream) {
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n, threads.x));

    Add_Mul_kernel_I<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Sub an arrays with another array that multiply by an array
//  d_o = d_i + d_i1 * d_i2
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void Sub_Mul_kernel(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, uint n) {
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] = d_i[id] - d_i1[id] * d_i2[id];
}

template<typename T>
void GMemOpers<T>::Sub_Mul(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, size_t n, StreamT stream) {
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n, threads.x));

    Sub_Mul_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, d_i2, n);
}

/////////////////////////////////////////////////////////////////////////////////
// Sub an arrays with another array that multiply by an array
//  d_o = d_o + d_i * d_i1
/////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void Sub_Mul_kernel_I(T* d_o, const T* d_i, const T* d_i1, uint n) {
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n)
        d_o[id] -= d_i[id] * d_i1[id];
}

template<typename T>
void GMemOpers<T>::Sub_Mul_I(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT stream) {
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n, threads.x));

    Sub_Mul_kernel_I<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);
}


/////////////////////////////////////////////////////////////////////////////////
// Quadary function
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void MulC_Add_MulC_kernel(T* d_o, const T* d_i, T a, const T* d_i1, T b,  uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] = d_i[id] * a + d_i1[id] * b;
    }
}

template<class T>
__global__ void MulC_Add_MulC_kernel(T* d_o, const T* d_i, const T* d_i1,  uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] = d_i[id] * fetch_y(0, (T*)NULL) + d_i1[id] * fetch_z(0, (T*)NULL);
    }
}

template<typename T>
void GMemOpers<T>::
MulC_Add_MulC(T* d_o, const T* d_i, const T& a, const T* d_i1, const T& b, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n,threads.x));

    if (onDev) {
        cache_bind_y(&a);
        cache_bind_z(&b);
        MulC_Add_MulC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);
    } else {
        MulC_Add_MulC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, a, d_i1, b, n);
    }
}

template<typename T>
__global__ void MulC_Add_MulC_kernel_I(T* d_o, T a, const T* d_i, T b,  uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] = d_o[id] * a + d_i[id] * b;
    }
}
template<class T>
__global__ void MulC_Add_MulC_kernel_I(T* d_o, const T* d_i,  uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] = d_o[id] * fetch_y(0, (T*)NULL) + d_i[id] * fetch_z(0, (T*)NULL);
    }
}

template<typename T>
void GMemOpers<T>::
MulC_Add_MulC_I(T* d_o, const T& a, const T* d_i, const T& b, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n, threads.x));

    if (onDev)
    {
        cache_bind_y(&a);
        cache_bind_z(&b);
        MulC_Add_MulC_kernel_I<<<grids, threads, 0, stream>>>(d_o, d_i, n);

    } else
        MulC_Add_MulC_kernel_I<<<grids, threads, 0, stream>>>(d_o, a, d_i, b, n);
}

////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void Add_AddMulC_kernel(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, T c, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] = d_i[id] + (d_i1[id] + d_i2[id]) * c;
    }
}

template<class T>
__global__ void Add_AddMulC_kernel(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] = d_i[id] + (d_i1[id] + d_i2[id]) * fetch(0,(T*)NULL);
    }
}

template<typename T>
void GMemOpers<T>::Add_AddMulC(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, const T& c, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n, threads.x));

    if (onDev) {
        cache_bind(&c);
        Add_AddMulC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, d_i2, n);
    } else
        Add_AddMulC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, d_i2, c, n);
}

////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void Add_AddMulC_I_kernel(T* d_o, const T* d_i, const T* d_i1, T c, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] += (d_i[id] + d_i1[id]) * c;
    }
}

template<class T>
__global__ void Add_AddMulC_I_kernel(T* d_o, const T* d_i, const T* d_i1, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] += (d_i[id] + d_i1[id]) * fetch(0,(T*)NULL);
    }
}

template<typename T>
void GMemOpers<T>::Add_AddMulC_I(T* d_o, const T* d_i, const T* d_i1, const T& c, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n, threads.x));

    if (onDev) {
        cache_bind(&c);
        Add_AddMulC_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);
    } else
        Add_AddMulC_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, c, n);
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
template<class T>
__global__ void Add_SubMulC_kernel(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] = d_i[id] + (d_i1[id] - d_i2[id]) * fetch(0,(T*)NULL);
    }
}

template<typename T>
__global__ void Add_SubMulC_kernel(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, T c, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] = d_i[id] + (d_i1[id] - d_i2[id]) * c;
    }
}

template<typename T>
void GMemOpers<T>::Add_SubMulC(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, const T& c, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n, threads.x));

    if (onDev){
        cache_bind(&c);
        Add_SubMulC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, d_i2, n);
    } else
        Add_SubMulC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, d_i2, c, n);
}

////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void Add_SubMulC_I_kernel(T* d_o, const T* d_i, const T* d_i1, T c, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] += (d_i[id] - d_i1[id]) * c;
    }
}
template<class T>
__global__ void Add_SubMulC_I_kernel(T* d_o, const T* d_i, const T* d_i1, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] += (d_i[id] - d_i1[id]) * fetch(0,(T*)NULL);
    }
}

template<typename T>
void GMemOpers<T>::Add_SubMulC_I(T* d_o, const T* d_i, const T* d_i1, const T& c, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n, threads.x));

    if (onDev) {
        cache_bind(&c);
        Add_SubMulC_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);
    } else
        Add_SubMulC_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, c, n);
}

////////////////////////////////////////////////////////////////////////////////
//// d_o = d_i + (d_i1 * d_i2 * c)
////////////////////////////////////////////////////////////////////////////////
template<class T>
__global__ void Add_MulMulC_kernel(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] = d_i[id] + (d_i1[id] * d_i2[id]) * fetch(0,(T*)NULL);
    }
}

template<typename T>
__global__ void Add_MulMulC_kernel(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, T c, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] = d_i[id] + (d_i1[id] * d_i2[id]) * c;
    }
}

template<typename T>
void GMemOpers<T>::Add_MulMulC(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, const T& c, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n, threads.x));

    if (onDev)
    {
        cache_bind(&c);
        Add_MulMulC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, d_i2, n);
    }else
        Add_MulMulC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, d_i2, c, n);
}

////////////////////////////////////////////////////////////////////////////////
// d_o = d_o + (d_i * d_i1 * d)
////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void Add_MulMulC_I_kernel(T* d_o, const T* d_i, const T* d_i1, T c, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] += (d_i[id] * d_i1[id]) * c;
    }
}

template<class T>
__global__ void Add_MulMulC_I_kernel(T* d_o, const T* d_i, const T* d_i1, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] += (d_i[id] * d_i1[id]) * fetch(0,(T*)NULL);
    }
}

template<typename T>
void GMemOpers<T>::Add_MulMulC_I(T* d_o, const T* d_i, const T* d_i1, const T& c, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n, threads.x));

    if (onDev) {
        cache_bind(&c);
        Add_MulMulC_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);
    } else
        Add_MulMulC_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, c, n);
}

////////////////////////////////////////////////////////////////////////////////
// d_o = d_i - (d_i1 + d_i2) * c
////////////////////////////////////////////////////////////////////////////////
template<class T>
__global__ void Sub_AddMulC_kernel(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] = d_i[id] - (d_i1[id] + d_i2[id]) * fetch(0,(T*)NULL);
    }
}

template<typename T>
__global__ void Sub_AddMulC_kernel(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, T c, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] = d_i[id] - (d_i1[id] + d_i2[id]) * c;
    }
}

template<typename T>
void GMemOpers<T>::Sub_AddMulC(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, const T& c, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n, threads.x));

    if (onDev){
        cache_bind(&c);
        Sub_AddMulC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, d_i2, n);

    } else
        Sub_AddMulC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, d_i2, c, n);
}

////////////////////////////////////////////////////////////////////////////////
// d_o -= (d_i + d_i1) * c
////////////////////////////////////////////////////////////////////////////////
template<class T>
__global__ void Sub_AddMulC_I_kernel(T* d_o, const T* d_i, const T* d_i1, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] -= (d_i[id] + d_i1[id]) * fetch(0,(T*)NULL);
    }
}

template<typename T>
__global__ void Sub_AddMulC_I_kernel(T* d_o, const T* d_i, const T* d_i1, T c, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] -= (d_i[id] + d_i1[id]) * c;
    }
}

template<typename T>
void GMemOpers<T>::Sub_AddMulC_I(T* d_o, const T* d_i, const T* d_i1, const T& c, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n, threads.x));

    if( onDev){
        cache_bind(&c);
        Sub_AddMulC_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);

    } else
        Sub_AddMulC_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, c, n);
}
////////////////////////////////////////////////////////////////////////////////
// d_o = d_i - (d_i1 - d_i2) * c
////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void Sub_SubMulC_kernel(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, T c, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] = d_i[id] - (d_i1[id] - d_i2[id]) * c;
    }
}

template<class T>
__global__ void Sub_SubMulC_kernel(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] = d_i[id] - (d_i1[id] - d_i2[id]) * fetch(0,(T*)NULL);
    }
}

template<typename T>
void GMemOpers<T>::Sub_SubMulC(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, const T& c, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n, threads.x));

    if (onDev){
        cache_bind(&c);
        Sub_SubMulC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, d_i2, n);
    } else
        Sub_SubMulC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, d_i2, c, n);
}

////////////////////////////////////////////////////////////////////////////////
// d_o -= (d_i - d_i1) * c
////////////////////////////////////////////////////////////////////////////////
template<class T>
__global__ void Sub_SubMulC_I_kernel(T* d_o, const T* d_i, const T* d_i1, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] -= (d_i[id] - d_i1[id]) * fetch(0,(T*)NULL);
    }
}

template<typename T>
__global__ void Sub_SubMulC_I_kernel(T* d_o, const T* d_i, const T* d_i1, T c, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] -= (d_i[id] - d_i1[id]) * c;
    }
}

template<typename T>
void GMemOpers<T>::Sub_SubMulC_I(T* d_o, const T* d_i, const T* d_i1, const T& c, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n, threads.x));

    if (onDev){
        cache_bind(&c);
        Sub_SubMulC_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);

    } else
        Sub_SubMulC_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, c, n);
}

////////////////////////////////////////////////////////////////////////////////
//// d_o = d_i + (d_i1 * d_i2 * c)
////////////////////////////////////////////////////////////////////////////////
template<class T>
__global__ void Sub_MulMulC_kernel(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] = d_i[id] - (d_i1[id] * d_i2[id]) * fetch(0,(T*)NULL);
    }
}

template<typename T>
__global__ void Sub_MulMulC_kernel(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, T c, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] = d_i[id] - (d_i1[id] * d_i2[id]) * c;
    }
}

template<typename T>
void GMemOpers<T>::Sub_MulMulC(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, const T& c, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n, threads.x));

    if (onDev){
        cache_bind(&c);
        Sub_MulMulC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, d_i2, n);
    } else
        Sub_MulMulC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, d_i2, c, n);
}

////////////////////////////////////////////////////////////////////////////////
// d_o = d_o + (d_i * d_i1 * d)
////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void Sub_MulMulC_I_kernel(T* d_o, const T* d_i, const T* d_i1, T c, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] -= (d_i[id] * d_i1[id]) * c;
    }
}
template<class T>
__global__ void Sub_MulMulC_I_kernel(T* d_o, const T* d_i, const T* d_i1, uint n){
    uint blockId = get_blockID();
    uint      id = get_threadID(blockId);
    if (id < n){
        d_o[id] -= (d_i[id] * d_i1[id]) * fetch(0,(T*)NULL);
    }
}

template<typename T>
void GMemOpers<T>::Sub_MulMulC_I(T* d_o, const T* d_i, const T* d_i1, const T& c, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n, threads.x));

    if (onDev){
        cache_bind(&c);
        Sub_MulMulC_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, n);

    } else
        Sub_MulMulC_I_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, d_i1, c, n);
}

////////////////////////////////////////////////////////////////////////////////
// d_o = (d_i + a) * b + c
////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void AddCMulCAddC_kernel(T* d_o, const T* d_i, T a, T b, T c, uint n){
    uint blockId = get_blockID();
    uint id      = get_threadID(blockId);
    if (id < n)
        d_o[id] = (d_i[id] + a) * b + c;
}
template<class T>
__global__ void AddCMulCAddC_kernel(T* d_o, const T* d_i, uint n){
    uint blockId = get_blockID();
    uint id      = get_threadID(blockId);
    if (id < n)
        d_o[id] = d_i[id] + fetch_y(0,(T*)NULL) * fetch_z(0,(T*)NULL) + fetch(0,(T*)NULL);
}

template<typename T>
void GMemOpers<T>::
AddCMulCAddC(T* d_o, const T* d_i, const T& a, const T& b, const T& c, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n, threads.x));

    if (onDev) {
        cache_bind(&c);
        cache_bind_y(&a);
        cache_bind_z(&b);

        AddCMulCAddC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, n);

    } else
        AddCMulCAddC_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, a, b, c, n);
}

////////////////////////////////////////////////////////////////////////////////
// d_o = (d_o + a) * b + c
////////////////////////////////////////////////////////////////////////////////
template<typename T>
__global__ void AddCMulCAddC_I_kernel(T* d_o, T a, T b, T c, uint n){
    uint blockId = get_blockID();
    uint id      = get_threadID(blockId);
    if (id < n)
        d_o[id] = (d_o[id] + a) * b + c;
}

template<class T>
__global__ void AddCMulCAddC_I_kernel(T* d_o, uint n){
    uint blockId = get_blockID();
    uint id      = get_threadID(blockId);
    if (id < n)
        d_o[id] = d_o[id] + fetch_y(0,(T*)NULL) * fetch_z(0,(T*)NULL) + fetch(0,(T*)NULL);
}


template<typename T>
void GMemOpers<T>::AddCMulCAddC_I(T* d_o, const T& a, const T& b, const T& c, size_t n, StreamT stream, bool onDev){
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n, threads.x));

    if (onDev){
        cache_bind(&c);
        cache_bind_y(&a);
        cache_bind_z(&b);

        AddCMulCAddC_I_kernel<<<grids, threads, 0, stream>>>(d_o, n);

    } else
        AddCMulCAddC_I_kernel<<<grids, threads, 0, stream>>>(d_o, a, b, c, n);
}

template<typename T>
__global__ void ReverseOrder_kernel(T* d_o, const T* d_i, uint n){
    uint blockId = get_blockID();
    uint id      = get_threadID(blockId);
    if (id < n)
        d_o[id] = (d_i[n -1 - id]);
}

template<typename T>
void GMemOpers<T>::ReverseOrder(T* d_o, const T* d_i, size_t n, StreamT stream)
{
    dim3 threads(REG_BLOCK_SIZE);
    dim3 grids=make_grid(iDivUp(n, threads.x));
    ReverseOrder_kernel<<<grids, threads, 0, stream>>>(d_o, d_i, n);
}


#define TILE_DIM 8
template<typename T>
__global__ void ShiftRightCoordinate_shared_kernel(T* d_o, const T* d_i,
                                                   uint sizeX, uint sizeY, uint sizeZ)
{
    __shared__ T sdata[TILE_DIM][TILE_DIM+1][TILE_DIM+1];

    const uint iPlaneSize = sizeX * sizeY;
    const uint oPlaneSize = sizeX * sizeZ;

    uint bx = blockIdx.x * TILE_DIM;
    uint by = blockIdx.y * TILE_DIM;

    uint tx = threadIdx.x;
    uint ty = threadIdx.y;

    uint iid = bx + tx + (by + ty) * sizeX;

    for (uint bz=0; bz < sizeZ; bz+= TILE_DIM){
        if ((bx + tx < sizeX) && (by + ty < sizeY))
            for (size_t tz=0; tz < TILE_DIM && (bz + tz < sizeZ); ++tz, iid+=iPlaneSize)
                sdata[tz][ty][tx] = d_i[iid];
        __syncthreads();

        uint oid = bz + ty + ((bx + tx) + by * sizeX) * sizeZ;
        if ((bz + ty < sizeZ) && (bx + tx < sizeX))
            for (uint tz = 0; tz < TILE_DIM && (by + tz < sizeY); ++tz, oid+=oPlaneSize)
                d_o[oid] = sdata[ty][tz][tx];
    }
}

template<typename T>
__global__ void ShiftLeftCoordinate_shared_kernel(T* d_o, const T* d_i,
                                                  uint sizeX, uint sizeY, uint sizeZ)
{
    __shared__ T sdata[TILE_DIM][TILE_DIM+1][TILE_DIM+1];

    const uint iPlaneSize = sizeX * sizeY;
    const uint oPlaneSize = sizeY * sizeZ;

    uint bx = blockIdx.x * TILE_DIM;
    uint by = blockIdx.y * TILE_DIM;

    uint tx = threadIdx.x;
    uint ty = threadIdx.y;

    uint bz = 0;
    uint iid = bx + tx + (by + ty) * sizeX;

    while (bz < sizeZ){
        if ((bx + tx < sizeX) && (by + ty < sizeY))
            for (uint tz = 0; tz < TILE_DIM && (bz + tz < sizeZ); ++tz, iid +=iPlaneSize)
                sdata[tz][ty][tx] = d_i[iid];

        __syncthreads();

        uint oid = by + bz * sizeY + bx * sizeY * sizeZ + (tx + ty * sizeY);

        if ((by + tx < sizeY) && (bz + ty < sizeZ))
            for (uint tz = 0; tz < TILE_DIM && (bx + tz < sizeX); ++tz, oid += oPlaneSize)
                d_o[oid] = sdata[ty][tx][tz];
        bz += TILE_DIM;
    }
}

template<typename T>
void GMemOpers<T>::
ShiftCoordinate(T* d_o, const T* d_i, size_t sizeX, size_t sizeY, size_t sizeZ, bool dir, StreamT stream){
    dim3 threads(8, 8);
    dim3 grids(iDivUp(sizeX, threads.x),iDivUp(sizeY, threads.y));
    if (dir) {
        ShiftRightCoordinate_shared_kernel<<<grids, threads, 0, stream>>>
            (d_o, d_i, sizeX, sizeY, sizeZ);
    }
    else
        ShiftLeftCoordinate_shared_kernel<<<grids, threads, 0, stream>>>
            (d_o, d_i, sizeX, sizeY, sizeZ);
}


template class GMemOpers<float>;
template class GMemOpers<int>;
template class GMemOpers<uint>;
} // end namespace PyCA

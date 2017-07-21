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

#include <CMemOpers.h>
#include <MOpers.h>
#include <algorithm>
#include "PyCAException.h"
#include "mem.h"

#include "MemOperDefs.h"

#define CMEM_UNARY_OPERS(OP) MEM_UNARY_OPERS(CMemOpers, OP)
#define CMEM_BINARY_OPERS(OP) MEM_BINARY_OPERS(CMemOpers, OP)

#define CMEM_BINARY_OPERC(OP) MEM_BINARY_OPERC(CMemOpers, OP)
#define CMEM_BINARY_OPERC_I(OP) MEM_BINARY_OPERC_I(CMemOpers, OP)
#define CMEM_BINARY_OPER_NOC(OP) MEM_BINARY_OPER_NOC(CMemOpers, OP)
#define CMEM_BINARY_OPER_NOC_I(OP) MEM_BINARY_OPER_NOC_I(CMemOpers, OP)
#define CMEM_BINARY_OPER(OP) MEM_BINARY_OPER(CMemOpers, OP)
#define CMEM_BINARY_OPER_I(OP) MEM_BINARY_OPER_I(CMemOpers, OP)

//#pragma GCC diagnostic ignored "-Wunknown-pragmas"

namespace PyCA {

template<typename T>
void CMemOpers<T>::Copy(T* h_o, const T* h_i, size_t n, StreamT stream){
   cpyArrayH2H(h_o, h_i, n);
}

template<typename T>
void CMemOpers<T>::Copy(T* h_o, const T* h_i, const T* h_mask,
			size_t n, StreamT stream) 
{
#pragma omp parallel for
    for (size_t id=0; id< n; ++id)
	if(h_mask[id])
	    h_o[id] = h_i[id];
}

template<typename T>
void CMemOpers<T>::SetMem(T* h_o, const T& c, size_t n, StreamT stream, bool onDev) {
    std::fill(h_o, h_o + n, c);
}

template<typename T>
void CMemOpers<T>::SetMem(T* h_o, const T& c, const T* h_mask, size_t n, 
			  StreamT stream, bool onDev) 
{
#pragma omp parallel for
    for (size_t id=0; id< n; ++id)
	if(h_mask[id])
	    h_o[id] = c;
}

/////////////////////////////////////////////////////////////////////////////////
//  SetLinear
//  Initiate the value of an unsigned array with a linear ram
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::SetLinear(T* h_o, size_t n, StreamT stream) {
#pragma omp parallel for
    for (size_t id=0; id< n; ++id)
        h_o[id] = id;
}


/////////////////////////////////////////////////////////////////////////////////
// SetLinearDown
//  Initiate the value of an unsigned array with a linear ram down
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::SetLinearDown(T* h_o,  size_t n, StreamT stream) {
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = n - 1 - id;
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Comp_unary
//  Return the result on a single operation on the input : abs, negative, sqrt ...
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
template<MATH_OPS op>
void CMemOpers<T>::Comp_unary(T* h_o, const T* h_i, size_t n, StreamT stream){
   if(!MOpers<T, op>::valid()){
      throw PyCAException(__FILE__, __LINE__, "Unimplemented math operation for specified type");
   }
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = MOpers<T, op>::op(h_i[id]);
    }
}

// Masked version
template<typename T>
template<MATH_OPS op>
void 
CMemOpers<T>::
Comp_unary(T* h_o, const T* h_i, const T* h_mask, 
	   size_t n, StreamT stream)
{
    if(!MOpers<T, op>::valid()){
	throw PyCAException(__FILE__, __LINE__, "Unimplemented math operation for specified type");
    }
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
	if(h_mask[id])
	    h_o[id] = MOpers<T, op>::op(h_i[id]);
    }
}


/////////////////////////////////////////////////////////////////////////////////
// Comp_unary Inplace version
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
template<MATH_OPS op>
void CMemOpers<T>::Comp_unary_I(T* h_o, size_t n, StreamT stream){
   if(!MOpers<T, op>::valid()){
      throw PyCAException(__FILE__, __LINE__, "Unimplemented math operation for specified type");
   }
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = MOpers<T, op>::op(h_o[id]);
    }
}

template<typename T>
template<MATH_OPS op>
void CMemOpers<T>::Comp_unary_I(T* h_o, const T* h_mask, 
				size_t n, StreamT stream)
{
    if(!MOpers<T, op>::valid()){
	throw PyCAException(__FILE__, __LINE__, "Unimplemented math operation for specified type");
    }
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
	if(h_mask[id])
	    h_o[id] = MOpers<T, op>::op(h_o[id]);
    }
}

CMEM_UNARY_OPERS(Abs)
CMEM_UNARY_OPERS(Cube)
CMEM_UNARY_OPERS(Exp)
CMEM_UNARY_OPERS(Log)
CMEM_UNARY_OPERS(Sqr)
CMEM_UNARY_OPERS(Neg)
CMEM_UNARY_OPERS(Ramp)
CMEM_UNARY_OPERS(Sgn)
CMEM_UNARY_OPERS(Step)
CMEM_UNARY_OPERS(Sqrt)
CMEM_UNARY_OPERS(Inv)
CMEM_UNARY_OPERS(Sin)
CMEM_UNARY_OPERS(Asin)
CMEM_UNARY_OPERS(Cos)
CMEM_UNARY_OPERS(Acos)
CMEM_UNARY_OPERS(Tan)
CMEM_UNARY_OPERS(Atan)
CMEM_UNARY_OPERS(Csc)
CMEM_UNARY_OPERS(Sec)
CMEM_UNARY_OPERS(Cot)
CMEM_UNARY_OPERS(Ceil)
CMEM_UNARY_OPERS(Floor)
CMEM_UNARY_OPERS(Round)

/////////////////////////////////////////////////////////////////////////////////
// binary function with constant
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
template<MATH_OPS op>
void CMemOpers<T>::binaryC(T* h_o, const T* h_i, const T& c, size_t n, StreamT stream, bool onDev){
   if(!MOpers<T, op>::valid()){
      throw PyCAException(__FILE__, __LINE__, "Unimplemented math operation for specified type");
   }
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = MOpers<T, op>::op(h_i[id], c);
    }
}

template<typename T>
template<MATH_OPS op>
void CMemOpers<T>::binaryC_I(T* h_o, const T& c, size_t n, StreamT stream, bool onDev){
   if(!MOpers<T, op>::valid()){
      throw PyCAException(__FILE__, __LINE__, "Unimplemented math operation for specified type");
   }
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        MOpers<T, op>::iop(h_o[id], c);
    }
}

// masked versions

template<typename T>
template<MATH_OPS op>
void CMemOpers<T>::
binaryC(T* h_o, const T* h_i, const T& c, const T* h_mask,
	size_t n, StreamT stream, bool onDev)
{
    if(!MOpers<T, op>::valid()){
	throw PyCAException(__FILE__, __LINE__, "Unimplemented math operation for specified type");
    }
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
	if(h_mask[id])
	    h_o[id] = MOpers<T, op>::op(h_i[id], c);
    }
}

template<typename T>
template<MATH_OPS op>
void CMemOpers<T>::
binaryC_I(T* h_o, const T& c, const T* h_mask,
	  size_t n, StreamT stream, bool onDev)
{
    if(!MOpers<T, op>::valid()){
	throw PyCAException(__FILE__, __LINE__, "Unimplemented math operation for specified type");
    }
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
	if(h_mask[id])
	    MOpers<T, op>::iop(h_o[id], c);
    }
}

/////////////////////////////////////////////////////////////////////////////////
// binary function with images
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
template<MATH_OPS op>
void CMemOpers<T>::binary(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT stream){
   if(!MOpers<T, op>::valid()){
      throw PyCAException(__FILE__, __LINE__, "Unimplemented math operation for specified type");
   }
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = MOpers<T, op>::op(h_i[id], h_i1[id]);
    }
}

template<typename T>
template<MATH_OPS op>
void CMemOpers<T>::binary_I(T* h_o, const T* h_i, size_t n, StreamT stream){
   if(!MOpers<T, op>::valid()){
      throw PyCAException(__FILE__, __LINE__, "Unimplemented math operation for specified type");
   }
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        MOpers<T, op>::iop(h_o[id], h_i[id]);
    }

}

// masked versions

template<typename T>
template<MATH_OPS op>
void CMemOpers<T>::
binary(T* h_o, const T* h_i, const T* h_i1, const T* h_mask,
       size_t n, StreamT stream)
{
   if(!MOpers<T, op>::valid()){
      throw PyCAException(__FILE__, __LINE__, "Unimplemented math operation for specified type");
   }
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
	if(h_mask[id])
	    h_o[id] = MOpers<T, op>::op(h_i[id], h_i1[id]);
    }
}

template<typename T>
template<MATH_OPS op>
void CMemOpers<T>::
binary_I(T* h_o, const T* h_i, const T* h_mask,
	 size_t n, StreamT stream)
{
   if(!MOpers<T, op>::valid()){
      throw PyCAException(__FILE__, __LINE__, "Unimplemented math operation for specified type");
   }
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
	if(h_mask[id])
	    MOpers<T, op>::iop(h_o[id], h_i[id]);
    }

}

CMEM_BINARY_OPERS(Add)
CMEM_BINARY_OPERS(Sub)
CMEM_BINARY_OPERS(Mul)
CMEM_BINARY_OPERS(Div)
CMEM_BINARY_OPERS(Max)
CMEM_BINARY_OPERS(Min)
CMEM_BINARY_OPERS(LT)
CMEM_BINARY_OPERS(LTE)
CMEM_BINARY_OPERS(EQ)
CMEM_BINARY_OPERS(NEQ)
CMEM_BINARY_OPERS(GT)
CMEM_BINARY_OPERS(GTE)
CMEM_BINARY_OPERS(Atan2)

// Pow doesn't have image version
CMEM_BINARY_OPERC(Pow)
CMEM_BINARY_OPERC_I(Pow)
CMEM_BINARY_OPER_NOC(Pow)
CMEM_BINARY_OPER_NOC_I(Pow)

/////////////////////////////////////////////////////////////////////////////////
// Absolute value of the difference
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::AbsDiff(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT stream){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = (h_i[id] >= h_i1[id]) ? (h_i[id] - h_i1[id]) : (h_i1[id] - h_i[id]);
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Absolute value of the difference
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::AbsDiff_I(T* h_o, const T* h_i, size_t n, StreamT stream){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = (h_o[id] >= h_i[id]) ? h_o[id] - h_i[id] : h_i[id] - h_o[id];
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Square value of the difference
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::SqrDiff(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT stream){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = (h_i[id] - h_i1[id]) * (h_i[id] - h_i1[id]);
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Square value of the difference
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::SqrDiff_I(T* h_o, const T* h_i, size_t n, StreamT stream){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = (h_o[id] - h_i[id]) * (h_o[id] - h_i[id]);
    }
}

/////////////////////////////////////////////////////////////////////////////////
// h_o = h_i * (h_i1 * c)
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::MulMulC(T* h_o, const T* h_i, const T* h_i1, const T& c, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = h_i[id] * (h_i1[id] * c);
    }
}

template<typename T>
void CMemOpers<T>::MulMulC_I(T* h_o, const T* h_i, const T& c, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] *= (h_i[id] * c);
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Multiply three arrays together
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::MulMul(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, size_t n, StreamT stream){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = h_i[id] * (h_i1[id] * h_i2[id]);
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Multiply three arrays together inplace version
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::MulMul_I(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT stream){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] *=  (h_i[id] * h_i1[id]);
    }
}


/////////////////////////////////////////////////////////////////////////////////
// Add two array together and divide by the third array
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::AddDiv(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, size_t n, StreamT stream){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = (h_i[id] + h_i1[id]) / h_i2[id];
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Add two array together and divide by the third array inplace
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::AddDiv_I(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT stream){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = (h_o[id] + h_i[id]) / h_i1[id];
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Sub two array together and divide by the third array
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::SubDiv(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, size_t n, StreamT stream){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = (h_i[id] - h_i1[id]) / h_i2[id];
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Sub two array together and divide by the third array inplace
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::SubDiv_I(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT stream){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = (h_o[id] - h_i[id]) / h_i1[id];
    }
}


/////////////////////////////////////////////////////////////////////////////////
// Multiply two array together and add the third one
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::MulAdd(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, size_t n, StreamT stream){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = h_i[id] * h_i1[id] + h_i2[id];
    }
}


/////////////////////////////////////////////////////////////////////////////////
// Multiply two array together and add the third one inplace version
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::MulAdd_I(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT stream){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = h_o[id] * h_i[id] + h_i1[id];
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Multiply two array together and sub the third one
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::MulSub(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, size_t n, StreamT stream){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = h_i[id] * h_i1[id] - h_i2[id];
    }
}


/////////////////////////////////////////////////////////////////////////////////
// Multiply two array together and sub the third one inplace version
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::MulSub_I(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT stream){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = h_o[id] * h_i[id] - h_i1[id];
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Add two arrays and multiply by the third one
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::AddMul(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, size_t n, StreamT stream){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = (h_i[id] + h_i1[id]) * h_i2[id];
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Add two array together and multiply by the third one inplace version
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::AddMul_I(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT stream){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = (h_o[id] + h_i[id]) * h_i1[id];
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Sub two arrays and multiply by the third one
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::SubMul(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, size_t n, StreamT stream){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = (h_i[id] - h_i1[id]) * h_i2[id];
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Sub two array together and multiply by the third one inplace version
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::SubMul_I(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT stream){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = (h_o[id] - h_i[id]) * h_i1[id];
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Add an array by a constant then multiply by other constant
// Used with normalized function
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::AddCMulC(T* h_o, const T* h_i, const T& a, const T& b, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = (h_i[id] + a) * b;
    }
}

template<typename T>
void CMemOpers<T>::AddCMulC_I(T* h_o, const T& a, const T& b, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = (h_o[id] + a) * b;
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Multiply an array by a constant and add the second one
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::MulCAdd(T* h_o, const T* h_i, const T& c, const T* h_i1, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = h_i[id] * c + h_i1[id];
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Multiply an array by a constant and add the second one inplace version
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::MulCAdd_I(T* h_o, const T& c, const T* h_i, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = h_o[id] * c + h_i[id];
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Multiply an array by a constant and sub the second one
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::MulCSub(T* h_o, const T* h_i, const T& c, const T* h_i1, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = h_i[id] * c - h_i1[id];
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Multiply an array by a constant and sub the second one inplace version
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::MulCSub_I(T* h_o, const T& c, const T* h_i, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = h_o[id] * c - h_i[id];
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Multiply an array by a constant and add the second constant
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::MulCAddC(T* h_o, const T* h_i, const T& a, const T& b, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = h_i[id] * a + b;
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Multiply an array by a constant and add the second constant inplace version
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::MulCAddC_I(T* h_o, const T& a, const T& b, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = h_o[id] * a + b;
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Add two array then multiply by a constant
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::AddMulC(T* h_o, const T* h_i, const T* h_i1, const T& c, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = (h_i[id] + h_i1[id]) * c;
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Add two array then multiply by a constant inplace version
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::AddMulC_I(T* h_o, const T* h_i, const T& c, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = (h_o[id] + h_i[id]) * c;
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Sub two array then multiply by a constant
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::SubMulC(T* h_o, const T* h_i, const T* h_i1, const T& c, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = (h_i[id] - h_i1[id]) * c;
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Sub two array then multiply by a constant inplace version
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::SubMulC_I(T* h_o, const T* h_i, const T& c, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = (h_o[id] - h_i[id]) * c;
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Add an arrays with another array that multiply by a constant
//  h_o = h_i + h_i1 * c
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::Add_MulC(T* h_o, const T* h_i, const T* h_i1, const T& c, size_t n, StreamT stream, bool onDev) {
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = h_i[id] + h_i1[id] * c;
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Add an arrays with another array that multiply by a constant
//  h_o = h_i + h_i1 * c
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::Add_MulC_I(T* h_o, const T* h_i, const T& c, size_t n, StreamT stream, bool onDev) {
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] += h_i[id] * c;
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Add an arrays with another array that multiply by an array
//  h_o = h_i + h_i1 * h_i2
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::Add_Mul(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, size_t n, StreamT stream) {
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = h_i[id] + h_i1[id] * h_i2[id];
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Add an arrays with another array that multiply by an array
//  h_o = h_o + h_i * h_i1
/////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::Add_Mul_I(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT stream) {
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] += h_i[id] * h_i1[id];
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Sub an arrays with another array that multiply by an array
//  h_o = h_i + h_i1 * h_i2
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::Sub_Mul(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, size_t n, StreamT stream) {
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = h_i[id] - h_i1[id] * h_i2[id];
    }
}

/////////////////////////////////////////////////////////////////////////////////
// Sub an arrays with another array that multiply by an array
//  h_o = h_o + h_i * h_i1
/////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::Sub_Mul_I(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT stream) {
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] -= h_i[id] * h_i1[id];
    }
}


/////////////////////////////////////////////////////////////////////////////////
// Quadary function
//////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::
MulC_Add_MulC(T* h_o, const T* h_i, const T& a, const T* h_i1, const T& b, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = h_i[id] * a + h_i1[id] * b;
    }
}

template<typename T>
void CMemOpers<T>::
MulC_Add_MulC_I(T* h_o, const T& a, const T* h_i, const T& b, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = h_o[id] * a + h_i[id] * b;
    }
}

////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::Add_AddMulC(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, const T& c, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = h_i[id] + (h_i1[id] + h_i2[id]) * c;
    }
}

////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::Add_AddMulC_I(T* h_o, const T* h_i, const T* h_i1, const T& c, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] += (h_i[id] + h_i1[id]) * c;
    }
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::Add_SubMulC(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, const T& c, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = h_i[id] + (h_i1[id] - h_i2[id]) * c;
    }
}

////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::Add_SubMulC_I(T* h_o, const T* h_i, const T* h_i1, const T& c, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] += (h_i[id] - h_i1[id]) * c;
    }
}

////////////////////////////////////////////////////////////////////////////////
//// h_o = h_i + (h_i1 * h_i2 * c)
////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::Add_MulMulC(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, const T& c, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = h_i[id] + h_i1[id] * (h_i2[id] * c);
    }
}

////////////////////////////////////////////////////////////////////////////////
// h_o = h_o + (h_i * h_i1 * d)
////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::Add_MulMulC_I(T* h_o, const T* h_i, const T* h_i1, const T& c, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] += h_i[id] * (h_i1[id] * c);
    }
}

////////////////////////////////////////////////////////////////////////////////
// h_o = h_i - (h_i1 + h_i2) *  c
////////////////////////////////////////////////////////////////////////////////

template<typename T>
void CMemOpers<T>::Sub_AddMulC(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, const T& c, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = h_i[id] - (h_i1[id] + h_i2[id]) * c;
    }
}

////////////////////////////////////////////////////////////////////////////////
// h_o -= (h_i + h_i1) * c
////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::Sub_AddMulC_I(T* h_o, const T* h_i, const T* h_i1, const T& c, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] -= (h_i[id] + h_i1[id]) * c;
    }
}
////////////////////////////////////////////////////////////////////////////////
// h_o = h_i - (h_i1 - h_i2) *  c
////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::Sub_SubMulC(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, const T& c, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = h_i[id] - (h_i1[id] - h_i2[id]) * c;
    }
}

////////////////////////////////////////////////////////////////////////////////
// h_o -= (h_i - h_i1) * c
////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::Sub_SubMulC_I(T* h_o, const T* h_i, const T* h_i1, const T& c, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] -= (h_i[id] - h_i1[id]) * c;
    }
}

////////////////////////////////////////////////////////////////////////////////
//// h_o = h_i + (h_i1 * h_i2 * c)
////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::Sub_MulMulC(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, const T& c, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = h_i[id] - h_i1[id] * (h_i2[id] * c);
    }
}

////////////////////////////////////////////////////////////////////////////////
// h_o = h_o + (h_i * h_i1 * c
////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::Sub_MulMulC_I(T* h_o, const T* h_i, const T* h_i1, const T& c, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] -= h_i[id] * (h_i1[id] * c);
    }
}

////////////////////////////////////////////////////////////////////////////////
// h_o = (h_i + a) * b + c
////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::AddCMulCAddC(T* h_o, const T* h_i, const T& a, const T& b, const T& c, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = (h_i[id] + a) * b + c;
    }
}

////////////////////////////////////////////////////////////////////////////////
// h_o = (h_o + a) * b + c
////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CMemOpers<T>::AddCMulCAddC_I(T* h_o, const T& a, const T& b, const T& c, size_t n, StreamT stream, bool onDev){
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = (h_o[id] + a) * b + c;
    }
}

template<typename T>
void CMemOpers<T>::ReverseOrder(T* h_o, const T* h_i, size_t n, StreamT st)
{
#pragma omp parallel for
    for (size_t id=0; id< n; ++id) {
        h_o[id] = h_i[n - 1 - id];
    }
}

template<typename T>
void CMemOpers<T>::ShiftCoordinate(T* a_o, const T* a_i, size_t sizeX, size_t sizeY, size_t sizeZ, bool dir, StreamT s)
{
    size_t id = 0;
    if (dir)
        for (size_t z=0; z< sizeZ; ++z)
            for (size_t y=0; y< sizeY; ++y)
                for (size_t x=0; x< sizeX; ++x, ++id)
                    a_o[ z + (x + y * sizeX) * sizeZ] = a_i[id];
    else
        for (size_t z=0; z< sizeZ; ++z)
            for (size_t y=0; y< sizeY; ++y)
                for (size_t x=0; x< sizeX; ++x, ++id)
                    a_o[ y + (z + x * sizeZ) * sizeY] = a_i[id];
}

template class CMemOpers<float>;
template class CMemOpers<int>;
template class CMemOpers<uint>;
} // end namespace PyCA

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

#include <ImageOpers.h>

#include "GImageOpers.h"
#include "CImageOpers.h"

#include <Image3D.h>
#include <MemOpers.h>

#define IMAGE_UNARY_ARRAY_OP_DEF(OP)                                       \
template<int mode>							   \
void ImageOpers<mode>							   \
::OP(Image3D& a_o, const Image3D& a_i, StreamT stream)			   \
{									   \
    PYCA_ASSERT(a_o.grid() == a_i.grid());				   \
    size_t nVox   = a_o.nVox();						   \
    MemOpers<mode,float>::OP(a_o.get(), a_i.get(), nVox, stream);	   \
}

#define IMAGE_UNARY_ARRAY_MASKED_OP_DEF(OP)                                \
template<int mode>							   \
void ImageOpers<mode>							   \
::OP(Image3D& a_o, const Image3D& a_i, const Image3D& a_mask, 		   \
     StreamT stream)							   \
{									   \
    PYCA_ASSERT(a_o.grid() == a_i.grid());				   \
    size_t nVox   = a_o.nVox();						   \
    MemOpers<mode,float>::OP(a_o.get(), a_i.get(), a_mask.get(),           \
			     nVox, stream);				   \
}

#define IMAGE_UNARY_ARRAY_OP_I_DEF(OP)					   \
template<int mode>							   \
void ImageOpers<mode>							   \
::OP##_I(Image3D& a_o, StreamT stream)					   \
{									   \
    size_t nVox   = a_o.nVox();						   \
    MemOpers<mode,float>::OP##_I(a_o.get(), nVox, stream);		   \
}

#define IMAGE_UNARY_ARRAY_MASKED_OP_I_DEF(OP)				   \
template<int mode>							   \
void ImageOpers<mode>							   \
::OP##_I(Image3D& a_o, const Image3D& a_mask, StreamT stream)		   \
{									   \
    size_t nVox   = a_o.nVox();						   \
    MemOpers<mode,float>::OP##_I(a_o.get(), a_mask.get(),		   \
				 nVox, stream);				   \
}

#define IMAGE_UNARY_ARRAY_OP_DEFS(OP)					   \
IMAGE_UNARY_ARRAY_OP_DEF(OP)						   \
IMAGE_UNARY_ARRAY_MASKED_OP_DEF(OP)					   \
IMAGE_UNARY_ARRAY_OP_I_DEF(OP)                                             \
IMAGE_UNARY_ARRAY_MASKED_OP_I_DEF(OP)

#define IMAGE_BINARY_ARRAY_OP_DEF(OP)  	       	       	       	       	   \
template<int mode>							   \
void ImageOpers<mode>							   \
::OP(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1,		   \
     StreamT stream)							   \
{									   \
    MK_CHECK3_SIZE(a_o, a_i, a_i1);					   \
    size_t nVox   = a_o.nVox();						   \
    MemOpers<mode,float>::OP(a_o.get(), a_i.get(), a_i1.get(), 		   \
			     nVox, stream);				   \
}

#define IMAGE_BINARY_ARRAY_MASKED_OP_DEF(OP)				   \
template<int mode>							   \
void ImageOpers<mode>							   \
::OP(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, 		   \
     const Image3D& a_mask, 						   \
     StreamT stream)							   \
{									   \
    MK_CHECK3_SIZE(a_o, a_i, a_i1);					   \
    size_t nVox   = a_o.nVox();						   \
    MemOpers<mode,float>::OP(a_o.get(), a_i.get(), a_i1.get(),		   \
			     a_mask.get(),				   \
			     nVox, stream);				   \
}

#define IMAGE_BINARY_ARRAY_OP_I_DEF(OP)					   \
template<int mode>							   \
void ImageOpers<mode>							   \
::OP##_I(Image3D& a_o, const Image3D& a_i, StreamT stream)		   \
{									   \
    MK_CHECK2_SIZE(a_o, a_i);						   \
    size_t nVox   = a_o.nVox();						   \
    MemOpers<mode,float>::OP##_I(a_o.get(), a_i.get(), 			   \
				 nVox, stream);				   \
}

#define IMAGE_BINARY_ARRAY_MASKED_OP_I_DEF(OP)				   \
template<int mode>							   \
void ImageOpers<mode>							   \
::OP##_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_mask, 	   \
	 StreamT stream)						   \
{									   \
    MK_CHECK2_SIZE(a_o, a_i);						   \
    size_t nVox   = a_o.nVox();						   \
    MemOpers<mode,float>::OP##_I(a_o.get(), a_i.get(), a_mask.get(),	   \
				 nVox, stream);				   \
}

#define IMAGE_BINARY_ARRAY_OPC_DEF(OP)					   \
template<int mode>							   \
void ImageOpers<mode>							   \
::OP##C(Image3D& a_o, const Image3D& a_i, const float& c,		   \
       StreamT stream, bool onDev)					   \
{									   \
    MK_CHECK2_SIZE(a_o, a_i);						   \
    size_t nVox   = a_o.nVox();						   \
    MemOpers<mode,float>::OP##C(a_o.get(), a_i.get(), c, 		   \
				nVox, stream, onDev);			   \
}

#define IMAGE_BINARY_ARRAY_MASKED_OPC_DEF(OP)				   \
template<int mode>							   \
void ImageOpers<mode>							   \
::OP##C(Image3D& a_o, const Image3D& a_i, const float& c,		   \
	const Image3D& a_mask, 						   \
	StreamT stream, bool onDev)					   \
{									   \
    MK_CHECK2_SIZE(a_o, a_i);						   \
    size_t nVox   = a_o.nVox();						   \
    MemOpers<mode,float>::OP##C(a_o.get(), a_i.get(), c,		   \
				a_mask.get(),				   \
				nVox, stream, onDev);			   \
}

#define IMAGE_BINARY_ARRAY_OPC_I_DEF(OP)       	       	       	       	   \
template<int mode>							   \
void ImageOpers<mode>							   \
::OP##C_I(Image3D& a_o, const float& c,					   \
	 StreamT stream, bool onDev)					   \
{									   \
    size_t nVox   = a_o.nVox();						   \
    MemOpers<mode,float>::OP##C_I(a_o.get(), c, 			   \
				  nVox, stream, onDev);			   \
}

#define IMAGE_BINARY_ARRAY_MASKED_OPC_I_DEF(OP)				   \
template<int mode>							   \
void ImageOpers<mode>							   \
::OP##C_I(Image3D& a_o, const float& c, const Image3D& a_mask, 		   \
	 StreamT stream, bool onDev)					   \
{									   \
    size_t nVox   = a_o.nVox();						   \
    MemOpers<mode,float>::OP##C_I(a_o.get(), c, a_mask.get(),		   \
				  nVox, stream, onDev);			   \
}

#define IMAGE_BINARY_ARRAY_OPC_DEF_NOC(OP)				   \
template<int mode>							   \
void ImageOpers<mode>							   \
::OP(Image3D& a_o, const Image3D& a_i, const float& c,		           \
       StreamT stream, bool onDev)					   \
{									   \
    MK_CHECK2_SIZE(a_o, a_i);						   \
    size_t nVox   = a_o.nVox();						   \
    MemOpers<mode,float>::OP(a_o.get(), a_i.get(), c, 		           \
				nVox, stream, onDev);			   \
}

#define IMAGE_BINARY_ARRAY_MASKED_OPC_DEF_NOC(OP)			   \
template<int mode>							   \
void ImageOpers<mode>							   \
::OP(Image3D& a_o, const Image3D& a_i, const float& c,			   \
     const Image3D& a_mask, 						   \
     StreamT stream, bool onDev)					   \
{									   \
    MK_CHECK2_SIZE(a_o, a_i);						   \
    size_t nVox   = a_o.nVox();						   \
    MemOpers<mode,float>::OP(a_o.get(), a_i.get(), c, a_mask.get(),	   \
			     nVox, stream, onDev);			   \
}

#define IMAGE_BINARY_ARRAY_OPC_I_DEF_NOC(OP)       	       	       	   \
template<int mode>							   \
void ImageOpers<mode>							   \
::OP##_I(Image3D& a_o, const float& c,					   \
	 StreamT stream, bool onDev)					   \
{									   \
    size_t nVox   = a_o.nVox();						   \
    MemOpers<mode,float>::OP##_I(a_o.get(), c, 			           \
				 nVox, stream, onDev);			   \
}

#define IMAGE_BINARY_ARRAY_MASKED_OPC_I_DEF_NOC(OP)			   \
template<int mode>							   \
void ImageOpers<mode>							   \
::OP##_I(Image3D& a_o, const float& c,					   \
	 const Image3D& a_mask, 					   \
	 StreamT stream, bool onDev)					   \
{									   \
    size_t nVox   = a_o.nVox();						   \
    MemOpers<mode,float>::OP##_I(a_o.get(), c, a_mask.get(),		   \
				 nVox, stream, onDev);			   \
}

#define IMAGE_BINARY_ARRAY_OP_DEFS(OP)					   \
IMAGE_BINARY_ARRAY_OP_DEF(OP)						   \
IMAGE_BINARY_ARRAY_MASKED_OP_DEF(OP)					   \
IMAGE_BINARY_ARRAY_OP_I_DEF(OP)						   \
IMAGE_BINARY_ARRAY_MASKED_OP_I_DEF(OP)					   \
IMAGE_BINARY_ARRAY_OPC_DEF(OP)						   \
IMAGE_BINARY_ARRAY_MASKED_OPC_DEF(OP)					   \
IMAGE_BINARY_ARRAY_OPC_I_DEF(OP)                                           \
IMAGE_BINARY_ARRAY_MASKED_OPC_I_DEF(OP)                                    \
IMAGE_BINARY_ARRAY_OPC_DEF_NOC(OP)                 			   \
IMAGE_BINARY_ARRAY_MASKED_OPC_DEF_NOC(OP)                 		   \
IMAGE_BINARY_ARRAY_OPC_I_DEF_NOC(OP)                                       \
IMAGE_BINARY_ARRAY_MASKED_OPC_I_DEF_NOC(OP)

namespace PyCA {

template<int mode>
void ImageOpers<mode>
::Copy(Image3D& a_o, const Image3D& a_i, const Image3D& a_mask, 
       StreamT stream)
{
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::Copy(a_o.get(), a_i.get(), a_mask.get(), 
			       nVox, stream);
}

template<int mode>
void ImageOpers<mode>
::SetMem(Image3D& a_o, const float& c, StreamT stream, bool onDev){
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::SetMem(a_o.get(), c, nVox, stream, onDev);
}

template<int mode>
void ImageOpers<mode>
::SetMem(Image3D& a_o, const float& c, const Image3D& a_mask, 
	 StreamT stream, bool onDev){
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::SetMem(a_o.get(), c, a_mask.get(), 
				 nVox, stream, onDev);
}

IMAGE_UNARY_ARRAY_OP_DEFS(Abs)
IMAGE_UNARY_ARRAY_OP_DEFS(Cube)
IMAGE_UNARY_ARRAY_OP_DEFS(Exp)
IMAGE_UNARY_ARRAY_OP_DEFS(Log)
IMAGE_UNARY_ARRAY_OP_DEFS(Neg)
IMAGE_UNARY_ARRAY_OP_DEFS(Ramp)
IMAGE_UNARY_ARRAY_OP_DEFS(Sgn)
IMAGE_UNARY_ARRAY_OP_DEFS(Step)
IMAGE_UNARY_ARRAY_OP_DEFS(Sqr)
IMAGE_UNARY_ARRAY_OP_DEFS(Sqrt)
IMAGE_UNARY_ARRAY_OP_DEFS(Inv)
IMAGE_UNARY_ARRAY_OP_DEFS(Ceil)
IMAGE_UNARY_ARRAY_OP_DEFS(Floor)
IMAGE_UNARY_ARRAY_OP_DEFS(Round)
IMAGE_UNARY_ARRAY_OP_DEFS(Sin)
IMAGE_UNARY_ARRAY_OP_DEFS(Asin)
IMAGE_UNARY_ARRAY_OP_DEFS(Cos)
IMAGE_UNARY_ARRAY_OP_DEFS(Acos)
IMAGE_UNARY_ARRAY_OP_DEFS(Tan)
IMAGE_UNARY_ARRAY_OP_DEFS(Atan)
IMAGE_UNARY_ARRAY_OP_DEFS(Csc)
IMAGE_UNARY_ARRAY_OP_DEFS(Sec)
IMAGE_UNARY_ARRAY_OP_DEFS(Cot)

///

IMAGE_BINARY_ARRAY_OP_DEFS(Add)
IMAGE_BINARY_ARRAY_OP_DEFS(Sub)
IMAGE_BINARY_ARRAY_OP_DEFS(Mul)
IMAGE_BINARY_ARRAY_OP_DEFS(Div)
IMAGE_BINARY_ARRAY_OP_DEFS(Min)
IMAGE_BINARY_ARRAY_OP_DEFS(Max)
IMAGE_BINARY_ARRAY_OP_DEFS(Atan2)
IMAGE_BINARY_ARRAY_OP_DEFS(GT)
IMAGE_BINARY_ARRAY_OP_DEFS(GTE)
IMAGE_BINARY_ARRAY_OP_DEFS(EQ)
IMAGE_BINARY_ARRAY_OP_DEFS(NEQ)
IMAGE_BINARY_ARRAY_OP_DEFS(LT)
IMAGE_BINARY_ARRAY_OP_DEFS(LTE)




template<int mode>
void ImageOpers<mode>
::PowC(Image3D& a_o, const Image3D& a_i, const float& c, StreamT stream, bool onDev){
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::PowC(a_o.get(), a_i.get(), c, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::PowC_I(Image3D& a_o, const float& c, StreamT stream, bool onDev){
   size_t nVox   = a_o.nVox();
   MemOpers<mode,float>::PowC_I(a_o.get(), c, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::AbsDiff(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::AbsDiff(a_o.get(), a_i.get(), a_i1.get(), nVox, stream);
}

template<int mode>
void ImageOpers<mode>
::SqrDiff(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::SqrDiff(a_o.get(), a_i.get(), a_i1.get(), nVox, stream);
}

template<int mode>
void ImageOpers<mode>
::AbsDiff_I(Image3D& a_o, const Image3D& a_i, StreamT stream){
    MK_CHECK2_SIZE(a_o, a_i);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::AbsDiff_I(a_o.get(), a_i.get(), nVox, stream);
}


template<int mode>
void ImageOpers<mode>
::SqrDiff_I(Image3D& a_o, const Image3D& a_i, StreamT stream){
    MK_CHECK2_SIZE(a_o, a_i);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::SqrDiff_I(a_o.get(), a_i.get(), nVox, stream);
}

///

template<int mode>
void ImageOpers<mode>
::MulMulC(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const float& c, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::MulMulC(a_o.get(), a_i.get(), a_i1.get(), c, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::MulMulC_I(Image3D& a_o, const Image3D& a_i, const float& c, StreamT stream, bool onDev){
    MK_CHECK2_SIZE(a_o, a_i);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::MulMulC_I(a_o.get(), a_i.get(), c, nVox, stream,onDev);
}



template<int mode>
void ImageOpers<mode>
::MulMul(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, StreamT stream){
    MK_CHECK4_SIZE(a_o, a_i, a_i1, a_i2);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::MulMul(a_o.get(), a_i.get(), a_i1.get(), a_i2.get(), nVox, stream);
}

template<int mode>
void ImageOpers<mode>
::MulMul_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::MulMul_I(a_o.get(), a_i.get(), a_i1.get(), nVox, stream);
}

template<int mode>
void ImageOpers<mode>
::MulAdd(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, StreamT stream){
    MK_CHECK4_SIZE(a_o, a_i, a_i1, a_i2);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::MulAdd(a_o.get(), a_i.get(), a_i1.get(), a_i2.get(), nVox, stream);
}

template<int mode>
void ImageOpers<mode>
::MulAdd_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t nVox   = a_o.nVox();
    MemOpers<mode, float>::MulAdd_I(a_o.get(), a_i.get(), a_i1.get(), nVox, stream);
}

template<int mode>
void ImageOpers<mode>
::MulSub(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, StreamT stream){
    MK_CHECK4_SIZE(a_o, a_i, a_i1, a_i2);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::MulSub(a_o.get(), a_i.get(), a_i1.get(), a_i2.get(), nVox, stream);
}

template<int mode>
void ImageOpers<mode>
::MulSub_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t nVox   = a_o.nVox();
    MemOpers<mode, float>::MulSub_I(a_o.get(), a_i.get(), a_i1.get(), nVox, stream);
}

template<int mode>
void ImageOpers<mode>
::AddMul(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, StreamT stream){
    MK_CHECK4_SIZE(a_o, a_i, a_i1, a_i2);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::AddMul(a_o.get(), a_i.get(), a_i1.get(), a_i2.get(), nVox, stream);
}

template<int mode>
void ImageOpers<mode>
::AddMul_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::AddMul_I(a_o.get(), a_i.get(), a_i1.get(), nVox, stream);
}

template<int mode>
void ImageOpers<mode>
::AddDiv(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, StreamT stream){
    MK_CHECK4_SIZE(a_o, a_i, a_i1, a_i2);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::AddDiv(a_o.get(), a_i.get(), a_i1.get(), a_i2.get(), nVox, stream);
}

template<int mode>
void ImageOpers<mode>
::AddDiv_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::AddDiv_I(a_o.get(), a_i.get(), a_i1.get(), nVox, stream);
}

template<int mode>
void ImageOpers<mode>
::SubMul(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, StreamT stream){
    MK_CHECK4_SIZE(a_o, a_i, a_i1, a_i2);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::SubMul(a_o.get(), a_i.get(), a_i1.get(), a_i2.get(), nVox, stream);
}

template<int mode>
void ImageOpers<mode>
::SubMul_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::SubMul_I(a_o.get(), a_i.get(), a_i1.get(), nVox, stream);
}

template<int mode>
void ImageOpers<mode>
::SubDiv(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, StreamT stream){
    MK_CHECK4_SIZE(a_o, a_i, a_i1, a_i2);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::SubDiv(a_o.get(), a_i.get(), a_i1.get(), a_i2.get(), nVox, stream);
}

template<int mode>
void ImageOpers<mode>
::SubDiv_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::SubDiv_I(a_o.get(), a_i.get(), a_i1.get(), nVox, stream);
}

template<int mode>
void ImageOpers<mode>
::MulCAdd(Image3D& a_o, const Image3D& a_i, const float& c, const Image3D& a_i1, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::MulCAdd(a_o.get(), a_i.get(), c, a_i1.get(), nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::MulCAdd_I(Image3D& a_o, const float& c, const Image3D& a_i, StreamT stream, bool onDev){
    MK_CHECK2_SIZE(a_o, a_i);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::MulCAdd_I(a_o.get(), c, a_i.get(), nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::MulCSub(Image3D& a_o, const Image3D& a_i, const float& c, const Image3D& a_i1, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::MulCSub(a_o.get(), a_i.get(), c, a_i1.get(), nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::MulCSub_I(Image3D& a_o, const float& c, const Image3D& a_i, StreamT stream, bool onDev){
    MK_CHECK2_SIZE(a_o, a_i);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::MulCSub_I(a_o.get(), c, a_i.get(), nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::MulCAddC(Image3D& a_o, const Image3D& a_i, const float& a, const float& b, StreamT stream, bool onDev){
    MK_CHECK2_SIZE(a_o, a_i);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::MulCAddC(a_o.get(), a_i.get(), a, b, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::MulCAddC_I(Image3D& a_o, const float& a, const float& b, StreamT stream, bool onDev){
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::MulCAddC_I(a_o.get(),a, b, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::AddMulC(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const float& c, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::AddMulC(a_o.get(), a_i.get(), a_i1.get(), c, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::AddMulC_I(Image3D& a_o, const Image3D& a_i, const float& c, StreamT stream, bool onDev){
    MK_CHECK2_SIZE(a_o, a_i);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::AddMulC_I(a_o.get(), a_i.get(), c, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::SubMulC(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const float& c, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::SubMulC(a_o.get(), a_i.get(), a_i1.get(), c, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::SubMulC_I(Image3D& a_o, const Image3D& a_i, const float& c, StreamT stream, bool onDev){
    MK_CHECK2_SIZE(a_o, a_i);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::SubMulC_I(a_o.get(), a_i.get(), c, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::AddCMulC(Image3D& a_o, const Image3D& a_i, const float& a, const float& b, StreamT stream, bool onDev){
    MK_CHECK2_SIZE(a_o, a_i);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::AddCMulC(a_o.get(), a_i.get(), a, b, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::AddCMulC_I(Image3D& a_o, const float& a, const float& b, StreamT stream, bool onDev){
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::AddCMulC_I(a_o.get(), a, b, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::Add_MulC(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const float& c, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::Add_MulC(a_o.get(), a_i.get(), a_i1.get(), c, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::Add_MulC_I(Image3D& a_o, const Image3D& a_i, const float& c, StreamT stream, bool onDev){
    MK_CHECK2_SIZE(a_o, a_i);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::Add_MulC_I(a_o.get(), a_i.get(), c, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::Add_Mul(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, StreamT stream){
    MK_CHECK4_SIZE(a_o, a_i, a_i1, a_i2);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::Add_Mul(a_o.get(), a_i.get(), a_i1.get(), a_i2.get(), nVox, stream);
}

template<int mode>
void ImageOpers<mode>
::Add_Mul_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::Add_Mul_I(a_o.get(), a_i.get(), a_i1.get(), nVox, stream);
}

template<int mode>
void ImageOpers<mode>
::Sub_Mul(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, StreamT stream){
    MK_CHECK4_SIZE(a_o, a_i, a_i1, a_i2);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::Sub_Mul(a_o.get(), a_i.get(), a_i1.get(), a_i2.get(), nVox, stream);
}

template<int mode>
void ImageOpers<mode>
::Sub_Mul_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::Sub_Mul_I(a_o.get(), a_i.get(), a_i1.get(), nVox, stream);
}

template<int mode>
void ImageOpers<mode>
::Add_AddMulC(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, const float& c, StreamT stream, bool onDev){
    MK_CHECK4_SIZE(a_o, a_i, a_i1, a_i2);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::Add_AddMulC(a_o.get(), a_i.get(), a_i1.get(), a_i2.get(), c, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::Add_SubMulC(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, const float& c, StreamT stream, bool onDev){
    MK_CHECK4_SIZE(a_o, a_i, a_i1, a_i2);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::Add_SubMulC(a_o.get(), a_i.get(), a_i1.get(), a_i2.get(), c, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::Add_MulMulC(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, const float& c, StreamT stream, bool onDev){
    MK_CHECK4_SIZE(a_o, a_i, a_i1, a_i2);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::Add_MulMulC(a_o.get(), a_i.get(), a_i1.get(), a_i2.get(), c, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::Add_AddMulC_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const float& c, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::Add_AddMulC_I(a_o.get(), a_i.get(), a_i1.get(), c, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::Add_SubMulC_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const float& c, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::Add_SubMulC_I(a_o.get(), a_i.get(), a_i1.get(), c, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::Add_MulMulC_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const float& c, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::Add_MulMulC_I(a_o.get(), a_i.get(), a_i1.get(), c, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::Sub_AddMulC(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, const float& c, StreamT stream, bool onDev){
    MK_CHECK4_SIZE(a_o, a_i, a_i1, a_i2);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::Sub_AddMulC(a_o.get(), a_i.get(), a_i1.get(), a_i2.get(), c, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::Sub_SubMulC(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, const float& c, StreamT stream, bool onDev){
    MK_CHECK4_SIZE(a_o, a_i, a_i1, a_i2);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::Sub_SubMulC(a_o.get(), a_i.get(), a_i1.get(), a_i2.get(), c, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::Sub_MulMulC(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, const float& c, StreamT stream, bool onDev){
    MK_CHECK4_SIZE(a_o, a_i, a_i1, a_i2);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::Sub_MulMulC(a_o.get(), a_i.get(), a_i1.get(), a_i2.get(), c, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::Sub_AddMulC_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const float& c, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::Sub_AddMulC_I(a_o.get(), a_i.get(), a_i1.get(), c, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::Sub_SubMulC_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const float& c, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::Sub_SubMulC_I(a_o.get(), a_i.get(), a_i1.get(), c, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::Sub_MulMulC_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const float& c, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::Sub_MulMulC_I(a_o.get(), a_i.get(), a_i1.get(), c, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::MulC_Add_MulC(Image3D& a_o, const Image3D& a_i, const float& a, const Image3D& a_i1, const float& b, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::MulC_Add_MulC(a_o.get(), a_i.get(), a, a_i1.get(), b, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::MulC_Add_MulC_I(Image3D& a_o, const float& a, const Image3D& a_i, const float& b, StreamT stream, bool onDev){
    MK_CHECK2_SIZE(a_o, a_i);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::MulC_Add_MulC_I(a_o.get(), a, a_i.get(), b, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::AddCMulCAddC(Image3D& a_o, const Image3D& a_i, const float& a, const float& b, const float& c, StreamT stream, bool onDev){
    MK_CHECK2_SIZE(a_o, a_i);
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::AddCMulCAddC(a_o.get(), a_i.get(), a, b, c, nVox, stream,onDev);
}

template<int mode>
void ImageOpers<mode>
::AddCMulCAddC_I(Image3D& a_o, const float& a, const float& b, const float& c, StreamT stream, bool onDev){
    size_t nVox   = a_o.nVox();
    MemOpers<mode,float>::AddCMulCAddC_I(a_o.get(), a, b, c, nVox, stream,onDev);
}


template<int mode>
void ImageOpers<mode>
::ShiftCoordinate(Image3D& a_o, const Image3D& a_i, bool dir, StreamT s)
{
    Vec3Di oSize = a_o.size();
    Vec3Di iSize = a_i.size();

    Vec3Df oSp = a_o.spacing();
    Vec3Df iSp = a_o.spacing();

    if (dir == true ) {
       PRECONDITION((oSize.z = iSize.x) && (oSize.x = iSize.y) && (oSize.y = iSize.z), "Input/output sizes incompatible");
       PRECONDITION((oSp.z = iSp.x) && (oSp.x = iSp.y) && (oSp.y = iSp.z), "Input/output spacings incompatible");
    }
    else {
       PRECONDITION((oSize.y = iSize.x) && (oSize.z = iSize.y) && (oSize.x = iSize.z), "Input/output sizes incompatible");
       PRECONDITION((oSp.y = iSp.x) && (oSp.z = iSp.y) && (oSp.x = iSp.z), "Input/output spacings incompatible");
    }
    MemOpers<mode,float>::ShiftCoordinate(a_o.get(), a_i.get(),
                                           iSize.x, iSize.y, iSize.z, dir, s);

}

template<int mode>
void ImageOpers<mode>::
SubVol(Image3D& a_o,const Image3D& a_i,
       const Vec3Di& start, StreamT st)
{
   Executer::SubVol(a_o, a_i, start, st);
}

template<int mode>
void ImageOpers<mode>::
SetSubVol_I(Image3D& a_o,const Image3D& a_i,
	    const Vec3Di& start, StreamT st)
{
   Executer::SetSubVol_I(a_o, a_i, start, st);
}

template<int mode>
void ImageOpers<mode>::
Shrink_I(Image3D& a_o, float c, StreamT st)
{
   Executer::Shrink_I(a_o, c, st);
}

template<int mode>
void ImageOpers<mode>::
Shrink(Image3D& a_o,const Image3D& a_i,
       float c, StreamT st)
{
   Executer::Shrink(a_o, a_i, c, st);
}

template<int mode>
void ImageOpers<mode>::
SoftAbs_I(Image3D& a_o, float eps, StreamT st)
{
   Executer::SoftAbs_I(a_o, eps, st);
}

template<int mode>
void ImageOpers<mode>::
SoftAbs(Image3D& a_o,const Image3D& a_i,
       float eps, StreamT st)
{
   Executer::SoftAbs(a_o, a_i, eps, st);
}

template<int mode>
void ImageOpers<mode>::
SoftSgn_I(Image3D& a_o, float eps, StreamT st)
{
   Executer::SoftSgn_I(a_o, eps, st);
}

template<int mode>
void ImageOpers<mode>::
SoftSgn(Image3D& a_o,const Image3D& a_i,
       float eps, StreamT st)
{
   Executer::SoftSgn(a_o, a_i, eps, st);
}

template<int mode>
template<BackgroundStrategy bg, InterpT interp, bool useOriginOffset>
void ImageOpers<mode>::
Resample(Image3D& a_o, const Image3D& a_i, StreamT stream){
    Executer::template Resample<bg, interp, useOriginOffset>
	(a_o, a_i, stream);
}

template<int mode>
template<BackgroundStrategy bg, InterpT interp>
void ImageOpers<mode>::
ResampleWorld(Image3D& a_o, const Image3D& a_i, StreamT stream){
    Executer::template ResampleWorld<bg, interp>
	(a_o, a_i, stream);
}

template<int mode>
template<BackgroundStrategy bg>
void ImageOpers<mode>::
SplatWorld(Image3D& a_o, const Image3D& a_i, StreamT stream){
   Executer::template SplatWorld<bg>(a_o, a_i, stream);
}

template<int mode>
template<BackgroundStrategy bg>
void ImageOpers<mode>::
SplatWorld(Image3D& a_o, const Image3D& a_i, Image3D& a_w, StreamT stream){
   Executer::template SplatWorld<bg>(a_o, a_i, a_w, stream);
}

template<int mode>
void ImageOpers<mode>::
Convolve(Image3D& a_o, const Image3D& a_i,
	 const Image3D& kernel, StreamT stream)
{
   Executer::Convolve(a_o, a_i, kernel, stream);
}


// template instantiation

#include "ImageOpers_inst.cxx"

} // end namespace PyCA

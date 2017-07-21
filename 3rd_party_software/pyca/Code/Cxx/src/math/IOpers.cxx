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
									    
#include <IOpers.h>							    
									    
#include "GImageOpers.h"						    
#include "CImageOpers.h"						    
									    
#include "MemOpers.h"							    
									    
#include <ImageOpers.h>							    
#include <Image3D.h>							    
									    
#define IMAGE_UNARY_ARRAY_OP_DEF(OP)                                        \
void Opers								    \
::OP(Image3D& a_o, const Image3D& a_i, StreamT stream)			    \
{									    \
    MK_CHECK2_ALL(a_o, a_i);						    \
    AUTO_EXEC(a_o.memType(), ImageOpers,				    \
	      OP(a_o, a_i, stream));					    \
}				       	       	       	       	       	     
									    
#define IMAGE_UNARY_ARRAY_MASKED_OP_DEF(OP)	       	       	       	    \
void Opers								    \
::OP(Image3D& a_o, const Image3D& a_i, const Image3D& a_mask,		    \
     StreamT stream)							    \
{									    \
    MK_CHECK3_ALL(a_o, a_i, a_mask);					    \
    AUTO_EXEC(a_o.memType(), ImageOpers,				    \
	      OP(a_o, a_i, a_mask, stream));				    \
}									     
									    
#define IMAGE_UNARY_ARRAY_OP_I_DEF(OP)					    \
void Opers								    \
 ::OP##_I(Image3D& a_o, StreamT stream)					    \
{									    \
    AUTO_EXEC(a_o.memType(), ImageOpers,				    \
	      OP##_I(a_o, stream));					    \
}					       	       	       	       	      
									    
#define IMAGE_UNARY_ARRAY_MASKED_OP_I_DEF(OP)				    \
void Opers								    \
::OP##_I(Image3D& a_o, const Image3D& a_mask, StreamT stream)		    \
{									    \
    MK_CHECK2_ALL(a_o, a_mask);						    \
    AUTO_EXEC(a_o.memType(), ImageOpers,				    \
	      OP##_I(a_o, a_mask, stream));	       	       	       	    \
}									     
									    
#define IMAGE_UNARY_ARRAY_OP_DEFS(OP)					    \
IMAGE_UNARY_ARRAY_OP_DEF(OP)						    \
IMAGE_UNARY_ARRAY_MASKED_OP_DEF(OP)					    \
IMAGE_UNARY_ARRAY_OP_I_DEF(OP)						    \
IMAGE_UNARY_ARRAY_MASKED_OP_I_DEF(OP)   	       	       	       	       	     
				    					    
									    
#define IMAGE_BINARY_ARRAY_OP_DEF(OP)  	       	       	       	       	    \
void Opers								    \
::OP(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, 		    \
     StreamT stream)							    \
{									    \
   MK_CHECK3_ALL(a_o, a_i, a_i1);					    \
   AUTO_EXEC(a_o.memType(), ImageOpers,					    \
	     OP(a_o, a_i, a_i1, stream));				    \
}									     
									    
#define IMAGE_BINARY_ARRAY_MASKED_OP_DEF(OP)				    \
void Opers								    \
::OP(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1,		    \
     const Image3D& a_mask, 						    \
     StreamT stream)							    \
{									    \
    MK_CHECK4_ALL(a_o, a_i, a_i1, a_mask);				    \
    AUTO_EXEC(a_o.memType(), ImageOpers,				    \
	      OP(a_o, a_i, a_i1, a_mask, stream));			    \
}									     
								       	     
#define IMAGE_BINARY_ARRAY_OP_I_DEF(OP)                	       	       	    \
void Opers								    \
::OP##_I(Image3D& a_o, const Image3D& a_i, 				    \
	 StreamT stream)						    \
{									    \
   MK_CHECK2_ALL(a_o, a_i);						    \
   AUTO_EXEC(a_o.memType(), ImageOpers,					    \
	     OP##_I(a_o, a_i, stream));					    \
}									     
									    
#define IMAGE_BINARY_ARRAY_MASKED_OP_I_DEF(OP)				    \
void Opers								    \
::OP##_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_mask, 	    \
	 StreamT stream)						    \
{									    \
    MK_CHECK3_ALL(a_o, a_i, a_mask);					    \
    AUTO_EXEC(a_o.memType(), ImageOpers,				    \
	      OP##_I(a_o, a_i, a_mask, stream));			    \
}									     
									     
#define IMAGE_BINARY_ARRAY_OPC_DEF(OP)         	       	       	       	    \
void Opers								    \
::OP##C(Image3D& a_o, const Image3D& a_i, const float& c, 		    \
	StreamT stream, bool onDev)					    \
{									    \
   MK_CHECK2_ALL(a_o, a_i);						    \
   AUTO_EXEC(a_o.memType(), ImageOpers,					    \
	     OP##C(a_o, a_i, c, stream, onDev));			    \
}									     
									     
#define IMAGE_BINARY_ARRAY_MASKED_OPC_DEF(OP)				    \
void Opers								    \
::OP##C(Image3D& a_o, const Image3D& a_i, const float& c,		    \
	const Image3D& a_mask, 						    \
	StreamT stream, bool onDev)					    \
{									    \
    MK_CHECK3_ALL(a_o, a_i, a_mask);					    \
    AUTO_EXEC(a_o.memType(), ImageOpers,				    \
	      OP##C(a_o, a_i, c, a_mask, stream, onDev));		    \
}									     
									    
#define IMAGE_BINARY_ARRAY_OPC_I_DEF(OP)       	       	       	       	    \
void Opers								    \
::OP##C_I(Image3D& a_o, const float& c, 				    \
	  StreamT stream, bool onDev)					    \
{								    	    \
   AUTO_EXEC(a_o.memType(), ImageOpers,				    	    \
	     OP##C_I(a_o, c, stream, onDev));			    	    \
}									     
									    
#define IMAGE_BINARY_ARRAY_MASKED_OPC_I_DEF(OP)				    \
void Opers								    \
::OP##C_I(Image3D& a_o, const float& c, const Image3D& a_mask, 		    \
	  StreamT stream, bool onDev)					    \
{									    \
    MK_CHECK2_ALL(a_o, a_mask);						    \
    AUTO_EXEC(a_o.memType(), ImageOpers,				    \
	      OP##C_I(a_o, c, a_mask, stream, onDev));			    \
}									    
									    
#define IMAGE_BINARY_ARRAY_OP_NOC_DEF(OP)         	       	       	    \
void Opers							    	    \
::OP(Image3D& a_o, const Image3D& a_i, const float& c, 		       	    \
	StreamT stream, bool onDev)				    	    \
{								    	    \
   MK_CHECK2_ALL(a_o, a_i);					    	    \
   AUTO_EXEC(a_o.memType(), ImageOpers,				    	    \
	     OP(a_o, a_i, c, stream, onDev));			       	    \
}									     
									    
#define IMAGE_BINARY_ARRAY_MASKED_OP_NOC_DEF(OP)			    \
void Opers								    \
::OP(Image3D& a_o, const Image3D& a_i, const float& c,			    \
     const Image3D& a_mask, 						    \
     StreamT stream, bool onDev)					    \
{									    \
    MK_CHECK3_ALL(a_o, a_i, a_mask);					    \
    AUTO_EXEC(a_o.memType(), ImageOpers,				    \
	      OP(a_o, a_i, c, a_mask, stream, onDev));			    \
}									     
									    
#define IMAGE_BINARY_ARRAY_OP_NOC_I_DEF(OP)       	       	       	    \
void Opers							   	    \
::OP##_I(Image3D& a_o, const float& c, 				       	    \
	  StreamT stream, bool onDev)				   	    \
{								   	    \
   AUTO_EXEC(a_o.memType(), ImageOpers,				   	    \
	     OP##_I(a_o, c, stream, onDev));			   	    \
}									     
									    
#define IMAGE_BINARY_ARRAY_MASKED_OP_NOC_I_DEF(OP)			    \
void Opers								    \
::OP##_I(Image3D& a_o, const float& c, const Image3D& a_mask, 		    \
	  StreamT stream, bool onDev)					    \
{									    \
    MK_CHECK2_ALL(a_o, a_mask);						    \
    AUTO_EXEC(a_o.memType(), ImageOpers,				    \
	      OP##_I(a_o, c, a_mask, stream, onDev));			    \
}									     
									    
#define IMAGE_BINARY_ARRAY_OP_DEFS(OP)					    \
IMAGE_BINARY_ARRAY_OPC_DEF(OP)					            \
IMAGE_BINARY_ARRAY_MASKED_OPC_DEF(OP)					    \
IMAGE_BINARY_ARRAY_OP_NOC_DEF(OP)					    \
IMAGE_BINARY_ARRAY_MASKED_OP_NOC_DEF(OP)				    \
IMAGE_BINARY_ARRAY_OPC_I_DEF(OP)					    \
IMAGE_BINARY_ARRAY_MASKED_OPC_I_DEF(OP)					    \
IMAGE_BINARY_ARRAY_OP_NOC_I_DEF(OP)					    \
IMAGE_BINARY_ARRAY_MASKED_OP_NOC_I_DEF(OP)				    \
IMAGE_BINARY_ARRAY_OP_DEF(OP)					            \
IMAGE_BINARY_ARRAY_MASKED_OP_DEF(OP)					    \
IMAGE_BINARY_ARRAY_OP_I_DEF(OP)                                             \
IMAGE_BINARY_ARRAY_MASKED_OP_I_DEF(OP) 	       	       	       	       	     


namespace PyCA {

void Opers
::Copy(Image3D& a_o, const Image3D& a_i, StreamT stream){
    MK_CHECK2_SIZE(a_o, a_i);
    Copy(a_o.getMemPool(), a_i.getMemPool(), a_o.nVox(), stream);

    bool on_host   = ((a_o.memType()!= MEM_DEVICE) && (a_i.memType()!= MEM_DEVICE));
    bool on_device = ((a_o.memType()== MEM_DEVICE) && (a_i.memType()== MEM_DEVICE));
    MK_CHECK2_SIZE(a_o, a_i);

    if (stream != NULL ) {
	PRECONDITION(on_device ||
		     a_o.memType() == MEM_HOST_PINNED ||
		     a_i.memType() == MEM_HOST_PINNED,
		     "Pinned memory is required for streaming between CPU and GPU");
    }

    size_t nVox   = a_o.nVox();
    if (on_host){
	cpyArrayH2H(a_o.get(), a_i.get(), nVox);
    } else if (on_device){
	acpyArrayD2D(a_o.get(), a_i.get(), nVox, stream);
    } else {
	if (a_o.memType() == MEM_DEVICE) {
	    (stream == NULL) ? cpyArrayH2D(a_o.get(), a_i.get(), nVox) :
		acpyArrayH2D(a_o.get(), a_i.get(), nVox, stream);
	}
	else {
	    (stream == NULL) ? cpyArrayD2H(a_o.get(), a_i.get(), nVox) :
		acpyArrayD2H(a_o.get(), a_i.get(), nVox, stream);
	}
    }
}

void Opers
::Copy(Image3D& a_o, const Image3D& a_i, const Image3D& a_mask,
       StreamT stream)
{
   AUTO_EXEC(a_o.memType(), ImageOpers,
	     Copy(a_o, a_i, a_mask, stream));
}

void Opers
::SetMem(Image3D& a_o, const float& c, StreamT stream, bool onDev){
   AUTO_EXEC(a_o.memType(), ImageOpers,
	     SetMem(a_o, c, stream, onDev));
}

void Opers
::SetMem(Image3D& a_o, const float& c, const Image3D& a_mask, 
	 StreamT stream, bool onDev){
   AUTO_EXEC(a_o.memType(), ImageOpers,
	     SetMem(a_o, c, a_mask, stream, onDev));
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

IMAGE_BINARY_ARRAY_OP_DEFS(Add)
IMAGE_BINARY_ARRAY_OP_DEFS(Sub)
IMAGE_BINARY_ARRAY_OP_DEFS(Mul)
IMAGE_BINARY_ARRAY_OP_DEFS(Div)
IMAGE_BINARY_ARRAY_OP_DEFS(Max)
IMAGE_BINARY_ARRAY_OP_DEFS(Min)
IMAGE_BINARY_ARRAY_OP_DEFS(Atan2)
IMAGE_BINARY_ARRAY_OP_DEFS(GT)
IMAGE_BINARY_ARRAY_OP_DEFS(GTE)
IMAGE_BINARY_ARRAY_OP_DEFS(EQ)
IMAGE_BINARY_ARRAY_OP_DEFS(NEQ)
IMAGE_BINARY_ARRAY_OP_DEFS(LT)
IMAGE_BINARY_ARRAY_OP_DEFS(LTE)

void Opers
::PowC(Image3D& a_o, const Image3D& a_i, const float& c, StreamT stream, bool onDev){
    MK_CHECK2_ALL(a_o, a_i);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      PowC(a_o, a_i, c, stream, onDev));
}

void Opers
::PowC_I(Image3D& a_o, const float& c, StreamT stream, bool onDev){
   AUTO_EXEC(a_o.memType(), ImageOpers,
	     PowC_I(a_o, c, stream, onDev));
}


void Opers
::AbsDiff(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_ALL(a_o, a_i, a_i1);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      AbsDiff(a_o, a_i, a_i1, stream));
}

void Opers
::SqrDiff(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_ALL(a_o, a_i, a_i1);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      SqrDiff(a_o, a_i, a_i1, stream));
}


void Opers
::AbsDiff_I(Image3D& a_o, const Image3D& a_i, StreamT stream){
    MK_CHECK2_ALL(a_o, a_i);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      AbsDiff_I(a_o, a_i, stream));
}


void Opers
::SqrDiff_I(Image3D& a_o, const Image3D& a_i, StreamT stream){
    MK_CHECK2_ALL(a_o, a_i);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      SqrDiff_I(a_o, a_i, stream));
}

///

void Opers
::MulMulC(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const float& c, StreamT stream, bool onDev){
    MK_CHECK3_ALL(a_o, a_i, a_i1);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      MulMulC(a_o, a_i, a_i1, c, stream, onDev));
}



void Opers
::MulMulC_I(Image3D& a_o, const Image3D& a_i, const float& c, StreamT stream, bool onDev){
    MK_CHECK2_ALL(a_o, a_i);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      MulMulC_I(a_o, a_i, c, stream, onDev));
}



void Opers
::MulMul(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, StreamT stream){
    MK_CHECK4_ALL(a_o, a_i, a_i1, a_i2);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      MulMul(a_o, a_i, a_i1, a_i2, stream));
}


void Opers
::MulMul_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_ALL(a_o, a_i, a_i1);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      MulMul_I(a_o, a_i, a_i1, stream));
}


void Opers
::MulAdd(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, StreamT stream){
    MK_CHECK4_ALL(a_o, a_i, a_i1, a_i2);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      MulAdd(a_o, a_i, a_i1, a_i2, stream));
}


void Opers
::MulAdd_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_ALL(a_o, a_i, a_i1);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      MulAdd_I(a_o, a_i, a_i1, stream));
}


void Opers
::MulSub(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, StreamT stream){
    MK_CHECK4_ALL(a_o, a_i, a_i1, a_i2);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      MulSub(a_o, a_i, a_i1, a_i2, stream));
}


void Opers
::MulSub_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_ALL(a_o, a_i, a_i1);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      MulSub_I(a_o, a_i, a_i1, stream));
}


void Opers
::AddMul(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, StreamT stream){
    MK_CHECK4_ALL(a_o, a_i, a_i1, a_i2);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      AddMul(a_o, a_i, a_i1, a_i2, stream));
}


void Opers
::AddMul_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_ALL(a_o, a_i, a_i1);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      AddMul_I(a_o, a_i, a_i1, stream));
}


void Opers
::AddDiv(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, StreamT stream){
    MK_CHECK4_ALL(a_o, a_i, a_i1, a_i2);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      AddDiv(a_o, a_i, a_i1, a_i2, stream));
}

void Opers
::AddDiv_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_ALL(a_o, a_i, a_i1);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      AddDiv_I(a_o, a_i, a_i1, stream));
}


void Opers
::SubMul(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, StreamT stream){
    MK_CHECK4_ALL(a_o, a_i, a_i1, a_i2);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      SubMul(a_o, a_i, a_i1, a_i2, stream));
}

void Opers
::SubMul_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_ALL(a_o, a_i, a_i1);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      SubMul_I(a_o, a_i, a_i1, stream));
}

void Opers
::SubDiv(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, StreamT stream){
    MK_CHECK4_ALL(a_o, a_i, a_i1, a_i2);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      SubDiv(a_o, a_i, a_i1, a_i2, stream));
}


void Opers
::SubDiv_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_ALL(a_o, a_i, a_i1);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      SubDiv_I(a_o, a_i, a_i1, stream));
}



void Opers
::MulCAdd(Image3D& a_o, const Image3D& a_i, const float& c, const Image3D& a_i1, StreamT stream, bool onDev){
    MK_CHECK3_ALL(a_o, a_i, a_i1);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      MulCAdd(a_o, a_i, c, a_i1, stream, onDev));
}


void Opers
::MulCAdd_I(Image3D& a_o, const float& c, const Image3D& a_i, StreamT stream, bool onDev){
    MK_CHECK2_ALL(a_o, a_i);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      MulCAdd_I(a_o, c, a_i, stream, onDev));
}


void Opers
::MulCSub(Image3D& a_o, const Image3D& a_i, const float& c, const Image3D& a_i1, StreamT stream, bool onDev){
    MK_CHECK3_ALL(a_o, a_i, a_i1);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      MulCSub(a_o, a_i, c, a_i1, stream, onDev));
}


void Opers
::MulCSub_I(Image3D& a_o, const float& c, const Image3D& a_i, StreamT stream, bool onDev){
    MK_CHECK2_ALL(a_o, a_i);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      MulCSub_I(a_o, c, a_i, stream, onDev));
}



void Opers
::MulCAddC(Image3D& a_o, const Image3D& a_i, const float& a, const float& b, StreamT stream, bool onDev){
    MK_CHECK2_ALL(a_o, a_i);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      MulCAddC(a_o, a_i, a, b, stream, onDev));
}


void Opers
::MulCAddC_I(Image3D& a_o, const float& a, const float& b, StreamT stream, bool onDev){
   AUTO_EXEC(a_o.memType(), ImageOpers,
	     MulCAddC_I(a_o, a, b, stream, onDev));
}



void Opers
::AddMulC(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const float& c, StreamT stream, bool onDev){
    MK_CHECK3_ALL(a_o, a_i, a_i1);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      AddMulC(a_o, a_i, a_i1, c, stream, onDev));
}


void Opers
::AddMulC_I(Image3D& a_o, const Image3D& a_i, const float& c, StreamT stream, bool onDev){
    MK_CHECK2_ALL(a_o, a_i);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      AddMulC_I(a_o, a_i, c, stream, onDev));
}


void Opers
::SubMulC(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const float& c, StreamT stream, bool onDev){
    MK_CHECK3_ALL(a_o, a_i, a_i1);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      SubMulC(a_o, a_i, a_i1, c, stream, onDev));
}


void Opers
::SubMulC_I(Image3D& a_o, const Image3D& a_i, const float& c, StreamT stream, bool onDev){
    MK_CHECK2_ALL(a_o, a_i);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      SubMulC_I(a_o, a_i, c, stream, onDev));
}


void Opers
::AddCMulC(Image3D& a_o, const Image3D& a_i, const float& a, const float& b, StreamT stream, bool onDev){
    MK_CHECK2_ALL(a_o, a_i);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      AddCMulC(a_o, a_i, a, b, stream, onDev));
}


void Opers
::AddCMulC_I(Image3D& a_o, const float& a, const float& b, StreamT stream, bool onDev){
   AUTO_EXEC(a_o.memType(), ImageOpers,
	     AddCMulC_I(a_o, a, b, stream, onDev));
}


void Opers
::Add_MulC(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const float& c, StreamT stream, bool onDev){
    MK_CHECK3_ALL(a_o, a_i, a_i1);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      Add_MulC(a_o, a_i, a_i1, c, stream, onDev));
}


void Opers
::Add_MulC_I(Image3D& a_o, const Image3D& a_i, const float& c, StreamT stream, bool onDev){
    MK_CHECK2_ALL(a_o, a_i);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      Add_MulC_I(a_o, a_i, c, stream, onDev));
}


void Opers
::Add_Mul(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, StreamT stream){
    MK_CHECK4_ALL(a_o, a_i, a_i1, a_i2);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      Add_Mul(a_o, a_i, a_i1, a_i2, stream));
}


void Opers
::Add_Mul_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_ALL(a_o, a_i, a_i1);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      Add_Mul_I(a_o, a_i, a_i1, stream));
}



void Opers
::Sub_Mul(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, StreamT stream){
    MK_CHECK4_ALL(a_o, a_i, a_i1, a_i2);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      Sub_Mul(a_o, a_i, a_i1, a_i2, stream));
}


void Opers
::Sub_Mul_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, StreamT stream){
    MK_CHECK3_ALL(a_o, a_i, a_i1);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      Sub_Mul_I(a_o, a_i, a_i1, stream));
}


void Opers
::Add_AddMulC(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, const float& c, StreamT stream, bool onDev){
    MK_CHECK4_ALL(a_o, a_i, a_i1, a_i2);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      Add_AddMulC(a_o, a_i, a_i1, a_i2, c, stream, onDev));
}



void Opers
::Add_SubMulC(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, const float& c, StreamT stream, bool onDev){
    MK_CHECK4_ALL(a_o, a_i, a_i1, a_i2);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      Add_SubMulC(a_o, a_i, a_i1, a_i2, c, stream, onDev));
}



void Opers
::Add_MulMulC(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, const float& c, StreamT stream, bool onDev){
    MK_CHECK4_ALL(a_o, a_i, a_i1, a_i2);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      Add_MulMulC(a_o, a_i, a_i1, a_i2, c, stream, onDev));
}


void Opers
::Add_AddMulC_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const float& c, StreamT stream, bool onDev){
    MK_CHECK3_ALL(a_o, a_i, a_i1);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      Add_AddMulC_I(a_o, a_i, a_i1, c, stream, onDev));
}


void Opers
::Add_SubMulC_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const float& c, StreamT stream, bool onDev){
    MK_CHECK3_ALL(a_o, a_i, a_i1);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      Add_SubMulC_I(a_o, a_i, a_i1, c, stream, onDev));
}


void Opers
::Add_MulMulC_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const float& c, StreamT stream, bool onDev){
    MK_CHECK3_ALL(a_o, a_i, a_i1);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      Add_MulMulC_I(a_o, a_i, a_i1, c, stream, onDev));
}


void Opers
::Sub_AddMulC(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, const float& c, StreamT stream, bool onDev){
    MK_CHECK4_ALL(a_o, a_i, a_i1, a_i2);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      Sub_AddMulC(a_o, a_i, a_i1, a_i2, c, stream, onDev));
}



void Opers
::Sub_SubMulC(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, const float& c, StreamT stream, bool onDev){
    MK_CHECK4_ALL(a_o, a_i, a_i1, a_i2);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      Sub_SubMulC(a_o, a_i, a_i1, a_i2, c, stream, onDev));
}



void Opers
::Sub_MulMulC(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const Image3D& a_i2, const float& c, StreamT stream, bool onDev){
    MK_CHECK4_ALL(a_o, a_i, a_i1, a_i2);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      Sub_MulMulC(a_o, a_i, a_i1, a_i2, c, stream, onDev));
}


void Opers
::Sub_AddMulC_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const float& c, StreamT stream, bool onDev){
    MK_CHECK3_ALL(a_o, a_i, a_i1);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      Sub_AddMulC_I(a_o, a_i, a_i1, c, stream, onDev));
}


void Opers
::Sub_SubMulC_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const float& c, StreamT stream, bool onDev){
    MK_CHECK3_ALL(a_o, a_i, a_i1);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      Sub_SubMulC_I(a_o, a_i, a_i1, c, stream, onDev));
}


void Opers
::Sub_MulMulC_I(Image3D& a_o, const Image3D& a_i, const Image3D& a_i1, const float& c, StreamT stream, bool onDev){
    MK_CHECK3_ALL(a_o, a_i, a_i1);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      Sub_MulMulC_I(a_o, a_i, a_i1, c, stream, onDev));
}


void Opers
::MulC_Add_MulC(Image3D& a_o, const Image3D& a_i, const float& a, const Image3D& a_i1, const float& b, StreamT stream, bool onDev){
    MK_CHECK3_ALL(a_o, a_i, a_i1);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      MulC_Add_MulC(a_o, a_i, a, a_i1, b, stream, onDev));
}


void Opers
::MulC_Add_MulC_I(Image3D& a_o, const float& a, const Image3D& a_i, const float& b, StreamT stream, bool onDev){
    MK_CHECK2_ALL(a_o, a_i);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      MulC_Add_MulC_I(a_o, a, a_i, b, stream, onDev));
}


void Opers
::AddCMulCAddC(Image3D& a_o, const Image3D& a_i, const float& a, const float& b, const float& c, StreamT stream, bool onDev){
    MK_CHECK2_ALL(a_o, a_i);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      AddCMulCAddC(a_o, a_i, a, b, c, stream, onDev));
}


void Opers
::AddCMulCAddC_I(Image3D& a_o, const float& a, const float& b, const float& c, StreamT stream, bool onDev){
   AUTO_EXEC(a_o.memType(), ImageOpers,
	     AddCMulCAddC_I(a_o, a, b, c, stream, onDev));
}


void Opers
::ShiftCoordinate(Image3D& a_o, const Image3D& a_i, bool dir, StreamT s)
{
   MK_CHECK2_MEM(a_o, a_i)
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
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      ShiftCoordinate(a_o, a_i, dir, s))
}

void Opers::
SubVol(Image3D& a_o, const Image3D& a_i,
       const Vec3Di& start, StreamT stream)
{
    MK_CHECK2_MEM(a_o, a_i);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      SubVol(a_o, a_i, start, stream));
}

void Opers::
SetSubVol_I(Image3D& a_o,const Image3D& a_i,
	    const Vec3Di& start, StreamT st)
{
   MK_CHECK2_MEM(a_o, a_i);
   AUTO_EXEC(a_o.memType(), ImageOpers,
	     SetSubVol_I(a_o, a_i, start, st));
}

void Opers::
Shrink_I(Image3D& a_o, float c, StreamT stream)
{
   AUTO_EXEC(a_o.memType(), ImageOpers,
	     Shrink_I(a_o, c, stream));
}

void Opers::
Shrink(Image3D& a_o, const Image3D& a_i,
       float c, StreamT stream)
{
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), ImageOpers,
	     Shrink(a_o, a_i, c, stream));
}

void Opers::
SoftAbs_I(Image3D& a_o, float eps, StreamT stream)
{
   AUTO_EXEC(a_o.memType(), ImageOpers,
	     SoftAbs_I(a_o, eps, stream));
}

void Opers::
SoftAbs(Image3D& a_o, const Image3D& a_i,
       float eps, StreamT stream)
{
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), ImageOpers,
	     SoftAbs(a_o, a_i, eps, stream));
}

void Opers::
SoftSgn_I(Image3D& a_o, float eps, StreamT stream)
{
   AUTO_EXEC(a_o.memType(), ImageOpers,
	     SoftSgn_I(a_o, eps, stream));
}

void Opers::
SoftSgn(Image3D& a_o, const Image3D& a_i,
       float eps, StreamT stream)
{
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), ImageOpers,
	     SoftSgn(a_o, a_i, eps, stream));
}

template<BackgroundStrategy bg, InterpT interp, bool useOriginOffset>
void Opers::
Resample(Image3D& a_o, const Image3D& a_i, StreamT stream){
   MK_CHECK2_MEM(a_o, a_i);
   AUTO_EXEC(a_o.memType(), ImageOpers,
	     Resample<bg COMMA interp COMMA useOriginOffset>(a_o, a_i, stream));
}

template<BackgroundStrategy bg, bool useOriginOffset>
void 
resample(Image3D& a_o,const Image3D& a_i,
	 InterpT interp,
	 StreamT stream)
{
   if(interp == INTERP_NN){
       Opers::Resample<bg, INTERP_NN, useOriginOffset>
	   (a_o, a_i, stream);
   }else if(interp == INTERP_LINEAR){
       Opers::Resample<bg, INTERP_LINEAR, useOriginOffset>
	   (a_o, a_i, stream);
   }else if(interp == INTERP_CUBIC){
       Opers::Resample<bg, INTERP_CUBIC, useOriginOffset>
	   (a_o, a_i, stream);
   }else{
       // throw exception
       throw PyCAException(__FILE__, __LINE__, "Unknown interpolation strategy");
   }
}

void Opers::
Resample(Image3D& a_o,const Image3D& a_i,
	 BackgroundStrategy bg, 
	 InterpT interp,
	 bool useOriginOffset,
	 StreamT stream)
{
    if(bg == BACKGROUND_STRATEGY_CLAMP){
	if(useOriginOffset){
	    resample<BACKGROUND_STRATEGY_CLAMP, true>
		(a_o, a_i, interp, stream);
	}else{
	    resample<BACKGROUND_STRATEGY_CLAMP, false>
		(a_o, a_i, interp, stream);
	}
    }else if(bg == BACKGROUND_STRATEGY_WRAP){
	if(useOriginOffset){
	    resample<BACKGROUND_STRATEGY_WRAP, true>
		(a_o, a_i, interp, stream);
	}else{
	    resample<BACKGROUND_STRATEGY_WRAP, false>
		(a_o, a_i, interp, stream);
	}
    }else if(bg == BACKGROUND_STRATEGY_ZERO){
	if(useOriginOffset){
	    resample<BACKGROUND_STRATEGY_ZERO, true>
		(a_o, a_i, interp, stream);
	}else{
	    resample<BACKGROUND_STRATEGY_ZERO, false>
		(a_o, a_i, interp, stream);
	}
    }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ZERO){
	if(useOriginOffset){
	    resample<BACKGROUND_STRATEGY_PARTIAL_ZERO, true>
		(a_o, a_i, interp, stream);
	}else{
	    resample<BACKGROUND_STRATEGY_PARTIAL_ZERO, false>
		(a_o, a_i, interp, stream);
	}
    }else{
	// throw exception
	throw PyCAException(__FILE__, __LINE__, "Unknown background strategy");
    }
}

template<BackgroundStrategy bg, InterpT interp>
void Opers::
ResampleWorld(Image3D& a_o, const Image3D& a_i, StreamT stream){
    MK_CHECK2_MEM(a_o, a_i);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      ResampleWorld<bg COMMA interp>(a_o, a_i, stream));
}

template<BackgroundStrategy bg>
void
resampleWorld(Image3D& a_o,const Image3D& a_i,
	      InterpT interp,
	      StreamT stream)
{
    if(interp == INTERP_NN){
	Opers::ResampleWorld<bg, INTERP_NN>
	    (a_o, a_i, stream);
    }else if(interp == INTERP_LINEAR){
	Opers::ResampleWorld<bg, INTERP_LINEAR>
	    (a_o, a_i, stream);
    }else if(interp == INTERP_CUBIC){
	Opers::ResampleWorld<bg, INTERP_CUBIC>
	    (a_o, a_i, stream);
    }else{
	// throw exception
	throw PyCAException(__FILE__, __LINE__, "Unknown background strategy");
    }
}

void Opers::
ResampleWorld(Image3D& a_o,const Image3D& a_i,
	      BackgroundStrategy bg,
	      InterpT interp,
	      StreamT stream)
{
    if(bg == BACKGROUND_STRATEGY_CLAMP){
	resampleWorld<BACKGROUND_STRATEGY_CLAMP>
	    (a_o, a_i, interp, stream);
    }else if(bg == BACKGROUND_STRATEGY_WRAP){
	resampleWorld<BACKGROUND_STRATEGY_WRAP>
	    (a_o, a_i, interp, stream);
    }else if(bg == BACKGROUND_STRATEGY_ZERO){
	resampleWorld<BACKGROUND_STRATEGY_ZERO>
	    (a_o, a_i, interp, stream);
    }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ZERO){
	resampleWorld<BACKGROUND_STRATEGY_PARTIAL_ZERO>
	    (a_o, a_i, interp, stream);
   }else{
      // throw exception
      throw PyCAException(__FILE__, __LINE__, "Unknown background strategy");
   }
}

template<BackgroundStrategy bg>
void Opers::
SplatWorld(Image3D& a_o, const Image3D& a_i, StreamT stream){
    MK_CHECK2_MEM(a_o, a_i);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      SplatWorld<bg>(a_o, a_i, stream));
}

void Opers::
SplatWorld(Image3D& a_o,const Image3D& a_i,
	 BackgroundStrategy bg,
	 StreamT stream)
{
   if(bg == BACKGROUND_STRATEGY_CLAMP){
      SplatWorld<BACKGROUND_STRATEGY_CLAMP>(a_o, a_i, stream);
   }else if(bg == BACKGROUND_STRATEGY_WRAP){
      SplatWorld<BACKGROUND_STRATEGY_WRAP>(a_o, a_i, stream);
   }else if(bg == BACKGROUND_STRATEGY_ZERO){
      SplatWorld<BACKGROUND_STRATEGY_ZERO>(a_o, a_i, stream);
   }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ZERO){
      SplatWorld<BACKGROUND_STRATEGY_PARTIAL_ZERO>(a_o, a_i, stream);
   }else{
      // throw exception
      throw PyCAException(__FILE__, __LINE__, "Unknown background strategy");
   }
}

template<BackgroundStrategy bg>
void Opers::
SplatWorld(Image3D& a_o, const Image3D& a_i, Image3D& a_w, StreamT stream){
    MK_CHECK3_MEM(a_o, a_i, a_w);
    MK_CHECK2_SIZE(a_o, a_w);
    AUTO_EXEC(a_o.memType(), ImageOpers,
	      SplatWorld<bg>(a_o, a_i, a_w, stream));
}

void Opers::
SplatWorld(Image3D& a_o,const Image3D& a_i, Image3D& a_w,
	 BackgroundStrategy bg,
	 StreamT stream)
{
   if(bg == BACKGROUND_STRATEGY_CLAMP){
      SplatWorld<BACKGROUND_STRATEGY_CLAMP>(a_o, a_i, a_w, stream);
   }else if(bg == BACKGROUND_STRATEGY_WRAP){
      SplatWorld<BACKGROUND_STRATEGY_WRAP>(a_o, a_i, a_w, stream);
   }else if(bg == BACKGROUND_STRATEGY_ZERO){
      SplatWorld<BACKGROUND_STRATEGY_ZERO>(a_o, a_i, a_w, stream);
   }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ZERO){
      SplatWorld<BACKGROUND_STRATEGY_PARTIAL_ZERO>(a_o, a_i, a_w, stream);
   }else{
      // throw exception
      throw PyCAException(__FILE__, __LINE__, "Unknown background strategy");
   }
}

void Opers::
Convolve(Image3D& a_o, const Image3D& a_i,
	 const Image3D& kernel, StreamT stream){
   MK_CHECK2_ALL(a_o, a_i);
   MK_CHECK3_MEM(a_o, a_i, kernel);
   AUTO_EXEC(a_o.memType(), ImageOpers,
	     Convolve(a_o, a_i, kernel, stream));
}

// template instantiations
#include "IOpers_inst.cxx"

} // end namespace PyCA

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

#ifndef __IMG_OPERS_H
#define __IMG_OPERS_H

#ifndef SWIG
#include <Selector.h>
#include <pycaConst.h>
#include <estream.h>
#include <Vec3D.h>
#endif // SWIG

#define IMAGE_UNARY_ARRAY_OP_DEC(OP)   	       	       	       	       	   \
void OP(Image3D& a_o,const Image3D& a_i,				   \
	StreamT st=NULL);

#define IMAGE_UNARY_ARRAY_MASKED_OP_DEC(OP)   	       	       	       	   \
    void OP(Image3D& a_o,const Image3D& a_i, const Image3D& a_mask,	   \
	StreamT st=NULL);

#define IMAGE_UNARY_ARRAY_OP_I_DEC(OP)					   \
void OP##_I(Image3D& a_o,StreamT st=NULL);

#define IMAGE_UNARY_ARRAY_MASKED_OP_I_DEC(OP)				   \
void OP##_I(Image3D& a_o, const Image3D& a_mask, StreamT st=NULL);

#define IMAGE_BINARY_ARRAY_OP_DEC(OP)					   \
void OP(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,		   \
	StreamT st=NULL);
#define IMAGE_BINARY_ARRAY_MASKED_OP_DEC(OP)				   \
void OP(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,               \
        const Image3D& a_mask,				                   \
	StreamT st=NULL);

#define IMAGE_BINARY_ARRAY_OP_I_DEC(OP)					   \
void OP##_I(Image3D& a_o,const Image3D& a_i,				   \
	    StreamT st=NULL);
#define IMAGE_BINARY_ARRAY_MASKED_OP_I_DEC(OP)				   \
void OP##_I(Image3D& a_o,const Image3D& a_i, const Image3D& a_mask,        \
	    StreamT st=NULL);

#define IMAGE_BINARY_ARRAY_OPC_DEC(OP)					   \
void OP##C(Image3D& a_o,const Image3D& a_i,const float& c,		   \
	   StreamT st=NULL,bool onDev=false);
#define IMAGE_BINARY_ARRAY_MASKED_OPC_DEC(OP)				   \
void OP##C(Image3D& a_o,const Image3D& a_i,const float& c,                 \
           const Image3D& a_mask,		                           \
	   StreamT st=NULL,bool onDev=false);

#define IMAGE_BINARY_ARRAY_OPC_I_DEC(OP)       				   \
void OP##C_I(Image3D& a_o,const float& c,      				   \
	     StreamT st=NULL,bool onDev=false);
#define IMAGE_BINARY_ARRAY_MASKED_OPC_I_DEC(OP)       			   \
void OP##C_I(Image3D& a_o,const float& c, const Image3D& a_mask,      	   \
	     StreamT st=NULL,bool onDev=false);

#define IMAGE_BINARY_ARRAY_OP_NOC_DEC(OP)				   \
void OP(Image3D& a_o,const Image3D& a_i,const float& c,		           \
	   StreamT st=NULL,bool onDev=false);
#define IMAGE_BINARY_ARRAY_MASKED_OP_NOC_DEC(OP)	                   \
void OP(Image3D& a_o,const Image3D& a_i,const float& c,	                   \
        const Image3D& a_mask,				                   \
	StreamT st=NULL,bool onDev=false);

#define IMAGE_BINARY_ARRAY_OP_NOC_I_DEC(OP)       			   \
void OP##_I(Image3D& a_o,const float& c,      				   \
	     StreamT st=NULL,bool onDev=false);
#define IMAGE_BINARY_ARRAY_MASKED_OP_NOC_I_DEC(OP)       		   \
void OP##_I(Image3D& a_o,const float& c, const Image3D& a_mask,		   \
	     StreamT st=NULL,bool onDev=false);

#define IMAGE_UNARY_ARRAY_OP_DECS(OP, MOD)                                 \
MOD IMAGE_UNARY_ARRAY_OP_DEC(OP)					   \
MOD IMAGE_UNARY_ARRAY_MASKED_OP_DEC(OP)

#define IMAGE_UNARY_ARRAY_OP_I_DECS(OP, MOD)                               \
MOD IMAGE_UNARY_ARRAY_OP_I_DEC(OP)					   \
MOD IMAGE_UNARY_ARRAY_MASKED_OP_I_DEC(OP)

#define IMAGE_BINARY_ARRAY_OP_DECS(OP, MOD)                                \
MOD IMAGE_BINARY_ARRAY_OP_DEC(OP)					   \
MOD IMAGE_BINARY_ARRAY_MASKED_OP_DEC(OP)

#define IMAGE_BINARY_ARRAY_OP_I_DECS(OP, MOD)				   \
MOD IMAGE_BINARY_ARRAY_OP_I_DEC(OP)					   \
MOD IMAGE_BINARY_ARRAY_MASKED_OP_I_DEC(OP)

#define IMAGE_BINARY_ARRAY_OPC_DECS(OP, MOD)				   \
MOD IMAGE_BINARY_ARRAY_OPC_DEC(OP)					   \
MOD IMAGE_BINARY_ARRAY_OP_NOC_DEC(OP)					   \
MOD IMAGE_BINARY_ARRAY_MASKED_OPC_DEC(OP)				   \
MOD IMAGE_BINARY_ARRAY_MASKED_OP_NOC_DEC(OP)

#define IMAGE_BINARY_ARRAY_OPC_I_DECS(OP, MOD)				   \
MOD IMAGE_BINARY_ARRAY_OPC_I_DEC(OP)					   \
MOD IMAGE_BINARY_ARRAY_OP_NOC_I_DEC(OP)					   \
MOD IMAGE_BINARY_ARRAY_MASKED_OPC_I_DEC(OP)				   \
MOD IMAGE_BINARY_ARRAY_MASKED_OP_NOC_I_DEC(OP)

namespace PyCA {

class Image3D;
class CImageOpers;
class GImageOpers;

template<int mode>
class ImageOpers {
public:
    typedef typename BinSelector<mode, CImageOpers, GImageOpers>::Result Executer;
    enum { exec_mode = mode };

    static void Copy(Image3D& a_o, const Image3D& a_i, const Image3D& a_mask, StreamT st=NULL);

    static void SetMem(Image3D& a_o,const float& c,
		       StreamT st=NULL,bool onDev=false);

    static void SetMem(Image3D& a_o,const float& c, const Image3D& a_mask, 
		       StreamT st=NULL,bool onDev=false);

    /** @brief a_o = |a_i| */
    IMAGE_UNARY_ARRAY_OP_DECS(Abs, static)
    /** @brief a_o = |a_o| */
    IMAGE_UNARY_ARRAY_OP_I_DECS(Abs, static)

    /** @brief a_o = a_i^3 */
    IMAGE_UNARY_ARRAY_OP_DECS(Cube, static)
    /** @brief a_o = a_o^3 */
    IMAGE_UNARY_ARRAY_OP_I_DECS(Cube, static)

    /** @brief a_o = exp(a_i) */
    IMAGE_UNARY_ARRAY_OP_DECS(Exp, static)
    /** @brief a_o = exp(a_o) */
    IMAGE_UNARY_ARRAY_OP_I_DECS(Exp, static)

    /** @brief a_o = exp(a_i) */
    IMAGE_UNARY_ARRAY_OP_DECS(Log, static)
    /** @brief a_o = exp(a_o) */
    IMAGE_UNARY_ARRAY_OP_I_DECS(Log, static)

    /** @brief a_o = -a_i */
    IMAGE_UNARY_ARRAY_OP_DECS(Neg, static)
    /** @brief a_o = -a_o */
    IMAGE_UNARY_ARRAY_OP_I_DECS(Neg, static)

    /** @brief a_o = a_i if a_i > 0, else 0*/
    IMAGE_UNARY_ARRAY_OP_DECS(Ramp, static)
    /** @brief a_o = a_o if a_o > 0, else 0*/
    IMAGE_UNARY_ARRAY_OP_I_DECS(Ramp, static)

    /** @brief a_o = -1 if a_i < 0, 1 if a_i > 0, 0 if a_i == 0 */
    IMAGE_UNARY_ARRAY_OP_DECS(Sgn, static)
    /** @brief a_o = -1 if a_o < 0, 1 if a_o > 0, 0 if a_i == 0 */
    IMAGE_UNARY_ARRAY_OP_I_DECS(Sgn, static)

    /** @brief a_o = 0 if a_i < 0, 1 if a_i > 0, 0.5 if a_i == 0 */
    IMAGE_UNARY_ARRAY_OP_DECS(Step, static)
    /** @brief a_o = 0 if a_o < 0, 1 if a_o > 0, 0.5 if a_i == 0 */
    IMAGE_UNARY_ARRAY_OP_I_DECS(Step, static)

    /** @brief a_o = a_i^2 */
    IMAGE_UNARY_ARRAY_OP_DECS(Sqr, static)
    /** @brief a_o = a_o^2 */
    IMAGE_UNARY_ARRAY_OP_I_DECS(Sqr, static)

    /** @brief a_o = sqrt(a_i) */
    IMAGE_UNARY_ARRAY_OP_DECS(Sqrt, static)
    /** @brief a_o = sqrt(a_o) */
    IMAGE_UNARY_ARRAY_OP_I_DECS(Sqrt, static)

    /** @brief a_o = 1.0/a_i */
    IMAGE_UNARY_ARRAY_OP_DECS(Inv, static)
    /** @brief a_o = 1.0/a_o */
    IMAGE_UNARY_ARRAY_OP_I_DECS(Inv, static)

    /** @brief a_o = ceil(a_i) */
    IMAGE_UNARY_ARRAY_OP_DECS(Ceil, static)
    /** @brief a_o = ceil(a_o) */
    IMAGE_UNARY_ARRAY_OP_I_DECS(Ceil, static)

    /** @brief a_o = floor(a_i) */
    IMAGE_UNARY_ARRAY_OP_DECS(Floor, static)
    /** @brief a_o = floor(a_o) */
    IMAGE_UNARY_ARRAY_OP_I_DECS(Floor, static)

    /** @brief a_o = round(a_i) */
    IMAGE_UNARY_ARRAY_OP_DECS(Round, static)
    /** @brief a_o = round(a_o) */
    IMAGE_UNARY_ARRAY_OP_I_DECS(Round, static)

    /** @brief a_o = sin(a_i) */
    IMAGE_UNARY_ARRAY_OP_DECS(Sin, static)
    /** @brief a_o = sin(a_o) */
    IMAGE_UNARY_ARRAY_OP_I_DECS(Sin, static)

    /** @brief a_o = asin(a_i) */
    IMAGE_UNARY_ARRAY_OP_DECS(Asin, static)
    /** @brief a_o = asin(a_o) */
    IMAGE_UNARY_ARRAY_OP_I_DECS(Asin, static)

    /** @brief a_o = cos(a_i) */
    IMAGE_UNARY_ARRAY_OP_DECS(Cos, static)
    /** @brief a_o = cos(a_o) */
    IMAGE_UNARY_ARRAY_OP_I_DECS(Cos, static)

    /** @brief a_o = acos(a_i) */
    IMAGE_UNARY_ARRAY_OP_DECS(Acos, static)
    /** @brief a_o = acos(a_o) */
    IMAGE_UNARY_ARRAY_OP_I_DECS(Acos, static)

    /** @brief a_o = tan(a_i) */
    IMAGE_UNARY_ARRAY_OP_DECS(Tan, static)
    /** @brief a_o = tan(a_o) */
    IMAGE_UNARY_ARRAY_OP_I_DECS(Tan, static)

    /** @brief a_o = atan(a_i) */
    IMAGE_UNARY_ARRAY_OP_DECS(Atan, static)
    /** @brief a_o = atan(a_o) */
    IMAGE_UNARY_ARRAY_OP_I_DECS(Atan, static)

    /** @brief a_o = csc(a_i) */
    IMAGE_UNARY_ARRAY_OP_DECS(Csc, static)
    /** @brief a_o = csc(a_o) */
    IMAGE_UNARY_ARRAY_OP_I_DECS(Csc, static)

    /** @brief a_o = sec(a_i) */
    IMAGE_UNARY_ARRAY_OP_DECS(Sec, static)
    /** @brief a_o = sec(a_o) */
    IMAGE_UNARY_ARRAY_OP_I_DECS(Sec, static)

    /** @brief a_o = cot(a_i) */
    IMAGE_UNARY_ARRAY_OP_DECS(Cot, static)
    /** @brief a_o = cot(a_o) */
    IMAGE_UNARY_ARRAY_OP_I_DECS(Cot, static)

    //
    // BINARY OPERATORS
    //

    /** @brief a_o = a_i + a_i1 */
    IMAGE_BINARY_ARRAY_OP_DECS(Add, static)
    /** @brief a_o += a_i1 */
    IMAGE_BINARY_ARRAY_OP_I_DECS(Add, static)
    /** @brief a_o = a_i + c */
    IMAGE_BINARY_ARRAY_OPC_DECS(Add, static)
    /** @brief a_o = a_o + c */
    IMAGE_BINARY_ARRAY_OPC_I_DECS(Add, static)

    /** @brief a_o = a_i - a_i1 */
    IMAGE_BINARY_ARRAY_OP_DECS(Sub, static)
    /** @brief a_o -= a_i1 */
    IMAGE_BINARY_ARRAY_OP_I_DECS(Sub, static)
    /** @brief a_o = a_i - c */
    IMAGE_BINARY_ARRAY_OPC_DECS(Sub, static)
    /** @brief a_o = a_o - c */
    IMAGE_BINARY_ARRAY_OPC_I_DECS(Sub, static)

    /** @brief a_o = a_i * a_i1 */
    IMAGE_BINARY_ARRAY_OP_DECS(Mul, static)
    /** @brief a_o *= a_i1 */
    IMAGE_BINARY_ARRAY_OP_I_DECS(Mul, static)
    /** @brief a_o = a_i * c */
    IMAGE_BINARY_ARRAY_OPC_DECS(Mul, static)
    /** @brief a_o = a_o * c */
    IMAGE_BINARY_ARRAY_OPC_I_DECS(Mul, static)

    /** @brief a_o = a_i / a_i1 */
    IMAGE_BINARY_ARRAY_OP_DECS(Div, static)
    /** @brief a_o /= a_i1 */
    IMAGE_BINARY_ARRAY_OP_I_DECS(Div, static)
    /** @brief a_o = a_i / c */
    IMAGE_BINARY_ARRAY_OPC_DECS(Div, static)
    /** @brief a_o = a_o / c */
    IMAGE_BINARY_ARRAY_OPC_I_DECS(Div, static)

    /** @brief a_o = max(a_i, a_i1) */
    IMAGE_BINARY_ARRAY_OP_DECS(Max, static)
    /** @brief a_o = max(a_o, a_i) */
    IMAGE_BINARY_ARRAY_OP_I_DECS(Max, static)
    /** @brief a_o = max(a_i, c) */
    IMAGE_BINARY_ARRAY_OPC_DECS(Max, static)
    /** @brief a_o = max(a_o, c) */
    IMAGE_BINARY_ARRAY_OPC_I_DECS(Max, static)

    /** @brief a_o = min(a_i, a_i1) */
    IMAGE_BINARY_ARRAY_OP_DECS(Min, static)
    /** @brief a_o = min(a_o, a_i) */
    IMAGE_BINARY_ARRAY_OP_I_DECS(Min, static)
    /** @brief a_o = min(a_i, c) */
    IMAGE_BINARY_ARRAY_OPC_DECS(Min, static)
    /** @brief a_o = min(a_o, c) */
    IMAGE_BINARY_ARRAY_OPC_I_DECS(Min, static)

    /** @brief a_o = atan2(a_i, a_i1) */
    IMAGE_BINARY_ARRAY_OP_DECS(Atan2, static)
    /** @brief a_o = atan2(a_o, a_i) */
    IMAGE_BINARY_ARRAY_OP_I_DECS(Atan2, static)
    /** @brief a_o = atan2(a_i, c) */
    IMAGE_BINARY_ARRAY_OPC_DECS(Atan2, static)
    /** @brief a_o = atan2(a_o, c) */
    IMAGE_BINARY_ARRAY_OPC_I_DECS(Atan2, static)

    /** @brief a_o = a_i > a_i1 */
    IMAGE_BINARY_ARRAY_OP_DECS(GT, static)
    /** @brief a_o = a_o > a_i */
    IMAGE_BINARY_ARRAY_OP_I_DECS(GT, static)
    /** @brief a_o = a_i > c) */
    IMAGE_BINARY_ARRAY_OPC_DECS(GT, static)
    /** @brief a_o = a_o > c) */
    IMAGE_BINARY_ARRAY_OPC_I_DECS(GT, static)

    /** @brief a_o = a_i >= a_i1 */
    IMAGE_BINARY_ARRAY_OP_DECS(GTE, static)
    /** @brief a_o = a_o >= a_i */
    IMAGE_BINARY_ARRAY_OP_I_DECS(GTE, static)
    /** @brief a_o = a_i >= c) */
    IMAGE_BINARY_ARRAY_OPC_DECS(GTE, static)
    /** @brief a_o = a_o >= c) */
    IMAGE_BINARY_ARRAY_OPC_I_DECS(GTE, static)

    /** @brief a_o = a_i == a_i1 */
    IMAGE_BINARY_ARRAY_OP_DECS(EQ, static)
    /** @brief a_o = a_o == a_i */
    IMAGE_BINARY_ARRAY_OP_I_DECS(EQ, static)
    /** @brief a_o = a_i == c) */
    IMAGE_BINARY_ARRAY_OPC_DECS(EQ, static)
    /** @brief a_o = a_o == c) */
    IMAGE_BINARY_ARRAY_OPC_I_DECS(EQ, static)

    /** @brief a_o = a_i != a_i1 */
    IMAGE_BINARY_ARRAY_OP_DECS(NEQ, static)
    /** @brief a_o = a_o != a_i */
    IMAGE_BINARY_ARRAY_OP_I_DECS(NEQ, static)
    /** @brief a_o = a_i != c) */
    IMAGE_BINARY_ARRAY_OPC_DECS(NEQ, static)
    /** @brief a_o = a_o != c) */
    IMAGE_BINARY_ARRAY_OPC_I_DECS(NEQ, static)

    /** @brief a_o = a_i < a_i1 */
    IMAGE_BINARY_ARRAY_OP_DECS(LT, static)
    /** @brief a_o = a_o < a_i */
    IMAGE_BINARY_ARRAY_OP_I_DECS(LT, static)
    /** @brief a_o = a_i < c) */
    IMAGE_BINARY_ARRAY_OPC_DECS(LT, static)
    /** @brief a_o = a_o < c) */
    IMAGE_BINARY_ARRAY_OPC_I_DECS(LT, static)

    /** @brief a_o = a_i <= a_i1 */
    IMAGE_BINARY_ARRAY_OP_DECS(LTE, static)
    /** @brief a_o = a_o <= a_i */
    IMAGE_BINARY_ARRAY_OP_I_DECS(LTE, static)
    /** @brief a_o = a_i <= c) */
    IMAGE_BINARY_ARRAY_OPC_DECS(LTE, static)
    /** @brief a_o = a_o <= c) */
    IMAGE_BINARY_ARRAY_OPC_I_DECS(LTE, static)


    static void PowC(Image3D& a_o,const Image3D& a_i,const float& c,StreamT st=NULL,bool onDev=false);
    static void PowC_I(Image3D& a_o,const float& c,StreamT st=NULL,bool onDev=false);

    static void AbsDiff(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,StreamT st=NULL);
    static void AbsDiff_I(Image3D& a_o,const Image3D& a_i,StreamT st=NULL);

    static void SqrDiff(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,StreamT st=NULL);
    static void SqrDiff_I(Image3D& a_o,const Image3D& a_i,StreamT st=NULL);

    static void MulMulC(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const float& c,StreamT st=NULL,bool onDev=false);
    static void MulMulC_I(Image3D& a_o,const Image3D& a_i,const float& c,StreamT st=NULL,bool onDev=false);

    static void MulMul(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,StreamT st=NULL);
    static void MulMul_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,StreamT st=NULL);
    static void MulAdd(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,StreamT st=NULL);
    static void MulAdd_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,StreamT st=NULL);
    static void MulSub(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,StreamT st=NULL);
    static void MulSub_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,StreamT st=NULL);


    static void AddMul(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,StreamT st=NULL);
    static void AddMul_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,StreamT st=NULL);
    static void AddDiv(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,StreamT st=NULL);
    static void AddDiv_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,StreamT st=NULL);

    static void SubMul(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,StreamT st=NULL);
    static void SubMul_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,StreamT st=NULL);
    static void SubDiv(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,StreamT st=NULL);
    static void SubDiv_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,StreamT st=NULL);

    static void MulCAdd(Image3D& a_o,const Image3D& a_i,const float& c,const Image3D& a_i1,StreamT st=NULL,bool onDev=false);
    static void MulCAdd_I(Image3D& a_o,const float& c,const Image3D& a_i,StreamT st=NULL,bool onDev=false);
    static void MulCSub(Image3D& a_o,const Image3D& a_i,const float& c,const Image3D& a_i1,StreamT st=NULL,bool onDev=false);
    static void MulCSub_I(Image3D& a_o,const float& c,const Image3D& a_i,StreamT st=NULL,bool onDev=false);

    static void MulCAddC(Image3D& a_o,const Image3D& a_i,const float& a,const float& b,StreamT st=NULL,bool onDev=false);
    static void MulCAddC_I(Image3D& a_o,const float& a,const float& b,StreamT st=NULL,bool onDev=false);

    static void AddMulC(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const float& c,StreamT st=NULL,bool onDev=false);
    static void AddMulC_I(Image3D& a_o,const Image3D& a_i,const float& c,StreamT st=NULL,bool onDev=false);
    static void SubMulC(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const float& c,StreamT st=NULL,bool onDev=false);
    static void SubMulC_I(Image3D& a_o,const Image3D& a_i,const float& c,StreamT st=NULL,bool onDev=false);


    static void AddCMulC(Image3D& a_o,const Image3D& a_i,const float& a,const float& b,StreamT st=NULL,bool onDev=false);
    static void AddCMulC_I(Image3D& a_o,const float& a,const float& b,StreamT st=NULL,bool onDev=false);

    static void Add_MulC(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const float& c,StreamT st=NULL,bool onDev=false);
    static void Add_MulC_I(Image3D& a_o,const Image3D& a_i,const float& c,StreamT st=NULL,bool onDev=false);
    static void Add_Mul(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,StreamT st=NULL);
    static void Add_Mul_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,StreamT st=NULL);

    static void Sub_Mul(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,StreamT st=NULL);
    static void Sub_Mul_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,StreamT st=NULL);


    static void Add_AddMulC(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,const float& c,StreamT st=NULL,bool onDev=false);
    static void Add_SubMulC(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,const float& c,StreamT st=NULL,bool onDev=false);
    static void Add_MulMulC(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,const float& c,StreamT st=NULL,bool onDev=false);
    static void Add_AddMulC_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const float& c,StreamT st=NULL,bool onDev=false);
    static void Add_SubMulC_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const float& c,StreamT st=NULL,bool onDev=false);
    static void Add_MulMulC_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const float& c,StreamT st=NULL,bool onDev=false);
static void Sub_AddMulC(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,const float& c,StreamT st=NULL,bool onDev=false);
    static void Sub_SubMulC(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,const float& c,StreamT st=NULL,bool onDev=false);
    static void Sub_MulMulC(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,const float& c,StreamT st=NULL,bool onDev=false);
    static void Sub_AddMulC_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const float& c,StreamT st=NULL,bool onDev=false);
    static void Sub_SubMulC_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const float& c,StreamT st=NULL,bool onDev=false);
    static void Sub_MulMulC_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const float& c,StreamT st=NULL,bool onDev=false);


    static void MulC_Add_MulC(Image3D& a_o,const Image3D& a_i,const float& a,const Image3D& a_i1,const float& b,StreamT st=NULL,bool onDev=false);
    static void MulC_Add_MulC_I(Image3D& a_o,const float& a,const Image3D& a_i,const float& b,StreamT st=NULL,bool onDev=false);
    static void AddCMulCAddC(Image3D& a_o,const Image3D& a_i,const float& a,const float& b,const float& c,StreamT st=NULL,bool onDev=false);
    static void AddCMulCAddC_I(Image3D& a_o,const float& a,const float& b,const float& c,StreamT st=NULL,bool onDev=false);

    static void ShiftCoordinate(Image3D& a_o,const Image3D& a_i,bool dir,StreamT s=NULL);

   static void SubVol(Image3D& a_o,const Image3D& a_i,
		      const Vec3Di& start, StreamT st=NULL);

   static void SetSubVol_I(Image3D& a_o,const Image3D& a_i,
			   const Vec3Di& start, StreamT st=NULL);

   static void Shrink_I(Image3D& a_o, float c, StreamT st = NULL);

   static void Shrink(Image3D& a_o, const Image3D& a_i,
		      float c, StreamT st = NULL);

   static void SoftAbs_I(Image3D& a_o,
			 float eps, StreamT st);
   static void SoftAbs(Image3D& a_o, const Image3D& a_i,
		       float eps, StreamT st);

   static void SoftSgn_I(Image3D& a_o,
			 float eps, StreamT st);
   static void SoftSgn(Image3D& a_o, const Image3D& a_i,
		       float eps, StreamT st);

    // Advance function
    template<BackgroundStrategy bg, InterpT interp, bool useOriginOffset>
    static void Resample(Image3D& a_o,const Image3D& a_i, StreamT s=NULL);

   template<BackgroundStrategy bg, InterpT interp>
   static void ResampleWorld(Image3D& a_o,const Image3D& a_i, StreamT s=NULL);

   template<BackgroundStrategy bg>
   static void SplatWorld(Image3D& a_o,const Image3D& a_i, StreamT s=NULL);
   template<BackgroundStrategy bg>
   static void SplatWorld(Image3D& a_o,const Image3D& a_i, Image3D& a_w, StreamT s=NULL);

   static void Convolve(Image3D& a_o,const Image3D& a_i,
			const Image3D& kernel, StreamT s=NULL);

    inline static void NormalizeImageRange(Image3D& a_o, const float& minV, const float& maxV, StreamT st=NULL, bool onDev=true) {
        if (onDev == false) {
            AddCMulC_I(a_o, -minV, 1.f / (maxV - minV), st, false);
        } else {
            //printErrorMessage("Have not implemented");
        }
    }
};

} // end namespace PyCA

#endif

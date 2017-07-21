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

#ifndef __I_OPERS_H
#define __I_OPERS_H

#ifndef SWIG
#include <estream.h>
#include <pycaConst.h>
#include <Vec3D.h>
#endif // SWIG

#include <ImageOpers.h>

namespace PyCA {

class Image3D;

#define NOMOD

namespace Opers {

void Copy(Image3D& a_o,const Image3D& a_i,StreamT st=NULL);

void Copy(Image3D& a_o,const Image3D& a_i, const Image3D& a_mask, 
	  StreamT st=NULL);

/** @brief Sets all pixels to c, where c is a float */
void SetMem(Image3D& a_o,const float& c,StreamT st=NULL,bool onDev=false);
void SetMem(Image3D& a_o,const float& c,const Image3D& a_mask, 
	    StreamT st=NULL,bool onDev=false);

/** @brief a_o = |a_i| */
IMAGE_UNARY_ARRAY_OP_DECS(Abs, NOMOD)
/** @brief a_o = |a_o| */
IMAGE_UNARY_ARRAY_OP_I_DECS(Abs, NOMOD)

/** @brief a_o = a_i^3 */
IMAGE_UNARY_ARRAY_OP_DECS(Cube, NOMOD)
/** @brief a_o = a_o^3 */
IMAGE_UNARY_ARRAY_OP_I_DECS(Cube, NOMOD)

/** @brief a_o = exp(a_i) */
IMAGE_UNARY_ARRAY_OP_DECS(Exp, NOMOD)
/** @brief a_o = exp(a_o) */
IMAGE_UNARY_ARRAY_OP_I_DECS(Exp, NOMOD)

/** @brief a_o = log(a_i) */
IMAGE_UNARY_ARRAY_OP_DECS(Log, NOMOD)
/** @brief a_o = log(a_o) */
IMAGE_UNARY_ARRAY_OP_I_DECS(Log, NOMOD)

/** @brief a_o = -a_i */
IMAGE_UNARY_ARRAY_OP_DECS(Neg, NOMOD)
/** @brief a_o = -a_o */
IMAGE_UNARY_ARRAY_OP_I_DECS(Neg, NOMOD)

/** @brief a_o = a_i if a_i > 0, else 0*/
IMAGE_UNARY_ARRAY_OP_DECS(Ramp, NOMOD)
/** @brief a_o = a_o if a_o > 0, else 0*/
IMAGE_UNARY_ARRAY_OP_I_DECS(Ramp, NOMOD)

/** @brief a_o = -1 if a_i < 0, 1 if a_i > 0, 0 if a_i == 0 */
IMAGE_UNARY_ARRAY_OP_DECS(Sgn, NOMOD)
/** @brief a_o = -1 if a_o < 0, 1 if a_o > 0, 0 if a_i == 0 */
IMAGE_UNARY_ARRAY_OP_I_DECS(Sgn, NOMOD)

/** @brief a_o = 0 if a_i < 0, 1 if a_i > 0, 0.5 if a_i == 0 */
IMAGE_UNARY_ARRAY_OP_DECS(Step, NOMOD)
/** @brief a_o = 0 if a_o < 0, 1 if a_o > 0, 0.5 if a_i == 0 */
IMAGE_UNARY_ARRAY_OP_I_DECS(Step, NOMOD)

/** @brief a_o = a_i^2 */
IMAGE_UNARY_ARRAY_OP_DECS(Sqr, NOMOD)
/** @brief a_o = a_o^2 */
IMAGE_UNARY_ARRAY_OP_I_DECS(Sqr, NOMOD)

/** @brief a_o = sqrt(a_i) */
IMAGE_UNARY_ARRAY_OP_DECS(Sqrt, NOMOD)
/** @brief a_o = sqrt(a_o) */
IMAGE_UNARY_ARRAY_OP_I_DECS(Sqrt, NOMOD)

/** @brief a_o = 1.0/a_i */
IMAGE_UNARY_ARRAY_OP_DECS(Inv, NOMOD)
/** @brief a_o = 1.0/a_o */
IMAGE_UNARY_ARRAY_OP_I_DECS(Inv, NOMOD)

/** @brief a_o = ceil(a_i) */
IMAGE_UNARY_ARRAY_OP_DECS(Ceil, NOMOD)
/** @brief a_o = ceil(a_o) */
IMAGE_UNARY_ARRAY_OP_I_DECS(Ceil, NOMOD)

/** @brief a_o = floor(a_i) */
IMAGE_UNARY_ARRAY_OP_DECS(Floor, NOMOD)
/** @brief a_o = floor(a_o) */
IMAGE_UNARY_ARRAY_OP_I_DECS(Floor, NOMOD)

/** @brief a_o = round(a_i) */
IMAGE_UNARY_ARRAY_OP_DECS(Round, NOMOD)
/** @brief a_o = round(a_o) */
IMAGE_UNARY_ARRAY_OP_I_DECS(Round, NOMOD)

/** @brief a_o = sin(a_i) */
IMAGE_UNARY_ARRAY_OP_DECS(Sin, NOMOD)
/** @brief a_o = sin(a_o) */
IMAGE_UNARY_ARRAY_OP_I_DECS(Sin, NOMOD)

/** @brief a_o = asin(a_i) */
IMAGE_UNARY_ARRAY_OP_DECS(Asin, NOMOD)
/** @brief a_o = asin(a_o) */
IMAGE_UNARY_ARRAY_OP_I_DECS(Asin, NOMOD)

/** @brief a_o = cos(a_i) */
IMAGE_UNARY_ARRAY_OP_DECS(Cos, NOMOD)
/** @brief a_o = cos(a_o) */
IMAGE_UNARY_ARRAY_OP_I_DECS(Cos, NOMOD)

/** @brief a_o = acos(a_i) */
IMAGE_UNARY_ARRAY_OP_DECS(Acos, NOMOD)
/** @brief a_o = acos(a_o) */
IMAGE_UNARY_ARRAY_OP_I_DECS(Acos, NOMOD)

/** @brief a_o = tan(a_i) */
IMAGE_UNARY_ARRAY_OP_DECS(Tan, NOMOD)
/** @brief a_o = tan(a_o) */
IMAGE_UNARY_ARRAY_OP_I_DECS(Tan, NOMOD)

/** @brief a_o = atan(a_i) */
IMAGE_UNARY_ARRAY_OP_DECS(Atan, NOMOD)
/** @brief a_o = atan(a_o) */
IMAGE_UNARY_ARRAY_OP_I_DECS(Atan, NOMOD)

/** @brief a_o = csc(a_i) */
IMAGE_UNARY_ARRAY_OP_DECS(Csc, NOMOD)
/** @brief a_o = csc(a_o) */
IMAGE_UNARY_ARRAY_OP_I_DECS(Csc, NOMOD)

/** @brief a_o = sec(a_i) */
IMAGE_UNARY_ARRAY_OP_DECS(Sec, NOMOD)
/** @brief a_o = sec(a_o) */
IMAGE_UNARY_ARRAY_OP_I_DECS(Sec, NOMOD)

/** @brief a_o = cot(a_i) */
IMAGE_UNARY_ARRAY_OP_DECS(Cot, NOMOD)
/** @brief a_o = cot(a_o) */
IMAGE_UNARY_ARRAY_OP_I_DECS(Cot, NOMOD)



/** @brief a_o = a_i + a_i1 */
IMAGE_BINARY_ARRAY_OP_DECS(Add, NOMOD)
/** @brief a_o += a_i1 */
IMAGE_BINARY_ARRAY_OP_I_DECS(Add, NOMOD)
/** @brief a_o = a_i + c */
IMAGE_BINARY_ARRAY_OPC_DECS(Add, NOMOD)
/** @brief a_o = a_o + c */
IMAGE_BINARY_ARRAY_OPC_I_DECS(Add, NOMOD)

/** @brief a_o = a_i - a_i1 */
IMAGE_BINARY_ARRAY_OP_DECS(Sub, NOMOD)
/** @brief a_o -= a_i1 */
IMAGE_BINARY_ARRAY_OP_I_DECS(Sub, NOMOD)
/** @brief a_o = a_i - c */
IMAGE_BINARY_ARRAY_OPC_DECS(Sub, NOMOD)
/** @brief a_o = a_o - c */
IMAGE_BINARY_ARRAY_OPC_I_DECS(Sub, NOMOD)

/** @brief a_o = a_i * a_i1 */
IMAGE_BINARY_ARRAY_OP_DECS(Mul, NOMOD)
/** @brief a_o *= a_i1 */
IMAGE_BINARY_ARRAY_OP_I_DECS(Mul, NOMOD)
/** @brief a_o = a_i * c */
IMAGE_BINARY_ARRAY_OPC_DECS(Mul, NOMOD)
/** @brief a_o = a_o * c */
IMAGE_BINARY_ARRAY_OPC_I_DECS(Mul, NOMOD)

/** @brief a_o = a_i / a_i1 */
IMAGE_BINARY_ARRAY_OP_DECS(Div, NOMOD)
/** @brief a_o /= a_i1 */
IMAGE_BINARY_ARRAY_OP_I_DECS(Div, NOMOD)
/** @brief a_o = a_i / c */
IMAGE_BINARY_ARRAY_OPC_DECS(Div, NOMOD)
/** @brief a_o = a_o / c */
IMAGE_BINARY_ARRAY_OPC_I_DECS(Div, NOMOD)

/** @brief a_o = max(a_i, a_i1) */
IMAGE_BINARY_ARRAY_OP_DECS(Max, NOMOD)
/** @brief a_o = max(a_o, a_i) */
IMAGE_BINARY_ARRAY_OP_I_DECS(Max, NOMOD)
/** @brief a_o = max(a_i, c) */
IMAGE_BINARY_ARRAY_OPC_DECS(Max, NOMOD)
/** @brief a_o = max(a_o, c) */
IMAGE_BINARY_ARRAY_OPC_I_DECS(Max, NOMOD)

/** @brief a_o = min(a_i, a_i1) */
IMAGE_BINARY_ARRAY_OP_DECS(Min, NOMOD)
/** @brief a_o = min(a_o, a_i) */
IMAGE_BINARY_ARRAY_OP_I_DECS(Min, NOMOD)
/** @brief a_o = min(a_i, c) */
IMAGE_BINARY_ARRAY_OPC_DECS(Min, NOMOD)
/** @brief a_o = min(a_o, c) */
IMAGE_BINARY_ARRAY_OPC_I_DECS(Min, NOMOD)

/** @brief a_o = atan2(a_i, a_i1) */
IMAGE_BINARY_ARRAY_OP_DECS(Atan2, NOMOD)
/** @brief a_o = atan2(a_o, a_i) */
IMAGE_BINARY_ARRAY_OP_I_DECS(Atan2, NOMOD)
/** @brief a_o = atan2(a_i, c) */
IMAGE_BINARY_ARRAY_OPC_DECS(Atan2, NOMOD)
/** @brief a_o = atan2(a_o, c) */
IMAGE_BINARY_ARRAY_OPC_I_DECS(Atan2, NOMOD)

/** @brief a_o = a_i > a_i1 */
IMAGE_BINARY_ARRAY_OP_DECS(GT, NOMOD)
/** @brief a_o = a_o > a_i */
IMAGE_BINARY_ARRAY_OP_I_DECS(GT, NOMOD)
/** @brief a_o = a_i > c) */
IMAGE_BINARY_ARRAY_OPC_DECS(GT, NOMOD)
/** @brief a_o = a_o > c) */
IMAGE_BINARY_ARRAY_OPC_I_DECS(GT, NOMOD)

/** @brief a_o = a_i >= a_i1 */
IMAGE_BINARY_ARRAY_OP_DECS(GTE, NOMOD)
/** @brief a_o = a_o >= a_i */
IMAGE_BINARY_ARRAY_OP_I_DECS(GTE, NOMOD)
/** @brief a_o = a_i >= c) */
IMAGE_BINARY_ARRAY_OPC_DECS(GTE, NOMOD)
/** @brief a_o = a_o >= c) */
IMAGE_BINARY_ARRAY_OPC_I_DECS(GTE, NOMOD)

/** @brief a_o = a_i == a_i1 */
IMAGE_BINARY_ARRAY_OP_DECS(EQ, NOMOD)
/** @brief a_o = a_o == a_i */
IMAGE_BINARY_ARRAY_OP_I_DECS(EQ, NOMOD)
/** @brief a_o = a_i == c) */
IMAGE_BINARY_ARRAY_OPC_DECS(EQ, NOMOD)
/** @brief a_o = a_o == c) */
IMAGE_BINARY_ARRAY_OPC_I_DECS(EQ, NOMOD)

/** @brief a_o = a_i != a_i1 */
IMAGE_BINARY_ARRAY_OP_DECS(NEQ, NOMOD)
/** @brief a_o = a_o != a_i */
IMAGE_BINARY_ARRAY_OP_I_DECS(NEQ, NOMOD)
/** @brief a_o = a_i != c) */
IMAGE_BINARY_ARRAY_OPC_DECS(NEQ, NOMOD)
/** @brief a_o = a_o != c) */
IMAGE_BINARY_ARRAY_OPC_I_DECS(NEQ, NOMOD)

/** @brief a_o = a_i < a_i1 */
IMAGE_BINARY_ARRAY_OP_DECS(LT, NOMOD)
/** @brief a_o = a_o < a_i */
IMAGE_BINARY_ARRAY_OP_I_DECS(LT, NOMOD)
/** @brief a_o = a_i < c) */
IMAGE_BINARY_ARRAY_OPC_DECS(LT, NOMOD)
/** @brief a_o = a_o < c) */
IMAGE_BINARY_ARRAY_OPC_I_DECS(LT, NOMOD)

/** @brief a_o = a_i <= a_i1 */
IMAGE_BINARY_ARRAY_OP_DECS(LTE, NOMOD)
/** @brief a_o = a_o <= a_i */
IMAGE_BINARY_ARRAY_OP_I_DECS(LTE, NOMOD)
/** @brief a_o = a_i <= c) */
IMAGE_BINARY_ARRAY_OPC_DECS(LTE, NOMOD)
/** @brief a_o = a_o <= c) */
IMAGE_BINARY_ARRAY_OPC_I_DECS(LTE, NOMOD)



/** @brief a_o = a_i^c */
    void PowC(Image3D& a_o,const Image3D& a_i,const float& c,StreamT st=NULL,bool onDev=false);
/** @brief a_o = a_o^c */
    void PowC_I(Image3D& a_o,const float& c,StreamT st=NULL,bool onDev=false);
/** @brief a_o = |a_i-a_i1| */
    void AbsDiff(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,StreamT st=NULL);
/** @brief a_o = |a_o-a_i| */
    void AbsDiff_I(Image3D& a_o,const Image3D& a_i,StreamT st=NULL);
/** @brief a_o = (a_i-a_i1)^2 */
    void SqrDiff(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,StreamT st=NULL);
/** @brief a_o = (a_o-a_i)^2 */
    void SqrDiff_I(Image3D& a_o,const Image3D& a_i,StreamT st=NULL);
/** @brief a_o = a_i*a_i1*c */
    void MulMulC(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const float& c,StreamT st=NULL,bool onDev=false);
/** @brief a_o = a_o*a_i*c */
    void MulMulC_I(Image3D& a_o,const Image3D& a_i,const float& c,StreamT st=NULL,bool onDev=false);
/** @brief a_o = a_i*a_i1*a_i2 */
    void MulMul(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,StreamT st=NULL);
/** @brief a_o = a_o*a_i*a_i1 */
    void MulMul_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,StreamT st=NULL);
/** @brief a_o = a_i*a_i1 + a_i2 */
    void MulAdd(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,StreamT st=NULL);
/** @brief a_o = a_o*a_i + a_i1 */
    void MulAdd_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,StreamT st=NULL);
/** @brief a_o = a_i*a_i1 - a_i2 */
    void MulSub(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,StreamT st=NULL);
/** @brief a_o = a_o*a_i - a_i1 */
    void MulSub_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,StreamT st=NULL);

/** @brief a_o = (a_i+a_i1)*a_i2 */
    void AddMul(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,StreamT st=NULL);
/** @brief a_o = (a_o+a_i)*a_i1 */
    void AddMul_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,StreamT st=NULL);
/** @brief a_o = (a_i+a_i1)/a_i2 */
    void AddDiv(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,StreamT st=NULL);
/** @brief a_o = (a_o+a_i)/a_i1 */
    void AddDiv_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,StreamT st=NULL);

/** @brief a_o = (a_i-a_i1)*a_i2 */
    void SubMul(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,StreamT st=NULL);
/** @brief a_o = (a_o-a_i)*a_i1 */
    void SubMul_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,StreamT st=NULL);
/** @brief a_o = (a_i-a_i1)/a_i2 */
    void SubDiv(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,StreamT st=NULL);
/** @brief a_o = (a_o+a_i)/a_i1 */
    void SubDiv_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,StreamT st=NULL);

/** @brief a_o = a_i*c + a_i1 */
    void MulCAdd(Image3D& a_o,const Image3D& a_i,const float& c,const Image3D& a_i1,StreamT st=NULL,bool onDev=false);
/** @brief a_o = a_o*c + a_i */
    void MulCAdd_I(Image3D& a_o,const float& c,const Image3D& a_i,StreamT st=NULL,bool onDev=false);
/** @brief a_o = a_i*c - a_i1 */
    void MulCSub(Image3D& a_o,const Image3D& a_i,const float& c,const Image3D& a_i1,StreamT st=NULL,bool onDev=false);
/** @brief a_o = a_o*c - a_i */
    void MulCSub_I(Image3D& a_o,const float& c,const Image3D& a_i,StreamT st=NULL,bool onDev=false);

/** @brief a_o = a_i*a + b */
    void MulCAddC(Image3D& a_o,const Image3D& a_i,const float& a,const float& b,StreamT st=NULL,bool onDev=false);
/** @brief a_o = a_o*a + b */
    void MulCAddC_I(Image3D& a_o,const float& a,const float& b,StreamT st=NULL,bool onDev=false);

/** @brief a_o = (a_i+a_i1)*c */
    void AddMulC(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const float& c,StreamT st=NULL,bool onDev=false);
/** @brief a_o = (a_o+a_i)*c */
    void AddMulC_I(Image3D& a_o,const Image3D& a_i,const float& c,StreamT st=NULL,bool onDev=false);
/** @brief a_o = (a_o-a_i)*c */
    void SubMulC(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const float& c,StreamT st=NULL,bool onDev=false);
/** @brief a_o = (a_o-a_i)*c */
    void SubMulC_I(Image3D& a_o,const Image3D& a_i,const float& c,StreamT st=NULL,bool onDev=false);

/** @brief a_o = (a_i+a)*b */
    void AddCMulC(Image3D& a_o,const Image3D& a_i,const float& a,const float& b,StreamT st=NULL,bool onDev=false);
/** @brief a_o = (a_o+a)*b */
    void AddCMulC_I(Image3D& a_o,const float& a,const float& b,StreamT st=NULL,bool onDev=false);

/** @brief a_o = a_i + a_i1*c */
    void Add_MulC(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const float& c,StreamT st=NULL,bool onDev=false);
/** @brief a_o += a_i*c */
    void Add_MulC_I(Image3D& a_o,const Image3D& a_i,const float& c,StreamT st=NULL,bool onDev=false);
/** @brief a_o = a_i + (a_i1*a_i2) */
    void Add_Mul(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,StreamT st=NULL);
/** @brief a_o = a_o + (a_i*a_i1) */
    void Add_Mul_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,StreamT st=NULL);

/** @brief a_o = a_i - (a_i1*a_i2) */
    void Sub_Mul(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,StreamT st=NULL);
/** @brief a_o = a_o - (a_i*a_i1) */
    void Sub_Mul_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,StreamT st=NULL);

/** @brief a_o = a_i + (a_i1+a_i2)*c */
    void Add_AddMulC(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,const float& c,StreamT st=NULL,bool onDev=false);
/** @brief a_o = a_o + (a_i+a_i1)*c */
    void Add_AddMulC_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const float& c,StreamT st=NULL,bool onDev=false);
/** @brief a_o = a_i + (a_i1-a_i2)*c */
    void Add_SubMulC(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,const float& c,StreamT st=NULL,bool onDev=false);
/** @brief a_o = a_o + (a_i-a_i1)*c */
    void Add_SubMulC_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const float& c,StreamT st=NULL,bool onDev=false);
/** @brief a_o = a_i + a_i1*a_i2*c */
    void Add_MulMulC(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,const float& c,StreamT st=NULL,bool onDev=false);
/** @brief a_o = a_o + a_i*a_i1*c */
    void Add_MulMulC_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const float& c,StreamT st=NULL,bool onDev=false);
/** @brief a_o = a_i - (a_i1+a_i2)*c */
    void Sub_AddMulC(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,const float& c,StreamT st=NULL,bool onDev=false);
/** @brief a_o = a_o - (a_i+a_i1)*c */
    void Sub_AddMulC_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const float& c,StreamT st=NULL,bool onDev=false);
/** @brief a_o = a_i - (a_i1-a_i2)*c */
    void Sub_SubMulC(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,const float& c,StreamT st=NULL,bool onDev=false);
/** @brief a_o = a_o + (a_i-a_i1)*c */
    void Sub_SubMulC_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const float& c,StreamT st=NULL,bool onDev=false);
/** @brief a_o = a_i - a_i1*a_i2*c */
    void Sub_MulMulC(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const Image3D& a_i2,const float& c,StreamT st=NULL,bool onDev=false);
/** @brief a_o = a_o - a_i*a_i1*c */
    void Sub_MulMulC_I(Image3D& a_o,const Image3D& a_i,const Image3D& a_i1,const float& c,StreamT st=NULL,bool onDev=false);

/** @brief a_o = a_i*a + a_i1*b */
    void MulC_Add_MulC(Image3D& a_o,const Image3D& a_i,const float& a,const Image3D& a_i1,const float& b,StreamT st=NULL,bool onDev=false);
/** @brief a_o = a_o*a + a_i*b */
    void MulC_Add_MulC_I(Image3D& a_o,const float& a,const Image3D& a_i,const float& b,StreamT st=NULL,bool onDev=false);
/** @brief a_o = (a_i + a) * b + c */
    void AddCMulCAddC(Image3D& a_o,const Image3D& a_i,const float& a,const float& b,const float& c,StreamT st=NULL,bool onDev=false);
/** @brief a_o = (a_o + a) * b + c */
    void AddCMulCAddC_I(Image3D& a_o,const float& a,const float& b,const float& c,StreamT st=NULL,bool onDev=false);

/**
 * @brief Permute image axes
 *
 * With dir=true, x=y, y=z, z=x.  With dir=false, x=z, y=x, z=y.
 */
void ShiftCoordinate(Image3D& a_o,const Image3D& a_i,bool dir,StreamT s=NULL);

   /**
    * a_o = a_i[start.x:start.x+a_o.size().x,
    *           start.y:start.y+a_o.size().y,
    *           start.z:start.z+a_o.size().z]
    */
   void SubVol(Image3D& a_o,const Image3D& a_i,
		      const Vec3Di& start, StreamT st=NULL);

   /**
    * a_o[start.x:start.x+a_i.size().x,
    *     start.y:start.y+a_i.size().y,
    *     start.z:start.z+a_i.size().z] = a_i
    */
   void SetSubVol_I(Image3D& a_o,const Image3D& a_i,
			   const Vec3Di& start, StreamT st=NULL);

  /**
   * @brief Pointwise 'shrink' algorithm for L1 penalty
   *
   * a_o = max(a_o-c,0) + min(a_o+c,0)
   */
   void Shrink_I(Image3D& a_o, float c, StreamT st = NULL);

  /**
   * @brief Pointwise 'shrink' algorithm for L1 penalty
   *
   * a_o = max(a_i-c,0) + min(a_i+c,0)
   */
   void Shrink(Image3D& a_o, const Image3D& a_i,
		      float c, StreamT st = NULL);

  /**
   * @brief absolute value with 'soft' shape near zero
   *
   * a_o = -a_o-(eps/2.0) if a_o < -eps
   *     =  a_o-(eps/2.0) if a_o > -eps
   *     = (a_o^2)/(2.0*eps) if -eps <= a_o <= eps
   */
   void SoftAbs_I(Image3D& a_o,
			 float eps, StreamT st=NULL);
  /**
   * @brief absolute value with 'soft' shape near zero
   *
   * a_o = -a_i-(eps/2.0) if a_i < -eps
   *     =  a_i-(eps/2.0) if a_i > eps
   *     = (a_i^2)/(2.0*eps) if -eps <= a_i <= eps
   */
   void SoftAbs(Image3D& a_o, const Image3D& a_i,
		       float eps, StreamT st=NULL);

  /**
   * @brief sign function with linear transition near zero
   *
   * a_o = -1 if a_o < -eps
   *     =  1 if a_o > eps
   *     = a_o/eps if -eps <= a_o <= eps
   */
   void SoftSgn_I(Image3D& a_o,
			 float eps, StreamT st=NULL);
  /**
   * @brief sign function with linear transition near zero
   *
   * a_o = -1 if a_i < -eps
   *     =  1 if a_i > eps
   *     = a_i/eps if -eps <= a_i <= eps
   */
   void SoftSgn(Image3D& a_o, const Image3D& a_i,
		       float eps, StreamT st=NULL);


  /**
   * @brief resample a_i onto grid of a_o without regard to origin or spacing
   */
    template<BackgroundStrategy bg, InterpT interp, bool useOriginOffset>
    void Resample(Image3D& a_o,const Image3D& a_i, StreamT s=NULL);

  // non-templated version
  /**
   * @brief resample a_i onto grid of a_o without regard to origin or spacing
   */
   void Resample(Image3D& a_o,const Image3D& a_i,
		 BackgroundStrategy bg=BACKGROUND_STRATEGY_CLAMP,
		 InterpT interp=INTERP_LINEAR,
		 bool useOriginOffset=true,
		 StreamT s=NULL);

  /**
   * @brief resample a_i onto grid of a_o with regard to origin or spacing
   */
    template<BackgroundStrategy bg, InterpT interp>
    void ResampleWorld(Image3D& a_o,const Image3D& a_i, StreamT s=NULL);

   // non-template version
  /**
   * @brief resample a_i onto grid of a_o with regard to origin or spacing
   */
   void ResampleWorld(Image3D& a_o,const Image3D& a_i,
		      BackgroundStrategy bg=BACKGROUND_STRATEGY_CLAMP,
		      InterpT interp=INTERP_LINEAR,
		      StreamT s=NULL);

  /**
   * @brief splat a_i onto grid of a_o with regard to origin or spacing
   */
    template<BackgroundStrategy bg>
    void SplatWorld(Image3D& a_o, const Image3D& a_i, StreamT s=NULL);
    template<BackgroundStrategy bg>
    void SplatWorld(Image3D& a_o, const Image3D& a_i, Image3D& a_w, StreamT s=NULL);

   // non-template version
  /**
   * @brief splat a_i onto grid of a_o with regard to origin or spacing
   */
   void SplatWorld(Image3D& a_o, const Image3D& a_i,
		      BackgroundStrategy bg=BACKGROUND_STRATEGY_CLAMP,
		      StreamT s=NULL);
   void SplatWorld(Image3D& a_o, const Image3D& a_i, Image3D& a_w,
		      BackgroundStrategy bg=BACKGROUND_STRATEGY_CLAMP,
		      StreamT s=NULL);

  /**
   * @brief convolution with clamped boundary conditions
   *
   * a_o = a_i (x) kernel, where (x) is convolution operator
   */
   void Convolve(Image3D& a_o,const Image3D& a_i,
		 const Image3D& kernel, StreamT s=NULL);

} // end namespace Opers

} // end namespace PyCA

#endif

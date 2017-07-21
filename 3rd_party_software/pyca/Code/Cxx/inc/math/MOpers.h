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

#ifndef __OPERS_H
#define __OPERS_H

#include <limits.h>
#include <float.h>
#include <math.h>

#include <algorithm>

#include "pycaConst.h"

#define PYCA_MIN(a,b) ((a)<(b)?(a):(b))
#define PYCA_MAX(a,b) ((a)>(b)?(a):(b))

typedef unsigned int uint;

namespace PyCA {

enum MATH_OPS {
  MATH_Add,
  MATH_And,
  MATH_Abs,
  MATH_Ceil,
  MATH_Exp,
  MATH_Log,
  MATH_Floor,
  MATH_Round,
  MATH_Cube,
  MATH_Div,
  MATH_Inv,
  MATH_Max,
  MATH_Min,
  MATH_Mul,
  MATH_Neg,
  MATH_Or,
  MATH_Pow,
  MATH_Ramp,
  MATH_Set,
  MATH_Sgn,
  MATH_Step,
  MATH_Sqrt,
  MATH_Sqr,
  MATH_Sub,
  MATH_Lerp,
  MATH_Xor,
  MATH_Clamp,
  MATH_Sin,
  MATH_Asin,
  MATH_Cos,
  MATH_Acos,
  MATH_Tan,
  MATH_Atan,
  MATH_Csc,
  MATH_Sec,
  MATH_Cot,
  MATH_Atan2,
  MATH_GT,
  MATH_GTE,
  MATH_EQ,
  MATH_NEQ,
  MATH_LT,
  MATH_LTE,
};

enum CMP_OPS {
  CMP_GREATER,
  CMP_LESS,
  CMP_EQUAL,
  CMP_LESS_OR_EQUAL,
  CMP_GREATER_OR_EQUAL
};

template <typename T, CMP_OPS oper>
class CmpOpers
{
public:

  //unary function
  static __HOSTDEVICE__ uint op(T a, T b){
    switch (oper){
    case CMP_GREATER:
      return (a > b);
    case CMP_LESS :
      return (a < b);
    case CMP_EQUAL:
      return (a == b);
    case CMP_LESS_OR_EQUAL:
      return (a <= b);
    case CMP_GREATER_OR_EQUAL:
      return (a >= b);
    default:
      return 0;
    }
  }
};

template <typename T, MATH_OPS oper>
class MOpers
{
public:
  // is this operator implemented in general template implementation?
  static __HOSTDEVICE__ bool valid(){
    switch (oper){
    case MATH_Add:
    case MATH_Cube:
    case MATH_Div:
    case MATH_Exp:
    case MATH_Log:
    case MATH_Max:
    case MATH_Min:
    case MATH_Mul:
    case MATH_Neg:
    case MATH_Ramp:
    case MATH_Set:
    case MATH_Sgn:
    case MATH_Step:
    case MATH_Sqr:
    case MATH_Sub:
    case MATH_Lerp:
    case MATH_Clamp:
    case MATH_GT:
    case MATH_GTE:
    case MATH_EQ:
    case MATH_NEQ:
    case MATH_LT:
    case MATH_LTE:
      return true;
    case MATH_And:
    case MATH_Abs:
    case MATH_Inv:
    case MATH_Or:
    case MATH_Pow:
    case MATH_Sqrt:
    case MATH_Xor:
    case MATH_Sin:
    case MATH_Asin:
    case MATH_Cos:
    case MATH_Acos:
    case MATH_Tan:
    case MATH_Atan:
    case MATH_Csc:
    case MATH_Sec:
    case MATH_Cot:
    case MATH_Atan2:
      return false;
    default:
      return false;
    }
  }

  //unary function
  static __HOSTDEVICE__ T op(T a){
    switch (oper){
    default:
      // error value
      return ((T)-42000);

    case MATH_Cube:
      return a * a * a;
    case MATH_Neg:
      return -a;
    case MATH_Set:
      return a;
    case MATH_Sqr:
      return a * a;
    case MATH_Ramp:
      return a > 0 ? a : 0;
    case MATH_Sgn:
      return a > ((T)0) ? ((T)1) :
        (a<((T)0) ? ((T)-1) : ((T)0));
    case MATH_Step:
      return a > ((T)0) ? ((T)1) :
        (a<((T)0) ? ((T)0) : ((T)0.5));
    }
  }

  // binary function
  static __HOSTDEVICE__ T op(T a, T b){
    switch (oper){
    default:
      // error value
      return ((T)-42000);

    case MATH_Add:
      return a + b;

    case MATH_Div:
      return a / b;

    case MATH_Mul:
      return a * b;

    case MATH_Min:
      return min(a, b);

    case MATH_Max:
      return max(a, b);

    case MATH_Sub:
      return a - b;

    case MATH_Atan2:
      return atan2f(a, b);

    case MATH_GT:
      return a > b ? ((T)1) : ((T)0);
    case MATH_GTE:
      return a >= b ? ((T)1) : ((T)0);
    case MATH_EQ:
      return a == b ? ((T)1) : ((T)0);
    case MATH_NEQ:
      return a != b ? ((T)1) : ((T)0);
    case MATH_LT:
      return a < b ? ((T)1) : ((T)0);
    case MATH_LTE:
      return a <= b ? ((T)1) : ((T)0);


    }
  }

  static __HOSTDEVICE__ T op(T a , T b, T t){
    switch (oper){
    default:
      // error value
      return ((T)-42000);

    case MATH_Lerp:
      return a + t * (b-a);
    case MATH_Clamp:
      return max(a, min(b, t));
    }
  }

  // inplace operation
  static __HOSTDEVICE__ void iop(T &a, T b){
    switch (oper){
    default:
      // error value
      a = ((T)-42000);
      return;

    case MATH_Add:
      a += b;
      return;

    case MATH_Div:
      a /= b;
      return;

    case MATH_Sub:
      a -= b;
      return;

    case MATH_Mul:
      a *= b;
      return;

    case MATH_Min:
      a = min(a, b);
      return;

    case MATH_Max:
      a = max(a, b);
      return;

    case MATH_Atan2:
      a = atan2f(a, b);
      return;

    case MATH_GT:
      a = a > b ? ((T)1) : ((T)0);
      return;
    case MATH_GTE:
      a = a >= b ? ((T)1) : ((T)0);
      return;
    case MATH_EQ:
      a = a == b ? ((T)1) : ((T)0);
      return;
    case MATH_NEQ:
      a = a != b ? ((T)1) : ((T)0);
      return;
    case MATH_LT:
      a = a < b ? ((T)1) : ((T)0);
      return;
    case MATH_LTE:
      a = a <= b ? ((T)1) : ((T)0);
      return;


    }
  }

  static __HOSTDEVICE__ T identity() {
    return 0;
  }
};


template <MATH_OPS oper>
class MOpers<float, oper>
{
public:

  // is this operator implemented in float specialization?
  static __HOSTDEVICE__ bool valid(){
    switch (oper){
    case MATH_Abs:
    case MATH_Add:
    case MATH_Ceil:
    case MATH_Clamp:
    case MATH_Cube:
    case MATH_Div:
    case MATH_Exp:
    case MATH_Log:
    case MATH_Floor:
    case MATH_Inv:
    case MATH_Lerp:
    case MATH_Max:
    case MATH_Min:
    case MATH_Mul:
    case MATH_Neg:
    case MATH_Pow:
    case MATH_Round:
    case MATH_Ramp:
    case MATH_Set:
    case MATH_Sgn:
    case MATH_Step:
    case MATH_Sqr:
    case MATH_Sqrt:
    case MATH_Sub:
    case MATH_Sin:
    case MATH_Asin:
    case MATH_Cos:
    case MATH_Acos:
    case MATH_Tan:
    case MATH_Atan:
    case MATH_Csc:
    case MATH_Sec:
    case MATH_Cot:
    case MATH_Atan2:
    case MATH_GT:
    case MATH_GTE:
    case MATH_EQ:
    case MATH_NEQ:
    case MATH_LT:
    case MATH_LTE:
      return true;
    case MATH_And:
    case MATH_Or:
    case MATH_Xor:
      return false;

    default:
      return false;
    }
  }

  //unary function
  static __HOSTDEVICE__ float op(float a){
    switch (oper){
    default:
      // error value
      return -42000.f;

    case MATH_Abs:
      return fabs(a);

    case MATH_Cube:
      return a * a * a;

    case MATH_Exp:
      return expf(a);

    case MATH_Log:
      return logf(a);

    case MATH_Set:
      return a;

    case MATH_Sqrt:
      return sqrtf(a);

    case MATH_Sqr:
      return a * a;

    case MATH_Inv:
      return 1.f / a;

    case MATH_Neg:
      return -a;

    case MATH_Ramp:
      return a > 0.f ? a : 0;

    case MATH_Sgn:
      return a > 0.f ? 1.f : (a<0.f ? -1.f : 0.f);

    case MATH_Step:
      return a > 0.f ? 1.f : (a<0.f ? 0.f : 0.5f);

    case MATH_Ceil:
      return ceil(a);

    case MATH_Floor:
      return floor(a);

    case MATH_Round:
      return (a > 0.0) ? floor(a + 0.5) : ceil(a - 0.5);

    case MATH_Sin:
      return sinf(a);

    case MATH_Asin:
      return asinf(a);

    case MATH_Cos:
      return cosf(a);

    case MATH_Acos:
      return acosf(a);

    case MATH_Tan:
      return tanf(a);

    case MATH_Atan:
      return atanf(a);

    case MATH_Csc:
      return 1 / sinf(a);

    case MATH_Sec:
      return 1 / cosf(a);

    case MATH_Cot:
      return 1 / tan(a);

    }
  }

  // binary function
  static __HOSTDEVICE__ float op(float a, float b){
    switch (oper){
    default:
      // error value
      return -42000.f;

    case MATH_Add:
      return a + b;

    case MATH_Div:
      return a / b;

    case MATH_Mul:
      return a * b;

    case MATH_Min:
      return PYCA_MIN(a, b);

    case MATH_Max:
      return PYCA_MAX(a, b);

    case MATH_Sub:
      return a - b;

    case MATH_Pow:
      return pow(a, b);

    case MATH_Atan2:
      return atan2f(a, b);

    case MATH_GT:
      return a > b ? 1.f : 0.f;
    case MATH_GTE:
      return a >= b ? 1.f : 0.f;
    case MATH_EQ:
      return a == b ? 1.f : 0.f;
    case MATH_NEQ:
      return a != b ? 1.f : 0.f;
    case MATH_LT:
      return a < b ? 1.f : 0.f;
    case MATH_LTE:
      return a <= b ? 1.f : 0.f;

    }
  }

  static __HOSTDEVICE__ float op(float a , float b, float t){
    switch (oper){
    default:
      // error value
      return -42000.f;

    case MATH_Lerp:
      return a + t * (b-a);
    case MATH_Clamp:
      return PYCA_MAX(a, PYCA_MIN(b, t));
    }
  }

  // inplace operation
  static __HOSTDEVICE__ void iop(float &a, float b){
    switch (oper){
    default:
      // error value
      a = -42000.f;
      return;

    case MATH_Add:
      a += b;
      return;

    case MATH_Div:
      a /= b;
      return;

    case MATH_Sub:
      a -= b;
      return;

    case MATH_Mul:
      a *= b;
      return;

    case MATH_Min:
      a = PYCA_MIN(a, b);
      return;

    case MATH_Max:
      a = PYCA_MAX(a, b);
      return;

    case MATH_Pow:
      a = pow(a, b);
      return;

    case MATH_Atan2:
      a = atan2f(a, b);
      return;

    case MATH_GT:
      a = a > b ? 1.f : 0.f;
      return;
    case MATH_GTE:
      a = a >= b ? 1.f : 0.f;
      return;
    case MATH_EQ:
      a = a == b ? 1.f : 0.f;
      return;
    case MATH_NEQ:
      a = a != b ? 1.f : 0.f;
      return;
    case MATH_LT:
      a = a < b ? 1.f : 0.f;
      return;
    case MATH_LTE:
      a = a <= b ? 1.f : 0.f;
      return;


    }
  }

  static __HOSTDEVICE__ float identity() {
    switch (oper){
    default:
    case MATH_Add:
      return 0;

    case MATH_Sub:
      return 0;

    case MATH_Mul:
      return 1;

    case MATH_Div:
      return 1;

    case MATH_Min:
      return FLT_MAX;

    case MATH_Max:
      return -FLT_MAX;

    case MATH_Pow:
      return 1;
    }
  }
};

template <MATH_OPS oper>
class MOpers<int, oper>
{
public:

  // is this operator implemented in int specialization?
  static __HOSTDEVICE__ bool valid(){
    switch (oper){
    case MATH_Abs:
    case MATH_Cube:
    case MATH_Exp:
    case MATH_Log:
    case MATH_Sqr:
    case MATH_Neg:
    case MATH_Ramp:
    case MATH_Set:
    case MATH_Sgn:
    case MATH_Step:
    case MATH_Add:
    case MATH_And:
    case MATH_Or:
    case MATH_Xor:
    case MATH_Div:
    case MATH_Mul:
    case MATH_Max:
    case MATH_Min:
    case MATH_Sub:
    case MATH_Clamp:
    case MATH_Round:
    case MATH_Floor:
    case MATH_Ceil:
    case MATH_GT:
    case MATH_GTE:
    case MATH_EQ:
    case MATH_NEQ:
    case MATH_LT:
    case MATH_LTE:
      return true;
    case MATH_Sqrt:
    case MATH_Inv:
    case MATH_Pow:
    case MATH_Lerp:
    case MATH_Sin:
    case MATH_Asin:
    case MATH_Cos:
    case MATH_Acos:
    case MATH_Tan:
    case MATH_Atan:
    case MATH_Csc:
    case MATH_Sec:
    case MATH_Cot:
    case MATH_Atan2:
      return false;
    default:
      return false;
    }
  }

  //unary function
  static __HOSTDEVICE__ int op(int a){
    switch (oper){
    default:
      // error value
      return -42000;

    case MATH_Abs:
      return (a >= 0) ? a : -a;

    case MATH_Cube:
      return a * a * a;

    case MATH_Set:
      return a;

    case MATH_Sqr:
      return a * a;

    case MATH_Neg:
      return -a;

    case MATH_Ramp:
      return a > 0 ? a : 0;

    case MATH_Sgn:
      return a > 0 ? 1 : (a<0 ? -1 : 0);

    case MATH_Step:
      return a > 0 ? 1 : 0;

    case MATH_Ceil:
      return a;

    case MATH_Floor:
      return a;

    case MATH_Round:
      return a;
    }
  }

  // binary function
  static __HOSTDEVICE__ int op(int a, int b){
    switch (oper){
    default:
      // error value
      return -42000;
	    case MATH_Ramp:
	       return a > 0 ? a : 0;

    case MATH_Add:
      return a + b;

    case MATH_And:
      return a & b;

    case MATH_Or:
      return a | b;

    case MATH_Xor:
      return a ^ b;

    case MATH_Div:
      return a / b;

    case MATH_Mul:
      return a * b;

    case MATH_Min:
      return (a < b) ? a : b;

    case MATH_Max:
      return (a > b) ? a : b;

    case MATH_Sub:
      return a - b;

    case MATH_GT:
      return a > b ? 1 : 0;
    case MATH_GTE:
      return a >= b ? 1 : 0;
    case MATH_EQ:
      return a == b ? 1 : 0;
    case MATH_NEQ:
      return a != b ? 1 : 0;
    case MATH_LT:
      return a < b ? 1 : 0;
    case MATH_LTE:
      return a <= b ? 1 : 0;


    }
  }

  // inplace operation
  static __HOSTDEVICE__ void iop(int &a, int b){
    switch (oper){
    default:
      // error value
      a = -42000;
      return;

    case MATH_Add:
      a += b;
      return;

    case MATH_And:
      a &= b;
      return;

    case MATH_Or:
      a |= b;
      return;

    case MATH_Xor:
      a ^= b;
      return;

    case MATH_Div:
      a /= b;
      return;

    case MATH_Sub:
      a -= b;
      return;

    case MATH_Mul:
      a *= b;
      return;

    case MATH_Min:
      a = (a < b) ? a : b;
      return;

    case MATH_Max:
      a = (a > b) ? a : b;
      return;

    case MATH_GT:
      a = a > b ? 1 : 0;
      return;
    case MATH_GTE:
      a = a >= b ? 1 : 0;
      return;
    case MATH_EQ:
      a = a == b ? 1 : 0;
      return;
    case MATH_NEQ:
      a = a != b ? 1 : 0;
      return;
    case MATH_LT:
      a = a < b ? 1 : 0;
      return;
    case MATH_LTE:
      a = a <= b ? 1 : 0;
      return;

    }
  }

  static __HOSTDEVICE__ int identity() {
    switch (oper){
    default:
    case MATH_And:
      return 0xFFFFFFFF;

    case MATH_Or:
      return 0;

    case MATH_Add:
      return 0;

    case MATH_Mul:
      return 1;

    case MATH_Min:
      return INT_MAX;

    case MATH_Max:
      return INT_MIN;
    }
  }

};

template <MATH_OPS oper>
class MOpers<unsigned int, oper>
{
public:

  // is this operator implemented in uint specialization?
  static __HOSTDEVICE__ bool valid(){
    switch (oper){
    case MATH_Abs:
    case MATH_Cube:
    case MATH_Exp:
    case MATH_Log:
    case MATH_Sqr:
    case MATH_Ramp:
    case MATH_Set:
    case MATH_Sgn:
    case MATH_Step:
    case MATH_Add:
    case MATH_And:
    case MATH_Or:
    case MATH_Xor:
    case MATH_Div:
    case MATH_Mul:
    case MATH_Max:
    case MATH_Min:
    case MATH_Clamp:
    case MATH_Floor:
    case MATH_Ceil:
    case MATH_Round:
    case MATH_GT:
    case MATH_GTE:
    case MATH_EQ:
    case MATH_NEQ:
    case MATH_LT:
    case MATH_LTE:
      return true;
    case MATH_Neg:
    case MATH_Sub:
    case MATH_Sqrt:
    case MATH_Inv:
    case MATH_Pow:
    case MATH_Lerp:
    case MATH_Sin:
    case MATH_Asin:
    case MATH_Cos:
    case MATH_Acos:
    case MATH_Tan:
    case MATH_Atan:
    case MATH_Csc:
    case MATH_Sec:
    case MATH_Cot:
    case MATH_Atan2:
      return false;
    default:
      return false;
    }
  }

  //unary function
  static __HOSTDEVICE__ uint op(uint a){
    switch (oper){
    default:
      // error value
      return 42000;

    case MATH_Abs:
      return a;

    case MATH_Cube:
      return a * a * a;

    case MATH_Set:
      return a;

    case MATH_Sqr:
      return a * a;

    case MATH_Ramp:
      return a;

    case MATH_Sgn:
      return a==0 ? 0 : 1;

    case MATH_Step:
      return 1;

    case MATH_Ceil:
      return a;

    case MATH_Floor:
      return a;

    case MATH_Round:
      return a;
    }
  }

  // binary function
  static __HOSTDEVICE__ uint op(uint a, uint b){
    switch (oper){
    default:
      // error value
      return 42000;

    case MATH_Add:
      return a + b;

    case MATH_And:
      return a & b;

    case MATH_Or:
      return a | b;

    case MATH_Xor:
      return a ^ b;

    case MATH_Div:
      return a / b;

    case MATH_Mul:
      return a * b;

    case MATH_Min:
      return (a <= b ) ? a : b;

    case MATH_Max:
      return (a >= b ) ? a : b;

    case MATH_GT:
      return a > b ? 1 : 0;
    case MATH_GTE:
      return a >= b ? 1 : 0;
    case MATH_EQ:
      return a == b ? 1 : 0;
    case MATH_NEQ:
      return a != b ? 1 : 0;
    case MATH_LT:
      return a < b ? 1 : 0;
    case MATH_LTE:
      return a <= b ? 1 : 0;


    }
  }

  // inplace operation
  static __HOSTDEVICE__ void iop(uint &a, uint b){
    switch (oper){
    default:
      // error value
      a=42000;
      return;

    case MATH_Add:
      a += b;
      return;

    case MATH_And:
      a &= b;
      return;

    case MATH_Or:
      a |= b;
      return;

    case MATH_Xor:
      a ^= b;
      return;

    case MATH_Div:
      a /= b;
      return;

    case MATH_Sub:
      a -= b;
      return;

    case MATH_Mul:
      a *= b;
      return;

    case MATH_Min:
      a = (a <= b) ?  a : b;
      return;

    case MATH_Max:
      a = (a >= b) ?  a : b;
      return;

    case MATH_GT:
      a = a > b ? 1 : 0;
      return;
    case MATH_GTE:
      a = a >= b ? 1 : 0;
      return;
    case MATH_EQ:
      a = a == b ? 1 : 0;
      return;
    case MATH_NEQ:
      a = a != b ? 1 : 0;
      return;
    case MATH_LT:
      a = a < b ? 1 : 0;
      return;
    case MATH_LTE:
      a = a <= b ? 1 : 0;
      return;

    }
  }

  static __HOSTDEVICE__ uint identity() {
    switch (oper){
    default:
    case MATH_Mul:
      return 1;
    case MATH_And:
      return 0xFFFFFFFF;
    case MATH_Or:
      return 0;
    case MATH_Add:
      return 0;
    case MATH_Min:
      return 0xFFFFFFFF;
    case MATH_Max:
      return 0;
    }
  }

};

} // end namespace PyCA

#endif

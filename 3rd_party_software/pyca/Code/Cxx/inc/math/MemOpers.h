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

#ifndef  __MEM_OPERATION__H
#define  __MEM_OPERATION__H

#include <estream.h>

#include <MemPool.h>

#include <Selector.h>
#include <CMemOpers.h>
#include <GMemOpers.h>

namespace PyCA {

namespace Opers {
    // MemPool Copy command
    template<typename T>
    void Copy(MemPool<T>& a_o, const MemPool<T>& a_i, int nVox,
              StreamT stream=NULL);
}

#define UNARY_OP_EXEC(OP)      	       	       	       	       	       	    \
template<int mode, typename T>                                              \
void 									    \
MemOpers<mode, T>::							    \
OP(T* a_o, const T* a_i, size_t n, StreamT st){			            \
    MemOpersExec::OP(a_o, a_i, n, st);					    \
}

#define UNARY_MASKED_OP_EXEC(OP)                                            \
template<int mode, typename T>						    \
void									    \
MemOpers<mode, T>::							    \
OP(T* a_o, const T* a_i, const T* a_mask, 				    \
   size_t n, StreamT st)						    \
{									    \
    MemOpersExec::OP(a_o, a_i, a_mask, n, st);				    \
}

#define UNARY_OP_I_EXEC(OP)						    \
template<int mode, typename T>                                              \
void 									    \
MemOpers<mode, T>::							    \
OP##_I(T* a_o, size_t n, StreamT st){				            \
    MemOpersExec::OP##_I(a_o, n, st);					    \
}

#define UNARY_MASKED_OP_I_EXEC(OP)					    \
template<int mode, typename T>						    \
void									    \
MemOpers<mode, T>::							    \
OP##_I(T* a_o, const T* a_mask, 					    \
       size_t n, StreamT st)						    \
{									    \
    MemOpersExec::OP##_I(a_o, a_mask, n, st);				    \
}

#define UNARY_OPS_EXEC(OP)						    \
UNARY_OP_EXEC(OP)							    \
UNARY_MASKED_OP_EXEC(OP)						    \
UNARY_OP_I_EXEC(OP)						            \
UNARY_MASKED_OP_I_EXEC(OP)

#define BINARY_OP_EXEC(OP)						    \
template<int mode, typename T>                                              \
void 									    \
MemOpers<mode, T>::							    \
OP(T* a_o, const T* a_i, const T* a_i1,					    \
		size_t n, StreamT st)				            \
{									    \
    MemOpersExec::OP(a_o, a_i, a_i1, n, st);				    \
}

#define BINARY_MASKED_OP_EXEC(OP)					    \
template<int mode, typename T>						    \
void									    \
MemOpers<mode, T>::							    \
OP(T* a_o, const T* a_i, const T* a_i1, const T* a_mask,		    \
   size_t n, StreamT st)						    \
{									    \
    MemOpersExec::OP(a_o, a_i, a_i1, a_mask, n, st);			    \
}

#define BINARY_OP_I_EXEC(OP)						    \
template<int mode, typename T>                                              \
void 									    \
MemOpers<mode, T>::							    \
OP##_I(T* a_o, const T* a_i,						    \
		  size_t n, StreamT st)				            \
{									    \
    MemOpersExec::OP##_I(a_o, a_i, n, st);				    \
}

#define BINARY_MASKED_OP_I_EXEC(OP)					    \
template<int mode, typename T>						    \
void									    \
MemOpers<mode, T>::							    \
OP##_I(T* a_o, const T* a_i, const T* a_mask,				    \
       size_t n, StreamT st)						    \
{									    \
    MemOpersExec::OP##_I(a_o, a_i, a_mask, n, st);			    \
}

#define BINARY_OPC_EXEC(OP)						    \
template<int mode, typename T>                                              \
void 									    \
MemOpers<mode, T>::							    \
OP##C(T* a_o, const T* a_i, const T& c,					    \
		  size_t n, StreamT st,bool onDev)		            \
{									    \
    MemOpersExec::OP##C(a_o, a_i, c, n, st, onDev);			    \
}

#define BINARY_MASKED_OPC_EXEC(OP)					    \
template<int mode, typename T>						    \
void									    \
MemOpers<mode, T>::							    \
OP##C(T* a_o, const T* a_i, const T& c, const T* a_mask,		    \
      size_t n, StreamT st,bool onDev)					    \
{									    \
    MemOpersExec::OP##C(a_o, a_i, c, a_mask, n, st, onDev);		    \
}

#define BINARY_OP_NOC_EXEC(OP)						    \
template<int mode, typename T>                                              \
void 									    \
MemOpers<mode, T>::							    \
OP(T* a_o, const T* a_i, const T& c,					    \
		  size_t n, StreamT st,bool onDev)		            \
{									    \
    MemOpersExec::OP(a_o, a_i, c, n, st, onDev);			    \
}

#define BINARY_MASKED_OP_NOC_EXEC(OP)					    \
template<int mode, typename T>						    \
void									    \
MemOpers<mode, T>::							    \
OP(T* a_o, const T* a_i, const T& c, const T* a_mask,			    \
   size_t n, StreamT st,bool onDev)					    \
{									    \
    MemOpersExec::OP(a_o, a_i, c, a_mask, n, st, onDev);		    \
}

#define BINARY_OPC_I_EXEC(OP)						    \
template<int mode, typename T>                                              \
void 									    \
MemOpers<mode, T>::							    \
OP##C_I(T* a_o, const T& c,						    \
		   size_t n, StreamT st,bool onDev)		            \
{									    \
    MemOpersExec::OP##C_I(a_o, c, n, st, onDev);			    \
}

#define BINARY_MASKED_OPC_I_EXEC(OP)					    \
template<int mode, typename T>						    \
void									    \
MemOpers<mode, T>::							    \
OP##C_I(T* a_o, const T& c, const T* a_mask,				    \
	size_t n, StreamT st, bool onDev)				    \
{									    \
    MemOpersExec::OP##C_I(a_o, c, a_mask, n, st, onDev);		    \
}

#define BINARY_OP_NOC_I_EXEC(OP)					    \
template<int mode, typename T>                                              \
void 									    \
MemOpers<mode, T>::							    \
OP##_I(T* a_o, const T& c,						    \
		   size_t n, StreamT st,bool onDev)		            \
{									    \
    MemOpersExec::OP##_I(a_o, c, n, st, onDev);				    \
}

#define BINARY_MASKED_OP_NOC_I_EXEC(OP)				            \
template<int mode, typename T>						    \
void									    \
MemOpers<mode, T>::							    \
OP##_I(T* a_o, const T& c, const T* a_mask,				    \
       size_t n, StreamT st,bool onDev)					    \
{									    \
    MemOpersExec::OP##_I(a_o, c, a_mask, n, st, onDev);			    \
}

#define BINARY_OPS_EXEC(OP)						    \
BINARY_OP_EXEC(OP)							    \
BINARY_MASKED_OP_EXEC(OP)						    \
BINARY_OP_I_EXEC(OP)							    \
BINARY_MASKED_OP_I_EXEC(OP)						    \
BINARY_OPC_EXEC(OP)							    \
BINARY_MASKED_OPC_EXEC(OP)						    \
BINARY_OP_NOC_EXEC(OP)							    \
BINARY_MASKED_OP_NOC_EXEC(OP)						    \
BINARY_OPC_I_EXEC(OP)							    \
BINARY_MASKED_OPC_I_EXEC(OP)						    \
BINARY_OP_NOC_I_EXEC(OP)                                                    \
BINARY_MASKED_OP_NOC_I_EXEC(OP)

// Template version of MemOpers
template<int mode, typename T>
class MemOpers{
public:
    enum { exec_mode = mode };
    typedef typename Selector<mode, CMemOpers<T>, GMemOpers<T>, GMemOpers<T> >::Result MemOpersExec;

    static void Copy(T* a_o, const T* a_i, size_t n, StreamT st=NULL) {
        MemOpersExec::Copy(a_o, a_i, n, st);
    }

    static void Copy(T* a_o, const T* a_i, const T* a_mask,
		     size_t n, StreamT st=NULL)
    {
        MemOpersExec::Copy(a_o, a_i, a_mask, n, st);
    }

/**
 * Fill array with constant value
 */
    static void SetMem(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::SetMem(a_o, c, n, st, onDev);
    }

    static void SetMem(T* a_o, const T& c, const T* a_mask, size_t n, 
		       StreamT st=NULL,bool onDev=false){
        MemOpersExec::SetMem(a_o, c, a_mask, n, st, onDev);
    }

/**
 * Set unsined array with the linear ramp up
 * d[i] = i
 */
    static void SetLinear(T* a_o,  size_t n, StreamT st=NULL){
        MemOpersExec::SetLinear(a_o, n, st);
    }

/**
 * Set unsined array with the linear ramp down
 * d[i] = n - 1 - i;
 */
    static void SetLinearDown(T* a_o,  size_t n, StreamT st=NULL){
        MemOpersExec::SetLinearDown(a_o, n, st);
    }

/////////////////////////////////////////////////////////////////////////////////
// Unary functions
//////////////////////////////////////////////////////////////////////////////////
/**
 * Absolute value of array (a_o = abs(a_i))
 */
    static void Abs(T* a_o, const T* a_i, size_t n, StreamT st=NULL);
    static void Abs(T* a_o, const T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * Absolute value of array, in-place (a_o = abs(a_o))
 */
    static void Abs_I(T* a_o, size_t n, StreamT st=NULL);
    static void Abs_I(T* a_o, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * Cube each value in array (a_o = a_i^3)
 */
    static void Cube(T* a_o, const T* a_i, size_t n, StreamT st=NULL);
    static void Cube(T* a_o, const T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * Cube each value in array, in-place (a_o = a_o^3)
 */
    static void Cube_I(T* a_o, size_t n, StreamT st=NULL);
    static void Cube_I(T* a_o, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * Exponential of each value in array (a_o = exp(a_i))
 */
    static void Exp(T* a_o, const  T* a_i, size_t n, StreamT st=NULL);
    static void Exp(T* a_o, const  T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * Exponential of each value in array (a_o = exp(a_o)) inplace version
 */
    static void Exp_I(T* a_o, size_t n, StreamT st=NULL);
    static void Exp_I(T* a_o, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * Logarithm of each value in array (a_o = log(a_i))
 */
    static void Log(T* a_o, const  T* a_i, size_t n, StreamT st=NULL);
    static void Log(T* a_o, const  T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * Logarithm of each value in array (a_o = log(a_o)) inplace version
 */
    static void Log_I(T* a_o, size_t n, StreamT st=NULL);
    static void Log_I(T* a_o, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * Negate each value in array (a_o = -a_i)
 */
    static void Neg(T* a_o, const T* a_i, size_t n, StreamT st=NULL);
    static void Neg(T* a_o, const T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * Negate each value in array, in-place (a_o = -a_o)
 */
    static void Neg_I(T* a_o, size_t n, StreamT st=NULL);
    static void Neg_I(T* a_o, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * Ramp function of each value in array (a_o = ramp(a_i))
 */
    static void Ramp(T* a_o, const T* a_i, size_t n, StreamT st=NULL);
    static void Ramp(T* a_o, const T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * Ramp function of each value in array in-place (a_o = ramp(a_o))
 */
    static void Ramp_I(T* a_o, size_t n, StreamT st=NULL);
    static void Ramp_I(T* a_o, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * Sign of each value in array (a_o = sgn(a_i))
 */
    static void Sgn(T* a_o, const T* a_i, size_t n, StreamT st=NULL);
    static void Sgn(T* a_o, const T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * Sign of each value in array in-place (a_o = sgn(a_o))
 */
    static void Sgn_I(T* a_o, size_t n, StreamT st=NULL);
    static void Sgn_I(T* a_o, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * Step function of each value in array (a_o = step(a_i))
 */
    static void Step(T* a_o, const T* a_i, size_t n, StreamT st=NULL);
    static void Step(T* a_o, const T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * Sign of each value in array in-place (a_o = step(a_o))
 */
    static void Step_I(T* a_o, size_t n, StreamT st=NULL);
    static void Step_I(T* a_o, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * Square of each value in array (a_o = sqr(a_i))
 */
    static void Sqr(T* a_o, const T* a_i, size_t n, StreamT st=NULL);
    static void Sqr(T* a_o, const T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * Square of each value in array, in-place (a_o = sqr(a_o))
 */
    static void Sqr_I(T* a_o, size_t n, StreamT st=NULL);
    static void Sqr_I(T* a_o, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * Square root of each value in array (a_o = sqrt(a_i))
 */
    static void Sqrt(T* a_o, const  T* a_i, size_t n, StreamT st=NULL);
    static void Sqrt(T* a_o, const  T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * Square root of each value in array (a_o = sqrt(a_o)) inplace version
 */
    static void Sqrt_I(T* a_o, size_t n, StreamT st=NULL);
    static void Sqrt_I(T* a_o, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * Inverse of each value in array (a_o = 1.0/a_i)
 */
    static void Inv(T* a_o, const T* a_i, size_t n, StreamT st=NULL);
    static void Inv(T* a_o, const T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * Inverse of each value in array (a_o = 1.0/a_o)
 */
    static void Inv_I(T* a_o, size_t n, StreamT st=NULL);
    static void Inv_I(T* a_o, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * floor
 */
    static void Floor(T* a_o, const T* a_i, size_t n, StreamT st=NULL);
    static void Floor(T* a_o, const T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * floor, inplace version
 */
    static void Floor_I(T* a_o, size_t n, StreamT st=NULL);
    static void Floor_I(T* a_o, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * ceil
 */
    static void Ceil(T* a_o, const T* a_i, size_t n, StreamT st=NULL);
    static void Ceil(T* a_o, const T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * ceil, inplace version
 */
    static void Ceil_I(T* a_o, size_t n, StreamT st=NULL);
    static void Ceil_I(T* a_o, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * round
 */
    static void Round(T* a_o, const T* a_i, size_t n, StreamT st=NULL);
    static void Round(T* a_o, const T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * round, inplace version
 */
    static void Round_I(T* a_o, size_t n, StreamT st=NULL);
    static void Round_I(T* a_o, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * sine in radians of each value in array (a_o = sin(a_i))
 */
    static void Sin(T* a_o, const  T* a_i, size_t n, StreamT st=NULL);
    static void Sin(T* a_o, const  T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * sine in radians of each value in array (a_o = sin(a_o)) inplace version
 */
    static void Sin_I(T* a_o, size_t n, StreamT st=NULL);
    static void Sin_I(T* a_o, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * arcsine in radians of each value in array (h_o = asin(h_i))
 */
    static void Asin(T* a_o, const  T* a_i, size_t n, StreamT st=NULL);
    static void Asin(T* a_o, const  T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * arcsine in radians of each value in array (h_o = asin(h_o)) inplace version
 */
    static void Asin_I(T* a_o, size_t n, StreamT st=NULL);
    static void Asin_I(T* a_o, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * cosine in radians of each value in array (h_o = cos(h_i))
 */
    static void Cos(T* a_o, const  T* a_i, size_t n, StreamT st=NULL);
    static void Cos(T* a_o, const  T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * cosine in radians of each value in array (h_o = cos(h_o)) inplace version
 */
    static void Cos_I(T* a_o, size_t n, StreamT st=NULL);
    static void Cos_I(T* a_o, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * arccosine in radians of each value in array (h_o = acos(h_i))
 */
    static void Acos(T* a_o, const  T* a_i, size_t n, StreamT st=NULL);
    static void Acos(T* a_o, const  T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * arccosine in radians of each value in array (h_o = acos(h_o)) inplace version
 */
    static void Acos_I(T* a_o, size_t n, StreamT st=NULL);
    static void Acos_I(T* a_o, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * tangent in radians of each value in array (h_o = tan(h_i))
 */
    static void Tan(T* a_o, const  T* a_i, size_t n, StreamT st=NULL);
    static void Tan(T* a_o, const  T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * tangent in radians of each value in array (h_o = tan(h_o)) inplace version
 */
    static void Tan_I(T* a_o, size_t n, StreamT st=NULL);
    static void Tan_I(T* a_o, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * arctangent in radians of each value in array (h_o = atan(h_i))
 */
    static void Atan(T* a_o, const  T* a_i, size_t n, StreamT st=NULL);
    static void Atan(T* a_o, const  T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * arctangent in radians of each value in array (h_o = atan(h_o)) inplace version
 */
    static void Atan_I(T* a_o, size_t n, StreamT st=NULL);
    static void Atan_I(T* a_o, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * cosecant in radians of each value in array (h_o = csc(h_i))
 */
    static void Csc(T* a_o, const  T* a_i, size_t n, StreamT st=NULL);
    static void Csc(T* a_o, const  T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * cosecant in radians of each value in array (h_o = csc(h_o)) inplace version
 */
    static void Csc_I(T* a_o, size_t n, StreamT st=NULL);
    static void Csc_I(T* a_o, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * secant in radians of each value in array (h_o = sec(h_i))
 */
    static void Sec(T* a_o, const  T* a_i, size_t n, StreamT st=NULL);
    static void Sec(T* a_o, const  T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * secant in radians of each value in array (h_o = sec(h_o)) inplace version
 */
    static void Sec_I(T* a_o, size_t n, StreamT st=NULL);
    static void Sec_I(T* a_o, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * cotangent in radians of each value in array (h_o = cot(h_i))
 */
    static void Cot(T* a_o, const  T* a_i, size_t n, StreamT st=NULL);
    static void Cot(T* a_o, const  T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * cotangent in radians of each value in array (h_o = cot(h_o)) inplace version
 */
    static void Cot_I(T* a_o, size_t n, StreamT st=NULL);
    static void Cot_I(T* a_o, const T* a_mask, size_t n, StreamT st=NULL);


/////////////////////////////////////////////////////////////////////////////////
// Binary functions
//////////////////////////////////////////////////////////////////////////////////
/**
 * Add constant to an array (a_o = a_i + c)
 */
    static void AddC(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void AddC(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void Add(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void Add(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

/**
 * Add constant to an array, in-place (a_o = a_o + c)
 */
    static void AddC_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void AddC_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void Add_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void Add_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

/**
 * Add two arrays, (a_o = a_i + a_i1)
 */
    static void Add(T* a_o, const T* a_i, const T* a_i1, size_t n, StreamT st=NULL);
    static void Add(T* a_o, const T* a_i, const T* a_i1, const T* a_mask, size_t n, StreamT st=NULL);

/**
 * Add two arrays, in-place (a_o = a_o + a_i)
 */
    static void Add_I(T* a_o, const T* a_i, size_t n, StreamT st=NULL);
    static void Add_I(T* a_o, const T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

    /**
     * Sub constant from an array (a_o = a_i - c)
     */
    static void SubC(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void SubC(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void Sub(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void Sub(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

/**
 * Sub constant from array, in-place (a_o = a_o - c)
 */
    static void SubC_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void SubC_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void Sub_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void Sub_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

/**
 * Sub two arrays, (a_o = a_i - a_i1)
 */
    static void Sub(T* a_o, const T* a_i, const T* a_i1, size_t n, StreamT st=NULL);
    static void Sub(T* a_o, const T* a_i, const T* a_i1, const T* a_mask, size_t n, StreamT st=NULL);
/**
 * Sub two arrays, in-place (a_o = a_o - a_i)
 */
    static void Sub_I(T* a_o, const T* a_i, size_t n, StreamT st=NULL);
    static void Sub_I(T* a_o, const T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

    /**
     * Mul constant from an array (a_o = a_i * c)
     */
    static void MulC(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void MulC(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void Mul(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void Mul(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

/**
 * Mul constant from array, in-place (a_o = a_o * c)
 */
    static void MulC_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void MulC_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void Mul_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void Mul_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

/**
 * Mul two arrays, (a_o = a_i * a_i1)
 */
    static void Mul(T* a_o, const T* a_i, const T* a_i1, size_t n, StreamT st=NULL);
    static void Mul(T* a_o, const T* a_i, const T* a_i1, const T* a_mask, size_t n, StreamT st=NULL);

    /**
     * Mul two arrays, in-place (a_o = a_o * a_i)
     */
    static void Mul_I(T* a_o, const T* a_i, size_t n, StreamT st=NULL);
    static void Mul_I(T* a_o, const T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

    /**
     * Div constant from an array (a_o = a_i / c)
     */
    static void DivC(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void DivC(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void Div(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void Div(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

/**
 * Div constant from array, in-place (a_o = a_o / c)
 */
    static void DivC_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void DivC_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void Div_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void Div_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

/**
 * Div two arrays, (a_o = a_i / a_i1)
 */
    static void Div(T* a_o, const T* a_i, const T* a_i1, size_t n, StreamT st=NULL);
    static void Div(T* a_o, const T* a_i, const T* a_i1, const T* a_mask, size_t n, StreamT st=NULL);

    /**
     * Div two arrays, in-place (a_o = a_o / a_i)
     */
    static void Div_I(T* a_o, const T* a_i, size_t n, StreamT st=NULL);
    static void Div_I(T* a_o, const T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

    /**
     * max of constant and array (a_o = max(a_i, c))
     */
    static void MaxC(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void MaxC(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void Max(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void Max(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

/**
 * max of constant and array, in-place (a_o = max(a_o, c))
 */
    static void MaxC_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void MaxC_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void Max_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void Max_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

/**
 * max of two arrays, (a_o = max(a_i, a_i1))
 */
    static void Max(T* a_o, const T* a_i, const T* a_i1, size_t n, StreamT st=NULL);
    static void Max(T* a_o, const T* a_i, const T* a_i1, const T* a_mask, size_t n, StreamT st=NULL);

    /**
     * max of two arrays, in-place (a_o = max(a_o, a_i))
     */
    static void Max_I(T* a_o, const T* a_i, size_t n, StreamT st=NULL);
    static void Max_I(T* a_o, const T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

    /**
     * max of constant and array (a_o = max(a_i, c))
     */
    static void MinC(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void MinC(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void Min(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void Min(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

/**
 * max of constant and array, in-place (a_o = max(a_o, c))
 */
    static void MinC_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void MinC_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void Min_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void Min_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

/**
 * max of two arrays, (a_o = max(a_i, a_i1))
 */
    static void Min(T* a_o, const T* a_i, const T* a_i1, size_t n, StreamT st=NULL);
    static void Min(T* a_o, const T* a_i, const T* a_i1, const T* a_mask, size_t n, StreamT st=NULL);

    /**
     * max of two arrays, in-place (a_o = max(a_o, a_i))
     */
    static void Min_I(T* a_o, const T* a_i, size_t n, StreamT st=NULL);
    static void Min_I(T* a_o, const T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

   /**
    * pow (exponentiation) ho(x) = hi(x)^c
    */
    static void PowC(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void PowC(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void Pow(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void Pow(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

   /**
    * pow (exponentiation), in-place (ho(x) = ho(x)^c)
    */
    static void PowC_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void PowC_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void Pow_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void Pow_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

    /**
     * arctangent given by y coordinate and x coordinate (a_o = a_i * a_i1)
     */
    static void Atan2(T* a_o, const T* a_i, const T* a_i1, size_t n, StreamT st=NULL);
    static void Atan2(T* a_o, const T* a_i, const T* a_i1, const T* a_mask, size_t n, StreamT st=NULL);

    /**
     * arctangent given by y coordinate and x coordinate, inplace (a_o = atan(a_o, a_i)
     */
    static void Atan2_I(T* a_o, const T* a_i, size_t n, StreamT st=NULL);
    static void Atan2_I(T* a_o, const T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

    /**
     * arctangent given by y coordinate and x coordinate (a_o = atan2(a_i, c))
     */
    static void Atan2C(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void Atan2C(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void Atan2(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void Atan2(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

    /**
     * arctangent given by y coordinate and x coordinate, inplace (a_o = atan(a_o, c)
     */
    static void Atan2C_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void Atan2C_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void Atan2_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void Atan2_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

    /**
     * Greater than of an array and constant (h_o = h_i > c)
     */
    static void GTC(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void GTC(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void GT(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void GT(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

    /**
     * Greater than of an array and constant, in-place (h_o = h_o > c)
     */
    static void GTC_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void GTC_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void GT_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void GT_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

    /**
     * Greater than of two arrays, (h_o = h_i > h_i1)
     */
    static void GT(T* a_o, const T* a_i, const T* a_i1, size_t n, StreamT st=NULL);
    static void GT(T* a_o, const T* a_i, const T* a_i1, const T* a_mask, size_t n, StreamT st=NULL);

    /**
     * Greater than of two arrays, in-place (h_o = h_o > h_i)
     */
    static void GT_I(T* a_o, const T* a_i, size_t n, StreamT st=NULL);
    static void GT_I(T* a_o, const T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

    /**
     * Greater than or equal of an array and constant (h_o = h_i >= c)
     */
    static void GTEC(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void GTEC(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void GTE(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void GTE(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

    /**
     * Greater than or equal of an array and constant, in-place (h_o = h_o >= c)
     */
    static void GTE_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void GTE_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void GTEC_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void GTEC_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

    /**
     * Greater than or equal of two arrays, (h_o = h_i >= h_i1)
     */
    static void GTE(T* a_o, const T* a_i, const T* a_i1, size_t n, StreamT st=NULL);
    static void GTE(T* a_o, const T* a_i, const T* a_i1, const T* a_mask, size_t n, StreamT st=NULL);

    /**
     * Greater than or equal of two arrays, in-place (h_o = h_o >= h_i)
     */
    static void GTE_I(T* a_o, const T* a_i, size_t n, StreamT st=NULL);
    static void GTE_I(T* a_o, const T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

    /**
     * Boolean Equal of an array and constant (h_o = h_i == c)
     */
    static void EQC(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void EQC(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void EQ(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void EQ(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

    /**
     * Boolean Equal of an array and constant, in-place (h_o = h_o == c)
     */
    static void EQC_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void EQC_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void EQ_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void EQ_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

    /**
     * Boolean Equal of two arrays, (h_o = h_i == h_i1)
     */
    static void EQ(T* a_o, const T* a_i, const T* a_i1, size_t n, StreamT st=NULL);
    static void EQ(T* a_o, const T* a_i, const T* a_i1, const T* a_mask, size_t n, StreamT st=NULL);

    /**
     * Boolean Equal of two arrays, in-place (h_o = h_o == h_i)
     */
    static void EQ_I(T* a_o, const T* a_i, size_t n, StreamT st=NULL);
    static void EQ_I(T* a_o, const T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

    /**
     * Boolean Not Equal of an array and constant (h_o = h_i != c)
     */
    static void NEQC(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void NEQC(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void NEQ(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void NEQ(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

    /**
     * Boolean Not Equal of an array and constant, in-place (h_o = h_o != c)
     */
    static void NEQC_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void NEQC_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void NEQ_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void NEQ_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

    /**
     * Boolean Not Equal of two arrays, (h_o = h_i != h_i1)
     */
    static void NEQ(T* a_o, const T* a_i, const T* a_i1, size_t n, StreamT st=NULL);
    static void NEQ(T* a_o, const T* a_i, const T* a_i1, const T* a_mask, size_t n, StreamT st=NULL);

    /**
     * Boolean Not Equal of two arrays, in-place (h_o = h_o != h_i)
     */
    static void NEQ_I(T* a_o, const T* a_i, size_t n, StreamT st=NULL);
    static void NEQ_I(T* a_o, const T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

    /**
     * Less than of an array and constant (h_o = h_i < c)
     */
    static void LTC(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void LTC(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void LT(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void LT(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

    /**
     * Less than of an array and constant, in-place (h_o = h_o < c)
     */
    static void LTC_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void LTC_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void LT_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void LT_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

    /**
     * Less than of two arrays, (h_o = h_i < h_i1)
     */
    static void LT(T* a_o, const T* a_i, const T* a_i1, size_t n, StreamT st=NULL);
    static void LT(T* a_o, const T* a_i, const T* a_i1, const T* a_mask, size_t n, StreamT st=NULL);

    /**
     * Less than of two arrays, in-place (h_o = h_o < h_i)
     */
    static void LT_I(T* a_o, const T* a_i, size_t n, StreamT st=NULL);
    static void LT_I(T* a_o, const T* a_i, const T* a_mask, size_t n, StreamT st=NULL);

    /**
     * Less than or equal of an array and constant (h_o = h_i <= c)
     */
    static void LTEC(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void LTEC(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void LTE(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void LTE(T* a_o, const T* a_i, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

    /**
     * Less than or equal of an array and constant, in-place (h_o = h_o <= c)
     */
    static void LTEC_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void LTEC_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);
    static void LTE_I(T* a_o, const T& c, size_t n, StreamT st=NULL,bool onDev=false);
    static void LTE_I(T* a_o, const T& c, const T* a_mask, size_t n, StreamT st=NULL,bool onDev=false);

    /**
     * Less than or equal of two arrays, (h_o = h_i <= h_i1)
     */
    static void LTE(T* a_o, const T* a_i, const T* a_i1, size_t n, StreamT st=NULL);
    static void LTE(T* a_o, const T* a_i, const T* a_i1, const T* a_mask, size_t n, StreamT st=NULL);

    /**
     * Less than or equal of two arrays, in-place (h_o = h_o <= h_i)
     */
    static void LTE_I(T* a_o, const T* a_i, size_t n, StreamT st=NULL);
    static void LTE_I(T* a_o, const T* a_i, const T* a_mask, size_t n, StreamT st=NULL);


/////////////////////////////////////////////////////////////////////////////////
// Triary functions
//////////////////////////////////////////////////////////////////////////////////

/**
 * Compute the absolution different between 2 vector
 * a_o = abs(a_i - a_i1)
 */
    static void AbsDiff(T* a_o, const T* a_i, const T* a_i1, size_t n, StreamT st=NULL){
        MemOpersExec::AbsDiff(a_o, a_i, a_i1, n, st);
    }

/**
 * Compute the absolution different between 2 vector in-place version
 * a_o = abs(a_o - a_i1)
 */
    static void AbsDiff_I(T* a_o, const T* a_i, size_t n, StreamT st=NULL){
        MemOpersExec::AbsDiff_I(a_o, a_i, n, st);
    }

/**
 * Compute the absolution different between 2 vector
 * a_o = sqr(a_i - a_i1)
 */

    static void SqrDiff(T* a_o, const T* a_i, const T* a_i1, size_t n, StreamT st=NULL){
        MemOpersExec::SqrDiff(a_o, a_i, a_i1, n, st);
    }

/**
 * Compute the absolution different between 2 vector in-place version
 * a_o = sqr(a_o - a_i1)
 */

    static void SqrDiff_I(T* a_o, const T* a_i, size_t n, StreamT st=NULL){
        MemOpersExec::SqrDiff_I(a_o, a_i, n, st);
    }

/**
 * a_o = a_i * (a_i1 * c)
 */

    static void MulMulC(T* a_o, const T* a_i, const T* a_i1, const T& c, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::MulMulC(a_o, a_i, a_i1, c, n, st, onDev);
    }

/**
 * a_o = a_o * (a_i1 * c)
 */

    static void MulMulC_I(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::MulMulC_I(a_o, a_i, c, n, st, onDev);
    }

/**
 * a_o = a_i * (a_i1 * a_i2)
 */

    static void MulMul(T* a_o, const T* a_i, const T* a_i1, const T* a_i2, size_t n, StreamT st=NULL){
        MemOpersExec::MulMul(a_o, a_i, a_i1, a_i2, n, st);
    }

/**
 * a_o = a_o * (a_i * a_i1)
 */

    static void MulMul_I(T* a_o, const T* a_i, const T* a_i1, size_t n, StreamT st=NULL){
        MemOpersExec::MulMul_I(a_o, a_i, a_i1, n, st);
    }

/**
 * a_o = a_i * a_i1 + a_i2
 */

    static void MulAdd(T* a_o, const T* a_i, const T* a_i1, const T* a_i2, size_t n, StreamT st=NULL){
        MemOpersExec::MulAdd(a_o, a_i, a_i1, a_i2, n, st);
    }

/**
 * a_o = a_o * a_i + a_i1
 */

    static void MulAdd_I(T* a_o, const T* a_i, const T* a_i1, size_t n, StreamT st=NULL){
        MemOpersExec::MulAdd_I(a_o, a_i, a_i1, n, st);
    }

/**
 * a_o = a_i * a_i1 - a_i2
 */

    static void MulSub(T* a_o, const T* a_i, const T* a_i1, const T* a_i2, size_t n, StreamT st=NULL){
        MemOpersExec::MulSub(a_o, a_i, a_i1, a_i2, n, st);
    }

/**
 * a_o = a_o * a_i - a_i1
 */

    static void MulSub_I(T* a_o, const T* a_i, const T* a_i1, size_t n, StreamT st=NULL){
        MemOpersExec::MulSub_I(a_o, a_i, a_i1, n, st);
    }


/**
 * a_o = (a_i + a_i1) * a_i2
 */

    static void AddMul(T* a_o, const T* a_i, const T* a_i1, const T* a_i2, size_t n, StreamT st=NULL){
        MemOpersExec::AddMul(a_o, a_i, a_i1, a_i2, n, st);
    }

/**
 * a_o = (a_o + a_i) * a_i1
 */

    static void AddMul_I(T* a_o, const T* a_i, const T* a_i1, size_t n, StreamT st=NULL){
        MemOpersExec::AddMul_I(a_o, a_i, a_i1, n, st);
    }

/**
 * a_o = (a_i + a_i1) / a_i2
 */

    static void AddDiv(T* a_o, const T* a_i, const T* a_i1, const T* a_i2, size_t n, StreamT st=NULL){
        MemOpersExec::AddDiv(a_o, a_i, a_i1, a_i2, n, st);
    }

/**
 * a_o = (a_o + a_i) / a_i1
 */

    static void AddDiv_I(T* a_o, const T* a_i, const T* a_i1, size_t n, StreamT st=NULL){
        MemOpersExec::AddDiv_I(a_o, a_i, a_i1, n, st);
    }

/**
 * a_o = (a_i - a_i1) * a_i2
 */

    static void SubMul(T* a_o, const T* a_i, const T* a_i1, const T* a_i2, size_t n, StreamT st=NULL){
        MemOpersExec::SubMul(a_o, a_i, a_i1, a_i2, n, st);
    }

/**
 * a_o = (a_o - a_i) * a_i1
 */

    static void SubMul_I(T* a_o, const T* a_i, const T* a_i1, size_t n, StreamT st=NULL){
        MemOpersExec::SubMul_I(a_o, a_i, a_i1, n, st);

    }

/**
 * a_o = (a_i - a_i1) / a_i2
 */

    static void SubDiv(T* a_o, const T* a_i, const T* a_i1, const T* a_i2, size_t n, StreamT st=NULL){
        MemOpersExec::SubDiv(a_o, a_i, a_i1, a_i2, n, st);
    }

/**
 * a_o = (a_o - a_i) / a_i1
 */

    static void SubDiv_I(T* a_o, const T* a_i, const T* a_i1, size_t n, StreamT st=NULL){
        MemOpersExec::SubDiv_I(a_o, a_i, a_i1, n, st);
    }

/**
 * a_o = (a_i * c) + a_i1
 */

    static void MulCAdd(T* a_o, const T* a_i, const T& c, const T* a_i1, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::MulCAdd(a_o, a_i, c, a_i1, n, st, onDev);
    }

/**
 * a_o = (a_o * c) + a_i
 */

    static void MulCAdd_I(T* a_o, const T& c, const T* a_i, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::MulCAdd_I(a_o, c, a_i, n, st, onDev);
    }

/**
 * a_o = (a_i * c) - a_i1
 */

    static void MulCSub(T* a_o, const T* a_i, const T& c, const T* a_i1, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::MulCSub(a_o, a_i, c, a_i1, n, st, onDev);
    }

/**
 * a_o = (a_o * c) - a_i
 */

    static void MulCSub_I(T* a_o, const T& c, const T* a_i, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::MulCSub_I(a_o, c, a_i, n, st, onDev);
    }

/**
 * a_o = (a_i * a) + b
 */

    static void MulCAddC(T* a_o, const T* a_i, const T& a, const T& b, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::MulCAddC(a_o, a_i, a, b, n, st, onDev);
    }

/**
 * a_o = (a_o * a) + b
 */
    static void MulCAddC_I(T* a_o, const T& a, const T& b, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::MulCAddC_I(a_o, a, b, n, st, onDev);
    }

/**
 * a_o = (a_i + a_i1) * c
 */
    static void AddMulC(T* a_o, const T* a_i, const T* a_i1, const T& c, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::AddMulC(a_o, a_i, a_i1, c, n, st, onDev);
    }

/**
 * a_o = (a_o + a_i) * c
 */
    static void AddMulC_I(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::AddMulC_I(a_o, a_i, c, n, st, onDev);
    }

/**
 * a_o = (a_i - a_i1) * c
 */
    static void SubMulC(T* a_o, const T* a_i, const T* a_i1, const T& c, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::SubMulC(a_o, a_i, a_i1, c, n, st, onDev);
    }

/**
 * a_o = (a_o - a_i) * c
 */
    static void SubMulC_I(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::SubMulC_I(a_o, a_i, c, n, st, onDev);
    }


/**
 * a_o = (a_i + a) * b
 */
    static void AddCMulC(T* a_o, const T* a_i, const T& a, const T& b, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::AddCMulC(a_o, a_i, a, b, n, st, onDev);
    }

/**
 * a_o = (a_o + a) * b
 */
    static void AddCMulC_I(T* a_o, const T& a, const T& b, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::AddCMulC_I(a_o, a, b, n, st, onDev);
    }

/**
 * a_o = a_i + a_i1 * c
 */
    static void Add_MulC(T* a_o, const T* a_i, const T* a_i1, const T& c, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::Add_MulC(a_o, a_i, a_i1, c, n, st, onDev);
    }

/**
 * a_o = a_o + a_i * c
 */
    static void Add_MulC_I(T* a_o, const T* a_i, const T& c, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::Add_MulC_I(a_o, a_i, c, n, st, onDev);
    }

/**
 * a_o = a_i + a_i1 * a_i2
 */
    static void Add_Mul(T* a_o, const T* a_i, const T* a_i1, const T* a_i2, size_t n, StreamT st=NULL)
        {
            MemOpersExec::Add_Mul(a_o, a_i, a_i1, a_i2, n, st);
        }
/**
 * a_o = a_o + a_i * a_i1
 */
    static void Add_Mul_I(T* a_o, const T* a_i, const T* a_i1, size_t n, StreamT st=NULL){
        MemOpersExec::Add_Mul_I(a_o, a_i, a_i1, n, st);
    }

    /**
     * a_o = a_i - a_i1 * a_i2
     */
    static void Sub_Mul(T* a_o, const T* a_i, const T* a_i1, const T* a_i2, size_t n, StreamT st=NULL){
        MemOpersExec::Sub_Mul(a_o, a_i, a_i1, a_i2, n, st);
    }
/**
 * a_o = a_o - a_i * a_i1
 */
    static void Sub_Mul_I(T* a_o, const T* a_i, const T* a_i1, size_t n, StreamT st=NULL){
        MemOpersExec::Sub_Mul_I(a_o, a_i, a_i1, n, st);
    }

/////////////////////////////////////////////////////////////////////////////////
// n-ary functions
//////////////////////////////////////////////////////////////////////////////////
/**
 * a_o = a_i + (a_i1 + a_i2) * c
 */
    static void Add_AddMulC(T* a_o, const T* a_i, const T* a_i1, const T* a_i2, const T& c, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::Add_AddMulC(a_o, a_i, a_i1, a_i2, c, n, st, onDev);
    }

/**
 * a_o = a_i + (a_i1 - a_i2) * c
 */
    static void Add_SubMulC(T* a_o, const T* a_i, const T* a_i1, const T* a_i2, const T& c, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::Add_SubMulC(a_o, a_i, a_i1, a_i2, c, n, st, onDev);
    }

/**
 * a_o = a_i + (a_i1 * a_i2) * c
 */
    static void Add_MulMulC(T* a_o, const T* a_i, const T* a_i1, const T* a_i2, const T& c, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::Add_MulMulC(a_o, a_i, a_i1, a_i2, c, n, st, onDev);
    }


/**
 * a_o += (a_i + a_i1) * c
 */
    static void Add_AddMulC_I(T* a_o, const T* a_i, const T* a_i1, const T& c, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::Add_AddMulC_I(a_o, a_i, a_i1, c, n, st, onDev);
    }

/**
 * a_o += (a_i - a_i1) * c
 */
    static void Add_SubMulC_I(T* a_o, const T* a_i, const T* a_i1, const T& c, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::Add_SubMulC_I(a_o, a_i, a_i1, c, n, st, onDev);
    }

/**
 * a_o += (a_i * a_i1) * c
 */
    static void Add_MulMulC_I(T* a_o, const T* a_i, const T* a_i1, const T& c, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::Add_MulMulC_I(a_o, a_i, a_i1, c, n, st, onDev);
    }

    /**
     * a_o = a_i - (a_i1 + a_i2) * c
     */
    static void Sub_AddMulC(T* a_o, const T* a_i, const T* a_i1, const T* a_i2, const T& c, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::Sub_AddMulC(a_o, a_i, a_i1, a_i2, c, n, st, onDev);
    }

/**
 * a_o = a_i - (a_i1 - a_i2) * c
 */
    static void Sub_SubMulC(T* a_o, const T* a_i, const T* a_i1, const T* a_i2, const T& c, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::Sub_SubMulC(a_o, a_i, a_i1, a_i2, c, n, st, onDev);
    }

/**
 * a_o = a_i - (a_i1 * a_i2) * c
 */
    static void Sub_MulMulC(T* a_o, const T* a_i, const T* a_i1, const T* a_i2, const T& c, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::Sub_MulMulC(a_o, a_i, a_i1, a_i2, c, n, st, onDev);
    }


/**
 * a_o -= (a_i + a_i1) * c
 */
    static void Sub_AddMulC_I(T* a_o, const T* a_i, const T* a_i1, const T& c, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::Sub_AddMulC_I(a_o, a_i, a_i1, c, n, st, onDev);
    }

/**
 * a_o -= (a_i - a_i1) * c
 */
    static void Sub_SubMulC_I(T* a_o, const T* a_i, const T* a_i1, const T& c, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::Sub_SubMulC_I(a_o, a_i, a_i1, c, n, st, onDev);
    }

/**
 * a_o -= (a_i * a_i1) * c
 */
    static void Sub_MulMulC_I(T* a_o, const T* a_i, const T* a_i1, const T& c, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::Sub_MulMulC_I(a_o, a_i, a_i1, c, n, st, onDev);
    }


// Functions with 4 inputs

/**
 * a_o = a_i * a + a_i1 * b
 */
    static void MulC_Add_MulC(T* a_o, const T* a_i, const T& a, const T* a_i1, const T& b, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::MulC_Add_MulC(a_o, a_i, a, a_i1, b, n, st, onDev);
    }

/**
 * a_o = a_o * a + a_i1 * b
 */
    static void MulC_Add_MulC_I(T* a_o, const T& a, const T* a_i, const T& b, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::MulC_Add_MulC_I(a_o, a, a_i, b, n, st, onDev);
    }

/**
 * a_o = (a_i + a) * b + c
 */
    static void AddCMulCAddC(T* a_o, const T* a_i, const T& a, const T& b, const T& c, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::AddCMulCAddC(a_o, a_i, a, b, c, n, st, onDev);
    }

/**
 * a_o = (a_o + a) * b + c
 */
    static void AddCMulCAddC_I(T* a_o, const T& a, const T& b, const T& c, size_t n, StreamT st=NULL,bool onDev=false){
        MemOpersExec::AddCMulCAddC_I(a_o, a, b, c, n, st, onDev);
    }

/**
 * a_o[i] = a_i[n-1 - i]
 */
    static void ReverseOrder(T* a_o, const T* a_i, size_t n, StreamT st=NULL){
        MemOpersExec::ReverseOrder(a_o, a_i, n, st);
    }

/**
 * Shift coordinates, so that x becomes y axis, y becomes z, and z becomes x (opposite direction if dir is false)
 */
    static void ShiftCoordinate(T* a_o, const T* a_i, size_t sizeX, size_t sizeY, size_t sizeZ, bool dir, StreamT st=NULL){
        MemOpersExec::ShiftCoordinate(a_o, a_i, sizeX, sizeY, sizeZ, dir, st);
    }
};

UNARY_OPS_EXEC(Abs)
UNARY_OPS_EXEC(Cube)
UNARY_OPS_EXEC(Exp)
UNARY_OPS_EXEC(Log)
UNARY_OPS_EXEC(Sqr)
UNARY_OPS_EXEC(Neg)
UNARY_OPS_EXEC(Ramp)
UNARY_OPS_EXEC(Sgn)
UNARY_OPS_EXEC(Step)
UNARY_OPS_EXEC(Sqrt)
UNARY_OPS_EXEC(Inv)
UNARY_OPS_EXEC(Sin)
UNARY_OPS_EXEC(Asin)
UNARY_OPS_EXEC(Cos)
UNARY_OPS_EXEC(Acos)
UNARY_OPS_EXEC(Tan)
UNARY_OPS_EXEC(Atan)
UNARY_OPS_EXEC(Csc)
UNARY_OPS_EXEC(Sec)
UNARY_OPS_EXEC(Cot)
UNARY_OPS_EXEC(Ceil)
UNARY_OPS_EXEC(Floor)
UNARY_OPS_EXEC(Round)

BINARY_OPS_EXEC(Add)
BINARY_OPS_EXEC(Sub)
BINARY_OPS_EXEC(Mul)
BINARY_OPS_EXEC(Div)
BINARY_OPS_EXEC(Max)
BINARY_OPS_EXEC(Min)
BINARY_OPS_EXEC(LT)
BINARY_OPS_EXEC(LTE)
BINARY_OPS_EXEC(EQ)
BINARY_OPS_EXEC(NEQ)
BINARY_OPS_EXEC(GT)
BINARY_OPS_EXEC(GTE)
BINARY_OPS_EXEC(Atan2)

// Pow doesn't have image version
BINARY_OPC_EXEC(Pow)
BINARY_OPC_I_EXEC(Pow)
BINARY_OP_NOC_EXEC(Pow)
BINARY_OP_NOC_I_EXEC(Pow)

} // end namespace PyCA

#endif

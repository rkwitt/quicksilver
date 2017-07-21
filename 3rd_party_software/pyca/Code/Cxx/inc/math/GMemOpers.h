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

#ifndef  __GPU_MEM_OPERATION__H
#define  __GPU_MEM_OPERATION__H

#include <pycaConst.h>
#include <MOpers.h>
#include <estream.h>

namespace PyCA {

template<typename T>
class GMemOpers{
private:

    template<MATH_OPS op>
    static void Comp_unary(T* d_o, const T* d_i, size_t n, StreamT st);
    template<MATH_OPS op>
    static void Comp_unary(T* d_o, const T* d_i, const T* d_mask, 
			   size_t n, StreamT st);

    template<MATH_OPS op>
    static void Comp_unary_I(T* d_o, size_t n, StreamT st);
    template<MATH_OPS op>
    static void Comp_unary_I(T* d_o, const T* d_mask, size_t n, StreamT st);

    template<MATH_OPS op>
    static void binaryC(T* d_o, const T* d_i, const T& c, 
			size_t n, StreamT st,bool onDev);
    template<MATH_OPS op>
    static void binaryC(T* d_o, const T* d_i, const T& c, const T* d_mask, 
			size_t n, StreamT st,bool onDev);

    template<MATH_OPS op>
    static void binaryC_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    template<MATH_OPS op>
    static void binaryC_I(T* d_o, const T& c, const T* d_mask, 
			  size_t n, StreamT st,bool onDev);

    template<MATH_OPS op>
    static void binary(T* d_o, const T* d_i, const T* d_i1, 
		       size_t n, StreamT st);
    template<MATH_OPS op>
    static void binary(T* d_o, const T* d_i, const T* d_i1, const T* d_mask, 
		       size_t n, StreamT st);

    template<MATH_OPS op>
    static void binary_I(T* d_o, const T* d_i, size_t n, StreamT st);
    template<MATH_OPS op>
    static void binary_I(T* d_o, const T* d_i, const T* d_mask, 
			 size_t n, StreamT st);

public:
    static void Copy(T* d_o, const T* d_i, size_t n, StreamT st);
    static void Copy(T* d_o, const T* d_i, const T* d_mask, 
		     size_t n, StreamT st);
/**
 * Fill array with constant value
 */
    static void SetMem(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void SetMem(T* d_o, const T& c, const T* d_mask, 
		       size_t n, StreamT st,bool onDev);

/**
 * Set unsined array with the linear ramp up
 * d[i] = i
 */
    static void SetLinear(T* d_o,  size_t n, StreamT st);

/**
 * Set unsined array with the linear ramp down
 * d[i] = n - 1 - i;
 */
    static void SetLinearDown(T* d_o,  size_t n, StreamT st);

/////////////////////////////////////////////////////////////////////////////////
// Unary functions
//////////////////////////////////////////////////////////////////////////////////
/**
 * Absolute value of array (d_o = abs(d_i))
 */
    static void Abs(T* d_o, const T* d_i, size_t n, StreamT st);
    static void Abs(T* d_o, const T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * Absolute value of array, in-place (d_o = abs(d_o))
 */
    static void Abs_I(T* d_o, size_t n, StreamT st);
    static void Abs_I(T* d_o, const T* d_mask, size_t n, StreamT st);

/**
 * Cube each value in array (d_o = d_i^3)
 */
    static void Cube(T* d_o, const T* d_i, size_t n, StreamT st);
    static void Cube(T* d_o, const T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * Cube each value in array, in-place (d_o = d_o^3)
 */
    static void Cube_I(T* d_o, size_t n, StreamT st);
    static void Cube_I(T* d_o, const T* d_mask, size_t n, StreamT st);

/**
 * Exp of each value in array (d_o = exp(d_i))
 */
    static void Exp(T* d_o, const  T* d_i, size_t n, StreamT st);
    static void Exp(T* d_o, const  T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * Exp of each value in array (d_o = exp(d_o)) inplace version
 */
    static void Exp_I(T* d_o, size_t n, StreamT st);
    static void Exp_I(T* d_o, const T* d_mask, size_t n, StreamT st);

/**
 * Log of each value in array (d_o = log(d_i))
 */
    static void Log(T* d_o, const  T* d_i, size_t n, StreamT st);
    static void Log(T* d_o, const  T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * Log of each value in array (d_o = log(d_o)) inplace version
 */
    static void Log_I(T* d_o, size_t n, StreamT st);
    static void Log_I(T* d_o, const T* d_mask, size_t n, StreamT st);

/**
 * Negate each value in array (d_o = -d_i)
 */
    static void Neg(T* d_o, const T* d_i, size_t n, StreamT st);
    static void Neg(T* d_o, const T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * Negate each value in array, in-place (d_o = -d_o)
 */
    static void Neg_I(T* d_o, size_t n, StreamT st);
    static void Neg_I(T* d_o, const T* d_mask, size_t n, StreamT st);

/**
 * Ramp function of each value in array (d_o = ramp(d_i))
 */
    static void Ramp(T* d_o, const T* d_i, size_t n, StreamT st);
    static void Ramp(T* d_o, const T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * Ramp function of each value in array in-place (d_o = sgn(d_o))
 */
    static void Ramp_I(T* d_o, size_t n, StreamT st);
    static void Ramp_I(T* d_o, const T* d_mask, size_t n, StreamT st);

/**
 * Sign of each value in array (d_o = sgn(d_i))
 */
    static void Sgn(T* d_o, const T* d_i, size_t n, StreamT st);
    static void Sgn(T* d_o, const T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * Sign of each value in array in-place (d_o = sgn(d_o))
 */
    static void Sgn_I(T* d_o, size_t n, StreamT st);
    static void Sgn_I(T* d_o, const T* d_mask, size_t n, StreamT st);

/**
 * Heaviside Step of each value in array (d_o = sgn(d_i))
 */
    static void Step(T* d_o, const T* d_i, size_t n, StreamT st);
    static void Step(T* d_o, const T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * Heaviside Step of each value in array in-place (d_o = sgn(d_o))
 */
    static void Step_I(T* d_o, size_t n, StreamT st);
    static void Step_I(T* d_o, const T* d_mask, size_t n, StreamT st);

/**
 * Square of each value in array (d_o = sqr(d_i))
 */
    static void Sqr(T* d_o, const T* d_i, size_t n, StreamT st);
    static void Sqr(T* d_o, const T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * Square of each value in array, in-place (d_o = sqr(d_o))
 */
    static void Sqr_I(T* d_o, size_t n, StreamT st);
    static void Sqr_I(T* d_o, const T* d_mask, size_t n, StreamT st);

/**
 * Square root of each value in array (d_o = sqrt(d_i))
 */
    static void Sqrt(T* d_o, const  T* d_i, size_t n, StreamT st);
    static void Sqrt(T* d_o, const  T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * Square root of each value in array (d_o = sqrt(d_o)) inplace version
 */
    static void Sqrt_I(T* d_o, size_t n, StreamT st);
    static void Sqrt_I(T* d_o, const T* d_mask, size_t n, StreamT st);

/**
 * Inverse of each value in array (d_o = 1.0/d_i)
 */
    static void Inv(T* d_o, const T* d_i, size_t n, StreamT st);
    static void Inv(T* d_o, const T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * Inverse of each value in array (d_o = 1.0/d_o)
 */
    static void Inv_I(T* d_o, size_t n, StreamT st);
    static void Inv_I(T* d_o, const T* d_mask, size_t n, StreamT st);

/**
 * ceiling of each value in array
 */
    static void Ceil(T* d_o, const T* d_i, size_t n, StreamT st);
    static void Ceil(T* d_o, const T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * ceiling of each value in array, inplace version
 */
    static void Ceil_I(T* d_o, size_t n, StreamT st);
    static void Ceil_I(T* d_o, const T* d_mask, size_t n, StreamT st);

/**
 * floor of each value in array
 */
    static void Floor(T* d_o, const T* d_i, size_t n, StreamT st);
    static void Floor(T* d_o, const T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * floor of each value in array, inplace version
 */
    static void Floor_I(T* d_o, size_t n, StreamT st);
    static void Floor_I(T* d_o, const T* d_mask, size_t n, StreamT st);

/**
 * round each value in array
 */
    static void Round(T* d_o, const T* d_i, size_t n, StreamT st);
    static void Round(T* d_o, const T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * round each value in array, inplace version
 */
    static void Round_I(T* d_o, size_t n, StreamT st);
    static void Round_I(T* d_o, const T* d_mask, size_t n, StreamT st);

/**
 * sine in radians of each value in array (d_o = sin(d_i))
 */
    static void Sin(T* d_o, const  T* d_i, size_t n, StreamT st);
    static void Sin(T* d_o, const  T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * sine in radians of each value in array (d_o = sin(d_o)) inplace version
 */
    static void Sin_I(T* d_o, size_t n, StreamT st);
    static void Sin_I(T* d_o, const T* d_mask, size_t n, StreamT st);

/**
 * arcsine in radians of each value in array (d_o = asin(d_i))
 */
    static void Asin(T* d_o, const  T* d_i, size_t n, StreamT st);
    static void Asin(T* d_o, const  T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * arcsine in radians of each value in array (d_o = asin(d_o)) inplace version
 */
    static void Asin_I(T* d_o, size_t n, StreamT st);
    static void Asin_I(T* d_o, const T* d_mask, size_t n, StreamT st);

/**
 * cosine in radians of each value in array (d_o = cos(d_i))
 */
    static void Cos(T* d_o, const  T* d_i, size_t n, StreamT st);
    static void Cos(T* d_o, const  T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * cosine in radians of each value in array (d_o = cos(d_o)) inplace version
 */
    static void Cos_I(T* d_o, size_t n, StreamT st);
    static void Cos_I(T* d_o, const T* d_mask, size_t n, StreamT st);

/**
 * arccosine in radians of each value in array (d_o = acos(d_i))
 */
    static void Acos(T* d_o, const  T* d_i, size_t n, StreamT st);
    static void Acos(T* d_o, const  T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * arccosine in radians of each value in array (d_o = acos(d_o)) inplace version
 */
    static void Acos_I(T* d_o, size_t n, StreamT st);
    static void Acos_I(T* d_o, const T* d_mask, size_t n, StreamT st);

/**
 * tangent in radians of each value in array (d_o = tan(d_i))
 */
    static void Tan(T* d_o, const  T* d_i, size_t n, StreamT st);
    static void Tan(T* d_o, const  T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * tangent in radians of each value in array (d_o = tan(d_o)) inplace version
 */
    static void Tan_I(T* d_o, size_t n, StreamT st);
    static void Tan_I(T* d_o, const T* d_mask, size_t n, StreamT st);

/**
 * arctangent in radians of each value in array (d_o = atan(d_i))
 */
    static void Atan(T* d_o, const  T* d_i, size_t n, StreamT st);
    static void Atan(T* d_o, const  T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * arctangent in radians of each value in array (d_o = atan(d_o)) inplace version
 */
    static void Atan_I(T* d_o, size_t n, StreamT st);
    static void Atan_I(T* d_o, const T* d_mask, size_t n, StreamT st);

/**
 * cosecant in radians of each value in array (d_o = csc(d_i))
 */
    static void Csc(T* d_o, const  T* d_i, size_t n, StreamT st);
    static void Csc(T* d_o, const  T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * cosecant in radians of each value in array (d_o = csc(d_o)) inplace version
 */
    static void Csc_I(T* d_o, size_t n, StreamT st);
    static void Csc_I(T* d_o, const T* d_mask, size_t n, StreamT st);

/**
 * secant in radians of each value in array (d_o = sec(d_i))
 */
    static void Sec(T* d_o, const  T* d_i, size_t n, StreamT st);
    static void Sec(T* d_o, const  T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * secant in radians of each value in array (d_o = sec(d_o)) inplace version
 */
    static void Sec_I(T* d_o, size_t n, StreamT st);
    static void Sec_I(T* d_o, const T* d_mask, size_t n, StreamT st);

/**
 * cotangent in radians of each value in array (d_o = cot(d_i))
 */
    static void Cot(T* d_o, const  T* d_i, size_t n, StreamT st);
    static void Cot(T* d_o, const  T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * cotangent in radians of each value in array (d_o = cot(d_o)) inplace version
 */
    static void Cot_I(T* d_o, size_t n, StreamT st);
    static void Cot_I(T* d_o, const T* d_mask, size_t n, StreamT st);




/////////////////////////////////////////////////////////////////////////////////
// Binary functions
//////////////////////////////////////////////////////////////////////////////////
/**
 * Add constant to an array (d_o = d_i + c)
 */
     static void AddC(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void AddC(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
     static void Add(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void Add(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Add constant to an array, in-place (d_o = d_o + c)
 */
     static void AddC_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void AddC_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
     static void Add_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void Add_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Add two arrays, (d_o = d_i + d_i1)
 */
     static void Add(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT st);
    static void Add(T* d_o, const T* d_i, const T* d_i1, const T* d_mask, size_t n, StreamT st);

/**
 * Add two arrays, in-place (d_o = d_o + d_i)
 */
    static void Add_I(T* d_o, const T* d_i, size_t n, StreamT st);
    static void Add_I(T* d_o, const T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * Subtract constant from an array (d_o = d_i - c)
 */
    static void SubC(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void SubC(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void Sub(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void Sub(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Subtract constant from an array, in-place (d_o = d_o - c)
 */
    static void SubC_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void SubC_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void Sub_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void Sub_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Subtract two arrays, (d_o = d_i - d_i1)
 */
    static void Sub(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT st);
    static void Sub(T* d_o, const T* d_i, const T* d_i1, const T* d_mask, size_t n, StreamT st);

/**
 * Subtract two arrays, in-place (d_o = d_o - d_i)
 */
    static void Sub_I(T* d_o, const T* d_i, size_t n, StreamT st);
    static void Sub_I(T* d_o, const T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * Multiply an array by a constant (scale array) (d_o = d_i * c)
 */
    static void MulC(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void MulC(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void Mul(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void Mul(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Multiply an array by a constant (scale array) in-place (d_o = d_o * c)
 */
    static void MulC_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void MulC_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void Mul_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void Mul_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Multiply two arrays (d_o = d_i * d_i1)
 */
    static void Mul(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT st);
    static void Mul(T* d_o, const T* d_i, const T* d_i1, const T* d_mask, size_t n, StreamT st);

/**
 * Multiply two arrays, in-place (d_o = d_o * d_i)
 */
    static void Mul_I(T* d_o, const T* d_i, size_t n, StreamT st);
    static void Mul_I(T* d_o, const T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * Divide an array by a constant (scale array) (d_o = d_i / c)
 */
     static void DivC(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void DivC(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void Div(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void Div(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Divide an array by a constant (scale array) in-place (d_o = d_o / c)
 */
    static void DivC_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void DivC_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void Div_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void Div_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Divide two arrays (d_o = d_i / d_i1)
 */
    static void Div(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT st);
    static void Div(T* d_o, const T* d_i, const T* d_i1, const T* d_mask, size_t n, StreamT st);

/**
 * Divide two arrays, in-place (d_o = d_o / d_i)
 */
    static void Div_I(T* d_o, const T* d_i, size_t n, StreamT st);
    static void Div_I(T* d_o, const T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * max of array and constant  (d_o = max(d_i, c))
 */
     static void MaxC(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void MaxC(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
     static void Max(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void Max(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * max of array and constant  (d_o = max(d_o, c))
 */
    static void MaxC_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void MaxC_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void Max_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void Max_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * component-wise max of two arrays (d_o = max(d_i, d_i1))
 */
    static void Max(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT st);
    static void Max(T* d_o, const T* d_i, const T* d_i1, const T* d_mask, size_t n, StreamT st);

/**
 * component-wise max of two arrays, in-place version (d_o = max(d_o, d_i))
 */
    static void Max_I(T* d_o, const T* d_i, size_t n, StreamT st);
    static void Max_I(T* d_o, const T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * min of array and constant  (d_o = min(d_i, c))
 */
     static void MinC(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void MinC(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
     static void Min(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void Min(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * min of array and constant  (d_o = min(d_o, c))
 */
    static void MinC_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void MinC_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void Min_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void Min_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * component-wise min of two arrays (d_o = min(d_i, d_i1))
 */
    static void Min(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT st);
    static void Min(T* d_o, const T* d_i, const T* d_i1, const T* d_mask, size_t n, StreamT st);

/**
 * component-wise min of two arrays, in-place version (d_o = min(d_o, d_i))
 */
    static void Min_I(T* d_o, const T* d_i, size_t n, StreamT st);
    static void Min_I(T* d_o, const T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * pow (exponentiation) do(x) = di(x)^c
 */
     static void PowC(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void PowC(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
     static void Pow(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void Pow(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * pow (exponentiation), in-place (do(x) = do(x)^c)
 */
    static void PowC_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void PowC_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void Pow_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void Pow_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * arctangent given by x coordinate and y coordinate (d_o = atan2(d_i1, d_i2)
 */
    static void Atan2(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT st);
    static void Atan2(T* d_o, const T* d_i, const T* d_i1, const T* d_mask, size_t n, StreamT st);

/**
 * arctangent given by x coordinate and y coordinate, inplace (d_o = atan2(d_o, d_i1)
 */
    static void Atan2_I(T* d_o, const T* d_i, size_t n, StreamT st);
    static void Atan2_I(T* d_o, const T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * arctangent given by x coordinate and y coordinate (d_o = atan2(d_i1, c)
 */
     static void Atan2C(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void Atan2C(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
     static void Atan2(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void Atan2(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * arctangent given by x coordinate and y coordinate, inplace (d_o = atan2(d_o, c)
 */
    static void Atan2C_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void Atan2C_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void Atan2_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void Atan2_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Greater than of an array and constant (d_o = d_i > c)
 */
    static void GTC(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void GTC(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void GT(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void GT(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Greater than of an array and constant, in-place (d_o = d_o > c)
 */
    static void GTC_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void GTC_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void GT_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void GT_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Greater than of two arrays, (d_o = d_i > d_i1)
 */
    static void GT(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT st);
    static void GT(T* d_o, const T* d_i, const T* d_i1, const T* d_mask, size_t n, StreamT st);

/**
 * Greater than of two arrays, in-place (d_o = d_o > d_i)
 */
    static void GT_I(T* d_o, const T* d_i, size_t n, StreamT st);
    static void GT_I(T* d_o, const T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * Greater than or equal of an array and constant (d_o = d_i >= c)
 */
    static void GTEC(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void GTEC(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void GTE(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void GTE(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Greater than or equal of an array and constant, in-place (d_o = d_o >= c)
 */
    static void GTEC_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void GTEC_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void GTE_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void GTE_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Greater than or equal of two arrays, (d_o = d_i >= d_i1)
 */
    static void GTE(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT st);
    static void GTE(T* d_o, const T* d_i, const T* d_i1, const T* d_mask, size_t n, StreamT st);

/**
 * Greater than or equal of two arrays, in-place (d_o = d_o >= d_i)
 */
    static void GTE_I(T* d_o, const T* d_i, size_t n, StreamT st);
    static void GTE_I(T* d_o, const T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * Boolean Equal of an array and constant (d_o = d_i == c)
 */
    static void EQC(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void EQC(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void EQ(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void EQ(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Boolean Equal of an array and constant, in-place (d_o = d_o == c)
 */
    static void EQC_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void EQC_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void EQ_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void EQ_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Boolean Equal of two arrays, (d_o = d_i == d_i1)
 */
    static void EQ(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT st);
    static void EQ(T* d_o, const T* d_i, const T* d_i1, const T* d_mask, size_t n, StreamT st);

/**
 * Boolean Equal of two arrays, in-place (d_o = d_o == d_i)
 */
    static void EQ_I(T* d_o, const T* d_i, size_t n, StreamT st);
    static void EQ_I(T* d_o, const T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * Boolean Not Equal of an array and constant (d_o = d_i != c)
 */
    static void NEQC(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void NEQC(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void NEQ(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void NEQ(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Boolean Not Equal of an array and constant, in-place (d_o = d_o != c)
 */
    static void NEQC_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void NEQC_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void NEQ_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void NEQ_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Boolean Not Equal of two arrays, (d_o = d_i != d_i1)
 */
    static void NEQ(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT st);
    static void NEQ(T* d_o, const T* d_i, const T* d_i1, const T* d_mask, size_t n, StreamT st);

/**
 * Boolean Not Equal of two arrays, in-place (d_o = d_o != d_i)
 */
    static void NEQ_I(T* d_o, const T* d_i, size_t n, StreamT st);
    static void NEQ_I(T* d_o, const T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * Less than of an array and constant (d_o = d_i < c)
 */
    static void LTC(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void LTC(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void LT(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void LT(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Less than of an array and constant, in-place (d_o = d_o < c)
 */
    static void LTC_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void LTC_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void LT_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void LT_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Less than of two arrays, (d_o = d_i < d_i1)
 */
    static void LT(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT st);
    static void LT(T* d_o, const T* d_i, const T* d_i1, const T* d_mask, size_t n, StreamT st);

/**
 * Less than of two arrays, in-place (d_o = d_o < d_i)
 */
    static void LT_I(T* d_o, const T* d_i, size_t n, StreamT st);
    static void LT_I(T* d_o, const T* d_i, const T* d_mask, size_t n, StreamT st);

/**
 * Less than or equal of an array and constant (d_o = d_i <= c)
 */
    static void LTEC(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void LTEC(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void LTE(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);
    static void LTE(T* d_o, const T* d_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Less than or equal of an array and constant, in-place (d_o = d_o <= c)
 */
    static void LTEC_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void LTEC_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void LTE_I(T* d_o, const T& c, size_t n, StreamT st,bool onDev);
    static void LTE_I(T* d_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Less than or equal of two arrays, (d_o = d_i <= d_i1)
 */
    static void LTE(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT st);
    static void LTE(T* d_o, const T* d_i, const T* d_i1, const T* d_mask, size_t n, StreamT st);

/**
 * Less than or equal of two arrays, in-place (d_o = d_o <= d_i)
 */
    static void LTE_I(T* d_o, const T* d_i, size_t n, StreamT st);
    static void LTE_I(T* d_o, const T* d_i, const T* d_mask, size_t n, StreamT st);


/////////////////////////////////////////////////////////////////////////////////
// Triary functions
//////////////////////////////////////////////////////////////////////////////////

/**
 * Compute the absolution different between 2 vector
 * d_o = abs(d_i - d_i1)
 */
    static void AbsDiff(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT st);

/**
 * Compute the absolution different between 2 vector in-place version
 * d_o = abs(d_o - d_i1)
 */
    static void AbsDiff_I(T* d_o, const T* d_i, size_t n, StreamT st);

/**
 * Compute the absolution different between 2 vector
 * d_o = sqr(d_i - d_i1)
 */

    static void SqrDiff(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT st);

/**
 * Compute the absolution different between 2 vector in-place version
 * d_o = sqr(d_o - d_i1)
 */

    static void SqrDiff_I(T* d_o, const T* d_i, size_t n, StreamT st);

/**
 * d_o = d_i * (d_i1 * c)
 */

    static void MulMulC(T* d_o, const T* d_i, const T* d_i1, const T& c, size_t n, StreamT st,bool onDev);

/**
 * d_o = d_o * (d_i1 * c)
 */

    static void MulMulC_I(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);

/**
 * d_o = d_i * (d_i1 * d_i2)
 */

    static void MulMul(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, size_t n, StreamT st);

/**
 * d_o = d_o * (d_i * d_i1)
 */

    static void MulMul_I(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT st);

/**
 * d_o = d_i * d_i1 + d_i2
 */

    static void MulAdd(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, size_t n, StreamT st);

/**
 * d_o = d_o * d_i + d_i1
 */

    static void MulAdd_I(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT st);

/**
 * d_o = d_i * d_i1 - d_i2
 */

    static void MulSub(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, size_t n, StreamT st);

/**
 * d_o = d_o * d_i - d_i1
 */

    static void MulSub_I(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT st);


/**
 * d_o = (d_i + d_i1) * d_i2
 */

    static void AddMul(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, size_t n, StreamT st);

/**
 * d_o = (d_o + d_i) * d_i1
 */

    static void AddMul_I(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT st);

/**
 * d_o = (d_i + d_i1) / d_i2
 */

    static void AddDiv(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, size_t n, StreamT st);

/**
 * d_o = (d_o + d_i) / d_i1
 */

    static void AddDiv_I(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT st);

/**
 * d_o = (d_i - d_i1) * d_i2
 */

    static void SubMul(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, size_t n, StreamT st);

/**
 * d_o = (d_o - d_i) * d_i1
 */

    static void SubMul_I(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT st);

/**
 * d_o = (d_i - d_i1) / d_i2
 */

    static void SubDiv(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, size_t n, StreamT st);

/**
 * d_o = (d_o - d_i) / d_i1
 */

    static void SubDiv_I(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT st);

/**
 * d_o = (d_i * c) + d_i1
 */

    static void MulCAdd(T* d_o, const T* d_i, const T& c, const T* d_i1, size_t n, StreamT st,bool onDev);

/**
 * d_o = (d_o * c) + d_i
 */

    static void MulCAdd_I(T* d_o, const T& c, const T* d_i, size_t n, StreamT st,bool onDev);

/**
 * d_o = (d_i * c) - d_i1
 */

    static void MulCSub(T* d_o, const T* d_i, const T& c, const T* d_i1, size_t n, StreamT st,bool onDev);

/**
 * d_o = (d_o * c) - d_i
 */

    static void MulCSub_I(T* d_o, const T& c, const T* d_i, size_t n, StreamT st,bool onDev);

/**
 * d_o = (d_i * a) + b
 */

    static void MulCAddC(T* d_o, const T* d_i, const T& a, const T& b, size_t n, StreamT st,bool onDev);

/**
 * d_o = (d_o * a) + b
 */
    static void MulCAddC_I(T* d_o, const T& a, const T& b, size_t n, StreamT st,bool onDev);

/**
 * d_o = (d_i + d_i1) * c
 */
    static void AddMulC(T* d_o, const T* d_i, const T* d_i1, const T& c, size_t n, StreamT st,bool onDev);

/**
 * d_o = (d_o + d_i) * c
 */
    static void AddMulC_I(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);

/**
 * d_o = (d_i - d_i1) * c
 */
    static void SubMulC(T* d_o, const T* d_i, const T* d_i1, const T& c, size_t n, StreamT st,bool onDev);

/**
 * d_o = (d_o - d_i) * c
 */
    static void SubMulC_I(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);


/**
 * d_o = (d_i + a) * b
 */
    static void AddCMulC(T* d_o, const T* d_i, const T& a, const T& b, size_t n, StreamT st,bool onDev);

/**
 * d_o = (d_o + a) * b
 */
    static void AddCMulC_I(T* d_o, const T& a, const T& b, size_t n, StreamT st,bool onDev);

/**
 * d_o = d_i + d_i1 * c
 */
    static void Add_MulC(T* d_o, const T* d_i, const T* d_i1, const T& c, size_t n, StreamT st,bool onDev);

/**
 * d_o = d_o + d_i * c
 */
    static void Add_MulC_I(T* d_o, const T* d_i, const T& c, size_t n, StreamT st,bool onDev);

/**
 * d_o = d_i + d_i1 * d_i2
 */
    static void Add_Mul(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, size_t n, StreamT st);
/**
 * d_o = d_o + d_i * d_i1
 */
    static void Add_Mul_I(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT st);

 /**
 * d_o = d_i - d_i1 * d_i2
 */
    static void Sub_Mul(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, size_t n, StreamT st);
/**
 * d_o = d_o - d_i * d_i1
 */
    static void Sub_Mul_I(T* d_o, const T* d_i, const T* d_i1, size_t n, StreamT st);

/////////////////////////////////////////////////////////////////////////////////
// n-ary functions
//////////////////////////////////////////////////////////////////////////////////
/**
 * d_o = d_i + (d_i1 + d_i2) * c
 */
    static void Add_AddMulC(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, const T& c, size_t n, StreamT st,bool onDev);

/**
 * d_o = d_i + (d_i1 - d_i2) * c
 */
    static void Add_SubMulC(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, const T& c, size_t n, StreamT st,bool onDev);

/**
 * d_o = d_i + (d_i1 * d_i2) * c
 */
    static void Add_MulMulC(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, const T& c, size_t n, StreamT st,bool onDev);


/**
 * d_o += (d_i + d_i1) * c
 */
    static void Add_AddMulC_I(T* d_o, const T* d_i, const T* d_i1, const T& c, size_t n, StreamT st,bool onDev);

/**
 * d_o += (d_i - d_i1) * c
 */
    static void Add_SubMulC_I(T* d_o, const T* d_i, const T* d_i1, const T& c, size_t n, StreamT st,bool onDev);

/**
 * d_o += (d_i * d_i1) * c
 */
    static void Add_MulMulC_I(T* d_o, const T* d_i, const T* d_i1, const T& c, size_t n, StreamT st,bool onDev);

    /**
 * d_o = d_i - (d_i1 + d_i2) * c
 */
    static void Sub_AddMulC(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, const T& c, size_t n, StreamT st,bool onDev);

/**
 * d_o = d_i - (d_i1 - d_i2) * c
 */
    static void Sub_SubMulC(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, const T& c, size_t n, StreamT st,bool onDev);

/**
 * d_o = d_i - (d_i1 * d_i2) * c
 */
    static void Sub_MulMulC(T* d_o, const T* d_i, const T* d_i1, const T* d_i2, const T& c, size_t n, StreamT st,bool onDev);


/**
 * d_o -= (d_i + d_i1) * c
 */
    static void Sub_AddMulC_I(T* d_o, const T* d_i, const T* d_i1, const T& c, size_t n, StreamT st,bool onDev);

/**
 * d_o -= (d_i - d_i1) * c
 */
    static void Sub_SubMulC_I(T* d_o, const T* d_i, const T* d_i1, const T& c, size_t n, StreamT st,bool onDev);

/**
 * d_o -= (d_i * d_i1) * c
 */
    static void Sub_MulMulC_I(T* d_o, const T* d_i, const T* d_i1, const T& c, size_t n, StreamT st,bool onDev);


// Functions with 4 inputs

/**
 * d_o = d_i * a + d_i1 * b
 */
    static void MulC_Add_MulC(T* d_o, const T* d_i, const T& a, const T* d_i1, const T& b, size_t n, StreamT st,bool onDev);

/**
 * d_o = d_o * a + d_i1 * b
 */
    static void MulC_Add_MulC_I(T* d_o, const T& a, const T* d_i, const T& b, size_t n, StreamT st,bool onDev);

/**
 * d_o = (d_i + a) * b + c
 */
    static void AddCMulCAddC(T* d_o, const T* d_i, const T& a, const T& b, const T& c, size_t n, StreamT st,bool onDev);

/**
 * d_o = (d_o + a) * b + c
 */
    static void AddCMulCAddC_I(T* d_o, const T& a, const T& b, const T& c, size_t n, StreamT st,bool onDev);

/**
 * d_o[i] = d_i[n-1 - i]
 */
    static void ReverseOrder(T* d_o, const T* d_i, size_t n, StreamT st);

    static void ShiftCoordinate(T* d_o, const T* d_i, size_t sizeX, size_t sizeY, size_t sizeZ, bool dir, StreamT s);
};

} // end namespace PyCA

#endif

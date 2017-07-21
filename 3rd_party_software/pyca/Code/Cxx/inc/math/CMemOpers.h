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

#ifndef  __CPU_MEM_OPERATION__H
#define  __CPU_MEM_OPERATION__H

#include <pycaConst.h>
#include <MOpers.h>
#include <estream.h>

namespace PyCA {

template<typename T>
class CMemOpers{
private:
    template<MATH_OPS op>
    static void Comp_unary(T* h_o, const T* h_i, size_t n, StreamT st);
    template<MATH_OPS op>
    static void Comp_unary(T* h_o, const T* h_i, const T* h_mask, 
			   size_t n, StreamT st);

    template<MATH_OPS op>
    static void Comp_unary_I(T* h_o, size_t n, StreamT st);
    template<MATH_OPS op>
    static void Comp_unary_I(T* h_o, const T* h_mask, size_t n, StreamT st);

    template<MATH_OPS op>
    static void binaryC(T* h_o, const T* h_i, const T& c, 
			size_t n, StreamT st,bool onDev);
    template<MATH_OPS op>
    static void binaryC(T* h_o, const T* h_i, const T& c, const T* h_mask, 
			size_t n, StreamT st,bool onDev);

    template<MATH_OPS op>
    static void binaryC_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    template<MATH_OPS op>
    static void binaryC_I(T* h_o, const T& c, const T* h_mask, 
			  size_t n, StreamT st,bool onDev);

    template<MATH_OPS op>
    static void binary(T* h_o, const T* h_i, const T* h_i1, 
		       size_t n, StreamT st);
    template<MATH_OPS op>
    static void binary(T* h_o, const T* h_i, const T* h_i1, const T* h_mask, 
		       size_t n, StreamT st);

    template<MATH_OPS op>
    static void binary_I(T* h_o, const T* h_i, size_t n, StreamT st);
    template<MATH_OPS op>
    static void binary_I(T* h_o, const T* h_i, const T* h_mask, 
			 size_t n, StreamT st);
public:
    static void Copy(T* h_o, const T* h_i, size_t n, StreamT st);
    static void Copy(T* h_o, const T* h_i, const T* d_mask,
		     size_t n, StreamT st);
/**
 * Fill array with constant value
 */
    static void SetMem(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void SetMem(T* h_o, const T& c, const T* h_mask, 
		       size_t n, StreamT st,bool onDev);

/**
 * Set unsined array with the linear ramp up
 * d[i] = i
 */
    static void SetLinear(T* h_o,  size_t n, StreamT st);

/**
 * Set unsined array with the linear ramp down
 * d[i] = n - 1 - i;
 */
    static void SetLinearDown(T* h_o,  size_t n, StreamT st);

/////////////////////////////////////////////////////////////////////////////////
// Unary functions
//////////////////////////////////////////////////////////////////////////////////
/**
 * Absolute value of array (h_o = abs(h_i))
 */
    static void Abs(T* h_o, const T* h_i, size_t n, StreamT st);
    static void Abs(T* h_o, const T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * Absolute value of array, in-place (h_o = abs(h_o))
 */
    static void Abs_I(T* h_o, size_t n, StreamT st);
    static void Abs_I(T* h_o, const T* d_mask, size_t n, StreamT st);

/**
 * Cube each value in array (h_o = h_i^3)
 */
    static void Cube(T* h_o, const T* h_i, size_t n, StreamT st);
    static void Cube(T* h_o, const T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * Cube each value in array, in-place (h_o = h_o^3)
 */
    static void Cube_I(T* h_o, size_t n, StreamT st);
    static void Cube_I(T* h_o, const T* d_mask, size_t n, StreamT st);

/**
 * Exponential of each value in array (h_o = exp(h_i))
 */
    static void Exp(T* h_o, const  T* h_i, size_t n, StreamT st);
    static void Exp(T* h_o, const  T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * Exponential of each value in array (h_o = exp(h_o)) inplace version
 */
    static void Exp_I(T* h_o, size_t n, StreamT st);
    static void Exp_I(T* h_o, const T* d_mask, size_t n, StreamT st);

/**
 * Logarithm of each value in array (h_o = log(h_i))
 */
    static void Log(T* h_o, const  T* h_i, size_t n, StreamT st);
    static void Log(T* h_o, const  T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * Logarithm of each value in array (h_o = log(h_o)) inplace version
 */
    static void Log_I(T* h_o, size_t n, StreamT st);
    static void Log_I(T* h_o, const T* d_mask, size_t n, StreamT st);

/**
 * Negate each value in array (h_o = -h_i)
 */
    static void Neg(T* h_o, const T* h_i, size_t n, StreamT st);
    static void Neg(T* h_o, const T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * Negate each value in array, in-place (h_o = -h_o)
 */
    static void Neg_I(T* h_o, size_t n, StreamT st);
    static void Neg_I(T* h_o, const T* d_mask, size_t n, StreamT st);

/**
 * Ramp function of each value in array (h_o = ramp(h_i))
 */
    static void Ramp(T* h_o, const T* h_i, size_t n, StreamT st);
    static void Ramp(T* h_o, const T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * Sign of each value in array in-place (h_o = ramp(h_o))
 */
    static void Ramp_I(T* h_o, size_t n, StreamT st);
    static void Ramp_I(T* h_o, const T* d_mask, size_t n, StreamT st);

/**
 * Sign of each value in array (h_o = sgn(h_i))
 */
    static void Sgn(T* h_o, const T* h_i, size_t n, StreamT st);
    static void Sgn(T* h_o, const T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * Sign of each value in array in-place (h_o = sgn(h_o))
 */
    static void Sgn_I(T* h_o, size_t n, StreamT st);
    static void Sgn_I(T* h_o, const T* d_mask, size_t n, StreamT st);

/**
 * Heaviside Step of each value in array (h_o = step(h_i))
 */
    static void Step(T* h_o, const T* h_i, size_t n, StreamT st);
    static void Step(T* h_o, const T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * Heaviside Step of each value in array in-place (h_o = sgn(h_o))
 */
    static void Step_I(T* h_o, size_t n, StreamT st);
    static void Step_I(T* h_o, const T* d_mask, size_t n, StreamT st);

/**
 * Square of each value in array (h_o = sqr(h_i))
 */
    static void Sqr(T* h_o, const T* h_i, size_t n, StreamT st);
    static void Sqr(T* h_o, const T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * Square of each value in array, in-place (h_o = sqr(h_o))
 */
    static void Sqr_I(T* h_o, size_t n, StreamT st);
    static void Sqr_I(T* h_o, const T* d_mask, size_t n, StreamT st);

/**
 * Square root of each value in array (h_o = sqrt(h_i))
 */
    static void Sqrt(T* h_o, const  T* h_i, size_t n, StreamT st);
    static void Sqrt(T* h_o, const  T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * Square root of each value in array (h_o = sqrt(h_o)) inplace version
 */
    static void Sqrt_I(T* h_o, size_t n, StreamT st);
    static void Sqrt_I(T* h_o, const T* d_mask, size_t n, StreamT st);

/**
 * Inverse of each value in array (h_o = 1.0/h_i)
 */
    static void Inv(T* h_o, const T* h_i, size_t n, StreamT st);
    static void Inv(T* h_o, const T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * Inverse of each value in array (h_o = 1.0/h_o)
 */
    static void Inv_I(T* h_o, size_t n, StreamT st);
    static void Inv_I(T* h_o, const T* d_mask, size_t n, StreamT st);

/**
 * ceiling of each value in array
 */
    static void Ceil(T* h_o, const T* h_i, size_t n, StreamT st);
    static void Ceil(T* h_o, const T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * ceiling of each value in array, inplace version
 */
    static void Ceil_I(T* h_o, size_t n, StreamT st);
    static void Ceil_I(T* h_o, const T* d_mask, size_t n, StreamT st);

/**
 * floor of each value in array
 */
    static void Floor(T* h_o, const T* h_i, size_t n, StreamT st);
    static void Floor(T* h_o, const T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * floor of each value in array, inplace version
 */
    static void Floor_I(T* h_o, size_t n, StreamT st);
    static void Floor_I(T* h_o, const T* d_mask, size_t n, StreamT st);

/**
 * round each value in array
 */
    static void Round(T* h_o, const T* h_i, size_t n, StreamT st);
    static void Round(T* h_o, const T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * round each value in array, inplace version
 */
    static void Round_I(T* h_o, size_t n, StreamT st);
    static void Round_I(T* h_o, const T* d_mask, size_t n, StreamT st);

/**
 * sine in radians of each value in array (h_o = sin(h_i))
 */
    static void Sin(T* h_o, const  T* h_i, size_t n, StreamT st);
    static void Sin(T* h_o, const  T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * sine in radians of each value in array (h_o = sin(h_o)) inplace version
 */
    static void Sin_I(T* h_o, size_t n, StreamT st);
    static void Sin_I(T* h_o, const T* d_mask, size_t n, StreamT st);

/**
 * arcsine in radians of each value in array (h_o = asin(h_i))
 */
    static void Asin(T* h_o, const  T* h_i, size_t n, StreamT st);
    static void Asin(T* h_o, const  T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * arcsine in radians of each value in array (h_o = asin(h_o)) inplace version
 */
    static void Asin_I(T* h_o, size_t n, StreamT st);
    static void Asin_I(T* h_o, const T* d_mask, size_t n, StreamT st);

/**
 * cosine in radians of each value in array (h_o = cos(h_i))
 */
    static void Cos(T* h_o, const  T* h_i, size_t n, StreamT st);
    static void Cos(T* h_o, const  T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * cosine in radians of each value in array (h_o = cos(h_o)) inplace version
 */
    static void Cos_I(T* h_o, size_t n, StreamT st);
    static void Cos_I(T* h_o, const T* d_mask, size_t n, StreamT st);

/**
 * arccosine in radians of each value in array (h_o = acos(h_i))
 */
    static void Acos(T* h_o, const  T* h_i, size_t n, StreamT st);
    static void Acos(T* h_o, const  T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * arccosine in radians of each value in array (h_o = acos(h_o)) inplace version
 */
    static void Acos_I(T* h_o, size_t n, StreamT st);
    static void Acos_I(T* h_o, const T* d_mask, size_t n, StreamT st);

/**
 * tangent in radians of each value in array (h_o = tan(h_i))
 */
    static void Tan(T* h_o, const  T* h_i, size_t n, StreamT st);
    static void Tan(T* h_o, const  T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * tangent in radians of each value in array (h_o = tan(h_o)) inplace version
 */
    static void Tan_I(T* h_o, size_t n, StreamT st);
    static void Tan_I(T* h_o, const T* d_mask, size_t n, StreamT st);

/**
 * arctangent in radians of each value in array (h_o = atan(h_i))
 */
    static void Atan(T* h_o, const  T* h_i, size_t n, StreamT st);
    static void Atan(T* h_o, const  T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * arctangent in radians of each value in array (h_o = atan(h_o)) inplace version
 */
    static void Atan_I(T* h_o, size_t n, StreamT st);
    static void Atan_I(T* h_o, const T* d_mask, size_t n, StreamT st);

/**
 * cosecant in radians of each value in array (h_o = csc(h_i))
 */
    static void Csc(T* h_o, const  T* h_i, size_t n, StreamT st);
    static void Csc(T* h_o, const  T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * cosecant in radians of each value in array (h_o = csc(h_o)) inplace version
 */
    static void Csc_I(T* h_o, size_t n, StreamT st);
    static void Csc_I(T* h_o, const T* d_mask, size_t n, StreamT st);

/**
 * secant in radians of each value in array (h_o = sec(h_i))
 */
    static void Sec(T* h_o, const  T* h_i, size_t n, StreamT st);
    static void Sec(T* h_o, const  T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * secant in radians of each value in array (h_o = sec(h_o)) inplace version
 */
    static void Sec_I(T* h_o, size_t n, StreamT st);
    static void Sec_I(T* h_o, const T* d_mask, size_t n, StreamT st);

/**
 * cotangent in radians of each value in array (h_o = cot(h_i))
 */
    static void Cot(T* h_o, const  T* h_i, size_t n, StreamT st);
    static void Cot(T* h_o, const  T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * cotangent in radians of each value in array (h_o = cot(h_o)) inplace version
 */
    static void Cot_I(T* h_o, size_t n, StreamT st);
    static void Cot_I(T* h_o, const T* d_mask, size_t n, StreamT st);



/////////////////////////////////////////////////////////////////////////////////
// Binary functions
//////////////////////////////////////////////////////////////////////////////////
/**
 * Add constant to an array (h_o = h_i + c)
 */
     static void AddC(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void AddC(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
     static void Add(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void Add(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Add constant to an array, in-place (h_o = h_o + c)
 */
     static void AddC_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void AddC_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
     static void Add_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void Add_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Add two arrays, (h_o = h_i + h_i1)
 */
     static void Add(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT st);
    static void Add(T* h_o, const T* h_i, const T* h_i1, const T* d_mask, size_t n, StreamT st);

/**
 * Add two arrays, in-place (h_o = h_o + h_i)
 */
    static void Add_I(T* h_o, const T* h_i, size_t n, StreamT st);
    static void Add_I(T* h_o, const T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * Subtract constant from an array (h_o = h_i - c)
 */
    static void SubC(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void SubC(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void Sub(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void Sub(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Subtract constant from an array, in-place (h_o = h_o - c)
 */
    static void SubC_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void SubC_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void Sub_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void Sub_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Subtract two arrays, (h_o = h_i - h_i1)
 */
    static void Sub(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT st);
    static void Sub(T* h_o, const T* h_i, const T* h_i1, const T* d_mask, size_t n, StreamT st);

/**
 * Subtract two arrays, in-place (h_o = h_o - h_i)
 */
    static void Sub_I(T* h_o, const T* h_i, size_t n, StreamT st);
    static void Sub_I(T* h_o, const T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * Multiply an array by a constant (scale array) (h_o = h_i * c)
 */
    static void MulC(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void MulC(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void Mul(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void Mul(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Multiply an array by a constant (scale array) in-place (h_o = h_o * c)
 */
    static void MulC_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void MulC_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void Mul_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void Mul_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Multiply two arrays (h_o = h_i * h_i1)
 */
    static void Mul(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT st);
    static void Mul(T* h_o, const T* h_i, const T* h_i1, const T* d_mask, size_t n, StreamT st);

/**
 * Multiply two arrays, in-place (h_o = h_o * h_i)
 */
    static void Mul_I(T* h_o, const T* h_i, size_t n, StreamT st);
    static void Mul_I(T* h_o, const T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * Divide an array by a constant (scale array) (h_o = h_i / c)
 */
     static void DivC(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void DivC(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
     static void Div(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void Div(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Divide an array by a constant (scale array) in-place (h_o = h_o / c)
 */
    static void DivC_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void DivC_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void Div_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void Div_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Divide two arrays (h_o = h_i / h_i1)
 */
    static void Div(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT st);
    static void Div(T* h_o, const T* h_i, const T* h_i1, const T* d_mask, size_t n, StreamT st);

/**
 * Divide two arrays, in-place (h_o = h_o / h_i)
 */
    static void Div_I(T* h_o, const T* h_i, size_t n, StreamT st);
    static void Div_I(T* h_o, const T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * max of array and constant  (h_o = max(h_i, c))
 */
     static void MaxC(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void MaxC(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void Max(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void Max(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * max of array and constant  (h_o = max(h_o, c))
 */
    static void MaxC_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void MaxC_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void Max_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void Max_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * component-wise max of two arrays (h_o = max(h_i, h_i1))
 */
    static void Max(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT st);
    static void Max(T* h_o, const T* h_i, const T* h_i1, const T* d_mask, size_t n, StreamT st);

/**
 * component-wise max of two arrays, in-place version (h_o = max(h_o, h_i))
 */
    static void Max_I(T* h_o, const T* h_i, size_t n, StreamT st);
    static void Max_I(T* h_o, const T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * min of array and constant  (h_o = min(h_i, c))
 */
     static void MinC(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void MinC(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
     static void Min(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void Min(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * min of array and constant  (h_o = min(h_o, c))
 */
    static void MinC_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void MinC_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void Min_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void Min_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * component-wise min of two arrays (h_o = min(h_i, h_i1))
 */
    static void Min(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT st);
    static void Min(T* h_o, const T* h_i, const T* h_i1, const T* d_mask, size_t n, StreamT st);

/**
 * component-wise min of two arrays, in-place version (h_o = min(h_o, h_i))
 */
    static void Min_I(T* h_o, const T* h_i, size_t n, StreamT st);
    static void Min_I(T* h_o, const T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * pow (exponentiation) ho(x) = hi(x)^c
 */
     static void PowC(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void PowC(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
     static void Pow(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void Pow(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * pow (exponentiation), in-place (ho(x) = ho(x)^c)
 */
    static void PowC_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void PowC_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void Pow_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void Pow_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * arctangent given by y coordinate and x coordinate (h_o = atan2(h_i1, h_i2)
 */
    static void Atan2(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT st);
    static void Atan2(T* h_o, const T* h_i, const T* h_i1, const T* d_mask, size_t n, StreamT st);

/**
 * arctangent given by y coordinate and x coordinate, inplace (h_o = atan2(h_o, h_i1)
 */
    static void Atan2_I(T* h_o, const T* h_i, size_t n, StreamT st);
    static void Atan2_I(T* h_o, const T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * arctangent given by y coordinate and x coordinate (h_o = atan2(h_i1, c)
 */
     static void Atan2C(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void Atan2C(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
     static void Atan2(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void Atan2(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * arctangent given by y coordinate and x coordinate, inplace (h_o = atan2(h_o, c)
 */
    static void Atan2C_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void Atan2C_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void Atan2_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void Atan2_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Greater than of an array and constant (h_o = h_i > c)
 */
    static void GTC(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void GTC(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void GT(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void GT(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Greater than of an array and constant, in-place (h_o = h_o > c)
 */
    static void GTC_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void GTC_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void GT_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void GT_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Greater than of two arrays, (h_o = h_i > h_i1)
 */
    static void GT(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT st);
    static void GT(T* h_o, const T* h_i, const T* h_i1, const T* d_mask, size_t n, StreamT st);

/**
 * Greater than of two arrays, in-place (h_o = h_o > h_i)
 */
    static void GT_I(T* h_o, const T* h_i, size_t n, StreamT st);
    static void GT_I(T* h_o, const T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * Greater than or equal of an array and constant (h_o = h_i >= c)
 */
    static void GTEC(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void GTEC(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void GTE(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void GTE(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Greater than or equal of an array and constant, in-place (h_o = h_o >= c)
 */
    static void GTEC_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void GTEC_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void GTE_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void GTE_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Greater than or equal of two arrays, (h_o = h_i >= h_i1)
 */
    static void GTE(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT st);
    static void GTE(T* h_o, const T* h_i, const T* h_i1, const T* d_mask, size_t n, StreamT st);

/**
 * Greater than or equal of two arrays, in-place (h_o = h_o >= h_i)
 */
    static void GTE_I(T* h_o, const T* h_i, size_t n, StreamT st);
    static void GTE_I(T* h_o, const T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * Boolean Equal of an array and constant (h_o = h_i == c)
 */
    static void EQC(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void EQC(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void EQ(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void EQ(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Boolean Equal of an array and constant, in-place (h_o = h_o == c)
 */
    static void EQC_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void EQC_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void EQ_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void EQ_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Boolean Equal of two arrays, (h_o = h_i == h_i1)
 */
    static void EQ(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT st);
    static void EQ(T* h_o, const T* h_i, const T* h_i1, const T* d_mask, size_t n, StreamT st);

/**
 * Boolean Equal of two arrays, in-place (h_o = h_o == h_i)
 */
    static void EQ_I(T* h_o, const T* h_i, size_t n, StreamT st);
    static void EQ_I(T* h_o, const T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * Boolean Not Equal of an array and constant (h_o = h_i != c)
 */
    static void NEQC(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void NEQC(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void NEQ(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void NEQ(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Boolean Not Equal of an array and constant, in-place (h_o = h_o != c)
 */
    static void NEQC_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void NEQC_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void NEQ_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void NEQ_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Boolean Not Equal of two arrays, (h_o = h_i != h_i1)
 */
    static void NEQ(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT st);
    static void NEQ(T* h_o, const T* h_i, const T* h_i1, const T* d_mask, size_t n, StreamT st);

/**
 * Boolean Not Equal of two arrays, in-place (h_o = h_o != h_i)
 */
    static void NEQ_I(T* h_o, const T* h_i, size_t n, StreamT st);
    static void NEQ_I(T* h_o, const T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * Less than of an array and constant (h_o = h_i < c)
 */
    static void LTC(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void LTC(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void LT(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void LT(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Less than of an array and constant, in-place (h_o = h_o < c)
 */
    static void LTC_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void LTC_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void LT_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void LT_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Less than of two arrays, (h_o = h_i < h_i1)
 */
    static void LT(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT st);
    static void LT(T* h_o, const T* h_i, const T* h_i1, const T* d_mask, size_t n, StreamT st);

/**
 * Less than of two arrays, in-place (h_o = h_o < h_i)
 */
    static void LT_I(T* h_o, const T* h_i, size_t n, StreamT st);
    static void LT_I(T* h_o, const T* h_i, const T* d_mask, size_t n, StreamT st);

/**
 * Less than or equal of an array and constant (h_o = h_i <= c)
 */
    static void LTEC(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void LTEC(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void LTE(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);
    static void LTE(T* h_o, const T* h_i, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Less than or equal of an array and constant, in-place (h_o = h_o <= c)
 */
    static void LTEC_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void LTEC_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);
    static void LTE_I(T* h_o, const T& c, size_t n, StreamT st,bool onDev);
    static void LTE_I(T* h_o, const T& c, const T* d_mask, size_t n, StreamT st,bool onDev);

/**
 * Less than or equal of two arrays, (h_o = h_i <= h_i1)
 */
    static void LTE(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT st);
    static void LTE(T* h_o, const T* h_i, const T* h_i1, const T* d_mask, size_t n, StreamT st);

/**
 * Less than or equal of two arrays, in-place (h_o = h_o <= h_i)
 */
    static void LTE_I(T* h_o, const T* h_i, size_t n, StreamT st);
    static void LTE_I(T* h_o, const T* h_i, const T* d_mask, size_t n, StreamT st);


/////////////////////////////////////////////////////////////////////////////////
// Triary functions
//////////////////////////////////////////////////////////////////////////////////

/**
 * Compute the absolution different between 2 vector
 * h_o = abs(h_i - h_i1)
 */
    static void AbsDiff(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT st);

/**
 * Compute the absolution different between 2 vector in-place version
 * h_o = abs(h_o - h_i1)
 */
    static void AbsDiff_I(T* h_o, const T* h_i, size_t n, StreamT st);

/**
 * Compute the absolution different between 2 vector
 * h_o = sqr(h_i - h_i1)
 */

    static void SqrDiff(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT st);

/**
 * Compute the absolution different between 2 vector in-place version
 * h_o = sqr(h_o - h_i1)
 */

    static void SqrDiff_I(T* h_o, const T* h_i, size_t n, StreamT st);

/**
 * h_o = h_i * (h_i1 * c)
 */

    static void MulMulC(T* h_o, const T* h_i, const T* h_i1, const T& c, size_t n, StreamT st,bool onDev);

/**
 * h_o = h_o * (h_i1 * c)
 */

    static void MulMulC_I(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);

/**
 * h_o = h_i * (h_i1 * h_i2)
 */

    static void MulMul(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, size_t n, StreamT st);

/**
 * h_o = h_o * (h_i * h_i1)
 */

    static void MulMul_I(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT st);

/**
 * h_o = h_i * h_i1 + h_i2
 */

    static void MulAdd(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, size_t n, StreamT st);

/**
 * h_o = h_o * h_i + h_i1
 */

    static void MulAdd_I(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT st);

/**
 * h_o = h_i * h_i1 - h_i2
 */

    static void MulSub(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, size_t n, StreamT st);

/**
 * h_o = h_o * h_i - h_i1
 */

    static void MulSub_I(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT st);


/**
 * h_o = (h_i + h_i1) * h_i2
 */

    static void AddMul(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, size_t n, StreamT st);

/**
 * h_o = (h_o + h_i) * h_i1
 */

    static void AddMul_I(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT st);

/**
 * h_o = (h_i + h_i1) / h_i2
 */

    static void AddDiv(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, size_t n, StreamT st);

/**
 * h_o = (h_o + h_i) / h_i1
 */

    static void AddDiv_I(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT st);

/**
 * h_o = (h_i - h_i1) * h_i2
 */

    static void SubMul(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, size_t n, StreamT st);

/**
 * h_o = (h_o - h_i) * h_i1
 */

    static void SubMul_I(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT st);

/**
 * h_o = (h_i - h_i1) / h_i2
 */

    static void SubDiv(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, size_t n, StreamT st);

/**
 * h_o = (h_o - h_i) / h_i1
 */

    static void SubDiv_I(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT st);

/**
 * h_o = (h_i * c) + h_i1
 */

    static void MulCAdd(T* h_o, const T* h_i, const T& c, const T* h_i1, size_t n, StreamT st,bool onDev);

/**
 * h_o = (h_o * c) + h_i
 */

    static void MulCAdd_I(T* h_o, const T& c, const T* h_i, size_t n, StreamT st,bool onDev);

/**
 * h_o = (h_i * c) - h_i1
 */

    static void MulCSub(T* h_o, const T* h_i, const T& c, const T* h_i1, size_t n, StreamT st,bool onDev);

/**
 * h_o = (h_o * c) - h_i
 */

    static void MulCSub_I(T* h_o, const T& c, const T* h_i, size_t n, StreamT st,bool onDev);

/**
 * h_o = (h_i * a) + b
 */

    static void MulCAddC(T* h_o, const T* h_i, const T& a, const T& b, size_t n, StreamT st,bool onDev);

/**
 * h_o = (h_o * a) + b
 */
    static void MulCAddC_I(T* h_o, const T& a, const T& b, size_t n, StreamT st,bool onDev);

/**
 * h_o = (h_i + h_i1) * c
 */
    static void AddMulC(T* h_o, const T* h_i, const T* h_i1, const T& c, size_t n, StreamT st,bool onDev);

/**
 * h_o = (h_o + h_i) * c
 */
    static void AddMulC_I(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);

/**
 * h_o = (h_i - h_i1) * c
 */
    static void SubMulC(T* h_o, const T* h_i, const T* h_i1, const T& c, size_t n, StreamT st,bool onDev);

/**
 * h_o = (h_o - h_i) * c
 */
    static void SubMulC_I(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);


/**
 * h_o = (h_i + a) * b
 */
    static void AddCMulC(T* h_o, const T* h_i, const T& a, const T& b, size_t n, StreamT st,bool onDev);

/**
 * h_o = (h_o + a) * b
 */
    static void AddCMulC_I(T* h_o, const T& a, const T& b, size_t n, StreamT st,bool onDev);

/**
 * h_o = h_i + h_i1 * c
 */
    static void Add_MulC(T* h_o, const T* h_i, const T* h_i1, const T& c, size_t n, StreamT st,bool onDev);

/**
 * h_o = h_o + h_i * c
 */
    static void Add_MulC_I(T* h_o, const T* h_i, const T& c, size_t n, StreamT st,bool onDev);

/**
 * h_o = h_i + h_i1 * h_i2
 */
    static void Add_Mul(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, size_t n, StreamT st);
/**
 * h_o = h_o + h_i * h_i1
 */
    static void Add_Mul_I(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT st);

 /**
 * h_o = h_i - h_i1 * h_i2
 */
    static void Sub_Mul(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, size_t n, StreamT st);
/**
 * h_o = h_o - h_i * h_i1
 */
    static void Sub_Mul_I(T* h_o, const T* h_i, const T* h_i1, size_t n, StreamT st);

/////////////////////////////////////////////////////////////////////////////////
// n-ary functions
//////////////////////////////////////////////////////////////////////////////////
/**
 * h_o = h_i + (h_i1 + h_i2) * c
 */
    static void Add_AddMulC(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, const T& c, size_t n, StreamT st,bool onDev);

/**
 * h_o = h_i + (h_i1 - h_i2) * c
 */
    static void Add_SubMulC(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, const T& c, size_t n, StreamT st,bool onDev);

/**
 * h_o = h_i + (h_i1 * h_i2) * c
 */
    static void Add_MulMulC(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, const T& c, size_t n, StreamT st,bool onDev);


/**
 * h_o += (h_i + h_i1) * c
 */
    static void Add_AddMulC_I(T* h_o, const T* h_i, const T* h_i1, const T& c, size_t n, StreamT st,bool onDev);

/**
 * h_o += (h_i - h_i1) * c
 */
    static void Add_SubMulC_I(T* h_o, const T* h_i, const T* h_i1, const T& c, size_t n, StreamT st,bool onDev);

/**
 * h_o += (h_i * h_i1) * c
 */
    static void Add_MulMulC_I(T* h_o, const T* h_i, const T* h_i1, const T& c, size_t n, StreamT st,bool onDev);

    /**
 * h_o = h_i - (h_i1 + h_i2) * c
 */
    static void Sub_AddMulC(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, const T& c, size_t n, StreamT st,bool onDev);

/**
 * h_o = h_i - (h_i1 - h_i2) * c
 */
    static void Sub_SubMulC(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, const T& c, size_t n, StreamT st,bool onDev);

/**
 * h_o = h_i - (h_i1 * h_i2) * c
 */
    static void Sub_MulMulC(T* h_o, const T* h_i, const T* h_i1, const T* h_i2, const T& c, size_t n, StreamT st,bool onDev);


/**
 * h_o -= (h_i + h_i1) * c
 */
    static void Sub_AddMulC_I(T* h_o, const T* h_i, const T* h_i1, const T& c, size_t n, StreamT st,bool onDev);

/**
 * h_o -= (h_i - h_i1) * c
 */
    static void Sub_SubMulC_I(T* h_o, const T* h_i, const T* h_i1, const T& c, size_t n, StreamT st,bool onDev);

/**
 * h_o -= (h_i * h_i1) * c
 */
    static void Sub_MulMulC_I(T* h_o, const T* h_i, const T* h_i1, const T& c, size_t n, StreamT st,bool onDev);


// Functions with 4 inputs

/**
 * h_o = h_i * a + h_i1 * b
 */
    static void MulC_Add_MulC(T* h_o, const T* h_i, const T& a, const T* h_i1, const T& b, size_t n, StreamT st,bool onDev);

/**
 * h_o = h_o * a + h_i1 * b
 */
    static void MulC_Add_MulC_I(T* h_o, const T& a, const T* h_i, const T& b, size_t n, StreamT st,bool onDev);

/**
 * h_o = (h_i + a) * b + c
 */
    static void AddCMulCAddC(T* h_o, const T* h_i, const T& a, const T& b, const T& c, size_t n, StreamT st,bool onDev);

/**
 * h_o = (h_o + a) * b + c
 */
    static void AddCMulCAddC_I(T* h_o, const T& a, const T& b, const T& c, size_t n, StreamT st,bool onDev);

/**
 * h_o[i] = h_i[n-1 - i]
 */
    static void ReverseOrder(T* h_o, const T* h_i, size_t n, StreamT st);

    static void ShiftCoordinate(T* h_o, const T* h_i, size_t sizeX, size_t sizeY, size_t sizeZ, bool dir, StreamT s);

};

} // end namespace PyCA

#endif

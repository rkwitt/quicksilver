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

#ifndef __CONST_DEVICE_H
#define __CONST_DEVICE_H

#ifndef SWIG
#include <cfloat>
#endif // !SWIG

#define CTA_SIZE              256
#define CUDA_DATA_BLOCK_ALIGN 128
#define BLOCK_ALIGN           128
#define REG_BLOCK_SIZE        256 

#define MAX_NUMBER_DEVICES    256

#define EXEC_CPU       0
#define EXEC_GPU       1
#define EXEC_GPU_PARAM 2

//namespace PyCA {

typedef unsigned int uint;

#define FIX_SCALE_20 1048576.f

// ================================================================
// FFT kernel constant
// ================================================================

#define MAX_FFT_TABLE_LENGTH 512

// ================================================================
// enums for interpolation method and background strategy
// ================================================================

#ifndef DEFAULT_INTERP_METHOD
#define DEFAULT_INTERP_METHOD INTERP_LINEAR
#endif // !DEFAULT_INTERP_METHOD

enum InterpT { INTERP_NN, 
	       INTERP_LINEAR, 
	       INTERP_CUBIC};

enum BackgroundStrategy { BACKGROUND_STRATEGY_PARTIAL_ID,
                          BACKGROUND_STRATEGY_ID,
                          BACKGROUND_STRATEGY_PARTIAL_ZERO,
                          BACKGROUND_STRATEGY_ZERO,
                          BACKGROUND_STRATEGY_CLAMP,
                          BACKGROUND_STRATEGY_WRAP,
			  BACKGROUND_STRATEGY_VAL};

// ================================================================
// defines/enums for FiniteDiff 
// ================================================================

#define WRAP_TRUE 1
#define WRAP_FALSE 0
#define ACCUM_TRUE 1
#define ACCUM_FALSE 0
#define SLICE_TRUE 1
#define SLICE_FALSE 0

// ================================================================
// remove need for __host__ and __device__ if not nvcc
// ================================================================

#ifdef __CUDACC__
#define __HOSTDEVICE__ __host__ __device__
#else
#define __HOSTDEVICE__
#endif

namespace PyCA {

const float MIN_VAL = FLT_MIN;
const float MAX_VAL = FLT_MAX;

// template params for specifying direction 
enum DimT {DIM_X, DIM_Y, DIM_Z};
// finite difference types
enum DiffT {DIFF_FORWARD, DIFF_BACKWARD, DIFF_CENTRAL};
// these only matter for central differences without wrapping --
// BC_APPROX takes the appropriate forward or backward difference at
// the boundary to approximate the central difference, BC_CLAMP
// assumes zero derivative at the boundary. For forward and backward
// differences, zero derivative assumed for either case since no
// reasonable 'approximation' exists.  The original code used the
// approximation boundary conditions, so this is the default
enum BoundaryCondT {BC_APPROX, BC_WRAP, BC_CLAMP};
// operation to perform before returning value.  OP_VAL returns the
// value (no additional operation) OP_SQR returns the square of the
// calculated value.
enum OpT{OP_VAL,OP_SQR};

} // end namespace PyCA

#endif

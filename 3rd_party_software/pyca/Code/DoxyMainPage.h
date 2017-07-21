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

/**
 * \mainpage
 *
 * For  documentation  of the  external  (python)  interface, see  the
 * following documentation sections:
 *   - Objects:
 *       - PyCA::Vec3D - single vector
 *       - PyCA::GridInfo - class representing the size/spacing/origin
 *           of an image or field
 *       - PyCA::Image3D - image object
 *       - PyCA::Field3D - field object
 *       - PyCA::FluidKernelFFT - the fluid kernel class
 *       - PyCA::MultiscaleManager - class to help manage multiple scale levels
 *       - PyCA::MultiscaleResampler - class which uses a
 *           PyCA::MultiscaleManager instance to help up/downsample
 *           PyCA::Image3D and PyCA::Field3D objects
 *       - PyCA::GaussianFilter - Gaussian filter implementation
 *       - PyCA::MemoryManager - Global memory manager for temporary
 *           PyCA::Image3D and PyCA::Field3D objects, currently only
 *           necessary when using PyCA::MultiscaleResampler
 *   - Functions:
 *       - PyCA::Opers - contains most functions in the external interface
 *       - PyCA::CudaUtils - extra functions for querying devices
 *       - PyCA::ITKFileIO - class containing IO functions for
 *           PyCA::Image3D and PyCA::Field3D objects using
 *           ITK-supported file formats (requires compilation with ITK
 *           support)
 */

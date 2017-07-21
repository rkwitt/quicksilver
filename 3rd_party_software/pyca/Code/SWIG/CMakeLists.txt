
# this contains SETUP_PyCA_SWIG_DEFS, which adds -DCUDA_ENABLED,
# -DENABLE_ITK to Core.i
INCLUDE(${PyCA_BINARY_DIR}/PyCAConfig.cmake)

SETUP_PyCA_SWIG_DEFS(Core.i)

# Generate docstring target dependencies

# docstrings for types.i
GEN_DOCSTRING(CLASS Vec3D)
GEN_DOCSTRING(CLASS GridInfo)
GEN_DOCSTRING(CLASS Mat3D)
GEN_DOCSTRING(CLASS Aff3D)
GEN_DOCSTRING(CLASS Image3D)
GEN_DOCSTRING(CLASS Field3D)

# docstrings for math.i
GEN_DOCSTRING(NAMESPACE Opers)

# docstrings for base.i
GEN_DOCSTRING(NAMESPACE CudaUtils)

# docstrings for alg.i
GEN_DOCSTRING(CLASS FluidKernelFFT)
GEN_DOCSTRING(CLASS MemoryManager)
GEN_DOCSTRING(CLASS GaussianFilter)
GEN_DOCSTRING(CLASS MultiscaleManager)
GEN_DOCSTRING(CLASS MultiscaleResampler)

# docstrings for io.i
IF(ENABLE_ITK)
  GEN_DOCSTRING(CLASS ITKFileIO)
ENDIF(ENABLE_ITK)

# set up swig wrapping
GEN_SWIG_MODULE("Core")


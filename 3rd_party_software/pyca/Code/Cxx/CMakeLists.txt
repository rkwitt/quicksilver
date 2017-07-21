
IF(CUDA_ENABLED)

  INCLUDE_DIRECTORIES(
    ${CUDA_TOOLKIT_INCLUDE}
    ${CUDA_SDK_ROOT_DIR}/inc
    ${CUDA_SDK_ROOT_DIR}/common/inc
    )

  SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -arch=sm_${CUDA_ARCH_VERSION} -Xcompiler -D__builtin_stdarg_start=__builtin_va_start)

ENDIF(CUDA_ENABLED)

ADD_SUBDIRECTORY(src)


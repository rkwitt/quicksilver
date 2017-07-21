#-----------------------------------------------------------------------------
# Top level CMakeLists.txt file for PyCA
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(PyCA)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

#-----------------------------------------------------------------------------
# set cmake module path
SET(PyCA_CMAKE_DIR "${PyCA_SOURCE_DIR}/CMake")
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PyCA_CMAKE_DIR}")

#-----------------------------------------------------------------------------
# include extra cmake functions
SET(PyCA_EXTRA_FUNCTIONS ${PyCA_SOURCE_DIR}/PyCAFunctions.cmake)
INCLUDE(${PyCA_EXTRA_FUNCTIONS})

#-----------------------------------------------------------------------------
# set directory for compiled libraries and executables
SET (LIBRARY_OUTPUT_PATH
  ${PROJECT_BINARY_DIR}/bin
  CACHE PATH
  "Single Directory for all libraries.")
SET (EXECUTABLE_OUTPUT_PATH
  ${PROJECT_BINARY_DIR}/bin
  CACHE PATH
  "Single Directory for all executables.")
MARK_AS_ADVANCED (
  LIBRARY_OUTPUT_PATH
  EXECUTABLE_OUTPUT_PATH
)

#-----------------------------------------------------------------------------
# set unix specific variables
IF(UNIX)
ENDIF(UNIX)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# -----------------------------------------------------------------------------
# set up a custom clean target that calls 'make clean', but can also
# have other dependencies
ADD_CUSTOM_TARGET(clean_all
  COMMAND ${CMAKE_MAKE_PROGRAM} clean
  WORKING_DIRECTORY ${CMAKE_CURRENT_DIR}
  )

#-----------------------------------------------------------------------------
# should shared libraries be built?
MARK_AS_ADVANCED(BUILD_SHARED_LIBS)
OPTION(BUILD_SHARED_LIBS "Build shared libraries?" ON)
IF(BUILD_SHARED_LIBS)
  SET(LIBTYPE SHARED)
ELSE(BUILD_SHARED_LIBS)
  SET(LIBTYPE STATIC)
ENDIF(BUILD_SHARED_LIBS)
# compute default library install dir

#-----------------------------------------------------------------------------
# Include debug messages?
IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
  # logDEBUG4 = 7
  ADD_DEFINITIONS(-DERRLOG_MAX_LEVEL=7)
ELSE(CMAKE_BUILD_TYPE STREQUAL "Debug")
  # logINFO = 2
  ADD_DEFINITIONS(-DERRLOG_MAX_LEVEL=2)
ENDIF(CMAKE_BUILD_TYPE STREQUAL "Debug")

#-----------------------------------------------------------------------------
# Define LIBRARY_INSTALL_DIR
set (_DEFAULT_LIBRARY_INSTALL_DIR lib)
if (EXISTS "${CMAKE_INSTALL_PREFIX}/lib32/" AND CMAKE_SIZEOF_VOID_P EQUAL 4)
  set (_DEFAULT_LIBRARY_INSTALL_DIR lib32)
elif (EXISTS "${CMAKE_INSTALL_PREFIX}/lib64/" AND CMAKE_SIZEOF_VOID_P EQUAL 8)
  set (_DEFAULT_LIBRARY_INSTALL_DIR lib64)
endif (EXISTS "${CMAKE_INSTALL_PREFIX}/lib32/" AND CMAKE_SIZEOF_VOID_P EQUAL 4)
# the library install dir
set(LIBRARY_INSTALL_DIR "${_DEFAULT_LIBRARY_INSTALL_DIR}" CACHE PATH "Installation directory for libraries")
# make the library install dir an absolute path (can be important e.g. when using CONFIGURE_FILE to embed
# the library installation directory into a file)
if(NOT IS_ABSOLUTE "${LIBRARY_INSTALL_DIR}")
  set(LIBRARY_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/${LIBRARY_INSTALL_DIR}")
endif(NOT IS_ABSOLUTE "${LIBRARY_INSTALL_DIR}")

# -----------------------------------------------------------------------------
# Insert PyCA version/commit info into
# Code/Python/Common/PyCAVersion.py.in (creating PyCAVersion.py)

# Get git executable
FIND_PACKAGE(Git)

# Force this to run every time 'make' is run
ADD_CUSTOM_TARGET(configure_pyca_version
  COMMAND ${CMAKE_COMMAND} 
  -DPyCA_SOURCE_DIR=${PyCA_SOURCE_DIR} 
  -DPyCA_BINARY_DIR=${PyCA_BINARY_DIR} 
  -P ${PyCA_SOURCE_DIR}/PyCAVersionConfig.cmake
  COMMENT "Configuring ${PyCA_SOURCE_DIR}/Code/Python/Common/PyCAVersion.py"
  VERBATIM)

#-----------------------------------------------------------------------------
# find python interpreter.  In newer versions of cmake FindPythonInterp should work, but for now we'll do it ourselves
find_program(PYTHON_EXECUTABLE 
  NAMES python2.7 python2
  DOC "Python interpreter (version 2.x required)"
  )
mark_as_advanced(PYTHON_EXECUTABLE)

#-----------------------------------------------------------------------------
# If OpenMP is available, always use it by default
include(CheckCCompilerFlag)
check_c_compiler_flag(-fopenmp HAS_FOPENMP)
check_c_compiler_flag(/openmp HAS_SLASHOPENMP)
# gcc style
if (HAS_FOPENMP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
# msvc style
elseif (HAS_SLASHOPENMP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
endif()

#-----------------------------------------------------------------------------
# Setup Doxygen
FIND_PACKAGE(Doxygen)
SET(DOXY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Docs")
SET(DOXY_CONFIG "${DOXY_DIR}/Doxyfile")
IF (DOXYGEN_FOUND)

  SET(DOXYFILE_IN ${PyCA_SOURCE_DIR}/Docs/Doxyfile)
  SET(DOXYFILE_GEN ${PyCA_BINARY_DIR}/Doxyfile)

  # set up a target so we can call 'make doc'
  ADD_CUSTOM_TARGET(doc
    COMMAND cp ${DOXYFILE_IN} ${DOXYFILE_GEN}
    COMMAND echo "INPUT = ${PyCA_SOURCE_DIR}/Code" >> ${DOXYFILE_GEN}
    COMMAND echo "HTML_OUTPUT = ${PyCA_SOURCE_DIR}/Docs/html" >> ${DOXYFILE_GEN}
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_GEN}
    VERBATIM
    )

  # set up a target to remove doxygen-generated files
  ADD_CUSTOM_TARGET(clean_doc
    COMMAND rm -rf ${PyCA_SOURCE_DIR}/Docs/html
    )

  ADD_DEPENDENCIES(clean_all clean_doc)
ENDIF (DOXYGEN_FOUND)

#-----------------------------------------------------------------------------
# targets for static code analysis (delinting)
SET(PYLINT_IGNORE_CODES C0103,R0913,R0914,R0915,W0511)
ADD_CUSTOM_TARGET(delint_python
    COMMAND find ${PyCA_SOURCE_DIR}/{Code,Testing,Examples} -name *.py | xargs pylint --include-ids yes --disable ${PYLINT_IGNORE_CODES}
    )
ADD_CUSTOM_TARGET(delint_cxx
    COMMAND find ${PyCA_SOURCE_DIR}/Code/Cxx -name *.h -o -name *.cxx -o -name *.txx -o -name *.tpl | xargs cpplint.py --counting=detailed
    )

ADD_CUSTOM_TARGET(delint
    DEPENDS delint_python
    DEPENDS delint_cxx)

#-----------------------------------------------------------------------------
# targets for autogenerating docstrings for python from doxygen xml

SET(DOCSTRING_DIR ${PyCA_BINARY_DIR}/docstring)
#MESSAGE(STATUS "Docstring Dir: ${DOCSTRING_DIR}")
LIST(APPEND PyCA_INCLUDE_DIRECTORIES ${DOCSTRING_DIR})

IF(DOXYGEN_FOUND)
  SET(DOXYFILE_XML ${PyCA_BINARY_DIR}/Doxyfile)
  SET(DOXYFILE_XML_DIR ${PyCA_BINARY_DIR}/xml)

  # target used to generate XML parsing of source
  ADD_CUSTOM_TARGET(doxy_xml
    COMMAND cp ${DOXYFILE_IN} ${DOXYFILE_XML}
    COMMAND echo "INPUT = ${PyCA_SOURCE_DIR}/Code" >> ${DOXYFILE_XML}
    COMMAND echo "GENERATE_HTML = NO" >> ${DOXYFILE_XML}
    COMMAND echo "GENERATE_XML = YES" >> ${DOXYFILE_XML}
    COMMAND echo "XML_OUTPUT = ${DOXYFILE_XML_DIR}" >> ${DOXYFILE_XML}
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_XML}
    VERBATIM
    )

  SET(DOXY2SWIG ${PyCA_SOURCE_DIR}/doxy2swig.py)

  # target to generate docstrings (in .i files) from xml
  ADD_CUSTOM_TARGET(docstring)

ENDIF(DOXYGEN_FOUND)

FUNCTION(GEN_DOCSTRING TYPE NAME)

# NOTE: doxygen generates xml file names where PyCA::ClassName is
# written as PyCA_1_1ClassName.xml.  This may break in the future --
# I tried initially with globbing expression, but this fails for
# MemoryManager since there's also a ThreadMemoryManager class.
IF(TYPE STREQUAL "CLASS")
  #SET(XML_FILE_PAT ${DOXYFILE_XML_DIR}/classPyCA*${NAME}.xml)
  SET(XML_FILE_PAT ${DOXYFILE_XML_DIR}/classPyCA_1_1${NAME}.xml)
ELSEIF(TYPE STREQUAL "NAMESPACE")
  #SET(XML_FILE_PAT ${DOXYFILE_XML_DIR}/namespacePyCA*${NAME}.xml)
  SET(XML_FILE_PAT ${DOXYFILE_XML_DIR}/namespacePyCA_1_1${NAME}.xml)
ELSE(TYPE STREQUAL "CLASS")
  MESSAGE(SEND_ERROR "GEN_DOCSTRING: Error, TYPE must be either CLASS or NAMESPACE (recieved ${TYPE})")
ENDIF(TYPE STREQUAL "CLASS")

SET(DOCSTRING_FILE ${DOCSTRING_DIR}/${NAME}.i)
SET(DOCSTRING_FILE_LIST ${DOCSTRING_FILE_LIST} ${DOCSTRING_FILE} PARENT_SCOPE)

IF(DOXYGEN_FOUND)
  SET(DOC_TARG doc_${TYPE}_${NAME}_Targ)
  
  ADD_CUSTOM_TARGET(${DOC_TARG}
    COMMAND mkdir -p ${DOCSTRING_DIR}
    COMMAND ${PYTHON_EXECUTABLE} ${DOXY2SWIG} ${XML_FILE_PAT} ${DOCSTRING_FILE}
    COMMENT "Generating ${DOCSTRING_FILE}"
    # COMMENT "Generating ${DOCSTRING_FILE} from input ${XML_FILE_PAT}"
    DEPENDS doxy_xml
    )
  
  ADD_DEPENDENCIES(docstring ${DOC_TARG})
ENDIF(DOXYGEN_FOUND)

ENDFUNCTION(GEN_DOCSTRING)

#-----------------------------------------------------------------------------
# should we use SWIG to generate python wrappers?
OPTION(WRAP_PYTHON "Generate python wrapping using SWIG?" ON)
IF(WRAP_PYTHON)
  FIND_PACKAGE(SWIG REQUIRED)
  INCLUDE(${SWIG_USE_FILE})

  LIST(APPEND PyCA_INCLUDE_DIRECTORIES ${PyCA_SOURCE_DIR}/Code/SWIG)
  LIST(APPEND PyCA_INCLUDE_DIRECTORIES ${PyCA_SOURCE_DIR}/Code/SWIG/extern)
  
  # NOTE: this may break in future versions of cmake if python version
  # 3.x is found by default
  FIND_PACKAGE(PythonLibs 2 REQUIRED)
  LIST(APPEND PyCA_INCLUDE_DIRECTORIES ${PYTHON_INCLUDE_PATH})
  #message(STATUS "PYTHON_INCLUDE_PATH is ${PYTHON_INCLUDE_PATH}")
  
  # is this still necessary? jsp 2013-03-03
  LIST(APPEND PyCA_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR})

  # find the numpy include dir
  EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c "import numpy;print numpy.get_include()"
      OUTPUT_VARIABLE NPINCIN)
  FILE(TO_CMAKE_PATH "${NPINCIN}" NUMPY_INCLUDE_DIR)

  LIST(APPEND PyCA_INCLUDE_DIRECTORIES ${NUMPY_INCLUDE_DIR})

  SET(CMAKE_SWIG_FLAGS "")

  EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c "import site;print(site.getsitepackages()[0])" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
  FILE(TO_CMAKE_PATH "${PYTHON_SITE_PACKAGES}" DEFAULT_PYTHON_INSTALL_DIR)
  SET(PYTHON_INSTALL_DIR
    ${DEFAULT_PYTHON_INSTALL_DIR}
    CACHE PATH
    "Location to install PyCA python module"
    )

  # ADD_CUSTOM_COMMAND(OUTPUT "${CMAKE_BINARY_DIR}/python_module/timestamp"
  #   COMMAND rm -rf ${CMAKE_BINARY_DIR}/python_module
  #   COMMAND sh ${PyCA_SOURCE_DIR}/GenPythonModule.sh ${CMAKE_BINARY_DIR}
  #   COMMENT "Running python packaging target"
  #   VERBATIM
  #   )

  # dependencies of this target are added by each module in
  # subdirectories
  ADD_CUSTOM_TARGET(py_package ALL
    ${PYTHON_EXECUTABLE} ${PyCA_SOURCE_DIR}/GenPythonModule.py ${PyCA_SOURCE_DIR} ${PyCA_BINARY_DIR} 
    COMMENT "Running package target"
    VERBATIM
    )

  ADD_DEPENDENCIES(py_package configure_pyca_version)

  ADD_CUSTOM_TARGET(clean_py_package
      COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/python_module
    COMMENT "cleaning python package dir"
    VERBATIM
    )
  ADD_DEPENDENCIES(clean_all clean_py_package)

  INSTALL(DIRECTORY ${CMAKE_BINARY_DIR}/python_module/PyCA DESTINATION ${PYTHON_INSTALL_DIR})

  # set up a target so we can call 'make test'
  ADD_CUSTOM_TARGET(test
    COMMAND ${PYTHON_EXECUTABLE} -m unittest discover -v -s ${PyCA_SOURCE_DIR}/Testing -p "*UnitTest.py"
    )

  # function to generate empty file if it doesn't exist.  this is used
  # for docstring files which can just be empty if no doxygen
  # docstrings have been generatedd via 'docstring' target
  FUNCTION(GEN_IF_NOT_EXIST TO_GEN)
  GET_FILENAME_COMPONENT(TO_GEN_ABS ${TO_GEN} ABSOLUTE)
  GET_FILENAME_COMPONENT(TO_GEN_PATH ${TO_GEN_ABS} PATH)
  ADD_CUSTOM_COMMAND(OUTPUT ${TO_GEN}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${TO_GEN_PATH}
    COMMAND ${CMAKE_COMMAND} -E touch ${TO_GEN_ABS}
    COMMENT "Generating ${TO_GEN}"
    )
  ENDFUNCTION(GEN_IF_NOT_EXIST)
  
  # 
  # Function for adding compilation target (_${MODNAME}) for module
  # ${MODNAME} defined in file ${MODNAME}.i
  #
  FUNCTION(GEN_SWIG_MODULE MODNAME)
  
    SET_SOURCE_FILES_PROPERTIES(${MODNAME}.i PROPERTIES CPLUSPLUS ON)
    # add extra options to any existing flags
    ADD_TO_SOURCE_FILE_PROPERTY(${MODNAME}.i SWIG_FLAGS "-includeall")
    # This is no longer needed since we name the module Core (used to
    # be PyCA_Core) 
    # set up proper dependency on target PyCA_${MODNAME}
    # SET_SOURCE_FILES_PROPERTIES(${MODNAME}.i PROPERTIES
    # SWIG_MODULE_NAME PyCA_${MODNAME} ) add extra dependencies
    SET(INC_PATH ${PyCA_SOURCE_DIR}/Code/Cxx/inc/${MODNAME})
    FILE(GLOB SWIG_${MODNAME}_INC_DEPS 
      "${INC_PATH}/*.h" 
      "${INC_PATH}/*.txx" 
      "${INC_PATH}/*.cxx" 
      "${INC_PATH}/*.cu"
      )
    SET(SRC_PATH ${PyCA_SOURCE_DIR}/Code/Cxx/src/${MODNAME})
    FILE(GLOB SWIG_${MODNAME}_SRC_DEPS 
      "${SRC_PATH}/*.h" 
      "${SRC_PATH}/*.txx" 
      "${SRC_PATH}/*.cxx" 
      "${SRC_PATH}/*.cu"
      )

    SET(SWIG_${MODNAME}_SWIGINC_DEPS 
      ${CMAKE_CURRENT_SOURCE_DIR}/types.i
      ${CMAKE_CURRENT_SOURCE_DIR}/math.i
      ${CMAKE_CURRENT_SOURCE_DIR}/base.i
      ${CMAKE_CURRENT_SOURCE_DIR}/io.i
      ${CMAKE_CURRENT_SOURCE_DIR}/alg.i
      ${CMAKE_CURRENT_SOURCE_DIR}/type_opers.i
      ${CMAKE_CURRENT_SOURCE_DIR}/usage_stats.i
      )

    # MESSAGE(STATUS "GEN_SWIG_MODULE ${MODNAME}: DOCSTRING_FILE_LIST: ${DOCSTRING_FILE_LIST}")
    FOREACH(DOCFILE IN LISTS DOCSTRING_FILE_LIST)
      GEN_IF_NOT_EXIST(${DOCFILE})
    ENDFOREACH(DOCFILE)

    SET(SWIG_MODULE_${MODNAME}_EXTRA_DEPS 
      ${SWIG_${MODNAME}_INC_DEPS}
      ${SWIG_${MODNAME}_SRC_DEPS}
      ${DOCSTRING_FILE_LIST}
      ${SWIG_${MODNAME}_SWIGINC_DEPS}
      )
    
    # add module target
    SWIG_ADD_MODULE(${MODNAME} python ${MODNAME}.i)
    SWIG_LINK_LIBRARIES(${MODNAME} 
      PyCA
      ${PYTHON_LIBRARIES}
      )

    # make the 'py_package' target dependent on _${MODNAME}
    ADD_DEPENDENCIES(py_package _${MODNAME})
    # make this module depend on 'check_docstring' target, which just
    # generates some empty .i files if the doxygen-based docstrings
    # haven't been generated
    ADD_DEPENDENCIES(_${MODNAME} check_docstring)

  ENDFUNCTION(GEN_SWIG_MODULE)

ENDIF(WRAP_PYTHON)


#-----------------------------------------------------------------------------
# Set up CUDA
OPTION(USE_CUDA "Attempt to compile CUDA code" ON)
IF(USE_CUDA)

  # Don't need this anymore, FindCUDA.cmake now ships with cmake
  #SET(CMAKE_MODULE_PATH "${PyCA_CMAKE_DIR}/cuda" ${CMAKE_MODULE_PATH})
  FIND_PACKAGE(CUDA)
  IF(NOT CUDA_FOUND)
    MESSAGE("Can't find CUDA, will not compile apps requiring CUDA")
  ENDIF(NOT CUDA_FOUND)
  # As of cuda 5.0, there is no longer a separate SDK, so we don't
  # generally need this
  MARK_AS_ADVANCED(CUDA_SDK_ROOT_DIR)

  # Add definition for preprocessor
  IF(CUDA_FOUND)
    SET(CUDA_ENABLED TRUE)
    ADD_DEFINITIONS(-DCUDA_ENABLED)
    
    SET(CUDA_LIBRARIES 
      ${CUDA_LIBRARIES}
      ${CUDA_CUDA_LIBRARY} 
      ${CUDA_CUBLAS_LIBRARIES}
      ${CUDA_CUFFT_LIBRARIES}
      )

    # Definitions for the cuda version
    IF(${CUDA_VERSION} VERSION_LESS "3.1")
      MESSAGE(SEND_ERROR
	"CUDA version 3.1 or greater required to build PyCA GPU code")
    ENDIF(${CUDA_VERSION} VERSION_LESS "3.1")

    MARK_AS_ADVANCED(CUDA_ARCH_VERSION)
    SET(CUDA_ARCH_VERSION
      "13"
      CACHE STRING
      "Cuda architecture version, ie -arch sm_??")

  ELSE(CUDA_FOUND)
    SET(CUDA_ENABLED FALSE)
  ENDIF(CUDA_FOUND)

ELSE(USE_CUDA)
  # Not using CUDA, make sure this is false
  SET(CUDA_ENABLED FALSE)
ENDIF(USE_CUDA)

#-----------------------------------------------------------------------------
# Boost config
ADD_DEFINITIONS(-DBOOST_ALL_DYN_LINK)
FIND_PACKAGE(Boost REQUIRED)
LIST(APPEND PyCA_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})
# this only needs to be set if boost isn't found in the system library
# locations, but it makes it look like boost wasn't found
MARK_AS_ADVANCED(Boost_DIR)

#-----------------------------------------------------------------------------
# Setup FFTW
# * include directories are added automatically
# * applications must add link libraries ${FFTWF_LIB} or ${FFTWF_THREADS_LIB}
MARK_AS_ADVANCED(USE_FFTWF)
OPTION(USE_FFTWF "Use single precision FFTW if found" ON)
MARK_AS_ADVANCED(USE_FFTWD)
OPTION(USE_FFTWD "Use single precision FFTW if found" ON)
FIND_PACKAGE(FFTW REQUIRED)
LIST(APPEND PyCA_INCLUDE_DIRECTORIES ${FFTW_INCLUDE})

#-----------------------------------------------------------------------------
# Setup ITK
# * include directories are added automatically
# * applications must add appropriate libraries to link libraries
#MARK_AS_ADVANCED(USE_ITK)
OPTION(USE_ITK "Use ITK if found" ON)
IF(USE_ITK)
  FIND_PACKAGE(ITK)
  IF(ITK_FOUND)
    INCLUDE (${ITK_USE_FILE})
    SET(ENABLE_ITK TRUE)
    ADD_DEFINITIONS(-DENABLE_ITK)
  ELSE(ITK_FOUND)
    MESSAGE("Could not find ITK directory.")
    SET(ENABLE_ITK FALSE)
  ENDIF(ITK_FOUND)
ELSE(USE_ITK)
    SET(ENABLE_ITK FALSE)
ENDIF(USE_ITK)

#-----------------------------------------------------------------------------
# add include dirs
IF(CUDA_ENABLED)
  CUDA_INCLUDE_DIRECTORIES(
    ${CUDA_CUT_INCLUDE_DIRS}
    ${CUDA_INCLUDE_DIRS}
    )
  LIST(APPEND PyCA_INCLUDE_DIRECTORIES
    ${CUDA_CUT_INCLUDE_DIRS} 
    ${CUDA_INCLUDE_DIRS}
    )
ENDIF(CUDA_ENABLED)

#-----------------------------------------------------------------------------
# add include directories
LIST(APPEND PyCA_INCLUDE_DIRECTORIES
  ${CMAKE_CURRENT_SOURCE_DIR}/Code/Cxx/inc/alg/
  ${CMAKE_CURRENT_SOURCE_DIR}/Code/Cxx/inc/base/
  ${CMAKE_CURRENT_SOURCE_DIR}/Code/Cxx/inc/const/
  ${CMAKE_CURRENT_SOURCE_DIR}/Code/Cxx/inc/io/
  ${CMAKE_CURRENT_SOURCE_DIR}/Code/Cxx/inc/math/
  ${CMAKE_CURRENT_SOURCE_DIR}/Code/Cxx/inc/thread/
  ${CMAKE_CURRENT_SOURCE_DIR}/Code/Cxx/inc/types/
)

#-----------------------------------------------------------------------------
# set CMake INCLUDE_DIRECTORIES
INCLUDE_DIRECTORIES(${PyCA_INCLUDE_DIRECTORIES})

#-----------------------------------------------------------------------------
# process PyCAConfig.cmake
CONFIGURE_FILE(
  ${PyCA_SOURCE_DIR}/PyCAConfig.cmake.in 
  ${PyCA_BINARY_DIR}/PyCAConfig.cmake 
  @ONLY)

#-----------------------------------------------------------------------------
# go to subdirectories
SUBDIRS(${PyCA_SOURCE_DIR}/Code/Cxx)
IF(WRAP_PYTHON)
  SUBDIRS(${PyCA_SOURCE_DIR}/Code/SWIG)
ENDIF(WRAP_PYTHON)

# hackish fix for linking problem on osx clang
IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    SET(CMAKE_CXX_FLAGS "-stdlib=libstdc++")
ENDIF()

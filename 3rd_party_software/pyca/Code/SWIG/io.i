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

#ifdef ENABLE_ITK
%{
#include "ITKFileIO.h"
%}
#endif

%feature("autodoc","1");

/**
 * Extra exception handling for ITK functions since they can throw
 * non-PyCAException exceptions
 */
%define ExceptionExtension(funcname)
%exception funcname {
  try {
    $action
  } catch (std::exception &e) {
    PyErr_SetString(PyExc_IOError, const_cast<char*>(e.what()));
    return NULL;
  }
}
%enddef

ExceptionExtension(ReadHeader)
ExceptionExtension(LoadImage)
ExceptionExtension(LoadField)
ExceptionExtension(LoadFieldComponents)
ExceptionExtension(SaveImage)
ExceptionExtension(SaveField)
ExceptionExtension(SaveFieldComponents)

#ifdef ENABLE_ITK

 //#import "itkImageIOBase.h"

%rename(_ITKFileIO) ITKFileIO;
// include doxygen docstrings
%include "ITKFileIO.i"

%apply int &OUTPUT{ int &DATA_TYPE_OUT };

%include "ITKFileIO.h"

#endif


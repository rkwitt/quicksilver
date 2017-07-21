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

%module(package="PyCA") Core

%include "std_string.i"
%include "std_vector.i"
%include "std_complex.i"

%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"

%init %{
   import_array();
%}

%feature("autodoc","1");

// print is python keyword
%rename(_print) print;

/**
 * wrap every swig-wrapped c++ function in a try-catch block that
 * catches a PyCaException and raises a python error
 */
%exception {
  try {
    $action
  } catch (PyCA::PyCAException &awe) {
    std::string errmsg("$decl:\n");
    errmsg += awe.what();
    PyErr_SetString(PyExc_RuntimeError, const_cast<char*>(errmsg.c_str()));
    return NULL;
  }
}

%include "types.i"
%include "math.i"
%include "base.i"
%include "alg.i"
%include "io.i"
// add in (limited) operator support
%include "type_opers.i"

// usage statistics reporting on first load of python module
%include "usage_stats.i"

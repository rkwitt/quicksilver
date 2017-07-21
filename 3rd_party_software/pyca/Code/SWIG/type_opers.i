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

 // Vec3D

%define Vec3DOp(ltpl, rtype, op, name)
%extend PyCA::Vec3D< ltpl > {

   Vec3D<ltpl>& PyCA::Vec3D<ltpl>::__##name##I(const rtype &rhs){
      *($self) op##= rhs;
      return *$self;
   }

   Vec3D<ltpl> PyCA::Vec3D<ltpl>::operator##op(const rtype &rhs){
      return *$self op rhs;
   }

};
%enddef

%define Vec3DPyIOpDef(tpl, name)
%extend PyCA::Vec3D< tpl > {
%pythoncode %{
def __i##name##__(self, rhs):
  self.__##name##I(rhs);
  return self;
%}
};
%enddef

%define Vec3DOpDefs(op, name)
Vec3DOp(float, PyCA::Vec3D<float>, op, name)
Vec3DOp(float, PyCA::Vec3D<int>, op, name)
Vec3DOp(float, float, op, name)
Vec3DOp(float, int, op, name)
Vec3DOp(int, PyCA::Vec3D<int>, op, name)
Vec3DOp(int, int, op, name)

Vec3DPyIOpDef(float, name)
Vec3DPyIOpDef(int, name)

%enddef

Vec3DOpDefs(+, add)
Vec3DOpDefs(-, sub)
Vec3DOpDefs(*, mul)
Vec3DOpDefs(/, div)

 // Image3D

 // this didn't work due to some strange object ownership issues
/* %define UnaryOpDef(slfT, rhsT, op, pycaOp) */
/* %extend PyCA::slfT { */
/*    slfT &PyCA::slfT::operator##op##=(const rhsT &rhs){ */
/*       PyCA::Opers::pycaOp##_I(*$self, rhs); */
/*       return *$self; */
/*    } */
/* }; */
/* %enddef */

 // due to problems with object ownership, just implement this in python.
 // First create overloaded member functions for doing addition, etc.
 // then call that function from python e.g. __add__ function

%define UnaryOpMemberDef(slfT, rhsT, memOp, pycaOp)
%extend PyCA::slfT {
   void PyCA::slfT::memOp(const rhsT &rhs){
      PyCA::Opers::pycaOp(*$self, rhs);
   }
};
%enddef

%define UnaryOpPythonDef(slfT, pyOp, pycaOp)
%extend PyCA::slfT {
%pythoncode %{
def pyOp(self, rhs):
  self.pycaOp(rhs);
  return self;
%}
};
%enddef

UnaryOpMemberDef(Image3D, PyCA::Image3D, __addI, Add_I)
UnaryOpMemberDef(Image3D, float, __addI, AddC_I)
UnaryOpMemberDef(Image3D, PyCA::Image3D, __subI, Sub_I)
UnaryOpMemberDef(Image3D, float, __subI, SubC_I)
UnaryOpMemberDef(Image3D, PyCA::Image3D, __mulI, Mul_I)
UnaryOpMemberDef(Image3D, float, __mulI, MulC_I)
UnaryOpMemberDef(Image3D, PyCA::Image3D, __divI, Div_I)
UnaryOpMemberDef(Image3D, float, __divI, DivC_I)

UnaryOpPythonDef(Image3D, __iadd__, __addI)
UnaryOpPythonDef(Image3D, __isub__, __subI)
UnaryOpPythonDef(Image3D, __imul__, __mulI)
UnaryOpPythonDef(Image3D, __idiv__, __divI)

UnaryOpMemberDef(Field3D, PyCA::Field3D, __addI, Add_I)
UnaryOpMemberDef(Field3D, PyCA::Vec3Df, __addI, AddC_I)
UnaryOpMemberDef(Field3D, float, __addI, AddC_I)
UnaryOpMemberDef(Field3D, PyCA::Field3D, __subI, Sub_I)
UnaryOpMemberDef(Field3D, PyCA::Vec3Df, __subI, SubC_I)
UnaryOpMemberDef(Field3D, float, __subI, SubC_I)
UnaryOpMemberDef(Field3D, PyCA::Image3D, __mulI, Mul_I)
UnaryOpMemberDef(Field3D, PyCA::Vec3Df, __mulI, MulC_I)
UnaryOpMemberDef(Field3D, float, __mulI, MulC_I)
UnaryOpMemberDef(Field3D, PyCA::Image3D, __divI, Div_I)
UnaryOpMemberDef(Field3D, PyCA::Vec3Df, __divI, DivC_I)
UnaryOpMemberDef(Field3D, float, __divI, DivC_I)

UnaryOpPythonDef(Field3D, __iadd__, __addI)
UnaryOpPythonDef(Field3D, __isub__, __subI)
UnaryOpPythonDef(Field3D, __imul__, __mulI)
UnaryOpPythonDef(Field3D, __idiv__, __divI)

// Memory usage has been tested, %newobject seems to be handling
// memory ownership properly
//
// note -- standard operator+ doesn't work since copy constructor is
// declared private.  We're directly implementing the python __add__
// function in c++ which maps to the python + operator
%define BinaryOpDef(slfT, rhsT, pyOp, pycaOp)
%extend PyCA::slfT {
   slfT *PyCA::slfT::pyOp(const rhsT &rhs) const {
      PyCA::slfT *rtn = new PyCA::slfT($self->grid(), $self->memType());
      PyCA::Opers::pycaOp(*rtn, *$self, rhs);
      return rtn;
   }
};
%enddef

// these declarations tell SWIG that the named functions return
// objects newly allocated from the heap whose memory python should be
// responsible for freeing.
%newobject __add__;
%newobject __sub__;
%newobject __mul__;
%newobject __div__;
// __truediv__ for Python3
%newobject __truediv__;
%newobject __pow__;

// Image3D

// + operator
BinaryOpDef(Image3D, PyCA::Image3D, __add__, Add)
BinaryOpDef(Image3D, float, __add__, AddC)
BinaryOpDef(Image3D, float, __radd__, AddC)
// - operator
BinaryOpDef(Image3D, PyCA::Image3D, __sub__, Sub)
BinaryOpDef(Image3D, float, __sub__, SubC)
// * operator
BinaryOpDef(Image3D, PyCA::Image3D, __mul__, Mul)
BinaryOpDef(Image3D, float, __mul__, MulC)
BinaryOpDef(Image3D, float, __rmul__, MulC)
// / (div) operator
BinaryOpDef(Image3D, PyCA::Image3D, __div__, Div)
BinaryOpDef(Image3D, float, __div__, DivC)
BinaryOpDef(Image3D, PyCA::Image3D, __truediv__, Div)
BinaryOpDef(Image3D, float, __truediv__, DivC)
// ** (pow)
BinaryOpDef(Image3D, float, __pow__, PowC)

// Field3D

// + operator
BinaryOpDef(Field3D, PyCA::Field3D, __add__, Add)
BinaryOpDef(Field3D, PyCA::Vec3Df, __add__, AddC)
BinaryOpDef(Field3D, float, __add__, AddC)
BinaryOpDef(Field3D, PyCA::Vec3Df, __radd__, AddC)
BinaryOpDef(Field3D, float, __radd__, AddC)
// - operator
BinaryOpDef(Field3D, PyCA::Field3D, __sub__, Sub)
BinaryOpDef(Field3D, PyCA::Vec3Df, __sub__, SubC)
BinaryOpDef(Field3D, float, __sub__, SubC)
// * operator
BinaryOpDef(Field3D, PyCA::Image3D, __mul__, Mul)
BinaryOpDef(Field3D, Vec3Df, __mul__, MulC)
BinaryOpDef(Field3D, float, __mul__, MulC)
BinaryOpDef(Field3D, Vec3Df, __rmul__, MulC)
BinaryOpDef(Field3D, float, __rmul__, MulC)
// / (div) operator
BinaryOpDef(Field3D, PyCA::Image3D, __div__, Div)
BinaryOpDef(Field3D, PyCA::Vec3Df, __div__, DivC)
BinaryOpDef(Field3D, float, __div__, DivC)
BinaryOpDef(Field3D, PyCA::Image3D, __truediv__, Div)
BinaryOpDef(Field3D, PyCA::Vec3Df, __truediv__, DivC)
BinaryOpDef(Field3D, float, __truediv__, DivC)

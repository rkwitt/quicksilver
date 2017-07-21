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

%{
#include <cfloat>
#include "mem.h"
#include "pycaConst.h"
#include "AccumType.h"
#include "MemoryManager.h"
#include "Vec3D.h"
#include "GridInfo.h"
#include "Mat3D.h"
#include "Aff3D.h"
#include "Field3D.h"
#include "Image3D.h"
#include "IOpers.h"
#include "IFOpers.h"
#include "FOpers.h"
%}

// just return false for == comparison
%define LogicalEqualReturnsFalseExtension(type)
%extend type {
%pythoncode %{
def __eq__(self, other):
  """returns false"""
  return False
%}
};
%enddef

// ################################################################
// streaming
// ################################################################

// defines StreamT
%include "estream.h"
// defines BackgroundStrategy
%include "pycaConst.h"

// ################################################################
// mem
// ################################################################

// mem.h defines MemType (MEM_HOST, MEM_DEVICE, etc.)
%include "mem.h"
%include "MemPool.h"

// ################################################################
// copy() functions return new objects whose memory needs to be
// managed by python
// ################################################################

%newobject copy;

// ################################################################
// grid
// ################################################################

// include doxygen docstrings
%include "GridInfo.i"

%include "GridInfo.h"

 // remember %newobject copy declaration, this memory should be
 // managed by python
%extend PyCA::GridInfo {
   PyCA::GridInfo* copy() {
      return new PyCA::GridInfo(*$self);
   }
 };

%extend PyCA::GridInfo {
   char *__str__() {
      static char tmp[1024];
      PyCA::Vec3Di size = $self->size();
      PyCA::Vec3Df origin = $self->origin();
      PyCA::Vec3Df spacing = $self->spacing();
      sprintf(tmp, "Grid: size (%d, %d, %d) origin (%g %g %g) spacing (%g %g %g)", 
	      size.x, size.y, size.z, 
	      origin.x, origin.y, origin.z,
	      spacing.x, spacing.y, spacing.z);
      return tmp;
   }
};

%extend PyCA::GridInfo {
   void clone(const PyCA::GridInfo &other) {
      (*$self) = other;
   }
 };


// ################################################################
// these extensions work for Vec3D, Mat3D, Aff3D
// ################################################################

 // remember %newobject copy declaration, this memory should be
 // managed by python
%define TemplateTypeCopyExtension(classname, typename)
%extend PyCA::classname<typename> {
   PyCA::classname<typename>* copy() {
      return new PyCA::classname<typename>(*$self);
   }
 };
%enddef

// clone functions
%define TemplateTypeClone(classname, typename)
%extend PyCA::classname< typename > {
   void clone(const PyCA::classname< typename > &other) {
      (*$self) = other;
   }
};
%enddef

// ################################################################
// Vec3D
// ################################################################

// AccumType needed for Vec3D
%include "AccumType.h"
%template(AccTi) PyCA::AccumType<int>;
%template(AccTu) PyCA::AccumType<unsigned int>;
%template(AccTf) PyCA::AccumType<float>;
%template(AccTd) PyCA::AccumType<double>;

// include doxygen docstrings
%include "Vec3D.i"

%include "Vec3D.h"

%define Vec3DPrintExtension(format, type)
%extend PyCA::Vec3D<type> {
  char *__str__() {
    static char tmp[1024];
    sprintf(tmp,"Vec3D(format,format,format)", $self->x,$self->y,$self->z);
    return tmp;
  }
};
%enddef

%define Vec3DEqualsExtension(type)
%extend PyCA::Vec3D<type> {
  bool __eq__(const PyCA::Vec3D<type> &other) const {
    return other == *$self;
  }
  bool __ne__(const PyCA::Vec3D<type> &other) const {
    return other != *$self;
  }
};
%enddef

%define Vec3DToList(type)
%extend PyCA::Vec3D< type > {
%pythoncode %{
def tolist(self):
  """Convert a Vec3D to a multidimensional list"""
  vec_list = [self.x, self.y, self.z];
  return vec_list
%}
};
%enddef

%define Vec3DFromList(type)
%extend PyCA::Vec3D< type > {
%pythoncode %{
def fromlist(self, lst):
  """Convert a list to a Vec3D"""
  self.x = lst[0];
  self.y = lst[1];
  self.z = lst[2];
%}
};
%enddef

%template(Vec3Dd) PyCA::Vec3D<double>;
TemplateTypeCopyExtension(Vec3D, double)
TemplateTypeClone(Vec3D, double)
Vec3DPrintExtension(%g,double)
Vec3DToList(double)
Vec3DFromList(double)
Vec3DEqualsExtension(double)

%template(Vec3Df) PyCA::Vec3D<float>;
TemplateTypeCopyExtension(Vec3D, float)
TemplateTypeClone(Vec3D, float)
Vec3DPrintExtension(%g,float)
Vec3DToList(float)
Vec3DFromList(float)
Vec3DEqualsExtension(float)

%template(Vec3Di) PyCA::Vec3D<int>;
TemplateTypeCopyExtension(Vec3D, int)
TemplateTypeClone(Vec3D, int)
Vec3DToList(int)
Vec3DFromList(int)
Vec3DPrintExtension(%d,int)
Vec3DEqualsExtension(int)

%template(Vec3Du) PyCA::Vec3D<unsigned int>;
TemplateTypeCopyExtension(Vec3D, unsigned int)
TemplateTypeClone(Vec3D, unsigned int)
Vec3DToList(unsigned int)
Vec3DFromList(unsigned int)
Vec3DPrintExtension(%d,unsigned int)
Vec3DEqualsExtension(unsigned int)

// ################################################################
// Mat3D
// ################################################################

// include doxygen docstrings
%include "Mat3D.i"

%include "Mat3D.h"

%define Mat3DPrintExtension(format, type)
%extend PyCA::Mat3D<type> {
  char *__str__() {
    static char tmp[1024];
    sprintf(tmp,
	    "Mat3D[format, format, format]\n"
	    "     [format, format, format]\n"
	    "     [format, format, format]\n",
	    $self->get(0,0),$self->get(0,1),$self->get(0,2),
	    $self->get(1,0),$self->get(1,1),$self->get(1,2),
	    $self->get(2,0),$self->get(2,1),$self->get(2,2)
	    );
    return tmp;
  }
};
%enddef

%define Mat3DToList(type)
%extend PyCA::Mat3D< type > {
%pythoncode %{
def tolist(self):
  """Convert a Mat3D to a multidimensional list"""
  mat_list = [[0, 0, 0] for _ in range(3)]
  for m in range(3):
    for n in range(3):
      mat_list[m][n] = self.get(m,n)
  return mat_list
%}
};
%enddef

%define Mat3DFromList(type)
%extend PyCA::Mat3D< type > {
%pythoncode %{
def fromlist(self, lst):
  """Convert a list to a Mat3D"""
  for m in range(3):
    for n in range(3):
      self.set(m,n,lst[m][n])
%}
};
%enddef

%template(Mat3Dd) PyCA::Mat3D<double>;
TemplateTypeCopyExtension(Mat3D, double)
TemplateTypeClone(Mat3D, double)
Mat3DPrintExtension(%g,double)
Mat3DToList(double)
Mat3DFromList(double)

%template(Mat3Df) PyCA::Mat3D<float>;
TemplateTypeCopyExtension(Mat3D, float)
TemplateTypeClone(Mat3D, float)
Mat3DPrintExtension(%g,float)
Mat3DToList(float)
Mat3DFromList(float)

// ################################################################
// Aff3D
// ################################################################

// include doxygen docstrings
%include "Aff3D.i"

%include "Aff3D.h"

%define Aff3DPrintExtension(format, type)
%extend PyCA::Aff3D<type> {
  char *__str__() {
    static char tmp[1024];
    PyCA::Aff3D<type>::MatrixType &m = $self->matrix;
    PyCA::Aff3D<type>::VectorType &v = $self->vector;
    sprintf(tmp,
	    "Aff3D[format, format, format, format]\n"
	    "     [format, format, format, format]\n"
	    "     [format, format, format, format]\n"
	    "     [format, format, format, format]\n",
	    m.get(0,0),m.get(0,1),m.get(0,2),v.get(0),
	    m.get(1,0),m.get(1,1),m.get(1,2),v.get(1),
	    m.get(2,0),m.get(2,1),m.get(2,2),v.get(2),
	    0.f,0.f,0.f,1.f);
    return tmp;
  }
};
%enddef

%define Aff3DToList(type)
%extend PyCA::Aff3D< type > {
%pythoncode %{
def tolist(self):
  """Convert a Aff3D to a multidimensional list"""
  aff_list = [[0, 0, 0, 0] for _ in range(4)]
  mat = self.matrix
  vec = self.vector
  for m in range(3):
    aff_list[m][3] = vec.get(m)
    aff_list[3][m] = 0.0
    for n in range(3):
      aff_list[m][n] = mat.get(m,n)
  aff_list[3][3] = 1.0
  return aff_list
%}
};
%enddef

%define Aff3DFromList(type)
%extend PyCA::Aff3D< type > {
%pythoncode %{
def fromlist(self, lst):
  """Convert a list to a Mat3D"""
  mat = self.matrix
  vec = self.vector
  for m in range(3):
    vec.set(m,lst[m][3])
    for n in range(3):
      mat.set(m,n,lst[m][n])
%}
};
%enddef

%template(Aff3Dd) PyCA::Aff3D<double>;
TemplateTypeCopyExtension(Aff3D, double)
TemplateTypeClone(Aff3D, double)
Aff3DPrintExtension(%g,double)
Aff3DToList(double)
Aff3DFromList(double)

%template(Aff3Df) PyCA::Aff3D<float>;
TemplateTypeCopyExtension(Mat3D, float)
TemplateTypeClone(Aff3D, float)
Aff3DPrintExtension(%g,float)
Aff3DToList(float)
Aff3DFromList(float)

// ################################################################
//  Image3D
// ################################################################

// include doxygen docstrings
%include "Image3D.i"

%include "Image3D.h"

LogicalEqualReturnsFalseExtension(PyCA::Image3D)

 // remember %newobject copy declaration, this memory should be
 // managed by python
%extend PyCA::Image3D {
   PyCA::Image3D* copy() {
      PyCA::Image3D *im = new PyCA::Image3D($self->grid(), $self->memType());
      PyCA::Opers::Copy(*im, *$self);
      return im;
   }
 };

//  print extension 
%extend PyCA::Image3D {
  char *__str__() {
     static const char* memUninitName = "MEM_UNINITIALIZED";
     static const char* memHostName = "MEM_HOST";
     static const char* memDeviceName = "MEM_DEVICE";
     static const char* memHostPinnedName = "MEM_HOST_PINNED";
     static const char* memUnknownName = "MEM_UNKNOWN";
     static char tmp[1024];
     PyCA::Vec3Di size = $self->size();
     PyCA::Vec3Df origin = $self->origin();
     PyCA::Vec3Df spacing = $self->spacing();
     const char* memstr;
     switch($self->memType()){
     case PyCA::MEM_UNINITIALIZED:
	memstr = memUninitName;
	break;
     case PyCA::MEM_HOST:
	memstr = memHostName;
	break;
     case PyCA::MEM_DEVICE:
	memstr = memDeviceName;
	break;
     case PyCA::MEM_HOST_PINNED:
	memstr = memHostPinnedName;
	break;
     default:
	memstr = memUnknownName;
     }
     sprintf(tmp, "Image3D: size (%d, %d, %d) origin (%g %g %g) spacing (%g %g %g) MemType: %s", 
	     size.x, size.y, size.z, 
	     origin.x, origin.y, origin.z,
	     spacing.x, spacing.y, spacing.z,
	     memstr);
     return tmp;
  }
};

// tolist()
%extend PyCA::Image3D {
%pythoncode %{
def tolist(self):
  """Convert a Image3D to a multidimensional list"""
  # transfer to host memory if necessary
  origType = self.memType()
  self.toType(MEM_HOST)
  # set up multidimensional array
  imsize = self.size()
  im_as_list = [0]*imsize.x
  for xi in range(imsize.x):
    im_as_list[xi] = [0]*imsize.y
    for yi in range(imsize.y):
      im_as_list[xi][yi] = [0]*imsize.z
  # fill array
  for zi in range(imsize.z):
    for yi in range(imsize.y):
      for xi in range(imsize.x):
        im_as_list[xi][yi][zi] = self.get(xi,yi,zi)
  # return to original memory
  self.toType(origType)
  return im_as_list
%}
};

// extension function using numpy.i typmaps to give direct access to 
// underlying data as a numpy array
%extend PyCA::Image3D {
   void __asnp(float **ARGOUTVIEW_FARRAY3,
	       int* DIM1, 
	       int* DIM2, 
	       int* DIM3) 
   {
      PyCA::Vec3Di sz = $self->size();
      *DIM1=sz.x;
      *DIM2=sz.y;
      *DIM3=sz.z;
      *ARGOUTVIEW_FARRAY3=$self->get();
   }
   
 };

%extend PyCA::Image3D {
%pythoncode %{
def asnp(self):
  if self.memType() == MEM_DEVICE:
    raise Exception("Cannot wrap device array as numpy array");
    return None
  return self.__asnp()
%}
 };

// fromlist()
%extend PyCA::Image3D {
%pythoncode %{
def fromlist(self, lst):
  """Convert a multidimensional list to a RealArray3D"""
  # convert to host memory if necessary
  origType = self.memType()
  self.toType(MEM_HOST)
  # resize if necessary
  imsize = Vec3Di(len(lst),len(lst[0]),len(lst[0][0]))
  if self.capacity() < imsize.prod():
      self.resize(GridInfo(imsize))
  self.setSize(imsize)
  # copy data
  for zi in range(imsize.z):
    for yi in range(imsize.y):
      for xi in range(imsize.x):
        self.set(xi,yi,zi,lst[xi][yi][zi])
  # return to original memory
  self.toType(origType)
%}
};

%extend PyCA::Image3D {
  float get(unsigned int xIndex,
	    unsigned int yIndex,
	    unsigned int zIndex) const
  {
    float val = $self->get(xIndex, yIndex, zIndex);
    return val;
  }

  float get(unsigned int elementIndex) const
  {
    float val = $self->get(elementIndex);
    return val;
  }
};

%extend PyCA::Image3D {
   unsigned long rawptr() 
   {
       return (unsigned long)$self->get();
   }
   
 };

// ################################################################
// Field3D
// ################################################################


// include doxygen docstrings
%include "Field3D.i"

%include "Field3D.h"

 // remember %newobject copy declaration, this memory should be
 // managed by python
%extend PyCA::Field3D {
   PyCA::Field3D* copy() {
      PyCA::Field3D *f = new PyCA::Field3D($self->grid(), $self->memType());
      PyCA::Opers::Copy(*f, *$self);
      return f;
   }
 };

//  print extension 
%extend PyCA::Field3D {
  char *__str__() {
     static const char* memHostName = "MEM_HOST";
     static const char* memDeviceName = "MEM_DEVICE";
     static const char* memHostPinnedName = "MEM_HOST_PINNED";
     static const char* memUnknownName = "MEM_UNKNOWN";
     static char tmp[1024];
     PyCA::Vec3Di size = $self->size();
     PyCA::Vec3Df origin = $self->origin();
     PyCA::Vec3Df spacing = $self->spacing();
     const char* memstr;
     switch($self->memType()){
     case PyCA::MEM_HOST:
	memstr = memHostName;
	break;
     case PyCA::MEM_DEVICE:
	memstr = memDeviceName;
	break;
     case PyCA::MEM_HOST_PINNED:
	memstr = memHostPinnedName;
	break;
     default:
	memstr = memUnknownName;
     }
     sprintf(tmp, "Field3D: size (%d, %d, %d) origin (%g %g %g) spacing (%g %g %g) MemType: %s", 
	     size.x, size.y, size.z, 
	     origin.x, origin.y, origin.z,
	     spacing.x, spacing.y, spacing.z,
	     memstr);
     return tmp;
  }
};

LogicalEqualReturnsFalseExtension(PyCA::Field3D)

%extend PyCA::Field3D {
   PyCA::Vec3Df   
      get(unsigned int xIndex,
	  unsigned int yIndex,
	  unsigned int zIndex) const
   {
      PyCA::Vec3Df v = $self->get(xIndex, yIndex, zIndex);
      return v;
   }

   PyCA::Vec3Df
      get(unsigned int index) const
   {
      PyCA::Vec3Df v = $self->get(index);
      return v;
   }
};

// tolist()
%extend  PyCA::Field3D {
%pythoncode %{
def tolist(self):
  """Convert a Field3D to a multidimensional list"""
  # convert to host memory if necessary
  origType = self.memType()
  self.toType(MEM_HOST)
  sz = self.size()
  # set up multidimensional array
  vf_as_list = [0]*sz.x
  for xi in range(sz.x):
    vf_as_list[xi] = [0]*sz.y
    for yi in range(sz.y):
      vf_as_list[xi][yi] = [0]*sz.z
      for zi in range(sz.z):
        vf_as_list[xi][yi][zi] = [0]*3
  # fill array
  for zi in range(sz.z):
    for yi in range(sz.y):
      for xi in range(sz.x):
        v = self(xi,yi,zi)
        vf_as_list[xi][yi][zi][0] = v.x
        vf_as_list[xi][yi][zi][1] = v.y
        vf_as_list[xi][yi][zi][2] = v.z
  # convert to original memory type
  self.toType(origType)
  return vf_as_list
%}
};

// extension function using numpy.i typmaps to give direct access to
// underlying data as a numpy arrays
%define Field3DNumpyWrapping(THE_DIM)
%extend PyCA::Field3D {
   void __##THE_DIM##_asnp(float **ARGOUTVIEW_FARRAY3,
			   int* DIM1,
			   int* DIM2,
			   int* DIM3)
   {
      PyCA::Vec3Di sz = $self->size();
      *DIM1=sz.x;
      *DIM2=sz.y;
      *DIM3=sz.z;
      *ARGOUTVIEW_FARRAY3=$self->THE_DIM;
   }
   
};

%extend PyCA::Field3D {
%pythoncode %{
def THE_DIM ## _asnp(self):
  if self.memType() == MEM_DEVICE:
    raise Exception("Cannot wrap device array as numpy array");
    return None
  return self.__##THE_DIM##_asnp()
%}
 };
%enddef

Field3DNumpyWrapping(x)
Field3DNumpyWrapping(y)
Field3DNumpyWrapping(z)

%extend PyCA::Field3D {
%pythoncode %{
def asnp(self):
  x = self.x_asnp()
  y = self.y_asnp()
  z = self.z_asnp()
  return (x,y,z)
%}
};

%extend PyCA::Field3D {
%pythoncode %{
def fromlist(self, l):
  """Convert a multidimensional list to a Array3D"""
  # convert to host memory if necessary
  origType = self.memType()
  self.toType(MEM_HOST)
  # resize if necessary
  sz = Vec3Di(len(l),len(l[0]),len(l[0][0]))
  if self.capacity() < sz.prod():
      self.resize(GridInfo(sz))
  self.grid().setSize(sz)
  for zi in range(sz.z):
    for yi in range(sz.y):
      for xi in range(sz.x):
        v = Vec3Df()
        v.x = l[xi][yi][zi][0]
        v.y = l[xi][yi][zi][1]
        v.z = l[xi][yi][zi][2]
        self.set(xi,yi,zi,v)
  # convert to original memory type
  self.toType(origType)
%}
};

%extend PyCA::Field3D {
   unsigned long rawptr_x() 
   {
       return (unsigned long)$self->getX();
   }
   unsigned long rawptr_y() 
   {
       return (unsigned long)$self->getY();
   }
   unsigned long rawptr_z() 
   {
       return (unsigned long)$self->getZ();
   }
   
 };

// ################################################################
//  MemoryManager
// ################################################################
/**
 * Need ThreadMemoryManager for MutiscaleResampler.  Note that in
 * order for ThreadMemoryManager singleton will only work properly in
 * the library it's wrapped in -- the static var in PyCA_alg will be
 * different from the one initialized if this was wrapped in another
 * library.  Compiling a single dynamic 'PyCA' shared library and then
 * linking it (dynamically) against all _PyCA_*.so python libs would
 * probably fix this, but for now just including this here works, as
 * only things in _PyCA_alg require the ThreadMemoryManager.
 */

 // include doxygen docstrings
%include "MemoryManager.i"
 // wrap file
%include "MemoryManager.h"



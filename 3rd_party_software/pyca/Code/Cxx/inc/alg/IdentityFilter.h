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

#ifndef __IDENTITY_FILTER_H__
#define __IDENTITY_FILTER_H__

#ifndef SWIG

#include <Vec3D.h>
#include <estream.h>
#include <Image3D.h>
#include <Field3D.h>

#endif // SWIG

#include <PyCAException.h>

namespace PyCA {

template<int mode>
class IdentityFilter 
{
public:
    enum { exec_mode = mode };

    IdentityFilter(){}

    virtual ~IdentityFilter(){}

    void updateParams(const Vec3Di& size, 
		      const Vec3Df& sig, 
		      const Vec3Di kRad)
    {
	mSize = size;
    }

   virtual void filter(Image3D &a_o, const Image3D &a_i, Image3D &a_tmp, StreamT stream=NULL)
   {
       MK_CHECK2_ALL(a_o, a_i);
       if(exec_mode == EXEC_CPU && a_o.memType() != MEM_HOST){
	   throw PyCAException(__FILE__,__LINE__,"Using EXEC_CPU object on non-host memory");
       }else if(exec_mode == EXEC_GPU && a_o.memType() == MEM_HOST){
	   throw PyCAException(__FILE__,__LINE__,"Using EXEC_GPU object on host memory");
       }
       this->filter(a_o.get(), a_i.get(), a_tmp.get(), stream);
   }
    
    virtual void filter(float *a_o, const float* a_i, float* a_t, StreamT stream=NULL)
    {
	size_t nVox = mSize.prod();
	if(exec_mode == EXEC_CPU){
	    cpyArrayH2H(a_o, a_i, nVox);
	}else{
	    acpyArrayD2D(a_o, a_i, nVox, stream);
	}
    }

    virtual void filter2(float *a_o, const float* a_i, float* a_t, StreamT stream=NULL)
    {
	size_t nVox = mSize.prod();

	if(exec_mode == EXEC_CPU){
	    cpyArrayH2H(a_o, a_i, nVox);
	}else{
	    acpyArrayD2D(a_o, a_i, nVox, stream);
	}
    }

    virtual void convolutionSingleAxis(float* a_o, const float *a_i,
                                       size_t sizeX, size_t sizeY, size_t sizeZ,
                                       int axis, StreamT stream)
    {
	size_t nVox = mSize.prod();

	if(exec_mode == EXEC_CPU){
	    cpyArrayH2H(a_o, a_i, nVox);
	}else{
	    acpyArrayD2D(a_o, a_i, nVox, stream);
	}
    }

protected:
    Vec3Di mSize;

};

} // end namespace PyCA

#endif // __IDENTITY_FILTER_H__

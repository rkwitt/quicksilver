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

#ifndef __GAUSSIAN_FILTER_H__
#define __GAUSSIAN_FILTER_H__

#include <GaussianFilterInterface.h>

#ifndef SWIG

#include <Vec3D.h>
#include <estream.h>
#include <Image3D.h>
#include <Field3D.h>

#include <GaussianFilterBase.h>
#include <GGaussianFilter.h>
// PyCA
#include <Selector.h>

#endif // SWIG

namespace PyCA {

template<int mode>
class GaussianFilter 
    : public GaussianFilterInterface
{
public:
    enum { exec_mode = mode };

    typedef 
    typename Selector<mode, 
			  GaussianFilterBase<EXEC_CPU>, 
			  GGaussianFilter, 
			  GGaussianFilter>::Result 
			  GaussianFilterExecT;
    
    GaussianFilter(){}

    virtual ~GaussianFilter(){}

    void updateParams(const Vec3Di& size, 
		      const Vec3Df& sig, 
		      const Vec3Di kRad)
    {
	mGaussianFilterExec.updateParams(size, sig, kRad);
    }

   virtual void filter(Image3D &a_o, const Image3D &a_i, Image3D &a_tmp, StreamT stream=NULL)
   {
      this->filter(a_o.get(), a_i.get(), a_tmp.get(), stream);
   }
    
    virtual void filter(float *a_o, const float* a_i, float* a_t, StreamT stream=NULL)
    {
	mGaussianFilterExec.filter(a_o, a_i, a_t, stream);
    }

    virtual void filter2(float *a_o, const float* a_i, float* a_t, StreamT stream=NULL)
    {
	mGaussianFilterExec.filter2(a_o, a_i, a_t, stream);
    }

    virtual void convolutionSingleAxis(float* a_o, const float *a_i,
                                       size_t sizeX, size_t sizeY, size_t sizeZ,
                                       int axis, StreamT stream)
    {
	mGaussianFilterExec.convolutionSingleAxis(a_o, a_i, 
						  sizeX, sizeY, sizeZ,
						  axis, stream);
    }

protected:
    GaussianFilterExecT mGaussianFilterExec;
    
};

} // end namespace PyCA

#endif // __GAUSSIAN_FILTER_H__

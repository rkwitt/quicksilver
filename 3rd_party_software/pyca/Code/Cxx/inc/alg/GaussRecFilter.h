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

#ifndef __GAUSS_REC_FILTER_H__
#define __GAUSS_REC_FILTER_H__

#include <GaussRecFilterInterface.h>

#ifndef SWIG

#include <Vec3D.h>
#include <estream.h>
#include <Image3D.h>
#include <Field3D.h>

#include <GaussRecFilterBase.h>
#include <GGaussRecFilter.h>
// PyCA
#include <Selector.h>

#endif // SWIG

namespace PyCA {

template<int mode>
class GaussRecFilter 
   : public GaussRecFilterInterface
{
public:
   enum { exec_mode = mode };

   typedef 
   typename Selector<mode, 
		     GaussRecFilterBase<EXEC_CPU>, 
		     GGaussRecFilter, 
		     GGaussRecFilter>::Result 
   GaussRecFilterExecT;
    
   GaussRecFilter(){}

   virtual ~GaussRecFilter(){}
   
   void updateParams(const Vec3Di& size, 
		     const Vec3Df& sig)
   {
      mGaussRecFilterExec.updateParams(size, sig);
   }
   
   virtual void filter(Image3D &a_o, const Image3D &a_i, Image3D &a_tmp, StreamT stream=NULL)
   {
      this->filter(a_o.get(), a_i.get(), a_tmp.get(), stream);
   }
   
   virtual void filter(float *a_o, const float* a_i, float* a_t, StreamT stream=NULL)
   {
      mGaussRecFilterExec.filter(a_o, a_i, a_t, stream);
   }
   
   virtual void convolutionSingleAxis(float* a_o, const float *a_i,
				      size_t sizeX, size_t sizeY, size_t sizeZ,
				      int axis, StreamT stream)
   {
      mGaussRecFilterExec.convolutionSingleAxis(a_o, a_i, 
						sizeX, sizeY, sizeZ,
						axis, stream);
   }
   
protected:
   GaussRecFilterExecT mGaussRecFilterExec;
   
};

} // end namespace PyCA

#endif // __GAUSS_REC_FILTER_H__



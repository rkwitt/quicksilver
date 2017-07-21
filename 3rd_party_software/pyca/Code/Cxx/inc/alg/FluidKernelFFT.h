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

#ifndef __FLUID_KERNEL_FFT_H__
#define __FLUID_KERNEL_FFT_H__

#include <FluidKernelFFTInterface.h>

#ifndef SWIG

#include <Vec3D.h>
#include <estream.h>

#include <CFluidKernelFFT.h>
#include <GFluidKernelFFT.h>

#include <Selector.h>

#endif // SWIG

namespace PyCA {

/**
 * DiffOper implements the differential operator:
 *
 * \f[
 * \mathbf{L} = -\alpha\nabla^2 - \beta\nabla(\nabla\cdot) + \gamma\mathbf{I}
 * \f]
 *
 * (where \f$\mathbf{I}\f$ is the identity matrix) used in fluid
 * registration.
 *
 * This operator comes from fluid mechanics, specifically a
 * simplification of the Navier-Stokes equation for compressible
 * fluids with a very low Reynolds number. The parameters \f$\alpha\f$
 * and \f$\beta\f$ control the viscosity of the fluid, while
 * \f$\gamma\f$ ensures that \f$\mathbf{L}\f$ is invertable.
 * Furthermore, the LPow parameter allows the operator
 * \f$\mathbf{L}^p\f$ to be used for values of p other than one --
 * theoretically this will cause greater smoothing of the deformation
 * fields.  If eigenvalue/vector lookup tables are used,
 * precomputation is done (SetUseEigenLUT()) gaining speed at the
 * expense of memory.
 *
 * This class is not safe for use with multple threads.
 *
 */
template<int mode>
class FluidKernelFFT :
    public FluidKernelFFTInterface
{
public:
    enum { exec_mode = mode };

    typedef 
    typename Selector<mode, 
			  CFluidKernelFFT, 
			  GFluidKernelFFT, 
			  GFluidKernelFFT>::Result 
    FluidKernelFFTExecT;

   void setAlpha(float alpha){ 
      mFluidKernelFFTExec.setAlpha(alpha);
   }
   float getAlpha(){ 
      return mFluidKernelFFTExec.getAlpha();
   }
   void setBeta(float beta){ 
      mFluidKernelFFTExec.setBeta(beta);
   }
   float getBeta(){ 
      return mFluidKernelFFTExec.getBeta();
   }
   void setGamma(float gamma){ 
      mFluidKernelFFTExec.setGamma(gamma);
   }
   float getGamma(){ 
      return mFluidKernelFFTExec.getGamma();
   }
   void setLPow(int lpow){ 
      mFluidKernelFFTExec.setLPow(lpow);
   }
   int getLPow(){ 
      return mFluidKernelFFTExec.getLPow();
   }
   void setDivergenceFree(bool divfree){ 
      mFluidKernelFFTExec.setDivergenceFree(divfree);
   }
   bool getDivergenceFree(){ 
      return mFluidKernelFFTExec.getDivergenceFree();
   }

    void setGrid(const GridInfo& grid)
    {
       mFluidKernelFFTExec.setGrid(grid);
    }

    MemoryType memType() const {
	return mFluidKernelFFTExec.memType();
    }
    
    
    /**
     * f = Lv
     * 
     * v field is overwritten in this operation
     */
    void applyOperator(Field3D& f, const Field3D& v, 
		       StreamT stream = NULL)
    {
	mFluidKernelFFTExec.applyOperator(f, v, stream);
    }


    /**
     * v = Kf
     * 
     * f field is overwritten in this operation
     */
    void applyInverseOperator(Field3D& v, const Field3D& f, 
			      StreamT stream = NULL)
    {
	mFluidKernelFFTExec.applyInverseOperator(v, f, stream);
    }

    
    /**
     * f = Lv
     * 
     * v field is overwritten in this operation (holds f).
     */
    void applyOperator(Field3D& f, StreamT stream = NULL)
    {
	mFluidKernelFFTExec.applyOperator(f, stream);
    }


    /**
     * v = Kf
     * 
     * f field is overwritten in this operation (holds v).
     */
    void applyInverseOperator(Field3D& v, StreamT stream = NULL)
    {
	mFluidKernelFFTExec.applyInverseOperator(v, stream);
    }
    

protected:
    FluidKernelFFTExecT mFluidKernelFFTExec;
    
};

} // end namespace PyCA

#endif // __FLUID_KERNEL_FFT_H__

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

#include <GaussRecFilterBase.h>
#include <conditionMacro.h>
#include "SeparableFilter.h"
#include <GaussUtils.h>

namespace PyCA {

template<int mode, int order> 
void computeGaussRecParams(typename GaussRecFilterBase<mode>::GaussRecParams& p, float sigma)
{
    float a0, a1, a2, a3;
    float b1, b2;
    float coefp, coefn;

    // compute filter coefficients
    const float  nsigma = sigma < 0.1f ? 0.1f : sigma;
    const float  alpha  = 1.695f / nsigma;
    const float  ema    = (float)std::exp(-alpha);
    const float  ema2   = (float)std::exp(-2*alpha);
    
    b1 = -2*ema, b2 = ema2;
    a0 = 0, a1 = 0, a2 = 0, a3 = 0, coefp = 0, coefn = 0;
    switch (order) {
        case 0: {
            const float k = (1-ema)*(1-ema)/(1+2*alpha*ema-ema2);
            a0 = k;
            a1 = k*(alpha-1)*ema;
            a2 = k*(alpha+1)*ema;
            a3 = -k*ema2;
        } break;

        case 1: {
            const float k = (1-ema)*(1-ema)/ema;
            a0 = k*ema;
            a1 = a3 = 0;
            a2 = -a0;
        } break;

        case 2: {
            const float
                ea = (float)std::exp(-alpha),
                k = -(ema2-1)/(2*alpha*ema),
                kn = (-2*(-1+3*ea-3*ea*ea+ea*ea*ea)/(3*ea+1+3*ea*ea+ea*ea*ea));
            a0 = kn;
            a1 = -kn*(1+k*alpha)*ema;
            a2 = kn*(1-k*alpha)*ema;
            a3 = -kn*ema2;
        } break;

        default:
            throw PyCAException("gaussianFilter: invalid order parameter!\n");
            return;
    }
    coefp = (a0+a1)/(1+b1+b2);
    coefn = (a2+a3)/(1+b1+b2);

    p.a0 = a0;p.a1 = a1;p.a2 = a2;p.a3 = a3;
    p.b1 = b1;p.b2 = b2;
    p.coefn=coefn;p.coefp=coefp;
}

template<int mode>
GaussRecFilterBase<mode>::
GaussRecFilterBase()
   : mSize(0,0,0), mSigma(0.f,0.f,0.f) 
{}

template<int mode>
void
GaussRecFilterBase<mode>::
updateParams(const Vec3Di& size, const Vec3Df& sigma){
   mSize = size;
   if (mSigma.x != sigma.x)
      computeGaussRecParams<exec_mode, 0>(mParams[0], sigma.x);
   if (mSigma.y != sigma.y)
      computeGaussRecParams<exec_mode, 0>(mParams[1], sigma.y);
   if (mSigma.z != sigma.z)
      computeGaussRecParams<exec_mode, 0>(mParams[2], sigma.z);
   mSigma = sigma;
}

template<int mode>
void
GaussRecFilterBase<mode>::
filter(float *a_o, const float* a_i, float* a_t, StreamT stream)
{
    SeparableFilter<GaussRecFilterBase, EXEC_CPU>(*this, a_o, a_i, a_t, mSize, stream);
}

template<int mode>
void
GaussRecFilterBase<mode>::
ConvolutionX3D(float* a_o, const float* a_i, const GaussRecParams& p,
               size_t sizeX, size_t sizeY, size_t sizeZ, StreamT stream)
{
    float a0 = p.a0, a1 = p.a1, a2 = p.a2, a3 = p.a3;
    float b1 = p.b1, b2 = p.b2;
    float coefn=p.coefn, coefp=p.coefp;

    const uint planeSize = sizeX * sizeY;
    bool clampToEdge = true;

    size_t id=0;
    
    for (size_t y=0; y < sizeY; ++y)
        for (size_t x=0; x < sizeX; ++x, ++id) {
            float xp = 0.f, yp = 0.f, yb = 0.f;

            size_t nId = id;
            
            if (clampToEdge) {
                xp = a_i[nId];  yb = coefp*xp;   yp = yb;
            }

            for (size_t z = 0; z < sizeZ; z++, nId +=planeSize){
                float xc = a_i[nId];
                float yc = a0*xc + a1*xp - b1*yp - b2*yb;
                a_o[nId] = yc;

                
                //shifting around input output 
                xp = xc; yb = yp; yp = yc;
            }

            // reset pointers to point to last element in column
            nId -= planeSize;
            
            // reverse pass
            // ensures response is symmetrical
            float xn = 0.0f, xa = 0.0f, yn = 0.0f, ya = 0.0f;
    
            if (clampToEdge){
                xn = xa = a_i[nId]; yn = coefn*xn; ya = yn;
            }

            for (int z = sizeZ-1; z >= 0; z--, nId -=planeSize) {
                float xc = a_i[nId];
                float yc = a2*xn + a3*xa - b1*yn - b2*ya;
                a_o[nId] = a_o[nId] + yc;

                //shifting around input output 
                xa = xn;
                xn = xc;
                ya = yn;
                yn = yc;
            }
        }
}

template<int mode>
void
GaussRecFilterBase<mode>::
convolutionSingleAxis(float* a_o, const float *a_i,
                      size_t sizeX, size_t sizeY, size_t sizeZ,
                      int axis, StreamT stream)
{
   PRECONDITION(sizeX * sizeY * sizeZ == (size_t)mSize.prod(), "Incompatible size");
    ConvolutionX3D(a_o, a_i, mParams[axis], sizeX, sizeY, sizeZ, stream);
}

// template instantiation
template class GaussRecFilterBase<EXEC_CPU>;
#ifdef CUDA_ENABLED
template class GaussRecFilterBase<EXEC_GPU>;
#endif

} // end namespace PyCA

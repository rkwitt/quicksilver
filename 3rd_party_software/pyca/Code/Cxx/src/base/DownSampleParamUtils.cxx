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

#include <DownSampleParamUtils.h>
#include <GaussParamUtils.h>

namespace PyCA {

namespace DownSampleParamUtils{
    //
    // int versions
    //
    float SigmaFromScale(size_t f){
        // return sqrt(((float)f)/2.f);
        return (float)f;
    }

    Vec3Df SigmaFromScale(const Vec3Di& f){
        // return Vec3Df(sqrt(((float)f.x)/2.f),
        //               sqrt(((float)f.y)/2.f),
        //               sqrt(((float)f.z)/2.f));
       return Vec3Df((float)f.x, (float)f.y, (float)f.z);
    }

    void SigmaFromScale(Vec3Df& sig, size_t scale){
	float s = SigmaFromScale(scale);
	sig = Vec3Df(s,s,s);
    }

    void SigmaRadiusFromScale(float& sigma, size_t& kRadius, size_t f){
        sigma   = SigmaFromScale(f);
        kRadius = GaussParamUtils::RadiusFromSigma(sigma);
    }
    
    void SigmaRadiusFromScale(Vec3Df& sigma, Vec3Di& kRadius, int f){
        float sig = SigmaFromScale((size_t)f);
        int   kR  = GaussParamUtils::RadiusFromSigma(sig);
        sigma   = Vec3Df(sig, sig, sig);
        kRadius = Vec3Di(kR, kR, kR);
    }

    void SigmaRadiusFromScale(Vec3Df& sigma, Vec3Di& kRadius, const Vec3Di& factor)
    {
        sigma   = SigmaFromScale(factor);
        kRadius = GaussParamUtils::RadiusFromSigma(sigma);
    }

    //
    // float versions
    //
    float SigmaFromScale(float f){
        return sqrt(f/2.f);
    }

    Vec3Df SigmaFromScale(const Vec3Df& f){
        return Vec3Df(sqrt(f.x/2.f),
                      sqrt(f.y/2.f),
                      sqrt(f.z/2.f));
    }

    void SigmaRadiusFromScale(float& sigma, size_t& kRadius, float f){
        sigma   = SigmaFromScale(f);
        kRadius = GaussParamUtils::RadiusFromSigma(sigma);
    }
    
    void SigmaRadiusFromScale(Vec3Df& sigma, Vec3Di& kRadius, float f){
        float sig = SigmaFromScale(f);
        int   kR  = GaussParamUtils::RadiusFromSigma(sig);
        sigma   = Vec3Df(sig, sig, sig);
        kRadius = Vec3Di(kR, kR, kR);
    }

    void SigmaRadiusFromScale(Vec3Df& sigma, Vec3Di& kRadius, const Vec3Df& factor)
    {
        sigma   = SigmaFromScale(factor);
        kRadius = GaussParamUtils::RadiusFromSigma(sigma);
    }
}
} // end namespace PyCA

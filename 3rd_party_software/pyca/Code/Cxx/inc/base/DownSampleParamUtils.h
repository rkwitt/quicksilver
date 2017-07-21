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

#ifndef __DOWNSAMPLE_PARAMS_UTILS_H
#define __DOWNSAMPLE_PARAMS_UTILS_H

#include <Vec3D.h>

namespace PyCA {

namespace DownSampleParamUtils{

// int versions
float SigmaFromScale(size_t f);
Vec3Df SigmaFromScale(const Vec3Di& f);
void SigmaFromScale(Vec3Df& f, size_t scale);

void SigmaRadiusFromScale(float& sigma, size_t& kRadius, size_t f);
void SigmaRadiusFromScale(Vec3Df& sigma, Vec3Di& kRadius, int f);
void SigmaRadiusFromScale(Vec3Df& sigma, Vec3Di& kRadius, const Vec3Di& factor);
// float versions
float SigmaFromScale(float f);
Vec3Df SigmaFromScale(const Vec3Df& f);
  
void SigmaRadiusFromScale(float& sigma, size_t& kRadius, float f);
void SigmaRadiusFromScale(Vec3Df& sigma, Vec3Di& kRadius, float f);
void SigmaRadiusFromScale(Vec3Df& sigma, Vec3Di& kRadius, const Vec3Df& factor);

} // end namespace DownSampleParamUtils

} // end namespace PyCA

#endif

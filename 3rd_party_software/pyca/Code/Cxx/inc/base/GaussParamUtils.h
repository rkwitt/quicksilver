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

#ifndef __GAUSS_PARAM_UTILS_H
#define __GAUSS_PARAM_UTILS_H

#include <Vec3D.h>

namespace PyCA {

namespace GaussParamUtils {

size_t RadiusFromSigma(float sigma);
Vec3Di RadiusFromSigma(const Vec3Df& sigma);

template<class ParamT>
ParamT GetParamFromSigma(float sigma);

}; // end namespace GaussParamUtils

} // end namespace PyCA

#endif

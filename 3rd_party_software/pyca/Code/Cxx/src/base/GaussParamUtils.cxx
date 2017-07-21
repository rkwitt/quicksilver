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

#include <GaussParamUtils.h>

namespace PyCA {

namespace GaussParamUtils {
    size_t RadiusFromSigma(float sigma){
        return 2*static_cast<int>(std::ceil(sigma));
    }

    Vec3Di RadiusFromSigma(const Vec3Df& sigma){
        return Vec3Di(2*static_cast<int>(std::ceil(sigma.x)),
                      2*static_cast<int>(std::ceil(sigma.y)),
                      2*static_cast<int>(std::ceil(sigma.z)));
    }
}
} // end namespace PyCA

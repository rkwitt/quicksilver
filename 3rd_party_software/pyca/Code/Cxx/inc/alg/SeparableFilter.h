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

#ifndef __SEPARABLE_FILTER_H
#define __SEPARABLE_FILTER_H

#include <pycaConst.h>
#include <estream.h>
#include <Vec3D.h>

namespace PyCA {

template<class FilterPolicy, int mode>
void SeparableFilter(FilterPolicy& impl,
                     float *a_o, const float* a_i, float* a_temp,
                     const Vec3Di& size, StreamT stream);

template<class FilterPolicy, int mode>
void SeparableFilter(FilterPolicy& impl,
                     float *a_o, float* a_temp,
                     const Vec3Di& size, StreamT stream);

} // end namespace PyCA

#include "SeparableFilter.txx"

#endif

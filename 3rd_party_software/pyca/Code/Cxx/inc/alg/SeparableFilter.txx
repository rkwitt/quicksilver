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

#include "SeparableFilter.h"

#include <MemOpers.h>

namespace PyCA {

template<class FilterPolicy, int mod>
void SeparableFilter(FilterPolicy& impl,
                     float *a_o, const float* a_i, float* a_temp,
                     const Vec3Di& size, StreamT stream){
    
    impl.convolutionSingleAxis(a_temp,a_i,size.x,size.y,size.z,0,stream);
    MemOpers<mod,float>::ShiftCoordinate(a_o,a_temp,size.x,size.y,size.z,true,stream);
    
    impl.convolutionSingleAxis(a_temp,a_o,size.z,size.x,size.y,1,stream);
    MemOpers<mod,float>::ShiftCoordinate(a_o,a_temp,size.z,size.x,size.y,true,stream);
    
    impl.convolutionSingleAxis(a_temp,a_o,size.y,size.z,size.x,2,stream);
    MemOpers<mod,float>::ShiftCoordinate(a_o,a_temp,size.y,size.z,size.x,true,stream);
}

template<class FilterPolicy, int mod>
void SeparableFilter(FilterPolicy& impl,
                     float *a_o, float* a_temp,
                     const Vec3Di& size, StreamT stream){
    
    impl.convolutionSingleAxis(a_temp,a_o,size.x,size.y,size.z,0,stream);
    MemOpers<mod,float>::ShiftCoordinate(a_o,a_temp,size.x,size.y,size.z,true,stream);
    
    impl.convolutionSingleAxis(a_temp,a_o,size.z,size.x,size.y,1,stream);
    MemOpers<mod,float>::ShiftCoordinate(a_o,a_temp,size.z,size.x,size.y,true,stream);
    
    impl.convolutionSingleAxis(a_temp,a_o,size.y,size.z,size.x,2,stream);
    MemOpers<mod,float>::ShiftCoordinate(a_o,a_temp,size.y,size.z,size.x,true,stream);
}

} // end namespace PyCA

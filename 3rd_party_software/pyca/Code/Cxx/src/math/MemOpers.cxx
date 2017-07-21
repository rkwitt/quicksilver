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

#include <MemOpers.h>

#include <MemPool.h>

namespace PyCA {

template<typename T>
void Opers 
::Copy(MemPool<T>& a_o, const MemPool<T>& a_i, int nVox, StreamT stream){
    bool on_host   = ((a_o.memType()!= MEM_DEVICE) && (a_i.memType()!= MEM_DEVICE));
    bool on_device = ((a_o.memType()== MEM_DEVICE) && (a_i.memType()== MEM_DEVICE));

    if (stream != NULL ) {
        PRECONDITION(on_device ||
                     a_o.memType() == MEM_HOST_PINNED ||
                     a_i.memType() == MEM_HOST_PINNED,
                     "Pinned memory is required for streaming between CPU and GPU");
    }
    
    if (nVox < 0)
        nVox   = std::min(a_i.capacity(), a_o.capacity());
    if (on_host){
       cpyArrayH2H(a_o.get(), a_i.get(), nVox);
    } else if (on_device){
       acpyArrayD2D(a_o.get(), a_i.get(), nVox, stream);
    } else {  // uploading or downloading
        if (a_o.memType() == MEM_DEVICE) {
            (stream == NULL) ? cpyArrayH2D(a_o.get(), a_i.get(), nVox) :
                               acpyArrayH2D(a_o.get(), a_i.get(), nVox, stream);
        }
        else {
            (stream == NULL) ? cpyArrayD2H(a_o.get(), a_i.get(), nVox) :
                               acpyArrayD2H(a_o.get(), a_i.get(), nVox, stream);
        }
    }
}

// instantiate this for the only actual voxel type we're likely to use (float)
template void Opers::Copy(MemPool<float>& a_o, const MemPool<float>& a_i,
                          int nVox, StreamT stream);

} // end namespace PyCA

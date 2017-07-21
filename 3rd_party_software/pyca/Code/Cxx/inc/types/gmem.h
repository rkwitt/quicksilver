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

#ifndef __GMEM_H
#define __GMEM_H

#include <estream.h>

namespace PyCA {

// Prototypes
/*
 * Create primiteve object in device memory and return a pointer to the object
 * with a initial value from the host
 */
template<typename T>
T* createGObj(const T& v=T());

/*
 * Setvalue of a primitive object with a initial value from the host
 */
template<typename T>
void setGObj(T* d_obj, const T& h_v);

/*
 * Read the value of the primitive object to hostmemory
 */
template<typename T>
T getGObj(T* d_obj);

/*
 * Perform the division with float object in the device (or host memory)
 * a_o[0] = a_v / a_o[0]
 */
// template<int mode>
// void IDiv(float* a_o, float a_v, StreamT s);

} // end namespace PyCA

#endif

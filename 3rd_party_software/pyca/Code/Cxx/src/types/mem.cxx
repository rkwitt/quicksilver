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

#include <mem.h>
#include <iostream>
#include <string.h>

namespace PyCA {

const char* MemTypeStr(MemoryType type) {
   switch (type){
       case 0:
           return "uninitialized memory";
       case 1:
           return "host memory";
       case 2:
           return "pinned host memory";
       case 3:
           return "gpu_device_memory";
       default:
           return "unknown memory type";
   }
}

} // end namespace PyCA

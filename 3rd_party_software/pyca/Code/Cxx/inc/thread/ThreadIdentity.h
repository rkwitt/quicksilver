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

#ifndef __THREAD_IDENDITY_H
#define __THREAD_IDENDITY_H

#include<stddef.h>

namespace PyCA {

/*
 * Function should be called to initialize the thread id at the beginning of
 * each thread
 */
void InitThreadId();
/*
 * Return the thread id of the current thread
 */
size_t GetThreadId();

} // end namespace PyCA

#endif

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

#ifndef __CONDITION_MACRO_H
#define __CONDITION_MACRO_H

#include <PyCAException.h>
#include <cassert>

#define PRECONDITION(cond, msg)						\
    do                                                                  \
    {                                                                   \
        if (!(cond))                                                    \
        {                                                               \
            throw PyCA::PyCAException(__FILE__, __LINE__, msg);	\
        }                                                               \
    } while(0);

#define PYCA_ASSERT(cond)						         \
    do                                                                   \
    {                                                                    \
        if (!(cond))                                                     \
        {                                                                \
            throw PyCA::PyCAException(__FILE__, __LINE__, #cond); \
        }                                                                \
    } while(0);

#define MK_CHECK2_SIZE(a, b) PRECONDITION((a.grid() == b.grid()), "Input grids are not equal")
#define MK_CHECK3_SIZE(a, b, c) PRECONDITION((a.grid() == b.grid())&&(a.grid() == c.grid()), "Input grids are not equal")
#define MK_CHECK4_SIZE(a, b, c, d) PRECONDITION((a.grid() == b.grid())&&(a.grid() == c.grid())&&(a.grid() == d.grid()), "Input grids are not equal")

#define MK_CHECK_REFSIZE_1(s, a) PRECONDITION(s == a.size(),"Inputs have different Size")
#define MK_CHECK_REFSIZE_2(s, a, b) PRECONDITION((s == a.size())&&(s == b.size()),"Inputs have different Size")
#define MK_CHECK_REFSIZE_3(s, a, b, c) PRECONDITION((s == a.size())&&(s == b.size())&&(s == c.size()),"Inputs have different Size")
#define MK_CHECK_REFSIZE_4(s, a, b, c, d) PRECONDITION((s == a.size())&&(s == b.size())&&(s == c.size())&&(s == d.size()),"Inputs have different Size")
#define MK_CHECK_REFSIZE_5(s, a, b, c, d, e) PRECONDITION((s == a.size())&&(s == b.size())&&(s == c.size())&&(s == d.size())&&(s == e.size()),"Inputs have different Size")

// note: 
//       BACKGROUND_STRATEGY_* are from base/pycaConst.h

#ifdef _MSC_VER
#define STATIC_ASSERT assert
#else
#define STATIC_ASSERT(x)                        \
	typedef char StaticAssert[(x) ? 1 : -1];
#endif

#define MK_CHECK_VFIELD_BACKGROUND(bg) STATIC_ASSERT(bg == BACKGROUND_STRATEGY_ZERO || bg == BACKGROUND_STRATEGY_PARTIAL_ZERO || bg == BACKGROUND_STRATEGY_CLAMP || bg == BACKGROUND_STRATEGY_WRAP)

#define MK_CHECK_HFIELD_BACKGROUND(bg) STATIC_ASSERT(bg == BACKGROUND_STRATEGY_ID || bg == BACKGROUND_STRATEGY_PARTIAL_ID || bg == BACKGROUND_STRATEGY_CLAMP)

#define MK_CHECK_IMAGE_BACKGROUND(bg) STATIC_ASSERT(bg == BACKGROUND_STRATEGY_ZERO || bg == BACKGROUND_STRATEGY_PARTIAL_ZERO || bg == BACKGROUND_STRATEGY_CLAMP || bg == BACKGROUND_STRATEGY_WRAP || bg == BACKGROUND_STRATEGY_VAL)

// note: 
//       EXEC_CPU, EXEC_GPU, EXEC_GPU_PARAM are from base/pycaConst.h
//       MEM_HOST, MEM_DEVICE, MEM_HOST_PINNED are from inc/types/mem.h

#define CHECK_VALID_MEM_EXEC(exec_mode, mem_type)			\
    PRECONDITION( (exec_mode == EXEC_CPU && mem_type == PyCA::MEM_DEVICE), \
		  "Cannot call CPU methods on device memory"));		\
    PRECONDITION( (exec_mode == EXEC_GPU && mem_type == PyCA::MEM_HOST),	\
		  "Cannot call GPU methods on host memory"));		\
    PRECONDITION( (exec_mode == EXEC_GPU_PARAM && mem_type == PyCA::MEM_HOST), \
		  "Cannot call GPU methods on host memory"));

#define MK_CHECK2_MEM(a, b)						\
    PRECONDITION((a.memType() == b.memType()),				\
		 "memory types are not compatible");
#define MK_CHECK3_MEM(a, b, c)						\
    MK_CHECK2_MEM(a,b)							\
    MK_CHECK2_MEM(b,c)							\
    MK_CHECK2_MEM(a,c)

#define MK_CHECK4_MEM(a, b, c, d)					\
    MK_CHECK2_MEM(a,b)							\
    MK_CHECK2_MEM(a,c)							\
    MK_CHECK2_MEM(a,d)							\
    MK_CHECK2_MEM(b,c)							\
    MK_CHECK2_MEM(b,d)							\
    MK_CHECK2_MEM(c,d)
    
#define MK_CHECK2_ALL(a, b)						      \
   MK_CHECK2_SIZE(a,b)                                                \
   MK_CHECK2_MEM(a,b)

#define MK_CHECK3_ALL(a, b, c)					      \
   MK_CHECK3_SIZE(a,b,c)						      \
   MK_CHECK3_MEM(a,b,c)

#define MK_CHECK4_ALL(a, b, c, d)					      \
   MK_CHECK4_SIZE(a,b,c,d)					      \
   MK_CHECK4_MEM(a,b,c,d)

#endif

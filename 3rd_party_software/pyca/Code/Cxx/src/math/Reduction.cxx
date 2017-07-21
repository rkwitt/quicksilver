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

#include <Reduction.h>
#include <ReduceOb.h>
#include <Image3D.h>
#include <Field3D.h>
#include <pycaConst.h>
#include <estream.h>

namespace PyCA {

#ifdef CUDA_ENABLED

// execute correct version of a function based on memory type and
// onDev pram
#define AUTO_EXEC_ONDEV(mem_type, the_class, the_function)	        \
    {									\
	if(mem_type == MEM_DEVICE ||				        \
	   mem_type == MEM_HOST_PINNED){				\
	   if(!onDev){                                                  \
              the_class<EXEC_GPU>::the_function;			\
           }else{                                                       \
              the_class<EXEC_GPU_PARAM>::the_function;			\
           }                                                            \
	}else{								\
	    the_class<EXEC_CPU>::the_function;			        \
	}								\
    }

#else // CUDA_ENABLED

// execute correct version of a function based on memory type and
// onDev pram
#define AUTO_EXEC_ONDEV(mem_type, the_class, the_function)	        \
    {									\
	if(mem_type == MEM_DEVICE ||				        \
	   mem_type == MEM_HOST_PINNED){				\
           throw PyCAException(__FILE__, __LINE__,                      \
                                  "Error, GPU code not compiled");      \
	}else{								\
	    the_class<EXEC_CPU>::the_function;			        \
	}								\
    }

#endif // CUDA_ENABLED

namespace Opers {

void Max(float& a_o,const Image3D& a_i,bool update,StreamT stream,bool onDev)
{
   AUTO_EXEC_ONDEV(a_i.memType(), ReduceOb, Max(a_o, a_i, update, stream, onDev));
}

void Min(float& a_o,const Image3D& a_i,bool update,StreamT stream,bool onDev)
{
   AUTO_EXEC_ONDEV(a_i.memType(), ReduceOb, Min(a_o, a_i, update, stream, onDev));
}
void Sum(float& a_o,const Image3D& a_i,bool update,StreamT stream,bool onDev)
{
   AUTO_EXEC_ONDEV(a_i.memType(), ReduceOb, Sum(a_o, a_i, update, stream, onDev));
}
void LInf(float& a_o,const Image3D& a_i,bool update,StreamT stream,bool onDev)
{
   AUTO_EXEC_ONDEV(a_i.memType(), ReduceOb, LInf(a_o, a_i, update, stream, onDev));
}
void L1(float& a_o,const Image3D& a_i,bool update,StreamT stream,bool onDev) 
{
   AUTO_EXEC_ONDEV(a_i.memType(), ReduceOb, L1(a_o, a_i, update, stream, onDev));
}
void Sum2(float& a_o,const Image3D& a_i,
	  bool update,StreamT stream,bool onDev) 
{
   AUTO_EXEC_ONDEV(a_i.memType(), ReduceOb, Sum2(a_o, a_i, update, stream, onDev));
}
void Dot(float& a_o,const Image3D& a_i,const Image3D& a_i1,
	 bool update,StreamT stream,bool onDev) 
{
   MK_CHECK2_ALL(a_i, a_i1)
   AUTO_EXEC_ONDEV(a_i.memType(), ReduceOb, Dot(a_o, a_i, a_i1, update, stream, onDev));
}
void MaxMin(Vec2D<float>& a_o,const Image3D& a_i,
	    bool update,StreamT stream,bool onDev)
{
   AUTO_EXEC_ONDEV(a_i.memType(), ReduceOb, MaxMin(a_o, a_i, update, stream, onDev));
}

void Max(float& a_o,const Field3D& a_i, bool update,StreamT stream,bool onDev)
{
   AUTO_EXEC_ONDEV(a_i.memType(), ReduceOb, Max(a_o, a_i, update, stream, onDev));
}

void Min(float& a_o,const Field3D& a_i, bool update,StreamT stream,bool onDev)
{
   AUTO_EXEC_ONDEV(a_i.memType(), ReduceOb, Min(a_o, a_i, update, stream, onDev));
}

void Sum(float& a_o,const Field3D& a_i, bool update,StreamT stream,bool onDev)
{
   AUTO_EXEC_ONDEV(a_i.memType(), ReduceOb, Sum(a_o, a_i, update, stream, onDev));
}

void Sum(Vec3D<float>& a_o,const Field3D& a_i, StreamT stream)
{
   AUTO_EXEC(a_i.memType(), ReduceOb, Sum(a_o, a_i, stream));
}

void LInf(float& a_o,const Field3D& a_i, bool update,StreamT stream,bool onDev)
{
   AUTO_EXEC_ONDEV(a_i.memType(), ReduceOb, LInf(a_o, a_i, update, stream, onDev));
}
void L1(float& a_o,const Field3D& a_i, bool update,StreamT stream,bool onDev)
{
   AUTO_EXEC_ONDEV(a_i.memType(), ReduceOb, L1(a_o, a_i, update, stream, onDev));
}

void Sum2(float& a_o,const Field3D& a_i, bool update,StreamT stream,bool onDev)
{
   AUTO_EXEC_ONDEV(a_i.memType(), ReduceOb, Sum2(a_o, a_i, update, stream, onDev));
}

void MaxMin(Vec2D<float>& a_o,const Field3D& a_i, bool update,StreamT stream,bool onDev)
{
   AUTO_EXEC_ONDEV(a_i.memType(), ReduceOb, MaxMin(a_o, a_i, update, stream, onDev));
}

void Dot(float& a_o,const Field3D& a_i,const Field3D& a_i1,bool update,StreamT stream,bool onDev)
{
   MK_CHECK2_ALL(a_i, a_i1)
   AUTO_EXEC_ONDEV(a_i.memType(), ReduceOb, Dot(a_o, a_i, a_i1, update, stream, onDev));
}

} // end namespace Opers
   
} // end namespace PyCA

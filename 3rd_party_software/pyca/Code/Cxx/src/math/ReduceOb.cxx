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

#include <ReduceOb.h>
#include <Reduce.h>
#include <Image3D.h>
#include <Field3D.h>

namespace PyCA {

template<int mode>
void ReduceOb<mode>::
Max(float& a_o,const Image3D& a_i,
    bool update,StreamT stream,bool onDev)
{
   if(onDev) PYCA_ASSERT(mode == EXEC_GPU_PARAM);
   size_t n = a_i.nVox();
   Reduce<mode>::Max(a_o, a_i.get(), n, update, stream);
}

template<int mode>
void ReduceOb<mode>::
Min(float& a_o,const Image3D& a_i,
    bool update,StreamT stream,bool onDev)
{
   if(onDev) PYCA_ASSERT(mode == EXEC_GPU_PARAM);
   size_t n = a_i.nVox();
   Reduce<mode>::Min(a_o, a_i.get(), n, update, stream);
}

template<int mode>
void ReduceOb<mode>::
Sum(float& a_o,const Image3D& a_i,
    bool update,StreamT stream,bool onDev)
{
   if(onDev) PYCA_ASSERT(mode == EXEC_GPU_PARAM);
   size_t n = a_i.nVox();
   Reduce<mode>::Sum(a_o, a_i.get(), n, update, stream);
}

template<int mode>
void ReduceOb<mode>::
LInf(float& a_o,const Image3D& a_i,
     bool update,StreamT stream,bool onDev)
{
   if(onDev) PYCA_ASSERT(mode == EXEC_GPU_PARAM);
   size_t n = a_i.nVox();
   Reduce<mode>::LInf(a_o, a_i.get(), n, update, stream);
}


template<int mode>
void ReduceOb<mode>::
L1(float& a_o,const Image3D& a_i,
   bool update,StreamT stream,bool onDev) 
{
   if(onDev) PYCA_ASSERT(mode == EXEC_GPU_PARAM);
   size_t n = a_i.nVox();
   Reduce<mode>::L1(a_o, a_i.get(), n, update, stream);
}

template<int mode>
void ReduceOb<mode>::
Sum2(float& a_o,const Image3D& a_i,
     bool update,StreamT stream,bool onDev) 
{
   if(onDev) PYCA_ASSERT(mode == EXEC_GPU_PARAM);
   size_t n = a_i.nVox();
   Reduce<mode>::Sum2(a_o, a_i.get(), n, update, stream);
}

template<int mode>
void ReduceOb<mode>::
Dot(float& a_o,const Image3D& a_i,const Image3D& a_i1,
    bool update,StreamT stream,bool onDev) 
{
   if(onDev) PYCA_ASSERT(mode == EXEC_GPU_PARAM);
   size_t n = a_i.nVox();
   Reduce<mode>::Dot(a_o, a_i.get(), a_i1.get(), n, update, stream);
}

template<int mode>
void ReduceOb<mode>::
MaxMin(Vec2D<float>& a_o,const Image3D& a_i, 
       bool update,StreamT stream,bool onDev)
{
   if(onDev) PYCA_ASSERT(mode == EXEC_GPU_PARAM);
   size_t n = a_i.nVox();
   Reduce<mode>::MaxMin(a_o, a_i.get(), n, update, stream);
}

template<int mode>
void ReduceOb<mode>::
Max(float& a_o,const Field3D& a_i, 
    bool update, StreamT stream,bool onDev)
{
   if(onDev) PYCA_ASSERT(mode == EXEC_GPU_PARAM);
   size_t n = a_i.nVox();
   if (a_i.isContinuous()) {
      Reduce<mode>::Max(a_o, a_i.x, 3 * n, update, stream);
   }else{
      Reduce<mode>::Max(a_o, a_i.x, n, update, stream);
      Reduce<mode>::Max(a_o, a_i.y, n, true, stream);
      Reduce<mode>::Max(a_o, a_i.z, n, true, stream);
   }
}

template<int mode>
void ReduceOb<mode>::
Min(float& a_o,const Field3D& a_i, 
    bool update, StreamT stream,bool onDev)
{
   if(onDev) PYCA_ASSERT(mode == EXEC_GPU_PARAM);
   size_t n = a_i.nVox();
   if (a_i.isContinuous()) {
      Reduce<mode>::Min(a_o, a_i.x, 3 * n, update, stream);
   }else{
      Reduce<mode>::Min(a_o, a_i.x, n, update, stream);
      Reduce<mode>::Min(a_o, a_i.y, n, true, stream);
      Reduce<mode>::Min(a_o, a_i.z, n, true, stream);
   }
}

template<int mode>
void ReduceOb<mode>::
Sum(float& a_o,const Field3D& a_i, 
    bool update, StreamT stream,bool onDev)
{
   if(onDev) PYCA_ASSERT(mode == EXEC_GPU_PARAM);
   size_t n = a_i.nVox();
   if (a_i.isContinuous()) {
      Reduce<mode>::Sum(a_o, a_i.x, 3 * n, update, stream);
   }else{
      Reduce<mode>::Sum(a_o, a_i.x, n, update, stream);
      Reduce<mode>::Sum(a_o, a_i.y, n, true, stream);
      Reduce<mode>::Sum(a_o, a_i.z, n, true, stream);
   }
}

template<int mode>
void ReduceOb<mode>::
Sum(Vec3D<float>& a_o,const Field3D& a_i, StreamT stream)
{
   size_t n = a_i.nVox();
   Reduce<mode>::Sum(a_o.x, a_i.x, n, false, stream);
   Reduce<mode>::Sum(a_o.y, a_i.y, n, false, stream);
   Reduce<mode>::Sum(a_o.z, a_i.z, n, false, stream);
}

template<int mode>
void ReduceOb<mode>::
LInf(float& a_o,const Field3D& a_i, 
     bool update, StreamT stream,bool onDev)
{
   if(onDev) PYCA_ASSERT(mode == EXEC_GPU_PARAM);
   size_t n = a_i.nVox();
   if (a_i.isContinuous()) {
      Reduce<mode>::LInf(a_o, a_i.x, 3 * n, update, stream);
   }else{
      Reduce<mode>::LInf(a_o, a_i.x, n, update, stream);
      Reduce<mode>::LInf(a_o, a_i.y, n, true, stream);
      Reduce<mode>::LInf(a_o, a_i.z, n, true, stream);
   }
}

template<int mode>
void ReduceOb<mode>::
L1(float& a_o,const Field3D& a_i, 
   bool update, StreamT stream,bool onDev)
{
   if(onDev) PYCA_ASSERT(mode == EXEC_GPU_PARAM);
   size_t n = a_i.nVox();
   if (a_i.isContinuous()) {
      Reduce<mode>::L1(a_o, a_i.x, 3 * n, update, stream);
   }else{
      Reduce<mode>::L1(a_o, a_i.x, n, update, stream);
      Reduce<mode>::L1(a_o, a_i.y, n, true, stream);
      Reduce<mode>::L1(a_o, a_i.z, n, true, stream);
   }
}

template<int mode>
void ReduceOb<mode>::
Sum2(float& a_o,const Field3D& a_i, 
     bool update, StreamT stream,bool onDev)
{
   if(onDev) PYCA_ASSERT(mode == EXEC_GPU_PARAM);
   size_t n = a_i.nVox();
   if (a_i.isContinuous()) {
      Reduce<mode>::Sum2(a_o, a_i.x, 3 * n, update, stream);
   }else{
      Reduce<mode>::Sum2(a_o, a_i.x, n, update, stream);
      Reduce<mode>::Sum2(a_o, a_i.y, n, true, stream);
      Reduce<mode>::Sum2(a_o, a_i.z, n, true, stream);
   }
}

template<int mode>
void ReduceOb<mode>::
MaxMin(Vec2D<float>& a_o,const Field3D& a_i, 
       bool update,StreamT stream,bool onDev)
{
   if(onDev) PYCA_ASSERT(mode == EXEC_GPU_PARAM);
   size_t n = a_i.nVox();
   if (a_i.isContinuous()) {
      Reduce<mode>::MaxMin(a_o, a_i.x, 3 * n, update, stream);
   }else{
      Reduce<mode>::MaxMin(a_o, a_i.x, n, update, stream);
      Reduce<mode>::MaxMin(a_o, a_i.y, n, true, stream);
      Reduce<mode>::MaxMin(a_o, a_i.z, n, true, stream);
   }
}


template<int mode>
void ReduceOb<mode>::
Dot(float& a_o, const Field3D& a_i, const Field3D& a_i1,
    bool update,StreamT stream,bool onDev)
{
   if(onDev) PYCA_ASSERT(mode == EXEC_GPU_PARAM);
   size_t n = a_i.nVox();
   if (a_i.isContinuous() && a_i1.isContinuous()) {
      Reduce<mode>::Dot(a_o, a_i.x, a_i1.x, 3 * n, update, stream);
   }else{
      Reduce<mode>::Dot(a_o, a_i.x, a_i1.x, n, update, stream);
      Reduce<mode>::Dot(a_o, a_i.y, a_i1.y, n, true, stream);
      Reduce<mode>::Dot(a_o, a_i.z, a_i1.z, n, true, stream);
   }
   
}

// template instantiation
#include "ReduceOb_inst.cxx"

} // end namespace PyCA

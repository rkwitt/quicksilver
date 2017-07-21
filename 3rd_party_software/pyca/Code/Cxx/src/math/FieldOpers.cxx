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

#include <FieldOpers.h>
#include "CFieldOpers.h"
#include "GFieldOpers.h"
#include <ImageFieldOpers.h>
#include <conditionMacro.h>
#include <Field3D.h>
#include <Image3D.h>

#define FIELD_BINARY_ARRAY_OPC_VEC_DEF(OP)     	       	       	       	   \
template<int mode>			  				   \
void FieldOpers<mode>::			  				   \
OP##C(Field3D& a_o, const Field3D& a_i, const Vec3Df& c,		   \
     StreamT stream, bool onDev)	  				   \
{					  				   \
    MK_CHECK2_SIZE(a_o, a_i);		  				   \
    MK_CHECK2_MEM(a_o, a_i);		  				   \
    size_t n = a_o.nVox();		  				   \
    MemOpers<mode, float>::OP##C(a_o.x, a_i.x, c.x, n, stream, onDev);	   \
    MemOpers<mode, float>::OP##C(a_o.y, a_i.y, c.y, n, stream, onDev);	   \
    MemOpers<mode, float>::OP##C(a_o.z, a_i.z, c.z, n, stream, onDev); 	   \
}

#define FIELD_BINARY_ARRAY_OPC_FLOAT_DEF(OP)   	       	       	       	   \
template<int mode>			  				   \
void FieldOpers<mode>::							   \
OP##C(Field3D& a_o, const Field3D& a_i, const float& c,			   \
     StreamT stream, bool onDev)					   \
{									   \
    MK_CHECK2_SIZE(a_o, a_i);						   \
    MK_CHECK2_MEM(a_o, a_i);						   \
    size_t n = a_o.nVox();						   \
    if (Is1D(a_o, a_i, n)){						   \
        MemOpers<mode, float>::OP##C(a_o.x, a_i.x, c, 3 * n,		   \
				    stream, onDev);			   \
        return;								   \
    }									   \
    MemOpers<mode, float>::OP##C(a_o.x, a_i.x, c, n, stream, onDev);	   \
    MemOpers<mode, float>::OP##C(a_o.y, a_i.y, c, n, stream, onDev);	   \
    MemOpers<mode, float>::OP##C(a_o.z, a_i.z, c, n, stream, onDev);	   \
}

#define FIELD_BINARY_ARRAY_OPC_I_VEC_DEF(OP)				   \
template<int mode>							   \
void FieldOpers<mode>::							   \
OP##C_I(Field3D& a_o, const Vec3Df& c,					   \
       StreamT stream, bool onDev)					   \
{									   \
    size_t n = a_o.nVox();						   \
    MemOpers<mode, float>::OP##C_I(a_o.x, c.x, n, stream, onDev);	   \
    MemOpers<mode, float>::OP##C_I(a_o.y, c.y, n, stream, onDev);	   \
    MemOpers<mode, float>::OP##C_I(a_o.z, c.z, n, stream, onDev);	   \
}

#define FIELD_BINARY_ARRAY_OPC_I_FLOAT_DEF(OP)				   \
template<int mode>							   \
void FieldOpers<mode>::							   \
OP##C_I(Field3D& a_o, const float& c,					   \
       StreamT stream, bool onDev)					   \
{									   \
    size_t n = a_o.nVox();						   \
    if (Is1D(a_o, n)){							   \
        MemOpers<mode, float>::OP##C_I(a_o.x, c, 3 * n,			   \
				      stream, onDev);			   \
        return;								   \
    }									   \
    MemOpers<mode, float>::OP##C_I(a_o.x, c, n, stream, onDev);		   \
    MemOpers<mode, float>::OP##C_I(a_o.y, c, n, stream, onDev);		   \
    MemOpers<mode, float>::OP##C_I(a_o.z, c, n, stream, onDev);		   \
}

#define FIELD_BINARY_ARRAY_OP_FIELD_DEF(OP)				   \
template<int mode>							   \
void FieldOpers<mode>::							   \
OP(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1,		   \
    StreamT stream)							   \
{									   \
    MK_CHECK3_SIZE(a_o, a_i, a_i1);					   \
    MK_CHECK3_MEM(a_o, a_i, a_i1);					   \
    size_t n = a_o.nVox();						   \
    if (Is1D(a_o, a_i, a_i1, n)){					   \
        MemOpers<mode, float>::OP(a_o.x, a_i.x, a_i1.x, 3 * n,		   \
				   stream);				   \
        return;								   \
    }									   \
    MemOpers<mode, float>::OP(a_o.x, a_i.x, a_i1.x, n, stream);		   \
    MemOpers<mode, float>::OP(a_o.y, a_i.y, a_i1.y, n, stream);		   \
    MemOpers<mode, float>::OP(a_o.z, a_i.z, a_i1.z, n, stream);		   \
}

#define FIELD_BINARY_ARRAY_OP_I_FIELD_DEF(OP)				   \
template<int mode>							   \
void FieldOpers<mode>::							   \
OP##_I(Field3D& a_o, const Field3D& a_i,				   \
      StreamT stream)							   \
{									   \
    MK_CHECK2_SIZE(a_o, a_i);						   \
    MK_CHECK2_MEM(a_o, a_i);						   \
    size_t n = a_o.nVox();						   \
									   \
    if (Is1D(a_o, a_i, n)){						   \
        MemOpers<mode, float>::OP##_I(a_o.x, a_i.x, 3 * n,		   \
				     stream);				   \
        return;								   \
    }									   \
    MemOpers<mode, float>::OP##_I(a_o.x, a_i.x, n, stream);		   \
    MemOpers<mode, float>::OP##_I(a_o.y, a_i.y, n, stream);		   \
    MemOpers<mode, float>::OP##_I(a_o.z, a_i.z, n, stream);		   \
}

#define FIELD_BINARY_ARRAY_OP_DEFS(OP)					   \
FIELD_BINARY_ARRAY_OPC_VEC_DEF(OP)					   \
FIELD_BINARY_ARRAY_OPC_FLOAT_DEF(OP)					   \
FIELD_BINARY_ARRAY_OPC_I_VEC_DEF(OP)					   \
FIELD_BINARY_ARRAY_OPC_I_FLOAT_DEF(OP)					   \
FIELD_BINARY_ARRAY_OP_FIELD_DEF(OP)					   \
FIELD_BINARY_ARRAY_OP_I_FIELD_DEF(OP)

namespace PyCA {

/**
 * ================ Non-Member Functions ================
 */
bool Is1D(const Field3D& d_i, size_t n) {
    return (d_i.isContinuous());
}

bool Is1D(const Field3D& d_i, const Field3D& d_i1, size_t n) {
    return Is1D(d_i, n) && Is1D(d_i1, n);
}

bool Is1D(const Field3D& d_i, const Field3D& d_i1, const Field3D& d_i2, size_t n) {
    return Is1D(d_i, n) && Is1D(d_i1, n) && Is1D(d_i2, n);
}

bool Is1D(const Field3D& d_i, const Field3D& d_i1, const Field3D& d_i2, const Field3D& d_i3, size_t n) {
    return Is1D(d_i, n) && Is1D(d_i1, n) && Is1D(d_i2, n) && Is1D(d_i3, n);
}

/**
 * ================ Member Functions ================
 */

/*
 * Set all the value of the Vector field with single value Vector3f(v)
 */
template<int mode>
void FieldOpers<mode>::
SetMem(Field3D& a_o, const Vec3Df& v,  StreamT stream, bool onDev){
    size_t n = a_o.nVox();
    MemOpers<mode, float>::SetMem(a_o.x, v.x, n, stream, onDev);
    MemOpers<mode, float>::SetMem(a_o.y, v.y, n, stream, onDev);
    MemOpers<mode, float>::SetMem(a_o.z, v.z, n, stream, onDev);
}

template<int mode>
void FieldOpers<mode>::
SetMem(Field3D& a_o, const Vec3Df& v, const Image3D& mask, StreamT stream, bool onDev){
    size_t n = a_o.nVox();
    MemOpers<mode, float>::SetMem(a_o.x, v.x, mask.get(), n, stream, onDev);
    MemOpers<mode, float>::SetMem(a_o.y, v.y, mask.get(), n, stream, onDev);
    MemOpers<mode, float>::SetMem(a_o.z, v.z, mask.get(), n, stream, onDev);
}
/**
 * Set all the value of the Vector field with single float value (normaly used to zero out the memory)
 */
template<int mode>
void FieldOpers<mode>::
SetMem(Field3D& a_o, const float& c,  StreamT stream, bool onDev){
    size_t n = a_o.nVox();
    if (Is1D(a_o, n)){
        MemOpers<mode, float>::SetMem(a_o.x, c, 3 * n, stream, onDev);
    } else {
        MemOpers<mode, float>::SetMem(a_o.x, c, n, stream, onDev);
        MemOpers<mode, float>::SetMem(a_o.y, c, n, stream, onDev);
        MemOpers<mode, float>::SetMem(a_o.z, c, n, stream, onDev);
    }
}

template<int mode>
void FieldOpers<mode>::
SetMem(Field3D& a_o, const float& c, const Image3D& mask, StreamT stream, bool onDev){
    size_t n = a_o.nVox();
    MemOpers<mode, float>::SetMem(a_o.x, c, mask.get(), n, stream, onDev);
    MemOpers<mode, float>::SetMem(a_o.y, c, mask.get(), n, stream, onDev);
    MemOpers<mode, float>::SetMem(a_o.z, c, mask.get(), n, stream, onDev);
}

FIELD_BINARY_ARRAY_OP_DEFS(Add)

FIELD_BINARY_ARRAY_OP_DEFS(Sub)

FIELD_BINARY_ARRAY_OP_DEFS(Mul)

FIELD_BINARY_ARRAY_OP_DEFS(Div)

FIELD_BINARY_ARRAY_OP_DEFS(Max)

FIELD_BINARY_ARRAY_OP_DEFS(Min)

/** @brief a_o = (a_i + a_i1) * c */
template<int mode>
void FieldOpers<mode>::
AddMulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const float& c,  StreamT stream, bool onDev){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    MK_CHECK3_MEM(a_o, a_i, a_i1);
    size_t n = a_o.nVox();
    if (Is1D(a_o, a_i, a_i1, n)){
        MemOpers<mode, float>::AddMulC(a_o.x, a_i.x, a_i1.x, c, 3 * n, stream, onDev);
        return;
    }
    MemOpers<mode, float>::AddMulC(a_o.x, a_i.x, a_i1.x, c, n, stream, onDev);
    MemOpers<mode, float>::AddMulC(a_o.y, a_i.y, a_i1.y, c, n, stream, onDev);
    MemOpers<mode, float>::AddMulC(a_o.z, a_i.z, a_i1.z, c, n, stream, onDev);
}

/** @brief a_o = (a_i - a_i1) * c */
template<int mode>
void FieldOpers<mode>::
SubMulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const float& c,  StreamT stream, bool onDev){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    MK_CHECK3_MEM(a_o, a_i, a_i1);
    size_t n = a_o.nVox();
    if (Is1D(a_o, a_i, a_i1, n)){
        MemOpers<mode, float>::SubMulC(a_o.x, a_i.x, a_i1.x, c, 3 * n, stream, onDev);
        return;
    }
    MemOpers<mode, float>::SubMulC(a_o.x, a_i.x, a_i1.x, c, n, stream, onDev);
    MemOpers<mode, float>::SubMulC(a_o.y, a_i.y, a_i1.y, c, n, stream, onDev);
    MemOpers<mode, float>::SubMulC(a_o.z, a_i.z, a_i1.z, c, n, stream, onDev);
}

/** @brief a_o = (a_i * a_i1) * c */
template<int mode>
void FieldOpers<mode>::
MulMulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const float& c,  StreamT stream, bool onDev){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    MK_CHECK3_MEM(a_o, a_i, a_i1);
    size_t n = a_o.nVox();
    if (Is1D(a_o, a_i, a_i1, n)){
        MemOpers<mode, float>::MulMulC(a_o.x, a_i.x, a_i1.x, c, 3 * n, stream, onDev);
        return;
    }
    MemOpers<mode, float>::MulMulC(a_o.x, a_i.x, a_i1.x, c, n, stream, onDev);
    MemOpers<mode, float>::MulMulC(a_o.y, a_i.y, a_i1.y, c, n, stream, onDev);
    MemOpers<mode, float>::MulMulC(a_o.z, a_i.z, a_i1.z, c, n, stream, onDev);
}

/** @brief a_o = (a_o + a_i) * c */
template<int mode>
void FieldOpers<mode>::
AddMulC_I(Field3D& a_o, const Field3D& a_i, const float& c,  StreamT stream, bool onDev){
    MK_CHECK2_SIZE(a_o, a_i);
    MK_CHECK2_MEM(a_o, a_i);
    size_t n = a_o.nVox();

    if (Is1D(a_o, a_i, n)){
        MemOpers<mode, float>::AddMulC_I(a_o.x, a_i.x, c, 3 * n, stream, onDev);
        return;
    }
    MemOpers<mode, float>::AddMulC_I(a_o.x, a_i.x, c, n, stream, onDev);
    MemOpers<mode, float>::AddMulC_I(a_o.y, a_i.y, c, n, stream, onDev);
    MemOpers<mode, float>::AddMulC_I(a_o.z, a_i.z, c, n, stream, onDev);
}

/** @brief a_o = (a_o -a_i) * c */
template<int mode>
void FieldOpers<mode>::
SubMulC_I(Field3D& a_o, const Field3D& a_i, const float& c,  StreamT stream, bool onDev){
    MK_CHECK2_SIZE(a_o, a_i);
    MK_CHECK2_MEM(a_o, a_i);
    size_t n = a_o.nVox();
    if (Is1D(a_o, a_i, n)){
        MemOpers<mode, float>::SubMulC_I(a_o.x, a_i.x, c, 3 * n, stream, onDev);
        return;
    }
    MemOpers<mode, float>::SubMulC_I(a_o.x, a_i.x, c, n, stream, onDev);
    MemOpers<mode, float>::SubMulC_I(a_o.y, a_i.y, c, n, stream, onDev);
    MemOpers<mode, float>::SubMulC_I(a_o.z, a_i.z, c, n, stream, onDev);
}

/** @brief a_o = (a_o * a_i) * c */
template<int mode>
void FieldOpers<mode>::
MulMulC_I(Field3D& a_o, const Field3D& a_i, const float& c,  StreamT stream, bool onDev){
    MK_CHECK2_SIZE(a_o, a_i);
    MK_CHECK2_MEM(a_o, a_i);
    size_t n = a_o.nVox();
    if (Is1D(a_o, a_i, n)){
        MemOpers<mode, float>::MulMulC_I(a_o.x, a_i.x, c, 3 * n, stream, onDev);
        return;
    }
    MemOpers<mode, float>::MulMulC_I(a_o.x, a_i.x, c, n, stream, onDev);
    MemOpers<mode, float>::MulMulC_I(a_o.y, a_i.y, c, n, stream, onDev);
    MemOpers<mode, float>::MulMulC_I(a_o.z, a_i.z, c, n, stream, onDev);
}



/** @brief d_o = a_i + a_i1 * c */
template<int mode>
void FieldOpers<mode>::
Add_MulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const float& c, StreamT stream, bool onDev)
{
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    MK_CHECK3_MEM(a_o, a_i, a_i1);
    size_t n = a_o.nVox();
    if (Is1D(a_o, a_i, a_i1, n)){
        MemOpers<mode, float>::Add_MulC(a_o.x, a_i.x, a_i1.x, c, 3 * n, stream, onDev);
        return;
    }
    MemOpers<mode, float>::Add_MulC(a_o.x, a_i.x, a_i1.x, c, n, stream, onDev);
    MemOpers<mode, float>::Add_MulC(a_o.y, a_i.y, a_i1.y, c, n, stream, onDev);
    MemOpers<mode, float>::Add_MulC(a_o.z, a_i.z, a_i1.z, c, n, stream, onDev);
}

/** @brief a_o = a_o + a_i * c */
template<int mode>
void FieldOpers<mode>::
Add_MulC_I(Field3D& a_o, const Field3D& a_i, const float& c, StreamT stream, bool onDev)
{
    MK_CHECK2_SIZE(a_o, a_i);
    MK_CHECK2_MEM(a_o, a_i);
    size_t n = a_o.nVox();
    if (Is1D(a_o, a_i, n)){
        MemOpers<mode, float>::Add_MulC_I(a_o.x, a_i.x, c, 3 * n, stream, onDev);
        return;
    }
    MemOpers<mode, float>::Add_MulC_I(a_o.x, a_i.x, c, n, stream, onDev);
    MemOpers<mode, float>::Add_MulC_I(a_o.y, a_i.y, c, n, stream, onDev);
    MemOpers<mode, float>::Add_MulC_I(a_o.z, a_i.z, c, n, stream, onDev);
}

/** @brief a_o = a_i * c + a_i1 */
template<int mode>
void FieldOpers<mode>::
MulCAdd(Field3D& a_o, const Field3D& a_i, const float& c, const Field3D& a_i1, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    MK_CHECK3_MEM(a_o, a_i, a_i1);
    size_t n = a_o.nVox();
    if (Is1D(a_o, a_i, a_i1, n)){
        MemOpers<mode, float>::MulCAdd(a_o.x, a_i.x, c, a_i1.x, 3 * n, stream, onDev);
        return;
    }
    MemOpers<mode, float>::MulCAdd(a_o.x, a_i.x, c, a_i1.x, n, stream, onDev);
    MemOpers<mode, float>::MulCAdd(a_o.y, a_i.y, c, a_i1.y, n, stream, onDev);
    MemOpers<mode, float>::MulCAdd(a_o.z, a_i.z, c, a_i1.z, n, stream, onDev);
}

/** @brief a_o = a_i * c - a_i1 */
template<int mode>
void FieldOpers<mode>::
MulCSub(Field3D& a_o, const Field3D& a_i, const float& c, const Field3D& a_i1, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    MK_CHECK3_MEM(a_o, a_i, a_i1);
    size_t n = a_o.nVox();
    if (Is1D(a_o, a_i, a_i1, n)){
        MemOpers<mode, float>::MulCSub(a_o.x, a_i.x, c, a_i1.x, 3 * n, stream, onDev);
        return;
    }
    MemOpers<mode, float>::MulCSub(a_o.x, a_i.x, c, a_i1.x, n, stream, onDev);
    MemOpers<mode, float>::MulCSub(a_o.y, a_i.y, c, a_i1.y, n, stream, onDev);
    MemOpers<mode, float>::MulCSub(a_o.z, a_i.z, c, a_i1.z, n, stream, onDev);
}


/** @brief a_o = a_o * c + a_i1 */
template<int mode>
void FieldOpers<mode>::
MulCAdd_I(Field3D& a_o, const float& c, const Field3D& a_i, StreamT stream, bool onDev){
    MK_CHECK2_SIZE(a_o, a_i);
    MK_CHECK2_MEM(a_o, a_i);
    size_t n = a_o.nVox();
    if (Is1D(a_o, a_i, n)){
        MemOpers<mode, float>::MulCAdd_I(a_o.x, c, a_i.x, 3 * n, stream, onDev);
        return;
    }
    MemOpers<mode, float>::MulCAdd_I(a_o.x, c, a_i.x, n, stream, onDev);
    MemOpers<mode, float>::MulCAdd_I(a_o.y, c, a_i.y, n, stream, onDev);
    MemOpers<mode, float>::MulCAdd_I(a_o.z, c, a_i.z, n, stream, onDev);
}

/** @brief a_o = a_o * c - a_i1 */
template<int mode>
void FieldOpers<mode>::
MulCSub_I(Field3D& a_o, const float& c, const Field3D& a_i, StreamT stream, bool onDev){
    MK_CHECK2_SIZE(a_o, a_i);
    MK_CHECK2_MEM(a_o, a_i);
    size_t n = a_o.nVox();
    if (Is1D(a_o, a_i, n)){
        MemOpers<mode, float>::MulCSub_I(a_o.x, c, a_i.x, 3 * n, stream, onDev);
        return;
    }
    MemOpers<mode, float>::MulCSub_I(a_o.x, c, a_i.x, n, stream, onDev);
    MemOpers<mode, float>::MulCSub_I(a_o.y, c, a_i.y, n, stream, onDev);
    MemOpers<mode, float>::MulCSub_I(a_o.z, c, a_i.z, n, stream, onDev);
}


/** @brief a_o = (a_i + a) * b */
template<int mode>
void FieldOpers<mode>::
AddCMulC(Field3D& a_o, const Field3D& a_i, const Vec3Df& a, const Vec3Df& b, StreamT stream, bool onDev){
    MK_CHECK2_SIZE(a_o, a_i);
    MK_CHECK2_MEM(a_o, a_i);
    size_t n = a_o.nVox();
    MemOpers<mode, float>::AddCMulC(a_o.x, a_i.x, a.x, b.x, n, stream, onDev);
    MemOpers<mode, float>::AddCMulC(a_o.y, a_i.y, a.y, b.y, n, stream, onDev);
    MemOpers<mode, float>::AddCMulC(a_o.z, a_i.z, a.z, b.z, n, stream, onDev);
}

/** @brief a_o = (a_i * a) - b */
template<int mode>
void FieldOpers<mode>::
MulCAddC(Field3D& a_o, const Field3D& a_i, const Vec3Df& a, const Vec3Df& b, StreamT stream, bool onDev){
    MK_CHECK2_SIZE(a_o, a_i);
    MK_CHECK2_MEM(a_o, a_i);
    size_t n = a_o.nVox();
    MemOpers<mode, float>::MulCAddC(a_o.x, a_i.x, a.x, b.x, n, stream, onDev);
    MemOpers<mode, float>::MulCAddC(a_o.y, a_i.y, a.y, b.y, n, stream, onDev);
    MemOpers<mode, float>::MulCAddC(a_o.z, a_i.z, a.z, b.z, n, stream, onDev);
}

/** @brief a_o = (a_o + a) * b */
template<int mode>
void FieldOpers<mode>::
AddCMulC_I(Field3D& a_o, const Vec3Df& a, const Vec3Df& b, StreamT stream, bool onDev){
    size_t n = a_o.nVox();
    MemOpers<mode, float>::AddCMulC_I(a_o.x, a.x, b.x, n, stream, onDev);
    MemOpers<mode, float>::AddCMulC_I(a_o.y, a.y, b.y, n, stream, onDev);
    MemOpers<mode, float>::AddCMulC_I(a_o.z, a.z, b.z, n, stream, onDev);
}

/** @brief a_o = (a_o * a) - b */
template<int mode>
void FieldOpers<mode>::
MulCAddC_I(Field3D& a_o, const Vec3Df& a, const Vec3Df& b, StreamT stream, bool onDev){
    size_t n = a_o.nVox();
    MemOpers<mode, float>::MulCAddC_I(a_o.x, a.x, b.x, n, stream, onDev);
    MemOpers<mode, float>::MulCAddC_I(a_o.y, a.y, b.y, n, stream, onDev);
    MemOpers<mode, float>::MulCAddC_I(a_o.z, a.z, b.z, n, stream, onDev);
}


/** @brief a_o = a_i * a + a_i1 * b */
template<int mode>
void FieldOpers<mode>::
MulC_Add_MulC(Field3D& a_o,const Field3D& a_i, const float& a,
              const Field3D& a_i1, const float& b, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    MK_CHECK3_MEM(a_o, a_i, a_i1);
    size_t n = a_o.nVox();
    if (Is1D(a_o, a_i, a_i1, n)){
        MemOpers<mode, float>::MulC_Add_MulC(a_o.x, a_i.x, a, a_i1.x, b, 3 * n, stream, onDev);
        return;
    }
    MemOpers<mode, float>::MulC_Add_MulC(a_o.x, a_i.x, a, a_i1.x, b, n, stream, onDev);
    MemOpers<mode, float>::MulC_Add_MulC(a_o.y, a_i.y, a, a_i1.y, b, n, stream, onDev);
    MemOpers<mode, float>::MulC_Add_MulC(a_o.z, a_i.z, a, a_i1.z, b, n, stream, onDev);
}

/** @brief a_o = a_o * a + a_i * b */
template<int mode>
void FieldOpers<mode>::
MulC_Add_MulC_I(Field3D& a_o, const float& a, const Field3D& a_i, const float& b, StreamT stream, bool onDev)
{
    MK_CHECK2_SIZE(a_o, a_i);
    MK_CHECK2_MEM(a_o, a_i);
    size_t n = a_o.nVox();
    if (Is1D(a_o, a_i, n)){
        MemOpers<mode, float>::MulC_Add_MulC_I(a_o.x, a, a_i.x, b, 3 * n, stream, onDev);
        return;
    }

    MemOpers<mode, float>::MulC_Add_MulC_I(a_o.x, a, a_i.x, b, n, stream, onDev);
    MemOpers<mode, float>::MulC_Add_MulC_I(a_o.y, a, a_i.y, b, n, stream, onDev);
    MemOpers<mode, float>::MulC_Add_MulC_I(a_o.z, a, a_i.z, b, n, stream, onDev);
}

/** @brief a_o = (a_i + a) * b + c */
template<int mode>
void FieldOpers<mode>::
AddCMulCAddC(Field3D& a_o, const Field3D& a_i, const float& a, const float& b, const float& c, StreamT stream, bool onDev){
    MK_CHECK2_SIZE(a_o, a_i);
    MK_CHECK2_MEM(a_o, a_i);
    size_t n = a_o.nVox();
    if (Is1D(a_o, a_i, n)){
        MemOpers<mode, float>::AddCMulCAddC(a_o.x, a_i.x, a, b, c, 3 * n, stream, onDev);
        return;
    }
    MemOpers<mode, float>::AddCMulCAddC(a_o.x, a_i.x, a, b, c, n, stream, onDev);
    MemOpers<mode, float>::AddCMulCAddC(a_o.y, a_i.y, a, b, c, n, stream, onDev);
    MemOpers<mode, float>::AddCMulCAddC(a_o.z, a_i.z, a, b, c, n, stream, onDev);
}

/** @brief a_o = (a_o + a) * b + c */
template<int mode>
void FieldOpers<mode>::
AddCMulCAddC_I(Field3D& a_o,  const float& a, const float& b, const float& c, StreamT stream, bool onDev){
    size_t n = a_o.nVox();
    if (Is1D(a_o, n)){
        MemOpers<mode, float>::AddCMulCAddC_I(a_o.x, a, b, c, n, stream, onDev);
        return;
    }
    MemOpers<mode, float>::AddCMulCAddC_I(a_o.x, a, b, c, n, stream, onDev);
    MemOpers<mode, float>::AddCMulCAddC_I(a_o.y, a, b, c, n, stream, onDev);
    MemOpers<mode, float>::AddCMulCAddC_I(a_o.z, a, b, c, n, stream, onDev);
}

/** @brief a_o = (a_i + a) * b + c */
template<int mode>
void FieldOpers<mode>::
AddCMulCAddC(Field3D& a_o, const Field3D& a_i, const Vec3Df& a, const Vec3Df& b, const Vec3Df& c, StreamT stream, bool onDev){
    MK_CHECK2_SIZE(a_o, a_i);
    MK_CHECK2_MEM(a_o, a_i);
    size_t n = a_o.nVox();
    MemOpers<mode, float>::AddCMulCAddC(a_o.x, a_i.x, a.x, b.x, c.x, n, stream, onDev);
    MemOpers<mode, float>::AddCMulCAddC(a_o.y, a_i.y, a.y, b.y, c.y, n, stream, onDev);
    MemOpers<mode, float>::AddCMulCAddC(a_o.z, a_i.z, a.z, b.z, c.z, n, stream, onDev);
}

/** @brief a_o = (a_o + a) * b + c */
template<int mode>
void FieldOpers<mode>::
AddCMulCAddC_I(Field3D& a_o,  const Vec3Df& a, const Vec3Df& b, const Vec3Df& c, StreamT stream, bool onDev){
    size_t n = a_o.nVox();
    MemOpers<mode, float>::AddCMulCAddC_I(a_o.x, a.x, b.x, c.x, n, stream, onDev);
    MemOpers<mode, float>::AddCMulCAddC_I(a_o.y, a.y, b.y, c.y, n, stream, onDev);
    MemOpers<mode, float>::AddCMulCAddC_I(a_o.z, a.z, b.z, c.z, n, stream, onDev);
}

/** @brief a_o = a_i + a_i1 * a_i2 * c */
template<int mode>
void FieldOpers<mode>::
Add_MulMulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1,
            const Field3D& a_i2, const float& c, StreamT stream, bool onDev)
{
    MK_CHECK4_SIZE(a_o, a_i, a_i1, a_i2);
    MK_CHECK4_MEM(a_o, a_i, a_i1, a_i2);
    size_t n = a_o.nVox();
    if (Is1D(a_o, a_i, a_i1, a_i2, n)){
        MemOpers<mode, float>::Add_MulMulC(a_o.x, a_i.x, a_i1.x, a_i2.x, c, 3 * n, stream, onDev);
        return;
    }
    MemOpers<mode, float>::Add_MulMulC(a_o.x, a_i.x, a_i1.x, a_i2.x, c, n, stream, onDev);
    MemOpers<mode, float>::Add_MulMulC(a_o.y, a_i.y, a_i1.y, a_i2.y, c, n, stream, onDev);
    MemOpers<mode, float>::Add_MulMulC(a_o.z, a_i.z, a_i1.z, a_i2.z, c, n, stream, onDev);
}
/** @brief a_o = a_o + a_i * a_i1 * d */
template<int mode>
void FieldOpers<mode>::
Add_MulMulC_I(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const float& c, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    MK_CHECK3_MEM(a_o, a_i, a_i1);
    size_t n = a_o.nVox();
    if (Is1D(a_o, a_i, a_i1, n)){
        MemOpers<mode, float>::Add_MulMulC_I(a_o.x, a_i.x, a_i1.x, c, 3 * n, stream, onDev);
        return;
    }
    MemOpers<mode, float>::Add_MulMulC_I(a_o.x, a_i.x, a_i1.x, c, n, stream, onDev);
    MemOpers<mode, float>::Add_MulMulC_I(a_o.y, a_i.y, a_i1.y, c, n, stream, onDev);
    MemOpers<mode, float>::Add_MulMulC_I(a_o.z, a_i.z, a_i1.z, c, n, stream, onDev);
}

/** @brief a_o = a_i - a_i1 * a_i2 * c */
template<int mode>
void FieldOpers<mode>::
Sub_MulMulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1,
            const Field3D& a_i2, const float& c, StreamT stream, bool onDev){
    MK_CHECK4_SIZE(a_o, a_i, a_i1, a_i2);
    MK_CHECK4_MEM(a_o, a_i, a_i1, a_i2);
    size_t n = a_o.nVox();
    if (Is1D(a_o, a_i, a_i1, a_i2, n)){
        MemOpers<mode, float>::Sub_MulMulC(a_o.x, a_i.x, a_i1.x, a_i2.x, c, 3 * n, stream, onDev);
        return;
    }
    MemOpers<mode, float>::Sub_MulMulC(a_o.x, a_i.x, a_i1.x, a_i2.x, c, n, stream, onDev);
    MemOpers<mode, float>::Sub_MulMulC(a_o.y, a_i.y, a_i1.y, a_i2.y, c, n, stream, onDev);
    MemOpers<mode, float>::Sub_MulMulC(a_o.z, a_i.z, a_i1.z, a_i2.z, c, n, stream, onDev);
}

/** @brief a_o = a_o - a_i * a_i1 * c*/
template<int mode>
void FieldOpers<mode>::
Sub_MulMulC_I(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const float& c, StreamT stream, bool onDev){
    MK_CHECK3_SIZE(a_o, a_i, a_i1);
    MK_CHECK3_MEM(a_o, a_i, a_i1);
    size_t n = a_o.nVox();
    if (Is1D(a_o, a_i, a_i1, n)){
        MemOpers<mode, float>::Sub_MulMulC_I(a_o.x, a_i.x, a_i1.x, c, 3 * n, stream, onDev);
        return;
    }
    MemOpers<mode, float>::Sub_MulMulC_I(a_o.x, a_i.x, a_i1.x, c, n, stream, onDev);
    MemOpers<mode, float>::Sub_MulMulC_I(a_o.y, a_i.y, a_i1.y, c, n, stream, onDev);
    MemOpers<mode, float>::Sub_MulMulC_I(a_o.z, a_i.z, a_i1.z, c, n, stream, onDev);
}

////////////////////////////////////////////////////////////////////////////
// compose a velocity and hfield to get an hfield
// h(x) = g(x) + delta * v(g(x))
////////////////////////////////////////////////////////////////////////////
template<int mode>
template<BackgroundStrategy bg>
void FieldOpers<mode>::
ComposeHH(Field3D& f, const Field3D& g, const Field3D& h,
          StreamT stream)
{
    Executer::template ComposeHH<bg>(f, g, h, stream);
}

////////////////////////////////////////////////////////////////////////////
// compose a velocity and hfield to get an hfield
// h(x) = g(x) + delta * v(g(x))
////////////////////////////////////////////////////////////////////////////
template<int mode>
template<BackgroundStrategy bg>
void FieldOpers<mode>::
ComposeVH(Field3D& h, const Field3D& v, const Field3D& g, const float& delta,
          StreamT stream, bool onDev) {
    Executer::template ComposeVH<bg>(h, v, g, delta, stream, onDev);
}

////////////////////////////////////////////////////////////////////////////
// compose an inverse velocity and hfield to get an hfield
// h(x) = g(x) - delta * v(g(x))
////////////////////////////////////////////////////////////////////////////
template<int mode>
template<BackgroundStrategy bg>
void FieldOpers<mode>::
ComposeVInvH(Field3D& h, const Field3D& v, const Field3D& g, const float& delta, StreamT stream, bool onDev) {
    Executer::template ComposeVInvH<bg>(h, v, g, delta, stream, onDev);
}

/**
 * compose a h field and a velocify field to get an hfield
 * h(x) = g(x+ delta * v(x))
 *
 * davisb 2007
 */
template<int mode>
template<BackgroundStrategy bg>
void FieldOpers<mode>::
ComposeHV(Field3D& h, const Field3D& g, const Field3D& v, const float& delta, StreamT stream, bool onDev){
    Executer::template ComposeHV<bg>(h, g, v, delta, stream, onDev);
}

template<int mode>
template<BackgroundStrategy bg>
void FieldOpers<mode>::
ComposeHVInv(Field3D& h, const Field3D& g, const Field3D& v, const float& delta, StreamT stream, bool onDev){
    Executer::template ComposeHVInv<bg>(h, g, v, delta, stream, onDev);
}


////////////////////////////////////////////////////////////////////////////
// compose field with translation
// creating a_o(x) = a_i(x + t)
////////////////////////////////////////////////////////////////////////////
template<int mode>
template<BackgroundStrategy bg>
void FieldOpers<mode>::
ComposeTranslation(Field3D& a_o, const Field3D& a_i, const Vec3Df& t, StreamT stream, bool onDev){
    Executer::template ComposeTranslation<bg>(a_o, a_i, t, stream, onDev);
}

template<int mode>
template<BackgroundStrategy bg>
void FieldOpers<mode>::
ApplyH(Field3D& a_o, const Field3D& a_i, const Field3D& h, StreamT stream) {
    Executer::template ApplyH<bg>(a_o, a_i, h, stream);
}

template<int mode>
template<BackgroundStrategy bg>
void FieldOpers<mode>::
ApplyV(Field3D& a_o, const Field3D& a_i, const Field3D& u, const float& delta, StreamT stream, bool onDev) {
    Executer::template ApplyV<bg>(a_o, a_i, u, delta, stream, onDev);
}

template<int mode>
template<BackgroundStrategy bg>
void FieldOpers<mode>::
ApplyVInv(Field3D& a_o, const Field3D& a_i, const Field3D& u, const float& delta, StreamT stream, bool onDev) {
    Executer::template ApplyVInv<bg>(a_o, a_i, u, delta, stream, onDev);
}

template<int mode>
void FieldOpers<mode>::
Jacobian(Field3D& d_Xg, Field3D& d_Yg, Field3D& d_Zg,
         const Field3D& d_h,
	 DiffT diffType, BoundaryCondT bc,
	 StreamT stream)
{
   ImageFieldOpers<mode>::Gradient(d_Xg, d_h.x, diffType, bc, stream);
   ImageFieldOpers<mode>::Gradient(d_Yg, d_h.y, diffType, bc, stream);
   ImageFieldOpers<mode>::Gradient(d_Zg, d_h.z, diffType, bc, stream);
}

template<int mode>
void FieldOpers<mode>::
SplatField(Field3D& a_o, const Field3D& a_i,
		  const Field3D& h, StreamT stream) {
   Executer::SplatField(a_o, a_i, h, stream);
}

template<int mode>
void FieldOpers<mode>::
SubVol(Field3D& a_o, const Field3D& a_i,
       const Vec3Di& start, StreamT stream)
{
   Executer::SubVol(a_o, a_i, start, stream);
}

template<int mode>
void FieldOpers<mode>::
SetSubVol_I(Field3D& a_o, const Field3D& a_i,
	    const Vec3Di& start, StreamT stream)
{
   Executer::SetSubVol_I(a_o, a_i, start, stream);
}

template<int mode>
template<BackgroundStrategy bg, bool rescaleVector>
void FieldOpers<mode>::
Resample(Field3D& a_o, const Field3D& a_i, StreamT stream){
    Executer::template Resample<bg, rescaleVector>(a_o, a_i, stream);
}

template<int mode>
void FieldOpers<mode>::
ReprojectToUnitVec(Field3D& a_o, StreamT stream){
    Executer::ReprojectToUnitVec(a_o, stream);
}

template<int mode>
void FieldOpers<mode>::
NormalizeSafe(Field3D& a_o, const Field3D& a_i, const float& eps,
	      StreamT stream)
{
   Executer::NormalizeSafe(a_o, a_i, eps, stream);
}

template<int mode>
void FieldOpers<mode>::
NormalizeSafe_I(Field3D& a_o, const float& eps, StreamT stream)
{
   Executer::NormalizeSafe_I(a_o, eps, stream);
}

template<int mode>
void FieldOpers<mode>::
Shrink(Field3D& a_o, const Field3D& a_i, const float& eps,
	      StreamT stream)
{
   Executer::Shrink(a_o, a_i, eps, stream);
}

template<int mode>
void FieldOpers<mode>::
Shrink_I(Field3D& a_o, const float& eps, StreamT stream)
{
   Executer::Shrink_I(a_o, eps, stream);
}

// Given a deformation g, a displacement vector field v, and an inverse ginv,
// of v, invert the deformation gt(x) = g(x) + dt v(g(x))
template<int mode>
template<BackgroundStrategy bg>
void FieldOpers<mode>::
FixedPointInverse(Field3D &ginv, const Field3D& g, unsigned int numIter, StreamT stream, bool onDev)
{
    throw PyCAException(__FILE__, __LINE__, "Inverse not yet implemented/tested");
    Executer::template FixedPointInverse<bg>(ginv, g, numIter, stream, onDev);
}

// Given a forward update on deformation g as w, and its inverse g^{-1},
// this function computes updated g^{-1} corresponding this forward update on g as per
// g_{t+1,0} = g_{t,0}\circ(Id - w\circ g_{t+1,0})
template<int mode>
void FieldOpers<mode>::
UpdateInverse(Field3D &ginv0t1,Field3D &scratchV, const Field3D& ginv0t, const Field3D& w, unsigned int numIter, StreamT stream, bool onDev)
{
  Executer::UpdateInverse(ginv0t1, scratchV, ginv0t, w, numIter, stream, onDev);
}
////////////////////////////////////////////////////////////////////////////
// Lie algebra operations
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
// adjoint action of Diff on its Lie algebra
// Z = Ad_g X = |Dg|\circ g^{-1} X\circ g^{-1}
////////////////////////////////////////////////////////////////////////////
template<int mode>
template<BackgroundStrategy bg>
void FieldOpers<mode>::
Ad(Field3D& Z, const Field3D& g, const Field3D& X,
          StreamT stream, bool onDev) {
    Executer::template Ad<bg>(Z, g, X, stream, onDev);
}
////////////////////////////////////////////////////////////////////////////
// infinitesimal adjoint action, equal to negative Jacobi-Lie bracket
// Z = ad_X Y = DX Y - DY X
////////////////////////////////////////////////////////////////////////////
template<int mode>
void FieldOpers<mode>::
AdInf(Field3D& Z, const Field3D& X, const Field3D& Y,
          StreamT stream, bool onDev) {
    Executer::AdInf(Z, X, Y, stream, onDev);
}

////////////////////////////////////////////////////////////////////////////
// coadjoint action of Diff on its Lie coalgebra
// n = Ad_g^* m = (Dg)^T m\circ g |Dg|
////////////////////////////////////////////////////////////////////////////
template<int mode>
template<BackgroundStrategy bg>
void FieldOpers<mode>::
CoAd(Field3D& n, const Field3D& g, const Field3D& m,
          StreamT stream, bool onDev) {
    Executer::template CoAd<bg>(n, g, m, stream, onDev);
}
////////////////////////////////////////////////////////////////////////////
// infinitesimal coadjoint action
// n = ad_X^* m = (DX)^T m + div(m \otimes X)
////////////////////////////////////////////////////////////////////////////
template<int mode>
void FieldOpers<mode>::
CoAdInf(Field3D& n, const Field3D& X, const Field3D& m,
          StreamT stream, bool onDev) {
    Executer::CoAdInf(n, X, m, stream, onDev);
}

////////////////////////////////////////////////////////////////////////////
// computes tensor divergence of outer product of two vector fields 
// Z = div(X \otimes Y)
////////////////////////////////////////////////////////////////////////////
template<int mode>
void FieldOpers<mode>::
DivergenceTensor(Field3D& Z, const Field3D& X, const Field3D& Y,
	StreamT stream, bool onDev){
  Executer::DivergenceTensor(Z, X, Y, stream, onDev);
}


////////////////////////////////////////////////////////////////////////////
// computes jacobian of X times Y
// Z = DX Y
////////////////////////////////////////////////////////////////////////////
template<int mode>
void FieldOpers<mode>::
JacobianXY(Field3D& Z, const Field3D& X, const Field3D& Y,
          StreamT stream, bool onDev) {
    Executer::JacobianXY(Z, X, Y, stream, onDev);
}

////////////////////////////////////////////////////////////////////////////
// computes jacobian of X transpose times Y
// Z = (DX)' Y
////////////////////////////////////////////////////////////////////////////
template<int mode>
void FieldOpers<mode>::
JacobianXtY(Field3D& Z, const Field3D& X, const Field3D& Y,
	    StreamT stream, bool onDev) {
    Executer::JacobianXtY(Z, X, Y, stream, onDev);
}

// TODO: Figure out how best to pass in a kernel (DiffOper) for the following -JDH

//////////////////////////////////////////////////////////////////////////////
//// adjoint transpose action, given metric
//// Z = Ad_g^T X = K Ad_g^* LX
//////////////////////////////////////////////////////////////////////////////
//template<int mode>
//template<BackgroundStrategy bg>
//void FieldOpers<mode>::
//AdTranspose(Field3D& Z, const Field3D& g, const Field3D& X,
          //StreamT stream, bool onDev) {
    //Executer::template AdTranspose<bg>(Z, g, X, stream, onDev);
//}
//////////////////////////////////////////////////////////////////////////////
//// transposed infinitesimal adjoint action
//// Z = ad_X^T Y = K ad_X^* LY
//////////////////////////////////////////////////////////////////////////////
//template<int mode>
//template<BackgroundStrategy bg>
//void FieldOpers<mode>::
//adtranspose(Field3D& Z, const Field3D& X, const Field3D& Y,
          //StreamT stream, bool onDev) {
    //Executer::template adtranspose<bg>(Z, X, Y, stream, onDev);
//}

//////////////////////////////////////////////////////////////////////////////
//// symmetric product, as defined in Bullo1995
//// Z = sym_X Y = -(ad_X^T Y + ad_Y^T X)
//////////////////////////////////////////////////////////////////////////////
//template<int mode>
//template<BackgroundStrategy bg>
//void FieldOpers<mode>::
//sym(Field3D& Z, const Field3D& X, const Field3D& Y,
          //StreamT stream, bool onDev) {
    //Executer::template sym<bg>(Z, X, Y, stream, onDev);
//}

// template instantiation
#include "FieldOpers_inst.cxx"

} // end namespace PyCA

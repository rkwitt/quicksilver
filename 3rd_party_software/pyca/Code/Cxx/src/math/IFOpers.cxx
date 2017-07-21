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

#include <IFOpers.h>

#include "ImageFieldOpers.h"

#include <algorithm>

#include <Image3D.h>
#include <Field3D.h>
#include <MemOpers.h>

namespace PyCA {

void Opers
::Copy(Image3D& a_o, const Field3D& a_i, int dim, StreamT stream){
    MK_CHECK2_SIZE(a_o, a_i);
    switch(dim){
    case 0:
        Copy(a_o.getMemPool(), a_i.getMemPoolX(), a_o.nVox(), stream);
        break;
    case 1:
        Copy(a_o.getMemPool(), a_i.getMemPoolY(), a_o.nVox(), stream);
        break;
    case 2:
        Copy(a_o.getMemPool(), a_i.getMemPoolZ(), a_o.nVox(), stream);
        break;
    default:
       throw PyCAException(__FILE__, __LINE__,"dim should be in range [0,2]");
    }
}

void Opers
::Copy(Field3D& a_o, const Image3D& a_i, int dim, StreamT stream){
    MK_CHECK2_SIZE(a_o, a_i);
    switch(dim){
    case 0:
        Copy(a_o.getMemPoolX(), a_i.getMemPool(), a_o.nVox(), stream);
        break;
    case 1:
        Copy(a_o.getMemPoolY(), a_i.getMemPool(), a_o.nVox(), stream);
        break;
    case 2:
        Copy(a_o.getMemPoolZ(), a_i.getMemPool(), a_o.nVox(), stream);
        break;
    default:
       throw PyCAException(__FILE__, __LINE__,"dim should be in range [0,2]");
    }
}

/////////////////////////////////////////////////////////////////////////////
// apply hField to an image
// defImage(x) = image(h(x))
/////////////////////////////////////////////////////////////////////////////
template<BackgroundStrategy bg, InterpT interp>
void Opers::
ApplyH(Image3D& a_o, const Image3D& a_i, const Field3D& h, StreamT s){
   MK_CHECK3_ALL(a_o, a_i, h);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     ApplyH<bg COMMA interp>(a_o, a_i, h, s));
}

template<BackgroundStrategy bg>
void
applyH(Image3D& a_o, const Image3D& a_i, const Field3D& a_h, 
       InterpT interp, StreamT s)
{
   if(interp == INTERP_NN){
       Opers::ApplyH<bg, INTERP_NN>(a_o, a_i, a_h, s);
   }else if(interp == INTERP_LINEAR){
       Opers::ApplyH<bg, INTERP_LINEAR>(a_o, a_i, a_h, s);
   }else if(interp == INTERP_CUBIC){
       Opers::ApplyH<bg, INTERP_CUBIC>(a_o, a_i, a_h, s);
   }else{
      // throw exception
      throw PyCAException(__FILE__, __LINE__, "Unknown/unimplemented interp type");
   }
}

void Opers::
ApplyH(Image3D& a_o, const Image3D& a_i, const Field3D& a_h, 
       BackgroundStrategy bg, InterpT interp, 
       StreamT s)
{
   if(bg == BACKGROUND_STRATEGY_CLAMP){
       applyH<BACKGROUND_STRATEGY_CLAMP>(a_o, a_i, a_h, interp, s);
   }else if(bg == BACKGROUND_STRATEGY_WRAP){
       applyH<BACKGROUND_STRATEGY_WRAP>(a_o, a_i, a_h, interp, s);
   }else if(bg == BACKGROUND_STRATEGY_ZERO){
       applyH<BACKGROUND_STRATEGY_ZERO>(a_o, a_i, a_h, interp, s);
   }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ZERO){
       applyH<BACKGROUND_STRATEGY_PARTIAL_ZERO>(a_o, a_i, a_h, interp, s);
   }else{
      // throw exception
      throw PyCAException(__FILE__, __LINE__, "Unknown/unimplemented background strategy");
   }
}

/////////////////////////////////////////////////////////////////////////////
// apply uField to an image
// defImage(x) = image(x + delta * u(x))
/////////////////////////////////////////////////////////////////////////////

template<BackgroundStrategy bg, InterpT interp>
void Opers::
ApplyV(Image3D& a_o, const Image3D& a_i, const Field3D& u,const float& delta, StreamT s, bool onDev)
{
   MK_CHECK3_ALL(a_o, a_i, u);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     ApplyV<bg COMMA interp>(a_o, a_i, u, delta, s, onDev));
}

template<BackgroundStrategy bg>
void
applyV(Image3D& a_o, const Image3D& a_i, const Field3D& u, 
       const float& delta, InterpT interp, StreamT s, bool onDev)
{
   if(interp == INTERP_NN){
       Opers::ApplyV<bg, INTERP_NN>
	   (a_o, a_i, u, delta, s, onDev);
   }else if(interp == INTERP_LINEAR){
       Opers::ApplyV<bg, INTERP_LINEAR>
	   (a_o, a_i, u, delta, s, onDev);
   }else if(interp == INTERP_CUBIC){
       Opers::ApplyV<bg, INTERP_CUBIC>
	   (a_o, a_i, u, delta, s, onDev);
   }else{
       // throw exception
       throw PyCAException(__FILE__, __LINE__, "Unknown/unimplemented interp type");
   }
}

void Opers::
ApplyV(Image3D& a_o, const Image3D& a_i, const Field3D& a_u, const float& delta, 
       BackgroundStrategy bg, InterpT interp, 
       StreamT s, bool onDev)
{
   if(bg == BACKGROUND_STRATEGY_CLAMP){
       applyV<BACKGROUND_STRATEGY_CLAMP>
	   (a_o, a_i, a_u, delta, interp, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_WRAP){
      applyV<BACKGROUND_STRATEGY_WRAP>
	   (a_o, a_i, a_u, delta, interp, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_ZERO){
      applyV<BACKGROUND_STRATEGY_ZERO>
	   (a_o, a_i, a_u, delta, interp, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ZERO){
      applyV<BACKGROUND_STRATEGY_PARTIAL_ZERO>
	   (a_o, a_i, a_u, delta, interp, s, onDev);
   }else{
      // throw exception
      throw PyCAException(__FILE__, __LINE__, "Unknown/unimplemented background strategy");
   }
}

/////////////////////////////////////////////////////////////////////////////
// apply uField to an image
// defImage(x) = image(x - delta * u(x))
/////////////////////////////////////////////////////////////////////////////

template<BackgroundStrategy bg, InterpT interp>
void Opers::
ApplyVInv(Image3D& a_o, const Image3D& a_i, const Field3D& u,const float& delta, StreamT s, bool onDev)
{
   MK_CHECK3_ALL(a_o, a_i, u);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     ApplyVInv<bg COMMA interp>(a_o, a_i, u, delta, s, onDev));
};

template<BackgroundStrategy bg>
void
applyVInv(Image3D& a_o, const Image3D& a_i, const Field3D& u, 
	  const float& delta, InterpT interp, StreamT s, bool onDev)
{
   if(interp == INTERP_NN){
       Opers::ApplyVInv<bg, INTERP_NN>
	   (a_o, a_i, u, delta, s, onDev);
   }else if(interp == INTERP_LINEAR){
       Opers::ApplyVInv<bg, INTERP_LINEAR>
	   (a_o, a_i, u, delta, s, onDev);
   }else if(interp == INTERP_CUBIC){
       Opers::ApplyVInv<bg, INTERP_CUBIC>
	   (a_o, a_i, u, delta, s, onDev);
   }else{
       // throw exception
       throw PyCAException(__FILE__, __LINE__, "Unknown/unimplemented interp type");
   }
}

void Opers::
ApplyVInv(Image3D& a_o, const Image3D& a_i, const Field3D& a_u, const float& delta, 
	  BackgroundStrategy bg, InterpT interp, 
	  StreamT s, bool onDev)
{
   if(bg == BACKGROUND_STRATEGY_CLAMP){
       applyVInv<BACKGROUND_STRATEGY_CLAMP>
	   (a_o, a_i, a_u, delta, interp, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_WRAP){
      applyVInv<BACKGROUND_STRATEGY_WRAP>
	  (a_o, a_i, a_u, delta, interp, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_ZERO){
      applyVInv<BACKGROUND_STRATEGY_ZERO>
	  (a_o, a_i, a_u, delta, interp, s, onDev);
   }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ZERO){
      applyVInv<BACKGROUND_STRATEGY_PARTIAL_ZERO>
	  (a_o, a_i, a_u, delta, interp, s, onDev);
   }else{
      // throw exception
      throw PyCAException(__FILE__, __LINE__, "Unknown/unimplemented background strategy");
   }
}

template<BackgroundStrategy bg, InterpT interp>
void Opers::
ComposeTranslation(Image3D& a_o, const Image3D& a_i, const Vec3Df& t, StreamT s, bool onDev){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     ComposeTranslation<bg COMMA interp>(a_o, a_i, t, s, onDev));
}

template<BackgroundStrategy bg>
void
composeTranslation(Image3D& a_o, const Image3D& a_i, const Vec3Df& t, 
		   InterpT interp, 
		   StreamT s, bool onDev)
{
   if(interp == INTERP_NN){
       Opers::ComposeTranslation<bg, INTERP_NN>
	   (a_o, a_i, t, s, onDev);
   }else if(interp == INTERP_LINEAR){
       Opers::ComposeTranslation<bg, INTERP_LINEAR>
	   (a_o, a_i, t, s, onDev);
   }else if(interp == INTERP_CUBIC){
       Opers::ComposeTranslation<bg, INTERP_CUBIC>
	   (a_o, a_i, t, s, onDev);
   }else{
       // throw exception
       throw PyCAException(__FILE__, __LINE__, "Unknown/unimplemented interp type");
   }
}

void Opers::
ComposeTranslation(Image3D& a_o, const Image3D& a_i, const Vec3Df& t, 
		   BackgroundStrategy bg, InterpT interp, 
		   StreamT s, bool onDev)
{
    if(bg == BACKGROUND_STRATEGY_CLAMP){
	composeTranslation<BACKGROUND_STRATEGY_CLAMP>
	    (a_o, a_i, t, interp, s, onDev);
    }else if(bg == BACKGROUND_STRATEGY_WRAP){
	composeTranslation<BACKGROUND_STRATEGY_WRAP>
	    (a_o, a_i, t, interp, s, onDev);
    }else if(bg == BACKGROUND_STRATEGY_ZERO){
	composeTranslation<BACKGROUND_STRATEGY_ZERO>
	    (a_o, a_i, t, interp, s, onDev);
    }else if(bg == BACKGROUND_STRATEGY_PARTIAL_ZERO){
	composeTranslation<BACKGROUND_STRATEGY_PARTIAL_ZERO>
	    (a_o, a_i, t, interp, s, onDev);
    }else{
	// throw exception
	throw PyCAException(__FILE__, __LINE__, "Unknown/unimplemented background strategy");
    }
}


void Opers
::Splat(Image3D& a_o, const Field3D& a_h, const Image3D& a_i,
               bool normalize, StreamT stream)
{
   MK_CHECK3_ALL(a_o, a_h, a_i);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     Splat(a_o, a_h, a_i, normalize, stream));
}


void Opers::
FiniteDiff(Image3D& a_o, const Image3D& a_i, 
	   DimT dim, DiffT diffType, 
	   enum BoundaryCondT bc, 
	   bool accum, OpT op,
	   StreamT stream)
{
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     FiniteDiff(a_o, a_i, dim, diffType, bc, accum, op, stream));
}

void Opers::
Gradient(Field3D& a_o, const float* a_i, DiffT diffType, BoundaryCondT bc, StreamT stream)
{
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     Gradient(a_o, a_i, diffType, bc, stream));
}

void Opers::
Gradient(Field3D& a_o, const Image3D& a_i, DiffT diffType, BoundaryCondT bc, StreamT stream){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     Gradient(a_o, a_i, diffType, bc, stream));
}

void Opers::
Gradient2(Field3D& a_o, const Image3D& a_i, DiffT diffType, BoundaryCondT bc, StreamT stream){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     Gradient2(a_o, a_i, diffType, bc, stream));
}

void Opers::
GradientMask(Field3D& a_o, 
	     const Image3D& a_i, const Image3D& a_mask, 
	     DiffT diffType, BoundaryCondT bc, 
	     StreamT stream)
{
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     GradientMask(a_o, a_i, a_mask, 
			  diffType, bc, stream));
}

void Opers::
GradFor(Field3D& a_o, const Image3D& a_i, BoundaryCondT bc, StreamT stream){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     GradFor(a_o, a_i, bc, stream));
}

void Opers::
GradForMag(Image3D& a_o, const Image3D& a_i, StreamT stream){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     GradForMag(a_o, a_i, stream));
}

void Opers::
GradientMag(Image3D& a_o, const Image3D& a_i, DiffT diffType, BoundaryCondT bc, StreamT stream){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     GradientMag(a_o, a_i, diffType, bc, stream));
}


void Opers::
UpwindDiff(Image3D& a_o, const Image3D& a_i, const Image3D& a_speed, 
	   DimT dim, StreamT stream)
{
   MK_CHECK3_ALL(a_o, a_i, a_speed);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     UpwindDiff(a_o, a_i, a_speed, dim, stream));
}

void Opers::
UpwindGradMag(Image3D& a_o, const Image3D& a_i,
	      const Image3D& a_speed, StreamT stream)
{
   MK_CHECK3_ALL(a_o, a_i, a_speed);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     UpwindGradMag(a_o, a_i, a_speed, stream));
}

void Opers::
Divergence(Image3D& a_o, const Field3D& a_i, DiffT diffType, BoundaryCondT bc, StreamT stream){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     Divergence(a_o, a_i, diffType, bc, stream));
}

void Opers::
DivBack(Image3D& a_o, const Field3D& a_i, BoundaryCondT bc, StreamT stream){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     DivBack(a_o, a_i, bc, stream));
}

/**
 * Compute the magnitude image
 * a_o[i] = sqrt(a_i[i].x^2 + a_i[i].y^2 + a_i[i].z^2) 
 */

void Opers::
Magnitude(Image3D& a_o, const Field3D& a_i, StreamT stream){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     Magnitude(a_o, a_i, stream));
}

/**
 * Compute the magnitude array
 * a_o[i] = a_i[i].x^2 + a_i[i].y^2 + a_i[i].z^2 
 */
void Opers::
SqrMagnitude(Image3D& a_o, const Field3D& a_i, StreamT stream){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     SqrMagnitude(a_o, a_i, stream));
}

/**
 * Compute dot product image
 * a_o[i] = a_i[i].x * a_i1[i].x + a_i[i].y * a_i1[i].y + a_i[i].z * a_i1[i].z
 */


void Opers::
ComponentDotProd(Image3D& a_o, const Field3D& a_i, const Field3D& a_i1, StreamT stream) {
   MK_CHECK3_ALL(a_o, a_i, a_i1);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     ComponentDotProd(a_o, a_i, a_i1, stream));
}

/** @brief a_o.x = a_i.x + a_i1, a_o.y = a_i.y + a_i1, a_o.z = a_i.z + a_i1 */


void Opers::
Add(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, StreamT stream){
   MK_CHECK3_ALL(a_o, a_i, a_i1);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     Add(a_o, a_i, a_i1, stream));
}

/** @brief a_o.x = a_i.x - a_i1, a_o.y = a_i.y - a_i1, a_o.z = a_i.z - a_i1 */


void Opers::
Sub(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, StreamT stream){
   MK_CHECK3_ALL(a_o, a_i, a_i1);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     Sub(a_o, a_i, a_i1, stream));
}
/** @brief a_o.x = a_i.x * a_i1, a_o.y = a_i.y * a_i1, a_o.z = a_i.z * a_i1 */


void Opers::
Mul(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, StreamT stream){
   MK_CHECK3_ALL(a_o, a_i, a_i1);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     Mul(a_o, a_i, a_i1, stream));
}

/** @brief a_o.x = a_i.x / a_i1, a_o.y = a_i.y / a_i1, a_o.z = a_i.z / a_i1 */

void Opers::
Div(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, StreamT stream){
   MK_CHECK3_ALL(a_o, a_i, a_i1);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     Div(a_o, a_i, a_i1, stream));
}

/** @brief a_o.x += a_i, a_o.y += a_i, a_o.y += a_i, */


void Opers::
Add_I(Field3D& a_o, const Image3D& a_i, StreamT stream){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     Add_I(a_o, a_i, stream));
}
/** @brief a_o.x -= a_i, a_o.y -= a_i, a_o.y -= a_i, */


void Opers::
Sub_I(Field3D& a_o, const Image3D& a_i, StreamT stream){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     Sub_I(a_o, a_i, stream));
}
/** @brief a_o.x *= a_i, a_o.y *= a_i, a_o.y *= a_i, */


void Opers::
Mul_I(Field3D& a_o, const Image3D& a_i, StreamT stream){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     Mul_I(a_o, a_i, stream));
}

/** @brief a_o.x /= a_i, a_o.y /= a_i, a_o.y /= a_i*/


void Opers::
Div_I(Field3D& a_o, const Image3D& a_i, StreamT stream){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     Div_I(a_o, a_i, stream));
}

/** @brief a_o = a_i + a_i1 * a_i2 */


void Opers::
Add_Mul(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const Image3D& a_i2, StreamT stream){
   MK_CHECK3_ALL(a_o, a_i, a_i1);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     Add_Mul(a_o, a_i, a_i1, a_i2, stream));
}
/** @brief a_o = a_o + a_i * a_i1 */


void Opers::
Add_Mul_I(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, StreamT stream){
   MK_CHECK3_ALL(a_o, a_i, a_i1);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     Add_Mul_I(a_o, a_i, a_i1, stream));
}

/** @brief a_o = a_i - a_i1 * a_i2 */


void Opers::
Sub_Mul(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const Image3D& a_i2, StreamT stream){
   MK_CHECK3_ALL(a_o, a_i, a_i1);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     Sub_Mul(a_o, a_i, a_i1, a_i2, stream));
}
/** @brief a_o = a_o - a_i * a_i1 */


void Opers::
Sub_Mul_I(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, StreamT stream){
   MK_CHECK3_ALL(a_o, a_i, a_i1);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     Sub_Mul_I(a_o, a_i, a_i1, stream));
}

/** @brief a_o = a_i * a_i1 * c (a_o.x = a_i.x * a_i1 * c)*/


void Opers::
MulMulC(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, const float& c, StreamT stream, bool onDev){
   MK_CHECK3_ALL(a_o, a_i, a_i1);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     MulMulC(a_o, a_i, a_i1, c, stream, onDev));
}
/** @brief a_o = a_o * a_i * c  (a_o.x = a_o.x * a_i * c)*/


void Opers::
MulMulC_I(Field3D& a_o, const Image3D& a_i, const float& c, StreamT stream, bool onDev){
   MK_CHECK2_ALL(a_o, a_i);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     MulMulC_I(a_o, a_i, c, stream, onDev));
}

/** @brief a_o = a_i + a_i1 * a_i2 * c */


void Opers::
Add_MulMulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const Image3D& a_i2, const float& c,  StreamT stream, bool onDev){
   MK_CHECK4_ALL(a_o, a_i, a_i1, a_i2);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     Add_MulMulC(a_o, a_i, a_i1, a_i2, c, stream, onDev));
}

/** @brief a_o = a_o + a_i * a_i1 * c */


void Opers::
Add_MulMulC_I(const Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, const float& c,  StreamT stream, bool onDev){
   MK_CHECK3_ALL(a_o, a_i, a_i1);
   AUTO_EXEC(a_o.memType(), ImageFieldOpers, 
	     Add_MulMulC_I(a_o, a_i, a_i1, c, stream, onDev));
}



void Opers::
JacDetH(Image3D& a_detJ,
        const Field3D& a_Xg, const Field3D& a_Yg, const Field3D& a_Zg,
        StreamT stream)
{
   MK_CHECK4_ALL(a_detJ, a_Xg, a_Yg, a_Zg);
   AUTO_EXEC(a_detJ.memType(), ImageFieldOpers, 
	     JacDetH(a_detJ, a_Xg, a_Yg, a_Zg, stream));
}


void Opers::
JacDetV(Image3D& a_detJ,
        const Field3D& a_Xg, const Field3D& a_Yg, const Field3D& a_Zg,
        StreamT stream)
{
   MK_CHECK4_ALL(a_detJ, a_Xg, a_Yg, a_Zg);
   AUTO_EXEC(a_detJ.memType(), ImageFieldOpers, 
	     JacDetV(a_detJ, a_Xg, a_Yg, a_Zg, stream));
}


void Opers::
JacDetVInv(Image3D& a_detJ,
        const Field3D& a_Xg, const Field3D& a_Yg, const Field3D& a_Zg,
        StreamT stream)
{
   MK_CHECK4_ALL(a_detJ, a_Xg, a_Yg, a_Zg);
   AUTO_EXEC(a_detJ.memType(), ImageFieldOpers, 
	     JacDetVInv(a_detJ, a_Xg, a_Yg, a_Zg, stream));
}

void Opers::
JacDetH(Image3D& a_jdet, const Field3D& a_h, 
	DiffT diffType, BoundaryCondT bc,
        StreamT stream)
{
   MK_CHECK2_ALL(a_jdet, a_h);
   AUTO_EXEC(a_jdet.memType(), ImageFieldOpers, 
	     JacDetH(a_jdet, a_h, diffType, bc, stream));
}

// template instantiations

#include "IFOpers_inst.cxx"

} // end namespace PyCA

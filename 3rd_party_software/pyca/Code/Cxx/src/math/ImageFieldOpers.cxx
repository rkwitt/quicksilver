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

#include <ImageFieldOpers.h>

#include <conditionMacro.h>
#include <Image3D.h>
#include <Field3D.h>

#include "CImageFieldOpers.h"
#include "GImageFieldOpers.h"

namespace PyCA {

/*
 * apply hField to an image
 *  defImage(x) = image(h(x))
 */
template<int mode>
template<BackgroundStrategy bg, InterpT interp>
void ImageFieldOpers<mode>::ApplyH(Image3D& a_o, const Image3D& a_i, const Field3D& a_h, StreamT s){
    MK_CHECK3_MEM(a_o, a_i, a_h);
    Executer::template ApplyH<bg, interp>
	(a_o, a_i, a_h, s);
}

/*
 *  apply uField to an image
 *  defImage(x) = image(x + delta * u(x))
 */
template<int mode>
template<BackgroundStrategy bg, InterpT interp>
void ImageFieldOpers<mode>::ApplyV(Image3D& a_o, const Image3D& a_i, const Field3D& a_u, const float& delta, StreamT s, bool onDev) {
    MK_CHECK3_MEM(a_o, a_i, a_u);
    Executer::template ApplyV<bg, interp>
	(a_o, a_i, a_u, delta, s, onDev);
}
    
/*
 *  apply uField to an image
 *  defImage(x) = image(x - delta * u(x))
 */
template<int mode>
template<BackgroundStrategy bg, InterpT interp>
void ImageFieldOpers<mode>:: ApplyVInv(Image3D& a_o, const Image3D& a_i, const Field3D& a_u, const float& delta, StreamT s, bool onDev) {
    MK_CHECK3_MEM(a_o, a_i, a_u);
    Executer::template ApplyVInv<bg, interp>
	(a_o, a_i, a_u, delta, s, onDev);
}

/*
 * a_o = a_i ( x + t)
 */
template<int mode>
template<BackgroundStrategy bg, InterpT interp>
void ImageFieldOpers<mode>:: ComposeTranslation(Image3D& a_o, const Image3D& a_i, const Vec3Df& t, StreamT s, bool onDev) {
    MK_CHECK2_MEM(a_o, a_i);
    Executer::template ComposeTranslation<bg, interp>
	(a_o, a_i, t, s, onDev);
}
    
template<int mode>
void ImageFieldOpers<mode>:: Splat(float* a_o, const Field3D& a_h, const float* a_i, bool normalize, StreamT s) {
    Executer::Splat(a_o, a_h, a_i, normalize, s);
}

template<int mode>
void ImageFieldOpers<mode>:: Splat(Image3D& a_o, const Field3D& a_h, const Image3D& a_i, bool normalize, StreamT s) {
    MK_CHECK3_MEM(a_o, a_h, a_i);
    Executer::Splat(a_o, a_h, a_i, normalize, s);
}
    
/**
 * Compute finite difference in the dimension indicated
 *  a_o = diff(a_i)
 */
template<int mode>
void 
ImageFieldOpers<mode>::
FiniteDiff(Image3D& a_o, const Image3D& a_i, 
	   DimT dim,  DiffT diffType, 
	   enum BoundaryCondT bc, 
	   bool accum, OpT op,
	   StreamT s)
{
   Executer::FiniteDiff(a_o, a_i, dim, diffType, bc, accum, op, s);
}

/**
 * Compute the gradient of an image
 *  d_o = grad (a_i)
 */
template<int mode>
void ImageFieldOpers<mode>:: Gradient(Field3D& a_o, const float* a_i, DiffT diffType, 
				      BoundaryCondT bc, StreamT s)
{
   Executer::Gradient(a_o, a_i, diffType, bc, s);
}

template<int mode>
void ImageFieldOpers<mode>:: Gradient(Field3D& a_o, const Image3D& a_i, DiffT diffType, 
				      BoundaryCondT bc, StreamT s) 
{
    MK_CHECK2_MEM(a_o, a_i);
    Executer::Gradient(a_o, a_i, diffType, bc, s);
}

template<int mode>
void ImageFieldOpers<mode>:: Gradient2(Field3D& a_o, const Image3D& a_i, DiffT diffType, 
				       BoundaryCondT bc, StreamT s) 
{
    MK_CHECK2_MEM(a_o, a_i);
    Executer::Gradient2(a_o, a_i, diffType, bc, s);
}

template<int mode>
void ImageFieldOpers<mode>:: 
GradientMask(Field3D& a_o, 
	     const Image3D& a_i, 
	     const Image3D& a_mask, 
	     DiffT diffType, 
	     BoundaryCondT bc, 
	     StreamT s) 
{
   MK_CHECK3_MEM(a_o, a_i, a_mask);
   Executer::GradientMask(a_o, a_i, a_mask, diffType, bc, s);
}

template<int mode>
void ImageFieldOpers<mode>::
GradFor(Field3D& a_o, const Image3D& a_i, BoundaryCondT bc, StreamT s){
   MK_CHECK2_ALL(a_o, a_i);
   Executer::GradFor(a_o, a_i, bc, s);
}

template<int mode>
void ImageFieldOpers<mode>:: GradForMag(Image3D& a_o, const Image3D& a_i, StreamT s) {
    MK_CHECK2_MEM(a_o, a_i);
    Executer::GradForMag(a_o, a_i, s);
}

template<int mode>
void ImageFieldOpers<mode>:: GradientMag(Image3D& a_o, const Image3D& a_i, DiffT diffType, 
					 BoundaryCondT bc, StreamT s) 
{
    MK_CHECK2_MEM(a_o, a_i);
    Executer::GradientMag(a_o, a_i, diffType, bc, s);
}

template<int mode>
void ImageFieldOpers<mode>::
UpwindDiff(Image3D& a_o, const Image3D& a_i, 
	   const Image3D& a_speed,
	   DimT dim,
	   StreamT s)
{
   MK_CHECK3_MEM(a_o, a_i, a_speed);
   Executer::UpwindDiff(a_o, a_i, a_speed, dim, s);
}

template<int mode>
void ImageFieldOpers<mode>::
UpwindGradMag(Image3D& a_o, const Image3D& a_i,
	      const Image3D& a_speed, StreamT s)
{
   MK_CHECK3_MEM(a_o, a_i, a_speed);
   Executer::UpwindGradMag(a_o, a_i, a_speed, s);
}

/**
 * Compute the diverence of a field
 *  d_o = div (a_i)
 */
template<int mode>
void ImageFieldOpers<mode>:: Divergence(Image3D& a_o, const Field3D& a_i, DiffT diffType, 
					BoundaryCondT bc, StreamT s) 
{
    MK_CHECK2_MEM(a_o, a_i);
    Executer::Divergence(a_o, a_i, diffType, bc, s);
}

/**
 * Compute the diverence of a field using backward differences
 *  d_o = div (a_i)
 */
template<int mode>
void ImageFieldOpers<mode>:: 
DivBack(Image3D& a_o, const Field3D& a_i, BoundaryCondT bc, StreamT s){
   MK_CHECK2_ALL(a_o, a_i);
   Executer::DivBack(a_o, a_i, bc, s);
}

/**
 * Compute the magnitude image
 * a_o[i] = sqrt(a_i[i].x^2 + a_i[i].y^2 + a_i[i].z^2) 
 */
template<int mode>
void ImageFieldOpers<mode>:: Magnitude(Image3D& a_o, const Field3D& a_i, StreamT s) {
    MK_CHECK2_MEM(a_o, a_i);
    Executer::Magnitude(a_o, a_i, s);
}

/**
 * Compute the magnitude array
 * a_o[i] = a_i[i].x^2 + a_i[i].y^2 + a_i[i].z^2 
 */
template<int mode>
void ImageFieldOpers<mode>:: SqrMagnitude(Image3D& a_o, const Field3D& a_i, StreamT s) {
    MK_CHECK2_MEM(a_o, a_i);
    Executer::SqrMagnitude(a_o, a_i, s);
}

/**
 * Compute dot product array 
 * d_o[i] = d_i[i].x * d_i1[i].x + d_i[i].y * d_i1[i].y + d_i[i].z * d_i1[i].z
 */
template<int mode>
void ImageFieldOpers<mode>:: ComponentDotProd(Image3D& a_o, const Field3D& a_i, const Field3D& a_i1, StreamT s) {
    MK_CHECK3_MEM(a_o, a_i, a_i1);
    Executer::ComponentDotProd(a_o, a_i, a_i1, s);
}
    
/** @brief a_o.x = a_i.x + a_i1, a_o.y = a_i.y + a_i1, a_o.z = a_i.z + a_i1 */
template<int mode>
void ImageFieldOpers<mode>:: Add(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, StreamT s) {
    MK_CHECK3_MEM(a_o, a_i, a_i1);
    Executer::Add(a_o, a_i, a_i1, s);
}

/** @brief a_o.x = a_i.x - a_i1, a_o.y = a_i.y - a_i1, a_o.z = a_i.z - a_i1 */
template<int mode>
void ImageFieldOpers<mode>:: Sub(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, StreamT s){
    MK_CHECK3_MEM(a_o, a_i, a_i1);
    Executer::Sub(a_o, a_i, a_i1, s);
}

/** @brief a_o.x = a_i.x * a_i1, a_o.y = a_i.y * a_i1, a_o.z = a_i.z * a_i1 */
template<int mode>
void ImageFieldOpers<mode>:: Mul(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, StreamT s){
    MK_CHECK3_MEM(a_o, a_i, a_i1);
    Executer::Mul(a_o, a_i, a_i1, s);
}

/** @brief a_o.x = a_i.x / a_i1, a_o.y = a_i.y / a_i1, a_o.z = a_i.z / a_i1 */
template<int mode>
void ImageFieldOpers<mode>:: Div(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, StreamT s){
    MK_CHECK3_MEM(a_o, a_i, a_i1);
    Executer::Div(a_o, a_i, a_i1, s);
}

/** @brief a_o.x += a_i, a_o.y += a_i, a_o.y += a_i, */
template<int mode>
void ImageFieldOpers<mode>:: Add_I(Field3D& a_o, const Image3D& a_i, StreamT s) {
    MK_CHECK2_MEM(a_o, a_i);
    Executer::Add_I(a_o, a_i, s);
}

/** @brief a_o.x -= a_i, a_o.y -= a_i, a_o.y -= a_i, */
template<int mode>
void ImageFieldOpers<mode>:: Sub_I(Field3D& a_o, const Image3D& a_i, StreamT s){
    MK_CHECK2_MEM(a_o, a_i);
    Executer::Sub_I(a_o, a_i, s);
}

/** @brief a_o.x *= a_i, a_o.y *= a_i, a_o.y *= a_i, */
template<int mode>
void ImageFieldOpers<mode>:: Mul_I(Field3D& a_o, const Image3D& a_i, StreamT s) {
    MK_CHECK2_MEM(a_o, a_i);
    Executer::Mul_I(a_o, a_i, s);
}

/** @brief a_o.x /= a_i, a_o.y = a_i.y / a_i, a_o.z = a_i.z / a_i */
template<int mode>
void ImageFieldOpers<mode>:: Div_I(Field3D& a_o, const Image3D& a_i, StreamT s) {
    MK_CHECK2_MEM(a_o, a_i);
    Executer::Div_I(a_o, a_i, s);
}

/** @brief a_o = a_i + a_i1 * a_i2 */
template<int mode>
void ImageFieldOpers<mode>:: Add_Mul(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const Image3D& a_i2, StreamT s){
    MK_CHECK4_MEM(a_o, a_i, a_i1, a_i2);
    Executer::Add_Mul(a_o, a_i, a_i1, a_i2, s);
}

/** @brief a_o = a_o + a_i * a_i1 */
template<int mode>
void ImageFieldOpers<mode>:: Add_Mul_I(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, StreamT s){
    MK_CHECK3_MEM(a_o, a_i, a_i1);
    Executer::Add_Mul_I(a_o, a_i, a_i1, s);
}

/** @brief a_o = a_i - a_i1 * a_i2 */
template<int mode>
void ImageFieldOpers<mode>:: Sub_Mul(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const Image3D& a_i2, StreamT s){
    MK_CHECK4_MEM(a_o, a_i, a_i1, a_i2);
    Executer::Sub_Mul(a_o, a_i, a_i1, a_i2, s);
}

/** @brief a_o = a_o - a_i * a_i1 */
template<int mode>
void ImageFieldOpers<mode>:: Sub_Mul_I(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, StreamT s){
    MK_CHECK3_MEM(a_o, a_i, a_i1);
    Executer::Sub_Mul_I(a_o, a_i, a_i1, s);
}

/** @brief a_o = a_i * a_i1 * c (a_o.x = a_i.x * a_i1 * c)*/
template<int mode>
void ImageFieldOpers<mode>:: MulMulC(Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, const float& c, StreamT s, bool onDev){
    MK_CHECK3_MEM(a_o, a_i, a_i1);
    Executer::MulMulC(a_o, a_i, a_i1, c, s, onDev);
}

/** @brief a_o = a_o * a_i * c  (a_o.x = a_o.x * a_i * c)*/
template<int mode>
void ImageFieldOpers<mode>:: MulMulC_I(Field3D& a_o, const Image3D& a_i, const float& c, StreamT s, bool onDev){
    MK_CHECK2_MEM(a_o, a_i);
    Executer::MulMulC_I(a_o, a_i, c, s, onDev);
}

/** @brief a_o = a_i + a_i1 * a_i2 * c */
template<int mode>
void ImageFieldOpers<mode>:: Add_MulMulC(Field3D& a_o, const Field3D& a_i, const Field3D& a_i1, const Image3D& a_i2, const float& c,  StreamT s, bool onDev){
    MK_CHECK4_MEM(a_o, a_i, a_i1, a_i2);
    Executer::Add_MulMulC(a_o, a_i, a_i1, a_i2, c, s, onDev);
}

/** @brief a_o = a_o + a_i * a_i1 * c */
template<int mode>
void ImageFieldOpers<mode>:: Add_MulMulC_I(const Field3D& a_o, const Field3D& a_i, const Image3D& a_i1, const float& c,  StreamT s, bool onDev){
    MK_CHECK3_MEM(a_o, a_i, a_i1);
    Executer::Add_MulMulC_I(a_o, a_i, a_i1, c, s, onDev);
}
    
template<int mode>
void ImageFieldOpers<mode>:: JacDetH(Image3D& a_detJ,const Field3D& a_Xg, const Field3D& a_Yg, const Field3D& a_Zg, StreamT s){
    MK_CHECK4_MEM(a_detJ, a_Xg, a_Yg, a_Zg);
    Executer::JacDetH(a_detJ, a_Xg, a_Yg, a_Zg, s);
}

template<int mode>
void ImageFieldOpers<mode>:: JacDetV(Image3D& a_detJ,const Field3D& a_Xg, const Field3D& a_Yg, const Field3D& a_Zg, StreamT s){
    MK_CHECK4_MEM(a_detJ, a_Xg, a_Yg, a_Zg);
    Executer::JacDetV(a_detJ, a_Xg, a_Yg, a_Zg, s);
}

template<int mode>
void ImageFieldOpers<mode>:: JacDetVInv(Image3D& a_detJ,const Field3D& a_Xg,const Field3D& a_Yg,const Field3D& a_Zg,StreamT s){
    MK_CHECK4_MEM(a_detJ, a_Xg, a_Yg, a_Zg);
    Executer::JacDetVInv(a_detJ, a_Xg, a_Yg, a_Zg, s);
}

template<int mode>
void ImageFieldOpers<mode>:: JacDetH(Image3D& a_jdet,const Field3D& a_h, DiffT diffType, BoundaryCondT bc, StreamT s){
   MK_CHECK2_MEM(a_jdet, a_h);
   Executer::JacDetH(a_jdet, a_h, diffType, bc, s);
}

// template instantiations
#include "ImageFieldOpers_inst.cxx"

} // end namespace PyCA

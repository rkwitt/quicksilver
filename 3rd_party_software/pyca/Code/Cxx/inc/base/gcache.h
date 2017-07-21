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

#ifndef __CUDA_TEXFETCH_H
#define __CUDA_TEXFETCH_H

namespace PyCA {

// Texture cache for array
texture<unsigned char, 1> com_tex_uchar_x;
texture<signed char, 1> com_tex_char_x;
texture<unsigned short, 1> com_tex_ushort_x;
texture<short, 1> com_tex_short_x;

texture<float, 1> com_tex_float_x;
texture<int, 1> com_tex_int_x;
texture<uint, 1> com_tex_uint_x;

texture<float, 1> com_tex_float_y;
texture<int, 1> com_tex_int_y;
texture<uint, 1> com_tex_uint_y;

texture<float, 1> com_tex_float_z;
texture<int, 1> com_tex_int_z;
texture<uint, 1> com_tex_uint_z;

texture<float2, 1> com_tex_float2_x;
texture<int2, 1> com_tex_int2_x;
texture<uint2, 1> com_tex_uint2_x;

texture<float2, 1> com_tex_float2_y;
texture<int2, 1> com_tex_int2_y;
texture<uint2, 1> com_tex_uint2_y;

texture<float2, 1> com_tex_float2_z;
texture<int2, 1> com_tex_int2_z;
texture<uint2, 1> com_tex_uint2_z;

texture<int2, 1> com_tex_double_x;
texture<int2, 1> com_tex_double_y;
texture<int2, 1> com_tex_double_z;

texture<float4, 1> com_tex_float4_x;
texture<int4, 1> com_tex_int4_x;
texture<uint4, 1> com_tex_uint4_x;

texture<float4, 1> com_tex_float4_y;
texture<int4, 1> com_tex_int4_y;
texture<uint4, 1> com_tex_uint4_y;

texture<float4, 1> com_tex_float4_z;
texture<int4, 1> com_tex_int4_z;
texture<uint4, 1> com_tex_uint4_z;

inline void cache_bind(const float* d_x, const float* d_y, const float* d_z)
{
    cudaBindTexture(NULL, com_tex_float_x, d_x);
    cudaBindTexture(NULL, com_tex_float_y, d_y);
    cudaBindTexture(NULL, com_tex_float_z, d_z);
}

template<typename T>
inline void cache_bind(const T* d_i){
    cudaBindTexture(NULL, com_tex_float_x, d_i);
}

template<>
inline void cache_bind(const unsigned char* d_i){
    cudaBindTexture(NULL, com_tex_uchar_x, d_i);
}

template<>
inline void cache_bind(const char* d_i){
    cudaBindTexture(NULL, com_tex_char_x, d_i);
}

template<>
inline void cache_bind(const unsigned short* d_i){
    cudaBindTexture(NULL, com_tex_ushort_x, d_i);
}

template<>
inline void cache_bind(const short* d_i){
    cudaBindTexture(NULL, com_tex_short_x, d_i);
}


template<>
inline void cache_bind(const float* d_i){
    cudaBindTexture(NULL, com_tex_float_x, d_i);
}

template<>
inline void cache_bind(const float2* d_i){
    cudaBindTexture(NULL, com_tex_float2_x, d_i);
}

template<>
inline void cache_bind(const float4*d_i){
    cudaBindTexture(NULL, com_tex_float4_x, d_i);
}

template<>
inline void cache_bind(const int* d_i){
    cudaBindTexture(NULL, com_tex_int_x, d_i);
}

template<>
inline void cache_bind(const int2* d_i){
    cudaBindTexture(NULL, com_tex_int2_x, d_i);
}

template<>
inline void cache_bind(const int4*d_i){
    cudaBindTexture(NULL, com_tex_int4_x, d_i);
}

template<>
inline void cache_bind(const uint* d_i){
    cudaBindTexture(NULL, com_tex_uint_x, d_i);
}

template<>
inline void cache_bind(const uint2* d_i){
    cudaBindTexture(NULL, com_tex_uint2_x, d_i);
}

template<>
inline void cache_bind(const uint4*d_i){
    cudaBindTexture(NULL, com_tex_uint4_x, d_i);
}

template<>
inline void cache_bind(const double*d_i){
    cudaBindTexture(NULL, com_tex_double_x, (int2*) d_i);
}


__inline__ __device__ float fetchf(uint i){
    return tex1Dfetch(com_tex_float_x, i);
}

__inline__ __device__ int fetchi(uint i){
    return tex1Dfetch(com_tex_int_x, i);
}

__inline__ __device__ uint fetchu(uint i){
    return  tex1Dfetch(com_tex_uint_x, i);
}

template<typename T>
__inline__ __device__ T fetch(uint i, const T* d_i){
    return (T) tex1Dfetch(com_tex_float_x, i);
}

template<>
__inline__ __device__ unsigned char fetch(uint i, const unsigned char* d_i){
    return tex1Dfetch(com_tex_uchar_x, i);
}
template<>
__inline__ __device__ char fetch(uint i, const char* d_i){
    return tex1Dfetch(com_tex_char_x, i);
}

template<>
__inline__ __device__ unsigned short fetch(uint i, const unsigned short* d_i){
    return tex1Dfetch(com_tex_ushort_x, i);
}
template<>
__inline__ __device__ short fetch(uint i, const short* d_i){
    return tex1Dfetch(com_tex_short_x, i);
}

template<>
__inline__ __device__ float fetch(uint i, const float* d_i){
    return tex1Dfetch(com_tex_float_x, i);
}

template<>
__inline__ __device__ float2 fetch(uint i, const float2* d_i){
    return tex1Dfetch(com_tex_float2_x, i);
}

template<>
__inline__ __device__ float4 fetch(uint i, const float4* d_i){
    return tex1Dfetch(com_tex_float4_x, i);
}

template<>
__inline__ __device__ uint fetch(uint i, const uint* d_i){
    return tex1Dfetch(com_tex_uint_x, i);
}

template<>
__inline__ __device__ uint2 fetch(uint i, const uint2* d_i){
    return tex1Dfetch(com_tex_uint2_x, i);
}

template<>
__inline__ __device__ uint4 fetch(uint i, const uint4* d_i){
    return tex1Dfetch(com_tex_uint4_x, i);
}

template<>
__inline__ __device__ int fetch(uint i, const int* d_i){
    return tex1Dfetch(com_tex_int_x, i);
}

template<>
__inline__ __device__ int2 fetch(uint i, const int2* d_i){
    return tex1Dfetch(com_tex_int2_x, i);
}

template<> __inline__ __device__ int4 fetch(uint i, const int4* d_i){
    return tex1Dfetch(com_tex_int4_x, i);
}

template<>
__inline__ __device__ double fetch(uint i, const double* d_i){
    int2 v =  tex1Dfetch(com_tex_double_x, i);
    return __hiloint2double(v.x, v.y);
}

#define cache_bind_x cache_bind
#define fetch_x      fetch

// Function for the second cache
template<typename T>
inline void cache_bind_y(const T* d_i){
    cudaBindTexture(NULL, com_tex_float_y, d_i);
}

template<>
inline void cache_bind_y(const float* d_i){
    cudaBindTexture(NULL, com_tex_float_y, d_i);
}
template<>
inline void cache_bind_y(const int* d_i){
    cudaBindTexture(NULL, com_tex_int_y, d_i);
}
template<>
inline void cache_bind_y(const uint* d_i){
    cudaBindTexture(NULL, com_tex_uint_y, d_i);
}

template<>
inline void cache_bind_y(const double* d_i){
    cudaBindTexture(NULL, com_tex_double_y, (int2*) d_i);
}

template<typename T> void cache_bind_xy(const T* d_x, const T* d_y)
{
    cache_bind_x(d_x);
    cache_bind_y(d_y);
}

////////////////////////////////////////////////////////////////////////////////
template<typename T>
__inline__ __device__ T fetch_y(uint i, const T* d_i){
    return (T) tex1Dfetch(com_tex_float_y, i);
}

template<>
__inline__ __device__ float fetch_y(uint i, const float* d_i){
    return tex1Dfetch(com_tex_float_y, i);
}

template<>
__inline__ __device__ int fetch_y(uint i, const int* d_i){
    return tex1Dfetch(com_tex_int_y, i);
}

template<>
__inline__ __device__ uint fetch_y(uint i, const uint* d_i){
    return tex1Dfetch(com_tex_uint_y, i);
}

template<>
__inline__ __device__ float2 fetch_y(uint i, const float2* d_i){
    return tex1Dfetch(com_tex_float2_y, i);
}

template<>
__inline__ __device__ int2 fetch_y(uint i, const int2* d_i){
    return tex1Dfetch(com_tex_int2_y, i);
}

template<>
__inline__ __device__ uint2 fetch_y(uint i, const uint2* d_i){
    return tex1Dfetch(com_tex_uint2_y, i);
}

template<>
__inline__ __device__ double fetch_y(uint i, const double* d_i){
    int2 v =  tex1Dfetch(com_tex_double_y, i);
    return __hiloint2double(v.x, v.y);
}


// Function for the third cache
template<typename T>
inline void cache_bind_z(const T* d_i){
    cudaBindTexture(NULL, com_tex_float_z, d_i);
}

template<>
inline void cache_bind_z(const float* d_i){
    cudaBindTexture(NULL, com_tex_float_z, d_i);
}

template<>
inline void cache_bind_z(const int* d_i){
    cudaBindTexture(NULL, com_tex_int_z, d_i);
}
template<>
inline void cache_bind_z(const uint* d_i){
    cudaBindTexture(NULL, com_tex_uint_z, d_i);
}

template<>
inline void cache_bind_z(const double* d_i){
    cudaBindTexture(NULL, com_tex_double_z, (int2*) d_i);
}

template<typename T>
__inline__ __device__ T fetch_z(uint i, const T* d_i){
    return (T) tex1Dfetch(com_tex_float_z, i);
}

template<>
__inline__ __device__ float fetch_z(uint i, const float* d_i){
    return tex1Dfetch(com_tex_float_z, i);
}

template<>
__inline__ __device__ int fetch_z(uint i, const int* d_i){
    return tex1Dfetch(com_tex_int_z, i);
}

template<>
__inline__ __device__ uint fetch_z(uint i, const uint* d_i){
    return tex1Dfetch(com_tex_uint_z, i);
}

template<>
__inline__ __device__ double fetch_z(uint i, const double* d_i){
    int2 v =  tex1Dfetch(com_tex_double_z, i);
    return __hiloint2double(v.x, v.y);
}


//Fetching with template
template<typename T>
__inline__ __device__ T fetch(int i){
    return (T) tex1Dfetch(com_tex_float_x, i);
}

template<>
__inline__ __device__ char fetch<char>(int i) {
    return tex1Dfetch(com_tex_char_x, i);
}

template<>
__inline__ __device__ unsigned char fetch<unsigned char>(int i) {
    return tex1Dfetch(com_tex_uchar_x, i);
}

template<>
__inline__ __device__ short fetch<short>(int i) {
    return tex1Dfetch(com_tex_short_x, i);
}

template<>
__inline__ __device__ unsigned short fetch<unsigned short>(int i) {
    return tex1Dfetch(com_tex_ushort_x, i);
}

template<>
__inline__ __device__ float fetch<float>(int i) {
    return tex1Dfetch(com_tex_float_x, i);
}

template<>
__inline__ __device__ float2 fetch<float2>(int i) {
    return tex1Dfetch(com_tex_float2_x, i);
}

template<>
__inline__ __device__ float4 fetch<float4>(int i) {
    return tex1Dfetch(com_tex_float4_x, i);
}

template<>
__inline__ __device__ int fetch<int>(int i) {
    return tex1Dfetch(com_tex_int_x, i);
}

template<>
__inline__ __device__ int2 fetch<int2>(int i) {
    return tex1Dfetch(com_tex_int2_x, i);
}

template<>
__inline__ __device__ int4 fetch<int4>(int i) {
    return tex1Dfetch(com_tex_int4_x, i);
}

template<>
__inline__ __device__ uint fetch<uint>(int i) {
    return tex1Dfetch(com_tex_uint_x, i);
}

template<>
__inline__ __device__ uint2 fetch<uint2>(int i) {
    return tex1Dfetch(com_tex_uint2_x, i);
}

template<>
__inline__ __device__ uint4 fetch<uint4>(int i) {
    return tex1Dfetch(com_tex_uint4_x, i);
}

template<>
__inline__ __device__ double fetch<double>(int i) {
    int2 v =  tex1Dfetch(com_tex_double_x, i);
    return __hiloint2double(v.x, v.y);
}

} // end namespace PyCA

#endif

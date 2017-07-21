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

#include <cstdio>

#include "GReduceStream.h"
#include "ReduceStreamKernel.cu"
#include <pycaUtils.h>
#include <mem.h>
#include <gmem.h>
#include <Vec2D.h>

// TEST make sure boost isn't included in nvcc code
#if defined(BOOST_COMPILER)
int bla[-1];
#endif

namespace PyCA {

GReduceStream::GReduceStream(){
    int *d_buf;
    dmemAlloc(d_buf, REDUCE_STREAM_SIZE * MAX_NUMBER_OF_REDUCE_STREAM);
    mdBuf = d_buf;
}


GReduceStream::~GReduceStream() {
    if (mdBuf != NULL) {
        dmemFree(mdBuf);
        mdBuf = NULL;
    }
}

////////////////////////////////////////////////////////////////////////////////
// Single operator single output function
// d_o = reduce(oper, d_i);  
////////////////////////////////////////////////////////////////////////////////
template<typename T, MATH_OPS op>
void GReduceStream::reduce(T* d_o, const T* d_i, size_t n, bool update, StreamT stream){
    reduce_L1_kernel<T, MOpers<T, op> >
        <<<M, N, 0, stream>>>((T*)mdBuf, d_i, n);
    if (update) {
        reduce_L2_kernel<T, MOpers<T, op>, true><<<1, N, 0, stream>>>(d_o, (T*) mdBuf);
    } else {
        reduce_L2_kernel<T, MOpers<T, op>, false><<<1, N, 0, stream>>>(d_o, (T*) mdBuf);
    }
}

////////////////////////////////////////////////////////////////////////////////
// Double operator single input single output function
// d_o = reduce(oper, oper1(d_i))
////////////////////////////////////////////////////////////////////////////////

template <typename T, MATH_OPS op, MATH_OPS op1>
void GReduceStream::compReduce(T* d_o, const T* d_i, size_t n, bool update, StreamT stream)
{
    compReduce_L1_kernel<T, MOpers<T, op>, MOpers<T, op1> >
        <<<M, N, 0, stream>>>((T*)mdBuf, d_i, n);
    if (update)
        reduce_L2_kernel<T, MOpers<T, op>, true><<<1, N, 0, stream>>>(d_o, (T*)mdBuf);
    else
        reduce_L2_kernel<T, MOpers<T, op>, false><<<1, N, 0, stream>>>(d_o, (T*)mdBuf);
}


////////////////////////////////////////////////////////////////////////////////
// Double operator Double output single input function
//      d_o[0] = reduce(op, d_i)
//      d_o[1] = reduce(op1, d_i)
////////////////////////////////////////////////////////////////////////////////

template <typename T, MATH_OPS op, MATH_OPS op1>
void GReduceStream::bireduce(T* d_o, const T* d_i, size_t n, bool update, StreamT stream)
{
   PRECONDITION(d_i!= NULL, "null pointer");

    T* d_buf  = (T*) mdBuf;
    T* d_buf1 = d_buf + REDUCE_STREAM_SIZE;

    bireduce_L1_kernel<T, MOpers<T,op>, MOpers<T,op1> >
        <<<M, N, 0, stream>>>(d_buf, d_buf1, d_i, n);
    
    if (update)
        reduce_ip2_op2_L2_kernel<T, MOpers<T,op>, MOpers<T,op1>, true>
            <<<1, N, 0, stream>>>(d_o, d_buf, d_buf1);
    else
        reduce_ip2_op2_L2_kernel<T, MOpers<T,op>, MOpers<T,op1>, false>
            <<<1, N, 0, stream>>>(d_o, d_buf, d_buf1);
}

////////////////////////////////////////////////////////////////////////////////
// Double operator single output double input function
//      d_o[0] = reduce(op, op1(d_i, d_i1))
////////////////////////////////////////////////////////////////////////////////
template <typename T, MATH_OPS op, MATH_OPS op1>
void GReduceStream::product(T* d_o, const T* d_i0, const T* d_i1, size_t n, bool update, StreamT stream)
{
   PRECONDITION((d_i0!= NULL) && (d_i1!= NULL), "null pointer");
    product_L1_kernel<T, MOpers<T,op>, MOpers<T,op1> >
        <<<M, N, 0, stream>>>( (T*)mdBuf, d_i0, d_i1, n);
    if (update)
        reduce_L2_kernel<T, MOpers<T,op>, true>
            <<<1, N, 0, stream>>>(d_o, (T*)mdBuf);
    else
        reduce_L2_kernel<T, MOpers<T,op>, false>
            <<<1, N, 0, stream>>>(d_o, (T*)mdBuf);
}


////////////////////////////////////////////////////////////////////////////////
// Instantiate for implementation
////////////////////////////////////////////////////////////////////////////////
template<typename T>
void GReduceStream::Max(T& d_o, const T* d_i, size_t n, bool update, StreamT stream){
    reduce<T, MATH_Max>((T*)&d_o, d_i, n, update, stream);
}

template void GReduceStream::Max(float& d_o, const float* d_i, size_t n, bool update, StreamT stream);
template void GReduceStream::Max(int& d_o, const int* d_i, size_t n, bool update, StreamT stream);
template void GReduceStream::Max(uint& d_o, const uint* d_i, size_t n, bool update, StreamT stream);

////////////////////////////////////////////////////////////////////////////////
template<typename T>
void GReduceStream::Min(T& d_o, const T* d_i, size_t n, bool update, StreamT stream){
    reduce<T, MATH_Min>((T*)&d_o, d_i, n, update, stream);
}

template void GReduceStream::Min(float& d_o, const float* d_i, size_t n, bool update, StreamT stream);
template void GReduceStream::Min(int& d_o, const int* d_i, size_t n, bool update, StreamT stream);
template void GReduceStream::Min(uint& d_o, const uint* d_i, size_t n, bool update, StreamT stream);

////////////////////////////////////////////////////////////////////////////////
template<typename T>
void GReduceStream::Sum(T& d_o, const T* d_i, size_t n, bool update, StreamT stream){
    reduce<T, MATH_Add>((T*)&d_o, d_i, n, update, stream);
}

template void GReduceStream::Sum(float& d_o, const float* d_i, size_t n, bool update, StreamT stream);
template void GReduceStream::Sum(int& d_o, const int* d_i, size_t n, bool update, StreamT stream);
template void GReduceStream::Sum(uint& d_o, const uint* d_i, size_t n, bool update, StreamT stream);


template<typename T>
void GReduceStream::LInf(T& d_o, const T* d_i, size_t n, bool update, StreamT stream){
    compReduce<T, MATH_Max, MATH_Abs>((T*)&d_o, d_i, n, update, stream);
};

template void GReduceStream::LInf(float& d_o, const float* d_i, size_t n, bool update, StreamT stream);
template void GReduceStream::LInf(int& d_o, const int* d_i, size_t n, bool update, StreamT stream);
template void GReduceStream::LInf(uint& d_o, const uint* d_i, size_t n, bool update, StreamT stream);

template<typename T>
void GReduceStream::L1(T& d_o, const T* d_i, size_t n, bool update, StreamT stream){
    compReduce<T, MATH_Add, MATH_Abs>((T*)&d_o, d_i, n, update, stream);
};

template void GReduceStream::L1(float& d_o, const float* d_i, size_t n, bool update, StreamT stream);
template void GReduceStream::L1(int& d_o, const int* d_i, size_t n, bool update, StreamT stream);
template void GReduceStream::L1(uint& d_o, const uint* d_i, size_t n, bool update, StreamT stream);

template<typename T>
void GReduceStream::Sum2(T& d_o, const T* d_i, size_t n, bool update, StreamT stream){
    compReduce<T, MATH_Add, MATH_Sqr>((T*)&d_o, d_i, n, update, stream);
};

template void GReduceStream::Sum2(float& d_o, const float* d_i, size_t n, bool update, StreamT stream);
template void GReduceStream::Sum2(int& d_o, const int* d_i, size_t n, bool update, StreamT stream);
template void GReduceStream::Sum2(uint& d_o, const uint* d_i, size_t n, bool update, StreamT stream);

template<typename T>
void GReduceStream::MaxMin(Vec2D<T>&d_o, const T* d_i, size_t n, bool update, StreamT stream){
    bireduce<T, MATH_Max, MATH_Min>((T*)&d_o.x, d_i, n, update, stream);
}

template void GReduceStream::MaxMin(Vec2D<float>& d_o, const float* d_i, size_t n, bool update, StreamT stream);
template void GReduceStream::MaxMin(Vec2D<int>& d_o, const int* d_i, size_t n, bool update, StreamT stream);
template void GReduceStream::MaxMin(Vec2D<uint>& d_o, const uint* d_i, size_t n, bool update, StreamT stream);

template<typename T>
void GReduceStream::Dot(T& d_o, const T* d_i, const T* d_i1, size_t n, bool update, StreamT stream){
    product<T, MATH_Add, MATH_Mul>((T*)&d_o, d_i, d_i1, n, update, stream);
}

template
void GReduceStream::Dot(float& d_o, const float* d_i, const float* d_i1, size_t n, bool update, StreamT stream);
template
void GReduceStream::Dot(int& d_o, const int* d_i, const int* d_i1, size_t n, bool update, StreamT stream);
template
void GReduceStream::Dot(uint& d_o, const uint* d_i, const uint* d_i1, size_t n, bool update, StreamT stream);

bool GReduceStream::selfTest(size_t n){
    int test = true;
    int* h_i  = new int [n];
    int* h_i1 = new int [n];

    for (size_t j=0; j< n; ++j) h_i[j] = (rand() & 0xff);
    for (size_t j=0; j< n; ++j) h_i1[j] = (rand() & 0xff);
    
    int *d_i;
    dmemAlloc(d_i, n);
    cpyArrayH2D(d_i, h_i, n);
    
    int *d_i1;
    dmemAlloc(d_i1,n);
    cpyArrayH2D(d_i1,h_i1,n);

    int *d_o;
    dmemAlloc(d_o, 256);
    
    int h_max    = -INT_MAX, h_min = INT_MAX;
    int h_LInf   = 0;
    int h_sum2   = 0;
    int h_sum    = 0;
    int h_L1     = 0;
    int h_dot    = 0;
    
    for (size_t i=0; i< n; ++i)
    {
        h_max  = fmaxf(h_max, h_i[i]);
        h_min  = fminf(h_min, h_i[i]);
        h_sum += h_i[i];

        h_LInf = fmaxf(h_LInf, h_i[i]);
        h_L1  += fabsf(h_i[i]);
        h_sum2+= h_i[i]*h_i[i];
        
        h_dot += h_i1[i] * h_i[i];
    }

    int d_max = -INT_MAX, d_min = INT_MAX;
    int d_LInf= 0;
    int d_L1= 0;
    int d_sum2= 0;
    int d_sum = 0;
    int d_dot = 0;
    Vec2D<int> d_pair;
    Vec2D<int>* d_maxmin = createGObj(Vec2Di(0,0));
    
    this->Sum(d_o[0], d_i, n);cpyArrayD2H(&d_sum, d_o, 1);
    this->Max(d_o[0], d_i, n);cpyArrayD2H(&d_max, d_o, 1);
    this->Min(d_o[0], d_i, n);cpyArrayD2H(&d_min, d_o, 1);
    this->LInf(d_o[0],d_i, n);cpyArrayD2H(&d_LInf, d_o, 1);
    this->L1(d_o[0], d_i, n);cpyArrayD2H(&d_L1, d_o, 1);
    this->Sum2(d_o[0],d_i, n);cpyArrayD2H(&d_sum2, d_o, 1);
    this->Dot(d_o[0],d_i,d_i1, n);cpyArrayD2H(&d_dot, d_o, 1);
    this->MaxMin(*d_maxmin, d_i, n);d_pair = getGObj(d_maxmin);

    fprintf(stderr, "Maximum value from CPU %d from GPU %d\n",h_max, d_max);
    fprintf(stderr, "Minumum value from CPU %d from GPU %d\n",h_min, d_min);
    fprintf(stderr, "Total value from CPU %d from GPU %d\n",h_sum, d_sum);
    
    fprintf(stderr, "Maximum abosulte value from CPU %d from GPU %d\n",h_LInf, d_LInf);
    fprintf(stderr, "Total square value from CPU %d from GPU %d\n",h_sum2, d_sum2);
    fprintf(stderr, "Dot product value from CPU %d from GPU %d\n",h_dot, d_dot);
    fprintf(stderr, "Max Min value from CPU %d %d from GPU %d %d\n",h_max, h_min, d_pair.x, d_pair.y);
    
    //Extensive test
    h_max = -INT_MAX, h_min = INT_MAX;
    h_LInf = 0;
    h_sum2 = 0;
    h_sum = 0;
    h_dot = 0;
    h_L1 = 0;
    
    for (int l=1; l < 10001;++l){
        h_max = fmaxf(h_max, h_i[l-1]);
        h_LInf = fmaxf(h_LInf, h_i[l-1]);
        h_min = fminf(h_min, h_i[l-1]);
        h_sum += h_i[l-1];
        h_sum2 += h_i[l-1]*h_i[l-1];
        h_dot += h_i1[l-1] * h_i[l-1];
        h_L1  += fabsf(h_i[l-1]);
        
        this->Sum(d_o[0], d_i, l);cpyArrayD2H(&d_sum, d_o, 1);
        this->Max(d_o[0], d_i, l);cpyArrayD2H(&d_max, d_o, 1);
        this->Min(d_o[0], d_i, l);cpyArrayD2H(&d_min, d_o, 1);
        this->LInf(d_o[0],d_i, l);cpyArrayD2H(&d_LInf, d_o, 1);
        this->L1(d_o[0], d_i, l);cpyArrayD2H(&d_L1, d_o, 1);
        this->Sum2(d_o[0],d_i, l);cpyArrayD2H(&d_sum2, d_o, 1);
        this->Dot(d_o[0],d_i,d_i1, l);cpyArrayD2H(&d_dot, d_o, 1);
        this->MaxMin(*d_maxmin, d_i, l);d_pair = getGObj(d_maxmin);
        
        
        if (d_max != h_max){
            fprintf(stderr, "Max Test FAILED at %d GPU %d CPU %d\n", l, d_max, h_max );
            test = false;
        }

        if (d_min != h_min){
            fprintf(stderr, "Min Test FAILED at %d GPU %d CPU %d\n", l, d_min, h_min );
            test = false;
        }

        if (d_LInf!= h_LInf){
            fprintf(stderr, "MaxAbs Test FAILED at %d GPU %d CPU %d\n", l, d_LInf, h_LInf );
            test = false;
        }

        if (d_sum!= h_sum){
            fprintf(stderr, "Sum Test FAILED at %d GPU %d CPU %d\n", l, d_sum, h_sum );
            test = false;
        }

        if (d_sum2!= h_sum2){
            fprintf(stderr, "Sum SQR Test FAILED at %d GPU %d CPU %d\n", l, d_sum2, h_sum2 );
            test = false;
        }

        if (d_dot!= h_dot){
            fprintf(stderr, "Dot Test FAILED at %d GPU %d CPU %d\n", l, d_dot, h_dot );
            test = false;
        }
        
        if ( d_pair.x != h_max || d_pair.y != h_min){
            fprintf(stderr, "Max Min Test FAILED at %d GPU %d %d CPU %d %d\n",
                    l, d_pair.x, d_pair.y, h_max, h_min);
            test = false;
        }
        if (test == false)
            break;
    }

    if (test)
        fprintf(stderr, "Test PASSED  \n");
    
    delete []h_i1;
    delete []h_i;
    cudaFree(d_i);
    cudaFree(d_i1);
    return test;
}


} // end namespace PyCA

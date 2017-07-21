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

#include "CReduce.h"

#include <Vec2D.h>
#include <stdio.h>

#include <cstdlib>

namespace PyCA {

CReduce::CReduce(){};

template<typename T, MATH_OPS op>
void CReduce::reduce(T* h_o, const T* h_i, size_t n, bool update, StreamT stream){
    typedef class MOpers<T, op> oper;

    T res = oper::identity();
    for (size_t i=0; i< n; ++i) {
        oper::iop(res, h_i[i]);
    }

    h_o[0] = update ? oper::op(h_o[0], res) : res;
}

template <typename T, MATH_OPS op, MATH_OPS op1>
void CReduce::compReduce(T* h_o, const T* h_i, size_t n, bool update, StreamT stream){
    typedef class MOpers<T, op> oper;
    typedef class MOpers<T, op1> oper1;

    T res = oper::identity();
    for (size_t i=0; i< n; ++i) {
        oper::iop(res, oper1::op(h_i[i]));
    }

    h_o[0] = update ? oper::op(h_o[0], res) : res;
}

template <typename T, MATH_OPS op, MATH_OPS op1>
void CReduce::bireduce(T* h_o, const T* h_i, size_t n, bool update, StreamT stream){
    typedef class MOpers<T, op> oper;
    typedef class MOpers<T, op1> oper1;

    T res  = oper::identity();
    T res1 = oper1::identity();
    
    for (size_t i=0; i< n; ++i) {
        oper::iop(res, h_i[i]);
        oper1::iop(res1, h_i[i]);
    }

    h_o[0] = update ? oper::op(h_o[0], res) : res;
    h_o[1] = update ? oper1::op(h_o[1], res1) : res1;

}
    
template <typename T, MATH_OPS op, MATH_OPS op1>
void CReduce::product(T* h_o, const T* h_i0, const T* h_i1, size_t n, bool update, StreamT stream){
    typedef class MOpers<T, op> oper;
    typedef class MOpers<T, op1> oper1;

    T res = oper::identity();
    for (size_t i=0; i< n; ++i) {
        oper::iop(res, oper1::op(h_i0[i], h_i1[i]));
    }
    h_o[0] = update ? oper::op(h_o[0], res) : res;
}

////////////////////////////////////////////////////////////////////////////////
// Instantiate for implementation
////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CReduce::Max(T& h_o, const T* h_i, size_t n, bool update, StreamT stream){
    reduce<T, MATH_Max>((T*)&h_o, h_i, n, update, stream);
}

template void CReduce::Max(float& h_o, const float* h_i, size_t n, bool update, StreamT stream);
template void CReduce::Max(int& h_o, const int* h_i, size_t n, bool update, StreamT stream);
template void CReduce::Max(uint& h_o, const uint* h_i, size_t n, bool update, StreamT stream);

////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CReduce::Min(T& h_o, const T* h_i, size_t n, bool update, StreamT stream){
    reduce<T, MATH_Min>((T*)&h_o, h_i, n, update, stream);
}

template void CReduce::Min(float& h_o, const float* h_i, size_t n, bool update, StreamT stream);
template void CReduce::Min(int& h_o, const int* h_i, size_t n, bool update, StreamT stream);
template void CReduce::Min(uint& h_o, const uint* h_i, size_t n, bool update, StreamT stream);

////////////////////////////////////////////////////////////////////////////////
template<typename T>
void CReduce::Sum(T& h_o, const T* h_i, size_t n, bool update, StreamT stream){
    reduce<T, MATH_Add>((T*)&h_o, h_i, n, update, stream);
}

template void CReduce::Sum(float& h_o, const float* h_i, size_t n, bool update, StreamT stream);
template void CReduce::Sum(int& h_o, const int* h_i, size_t n, bool update, StreamT stream);
template void CReduce::Sum(uint& h_o, const uint* h_i, size_t n, bool update, StreamT stream);


template<typename T>
void CReduce::LInf(T& h_o, const T* h_i, size_t n, bool update, StreamT stream){
    compReduce<T, MATH_Max, MATH_Abs>((T*)&h_o, h_i, n, update, stream);
};

template void CReduce::LInf(float& h_o, const float* h_i, size_t n, bool update, StreamT stream);
template void CReduce::LInf(int& h_o, const int* h_i, size_t n, bool update, StreamT stream);
template void CReduce::LInf(uint& h_o, const uint* h_i, size_t n, bool update, StreamT stream);

template<typename T>
void CReduce::L1(T& h_o, const T* h_i, size_t n, bool update, StreamT stream){
    compReduce<T, MATH_Add, MATH_Abs>((T*)&h_o, h_i, n, update, stream);
};

template void CReduce::L1(float& h_o, const float* h_i, size_t n, bool update, StreamT stream);
template void CReduce::L1(int& h_o, const int* h_i, size_t n, bool update, StreamT stream);
template void CReduce::L1(uint& h_o, const uint* h_i, size_t n, bool update, StreamT stream);

template<typename T>
void CReduce::Sum2(T& h_o, const T* h_i, size_t n, bool update, StreamT stream){
    compReduce<T, MATH_Add, MATH_Sqr>((T*)&h_o, h_i, n, update, stream);
};

template void CReduce::Sum2(float& h_o, const float* h_i, size_t n, bool update, StreamT stream);
template void CReduce::Sum2(int& h_o, const int* h_i, size_t n, bool update, StreamT stream);
template void CReduce::Sum2(uint& h_o, const uint* h_i, size_t n, bool update, StreamT stream);

template<typename T>
void CReduce::MaxMin(Vec2D<T>&h_o, const T* h_i, size_t n, bool update, StreamT stream){
    bireduce<T, MATH_Max, MATH_Min>((T*)&h_o.x, h_i, n, update, stream);
}

template void CReduce::MaxMin(Vec2D<float>& h_o, const float* h_i, size_t n, bool update, StreamT stream);
template void CReduce::MaxMin(Vec2D<int>& h_o, const int* h_i, size_t n, bool update, StreamT stream);
template void CReduce::MaxMin(Vec2D<uint>& h_o, const uint* h_i, size_t n, bool update, StreamT stream);

template<typename T>
void CReduce::Dot(T& h_o, const T* h_i, const T* h_i1, size_t n, bool update, StreamT stream){
    product<T, MATH_Add, MATH_Mul>((T*)&h_o, h_i, h_i1, n, update, stream);
}

template
void CReduce::Dot(float& h_o, const float* h_i, const float* h_i1, size_t n, bool update, StreamT stream);
template
void CReduce::Dot(int& h_o, const int* h_i, const int* h_i1, size_t n, bool update, StreamT stream);
template
void CReduce::Dot(uint& h_o, const uint* h_i, const uint* h_i1, size_t n, bool update, StreamT stream);

bool CReduce::selfTest(size_t n){
    int test = true;
    int* h_i  = new int [n];
    int* h_i1 = new int [n];

    for (size_t j=0; j< n; ++j) h_i[j] = (rand() & 0xff);
    for (size_t j=0; j< n; ++j) h_i1[j] = (rand() & 0xff);
    
    int c_o[2];
    
    int h_max    = -INT_MAX, h_min = INT_MAX;
    int h_LInf   = 0;
    int h_sum2   = 0;
    int h_sum    = 0;
    int h_L1     = 0;
    int h_dot    = 0;
    
    for (size_t i=0; i< n; ++i)
    {
        h_max  = std::max(h_max, h_i[i]);
        h_min  = std::min(h_min, h_i[i]);
        h_sum += h_i[i];

        h_LInf = std::max(h_LInf, h_i[i]);
        h_L1  += abs(h_i[i]);
        h_sum2+= h_i[i]*h_i[i];
        
        h_dot += h_i1[i] * h_i[i];
    }

    int c_max = -INT_MAX, c_min = INT_MAX;
    int c_LInf= 0;
    int c_L1= 0;
    int c_sum2= 0;
    int c_sum = 0;
    int c_dot = 0;
    
    this->Sum(c_sum, h_i, n);
    this->Max(c_max, h_i, n);
    this->Min(c_min, h_i, n);
    this->LInf(c_LInf,h_i, n);
    this->L1(c_L1, h_i, n);
    this->Sum2(c_sum2,h_i, n);
    this->Dot(c_dot,h_i,h_i1, n);
    this->MaxMin(*(Vec2D<int>*)c_o, h_i, n);

    fprintf(stderr, "Maximum value from CPU %d from CPU2 %d\n",h_max, c_max);
    fprintf(stderr, "Minumum value from CPU %d from CPU2 %d\n",h_min, c_min);
    fprintf(stderr, "Total value from CPU %d from CPU2 %d\n",h_sum, c_sum);
    
    fprintf(stderr, "Maximum abosulte value from CPU %d from CPU2 %d\n",h_LInf, c_LInf);
    fprintf(stderr, "Total square value from CPU %d from CPU2 %d\n",h_sum2, c_sum2);
    fprintf(stderr, "Dot product value from CPU %d from CPU2 %d\n",h_dot, c_dot);
    fprintf(stderr, "Max Min value from CPU %d %d from CPU2 %d %d\n",h_max, h_min,
            c_o[0], c_o[1]);
    
    //Extensive test
    h_max = -INT_MAX, h_min = INT_MAX;
    h_LInf= 0;
    h_sum2= 0;
    h_sum = 0;
    h_dot = 0;
    h_L1  = 0;
    
    for (int l=1; l < 10001;++l){
        h_max = std::max(h_max, h_i[l-1]);
        h_LInf = std::max(h_LInf, h_i[l-1]);
        h_min = std::min(h_min, h_i[l-1]);
        h_sum += h_i[l-1];
        h_sum2 += h_i[l-1]*h_i[l-1];
        h_dot += h_i1[l-1] * h_i[l-1];
        h_L1  += abs(h_i[l-1]);

        this->Sum(c_sum, h_i, l);
        this->Max(c_max, h_i, l);
        this->Min(c_min, h_i, l);
        this->LInf(c_LInf,h_i,l);
        this->L1(c_L1, h_i, l);
        this->Sum2(c_sum2,h_i, l);
        this->Dot(c_dot,h_i,h_i1, l);
        this->MaxMin(*(Vec2D<int>*)c_o, h_i, l);
        
        if (c_max != h_max){
            fprintf(stderr, "Max Test FAILED at %d CPU2 %d CPU %d\n", l, c_max, h_max );
            test = false;
        }

        if (c_min != h_min){
            fprintf(stderr, "Min Test FAILED at %d CPU2 %d CPU %d\n", l, c_min, h_min );
            test = false;
        }

        if (c_LInf!= h_LInf){
            fprintf(stderr, "MaxAbs Test FAILED at %d CPU2 %d CPU %d\n", l, c_LInf, h_LInf );
            test = false;
        }

        if (c_sum!= h_sum){
            fprintf(stderr, "Sum Test FAILED at %d CPU2 %d CPU %d\n", l, c_sum, h_sum );
            test = false;
        }

        if (c_sum2!= h_sum2){
            fprintf(stderr, "Sum SQR Test FAILED at %d CPU2 %d CPU %d\n", l, c_sum2, h_sum2 );
            test = false;
        }

        if (c_dot!= h_dot){
            fprintf(stderr, "Dot Test FAILED at %d CPU2 %d CPU %d\n", l, c_dot, h_dot );
            test = false;
        }
        
        if ( c_o[0] != h_max || c_o[1] != h_min){
            fprintf(stderr, "Max Min Test FAILED at %d CPU2 %d %d CPU %d %d\n",
                    l, c_o[0], c_o[1], h_max, h_min);
            test = false;
        }
        if (test == false)
            break;
    }

    if (test)
        fprintf(stderr, "Test PASSED  \n");
    
    delete []h_i1;
    delete []h_i;
    return test;
}


} // end namespace PyCA

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

#include <FluidKernelFFTBase.h>

#include <FOpers.h>
#include <CudaUtils.h>
#include <boost_mem.h>

// MSVC lacks a few constants and math functions
#ifdef _MSC_VER
#define M_PI 3.141592653589793
#endif

namespace PyCA {

// ================================================================
// FFTLookupTable3D
// ================================================================

template<class T>
void 
FFTLookupTable3D<T>::
InitTable(){

    //
    // Test that MAX_FFT_TABLE_LENGTH is big enough
    //
    if(mSize.x > MAX_FFT_TABLE_LENGTH || 
       mSize.y > MAX_FFT_TABLE_LENGTH || 
       mSize.z > MAX_FFT_TABLE_LENGTH)
    {
       throw PyCAException(__FILE__,__LINE__,"Error, not enough static memory allocated!");
    }

    //
    // precompute some values
    //
    double sX = 2.0 * M_PI / mSize.x; 
    double sY = 2.0 * M_PI / mSize.y; 
    double sZ = 2.0 * M_PI / mSize.z; 

    double deltaXSq = mSp.x * mSp.x;
    double deltaYSq = mSp.y * mSp.y;
    double deltaZSq = mSp.z * mSp.z;

    //
    // fill in luts
    //
    for (size_t x = 0; x < (size_t)mSize.x; ++x) {
        mCosWX[x] = (2.0 * cos(sX * static_cast<float>(x)) - 2.0) / deltaXSq;
        mSinWX[x] = sin(sX * static_cast<float>(x)) / mSp.x;
    }

    for (size_t y = 0; y < (size_t)mSize.y; ++y) {
        mCosWY[y] = (2.0 * cos(sY * static_cast<float>(y)) - 2.0) / deltaYSq;
        mSinWY[y] = sin(sY * static_cast<float>(y)) / mSp.y;
    }

    for (size_t z = 0; z < (size_t)mSize.z; ++z) {
        mCosWZ[z] = (2.0 * cos(sZ * static_cast<float>(z)) - 2.0) / deltaZSq;
        mSinWZ[z] = sin(sZ * static_cast<float>(z)) / mSp.z;
    }
}
    
template<class T>
void 
FFTLookupTable3D<T>::
Allocate(){
    mCosWX = new T [MAX_FFT_TABLE_LENGTH];
    mCosWY = new T [MAX_FFT_TABLE_LENGTH];
    mCosWZ = new T [MAX_FFT_TABLE_LENGTH];
  
    mSinWX = new T [MAX_FFT_TABLE_LENGTH];
    mSinWY = new T [MAX_FFT_TABLE_LENGTH];
    mSinWZ = new T [MAX_FFT_TABLE_LENGTH];
}

template<class T>
void 
FFTLookupTable3D<T>::
setSize(const Vec3Di& size, 
	const Vec3Df& spacing)
{
    mSize = size;
    mSp   = spacing;
    InitTable();
}

template<class T>
void 
FFTLookupTable3D<T>::
Clear(){
    delete []mCosWX;
    delete []mCosWY;
    delete []mCosWZ;

    delete []mSinWX;
    delete []mSinWY;
    delete []mSinWZ;
}

// ================================================================
// FluidKernelFFTBase
// ================================================================

template<int mode, MemoryType mType, class T>
FluidKernelFFTBase<mode, mType, T>::
FluidKernelFFTBase()
  :mAlpha(1.0), mBeta(0.0), mGamma(0.0), mLPow(1.0),
   mDivergenceFree(false),
   mHasFFTPlan(false),
   mSize(0,0,0), mSpacing(0,0,0),
   mComplexSize(0),
   mAllocateSize(0),
   mFFTArrayX(NULL),
   mFFTArrayY(NULL),
   mFFTArrayZ(NULL)
{
    
}


template<int mode, MemoryType mType, class T>
FluidKernelFFTBase<mode, mType, T>::
~FluidKernelFFTBase()
{
}

template<int mode, MemoryType mType, class T>
void 
FluidKernelFFTBase<mode, mType, T>::
setGrid(const GridInfo& grid)
{
    this->setSize(grid);
}

template<int mode, MemoryType mType, class T>
void
FluidKernelFFTBase<mode, mType, T>::setSize(const GridInfo& g)
{
  
    Vec3Di newSize   = g.size();
    Vec3Df newSpacing= g.spacing();

    bool changeSize    = (mSize    != newSize);
    bool changeSpacing = (mSpacing != newSpacing);

    // TEST
    if(mType != MEM_HOST){
	CudaUtils::CheckCUDAError(__FILE__, __LINE__);
    }
    // END TEST

    if (changeSize) {
        mSize = newSize;

        mComplexSize   = mSize;
        mComplexSize.x = mSize.x/2+1;

        size_t n = mComplexSize.prod();
        if ( n > mAllocateSize)
            Alloc(n);
    }

    if (changeSpacing) {
        mSpacing  = newSpacing;
    }
    
    // recompute the lookup table
    if (changeSize || changeSpacing)
        mLookupTable.setSize(mSize, mSpacing);

    // TEST
    if(mType != MEM_HOST){
	CudaUtils::CheckCUDAError(__FILE__, __LINE__);
    }
    // END TEST

    // create FFT plan
    if (changeSize) {
        if (mHasFFTPlan)
            DestroyFFTPlan();
        CreateFFTPlan();
    }
}

template<int mode, MemoryType mType, class T>
void 
FluidKernelFFTBase<mode, mType, T>::
Alloc(size_t n){
    mFFTArray = CreateSharedArray<mType, ComplexT<T> >(3*n);
    mFFTArrayX = &(mFFTArray.get()[0]);
    mFFTArrayY = &(mFFTArray.get()[n]);
    mFFTArrayZ = &(mFFTArray.get()[2*n]);
    mAllocateSize = n;
}

template<int mode, MemoryType mType, class T>
void 
FluidKernelFFTBase<mode, mType, T>::
ApplyImpl(T* dataX, T* dataY, T* dataZ, bool inverseOp, StreamT stream)
{
    //1. Compute FFT of the data, store result in mFFTArray buffers
    // convert the input from real to frequency(complex image)
    toFrequencyDomain(dataX, dataY, dataZ, stream);

    //2. Solve system
    if (mLPow == floor(mLPow)) {
        size_t nPow = (size_t) mLPow;
        for (size_t i=0; i < nPow; ++i)
            frequencyDomainApply(inverseOp, stream);
    }

    //3. convert the output back to time domain
    toSpaceDomain(dataX, dataY, dataZ, stream);

};

template<int mode, MemoryType mType, class T>
void 
FluidKernelFFTBase<mode, mType, T>::
Apply(Field3D& out, const Field3D& in, 
      bool inverseOp, StreamT stream)
{
    MK_CHECK_REFSIZE_2(mSize, out, in);
    //size_t nElems = mSize.prod();
    //FieldOpers<mode>::MulC(out, in, 1.f / nElems, stream);
    Opers::Copy(out, in, NULL);
    ApplyImpl(out.x, out.y, out.z, inverseOp, stream);
}

template<int mode, MemoryType mType, class T>
void 
FluidKernelFFTBase<mode, mType, T>::
Apply(Field3D& out, bool inverseOp, StreamT stream){
    MK_CHECK_REFSIZE_1(mSize, out);
    //size_t nElems = mSize.prod();
    //FieldOpers<mode>::MulC_I(out, 1.f / nElems, stream);
    ApplyImpl(out.x, out.y, out.z, inverseOp, stream);
}

//--------------------------------------------------------------------------------
/**
 * f = Lv
 * 
 * v field is overwritten in this operation
 */
template<int mode, MemoryType mType, class T>
void 
FluidKernelFFTBase<mode, mType, T>::
applyOperator(Field3D& f, const Field3D& v, StreamT stream){
    Apply(f, v, false, stream);
}

/**
 * v = Kf
 * 
 * f field is overwritten in this operation
 */
template<int mode, MemoryType mType, class T>
void 
FluidKernelFFTBase<mode, mType, T>::
applyInverseOperator(Field3D& v, const Field3D& f, StreamT stream){
    Apply(v, f, true, stream);
}
    

/**
 * f = Lv
 * 
 * v field is overwritten in this operation (holds f).
 */
template<int mode, MemoryType mType, class T>
void 
FluidKernelFFTBase<mode, mType, T>::
applyOperator(Field3D& f, StreamT stream){
    Apply(f, false, stream);
}

/**
 * v = Kf
 * 
 * f field is overwritten in this operation (holds v).
 */
template<int mode, MemoryType mType, class T>
void 
FluidKernelFFTBase<mode, mType, T>::
applyInverseOperator(Field3D& v, StreamT stream){
    Apply(v, true, stream);
}

// template instantiations
template class FluidKernelFFTBase<EXEC_CPU, MEM_HOST, float>;
#ifdef CUDA_ENABLED
template class FluidKernelFFTBase<EXEC_GPU, MEM_DEVICE, float>;
#endif

}; // end namespace PyCA

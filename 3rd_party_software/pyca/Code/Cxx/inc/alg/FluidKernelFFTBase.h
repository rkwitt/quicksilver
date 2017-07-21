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

#ifndef __FLUID_KERNEL_FFT_BASE_H__
#define __FLUID_KERNEL_FFT_BASE_H__

#include <FluidKernelFFTInterface.h>

// if REPLACE_ZERO_FREQ is defined, FluidKernelFFT will replace the
// zero frequency when gamma is zero (ie L is non-invertible),
// otherwise the zero frequency is set to zero
// #define REPLACE_ZERO_FREQ

namespace PyCA {

// forward declaration
class Field3D;

// ================================================================
// FFTLookupTable3D
// ================================================================

template<class T>
class FFTLookupTable3D {
public:
    FFTLookupTable3D():mCosWX(NULL), mCosWY(NULL), mCosWZ(NULL),
		       mSinWX(NULL), mSinWY(NULL), mSinWZ(NULL)
    {
	Allocate();
    };
    
    ~FFTLookupTable3D(){
        Clear();
    }

    void setSize(const Vec3Di& size, 
		 const Vec3Df& sp=Vec3Df(1.f, 1.f, 1.f));

    const T* CosWX() const { return mCosWX; }
    const T* CosWY() const { return mCosWY; }
    const T* CosWZ() const { return mCosWZ; }
    const T* SinWX() const { return mSinWX; }
    const T* SinWY() const { return mSinWY; }
    const T* SinWZ() const { return mSinWZ; }
    Vec3Di size() const { return mSize; }
    Vec3Df spacing() const { return mSp; }
    
protected:
    
    void Allocate();
    void InitTable();
    void Clear();
    
    // host data 
    T* mCosWX,* mCosWY,* mCosWZ;
    T* mSinWX,* mSinWY,* mSinWZ;
    
    Vec3Di mSize;
    Vec3Df mSp;
}; // end FFTLookupTable3D

// ================================================================
// FluidKernelFFTBase
// ================================================================

template<int mode, MemoryType mType, class T>
class FluidKernelFFTBase :
    public FluidKernelFFTInterface
{
public:

    FluidKernelFFTBase();
    ~FluidKernelFFTBase();

   // Implementation of FluidKernelFFTInterface

   void setAlpha(float alpha){ mAlpha = alpha; }
   float getAlpha(){ return mAlpha; }
   void setBeta(float beta){ mBeta = beta; }
   float getBeta(){ return mBeta; }
   void setGamma(float gamma){ mGamma = gamma; }
   float getGamma(){ return mGamma; }
   void setLPow(int lpow){ mLPow = lpow; }
   int getLPow(){ return mLPow; }
   void setDivergenceFree(bool divfree){ mDivergenceFree = divfree; }
   bool getDivergenceFree(){ return mDivergenceFree; }

   void setGrid(const GridInfo& grid);

   MemoryType memType() const { return mType; }
    

    /**
     * Getters for accessing internal FFT arrays, just for c++
     * interface, not accessible from wrapped code
     */
    ComplexT<T>* getFFTArrayX(){ return mFFTArrayX; }
    ComplexT<T>* getFFTArrayY(){ return mFFTArrayY; }
    ComplexT<T>* getFFTArrayZ(){ return mFFTArrayZ; }

    
    /**
     * f = Lv
     * 
     * v field is overwritten in this operation
     */
    void applyOperator(Field3D& f, const Field3D& v, 
		       StreamT stream);

    /**
     * v = Kf
     * 
     * f field is overwritten in this operation
     */
    void applyInverseOperator(Field3D& v, const Field3D& f, 
			      StreamT stream);
    
    /**
     * f = Lv
     * 
     * v field is overwritten in this operation (holds f).
     */
    void applyOperator(Field3D& f, StreamT stream);

    /**
     * v = Kf
     * 
     * f field is overwritten in this operation (holds v).
     */
    void applyInverseOperator(Field3D& v, StreamT stream);

protected:

    virtual void setSize(const GridInfo& g);

    virtual void Apply(Field3D& out, const Field3D& in, 
		       bool inverseOp, StreamT stream);
    virtual void Apply(Field3D& out, bool inverseOp, StreamT stream);
    virtual void ApplyImpl(T* dataX, T* dataY, T* dataZ, 
			   bool inverseOp, StreamT stream);
            
    virtual void CreateFFTPlan() = 0;
    virtual void DestroyFFTPlan() = 0;
    
    /**
     * Perform FFT on dX, dY, dZ, placing the results in the internal
     * complex arrays mFFTArrayX, mFFTArrayY, mFFTArrayZ.
     */
    virtual void toFrequencyDomain(const T* x, const T* y, const T* z,
				   StreamT s = NULL) = 0;
    /**
     * compute inverse FFT on the internal complex arrays mFFTArrayX,
     * mFFTArrayY, mFFTArrayZ, place the result in the real arrays dX,
     * dY, dZ.
     */
    virtual void toSpaceDomain(T* x, T* y, T* z,
			      StreamT s = NULL) = 0;

    /**
     * Apply the forward or inverse differential operator to the
     * frequency domain data in the internal complex arrays
     * mFFTArrayX, mFFTArrayY, mFFTArrayZ.
     */
    virtual void frequencyDomainApply(bool inverseOp, 
				      StreamT s = NULL) = 0;
    
    virtual void Alloc(size_t n);
    
    T mAlpha, mBeta, mGamma;
    int mLPow;
    bool mDivergenceFree;

    bool mHasFFTPlan;
    
    // ===== Member Data =====
    Vec3Di mSize;
    Vec3Df mSpacing;
    Vec3Di mComplexSize;
    
    // Lookup table
    FFTLookupTable3D<T> mLookupTable;
  
    unsigned int mAllocateSize;
    
    // FFT array pointers
    boost::shared_ptr<ComplexT<T> > mFFTArray;
    // These point to locations in mFFTArray
    ComplexT<T>* mFFTArrayX;
    ComplexT<T>* mFFTArrayY;
    ComplexT<T>* mFFTArrayZ;

}; // end FluidKernelFFTBase

} // end namespace PyCA

#endif // __FLUID_KERNEL_FFT_BASE_H__

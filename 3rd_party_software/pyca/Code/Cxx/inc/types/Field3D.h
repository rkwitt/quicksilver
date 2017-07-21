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

#ifndef __FIELD3D_H
#define __FIELD3D_H

#ifndef SWIG
#include <MemPool.h>
#include <mem.h>
#include <GridInfo.h>
#endif // SWIG


namespace PyCA {

class Field3D{
private:
    // Will not allow copy to write 
    Field3D(const Field3D& other);
    Field3D& operator=(const Field3D& other);
    void init();
protected:
    MemoryType     mType;
    GridInfo       mGrid;
    size_t         mCap;
    MemPool<float> mDataX, mDataY, mDataZ;
public:
    float *x, *y, *z;
    Field3D(MemoryType type);
    Field3D(const GridInfo& g, MemoryType type);
    Field3D(const GridInfo& g,
            const MemPool<float>& dataX,
            const MemPool<float>& dataY,
            const MemPool<float>& dataZ);

    // convenience constructor (usually for python)
    Field3D(int sx, int sy, int sz, MemoryType type=MEM_HOST);

    virtual ~Field3D();

    void print(std::ostream& oss) const;
    // get memory type
    MemoryType memType() const { return mType; };
    
    // get the basic properties
    const GridInfo& grid() const { return mGrid; };
   void setGrid(const GridInfo& g);

    // get the raw pointer to data
    const float* getX() const { return mDataX.get(); };
    float* getX()             { return mDataX.get(); };
    const float* getY() const { return mDataY.get(); };
    float* getY()             { return mDataY.get(); };
    const float* getZ() const { return mDataZ.get(); };
    float* getZ()             { return mDataZ.get(); };

    // ================================================================
    // convenience functions, only work for CPU objects.  Very
    // inefficient, don't use in production code.  Used in python
    // wrapping

    void set(size_t xIndex, 
	     size_t yIndex,
	     size_t zIndex, 
	     const Vec3Df& item)
   {
      PYCA_ASSERT(mType == MEM_HOST);
      size_t idx = mGrid.size().x*(mGrid.size().y*zIndex + yIndex) + xIndex;
      this->x[idx] = item.x;
      this->y[idx] = item.y;
      this->z[idx] = item.z;
   }
  
    Vec3Df get(size_t xIndex, 
	       size_t yIndex,
	       size_t zIndex) const
   {
      PYCA_ASSERT(mType == MEM_HOST);
      size_t idx = mGrid.size().x*(mGrid.size().y*zIndex + yIndex) + xIndex;
      Vec3Df v;
      v.x = this->x[idx];
      v.y = this->y[idx];
      v.z = this->z[idx];
      return v;
   }
  
    Vec3Df operator()(size_t xIndex, 
		      size_t yIndex,
		      size_t zIndex) const
    {
      return this->get(xIndex, yIndex, zIndex);
    }

    void set(size_t index, 
	     const Vec3Df& item)
   {
      PYCA_ASSERT(mType == MEM_HOST);
      this->x[index] = item.x;
      this->y[index] = item.y;
      this->z[index] = item.z;
   }
  
   Vec3Df get(size_t index) const
   {
      PYCA_ASSERT(mType == MEM_HOST);
      Vec3Df v;
      v.x = this->x[index];
      v.y = this->y[index];
      v.z = this->z[index];
      return v;
   }
  
    Vec3Df operator()(size_t index) const
    {
      return this->get(index);
    }

    // ================================================================

    const MemPool<float>& getMemPoolX() const { return mDataX; };
    MemPool<float>&       getMemPoolX()       { return mDataX; };
    const MemPool<float>& getMemPoolY() const { return mDataY; };
    MemPool<float>&       getMemPoolY()       { return mDataY; };
    const MemPool<float>& getMemPoolZ() const { return mDataZ; };
    MemPool<float>&       getMemPoolZ()       { return mDataZ; };

    void swap(Field3D& rhs);
    
    // convenient functions
    const Vec3Di& size() const { return mGrid.size(); };
    const Vec3Df& spacing() const { return mGrid.spacing(); };
    const Vec3Df& origin() const { return mGrid.origin(); };

   void clone(Field3D& other);
   /**
    * Set the given grid, resizing the field's capacity to match
    */
   void resize(const GridInfo &grid, bool preserveData = false, StreamT stream = NULL);
   
   /**
    * change the capacity of this field
    */
   void setCap(const Vec3Di &size, bool preserveData = false, StreamT stream = NULL);
   void setCap(size_t nVox, bool preserveData = false, StreamT stream = NULL);

   /**
    * Change the memory type of this field, preserving contents
    */ 
   void toType(MemoryType memType, StreamT stream = NULL);

    size_t nVox()      const { return mGrid.nVox(); };
    float  voxelVol()  const { return mGrid.voxelVol(); };
    size_t capacity()  const { return mCap; };
    
    // Is the memory for the three components actually contiguous?
    bool isContinuous() const { return y == x+mGrid.nVox()
                                    && z == y+mGrid.nVox(); };
    bool onHost() const { return mType != MEM_DEVICE; }
};

} // end namespace PyCA

#endif

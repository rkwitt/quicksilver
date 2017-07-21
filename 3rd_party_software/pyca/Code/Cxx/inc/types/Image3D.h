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

#ifndef __IMAGE3D_H
#define __IMAGE3D_H

#ifndef SWIG
#include <MemPool.h>
#include <mem.h>
#include <GridInfo.h>
#endif // SWIG

namespace PyCA {

class Image3D{
private:
    // Will not allow copy to write 
    Image3D();
    Image3D(const Image3D& other);
    Image3D& operator=(const Image3D& other);
public:
    Image3D(MemoryType type);
    Image3D(const GridInfo& g, MemoryType type);
    Image3D(const GridInfo& g, const MemPool<float>& data);

    // convenience constructor (usually for python)
    Image3D(int sx, int sy, int sz, MemoryType type=MEM_HOST);

    // Virtual destructor, overridden in ManagedImage3D
    virtual ~Image3D();

    void print(std::ostream& oss) const;

    // get memory type
    MemoryType memType() const { return mType; };
    bool onDev() const { return mType == MEM_DEVICE; };
    bool onHost() const { return !onDev(); };
    
    // get the basic properties
    const GridInfo& grid() const { return mGrid; };
   // setGrid does not resize data -- will error if grid is greater
   // than capacity.  use resize to change allocated array size.
    void setGrid(const GridInfo& g);

    // get the raw pointer to data
    const float* get() const { return mData.get(); };
    float* get()             { return mData.get(); };

    const MemPool<float>& getMemPool() const { return mData; };
    MemPool<float>&       getMemPool()       { return mData; };

    void swap(Image3D& rhs);

   // ================================================================
   // convenience functions, only work for CPU objects.  Used in
   // python wrapping.

   const float& get(size_t idx) const {
      PYCA_ASSERT(mType != MEM_DEVICE); // Should we allow MEM_HOST_PINNED?
      return mData[idx];
   }
    
   float& get(size_t idx) {
      PYCA_ASSERT(mType != MEM_DEVICE); // Should we allow MEM_HOST_PINNED?
      return mData[idx];
   }

   const float& get(size_t xIndex, 
	     size_t yIndex,
	     size_t zIndex) const
   {
      PYCA_ASSERT(mType != MEM_DEVICE); // Should we allow MEM_HOST_PINNED?
      size_t idx = mGrid.size().x*(mGrid.size().y*zIndex + yIndex) + xIndex;
      return mData[idx];
   }

   float& get(size_t xIndex, 
	      size_t yIndex,
	      size_t zIndex)
   {
      PYCA_ASSERT(mType != MEM_DEVICE); // Should we allow MEM_HOST_PINNED?
      size_t idx = mGrid.size().x*(mGrid.size().y*zIndex + yIndex) + xIndex;
      return mData[idx];
   }

#ifndef SWIG   
   const float& operator[] (size_t idx) const {
      return this->get(idx);
   }
    
   float& operator[] (size_t idx) {
      return this->get(idx);
   }
#endif // ! SWIG   

   const float& operator()(size_t idx) const
   {
      return this->get(idx);
   }
   
   float& operator()(size_t idx)
   {
      return this->get(idx);
   }

   const float& operator()(size_t xIndex, 
			   size_t yIndex,
			   size_t zIndex) const
   {
      return this->get(xIndex, yIndex, zIndex);
   }
   
   float& operator()(size_t xIndex, 
		     size_t yIndex,
		     size_t zIndex)
   {
      return this->get(xIndex, yIndex, zIndex);
   }

   void set(size_t idx, float item) {
      this->get(idx) = item;
   }
    
   void set(size_t xIndex,
	    size_t yIndex,
	    size_t zIndex, 
	    float item) 
   {
      this->get(xIndex, yIndex, zIndex) = item;
   }
    
   // ================================================================

    const Vec3Di& size() const { return mGrid.size(); };
    const Vec3Df& spacing() const { return mGrid.spacing(); };
    const Vec3Df& origin() const { return mGrid.origin(); };

    void setSize(const Vec3Di &size);
    void setSpacing(const Vec3Df &sp){ mGrid.setSpacing(sp); }
    void setOrigin(const Vec3Df &org){ mGrid.setOrigin(org); }

    
   void clone(const Image3D &other);

   /**
    * Set the given grid, resizing the field's capacity to match
    */
   void resize(const GridInfo &grid, bool preserveData = false, StreamT stream = NULL);
   
   /**
    * change the capacity of this field
    */
   void setCap(const Vec3Di &size, bool preserveData = false, StreamT stream = NULL);
   void setCap(size_t nVox, bool preserveData = false, StreamT stream = NULL);

   void toType(MemoryType memType, StreamT stream = NULL);

    size_t nVox()     const { return mGrid.nVox(); };
    float  voxelVol() const { return mGrid.voxelVol(); };
    size_t capacity() const { return mData.capacity(); };
    
protected:
    MemoryType     mType;
    GridInfo       mGrid;
    MemPool<float> mData;
};

} // end namespace PyCA

#endif

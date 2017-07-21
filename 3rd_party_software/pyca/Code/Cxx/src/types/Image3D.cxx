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

#include <Image3D.h>
#include <GridInfo.h>

namespace PyCA {

Image3D::Image3D(MemoryType type):mType(type), mGrid(), mData(type){
    
}

Image3D::Image3D(const GridInfo& g, MemoryType type):
    mType(type), mGrid(g), mData(g.nVox(), type)
{
    
}

Image3D::Image3D(const GridInfo& g, const MemPool<float>& data):
    mType(data.memType()), mGrid(g), mData(data)
{
    // Check to see if the size of the image is less
    // than the memory available
   PRECONDITION(g.nVox() <= (int)mData.capacity(),"Data too large for given MemPool");
}

Image3D::~Image3D()
{
   mType = MEM_UNINITIALIZED;
}

Image3D::Image3D(int sx, int sy, int sz, MemoryType type)
    : mType(type), 
      mGrid(Vec3Di(sx,sy,sz)), 
      mData(mGrid.nVox(), type)
{

}

void Image3D::print(std::ostream& oss) const{
    oss << "Grid " << mGrid << std::endl;
    oss << "Number of elements " << mGrid.nVox() << " Capacity " << mData.capacity() << std::endl;
    oss << "Address " << this->get() << std::endl;
}

void Image3D::setGrid(const GridInfo& g) {
   PRECONDITION(g.nVox() <= (int)mData.capacity(),"grid too large for current image");
    mGrid = g;
}

void Image3D::setSize(const Vec3Di &size){ 
   mGrid.setSize(size);  
   PRECONDITION(mGrid.nVox() <= (int)mData.capacity(),"size too large for current image");
}

void Image3D::swap(Image3D& rhs){
   MK_CHECK2_MEM((*this),rhs);
   
   std::swap(mGrid, rhs.mGrid);
   mData.swap(rhs.mData);
}

void 
Image3D::
clone(const Image3D &other)
{
   mData.clone(other.mData);
   mGrid = other.mGrid;
   mType = other.mType;
}

void 
Image3D::
resize(const GridInfo &grid, bool preserveData, StreamT stream)
{
   this->setCap(grid.size(), preserveData, stream);
   this->setGrid(grid);
}

void 
Image3D::
setCap(const Vec3Di &size, bool preserveData, StreamT stream)
{
   size_t nVox = size.prod();
   setCap(nVox, preserveData, stream);
}

void 
Image3D::
setCap(size_t nVox, bool preserveData, StreamT stream)
{
   if(mGrid.nVox() == (int)nVox) return;
   mData.resize(nVox, preserveData, stream);
}

void 
Image3D::
toType(MemoryType memType, StreamT stream)
{
    if(memType == mType) return;

    mData.toType(memType, stream);
    mType = memType;
}

} // end namespace PyCA


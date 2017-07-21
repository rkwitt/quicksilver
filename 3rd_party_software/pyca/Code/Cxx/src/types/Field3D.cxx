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

#include <Field3D.h>
#include <pycaUtils.h>
#include <pycaConst.h>

namespace PyCA {

Field3D::Field3D(MemoryType type):
   mType(type), 
   mGrid(), 
   mCap(0),
   mDataX(type),
   mDataY(type),
   mDataZ(type)
{
    
}

Field3D::Field3D(const GridInfo& g, MemoryType type):
    mType(type),
    mGrid(g),
    mCap(g.nVox()),
    mDataX(g.nVox(), type),
    mDataY(g.nVox(), type),
    mDataZ(g.nVox(), type)
{
    init();
}

Field3D::Field3D(const GridInfo& g,
                 const MemPool<float>& dataX,
                 const MemPool<float>& dataY,
                 const MemPool<float>& dataZ):
    mType(dataX.memType()),
    mGrid(g),
    mCap(g.nVox()),
    mDataX(dataX),
    mDataY(dataY),
    mDataZ(dataZ)
{
    PRECONDITION(dataX.capacity() == dataY.capacity() 
                 && dataY.capacity() == dataZ.capacity(),
                 "Provided mempools are different sizes");
    PRECONDITION(dataX.capacity() >= g.nVox(),"Grid is larger than pre-allocated memory");
    init();
}

Field3D::Field3D(int sx, int sy, int sz, MemoryType type)
    : mType(type), 
      mGrid(Vec3Di(sx,sy,sz)), 
      mCap(mGrid.nVox()),
      mDataX(mGrid.nVox(), type),
      mDataY(mGrid.nVox(), type),
      mDataZ(mGrid.nVox(), type)
{
    init();
}

Field3D::~Field3D() {
    mType = MEM_UNINITIALIZED;
}

void Field3D::print(std::ostream& oss) const{
    oss << "Grid " << mGrid << std::endl;
    oss << "Number of elements " << mGrid.nVox() << std::endl;
    oss << "Address X " << this->x << std::endl;
    oss << "Address Y " << this->y << std::endl;
    oss << "Address Z " << this->z << std::endl;
}

void Field3D::setGrid(const GridInfo& g) {
    mGrid   = g;
    PRECONDITION(mDataX.capacity() >= g.nVox(),"Grid is larger than pre-allocated memory");
    init();
}

void Field3D::init(){
    x = mDataX.get();
    y = mDataY.get();
    z = mDataZ.get();
}

void Field3D::swap(Field3D& rhs){
    PYCA_ASSERT(mType == rhs.mType);
    std::swap(mGrid  , rhs.mGrid);
    std::swap(mCap, rhs.mCap);
    std::swap(x, rhs.x);
    std::swap(y, rhs.y);
    std::swap(z, rhs.z);
    mDataX.swap(rhs.mDataX);
    mDataY.swap(rhs.mDataY);
    mDataZ.swap(rhs.mDataZ);
}

void Field3D::clone(Field3D& other)
{
   mDataX.clone(other.mDataX);
   mDataY.clone(other.mDataY);
   mDataZ.clone(other.mDataZ);
   mType = other.mType;
   mGrid = other.mGrid;
   mCap = other.mCap;
   init();
}

void 
Field3D::
resize(const GridInfo &grid, bool preserveData, StreamT stream)
{
   this->setCap(grid.size(), preserveData, stream);
   this->setGrid(grid);
}

void 
Field3D::
setCap(const Vec3Di &size, bool preserveData, StreamT stream)
{
   size_t nVox = size.prod();
   setCap(nVox, preserveData, stream);
}

void 
Field3D::
setCap(size_t nVox, bool preserveData, StreamT stream)
{
   if(mGrid.nVox() == (int)nVox) return;
   size_t newCap = nVox;
   mDataX.resize(newCap, preserveData, stream);
   mDataY.resize(newCap, preserveData, stream);
   mDataZ.resize(newCap, preserveData, stream);
}


void 
Field3D::
toType(MemoryType memType, StreamT stream)
{
    if(memType == mType) return;
    
    mDataX.toType(memType, stream);
    mDataY.toType(memType, stream);
    mDataZ.toType(memType, stream);
    mType = memType;
    init();
}

} // end namespace PyCA

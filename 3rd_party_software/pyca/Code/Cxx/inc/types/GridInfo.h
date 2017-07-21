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

#ifndef __GRID_H
#define __GRID_H

#include <Vec3D.h>

#define GRIDINFO_EPS 1e-6

namespace PyCA {

class GridInfo{
public:
    

    typedef Vec3Di size_type;
    typedef Vec3Df spacing_type;
    typedef Vec3Df point_type;
    
    GridInfo(): mSize(0,0,0), mSp(1.f,1.f,1.f), mOrg(0.f, 0.f, 0.f){};
    explicit GridInfo(const size_type& size,
                      const spacing_type& sp = spacing_type(1.f, 1.f, 1.f),
                      const point_type& org= point_type(0.f, 0.f, 0.f))
        :mSize(size), mSp(sp), mOrg(org)
        {
        }

   GridInfo(const GridInfo &grid)
      :mSize(grid.size()), mSp(grid.spacing()), mOrg(grid.origin())
   {
   }
    
    void set(const size_type& nSize,
             const spacing_type& nSp  = spacing_type(1.f, 1.f, 1.f),
             const point_type& nOrg = point_type(0.f, 0.f, 0.f))
        {
            mSize = nSize;
            mSp   = nSp;
            mOrg  = nOrg;
        }
    
    const size_type& size()    const { return mSize; };
    const spacing_type& spacing() const { return mSp;   };
    const point_type& origin()  const { return mOrg;  };

    void setSize(const size_type& size)
    {
	mSize = size;
    }
    
    void setSpacing(const spacing_type& sp)
    {
	mSp = sp;
    }
    
    void setOrigin(const point_type& org)
    {
	mOrg = org;
    }
    
    void read(std::istream& iss) {
        std::string temp;
        iss >> temp >> mSize >> temp;
        iss >> temp >> mSp >> temp;
        iss >> temp >> mOrg;
    }

#ifndef SWIG
   // types.i defines copy() and clone(...) in stead
   GridInfo& operator=(const GridInfo& rhs) {
      mSize = rhs.size();
      mSp = rhs.spacing();
      mOrg = rhs.origin();
      return *this;
   }
#endif // !SWIG

    bool operator==(const GridInfo& rhs) const{
        return ((mSize == rhs.mSize) && 
		((mSp-rhs.mSp).length() < GRIDINFO_EPS) && 
		((mOrg-rhs.mOrg).length() < GRIDINFO_EPS));
    }
        
    bool operator!=(const GridInfo& rhs) const{
        return !(*this == rhs);
    }
    
    void print(std::ostream& oss) const {
        oss << "Size= "    << mSize    <<", ";
        oss << "Spacing= " << mSp << ", ";
        oss << "Origin= "  << mOrg;
    }

    size_t    nVox() const { return mSize.prod(); };
    float voxelVol() const { return (float)mSp.prod(); };
private:

    // basic grid information
    size_type    mSize;
    spacing_type mSp;
    point_type   mOrg;
};

inline std::ostream& operator<<(std::ostream& oss, const GridInfo& p){
    p.print(oss);
    return oss;
}

inline std::istream& operator>>(std::istream& iss, GridInfo& p){
    p.read(iss);
    return iss;
}

} // end namespace PyCA

#endif

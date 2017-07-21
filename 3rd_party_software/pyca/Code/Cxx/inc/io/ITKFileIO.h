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

#ifndef __ITK_FILEIO_H
#define __ITK_FILEIO_H

#ifndef SWIG
#include "itkImageIOBase.h"
#include "itkImageIOFactory.h"
#endif // SWIG

#include <GridInfo.h>
#include <Aff3D.h>


namespace PyCA {

class Image3D;
class Field3D;

class ITKFileIO{
public:

#ifndef SWIG
    typedef itk::ImageIOBase::IOComponentType  ComponentType;
#else
    // fake typedef for SWIG
    typedef int ComponentType;
#endif

    /**
     * Read the grid and dataType (itk::ImageIOBase::IOComponentType)
     * from the file given by `name`
     */
    static void ReadHeader(GridInfo& gridOut, ComponentType& DATA_TYPE_OUT, const char* name);

    /**
     * Load the image in file `name` into `image`.  The capacity of
     * image will be changed to accommodate the image, and grid will
     * be set with correct values.
     */
    static void LoadImage(Image3D& image, const char* name);

    /**
     * Load the field in file `name` into `h`.  The capacity of h will
     * be changed to accommodate the image, and grid will be set with
     * correct values.
     */
    static void LoadField(Field3D& h, const char* name);

    /**
     * Load a Field3D from three images, with names prefix_x.mha,
     * prefix_y.mha, and prefix_z.mha
     */
    static void LoadFieldComponents(Field3D& h, const char* prefix);

    /**
     * Save image, optionally using compression (if available for the
     * given type).  Type is inferred from the filename extension.
     */
    static void SaveImage(const Image3D& image, const char* name, bool useCompression=false);

    /**
     * Save field, optionally using compression (if available for the
     * given type).  Type is inferred from the filename extension.
     */
    static void SaveField(const Field3D& h, const char* name, bool useCompression=false);

    /**
     * Save a Field3D to three images, with names prefix_x.mha,
     * prefix_y.mha, and prefix_z.mha, and optionally using compression.
     */
    static void SaveFieldComponents(const Field3D& h, const char* prefix, bool useCompression=false);

    /**
     * Load an affine transformation from an ITK AffineTransform file
     */
    template<typename T>
    void LoadAff3D(Aff3D<T> &aff, const char* filename);
};

} // end namespace PyCA

#endif

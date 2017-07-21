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

#include <ITKFileIO.h>

#include <Image3D.h>
#include <Field3D.h>

#include <boost_mem.h>

#include <itkImage.h>
#include <itkImportImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>

#include <itkTransformFileReader.h>
#include <itkTransformBase.h>
#include <itkAffineTransform.h>

#if ITK_VERSION_MAJOR < 4
#include <itkCompose3DCovariantVectorImageFilter.h>
#else
#include <itkComposeImageFilter.h>
#endif

namespace PyCA {

template<typename inT, typename outT>
void CopyFromITKBuffer(const char* fName, outT* oBuf, size_t n){

    typedef itk::Image<inT, 3>              ImageType;
    typedef itk::ImageFileReader<ImageType> ReaderType;

    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(fName);
    typename ImageType::Pointer  imagePtr = reader->GetOutput();
    try
    {
        reader->Update();
    }
    catch(itk::ImageFileReaderException &exc)
    {
        throw PyCAException(__FILE__, __LINE__, "ITK File Reader Exception", &exc);
    }
    
    catch(itk::ExceptionObject &exc)
    {
        throw PyCAException(__FILE__, __LINE__, "ITK Exception", &exc);
    }
    
    catch(...)
    {
        throw PyCAException(__FILE__, __LINE__, "Unknown exception loading");
    }

    std::copy(imagePtr->GetBufferPointer(), imagePtr->GetBufferPointer() + n, oBuf);
}

template<class T>
typename itk::ImportImageFilter<T, 3>::Pointer CreateImportFilter(const GridInfo& grid) {
    typedef itk::Image<T, 3>    ImageType;
    typedef itk::ImportImageFilter<T, 3> ImportFilterType;

    Vec3Di ImSize = grid.size();
    Vec3Df ImOrg  = grid.origin();
    Vec3Df ImSp   = grid.spacing();
    typename ImportFilterType::Pointer  importFilter = ImportFilterType::New();
    
    typename ImageType::PointType origin;
    origin[0] = ImOrg.x; origin[1] = ImOrg.y; origin[2] = ImOrg.z;
    importFilter->SetOrigin(origin);
        
    typename ImageType::SpacingType spacing;
    spacing[0]= ImSp.x; spacing[1] = ImSp.y; spacing[2] = ImSp.z;
    importFilter->SetSpacing(spacing);
    
    typename ImageType::IndexType  start;
    typename ImageType::SizeType   size;
    typename ImageType::RegionType region;
    start[0]  = start[1] = start[2] = 0;
    size[0]   = ImSize.x;  size[1] = ImSize.y; size[2]  = ImSize.z;
    region.SetSize(size);
    region.SetIndex(start);
    importFilter->SetRegion(region);
    
    return importFilter;
}

////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
void ITKFileIO::ReadHeader(GridInfo& grid, ComponentType& dataType, const char* fName)
{
    itk::ImageIOBase::Pointer imageIO =
        itk::ImageIOFactory::CreateImageIO(fName, itk::ImageIOFactory::ReadMode);
    
    if ( !imageIO )  {
       throw PyCAException(__FILE__, __LINE__, std::string("Create ITKImageIO failed: ")+fName);
    }
    imageIO->SetFileName(fName);
    imageIO->ReadImageInformation();
    
    size_t nDims = imageIO->GetNumberOfDimensions();
    if(nDims < 2 || nDims > 3){
        throw PyCAException(__FILE__,__LINE__,
			    "Error, number of dimensions not 2 or 3");
    }
    Vec3Di imSize;
    Vec3Df imSp, imOrg;
    for (size_t i=0; i< nDims; ++i) {
        imSize[i] = imageIO->GetDimensions(i);
        imOrg[i]  = imageIO->GetOrigin(i);
        imSp[i]   = imageIO->GetSpacing(i);
    }
    if(nDims == 2){
       imSize[2] = 1;
       imOrg[2] = 0.f;
       imSp[2] = 1.f;
    }
    grid = GridInfo(imSize, imSp, imOrg);

    dataType = imageIO->GetComponentType();
}

void ITKFileIO::LoadImage(Image3D& image, const char* name)
{
    GridInfo grid;
    ComponentType type;

    std::cerr << "Loading " << name << "...";

    ReadHeader(grid, type, name);
    if(type != itk::ImageIOBase::FLOAT){
       std::cerr << "converting from " 
		 << itk::ImageIOBase::GetComponentTypeAsString(type)
		 << " to float.";
    }
    // if (type != itk::ImageIOBase::FLOAT) {
    //     throw PyCAException(__FILE__,__LINE__,"We do not support type other than float at this time");
    // }

    size_t nVox = grid.nVox();
    Image3D temp(grid, image.memType());
    if (image.memType() != MEM_DEVICE) {
        CopyFromITKBuffer<float, float>(name, temp.get(), nVox);
    } else {
        boost::shared_ptr<float> tempH = CreateSharedArray<MEM_HOST_PINNED,float>(nVox);
        CopyFromITKBuffer<float, float>(name, tempH.get(), nVox);
        cpyArrayH2D(temp.get(), tempH.get(), nVox);
    }
    image.swap(temp);

    std::cerr << "done " << std::endl;
}


void ITKFileIO::LoadField(Field3D& h, const char* name)
{
    GridInfo grid;
    ComponentType type;

    typedef itk::CovariantVector< float, 3> VectorPixelType;
    typedef itk::Image<VectorPixelType, 3> VectorImageType;
    typedef itk::Image<float, 3> ScalarImageType;

    ReadHeader(grid, type, name);

    size_t nVox = grid.nVox();
    h.resize(grid);

    boost::shared_ptr<float> tmp_buff[3];
    float *buff[3];
    if (h.memType() != MEM_DEVICE) {
       buff[0] = h.x;
       buff[1] = h.y;
       buff[2] = h.z;
    } else {
        tmp_buff[0] = CreateSharedArray<MEM_HOST_PINNED,float>(nVox);
	buff[0] = tmp_buff[0].get();
        tmp_buff[1] = CreateSharedArray<MEM_HOST_PINNED,float>(nVox);
	buff[1] = tmp_buff[1].get();
        tmp_buff[2] = CreateSharedArray<MEM_HOST_PINNED,float>(nVox);
	buff[2] = tmp_buff[2].get();
    }

    typedef itk::ImageFileReader<VectorImageType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(name);
    VectorImageType::Pointer  imagePtr = reader->GetOutput();
    try
    {
        reader->Update();
    }
    catch(itk::ImageFileReaderException &exc)
    {
        throw PyCAException(__FILE__, __LINE__, "ITK File Reader Exception", &exc);
    }
    
    catch(itk::ExceptionObject &exc)
    {
        throw PyCAException(__FILE__, __LINE__, "ITK Exception", &exc);
    }
    
    catch(...)
    {
        throw PyCAException(__FILE__, __LINE__, "Unknown exception loading");
    }

    typedef 
       itk::VectorIndexSelectionCastImageFilter<VectorImageType, ScalarImageType> 
       ExtractComponentFilterType;
    ExtractComponentFilterType::Pointer extractFilter = ExtractComponentFilterType::New();
    
    extractFilter->SetInput(reader->GetOutput());

    for(int dim=0; dim<3; ++dim){
       extractFilter->SetIndex(dim);
       extractFilter->Update();
       std::copy(extractFilter->GetOutput()->GetBufferPointer(), 
		 extractFilter->GetOutput()->GetBufferPointer() + nVox,
		 buff[dim]);
    }

    if(h.memType() == MEM_DEVICE){
       cpyArrayH2D(h.x, buff[0], nVox);
       cpyArrayH2D(h.y, buff[1], nVox);
       cpyArrayH2D(h.z, buff[2], nVox);
    }

}

void ITKFileIO::LoadFieldComponents(Field3D& h, const char* prefix)
{
    std::string oPrefix(prefix);
    std::string name;

    GridInfo grid;
    ComponentType type;
    name = oPrefix + "_x.mha";
    ReadHeader(grid, type, name.c_str());
    if (type != itk::ImageIOBase::FLOAT) {
        throw PyCAException(__FILE__,__LINE__,"We do not support type other than float at this time");
    }

    Field3D temp(grid, h.memType());
    size_t nVox = grid.nVox();

    boost::shared_ptr<float> tempH = CreateSharedArray<MEM_HOST_PINNED,float>(nVox);
    if (h.memType() != MEM_DEVICE) {
        CopyFromITKBuffer<float, float>(name.c_str(), temp.x, nVox);
        name = oPrefix + "_y.mha";
        CopyFromITKBuffer<float, float>(name.c_str(), temp.y, nVox);
        name = oPrefix + "_z.mha";
        CopyFromITKBuffer<float, float>(name.c_str(), temp.z, nVox);
    } else {
        CopyFromITKBuffer<float, float>(name.c_str(), tempH.get(), nVox);
        cpyArrayH2D(temp.x, tempH.get(), nVox);
        name = oPrefix + "_y.mha";
        CopyFromITKBuffer<float, float>(name.c_str(), tempH.get(), nVox);
        cpyArrayH2D(temp.y, tempH.get(), nVox);
        name = oPrefix + "_z.mha";
        CopyFromITKBuffer<float, float>(name.c_str(), tempH.get(), nVox);
        cpyArrayH2D(temp.z, tempH.get(), nVox);
    }
    h.swap(temp);
}

/** internal use only */
template<class ITKFilterType>
void SaveITKImage(typename ITKFilterType::Pointer& filter,
		  const char* name, 
		  bool useCompression)
{
    typedef typename ITKFilterType::OutputImageType ImageType;
    typedef typename itk::ImageFileWriter<ImageType> VolumeWriterType;
    //
    // create writer
    //
    typename VolumeWriterType::Pointer writer = VolumeWriterType::New();
    writer->SetFileName(name);
    writer->SetUseCompression(useCompression);
    writer->SetInput(filter->GetOutput());
    
    //
    // write image
    //
    try
    {
        writer->Update();
    }
    catch(itk::ImageFileWriterException &exc)
    {
        throw PyCAException(__FILE__, __LINE__, "ITK File Writer Exception", &exc);
    }
    catch(itk::ExceptionObject &exc)
    {
        throw PyCAException(__FILE__, __LINE__, "ITK Exception", &exc);
    }
    catch(...)
    {
        throw PyCAException(__FILE__, __LINE__, "Unknown exception");
    }
}

/** internal use only */
template<typename T>
void SaveImageBuffer(const char* name, 
		     const T* iBuf,
                     typename itk::ImportImageFilter<T, 3>::Pointer& importFilter,
                     size_t nVox, bool useCompression){

    typedef itk::Image<T, 3>                ImageType;
    const bool importImageFilterWillOwnTheBuffer = false;

    importFilter->SetImportPointer(const_cast<T*>(iBuf), nVox,
                                   importImageFilterWillOwnTheBuffer);
    
    SaveITKImage<itk::ImportImageFilter<T, 3> >(importFilter, name, useCompression);
}

void ITKFileIO::SaveImage(const Image3D& image, const char* name, bool useCompression){
    std::cerr << "Saving " << name << "...";
    size_t nVox = image.nVox();
    itk::ImportImageFilter<float, 3>::Pointer
        importFilter= CreateImportFilter<float>(image.grid());
        
    if (image.memType() != MEM_DEVICE) {
        SaveImageBuffer<float>(name, image.get(), importFilter, nVox, useCompression);
    } else {
        boost::shared_ptr<float> tempH = CreateSharedArray<MEM_HOST_PINNED,float>(nVox);
        cpyArrayD2H(tempH.get(), image.get(), nVox);
        SaveImageBuffer<float>(name, tempH.get(), importFilter, nVox, useCompression);
    }
    std::cerr << "done" << std::endl;
}

void ITKFileIO::SaveField(const Field3D& h, const char* name, bool useCompression)
{
    typedef itk::CovariantVector< float, 3> VectorPixelType;
    typedef itk::Image<VectorPixelType, 3> VectorImageType;
    typedef itk::Image<float, 3> ScalarImageType;

    size_t nVox = h.nVox();
    itk::ImportImageFilter<float, 3>::Pointer importFilterX = CreateImportFilter<float>(h.grid());
    itk::ImportImageFilter<float, 3>::Pointer importFilterY = CreateImportFilter<float>(h.grid());
    itk::ImportImageFilter<float, 3>::Pointer importFilterZ = CreateImportFilter<float>(h.grid());

#if ITK_VERSION_MAJOR < 4
    typedef itk::Compose3DCovariantVectorImageFilter
	<ScalarImageType,VectorImageType> 
	ComposeCovariantVectorImageFilterType;
#else
    typedef itk::ComposeImageFilter<ScalarImageType,VectorImageType> 
	ComposeCovariantVectorImageFilterType;
#endif
    
    const bool importImageFilterWillOwnTheBuffer = false;
    
    boost::shared_ptr<float> tempX, tempY, tempZ;
    float* buffX;
    float* buffY;
    float* buffZ;
    
    if (h.memType() != MEM_DEVICE) {
	buffX = h.x;
	buffY = h.y;
	buffZ = h.z;
    } else {
        tempX = CreateSharedArray<MEM_HOST_PINNED,float>(nVox);
        tempY = CreateSharedArray<MEM_HOST_PINNED,float>(nVox);
        tempZ = CreateSharedArray<MEM_HOST_PINNED,float>(nVox);

        cpyArrayD2H(tempX.get(), h.x, nVox);
        cpyArrayD2H(tempY.get(), h.y, nVox);
        cpyArrayD2H(tempZ.get(), h.z, nVox);

	buffX = tempX.get();
	buffY = tempY.get();
	buffZ = tempZ.get();
    }

    importFilterX->SetImportPointer(buffX, nVox,
				    importImageFilterWillOwnTheBuffer);
    importFilterY->SetImportPointer(buffY, nVox,
				    importImageFilterWillOwnTheBuffer);
    importFilterZ->SetImportPointer(buffZ, nVox,
				    importImageFilterWillOwnTheBuffer);

    ComposeCovariantVectorImageFilterType::Pointer composeFilter = 
	ComposeCovariantVectorImageFilterType::New();
    composeFilter->SetInput1(importFilterX->GetOutput());
    composeFilter->SetInput2(importFilterY->GetOutput());
    composeFilter->SetInput3(importFilterZ->GetOutput());
    composeFilter->Update();

    SaveITKImage< ComposeCovariantVectorImageFilterType >(composeFilter, name, useCompression);

}

void ITKFileIO::SaveFieldComponents(const Field3D& h, const char* prefix, bool useCompression)
{
    std::string oPrefix(prefix);
    std::string name;
    
    size_t nVox = h.nVox();
    itk::ImportImageFilter<float, 3>::Pointer importFilter= CreateImportFilter<float>(h.grid());
    
    if (h.memType() != MEM_DEVICE) {
        name = oPrefix + "_x.mha";
        SaveImageBuffer<float>(name.c_str(), h.x, importFilter, nVox, useCompression);
        name = oPrefix + "_y.mha";
        SaveImageBuffer<float>(name.c_str(), h.y, importFilter, nVox, useCompression);
        name = oPrefix + "_z.mha";
        SaveImageBuffer<float>(name.c_str(), h.z, importFilter, nVox, useCompression);
    } else {
        boost::shared_ptr<float> tempH = CreateSharedArray<MEM_HOST_PINNED,float>(nVox);

        name = oPrefix + "_x.mha";
        cpyArrayD2H(tempH.get(), h.x, nVox);
        SaveImageBuffer<float>(name.c_str(), tempH.get(), importFilter, nVox, useCompression);
        name = oPrefix + "_y.mha";
        cpyArrayD2H(tempH.get(), h.y, nVox);
        SaveImageBuffer<float>(name.c_str(), tempH.get(), importFilter, nVox, useCompression);
        name = oPrefix + "_z.mha";
        cpyArrayD2H(tempH.get(), h.z, nVox);
        SaveImageBuffer<float>(name.c_str(), tempH.get(), importFilter, nVox, useCompression);
    }
}

template <class T>
void 
ITKFileIO::
LoadAff3D(Aff3D<T> &aff, const char* filename)
{
  itk::TransformFileReader::Pointer transformreader = itk::TransformFileReader::New();
  // read ITK-style transform
  transformreader->SetFileName(filename);
  try{
    transformreader->Update();
    itk::AffineTransform<float,3>::Pointer itkTransform = 
      dynamic_cast<itk::AffineTransform<float,3>*>((*transformreader->GetTransformList()->begin()).GetPointer());
    if(!itkTransform)
      {
	throw PyCAException(__FILE__, __LINE__, 
			       "Could not read transform into an itk::AffineTransform<float, 3>");
      }

    // set the center to zero
//     itk::AffineTransform<float,3>::InputPointType center;
//     center[0] = center[1] = center[2] = 0.0f;
//     itkTransform->SetCenter(center);

    // convert to RealAffineTransform format
    
    itk::Matrix<float,3,3> matrix = itkTransform->GetMatrix();
    itk::FixedArray<float, 3> vector = itkTransform->GetOffset();
    for(int row = 0; row<3; row++){
      for(int col = 0; col<3; col++){
	aff.matrix(row, col) = matrix(row, col);
      }
    }
    
    aff.vector.x = vector[0];
    aff.vector.y = vector[1];
    aff.vector.z = vector[2];
    
  }catch (std::exception* e){
    throw PyCAException(__FILE__, __LINE__, "Error reading ITK transform", e);
  }
  
}

// template instantiation
template void ITKFileIO::LoadAff3D(Aff3D<float> &aff, const char* filename);
template void ITKFileIO::LoadAff3D(Aff3D<double> &aff, const char* filename);


} // end namespace PyCA

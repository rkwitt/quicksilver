import numpy as np
import vtk


def WrapNPAsVTKImageData(nparr,
                         spacing=[1.0,1.0,1.0],
                         origin=[1.0,1.0,1.0],
                         AutoCopy=False):

    vtkArr = vtk.vtkFloatArray()

    # the SetVoidArray function uses the python 'buffer object'
    # interface, which is supported by numpy.  Only single-segment
    # memory objects are supported, however, so if AutoCopy is true
    # and nparr is not single-segment, a copy will be made and the
    # result passed back (ascontiguousarray returns original array if
    # already contiguous).  If AutoCopy is false and nparr is not
    # single-segment, an exception will be raised.
    if AutoCopy:
        nparr = np.ascontiguousarray(nparr)
    vtkArr.SetVoidArray(nparr, np.prod(nparr.shape), 1)

    vtkIm = vtk.vtkImageData()
    vtkIm.GetPointData().SetScalars(vtkArr)
    #vtkIm.SetDimensions(sz.x, sz.y, sz.z)
    vtkIm.SetDimensions(nparr.shape[1], nparr.shape[0], nparr.shape[2])
    vtkIm.SetScalarType(vtk.VTK_FLOAT)
    vtkIm.SetSpacing(spacing[1], spacing[0], spacing[2])
    vtkIm.SetOrigin(origin[1], origin[0], origin[2])

    if AutoCopy:
        return (vtkIm, nparr)
    else:
        return vtkIm

    
def WrapPyCAAsVTKImageData(im, AutoCopy=False):

    # have to go through the numpy interface as it has 'buffer object'
    # python interface
    return WrapNPAsVTKImageData(im.asnp(),
                                spacing=im.spacing().tolist(),
                                origin=im.origin().tolist(),
                                AutoCopy=AutoCopy)

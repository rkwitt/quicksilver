"""File IO methods"""

import PyCA.Core as core
import PyCACommon as common

import numpy as np
import os.path


def SaveImageNPY(Im, fname):
    """Save an image in *.npy format (no origin/spacing info)"""
    np.save(fname, Im.asnp())


def LoadImageNPY(fname, mType=core.MEM_HOST):
    """Load an image in *.npy format (no origin/spacing info)"""
    npy = np.load(fname)
    I = common.ImFromNPArr(npy, mType=mType)
    del npy.f
    npy.close()


def SaveImageNPZ(Im, fname, useCompression=True):
    """Save an image in *.npz format (records origin/spacing info)"""
    origin = Im.grid().origin().tolist()
    spacing = Im.grid().spacing().tolist()
    if useCompression:
        np.savez_compressed(fname, data=Im.asnp(), origin=origin,
                spacing=spacing)
    else:
        np.savez(f, data=Im.asnp(), origin=origin, spacing=spacing)


def LoadImageNPZ(fname, mType=core.MEM_HOST):
    """Load an image in *.npz format (with proper origin/spacing info)"""
    # TODO: put some format validation in here
    npz = np.load(fname)
    spacing = core.Vec3Df()
    spacing.fromlist(npz['spacing'])
    origin = core.Vec3Df()
    origin.fromlist(npz['origin'])
    I = common.ImFromNPArr(npz['data'], mType=mType, sp=spacing, orig=origin)
    del npz.f
    npz.close()
    return I


def SaveFieldNPZ(v, fname, useCompression=True):
    """Save a field in *.npz format (records origin/spacing info)"""
    origin = v.grid().origin().tolist()
    spacing = v.grid().spacing().tolist()
    # we will save in X by Y by Z by 3 format
    arr = np.transpose(v.asnp(), (1,2,3,0))
    if useCompression:
        np.savez_compressed(fname, data=arr, origin=origin,
                spacing=spacing)
    else:
        np.savez(fname, data=arr, origin=origin, spacing=spacing)


def LoadFieldNPZ(fname, mType=core.MEM_HOST):
    """Load a field in *.npz format (with proper origin/spacing info)"""
    # TODO: put some format validation in here
    npz = np.load(fname)
    spacing = core.Vec3Df()
    spacing.fromlist(npz['spacing'])
    origin = core.Vec3Df()
    origin.fromlist(npz['origin'])
    f = common.FieldFromNPArr(npz['data'], mType=mType, sp=spacing, orig=origin)
    del npz.f
    npz.close()
    return f


def LoadImageDCM(fname, mType=core.MEM_HOST):
    """Load a DICOM image in *.dcm format (with proper origin/spacing info)"""
    try:
        import dicom
    except:
        raise Exception('PyCA DICOM support requires pydicom, not found.')

    ds = dicom.read_file(fname)
    pixels = ds.pixel_array
    origin = core.Vec3Df()
    spacing = core.Vec3Df()
    if 'ImagePositionPatient' in ds.__dict__:
        origin.fromlist([ds.ImagePositionPatient[0], ds.ImagePositionPatient[1], 0])
    else:
        origin.fromlist([0,0,0])
    if 'PixelSpacing' in ds.__dict__:
        spacing.fromlist([ds.PixelSpacing[0], ds.PixelSpacing[1], 1])
    else:
        spacing.fromlist([1,1,1])

    return common.ImFromNPArr(pixels, mType=mType, sp=spacing, orig=origin)


def SaveImage(Im, fname, useCompression=True):
    """Save an image in any known format"""
    # get file extension
    ext = os.path.splitext(fname)[1].lower()

    # dispatch based on file extension
    if ext == '.npy':
        SaveImageNPY(Im, fname)
    elif ext == '.npz':
        SaveImageNPZ(Im, fname, useCompression=useCompression)
    elif ext == '.png':
        common.SavePNGImage(Im, fname)
    else:
        try:
            common.SaveITKImage(Im, fname, useCompression=useCompression)
        except IOError:
            raise Exception('File extension "'+ext+'" unknown.')


def LoadImage(fname, mType=core.MEM_HOST):
    """Load an image"""
    # get file extension
    ext = os.path.splitext(fname)[1].lower()

    # dispatch based on file extension
    if ext == '.npy':
        return LoadImageNPY(fname, mType=mType)
    elif ext == '.npz':
        return LoadImageNPZ(fname, mType=mType)
    elif ext == '.png':
        return common.LoadPNGImage(fname, mType=mType)
    elif ext == '.dcm':
        return LoadImageDCM(fname, mType=mType)
    else:
        try:
            return common.LoadITKImage(fname, mType=mType)
        except IOError:
            raise Exception('Could not load file.  Either the file extension "'
                +ext+'" is unknown or the file "'+fname+'" is missing or '
                +'corrupted.')


def SaveField(v, fname, useCompression=True):
    """Save a field in any known format"""
    # get file extension
    ext = os.path.splitext(fname)[1].lower()

    # dispatch based on file extension
    if ext == '.npy':
        raise Exception('Numpy field format not implemented')
    elif ext == '.npz':
        SaveFieldNPZ(v, fname, useCompression=useCompression)
    else:
        try:
            common.SaveITKField(v, fname, useCompression=useCompression)
        except IOError:
            raise Exception('File extension "'+ext+'" unknown.')


def LoadField(fname, mType=core.MEM_HOST):
    """Load a field in any known format"""
    # get file extension
    ext = os.path.splitext(fname)[1].lower()

    # dispatch based on file extension
    if ext == '.npy':
        raise Exception('Numpy field format not implemented')
    elif ext == '.npz':
        return LoadFieldNPZ(fname, mType=mType)
    else:
        try:
            return common.LoadITKField(fname, mType=mType)
        except IOError:
            raise Exception('File extension "'+ext+'" unknown.')

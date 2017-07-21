import PyCA.Core as core
#
# import numpy and scipy (arrays and linear algebra goodies)
#

import numpy as np
import os, errno, re
import math

try:
    import scipy.ndimage as ndimage
    ndimage_loaded = True
except ImportError:
    print "Warning: scipy.ndimage not found, some functionality disabled"

try:
    import matplotlib.pyplot as plt
    plt_loaded = True
except:
    print "Warning: matplotlib import failed, some functionality disabled"


def TukeyWindow(imSz, alpha):
    tx = TukeyWindow1D(imSz[0], alpha)
    if len(imSz) < 2:
        tw = tx
        return tw
    ty = TukeyWindow1D(imSz[1], alpha)
    tw = np.dot(tx.T, ty)
    if len(imSz) > 2:
        tz = TukeyWindow1D(imSz[2], alpha)
        tz = np.expand_dims(tz, 2)
        tz = np.transpose(tz, (2, 0, 1))
        tw = np.expand_dims(tw, 2)
        tw = np.tile(tw, (1, 1, imSz[2])) * np.tile(tz, (imSz[0], imSz[1], 1))
    return tw


def TukeyWindow1D(N, alpha):
    tw = np.ones([1, N])
    b = 0.5-0.5*np.cos(np.pi*np.array(range(0, alpha))/(alpha-1.0))
    tw[0, 0:alpha] = b
    tw[0, -1:-1-alpha:-1] = b
    return tw


def AddWhiteNoise(im, sig):
    if sig > 0:
        im = im + np.random.normal(0.0, sig, im.shape)
    return im


def GaussianBlur(im, sig, mode='nearest'):
    im = ndimage.filters.gaussian_filter(im, sig, mode=mode)
    return im


def ImFromNPArr(arr, mType=core.MEM_HOST,
                sp=core.Vec3Df(1.0, 1.0, 1.0),
                orig=core.Vec3Df(0.0, 0.0, 0.0)):
    arr3d = np.atleast_3d(arr)
    imSz = core.Vec3Di()
    imSz.fromlist(arr3d.shape)
    grid = core.GridInfo(imSz, sp, orig)
    im = core.Image3D(grid, core.MEM_HOST)
    im.asnp()[:,:,:] = arr3d
    im.toType(mType)
    return im

def FieldFromNPArr(arr, mType=core.MEM_HOST,
                   sp=core.Vec3Df(1.0, 1.0, 1.0),
                   orig=core.Vec3Df(0.0, 0.0, 0.0)):
    """
    Expects X-by-Y-by-Z-by-3 or X-by-Y-by-2 array
    """
    arrx, arry, arrz = _ComponentArraysFromField(arr)
    fSz = core.Vec3Di()
    fSz.fromlist(arrx.shape)

    grid = core.GridInfo(fSz, sp, orig)
    f = core.Field3D(grid, core.MEM_HOST)
    (fArrX, fArrY, fArrZ) = f.asnp()
    fArrX[:,:,:] = arrx
    fArrY[:,:,:] = arry
    fArrZ[:,:,:] = arrz
    f.toType(mType)

    return f

def IsImage3D(ob):
    return issubclass(type(ob), core.Image3D)

def IsField3D(ob):
    return issubclass(type(ob), core.Field3D)

def AsNPCopy(indata):

    origType = indata.memType()
    indata.toType(core.MEM_HOST)

    if IsImage3D(indata):
        out = indata.asnp().copy()
    elif IsField3D(indata):
        core.Field3D(indata.grid(), core.MEM_HOST)
        out_x, out_y, out_z = indata.asnp()
        sz = indata.grid().size().tolist()
        sz.append(3)
        out = np.zeros(sz)
        out[:,:,:,0] = out_x
        out[:,:,:,1] = out_y
        out[:,:,:,2] = out_z
    else:
        raise Exception('Expected `out` to be Image3D or Field3D')

    indata.toType(origType)

    return out

def CopyDataFromNPArr(out, arrIn):
    if IsImage3D(out):
        _CopyImDataFromNPArr(out, arrIn)
    elif IsField3D(out):
        _CopyFieldDataFromNPArr(out, arrIn)
    else:
        raise Exception('Expected `out` to be Image3D or Field3D')

def _CopyImDataFromNPArr(imOut, arrIn):
    arr3d = np.atleast_3d(arrIn)
    imSz = core.Vec3Di()
    imSz.fromlist(arr3d.shape)
    if imSz != imOut.size():
        raise Exception('imOut and arrIn must have compatible sizes: ' +
                        str(imSz) + ' != ' + str(imOut.size()))
    origType = imOut.memType()
    imOut.toType(core.MEM_HOST)
    imOut.asnp()[:,:,:] = arr3d
    imOut.toType(origType)

def _CopyFieldDataFromNPArr(fOut, arrIn):
    """
    Expects X-by-Y-by-Z-by-3 or X-by-Y-by-2 array
    """
    arrx, arry, arrz = _ComponentArraysFromField(arrIn)
    fSz = core.Vec3Di()
    fSz.fromlist(arrx.shape[:3])

    if fSz != fOut.size():
        raise Exception('fOut and arrIn must have compatible sizes: ' +
                        str(fSz) + ' != ' + str(fOut.size()))

    origType = fOut.memType()
    fOut.toType(core.MEM_HOST)
    (fArrX, fArrY, fArrZ) = fOut.asnp()
    fArrX[:,:,:] = arrx
    fArrY[:,:,:] = arry
    fArrZ[:,:,:] = arrz
    fOut.toType(origType)

def _ComponentArraysFromField(arr):
    if len(arr.shape) == 4:
        if arr.shape[3] != 3:
            raise Exception('Expected X-by-Y-by-Z-by-3 array')
        arrx = arr[:,:,:,0]
        arry = arr[:,:,:,1]
        arrz = arr[:,:,:,2]
    elif len(arr.shape) == 3:
        if arr.shape[2] == 2:
            arrx = np.atleast_3d(arr[:,:,0])
            arry = np.atleast_3d(arr[:,:,1])
            arrz = np.zeros(arrx.shape)
        elif arr.shape[2] == 3:
            arrx = np.atleast_3d(arr[:,:,0])
            arry = np.atleast_3d(arr[:,:,1])
            arrz = np.atleast_3d(arr[:,:,2])
        else:
            raise Exception('Expected X-by-Y-by-Z-by-3 or X-by-Y-by-2 array')
    else:
        raise Exception('Expected X-by-Y-by-Z-by-3 or X-by-Y-by-2 array')
    return arrx, arry, arrz


def LoadITKImage(fname, mType=core.MEM_HOST):
    im = core.Image3D(mType)
    core._ITKFileIO.LoadImage(im, fname)
    return im

def LoadITKField(fname, mType=core.MEM_HOST):
    f = core.Field3D(mType)
    core._ITKFileIO.LoadField(f, fname)
    return f

def SaveITKImage(im, fname, useCompression=False):
    core._ITKFileIO.SaveImage(im, fname, useCompression)

def SaveITKField(f, fname, useCompression=False):
    core._ITKFileIO.SaveField(f, fname, useCompression)

def ReadGrid(fname):
    g = core.GridInfo()
    dataType = core._ITKFileIO.ReadHeader(g, fname)
    return g

def LoadPNGImage(fname, mType=core.MEM_HOST,
                 sp=core.Vec3Df(1.0, 1.0, 1.0),
                 orig=core.Vec3Df(0.0, 0.0, 0.0)):
    imArr = plt.imread(fname)
    im = ImFromNPArr(imArr, mType=mType, sp=sp, orig=orig)
    return im

def SavePNGImage(im, fname, rng=None, cmap='gray'):
    imArr = np.squeeze(AsNPCopy(im))
    if rng is None:
        vmax=None
        vmin=None
    else:
        vmin=rng[0]
        vmax=rng[1]
    plt.imsave(fname, imArr, cmap=cmap, vmin=vmin, vmax=vmax)

def ExtractROI(im, extent):

    roiStart = core.Vec3Di(extent[0], extent[2], extent[4])
    roiSize = core.Vec3Di(extent[1]-extent[0]+1,
                          extent[3]-extent[2]+1,
                          extent[5]-extent[4]+1)
    roiGrid = core.GridInfo(im.grid().size(),
                            im.grid().spacing(),
                            im.grid().origin())
    roiGrid.setSize(roiSize)
    roi = core.Image3D(roiGrid, im.memType())
    core.SubVol(roi, im, roiStart)
    return roi

DIMMAP = {0: 0, 'x': 0, 1: 1, 'y': 1, 2: 2, 'z': 2}


def ExtractSliceArr(arr, sliceIdx=None, dim='z'):
    """
    Given a numpy array 'arr', extract the given slice along the given
    dimension.  If sliceIdx is None, extract the middle slice.
    """
    arr = np.atleast_3d(arr)
    sz = arr.shape
    # take mid slice if none specified
    if sliceIdx is None:
        sliceIdx = sz[DIMMAP[dim]]//2
    roiStart = np.array([0, 0, 0])
    roiSize = np.array([sz[0], sz[1], sz[2]])
    roiStart[DIMMAP[dim]] = sliceIdx
    roiSize[DIMMAP[dim]] = 1
    roiEnd = roiStart+roiSize
    sliceIm = arr[roiStart[0]:roiEnd[0],
                  roiStart[1]:roiEnd[1],
                  roiStart[2]:roiEnd[2]]
    return np.squeeze(sliceIm)

def ExtractSliceIm(im, sliceIdx=None, dim='z'):
    """
    Given an Image3D 'im', extract the given slice along the given dimension.
    If sliceIdx is None, extract the middle slice.
    """
    sz = im.size().tolist()
    # take mid slice if none specified
    if sliceIdx is None:
        sliceIdx = sz[DIMMAP[dim]]//2
    roiStart = core.Vec3Di(0, 0, 0)
    roiSize = core.Vec3Di(sz[0], sz[1], sz[2])
    roiStart.set(DIMMAP[dim], sliceIdx)
    roiSize.set(DIMMAP[dim], 1)
    roiGrid = core.GridInfo(im.grid().size(),
                            im.grid().spacing(), im.grid().origin())
    roiGrid.setSize(roiSize)
    sliceIm = core.Image3D(roiGrid, im.memType())
    core.SubVol(sliceIm, im, roiStart)
    return sliceIm


def ExtractSliceArrIm(im, sliceIdx=None, dim='z'):
    """
    Given an Image3D 'im', extract the given slice along the given dimension
    and return a copy of it as a numpy array.  If sliceIdx is None, extract
    the middle slice.
    """
    sliceIm = ExtractSliceIm(im, sliceIdx, dim)
    sliceArr = AsNPCopy(sliceIm)
    return sliceArr


def ExtractSliceVF(vf, sliceIdx=None, dim='z'):
    """
    Given a Field3D 'vf', extract the given slice along the given dimension.
    If sliceIdx is None, extract the middle slice.
    """
    sz = vf.size().tolist()
    # take mid slice if none specified
    if sliceIdx is None:
        sliceIdx = sz[DIMMAP[dim]]//2
    roiStart = core.Vec3Di(0, 0, 0)
    roiSize = core.Vec3Di(sz[0], sz[1], sz[2])
    roiStart.set(DIMMAP[dim], sliceIdx)
    roiSize.set(DIMMAP[dim], 1)
    roiGrid = core.GridInfo(vf.grid().size(),
                            vf.grid().spacing(), vf.grid().origin())
    roiGrid.setSize(roiSize)
    sliceVF = core.Field3D(roiGrid, vf.memType())
    core.SubVol(sliceVF, vf, roiStart)
    return sliceVF


def ExtractSliceArrVF(vf, sliceIdx=None, dim='z'):
    """
    Given a Field3D 'vf', extract the given slice along the given dimension
    and return a copy of it as a numpy array.  If sliceIdx is None, extract
    the middle slice.
    """
    sliceVF = ExtractSliceVF(vf, sliceIdx, dim)
    sliceArr = AsNPCopy(sliceVF)
    return sliceArr


def ComposeDef(V, t, asVField=False, inverse=False,
               scratchV1=None, scratchV2=None):
    """
    Takes an array of Field3Ds and returns a Field3Ds
    containting the vectors Composed to non-integer time t
    """
    vlen = len(V)
    grid = V[0].grid()
    mType = V[0].memType()
    # just clamp to final time
    if t > vlen:
        t = vlen
    t_int = int(math.floor(t))
    t_frac = t-t_int
    h = core.Field3D(grid, mType)
    if scratchV1 is None:
        scratchV1 = core.Field3D(grid, mType)
    core.SetToIdentity(h)
    for s in range(t_int):
        if inverse:
            core.ComposeHVInv(scratchV1, h, V[s])
        else:
            core.ComposeVH(scratchV1, V[s], h)
        h.swap(scratchV1)
    if t_frac != 0.0:
        if scratchV2 is None:
            scratchV2 = core.Field3D(grid, mType)
        core.Copy(scratchV2, V[t_int])
        core.MulC_I(scratchV2, core.Vec3Df(t_frac, t_frac, t_frac))
        if inverse:
            core.ComposeHVInv(scratchV1, h, scratchV2)
        else:
            core.ComposeVH(scratchV1, scratchV2, h)
        core.Copy(h, scratchV1)
    if asVField:
        core.SetToIdentity(scratchV1)
        core.Sub_I(h, scratchV1)
    return h


def ReadVolWriteFrames(volName, frameFmt,
                       mType=core.MEM_HOST,
                       frameRange=None, t=False):
    vol = core.Image3D(mType)
    core._ITKFileIO.LoadImage(vol, volName)
    imSz = vol.size()
    sz = imSz.tolist()
    if frameRange is None:
        frameRange = [0, sz[2]]
    volArr = vol.asnp()
    for z in range(frameRange[0], frameRange[1]):
        sliceArr = volArr[:,:,z]
        if t:
            sliceArr = sliceArr.T
        plt.imsave(frameFmt % z, sliceArr, cmap='gray')


def SaveSliceStack(fname, imList, sliceIdx=None, dim='z'):
    nIm = len(imList)
    imSz = imList[0].size()
    sz = imSz.tolist()
    mType = imList[0].memType()
    volSz = core.Vec3Di(imSz)
    volSz.set(DIMMAP[dim], nIm)
    volGrid = core.GridInfo(volSz)
    vol = core.Image3D(volGrid, mType)
    start = core.Vec3Di(0,0,0)
    for imIdx in range(nIm):
        # if this volume is not a slice already, extract a slice
        im = imList[imIdx]
        if sz[DIMMAP[dim]] > 1:
            im = ExtractSliceArrIm(im, sliceIdx, dim)
        start.set(DIMMAP[dim], imIdx)
        core.SetSubVol_I(vol, im, start)
    core._ITKFileIO.SaveImage(vol, fname)


def Distribute(nEl, size, rank):
    """
    Distribute 'nEl' elements among 'size' containers, return the list of
    elements in container 'rank'
    nEl is the number of elements to distribute
    size is total number of process
    rank is the ID of this process
    returns list of indices for this process to handle
    """
    if size < 1:
        raise Exception('size cannot be zero')
    if size == 1:
        return range(nEl)
    else:
        if nEl < size:
            if rank < nEl:
                return [rank]
            else:
                return []
        else:
            # minimum number of elements given to any process
            nProcEl = int(math.floor(nEl/size))
            # number of remaining elments to split up
            remaining = nEl % size

            if rank < remaining:
                # one extra element in this bin
                nProcEl += 1
                startEl = nProcEl*rank
            else:
                startEl = (nProcEl+1)*remaining + nProcEl*(rank-remaining)
            return range(startEl, startEl+nProcEl)
    raise Exception('Should not get here')


def Mkdir_p(dir):
    try:
        print 'creating dir: ', dir
        os.makedirs(dir)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            print 'already exists'
        else:
            raise


def ReadPLUNC(fname):
    f = open(fname, 'r')
    # throw away first line
    f.readline()
    # create the matrix to return
    mat = np.zeros([4, 4])
    for m in range(4):
        matline = f.readline()
        # split on space or tab
        splitchars = re.compile('[ \t\n]+')
        matentries = splitchars.split(matline)
        matentries = filter(lambda s: len(s) > 0, matentries)
        if len(matentries) != 4:
            raise Exception('error, wrong number of entries')
        for n in range(4):
            mat[m, n] = float(matentries[n])
    return mat.T

#def ReadFramesWriteVol(frameFmt, volName):

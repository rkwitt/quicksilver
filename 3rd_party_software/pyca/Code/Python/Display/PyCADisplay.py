import PyCA.Core as core
import PyCA.Common as common

import numpy as np

# pyplot contains most of the plotting functionality
try:
    import matplotlib.pyplot as plt
    import matplotlib.colors
    plt_loaded = True
    # tell it to use interactive mode -- see results immediately
    plt.ion()
except:
    print "Warning: matplotlib import failed, some functionality disabled"


def EnergyPlot(energy, trim=0.05, legend=None):
    # trim
    ttlEn = np.array(energy[-1])
    nEl = len(ttlEn)
    ascendingEn = sorted(ttlEn)
    nToTrim = int(nEl*trim)
    trimVal = ascendingEn[-(1+nToTrim)]
    en = np.array(energy)
    en[en > trimVal] = float('NaN')
    # plot
    colorList = ['r', 'g', 'b', 'k', 'm', 'c', 'y']
    for idx in range(len(en)):
        plt.plot(en[idx][:], colorList[idx])
        if idx == 0:
            plt.hold(True)
    plt.hold(False)
    if legend is not None:
        plt.legend(legend)
    plt.draw()

    
def Quiver(vf, sliceIdx=None, dim='z'):
    sliceArr = common.ExtractSliceArrVF(vf, sliceIdx, dim)
    plt.quiver(np.squeeze(sliceArr[:,:,0]).T,
               np.squeeze(sliceArr[:,:,1]).T,
               scale=None)
    # change axes to match image axes
    if not plt.gca().yaxis_inverted():
        plt.gca().invert_yaxis()
        # force redraw
        plt.draw()

        
def GridPlot(vf, sliceIdx=None, dim='z', every=1,
             isVF=True, color='g', plotBase=True, colorbase='#A0A0FF'):
    sliceArr = np.squeeze(common.ExtractSliceArrVF(vf, sliceIdx, dim))
    sz = sliceArr.shape
    hID = np.mgrid[0:sz[0], 0:sz[1]]
    d1 = np.squeeze(hID[1, ::every, ::every])
    d2 = np.squeeze(hID[0, ::every, ::every])
    sliceArr = sliceArr[::every, ::every, :]
    if plotBase:
        plt.plot(d1, d2, colorbase)
        plt.hold(True)
        plt.plot(d1.T, d2.T, colorbase)
    if not isVF:
        d1 = np.zeros(d1.shape)
        d2 = np.zeros(d2.shape)
    plt.plot(d1+np.squeeze(sliceArr[:,:,1]),
             d2+np.squeeze(sliceArr[:,:,0]), color)
    plt.hold(True)
    plt.plot((d1+np.squeeze(sliceArr[:,:,1])).T,
             (d2+np.squeeze(sliceArr[:,:,0])).T, color)
    plt.hold(False)
    # change axes to match image axes
    if not plt.gca().yaxis_inverted():
        plt.gca().invert_yaxis()
        # force redraw
        plt.draw()

def JacDetCMap(jd_max=10.0, cmap='PRGn',
               nonpos_clr=(1.0, 0.0, 0.0, 1.0)):
    """
    Create a colormap for displaying jacobian determinant that clamps
    over and under values to colormap boundary colors, but sets 'bad'
    values to nonpos_clr since logmapped zero or negative numbers will
    be -inf or nan
    """
    jacdet_cmap = plt.cm.__dict__[cmap]
    jacdet_cmap.set_over(color=jacdet_cmap(1.0))
    jacdet_cmap.set_under(color=jacdet_cmap(0.0))
    jacdet_cmap.set_bad(color=nonpos_clr)
    return jacdet_cmap

def JacDetPlot(vf, title='Jac. Det.',
               jd_max=10.0, cmap='PRGn',
               nonpos_clr=(1.0, 0.0, 0.0, 1.0),
               sliceIdx=None, dim='z',
               isVF=True,
               newFig=False):
    """
    Plot the jacobian determinant using logmapped colors, and setting
    zero or negative values to 'nonpos_clr'.  jd_max is the maximum
    value to display without clamping, and also defines the min value
    as 1.0/jd_max to assure 1.0 is centered in the colormap.  If vf is
    a vector field, compute the jacobian determinant.  If it is an
    Image3D, assume it is the jacobian determinant to be plotted.
    """

    if common.IsField3D(vf):
        grid = vf.grid()
        mType = vf.memType()
        h = core.ManagedField3D(grid, mType)
        jacdet = core.ManagedImage3D(grid, mType)
        core.Copy(h, vf)
        if isVF:
            core.VtoH_I(h)
        core.JacDetH(jacdet, h)
    elif common.IsImage3D(vf):
        jacdet = vf
    else:
        raise Exception('unknown input type to JacDetPlot, %s'%\
                        str(type(vf)))
    jd_cmap = JacDetCMap(jd_max=jd_max, cmap=cmap,
                               nonpos_clr=nonpos_clr)
    DispImage(jacdet, title=title,
              sliceIdx=sliceIdx, dim=dim,
              rng=[1.0/jd_max, jd_max],
              cmap=jd_cmap,
              log=True,
              newFig=newFig)
        
def SaveSlice(fname, im, sliceIdx=None, dim='z',
              cmap='gray', rng=None, t=True):
    sliceIm = common.ExtractSliceIm(im, sliceIdx, dim)
    sliceIm.toType(core.MEM_HOST)
    if rng is None:
        vmin = None
        vmax = None
    else:
        vmin = rng[0]
        vmax = rng[1]
    if t:
        plt.imsave(fname, np.squeeze(sliceIm.asnp()).T,
                   cmap=cmap, vmin=vmin, vmax=vmax)
    else:
        plt.imsave(fname, np.squeeze(sliceIm.asnp()),
                   cmap=cmap, vmin=vmin, vmax=vmax)

        
def SaveThreePaneArr(npim, outname, cmap='gray', rng=None):

    sz = npim.shape
    xSlice = np.squeeze(npim[sz[0]/2, :, :])
    ySlice = np.squeeze(npim[:, sz[1]/2, :])
    zSlice = np.squeeze(npim[:, :, sz[2]/2])
    sumx = xSlice.shape[1] + ySlice.shape[1] + zSlice.shape[1]
    maxy = max([xSlice.shape[0], ySlice.shape[0], zSlice.shape[0]])
    threePane = np.zeros([maxy, sumx])
    starty = (maxy-xSlice.shape[0])/2
    threePane[starty:starty+xSlice.shape[0], 0:xSlice.shape[1]] = xSlice
    startx = xSlice.shape[1]
    starty = (maxy-ySlice.shape[0])/2
    threePane[starty:starty+ySlice.shape[0],
              startx:startx+ySlice.shape[1]] = ySlice
    startx = startx+ySlice.shape[1]
    starty = (maxy-zSlice.shape[0])/2
    threePane[starty:starty+zSlice.shape[0],
              startx:startx+zSlice.shape[1]] = zSlice
    
    if rng is None:
        vmin = None
        vmax = None
    else:
        vmin = rng[0]
        vmax = rng[1]

    plt.imsave(outname,
               threePane, cmap=cmap,
               vmin=vmin, vmax=vmax)


#
# Takes an array of Field3Ds, an Image3D, and a series of
# timepoints, and a function and optional args.  Calls function as:
# func(Idef, t, *args)
# for deformed image at each timepoint, and appends the result to a
# list that is returned from this function
#
def DefSeriesIter(I, V, t, func, args):
    grid = I.grid()
    mType = I.memType()
    tlen = len(t)
    h = core.Field3D(grid, mType)
    IDef = core.Image3D(grid, mType)
    core.SetToIdentity(h)
    scratchV1 = core.Field3D(grid, mType)
    scratchV2 = core.Field3D(grid, mType)
    rtnarr = []
    for tidx in range(tlen):
        curt = t[tidx]
        h = common.ComposeDef(V, curt, inverse=True,
                              asVField=False,
                              scratchV1=scratchV1,
                              scratchV2=scratchV2)
        core.ApplyH(IDef, I, h)
        r = func(IDef, curt, *args)
        rtnarr.append(r)
    return rtnarr


#
# Takes an array of Field3Ds, an Image3D, and a series of
# timepoints and returns a volume containing the frames at these times
#
def ExtractSliceFunc(IDef, t):
    sliceIm = common.ExtractSliceIm(IDef)
    return sliceIm


def GenDefFrames(I, V, t):
    """
    Take an image, a series of velocity fields, and a set of
    timepoints, and return an image that is a stack of slices of I
    deformed to each of the timepoints t.
    """
    grid = I.grid()
    mType = I.memType()
    
    sliceArr = \
        DefSeriesIter(I, V, t, func=ExtractSliceFunc, args=())

    nSlice = len(sliceArr)

    outGrid = grid.copy()
    outGrid.size().z = nSlice
    sliceStack = core.Image3D(outGrid, mType)
    start = core.Vec3Di(0,0,0)

    for sIdx in range(nSlice):
        start.z = sIdx
        core.SetSubVol_I(sliceStack, sliceArr[sIdx], start)
    return sliceStack


# display an Image3D slice
def DispImage(im, title=None,
              sliceIdx=None, dim='z',
              cmap='gray', newFig=True,
              rng=None, t=False, log=False):

    # if this volume is not a slice already, extract a slice
    sz = im.size().tolist()
    if sz[common.DIMMAP[dim]] > 1:
        im = common.ExtractSliceIm(im, sliceIdx, dim)
        im.toType(core.MEM_HOST)

    # transfer to host memory if necessary
    if im.memType() == core.MEM_DEVICE:
        tmp = core.Image3D(im.grid(), core.MEM_HOST)
        core.Copy(tmp, im)
        im = tmp
        
    # create the figure if requested
    if newFig:
        if title is None:
            plt.figure()
        else:
            plt.figure(title)
        plt.clf()
            
    # set the range
    if rng is None:
        vmin = None
        vmax = None
    else:
        vmin = rng[0]
        vmax = rng[1]

    if log:
        norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

    # convert to numpy array
    arr = np.squeeze(im.asnp().copy())

    # transpose if requested
    if t:
        arr = arr.T

    # display
    plt.imshow(arr, cmap=cmap, vmin=vmin, vmax=vmax, norm=norm, interpolation='nearest')
    plt.axis('tight')
    plt.axis('image')
    if title is not None:
        plt.title(title)
    plt.xticks([])
    plt.yticks([])
    plt.draw()

def SampleCMap(nSamples, cmap=matplotlib.cm.jet):
    """
    Evenly sample a colormap.  Used to create easy index-to-color
    mappings.
    """
    norm = matplotlib.colors.Normalize(vmin=0, vmax=nSamples)
    sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    colors = sm.to_rgba(np.arange(nSamples))
    return colors

    

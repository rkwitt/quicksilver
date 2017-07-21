import PyCA.Core as core
import math
import numpy as np

import PyCACommon as common
reload(common)


def DrawCircle(sz, center, rad):
    grid = np.mgrid[0:sz[0], 0:sz[1]]
    coords = grid - np.transpose(np.atleast_3d(center), (1,0,2))
    d = np.squeeze(np.sqrt(np.sum(coords**2, 0)))
    circ = np.zeros(sz)
    circ[d < rad] = 1.0
    return circ


def DrawEllipse(sz, center, radx, rady):
    grid = np.mgrid[0:sz[0], 0:sz[1]]
    coords = grid - np.transpose(np.atleast_3d(center), (1,0,2))
    coords = coords / np.transpose(
        np.atleast_3d(np.array([float(radx), float(rady)])), (1, 0, 2))
    d = np.squeeze(np.sum(coords**2, 0))
    circ = np.zeros(sz)
    circ[d < 1.0] = 1.0
    return circ


def DrawRect(sz, upleft, wh):
    rect = np.zeros(sz)
    rect[upleft[0]:upleft[0]+wh[0], upleft[1]:upleft[1]+wh[1]] = 1.0
    return rect


def DrawRectFromCenter(sz, center, wh):
    rect = np.zeros(sz)
    r = np.array(wh)/2
    rect[center[0]-r[0]:center[0]+(wh[0]-r[0]),
         center[1]-r[1]:center[1]+(wh[1]-r[1])] = 1.0
    return rect


def DrawSin(sz, nperiods):
    grid = np.meshgrid(np.linspace(0, 2*np.pi*nperiods, sz[1]),
                       np.linspace(0, 2*np.pi*nperiods, sz[0]))
    bg = 0.5*((np.sin(grid[1])*np.sin(grid[0]))+1.0)
    return bg


def DrawLine(sz, slope, intercept, width):
    # normalize
    A = 1.0/math.sqrt(1+slope**2)
    B = -slope*A
    C = -intercept
    #print A, B, C
    grid = np.mgrid[0:sz[0], 0:sz[1]]
    grid[0,:,:] = grid[0,:,:]*A
    grid[1,:,:] = grid[1,:,:]*B
    d = np.abs(np.squeeze(np.sum(grid, 0)+C))
    line = np.zeros(sz)
    line[d < width/2] = 1.0
    return line


def RandImage(sz, nSig=1.0, gSig=0.0, NonNeg=False,
              mType=core.MEM_HOST,
              sp=core.Vec3Df(1.0,1.0,1.0),
              orig=core.Vec3Df(0.0,0.0,0.0)):
    imArr = nSig*np.random.randn(*sz)
    if gSig > 0:
        imArr = common.GaussianBlur(imArr, gSig)
    if NonNeg:
        imArr = np.abs(imArr)
    im = common.ImFromNPArr(imArr, mType, sp=sp, orig=orig)
    return im


def RandUnifImage(sz, lbound=0, ubound=1,
                  mType=core.MEM_HOST,
                  sp=core.Vec3Df(1.0,1.0,1.0),
                  orig=core.Vec3Df(0.0,0.0,0.0)):
    imArr = np.random.rand(*sz)*(ubound - lbound) + lbound
    im = common.ImFromNPArr(imArr, mType, sp=sp, orig=orig)
    return im


def RandField(sz, nSig=1.0, gSig=0.0, NonNeg=False,
              mType=core.MEM_HOST,
              sp=core.Vec3Df(1.0,1.0,1.0),
              orig=core.Vec3Df(0.0,0.0,0.0)):
    vSz = list(sz)
    dim = len(vSz)
    assert dim == 2 or dim == 3
    if dim == 3 and vSz[2] == 1:
        dim = 2
    if dim == 2:
        vSz.append(1)
    vSz.append(3)
    vArr = np.zeros(vSz)
    for d in range(dim):
        vDimArr = np.atleast_3d(nSig*np.random.randn(*vSz[0:2]))
        if gSig > 0:
            vDimArr = common.GaussianBlur(vDimArr, gSig)
        vArr[:,:,:,d] = vDimArr
    if NonNeg:
        vArr = np.abs(vArr)
    v = common.FieldFromNPArr(vArr, mType, sp=sp, orig=orig)
    return v

def DrawChecker(sz, sqrSz):
    gridy, gridx = np.meshgrid(np.arange(0, sz[1]),
                               np.arange(0, sz[0]))
    patx = (gridx//sqrSz)%2
    paty = (gridy//sqrSz)%2
    bg = np.zeros(sz)
    bg[patx+paty == 1] = 1.0
    return bg

def DrawSlantyDiamonds(sz, zagSz):
    gridy, gridx = np.meshgrid(np.arange(0, sz[1]),
                               np.arange(0, sz[0]))
    gridx = gridx + gridy%zagSz
    bg = (gridx//zagSz)%2
    return bg

def DrawDiamonds(sz, zagSz):
    gridy, gridx = np.meshgrid(np.arange(0, sz[1]),
                               np.arange(0, sz[0]))
    gridA = ((gridx + gridy)//zagSz)%2
    gridB = ((gridx - gridy + sz[1])//zagSz)%2
    bg = np.zeros(sz)
    bg[gridA+gridB != 1] = 1.0
    return bg

def DrawSlantStripes(sz, zagSz):
    gridy, gridx = np.meshgrid(np.arange(0, sz[1]),
                               np.arange(0, sz[0]))
    bg = ((gridx + gridy)//zagSz)%2
    return bg

def DrawWavyStripes(sz, nWaves=5, waveSz=5.0, stripeWidth=12):
    gridy, gridx = np.meshgrid(np.arange(0, sz[1]),
                               np.arange(0, sz[0]))
    off = waveSz*np.sin(np.pi*nWaves*gridy/float(sz[1]))
    bg = ((gridx + off)//stripeWidth)%2
    return bg

def DrawWavySinStripes(sz, nWaves=5, waveSz=5.0, nStripes=12):
    gridy, gridx = np.meshgrid(np.arange(0, sz[1]),
                               np.arange(0, sz[0]))
    off = (waveSz/np.pi)*np.sin(np.pi*nWaves*gridy/float(sz[1]))
    bg = np.sin(2*np.pi*nStripes*(gridx + off)/float(sz[0]))
    # scale to [0,1]
    bg = 0.5*(bg + 1.0)
    return bg

def DrawSinSlantStripes(sz, zagSz):
    gridy, gridx = np.meshgrid(np.arange(0, sz[1]),
                               np.arange(0, sz[0]))
    bg = np.sin((gridx + gridy)*(np.pi/float(zagSz)))
    bg = bg/2.0 + 0.5
    return bg

def WavyDef(sz, nWaves=4, waveAmp=10, waveDim=0,
            mType=core.MEM_HOST, deformation=True):
    """
    Generate and return a 'wavy' vector field.  If deformation is
    True, return a deformation, otherwise return a displacement.
    """
    gridy, gridx = np.meshgrid(np.arange(0, sz[1]),
                               np.arange(0, sz[0]))
    if waveDim == 0:
        wGrid = gridy
    elif waveDim == 1:
        wGrid = gridx
    else:
        raise Exception('waveDim should be 0 or 1')
    
    vfArr = np.zeros((sz[0], sz[1], 2))
    vfArr[:,:,waveDim] = \
        waveAmp*np.sin(2*np.pi*nWaves*wGrid/float(sz[1]))
    
    vf = common.FieldFromNPArr(vfArr, mType=mType)

    if deformation:
        core.VtoH_I(vf)
        
    return vf

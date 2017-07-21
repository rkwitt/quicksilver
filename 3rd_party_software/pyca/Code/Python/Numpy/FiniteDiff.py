"""
Finite differences implemented in numpy.  Meant as reference
implementation to test PyCA functions against, not as optimized
implementations.
"""

import numpy as np

from PyCA.Core import Vec3Df
# need constants from PyCA
from PyCA.Core import DIM_X, DIM_Y, DIM_Z
from PyCA.Core import DIFF_FORWARD, DIFF_BACKWARD, DIFF_CENTRAL
from PyCA.Core import BC_APPROX, BC_CLAMP, BC_WRAP

#
# Numpy finite difference implementation to test against
#
def FiniteDiff(arr, dim, diffType, bc, sp=Vec3Df(1.0, 1.0, 1.0)):

    diff = np.zeros(arr.shape)
    if dim == DIM_X:
        if diffType == DIFF_FORWARD:
            diff[:-1,:] = (arr[1:,:]-arr[:-1,:])/sp.x
            if bc == BC_APPROX:
                pass
            elif bc == BC_CLAMP:
                pass
            elif bc == BC_WRAP:
                diff[-1,:] = (arr[0,:]-arr[-1,:])/sp.x
            else:
                raise Exception('Unknwn boundary condition')
        elif diffType == DIFF_BACKWARD:
            diff[1:,:] = (arr[1:,:]-arr[:-1,:])/sp.x
            if bc == BC_APPROX:
                pass
            elif bc == BC_CLAMP:
                pass
            elif bc == BC_WRAP:
                diff[0,:] = (arr[0,:]-arr[-1,:])/sp.x
            else:
                raise Exception('Unknwn boundary condition')
        elif diffType == DIFF_CENTRAL:
            diff[1:-1,:] = (arr[2:,:]-arr[:-2,:])/(2.0*sp.x)
            if bc == BC_APPROX:
                diff[0,:] = (arr[1,:]-arr[0,:])/sp.x
                diff[-1,:] = (arr[-1,:]-arr[-2,:])/sp.x
            elif bc == BC_CLAMP:
                diff[0,:] = (arr[1,:]-arr[0,:])/(2.0*sp.x)
                diff[-1,:] = (arr[-1,:]-arr[-2,:])/(2.0*sp.x)
            elif bc == BC_WRAP:
                diff[0,:] = (arr[1,:]-arr[-1,:])/(2.0*sp.x)
                diff[-1,:] = (arr[0,:]-arr[-2,:])/(2.0*sp.x)
        else:
            raise Exception('Unknwn diff type')
    elif dim == DIM_Y:
        if diffType == DIFF_FORWARD:
            diff[:,:-1] = (arr[:,1:]-arr[:,:-1])/sp.y
            if bc == BC_APPROX:
                pass
            elif bc == BC_CLAMP:
                pass
            elif bc == BC_WRAP:
                diff[:,-1] = (arr[:,0]-arr[:,-1])/sp.y
            else:
                raise Exception('Unknwn boundary condition')
        elif diffType == DIFF_BACKWARD:
            diff[:,1:] = (arr[:,1:]-arr[:,:-1])/sp.y
            if bc == BC_APPROX:
                pass
            elif bc == BC_CLAMP:
                pass
            elif bc == BC_WRAP:
                diff[:,0] = (arr[:,0]-arr[:,-1])/sp.y
            else:
                raise Exception('Unknwn boundary condition')
        elif diffType == DIFF_CENTRAL:
            diff[:,1:-1] = (arr[:,2:]-arr[:,:-2])/(2.0*sp.y)
            if bc == BC_APPROX:
                diff[:,0] = (arr[:,1]-arr[:,0])/sp.y
                diff[:,-1] = (arr[:,-1]-arr[:,-2])/sp.y
            elif bc == BC_CLAMP:
                diff[:,0] = (arr[:,1]-arr[:,0])/(2.0*sp.y)
                diff[:,-1] = (arr[:,-1]-arr[:,-2])/(2.0*sp.y)
            elif bc == BC_WRAP:
                diff[:,0] = (arr[:,1]-arr[:,-1])/(2.0*sp.y)
                diff[:,-1] = (arr[:,0]-arr[:,-2])/(2.0*sp.y)
        else:
            raise Exception('Unknwn diff type')
    elif dim == DIM_Z:
        raise Exception('DIM_Z unimplemented')
    else:
        raise Exception('Unknwn dimension')
    return diff
# end FiniteDiff

#
# Numpy finite difference implementation to test against
#
def Grad(arr, diffType, bc, sp=Vec3Df(1.0, 1.0, 1.0)):
    dx = FiniteDiff(arr, DIM_X, diffType, bc, sp)
    dy = FiniteDiff(arr, DIM_Y, diffType, bc, sp)
    g = np.concatenate((np.atleast_3d(dx), np.atleast_3d(dy)),axis=2)
    return g

#
# Numpy implementation of 2D jacobian determinant
#
def JacDet(varr, diffType, bc, sp=Vec3Df(1.0, 1.0, 1.0)):
    gx = Grad(varr[:,:,0], diffType, bc, sp)
    gy = Grad(varr[:,:,1], diffType, bc, sp)
    jacDet = gx[:,:,0]*gy[:,:,1] - gx[:,:,1]*gy[:,:,0]
    return jacDet

#
# Numpy masked finite difference implementation to test against
#
def FiniteDiffMask(arr, mask, dim, diffType, bc, sp=Vec3Df(1.0, 1.0, 1.0)):

    diff = np.zeros(arr.shape)
    if dim == DIM_X:
        if diffType == DIFF_FORWARD:
            v = arr[:-1,:]
            vn = arr[1:,:].copy()
            vn[mask[1:,:] == 0] = v[mask[1:,:]==0]
            vp = arr[:-1,:]
            diff[:-1,:] = (vn-vp)/sp.x
            if bc == BC_APPROX:
                pass
            elif bc == BC_CLAMP:
                pass
            elif bc == BC_WRAP:
                raise Exception('wrap bc unimplemented for masked fd')
            else:
                raise Exception('Unknwn boundary condition')
        elif diffType == DIFF_BACKWARD:
            v = arr[1:,:]
            vn = arr[1:,:]
            vp = arr[:-1,:].copy()
            vp[mask[:-1,:] == 0] = v[mask[:-1,:] == 0]
            diff[1:,:] = (vn-vp)/sp.x
            if bc == BC_APPROX:
                pass
            elif bc == BC_CLAMP:
                pass
            elif bc == BC_WRAP:
                raise Exception('wrap bc unimplemented for masked fd')
            else:
                raise Exception('Unknwn boundary condition')
        elif diffType == DIFF_CENTRAL:
            v = arr[1:-1,:]
            vn = arr[2:,:].copy()
            vn[mask[2:,:] == 0] = v[mask[2:,:] == 0]
            vp = arr[:-2,:].copy()
            vp[mask[:-2,:] == 0] = v[mask[:-2,:] == 0]
            if bc == BC_APPROX:
                diff[1:-1,:] = (vn-vp)/(sp.x)
            elif bc == BC_CLAMP:
                diff[1:-1,:] = (vn-vp)/(2.0*sp.x)
            elif bc == BC_WRAP:
                raise Exception('wrap bc unimplemented for masked fd')
        else:
            raise Exception('Unknwn diff type')
    elif dim == DIM_Y:
        if diffType == DIFF_FORWARD:
            v = arr[:,:-1]
            vn = arr[:,1:].copy()
            vn[mask[:,1:] == 0] = v[mask[:,1:] == 0]
            vp = arr[:,:-1]
            diff[:,:-1] = (vn-vp)/sp.y
            if bc == BC_APPROX:
                pass
            elif bc == BC_CLAMP:
                pass
            elif bc == BC_WRAP:
                raise Exception('wrap bc unimplemented for masked fd')
            else:
                raise Exception('Unknwn boundary condition')
        elif diffType == DIFF_BACKWARD:
            v = arr[:,1:]
            vn = arr[:,1:]
            vp = arr[:,:-1].copy()
            vp[mask[:,:-1] == 0] = v[mask[:,:-1] == 0]
            diff[:,1:] = (vn-vp)/sp.y
            if bc == BC_APPROX:
                pass
            elif bc == BC_CLAMP:
                pass
            elif bc == BC_WRAP:
                raise Exception('wrap bc unimplemented for masked fd')
            else:
                raise Exception('Unknwn boundary condition')
        elif diffType == DIFF_CENTRAL:
            v = arr[:,1:-1]
            vn = arr[:,2:].copy()
            vn[mask[:,2:] == 0] = v[mask[:,2:] == 0]
            vp = arr[:,:-2].copy()
            vp[mask[:,:-2] == 0] = v[mask[:,:-2] == 0]
            if bc == BC_APPROX:
                diff[:,1:-1] = (vn-vp)/(sp.y)
            elif bc == BC_CLAMP:
                diff[:,1:-1] = (vn-vp)/(2.0*sp.y)
            elif bc == BC_WRAP:
                raise Exception('wrap bc unimplemented for masked fd')
        else:
            raise Exception('Unknwn diff type')
    elif dim == DIM_Z:
        raise Exception('DIM_Z unimplemented')
    else:
        raise Exception('Unknwn dimension')

    diff[mask == 0] = 0

    return diff

#
# Numpy masked gradient finite difference implementation to test
# against
#
def GradMask(arr, mask, diffType, bc, sp=Vec3Df(1.0, 1.0, 1.0)):
    dx = FiniteDiffMask(arr, mask, DIM_X, diffType, bc, sp)
    dy = FiniteDiffMask(arr, mask, DIM_Y, diffType, bc, sp)
    g = np.concatenate((np.atleast_3d(dx), np.atleast_3d(dy)),axis=2)
    return g


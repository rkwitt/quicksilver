#
# This file contains testing where PyCA results are compared to
# results from numpy.  All tests can be run from the command line by:
#
# > python -m unittest discover -v -p '*UnitTest.py'
#
# To run an individual test with graphical output from ipython:
#
# import CpuGpuUnitTest as cgtest
# cgtc = cgtest.CpuGpuTestCase()
# cgtc.test_Exp(disp=True)
#
import sys
import unittest
import PyCATest

from PyCA.Core import *

import PyCA.Common as common
reload(common)

import numpy as np

try:
    import matplotlib.pyplot as plt
    plt.ion()
except ImportError:
    print "Warning: matplotlib.pyplot not found, some functionality disabled"

    
def CheckIm(hIm, dIm, name, disp):
    hdIm = dIm.copy()
    hdIm.toType(MEM_HOST)
    dImArr = np.squeeze(hdIm.asnp())
    hImArr = np.squeeze(hIm.asnp())
    diff = hImArr-dImArr
    if disp:
        title = name
        plt.figure(title)
        plt.subplot(1,3,1)
        plt.imshow(hImArr)
        plt.colorbar();
        plt.title('host')
        plt.draw()
        plt.subplot(1,3,2)
        plt.imshow(dImArr)
        plt.colorbar();
        plt.title('device')
        plt.draw()
        plt.subplot(1,3,3)
        plt.imshow(diff)
        plt.colorbar();
        plt.title('diff')
        plt.draw()
        plt.show()
    diffMax = np.max(np.abs(diff))
    diffAv = np.sum(np.abs(diff))/np.prod(diff.shape)
    return (diffAv, diffMax)

def CheckField(hF,dF,name,disp):
    hdF = dF.copy()
    hdF.toType(MEM_HOST)
    dFArr_x, dFArr_y, dFArr_z = hdF.asnp()
    dFArr_x = np.squeeze(dFArr_x)
    dFArr_y = np.squeeze(dFArr_y)
    hFArr_x, hFArr_y, hFArr_z = hF.asnp()
    hFArr_x = np.squeeze(hFArr_x)
    hFArr_y = np.squeeze(hFArr_y)
    diff_x = hFArr_x-dFArr_x
    diff_y = hFArr_y-dFArr_y
    if disp:
        title = name
        plt.figure(title)

        plt.subplot(2,3,1)
        plt.imshow(np.squeeze(hFArr_x))
        plt.colorbar();
        plt.title('host x')
        plt.draw()
        plt.subplot(2,3,2)
        plt.imshow(np.squeeze(dFArr_x))
        plt.colorbar();
        plt.title('device x')
        plt.draw()
        plt.subplot(2,3,3)
        plt.imshow(np.squeeze(diff_x))
        plt.colorbar();
        plt.title('diff x')
        plt.draw()

        plt.subplot(2,3,4)
        plt.imshow(np.squeeze(hFArr_y))
        plt.colorbar();
        plt.title('host y')
        plt.draw()
        plt.subplot(2,3,5)
        plt.imshow(np.squeeze(dFArr_y))
        plt.colorbar();
        plt.title('device y')
        plt.draw()
        plt.subplot(2,3,6)
        plt.imshow(np.squeeze(diff_y))
        plt.colorbar();
        plt.title('diff y')
        plt.draw()
        plt.show()
    diffMax = max(np.max(np.abs(diff_x)), np.max(np.abs(diff_y)))
    diffAv = (np.sum(np.abs(diff_x)) + np.sum(np.abs(diff_y))) \
        / (2*np.prod(diff_x.shape))
    return (diffAv, diffMax)

#
# Test Class
#
class CpuGpuTestCase(unittest.TestCase):

    def __init__(self, methodName='runTest'):
        super(CpuGpuTestCase, self).__init__(methodName)
        
        self.cudaEnabled = (GetNumberOfCUDADevices() > 0)

        if self.cudaEnabled:
            # allowable average abs. diff
            self.AvEps = 1e-6
            # allowable max abs. diff
            self.MaxEps = 1e-4
            # image size
            self.sz = np.array([127, 119])
            # spacing
            self.sp = np.array([1.5, 2.1])
            # fluid parameters
            self.fluidParams = [1.0, 1.0, 0.0]

            self.vsz = np.append(self.sz, 2)
            self.imSz = Vec3Di(int(self.sz[0]), int(self.sz[1]), 1)
            self.imSp = Vec3Df(float(self.sp[0]), float(self.sp[1]), 1.0)
            # set up grid
            self.grid = GridInfo(self.imSz, self.imSp)

            # set up host / device images
            self.I0Arr = common.DrawEllipse(self.sz, self.sz/2,
                                            self.sz[0]/4, self.sz[1]/3)
            self.I1Arr = common.DrawEllipse(self.sz, self.sz/2,
                                            self.sz[0]/3, self.sz[1]/4)

            self.I0Arr = common.GaussianBlur(self.I0Arr,1.5)
            self.I1Arr = common.GaussianBlur(self.I1Arr,1.5)

            self.hI0Orig = common.ImFromNPArr(self.I0Arr, 
                                              mType=MEM_HOST, 
                                              sp=self.imSp)
            self.hI1Orig = common.ImFromNPArr(self.I1Arr, 
                                              mType=MEM_HOST, 
                                              sp=self.imSp)
            self.dI0Orig = common.ImFromNPArr(self.I0Arr, 
                                              mType=MEM_DEVICE, 
                                              sp=self.imSp)
            self.dI1Orig = common.ImFromNPArr(self.I1Arr, 
                                              mType=MEM_DEVICE, 
                                              sp=self.imSp)

    # automatically called before test functions
    def setUp(self):
        if not self.cudaEnabled:
            self.skipTest('Cannot run test, no CUDA device found or CUDA support not compiled')

    def TestIm(self, pycaI, npI, name, disp, avEps=None, maxEps=None):

        if avEps is None:
            avEps = self.AvEps
        if maxEps is None:
            maxEps = self.MaxEps

        diffAv, diffMax = CheckIm(pycaI, npI, name=name, disp=disp)

        self.assertLess(diffAv, avEps)
        self.assertLess(diffAv, self.MaxEps)

    def TestField(self, pycaF, npF, name, disp, avEps=None, maxEps=None):

        if avEps is None:
            avEps = self.AvEps
        if maxEps is None:
            maxEps = self.MaxEps

        diffAv, diffMax = CheckField(pycaF, npF, name=name, disp=disp)

        self.assertLess(diffAv, avEps)
        self.assertLess(diffAv, self.MaxEps)
        

    def setUpI0(self):
        self.hI0 = self.hI0Orig.copy()
        self.dI0 = self.dI0Orig.copy()

    def tearDownI0(self):
        self.hI0 = None
        self.dI0 = None

    def setUpI1(self):
        self.hI1 = self.hI0Orig.copy()
        self.dI1 = self.dI0Orig.copy()

    def tearDownI1(self):
        self.hI1 = None
        self.dI1 = None

    def setUpDiffOp(self):
        self.hDiffOp = FluidKernelFFTCPU()
        self.hDiffOp.setAlpha(self.fluidParams[0])
        self.hDiffOp.setBeta(self.fluidParams[1])
        self.hDiffOp.setGamma(self.fluidParams[2])
        self.hDiffOp.setGrid(self.grid)

        self.dDiffOp = FluidKernelFFTGPU()
        self.dDiffOp.setAlpha(self.fluidParams[0])
        self.dDiffOp.setBeta(self.fluidParams[1])
        self.dDiffOp.setGamma(self.fluidParams[2])
        self.dDiffOp.setGrid(self.grid)

    def tearDownDiffOp(self):
        self.hDiffOp = None
        self.dDiffOp = None

    def setUpGrad(self):
        self.hGrad = Field3D(self.grid, MEM_HOST)
        self.dGrad = Field3D(self.grid, MEM_DEVICE)
        Gradient(self.hGrad, self.hI0Orig)
        Gradient(self.dGrad, self.dI0Orig)

    def tearDownGrad(self):
        self.hGrad = None
        self.dGrad = None

    def randImSetUp(self):
        self.hRandIm = \
            common.RandImage(self.sz, nSig=1.0, gSig=0.0, 
                             mType = MEM_HOST, sp = self.imSp)
        self.dRandIm = self.hRandIm.copy()
        self.dRandIm.toType(MEM_DEVICE)

    def randImTearDown(self):
        self.hRandIm = None
        self.dRandIm = None

    def randVPair(self):
        hV = common.RandField(self.sz, nSig=5.0, gSig=4.0, 
                              mType = MEM_HOST, sp = self.imSp)
        dV = hV.copy()
        dV.toType(MEM_DEVICE)
        return hV, dV

    def randVSetUp(self):
        hV, dV = self.randVPair()
        self.hRandV = hV
        self.dRandV = dV

    def randVTearDown(self):
        self.hRandV = None
        self.dRandV = None

    def randHSetUp(self):
        self.hRandH = \
            common.RandField(self.sz, nSig=5.0, gSig=4.0, 
                             mType = MEM_HOST, sp = self.imSp)
        VtoH_I(self.hRandH)
        self.dRandH = self.hRandH.copy()
        self.dRandH.toType(MEM_DEVICE)

    def randHTearDown(self):
        self.hRandH = None
        self.dRandH = None

    def resultImSetUp(self):
        self.hIm = Image3D(self.grid, MEM_HOST)
        self.dIm = Image3D(self.grid, MEM_DEVICE)

    def resultImTearDown(self):
        self.hIm = None
        self.dIm = None

    def resultFieldSetUp(self):
        self.hField = Field3D(self.grid, MEM_HOST)
        self.dField = Field3D(self.grid, MEM_DEVICE)

    def resultFieldTearDown(self):
        self.hField = None
        self.dField = None

    def randMaskSetUp(self):
        randArr = np.random.rand(self.sz[0], self.sz[1])
        maskArr = np.zeros(randArr.shape)
        maskArr[randArr > 0.5] = 1.0
        self.hRandMask = common.ImFromNPArr(maskArr,
                                            mType=MEM_HOST,
                                            sp=self.imSp)
        self.dRandMask = self.hRandMask.copy()
        self.dRandMask.toType(MEM_DEVICE)

    def randMaskTearDown(self):
        self.hRandMask = None
        self.dRandMask = None


    ################################################################
    #
    # Begin Tests
    #
    ################################################################

    
    #
    # Check image generation
    #
    @PyCATest.AddSetUp(setUpI0, tearDownI0)
    @PyCATest.AddSetUp(setUpI1, tearDownI1)
    def test_ImageGen(self, disp=False):
        self.TestIm(self.hI0,self.dI0,
                    name='I0',disp=disp)
        self.TestIm(self.hI1,self.dI1,
                    name='I1',disp=disp)

    #
    # Gradient
    #
    @PyCATest.AddSetUp(setUpGrad, tearDownGrad)
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(resultFieldSetUp, resultFieldTearDown)
    def test_Gradient(self, disp=False):
        self.TestField(self.hGrad,self.dGrad,
                        name='Gradient',disp=disp)
        Gradient(self.hField, self.hRandIm)
        Gradient(self.dField, self.dRandIm)
        self.TestField(self.hField, self.dField,
                       name='Gradient', disp=disp)

    #
    # Gradient2 (note: only CPU version implemented, should probably
    # just delete this function)
    #
    @PyCATest.AddSetUp(setUpI0, tearDownI0)
    @PyCATest.AddSetUp(setUpGrad, tearDownGrad)
    def test_Gradient2(self, disp=False):
        hGrad2 = Field3D(self.grid, MEM_HOST)
        Gradient2(hGrad2, self.hI0)
        self.TestField(hGrad2,self.dGrad,
                       name='Grad2 v Grad',disp=disp)

    #
    # FluidKernelFFT
    #
    @PyCATest.AddSetUp(setUpDiffOp, tearDownDiffOp)
    @PyCATest.AddSetUp(setUpGrad, tearDownGrad)
    def test_DiffOp(self, disp=False):
        hKGrad = Field3D(self.grid, MEM_HOST)
        dKGrad = Field3D(self.grid, MEM_DEVICE)
        self.hDiffOp.applyInverseOperator(hKGrad, self.hGrad)
        self.dDiffOp.applyInverseOperator(dKGrad, self.dGrad)
        self.TestField(hKGrad, dKGrad,
                       name='Kernel Inverse Op', disp=disp, avEps=2e-6)
        hLKGrad = Field3D(self.grid, MEM_HOST)
        dLKGrad = Field3D(self.grid, MEM_DEVICE)
        self.hDiffOp.applyOperator(hLKGrad, hKGrad)
        self.dDiffOp.applyOperator(dLKGrad, dKGrad)
        self.TestField(hLKGrad, dLKGrad,
                       name='Kernel Forward Op', disp=disp, avEps=3e-6)

        self.TestField(self.hGrad, dLKGrad,
                       name='Diff Op LK = Id', disp=disp, avEps=2e-6)
        self.sz = np.array([127, 121])
        self.TestField(self.hGrad, dLKGrad,
                       name='Diff Op LK = Id', disp=disp, avEps=2e-6)
        
        

    #
    # FiniteDiff
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(resultImSetUp, resultImTearDown)
    def test_FiniteDiff(self, disp=False):
        for dim in [DIM_X, DIM_Y, DIM_Z]:
            for diffType in [DIFF_FORWARD, DIFF_BACKWARD, DIFF_CENTRAL]:
                for bc in [BC_CLAMP, BC_WRAP, BC_APPROX]:
                    FiniteDiff(self.hIm, self.hRandIm, dim, diffType, bc)
                    FiniteDiff(self.dIm, self.dRandIm, dim, diffType, bc)
                    pltname = '%s %s %s'%\
                        (PyCATest.DIMNAMES[dim], 
                         PyCATest.DIFFTNAMES[diffType], 
                         PyCATest.BCNAMES[bc])
                    self.TestIm(self.hIm,self.dIm, 
                                name=pltname,disp=disp)

    #
    # Jacobian Determinant
    #
    @PyCATest.AddSetUp(resultImSetUp, resultImTearDown)
    @PyCATest.AddSetUp(randHSetUp, randHTearDown)
    def test_JacDet(self, disp=False):
        hVx = Field3D(self.grid, MEM_HOST)
        hVy = Field3D(self.grid, MEM_HOST)
        hVz = Field3D(self.grid, MEM_HOST)
        dVx = Field3D(self.grid, MEM_DEVICE)
        dVy = Field3D(self.grid, MEM_DEVICE)
        dVz = Field3D(self.grid, MEM_DEVICE)
        Jacobian(hVx, hVy, hVz, self.hRandH)
        JacDetH(self.hIm, hVx, hVy, hVz)
        Jacobian(dVx, dVy, dVz, self.dRandH)
        JacDetH(self.dIm, dVx, dVy, dVz)

        self.TestField(hVx, dVx, name='Jacobian Vx', disp=disp)
        self.TestField(hVy, dVy, name='Jacobian Vy', disp=disp)
        self.TestField(hVz, dVz, name='Jacobian Vz', disp=disp)

        self.TestIm(self.hIm,self.dIm, 
                    name='jac. det. (full)', disp=disp)

        hJD2 = Image3D(self.grid, MEM_HOST)
        dJD2 = Image3D(self.grid, MEM_DEVICE)

        JacDetH(hJD2, self.hRandH)
        JacDetH(dJD2, self.dRandH)

        self.TestIm(hJD2,dJD2, 
                    name='jac. det. (pointwise)', disp=disp)

        self.TestIm(self.hIm, dJD2, 
                    name='jac. det. (full vs. pointwise)', disp=disp)


    #
    # UpwindGradMag
    #
    @PyCATest.AddSetUp(setUpI0, tearDownI0)
    @PyCATest.AddSetUp(resultImSetUp, resultImTearDown)
    def test_UpwindGradMag(self, disp=False):
        UpwindGradMag(self.hIm, self.hI0, self.hI0)
        UpwindGradMag(self.dIm, self.dI0, self.dI0)
        self.TestIm(self.hIm,self.dIm, 
                    name='UpwindGradMag', disp=disp)

    #
    # Splat (splatting)
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(randHSetUp, randHTearDown)
    @PyCATest.AddSetUp(resultImSetUp, resultImTearDown)
    def test_Splat(self, disp=False):
        Splat(self.hIm, self.hRandH, self.hRandIm)
        Splat(self.dIm, self.dRandH, self.dRandIm)
        self.TestIm(self.hIm,self.dIm, 
                    name='Splat', disp=disp, avEps=2e-6)

    #
    # SplatLargeRange (splatting images with large intensity values > 2048)
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(randHSetUp, randHTearDown)
    @PyCATest.AddSetUp(resultImSetUp, resultImTearDown)
    def test_SplatLargeRange(self, disp=False):
        MulC_I(self.hRandIm, 1000.00)
        MulC_I(self.dRandIm, 1000.00)
        Splat(self.hIm, self.hRandH, self.hRandIm)
        Splat(self.dIm, self.dRandH, self.dRandIm)
        self.TestIm(self.hIm,self.dIm, 
                    name='SplatLargeRange', disp=disp)

    #
    # Ad (Adjoint representation operator)
    #
    @PyCATest.AddSetUp(randVSetUp, randVTearDown)
    @PyCATest.AddSetUp(randHSetUp, randHTearDown)
    @PyCATest.AddSetUp(resultFieldSetUp, resultFieldTearDown)
    def test_Ad(self, disp=False):
        Ad(self.hField, self.hRandH, self.hRandV)
        Ad(self.dField, self.dRandH, self.dRandV)
        self.TestField(self.hField,self.dField, 
                       name='Ad', disp=disp, avEps=2e-6)

    #
    # CoAd (CoAdjoint  operator)
    #
    @PyCATest.AddSetUp(randVSetUp, randVTearDown)
    @PyCATest.AddSetUp(randHSetUp, randHTearDown)
    @PyCATest.AddSetUp(resultFieldSetUp, resultFieldTearDown)
    def test_CoAd(self, disp=False):
        CoAd(self.hField, self.hRandH, self.hRandV)
        CoAd(self.dField, self.dRandH, self.dRandV)
        self.TestField(self.hField,self.dField, 
                       name='CoAd', disp=disp)


    #
    # UpwindDiff
    #
    @PyCATest.AddSetUp(setUpI0, tearDownI0)
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(resultImSetUp, resultImTearDown)
    def test_UpwindDiff(self, disp=False):
        for dim in [DIM_X, DIM_Y]:
            pltname = 'Upwind Diff %s'%PyCATest.DIMNAMES[dim]
            UpwindDiff(self.hIm, self.hI0, self.hRandIm, dim)
            UpwindDiff(self.dIm, self.dI0, self.dRandIm, dim)
        self.TestIm(self.hIm,self.dIm, 
                    name=pltname, disp=disp)

    #
    # Convolve
    #
    @PyCATest.AddSetUp(setUpI0, tearDownI0)
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(resultImSetUp, resultImTearDown)
    def test_Convolve(self, disp=False):

        # create convolution kernel (laplacian stencil)
        kArr = np.zeros([3,3])
        kArr[1,:] = 1.0
        kArr[:,1] = 1.0
        kArr[1,1] = -4.0
        kArr = np.atleast_3d(kArr)

        hK = Image3D(MEM_HOST)
        hK.fromlist(kArr.tolist())
        dK = Image3D(MEM_DEVICE)
        dK.fromlist(kArr.tolist())

        Convolve(self.hIm, self.hI0, hK)
        Convolve(self.dIm, self.dI0, dK)
        self.TestIm(self.hIm,self.dIm, 
                    name='Convolve', disp=disp)

    #
    # SubVol
    #
    @PyCATest.AddSetUp(setUpI0, tearDownI0)
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(resultImSetUp, resultImTearDown)
    def test_SubVol(self, disp=False):

        newSz = (0.66*self.sz).astype('int')
        newSz = Vec3Di(newSz[0], newSz[1], 1)
        newGrid = GridInfo(newSz)
        start = Vec3Di(5,7,0)
        hNewIm = Image3D(newGrid, MEM_HOST)
        dNewIm = Image3D(newGrid, MEM_DEVICE)
        SubVol(hNewIm, self.hRandIm, start)
        SubVol(dNewIm, self.dRandIm, start)
        self.TestIm(hNewIm, dNewIm, 
                    name='SubVol', disp=disp)
        
    #
    # SetSubVol_I
    #
    @PyCATest.AddSetUp(setUpI0, tearDownI0)
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(resultImSetUp, resultImTearDown)
    def test_SetSubVol_I(self, disp=False):

        newSz = (0.66*self.sz).astype('int')
        newSz = Vec3Di(newSz[0], newSz[1], 1)
        newGrid = GridInfo(newSz)
        start = Vec3Di(5,7,0)
        hNewIm = Image3D(newGrid, MEM_HOST)
        SetMem(hNewIm, 4.0)
        dNewIm = Image3D(newGrid, MEM_DEVICE)
        SetMem(dNewIm, 4.0)
        SetSubVol_I(self.hRandIm, hNewIm, start)
        SetSubVol_I(self.dRandIm, dNewIm, start)
        self.TestIm(self.hRandIm, self.dRandIm, 
                    name='SetSubVol_I', disp=disp)
        
    #
    # Resample
    #
    @PyCATest.AddSetUp(setUpI0, tearDownI0)
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(resultImSetUp, resultImTearDown)
    def test_Resample(self, disp=False):

        newSz = (0.66*self.sz).astype('int')
        newSz = Vec3Di(newSz[0], newSz[1], 1)
        newGrid = GridInfo(newSz)
        hNewIm = Image3D(newGrid, MEM_HOST)
        dNewIm = Image3D(newGrid, MEM_DEVICE)
        Resample(hNewIm, self.hRandIm)
        Resample(dNewIm, self.dRandIm)
        self.TestIm(hNewIm, dNewIm, 
                    name='Resample', disp=disp,
                    avEps=4e-6, maxEps=4e-6)

    #
    # ResampleWorld
    #
    @PyCATest.AddSetUp(setUpI0, tearDownI0)
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(resultImSetUp, resultImTearDown)
    def test_ResampleWorld(self, disp=False):

        fac = 0.66
        newSz = (fac*self.sz).astype('int')
        newSp = (1.0/fac)*self.sp
        newSz = Vec3Di(newSz[0], newSz[1], 1)
        newSp = Vec3Df(newSp[0], newSp[1], 1)
        newOr = Vec3Df(1.0, -2.0, 0.0)
        newGrid = GridInfo(newSz, newSp, newOr)
        hNewIm = Image3D(newGrid, MEM_HOST)
        dNewIm = Image3D(newGrid, MEM_DEVICE)
        ResampleWorld(hNewIm, self.hRandIm)
        ResampleWorld(dNewIm, self.dRandIm)
        self.TestIm(hNewIm, dNewIm, 
                    name='ResampleWorld', disp=disp,
                    avEps=4e-6, maxEps=4e-6)

    #
    # ApplyH
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(randHSetUp, randHTearDown)
    @PyCATest.AddSetUp(resultImSetUp, resultImTearDown)
    def test_ApplyH(self, disp=False):
        ApplyH(self.hIm, self.hRandIm, self.hRandH)
        ApplyH(self.dIm, self.dRandIm, self.dRandH)
        self.TestIm(self.hIm,self.dIm, 
                    name='ApplyH', disp=disp)

    #
    # ApplyV
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(randVSetUp, randVTearDown)
    @PyCATest.AddSetUp(resultImSetUp, resultImTearDown)
    def test_ApplyV(self, disp=False):
        ApplyV(self.hIm, self.hRandIm, self.hRandV)
        ApplyV(self.dIm, self.dRandIm, self.dRandV)
        self.TestIm(self.hIm,self.dIm, 
                    name='ApplyV', disp=disp)

    #
    # ApplyVInv
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(randVSetUp, randVTearDown)
    @PyCATest.AddSetUp(resultImSetUp, resultImTearDown)
    def test_ApplyVInv(self, disp=False):
        ApplyVInv(self.hIm, self.hRandIm, self.hRandV)
        ApplyVInv(self.dIm, self.dRandIm, self.dRandV)
        self.TestIm(self.hIm,self.dIm, 
                    name='ApplyVInv', disp=disp)

    #
    # ComposeHH
    #
    @PyCATest.AddSetUp(setUpI0, tearDownI0)
    @PyCATest.AddSetUp(randHSetUp, randHTearDown)
    @PyCATest.AddSetUp(setUpGrad, tearDownGrad)
    @PyCATest.AddSetUp(resultFieldSetUp, resultFieldTearDown)
    def test_ComposeHV(self, disp=False):
        VtoH_I(self.hGrad)
        VtoH_I(self.dGrad)
        ComposeHH(self.hField, self.hRandH, self.hGrad)
        ComposeHH(self.dField, self.dRandH, self.dGrad)
        self.TestField(self.hField,self.dField, 
                       name='ComposeHH', disp=disp)
        hF2 = Field3D(self.grid, MEM_HOST)
        dF2 = Field3D(self.grid, MEM_DEVICE)
        bg = BACKGROUND_STRATEGY_PARTIAL_ID
        ApplyH(hF2, self.hRandH, self.hGrad, bg)
        ApplyH(dF2, self.dRandH, self.dGrad, bg)
        self.TestField(hF2, dF2, 
                        name='ApplyH', disp=disp)
        self.TestField(self.hField, dF2, 
                       name='ComposeHH vs. ApplyH', disp=disp)
    #
    # ComposeHV
    #
    @PyCATest.AddSetUp(setUpI0, tearDownI0)
    @PyCATest.AddSetUp(randHSetUp, randHTearDown)
    @PyCATest.AddSetUp(setUpGrad, tearDownGrad)
    @PyCATest.AddSetUp(resultFieldSetUp, resultFieldTearDown)
    def test_ComposeHV(self, disp=False):
        ComposeHV(self.hField, self.hRandH, self.hGrad)
        ComposeHV(self.dField, self.dRandH, self.dGrad)
        self.TestField(self.hField,self.dField, 
                       name='ComposeHV', disp=disp)
        hF2 = Field3D(self.grid, MEM_HOST)
        dF2 = Field3D(self.grid, MEM_DEVICE)
        bg = BACKGROUND_STRATEGY_PARTIAL_ID
        ApplyV(hF2, self.hRandH, self.hGrad, bg)
        ApplyV(dF2, self.dRandH, self.dGrad, bg)
        self.TestField(hF2, dF2, 
                       name='ApplyV', disp=disp)
        self.TestField(self.hField, dF2, 
                       name='ComposeHV vs. ApplyV', disp=disp)

    #
    # ComposeHVInv
    #
    @PyCATest.AddSetUp(setUpI0, tearDownI0)
    @PyCATest.AddSetUp(randHSetUp, randHTearDown)
    @PyCATest.AddSetUp(setUpGrad, tearDownGrad)
    @PyCATest.AddSetUp(resultFieldSetUp, resultFieldTearDown)
    def test_ComposeHVInv(self, disp=False):
        ComposeHVInv(self.hField, self.hRandH, self.hGrad)
        ComposeHVInv(self.dField, self.dRandH, self.dGrad)
        self.TestField(self.hField,self.dField, 
                       name='ComposeHVInv', disp=disp)
        hF2 = Field3D(self.grid, MEM_HOST)
        dF2 = Field3D(self.grid, MEM_DEVICE)
        bg = BACKGROUND_STRATEGY_PARTIAL_ID
        ApplyVInv(hF2, self.hRandH, self.hGrad, bg)
        ApplyVInv(dF2, self.dRandH, self.dGrad, bg)
        self.TestField(hF2, dF2, 
                       name='ApplyVInv', disp=disp)
        self.TestField(self.hField, dF2, 
                       name='ComposeHVInv vs. ApplyVInv', disp=disp)
        
    #
    # ComposeVH
    #
    @PyCATest.AddSetUp(setUpI0, tearDownI0)
    @PyCATest.AddSetUp(randHSetUp, randHTearDown)
    @PyCATest.AddSetUp(setUpGrad, tearDownGrad)
    @PyCATest.AddSetUp(resultFieldSetUp, resultFieldTearDown)
    def test_ComposeVH(self, disp=False):
        ComposeVH(self.hField, self.hGrad, self.hRandH)
        ComposeVH(self.dField, self.dGrad, self.dRandH)
        self.TestField(self.hField,self.dField, 
                       name='ComposeVH', disp=disp)
        hF2 = Field3D(self.grid, MEM_HOST)
        dF2 = Field3D(self.grid, MEM_DEVICE)
        bg = BACKGROUND_STRATEGY_PARTIAL_ZERO
        ApplyH(hF2, self.hGrad, self.hRandH, bg)
        ApplyH(dF2, self.dGrad, self.dRandH, bg)
        self.TestField(hF2, dF2, 
                       name='ApplyH', disp=disp)
        invSp = self.imSp.inverse()
        MulC_I(hF2, invSp)
        MulC_I(dF2, invSp)
        Add_I(hF2, self.hRandH)
        Add_I(dF2, self.dRandH)
        self.TestField(hF2, self.dField, 
                       name='ApplyH vs. ComposeVH', disp=disp)
        self.TestField(self.hField, dF2, 
                       name='ComposeVH vs. ApplyH', disp=disp)

    #
    # ComposeVInvH
    #
    @PyCATest.AddSetUp(setUpI0, tearDownI0)
    @PyCATest.AddSetUp(randHSetUp, randHTearDown)
    @PyCATest.AddSetUp(setUpGrad, tearDownGrad)
    @PyCATest.AddSetUp(resultFieldSetUp, resultFieldTearDown)
    def test_ComposeVInvH(self, disp=False):
        ComposeVInvH(self.hField, self.hGrad, self.hRandH)
        ComposeVInvH(self.dField, self.dGrad, self.dRandH)
        self.TestField(self.hField,self.dField, 
                       name='ComposeVInvH', disp=disp)
        hF2 = Field3D(self.grid, MEM_HOST)
        dF2 = Field3D(self.grid, MEM_DEVICE)
        bg = BACKGROUND_STRATEGY_PARTIAL_ZERO
        ApplyH(hF2, self.hGrad, self.hRandH, bg)
        ApplyH(dF2, self.dGrad, self.dRandH, bg)
        self.TestField(hF2, dF2, 
                       name='ApplyH', disp=disp)
        invSp = self.imSp.inverse()
        MulC_I(hF2, invSp)
        MulC_I(dF2, invSp)
        Sub(hF2, self.hRandH, hF2)
        Sub(dF2, self.dRandH, dF2)
        self.TestField(self.hField, dF2, 
                       name='ComposeVInvH vs. ApplyH', disp=disp)

    @PyCATest.AddSetUp(randMaskSetUp, randMaskTearDown)
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(resultImSetUp, resultImTearDown)
    def test_CopyMasked(self, disp=False):
        SetMem(self.hIm, 0.0)
        SetMem(self.dIm, 0.0)
        Copy(self.hIm, self.hRandIm, self.hRandMask)
        Copy(self.dIm, self.dRandIm, self.dRandMask)
        self.TestIm(self.hIm,self.dIm, 
                    name='CopyMasked', disp=disp)
        
    @PyCATest.AddSetUp(randMaskSetUp, randMaskTearDown)
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(resultImSetUp, resultImTearDown)
    def test_AbsMasked(self, disp=False):
        SetMem(self.hIm, 0.0)
        SetMem(self.dIm, 0.0)
        Abs(self.hIm, self.hRandIm, self.hRandMask)
        Abs(self.dIm, self.dRandIm, self.dRandMask)
        self.TestIm(self.hIm,self.dIm, 
                    name='AbsMasked', disp=disp)
        
    @PyCATest.AddSetUp(randMaskSetUp, randMaskTearDown)
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(resultImSetUp, resultImTearDown)
    def test_AddMasked(self, disp=False):
        SetMem(self.hIm, 0.0)
        SetMem(self.dIm, 0.0)
        hIm2 = common.RandImage(self.sz, nSig=1.0, gSig=0.0, 
                                mType = MEM_HOST, sp = self.imSp)
        dIm2 = hIm2.copy()
        dIm2.toType(MEM_DEVICE)
        Add(self.hIm, self.hRandIm, hIm2, self.hRandMask)
        Add(self.dIm, self.dRandIm, dIm2, self.dRandMask)
        self.TestIm(self.hIm,self.dIm, 
                    name='AddMasked', disp=disp)


    #
    # ImageField Opers
    #
    @PyCATest.AddSetUp(randVSetUp, randVTearDown)
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(resultFieldSetUp, resultFieldTearDown)
    def test_IF_BinaryOps(self, disp=False):
        Add(self.hField, self.hRandV, self.hRandIm)
        Add(self.dField, self.dRandV, self.dRandIm)
        self.TestField(self.hField, self.dField, 
                       name='IF_Add', disp=disp)
        Sub(self.hField, self.hRandV, self.hRandIm)
        Sub(self.dField, self.dRandV, self.dRandIm)
        self.TestField(self.hField, self.dField, 
                       name='IF_Sub', disp=disp)
        Mul(self.hField, self.hRandV, self.hRandIm)
        Mul(self.dField, self.dRandV, self.dRandIm)
        self.TestField(self.hField, self.dField, 
                       name='IF_Mul', disp=disp)
        Abs_I(self.hRandIm)
        Add_I(self.hRandIm, 1.0)
        Abs_I(self.dRandIm)
        Add_I(self.dRandIm, 1.0)
        Div(self.hField, self.hRandV, self.hRandIm)
        Div(self.dField, self.dRandV, self.dRandIm)
        self.TestField(self.hField, self.dField, 
                       name='IF_Div', disp=disp)


    @PyCATest.AddSetUp(randVSetUp, randVTearDown)
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    def test_IF_Add_I(self, disp=False):
        Add_I(self.hRandV, self.hRandIm)
        Add_I(self.dRandV, self.dRandIm)
        self.TestField(self.hRandV, self.dRandV, 
                       name='IF_Add_I', disp=disp)

    @PyCATest.AddSetUp(randVSetUp, randVTearDown)
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    def test_IF_Sub_I(self, disp=False):
        Sub_I(self.hRandV, self.hRandIm)
        Sub_I(self.dRandV, self.dRandIm)
        self.TestField(self.hRandV, self.dRandV, 
                       name='IF_Add_I', disp=disp)

    @PyCATest.AddSetUp(randVSetUp, randVTearDown)
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    def test_IF_Mul_I(self, disp=False):
        Mul_I(self.hRandV, self.hRandIm)
        Mul_I(self.dRandV, self.dRandIm)
        self.TestField(self.hRandV, self.dRandV, 
                       name='IF_Add_I', disp=disp)

    @PyCATest.AddSetUp(randVSetUp, randVTearDown)
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    def test_IF_Div_I(self, disp=False):
        Abs_I(self.hRandIm)
        Add_I(self.hRandIm, 1.0)
        Abs_I(self.dRandIm)
        Add_I(self.dRandIm, 1.0)
        Div_I(self.hRandV, self.hRandIm)
        Div_I(self.dRandV, self.dRandIm)
        self.TestField(self.hRandV, self.dRandV, 
                       name='IF_Add_I', disp=disp)


        

    @PyCATest.AddSetUp(randVSetUp, randVTearDown)
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(resultFieldSetUp, resultFieldTearDown)
    def test_IF_Add_Mul(self, disp=False):
        hV1, dV1 = self.randVPair()
        Add_Mul(self.hField, self.hRandV, hV1, self.hRandIm)
        Add_Mul(self.dField, self.dRandV, dV1, self.dRandIm)
        self.TestField(self.hField, self.dField, 
                       name='IF_Add_Mul', disp=disp)
        Add_Mul_I(self.hRandV, hV1, self.hRandIm)
        Add_Mul_I(self.dRandV, dV1, self.dRandIm)
        self.TestField(self.hField, self.dRandV, 
                       name='IF_Add_Mul_I', disp=disp)
        self.TestField(self.hRandV, self.dField, 
                       name='IF_Add_Mul_I', disp=disp)

    @PyCATest.AddSetUp(randVSetUp, randVTearDown)
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(resultFieldSetUp, resultFieldTearDown)
    def test_IF_Sub_Mul(self, disp=False):
        hV1, dV1 = self.randVPair()
        Sub_Mul(self.hField, self.hRandV, hV1, self.hRandIm)
        Sub_Mul(self.dField, self.dRandV, dV1, self.dRandIm)
        self.TestField(self.hField, self.dField, 
                       name='IF_Sub_Mul', disp=disp)
        Sub_Mul_I(self.hRandV, hV1, self.hRandIm)
        Sub_Mul_I(self.dRandV, dV1, self.dRandIm)
        self.TestField(self.hField, self.dRandV, 
                       name='IF_Sub_Mul_I', disp=disp)
        self.TestField(self.hRandV, self.dField, 
                       name='IF_Sub_Mul_I', disp=disp)

    @PyCATest.AddSetUp(randVSetUp, randVTearDown)
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(resultFieldSetUp, resultFieldTearDown)
    def test_IF_MulMulC(self, disp=False):
        randNum = np.random.rand()
        MulMulC(self.hField, self.hRandV, self.hRandIm, randNum)
        MulMulC(self.dField, self.dRandV, self.dRandIm, randNum)
        self.TestField(self.hField, self.dField, 
                       name='IF_MulMulC', disp=disp)
        MulMulC_I(self.hRandV, self.hRandIm, randNum)
        MulMulC_I(self.dRandV, self.dRandIm, randNum)
        self.TestField(self.hField, self.dRandV, 
                       name='IF_MulMulC_I', disp=disp)
        self.TestField(self.hRandV, self.dField, 
                       name='IF_MulMulC_I', disp=disp)
        
    @PyCATest.AddSetUp(randVSetUp, randVTearDown)
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(resultFieldSetUp, resultFieldTearDown)
    def test_IF_Add_MulMulC(self, disp=False):
        hV1, dV1 = self.randVPair()
        randNum = np.random.rand()
        Add_MulMulC(self.hField, self.hRandV, hV1, self.hRandIm, randNum)
        Add_MulMulC(self.dField, self.dRandV, dV1, self.dRandIm, randNum)
        self.TestField(self.hField, self.dField, 
                       name='IF_Add_MulMulC', disp=disp)
        Add_MulMulC_I(self.hRandV, hV1, self.hRandIm, randNum)
        Add_MulMulC_I(self.dRandV, dV1, self.dRandIm, randNum)
        self.TestField(self.hField, self.dRandV, 
                       name='IF_Add_MulMulC_I', disp=disp)
        self.TestField(self.hRandV, self.dField, 
                       name='IF_Add_MulMulC_I', disp=disp)

        
    # runTest is only added so that the class can be instantiated
    # directly in order to call individual tests
    def runTest():
        print 'No tests to run directly, all are member functions'

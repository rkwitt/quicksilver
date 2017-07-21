"""
 This file contains testing that adjoint operators are in fact numerically
 adjoint.  All tests can be run from the command line by:

 > python -m unittest discover -v -p '*UnitTest.py'

 To run an individual test with graphical output from ipython:

 import AdjointUnitTest as cgtest
 cgtc = cgtest.MiscTestCase()
 cgtc.test_MemTypeException(disp=True)

"""
#import sys
import unittest
#import PyCATest

import PyCA.Core as ca

import PyCA.Common as common
reload(common)

import numpy as np

# try:
#     import matplotlib.pyplot as plt
#     plt.ion()
# except ImportError:
#     print "Warning: matplotlib.pyplot not found, some functionality disabled"


#
# Test Class
#
class MiscTestCase(unittest.TestCase):

    def __init__(self, methodName='runTest'):
        super(MiscTestCase, self).__init__(methodName)

        self.cudaEnabled = (ca.GetNumberOfCUDADevices() > 0)

        # image size
        self.sz = np.array([128, 120])
        # spacing
        self.sp = [1.5, 2.1]

        # Vec3D versions
        self.imSz = ca.Vec3Di(int(self.sz[0]), int(self.sz[1]), 1)
        self.imSp = ca.Vec3Df(float(self.sp[0]), float(self.sp[1]), 1.0)
        # set up grid
        self.grid = ca.GridInfo(self.imSz, self.imSp)

    def skipIfNoCUDA(self):
        if not self.cudaEnabled:
            self.skipTest('Cannot run test, no CUDA device found or' +
                          'CUDA support not compiled')

    ################################################################
    #
    # Begin Tests
    #
    ################################################################

    #
    # Check that interpolation (via applyH) is adjoint to splatting
    #
    def test_SplatAdjoint(self, disp=False):
        hI = common.RandImage(self.sz, nSig=1.0, gSig=0.0,
                              mType=ca.MEM_HOST, sp=self.imSp)
        hJ = common.RandImage(self.sz, nSig=1.0, gSig=0.0,
                              mType=ca.MEM_HOST, sp=self.imSp)
        hPhi = common.RandField(self.sz, nSig=5.0, gSig=4.0,
                                mType=ca.MEM_HOST, sp=self.imSp)
        tmp = ca.Image3D(self.grid, ca.MEM_HOST)
        # compute < I(Phi(x)), J(x) >
        ca.ApplyH(tmp, hI, hPhi)
        phiIdotJ = ca.Dot(tmp, hJ)
        # compute < I(y), |DPhi^{-1}(y)| J(Phi^{-1}(y)) >
        ca.Splat(tmp, hPhi, hJ)
        IdotphiJ = ca.Dot(tmp, hI)
        #print "a=%f b=%f" % (phiIdotJ, IdotphiJ)
        self.assertLess(abs(phiIdotJ-IdotphiJ), 2e-6)
    #
    # Check that ResampleWorld is adjoint to SplatWorld
    #
    def test_SplatWorldAdjoint(self, disp=False):
        # Random input images
        reg = common.RandImage(self.sz, nSig=1.0, gSig=0.0,
                              mType=ca.MEM_HOST, sp=self.imSp)
        # adjust hJ's grid so that voxels don't align perfectly
        spfactor = 0.77382
        small = common.RandImage(self.sz, nSig=1.0, gSig=0.0,
                              mType=ca.MEM_HOST, sp=self.imSp*spfactor)
        tmpsmall = ca.Image3D(small.grid(), ca.MEM_HOST)
        ca.SetMem(tmpsmall, 0)
        # compute < I(Phi(x)), J(x) >
        ca.ResampleWorld(tmpsmall, reg)
        smallIP = ca.Dot(tmpsmall, small)
        # compute < I(y), |DPhi^{-1}(y)| J(Phi^{-1}(y)) >
        tmpreg = ca.Image3D(reg.grid(), ca.MEM_HOST)
        ca.SetMem(tmpreg, 0)
        ca.SplatWorld(tmpreg, small)
        regIP = ca.Dot(tmpreg, reg)
        #print "a=%f b=%f" % (phiIdotJ, IdotphiJ)
        self.assertLess(abs(smallIP-regIP), 2e-6)

    #
    # Test (co)Adjoint action is numerical adjoint of Adjoint action
    #

    def test_GroupAdjointAction(self, disp=False):
        hM = common.RandField(self.sz, nSig=5.0, gSig=4.0,
                                mType=ca.MEM_HOST, sp=self.imSp)
        hV = common.RandField(self.sz, nSig=5.0, gSig=4.0,
                                mType=ca.MEM_HOST, sp=self.imSp)
        hPhi = common.RandField(self.sz, nSig=5.0, gSig=4.0,
                                mType=ca.MEM_HOST, sp=self.imSp)
        tmp = ca.Field3D(self.grid, ca.MEM_HOST)
        # compute < m, Ad_\phi v >
        ca.Ad(tmp, hPhi, hV)
        rhs = ca.Dot(tmp, hM)
        # compute < Ad^*_\phi m,  v >
        ca.CoAd(tmp, hPhi, hM)
        lhs = ca.Dot(tmp, hV)
        #print "a=%f b=%f" % (rhs, lhs)
        self.assertLess(abs(rhs-lhs), 2e-6)

    #
    # TODO: Test that diffOp is self-adjoint
    #

    # runTest is only added so that the class can be instantiated
    # directly in order to call individual tests
    def runTest():
        print 'No tests to run directly, all are member functions'

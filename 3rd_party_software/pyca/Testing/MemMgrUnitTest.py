#
# This file contains testing of the MemoryManager.
# All tests can be run from the command line by:
#
# > python -m unittest discover -v -p '*UnitTest.py'
#
# To run an individual test with graphical output from ipython:
#
# import MemMgrUnitTest as cgtest
# cgtc = cgtest.MiscTestCase()
# cgtc.test_MemTypeException(disp=True)
#
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
            self.skipTest('Cannot run test, no CUDA device found or CUDA ' +
                          'support not compiled')

    ################################################################
    #
    # Begin Tests
    #
    ################################################################

    def InitMgr(self):
        if not ca.ThreadMemoryManager.isInitialized():
            ca.ThreadMemoryManager.init(self.grid, ca.MEM_HOST, 0)

    # test allocation
    def MemManAlloc(self):
        # this should just create an unmanaged image
        assert(not ca.ThreadMemoryManager.isInitialized())
        I = ca.ManagedImage3D(self.grid, ca.MEM_HOST)
        J = ca.ManagedImage3D(self.grid, ca.MEM_HOST)
        ca.MulC(I, J, 2.0) # test that we can do something with these things
        # Finally, make sure the TMM stays uninitialized
        assert(not ca.ThreadMemoryManager.isInitialized())

    def MemManDeAlloc(self):
        # test that the memory manager releases memory when deleting a managed
        # image
        ca.ThreadMemoryManager.init(self.grid, ca.MEM_HOST, 0)
        tmm = ca.ThreadMemoryManager.instance()
        #nInitialIm = tmm.getNumPools()
        I = ca.ManagedImage3D(self.grid, ca.MEM_HOST)
        nImAfterFirstAlloc = tmm.getNumPools()
        del I
        J = ca.ManagedImage3D(self.grid, ca.MEM_HOST)
        nFinalIm = tmm.getNumPools()
        del J
        # Test that we actually freed I, so that J now resides where I was
        assert(nImAfterFirstAlloc == nFinalIm)

    def test_ThreadMemMan(self):
        # these have to be run in this order
        self.MemManAlloc()
        self.InitMgr()
        self.MemManDeAlloc()

        mType = ca.MEM_HOST
        
        self.assertTrue(ca.ThreadMemoryManager.isInitialized())
            
        tmm = ca.ThreadMemoryManager.instance()

        NImages = 50
        NFields = 50

        for k in range(20):
            
            # print 'k=',k
            
            imList = []
            fieldList = []
            
            for i in range(NImages):
                imList.append(ca.ManagedImage3D(self.grid, mType))
                # print len(imList)
                # print tmm.getNumPools()
                              
            for i in range(NFields):
                fieldList.append(ca.ManagedField3D(self.grid, mType))
                # print len(fieldList)
                # print tmm.getNumPools()

            self.assertEqual(tmm.getNumPools(), NFields*3+NImages)
       
    # runTest is only added so that the class can be instantiated
    # directly in order to call individual tests
    def runTest():
        print 'No tests to run directly, all are member functions'

#
# This file contains testing where PyCA results are compared to
# results from numpy.  All tests can be run from the command line by:
#
# > python -m unittest discover -v -p '*UnitTest.py'
#
# To run an individual test with graphical output from ipython:
#
# import MiscUnitTest as cgtest
# cgtc = cgtest.MiscTestCase()
# cgtc.test_MemTypeException(disp=True)
#
import sys
import unittest
import PyCATest

from PyCA.Core import *

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

        self.cudaEnabled = (GetNumberOfCUDADevices() > 0)
        
        # image size
        self.sz = np.array([128,120])
        # spacing
        self.sp = [1.5, 2.1]

        # Vec3D versions
        self.imSz = Vec3Di(int(self.sz[0]), int(self.sz[1]), 1)
        self.imSp = Vec3Df(float(self.sp[0]), float(self.sp[1]), 1.0)
        # set up grid
        self.grid = GridInfo(self.imSz, self.imSp)

    def skipIfNoCUDA(self):
        if not self.cudaEnabled:
            self.skipTest('Cannot run test, no CUDA device found or CUDA support not compiled')
        
    ################################################################
    #
    # Begin Tests
    #
    ################################################################

    #
    # Check exception wrapping
    #
    @PyCATest.AddSetUp(skipIfNoCUDA)
    def test_MemTypeException(self, disp=False):
        hIm = Image3D(self.grid, MEM_HOST)
        dIm = Image3D(self.grid, MEM_DEVICE)
        with self.assertRaises(Exception):
            MulC(dIm, hIm, 2.0)

    #
    # Check exception wrapping
    #
    @PyCATest.AddSetUp(skipIfNoCUDA)
    def test_MemSizeException(self, disp=False):
        im = Image3D(self.grid, MEM_HOST)
        # set up different grid
        imSz2 = Vec3Di(100, 100, 1)
        grid2 = GridInfo(imSz2, self.imSp)
        im2 = Image3D(self.grid, MEM_DEVICE)
        with self.assertRaises(Exception):
            MulC(im, im2, 2.0)

    # runTest is only added so that the class can be instantiated
    # directly in order to call individual tests
    def runTest():
        print 'No tests to run directly, all are member functions'

#
# This file contains testing where PyCA results are compared to
# results from numpy.  All tests can be run from the command line by:
#
# > python -m unittest discover -v -p '*UnitTest.py'
#
# To run an individual test with graphical output from ipython:
#
# import InterpUnitTest as cgtest
# cgtc = cgtest.InterpTestCase()
# cgtc.test_ResampleInterp(disp=True)
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
class InterpTestCase(unittest.TestCase):

    def __init__(self, methodName='runTest'):
        super(InterpTestCase, self).__init__(methodName)

        self.cudaEnabled = (GetNumberOfCUDADevices() > 0)

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
    def test_ResampleInterp(self, disp=False):
        # generate small integer-valued image
        initMax = 5
        randArrSmall = (np.random.rand(10,10)*initMax).astype(int)
        randImSmall = common.ImFromNPArr(randArrSmall)
        imLarge = Image3D(50,50,1)
        Resample(imLarge, randImSmall,
                 BACKGROUND_STRATEGY_CLAMP,
                 INTERP_NN)
        nUnique = len(np.unique(imLarge.asnp()))
        self.assertEqual(nUnique,initMax)
        

    # runTest is only added so that the class can be instantiated
    # directly in order to call individual tests
    def runTest():
        print 'No tests to run directly, all are member functions'

if __name__ == '__main__':
    """
    Run a test showing different interpolation methods used for
    upsampling and deformation.
    """
    
    import PyCA.Core as ca
    import PyCA.Common as common
    import PyCA.Display as display

    import numpy as np
    import matplotlib.pyplot as plt

    plt.ion()
    
    initMax = 5
    randArrSmall = (np.random.rand(10,10)*initMax).astype(int)
    imSmall = common.ImFromNPArr(randArrSmall)
    imLargeNN = ca.Image3D(50,50,1)
    imLargeLinear = ca.Image3D(50,50,1)
    imLargeCubic = ca.Image3D(50,50,1)
    ca.Resample(imLargeNN, imSmall,
                ca.BACKGROUND_STRATEGY_CLAMP,
                ca.INTERP_NN)
    ca.Resample(imLargeLinear, imSmall,
                ca.BACKGROUND_STRATEGY_CLAMP,
                ca.INTERP_LINEAR)
    ca.Resample(imLargeCubic, imSmall,
                ca.BACKGROUND_STRATEGY_CLAMP,
                ca.INTERP_CUBIC)
    plt.figure('interp test')
    plt.subplot(2,3,1)
    display.DispImage(imLargeNN, 'NN', newFig=False)
    plt.subplot(2,3,2)
    display.DispImage(imLargeLinear, 'Linear', newFig=False)
    plt.subplot(2,3,3)
    display.DispImage(imLargeCubic, 'Cubic', newFig=False)
    plt.subplot(2,3,5)
    display.DispImage(imSmall, 'small', newFig=False)
    plt.show()


    h = common.WavyDef([50,50], nWaves=1, waveAmp=10, waveDim=0,
                       mType=ca.MEM_HOST, deformation=True)
    imDefNN = imLargeNN.copy()
    ca.ApplyH(imDefNN, imLargeNN, h, ca.BACKGROUND_STRATEGY_CLAMP, ca.INTERP_NN)
    imDefLinear = imLargeNN.copy()
    ca.ApplyH(imDefLinear, imLargeNN, h, ca.BACKGROUND_STRATEGY_CLAMP, ca.INTERP_LINEAR)
    imDefCubic = imLargeNN.copy()
    ca.ApplyH(imDefCubic, imLargeNN, h, ca.BACKGROUND_STRATEGY_CLAMP, ca.INTERP_CUBIC)
    plt.figure('interp def test')
    plt.subplot(2,3,1)
    display.DispImage(imDefNN, 'NN', newFig=False)
    plt.subplot(2,3,2)
    display.DispImage(imDefLinear, 'Linear', newFig=False)
    plt.subplot(2,3,3)
    display.DispImage(imDefCubic, 'Cubic', newFig=False)
    plt.subplot(2,3,5)
    display.DispImage(imLargeNN, 'orig', newFig=False)
    plt.show()

    # make sure plots don't close on exit
    plt.ioff()
    plt.show()

#
# This file contains testing where PyCA results are compared to
# results from numpy.  All tests can be run from the command line by:
#
# > python -m unittest discover -v -p '*UnitTest.py'
#
# To run an individual test with graphical output from ipython:
#
# import NumpyUnitTest as nptest
# ntc = nptest.NumpyTestCase()
# ntc.test_Exp(disp=True)
#
import sys
import unittest
import PyCATest

from PyCA.Core import *

import PyCA.Common as common
reload(common)

# numpy reference implementations of PyCA functions
import PyCA.Numpy as canp
reload(canp)

import numpy as np

try:
    import matplotlib.pyplot as plt
    plt.ion()
except ImportError:
    print "Warning: matplotlib.pyplot not found, some functionality disabled"

def CheckIm(im, arr, name, disp):
    him = im.copy()
    him.toType(MEM_HOST)
    imArr = np.squeeze(him.asnp())
    diff = imArr-np.squeeze(arr)
    if disp:
        title = name
        plt.figure(title)
        plt.clf()
        plt.subplot(1,3,1)
        plt.imshow(imArr)
        plt.colorbar();
        plt.title('AW')
        plt.draw()
        plt.subplot(1,3,2)
        plt.imshow(np.squeeze(arr))
        plt.colorbar();
        plt.title('np')
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

def CheckField(v, varr, name, disp):
    hv = v.copy()
    hv.toType(MEM_HOST)
    (vx,vy,vz) = hv.asnp()
    vx = np.squeeze(vx)
    vy = np.squeeze(vy)
    diff_x = vx-varr[:,:,0]
    diff_y = vy-varr[:,:,1]
    if disp:
        plt.figure(name)

        plt.subplot(2,3,1)
        plt.imshow(vx)
        plt.colorbar();
        plt.title('aw x')
        plt.draw()
        plt.subplot(2,3,2)
        plt.imshow(varr[:,:,0])
        plt.colorbar();
        plt.title('np x')
        plt.draw()
        plt.subplot(2,3,3)
        plt.imshow(diff_x)
        plt.colorbar();
        plt.title('diff x')
        plt.draw()

        plt.subplot(2,3,4)
        plt.imshow(vy)
        plt.colorbar();
        plt.title('aw y')
        plt.draw()
        plt.subplot(2,3,5)
        plt.imshow(varr[:,:,1])
        plt.colorbar();
        plt.title('np y')
        plt.draw()
        plt.subplot(2,3,6)
        plt.imshow(diff_y)
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
class NumpyTestCase(unittest.TestCase):

    def __init__(self, methodName='runTest'):
        super(NumpyTestCase, self).__init__(methodName)

        # memory type
        if GetNumberOfCUDADevices() > 0:
            self.mType = MEM_DEVICE
        else:
            self.mType = MEM_HOST

        # allowable average abs. diff
        self.AvEps = 1e-6
        # allowable max abs. diff
        self.MaxEps = 1e-4
        # image size
        self.sz = np.array([31,33])
        # spacing
        self.imSp = Vec3Df(1.5, 2.1, 1.0)

        self.vsz = np.append(self.sz, 2)
        self.imSz = Vec3Di(int(self.sz[0]),int(self.sz[1]),1)
        # set up grid
        self.grid = GridInfo(self.imSz, self.imSp)

        # minimum value if NonNeg=True
        self.minNonNeg = 1e-2

        # range of ints for testing
        self.minRandInt = -1024
        self.maxRandInt = 1024

        self.BinaryOperTests = [
            #
            #(arg1Type, arg2Type, outType, pycaOp, npOp, inPlace)
            #
            # add
            #
            ('Image3D', 'Image3D', 'Image3D', 'Add', '__add__', False),
            ('Image3D', 'Image3D', 'Image3D', '__add__', '__add__', False),
            ('Image3D', 'Image3D', 'Image3D', 'Add_I', '__iadd__', True),
            ('Image3D', 'Image3D', 'Image3D', '__iadd__', '__iadd__', True),
            ('Image3D', 'float', 'Image3D', 'AddC', '__add__', False),
            ('Image3D', 'float', 'Image3D', '__add__', '__add__', False),
            ('Image3D', 'float', 'Image3D', 'AddC_I', '__iadd__', True),
            ('Image3D', 'float', 'Image3D', '__iadd__', '__iadd__', True),
            ('Image3D', 'float', 'Image3D', '__radd__', '__radd__', False),
            ('Field3D', 'Field3D', 'Field3D', 'Add', '__add__', False),
            ('Field3D', 'Field3D', 'Field3D', '__add__', '__add__', False),
            ('Field3D', 'Field3D', 'Field3D', 'Add_I', '__iadd__', True),
            ('Field3D', 'Field3D', 'Field3D', '__iadd__', '__iadd__', True),
            ('Field3D', 'Vec3Df', 'Field3D', 'AddC', '__add__', False),
            ('Field3D', 'Vec3Df', 'Field3D', '__add__', '__add__', False),
            ('Field3D', 'Vec3Df', 'Field3D', 'AddC_I', '__iadd__', True),
            ('Field3D', 'Vec3Df', 'Field3D', '__iadd__', '__iadd__', True),
            ('Field3D', 'Vec3Df', 'Field3D', '__radd__', '__radd__', False),
            ('Field3D', 'float', 'Field3D', 'AddC', '__add__', False),
            ('Field3D', 'float', 'Field3D', '__add__', '__add__', False),
            ('Field3D', 'float', 'Field3D', 'AddC_I', '__iadd__', True),
            ('Field3D', 'float', 'Field3D', '__iadd__', '__iadd__', True),
            ('Field3D', 'float', 'Field3D', '__radd__', '__radd__', False),
            ('Vec3Df', 'Vec3Df', 'Vec3Df', '__add__', '__add__', False),
            ('Vec3Df', 'Vec3Df', 'Vec3Df', '__iadd__', '__iadd__', True),
            ('Vec3Df', 'Vec3Di', 'Vec3Df', '__add__', '__add__', False),
            ('Vec3Df', 'Vec3Di', 'Vec3Df', '__iadd__', '__iadd__', True),
            ('Vec3Df', 'float', 'Vec3Df', '__add__', '__add__', False),
            # not implemented, but probably should be
            # ('Vec3Df', 'float', 'Vec3Df', '__radd__', '__radd__', False),
            ('Vec3Df', 'float', 'Vec3Df', '__iadd__', '__iadd__', True),
            ('Vec3Df', 'int', 'Vec3Df', '__add__', '__add__', False),
            # not implemented, but probably should be
            # ('Vec3Df', 'int', 'Vec3Df', '__radd__', '__radd__', False),
            ('Vec3Df', 'int', 'Vec3Df', '__iadd__', '__iadd__', True),
            ('Vec3Di', 'Vec3Di', 'Vec3Di', '__add__', '__add__', False),
            ('Vec3Di', 'Vec3Di', 'Vec3Di', '__iadd__', '__iadd__', True),
            ('Vec3Di', 'int', 'Vec3Di', '__add__', '__add__', False),
            # not implemented, but probably should be
            # ('Vec3Di', 'int', 'Vec3Di', '__radd__', '__radd__', False),
            ('Vec3Di', 'int', 'Vec3Di', '__iadd__', '__iadd__', True),
            # not implemented
            # ('Vec3Di', 'Vec3Df', 'Vec3Df', '__radd__', '__radd__', False),

            #
            # sub
            #
            ('Image3D', 'Image3D', 'Image3D', 'Sub', '__sub__', False),
            ('Image3D', 'Image3D', 'Image3D', '__sub__', '__sub__', False),
            ('Image3D', 'Image3D', 'Image3D', 'Sub_I', '__isub__', True),
            ('Image3D', 'Image3D', 'Image3D', '__isub__', '__isub__', True),
            ('Image3D', 'float', 'Image3D', 'SubC', '__sub__', False),
            ('Image3D', 'float', 'Image3D', '__sub__', '__sub__', False),
            ('Image3D', 'float', 'Image3D', 'SubC_I', '__isub__', True),
            ('Image3D', 'float', 'Image3D', '__isub__', '__isub__', True),
            ('Field3D', 'Field3D', 'Field3D', 'Sub', '__sub__', False),
            ('Field3D', 'Field3D', 'Field3D', '__sub__', '__sub__', False),
            ('Field3D', 'Field3D', 'Field3D', 'Sub_I', '__isub__', True),
            ('Field3D', 'Field3D', 'Field3D', '__isub__', '__isub__', True),
            ('Field3D', 'Vec3Df', 'Field3D', 'SubC', '__sub__', False),
            ('Field3D', 'Vec3Df', 'Field3D', '__sub__', '__sub__', False),
            ('Field3D', 'Vec3Df', 'Field3D', 'SubC_I', '__isub__', True),
            ('Field3D', 'Vec3Df', 'Field3D', '__isub__', '__isub__', True),
            ('Field3D', 'float', 'Field3D', 'SubC', '__sub__', False),
            ('Field3D', 'float', 'Field3D', '__sub__', '__sub__', False),
            ('Field3D', 'float', 'Field3D', 'SubC_I', '__isub__', True),
            ('Field3D', 'float', 'Field3D', '__isub__', '__isub__', True),
            ('Vec3Df', 'Vec3Df', 'Vec3Df', '__sub__', '__sub__', False),
            ('Vec3Df', 'Vec3Df', 'Vec3Df', '__isub__', '__isub__', True),
            ('Vec3Df', 'Vec3Di', 'Vec3Df', '__sub__', '__sub__', False),
            ('Vec3Df', 'Vec3Di', 'Vec3Df', '__isub__', '__isub__', True),
            ('Vec3Df', 'float', 'Vec3Df', '__sub__', '__sub__', False),
            ('Vec3Df', 'float', 'Vec3Df', '__isub__', '__isub__', True),
            ('Vec3Df', 'int', 'Vec3Df', '__sub__', '__sub__', False),
            ('Vec3Df', 'int', 'Vec3Df', '__isub__', '__isub__', True),
            ('Vec3Di', 'Vec3Di', 'Vec3Di', '__sub__', '__sub__', False),
            ('Vec3Di', 'Vec3Di', 'Vec3Di', '__isub__', '__isub__', True),
            ('Vec3Di', 'int', 'Vec3Di', '__sub__', '__sub__', False),
            ('Vec3Di', 'int', 'Vec3Di', '__isub__', '__isub__', True),
            #
            # mul
            #
            ('Image3D', 'Image3D', 'Image3D', 'Mul', '__mul__', False),
            ('Image3D', 'Image3D', 'Image3D', '__mul__', '__mul__', False),
            ('Image3D', 'Image3D', 'Image3D', 'Mul_I', '__imul__', True),
            ('Image3D', 'Image3D', 'Image3D', '__imul__', '__imul__', True),
            ('Image3D', 'float', 'Image3D', 'MulC', '__mul__', False),
            ('Image3D', 'float', 'Image3D', '__mul__', '__mul__', False),
            ('Image3D', 'float', 'Image3D', '__rmul__', '__rmul__', False),
            ('Image3D', 'float', 'Image3D', 'MulC_I', '__imul__', True),
            ('Image3D', 'float', 'Image3D', '__imul__', '__imul__', True),
            ('Field3D', 'Image3D', 'Field3D', 'Mul', '__mul__', False),
            ('Field3D', 'Image3D', 'Field3D', '__mul__', '__mul__', False),
            ('Field3D', 'Image3D', 'Field3D', 'Mul_I', '__imul__', True),
            ('Field3D', 'Image3D', 'Field3D', '__imul__', '__imul__', True),
            ('Field3D', 'Vec3Df', 'Field3D', 'MulC', '__mul__', False),
            ('Field3D', 'Vec3Df', 'Field3D', '__mul__', '__mul__', False),
            ('Field3D', 'Vec3Df', 'Field3D', '__rmul__', '__rmul__', False),
            ('Field3D', 'Vec3Df', 'Field3D', 'MulC_I', '__imul__', True),
            ('Field3D', 'Vec3Df', 'Field3D', '__imul__', '__imul__', True),
            ('Field3D', 'float', 'Field3D', 'MulC', '__mul__', False),
            ('Field3D', 'float', 'Field3D', '__mul__', '__mul__', False),
            ('Field3D', 'float', 'Field3D', '__rmul__', '__rmul__', False),
            ('Field3D', 'float', 'Field3D', 'MulC_I', '__imul__', True),
            ('Field3D', 'float', 'Field3D', '__imul__', '__imul__', True),
            ('Vec3Df', 'Vec3Df', 'Vec3Df', '__mul__', '__mul__', False),
            ('Vec3Df', 'Vec3Df', 'Vec3Df', '__imul__', '__imul__', True),
            ('Vec3Df', 'Vec3Di', 'Vec3Df', '__mul__', '__mul__', False),
            ('Vec3Df', 'Vec3Di', 'Vec3Df', '__imul__', '__imul__', True),
            ('Vec3Df', 'float', 'Vec3Df', '__mul__', '__mul__', False),
            ('Vec3Df', 'float', 'Vec3Df', '__imul__', '__imul__', True),
            ('Vec3Df', 'int', 'Vec3Df', '__mul__', '__mul__', False),
            ('Vec3Df', 'int', 'Vec3Df', '__imul__', '__imul__', True),
            ('Vec3Di', 'Vec3Di', 'Vec3Di', '__mul__', '__mul__', False),
            ('Vec3Di', 'Vec3Di', 'Vec3Di', '__imul__', '__imul__', True),
            ('Vec3Di', 'int', 'Vec3Di', '__mul__', '__mul__', False),
            ('Vec3Di', 'int', 'Vec3Di', '__imul__', '__imul__', True),
            #
            # div
            #
            ('Image3D', 'Image3D', 'Image3D', 'Div', '__div__', False),
            ('Image3D', 'Image3D', 'Image3D', '__div__', '__div__', False),
            ('Image3D', 'Image3D', 'Image3D', 'Div_I', '__idiv__', True),
            ('Image3D', 'Image3D', 'Image3D', '__idiv__', '__idiv__', True),
            ('Image3D', 'float', 'Image3D', 'DivC', '__div__', False),
            ('Image3D', 'float', 'Image3D', '__div__', '__div__', False),
            ('Image3D', 'float', 'Image3D', 'DivC_I', '__idiv__', True),
            ('Image3D', 'float', 'Image3D', '__idiv__', '__idiv__', True),
            ('Field3D', 'Image3D', 'Field3D', 'Div', '__div__', False),
            ('Field3D', 'Image3D', 'Field3D', '__div__', '__div__', False),
            ('Field3D', 'Image3D', 'Field3D', 'Div_I', '__idiv__', True),
            ('Field3D', 'Image3D', 'Field3D', '__idiv__', '__idiv__', True),
            ('Field3D', 'Vec3Df', 'Field3D', 'DivC', '__div__', False),
            ('Field3D', 'Vec3Df', 'Field3D', '__div__', '__div__', False),
            ('Field3D', 'Vec3Df', 'Field3D', 'DivC_I', '__idiv__', True),
            ('Field3D', 'Vec3Df', 'Field3D', '__idiv__', '__idiv__', True),
            ('Field3D', 'float', 'Field3D', 'DivC', '__div__', False),
            ('Field3D', 'float', 'Field3D', '__div__', '__div__', False),
            ('Field3D', 'float', 'Field3D', 'DivC_I', '__idiv__', True),
            ('Field3D', 'float', 'Field3D', '__idiv__', '__idiv__', True),
            ('Vec3Df', 'Vec3Df', 'Vec3Df', '__div__', '__div__', False),
            ('Vec3Df', 'Vec3Df', 'Vec3Df', '__idiv__', '__idiv__', True),
            ('Vec3Df', 'Vec3Di', 'Vec3Df', '__div__', '__div__', False),
            ('Vec3Df', 'Vec3Di', 'Vec3Df', '__idiv__', '__idiv__', True),
            ('Vec3Df', 'float', 'Vec3Df', '__div__', '__div__', False),
            ('Vec3Df', 'float', 'Vec3Df', '__idiv__', '__idiv__', True),
            ('Vec3Df', 'int', 'Vec3Df', '__div__', '__div__', False),
            ('Vec3Df', 'int', 'Vec3Df', '__idiv__', '__idiv__', True),
            ('Vec3Di', 'Vec3Di', 'Vec3Di', '__div__', '__div__', False),
            ('Vec3Di', 'Vec3Di', 'Vec3Di', '__idiv__', '__idiv__', True),
            ('Vec3Di', 'int', 'Vec3Di', '__div__', '__div__', False),
            ('Vec3Di', 'int', 'Vec3Di', '__idiv__', '__idiv__', True),
            #
            # min/max
            #
            ('Image3D', 'Image3D', 'Image3D', 'Min', 'np.minimum', False),
            ('Image3D', 'float', 'Image3D', 'MinC', 'np.minimum', False),
            ('Image3D', 'Image3D', 'Image3D', 'Max', 'np.maximum', False),
            ('Image3D', 'float', 'Image3D', 'MaxC', 'np.maximum', False),
            #
            # Atan2
            #
            ('Image3D', 'Image3D', 'Image3D', 'Atan2', 'np.arctan2', False),
            ('Image3D', 'float', 'Image3D', 'Atan2', 'np.arctan2', False),
            #
            # Logic Operations
            #
            ('Image3D', 'Image3D', 'Image3D', 'GT', '__gt__', False),
            ('Image3D', 'float', 'Image3D', 'GT', '__gt__', False),
            ('Image3D', 'Image3D', 'Image3D', 'GTE', '__ge__', False),
            ('Image3D', 'float', 'Image3D', 'GTE', '__ge__', False),
            ('Image3D', 'Image3D', 'Image3D', 'EQ', '__eq__', False),
            ('Image3D', 'float', 'Image3D', 'EQ', '__eq__', False),
            ('Image3D', 'Image3D', 'Image3D', 'NEQ', '__ne__', False),
            ('Image3D', 'float', 'Image3D', 'NEQ', '__ne__', False),
            ('Image3D', 'Image3D', 'Image3D', 'LT', '__lt__', False),
            ('Image3D', 'float', 'Image3D', 'LT', '__lt__', False),
            ('Image3D', 'Image3D', 'Image3D', 'LTE', '__le__', False),
            ('Image3D', 'float', 'Image3D', 'LTE', '__le__', False)
            ]


    def TestIm(self, pycaI, npI, name, disp, avEps=None, maxEps=None):

        if avEps is None:
            avEps = self.AvEps
        if maxEps is None:
            maxEps = self.MaxEps

        diffAv, diffMax = CheckIm(pycaI, npI, name=name, disp=disp)

        self.assertLess(diffAv, avEps)
        self.assertLess(diffMax, maxEps)

    def TestField(self, pycaF, npF, name, disp, avEps=None, maxEps=None):

        if avEps is None:
            avEps = self.AvEps
        if maxEps is None:
            maxEps = self.MaxEps

        diffAv, diffMax = CheckField(pycaF, npF, name=name, disp=disp)

        self.assertLess(diffAv, avEps)
        self.assertLess(diffMax, self.MaxEps)


    def genRandImPair(self, NonNeg=False):
        randIm = common.RandImage(self.sz, \
                                  MEM_HOST, \
                                  sp=self.imSp, \
                                  NonNeg=NonNeg)
        if NonNeg:
            AddC_I(randIm, self.minNonNeg)

        randArr = randIm.asnp().copy()
        randIm.toType(self.mType)

        return (randArr, randIm)

    def genRandImUnifPair(self, lbound=-1, ubound=1):
        randIm = common.RandUnifImage(self.sz, lbound, ubound,
                                      MEM_HOST,
                                      sp=self.imSp)

        randArr = randIm.asnp().copy()
        randIm.toType(self.mType)

        return (randArr, randIm)

    def randImSetUp(self):
        (self.randArr, self.randIm) = self.genRandImPair()

    def randImUnifSetUp(self):
        (self.randArr, self.randIm) = self.genRandImUnifPair()

    def randImTearDown(self):
        self.randArr = None
        self.randIm = None

    def genRandFieldPair(self, NonNeg=False):
        randF = common.RandField(self.sz, \
                                 MEM_HOST, \
                                 sp=self.imSp, \
                                 NonNeg=NonNeg)
        if NonNeg:
            AddC_I(randF, self.minNonNeg)

        randFArr = np.zeros(self.vsz)
        (vx,vy,vz) = randF.asnp()
        randFArr[:,:,0] = np.squeeze(vx)
        randFArr[:,:,1] = np.squeeze(vy)
        randF.toType(self.mType)

        return (randFArr, randF)

    def randFieldSetUp(self):
        (self.randFArr, self.randF) = self.genRandFieldPair()

    def randFieldTearDown(self):
        self.randFArr = None
        self.randF = None

    def genRandVecfPair(self, NonNeg=False):
        randVecArr = np.random.rand(3)
        if NonNeg:
            randVecArr = abs(randVecArr) + self.minNonNeg
        randVec = Vec3Df(float(randVecArr[0]),
                         float(randVecArr[1]),
                         float(randVecArr[2]))
        return (randVecArr, randVec)

    def genRandVeciPair(self, NonNeg=False):
        randVecArr = np.random.random_integers(self.minRandInt,
                                               self.maxRandInt,
                                               3)
        if NonNeg:
            randVecArr = np.maximum(randVecArr,1)
        randVec = Vec3Di(int(randVecArr[0]),
                         int(randVecArr[1]),
                         int(randVecArr[2]))
        return (randVecArr, randVec)

    def resultImSetUp(self):
        self.im = Image3D(self.grid, self.mType)

    def resultImTearDown(self):
        self.im = None

    def resultFieldSetUp(self):
        self.field = Field3D(self.grid, self.mType)

    def resultFieldTearDown(self):
        self.field = None

    def maskSetUp(self):
        bdr=10
        self.maskArr = np.zeros(self.sz)
        self.maskArr[bdr:-bdr,bdr:-bdr] = 1
        self.mask = common.ImFromNPArr(self.maskArr, \
                                           mType=self.mType, \
                                           sp=self.imSp)

    def maskTearDown(self):
        self.maskArr = None
        self.mask = None


    def randMaskSetUp(self):
        randArr = np.random.rand(self.sz[0], self.sz[1])
        self.randMaskArr = np.zeros(randArr.shape)
        self.randMaskArr[randArr > 0.5] = 1.0
        self.randMask = common.ImFromNPArr(self.randMaskArr,
                                           mType=self.mType,
                                           sp=self.imSp)
    def randMaskTearDown(self):
        self.randMaskArr = None
        self.randMask = None

    ################################################################
    #
    # Begin Tests
    #
    ################################################################


    #
    # Sgn
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    def test_Sgn(self, disp=False):
        Sgn_I(self.randIm)
        arr = np.sign(self.randArr)
        self.TestIm(self.randIm, arr, \
                        name='Sgn', disp=disp)

    #
    # MinMax
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    def test_MinMax(self, disp=False):
        imMax = Max(self.randIm)
        arrMax = np.max(self.randArr)
        self.assertLess(np.abs(imMax-arrMax), self.AvEps)
        imMin = Min(self.randIm)
        arrMin = np.min(self.randArr)
        self.assertLess(np.abs(imMin-arrMin), self.AvEps)
        imMinMaxMin,imMinMaxMax = MinMax(self.randIm)
        self.assertLess(np.abs(imMinMaxMax-arrMax), self.AvEps)
        self.assertLess(np.abs(imMinMaxMin-arrMin), self.AvEps)

    #
    # Sum (Reduction)
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    def test_Sum(self, disp=False):
        imSum = Sum(self.randIm)
        arrSum = np.sum(self.randArr)
        eps = 1e-5*abs(arrSum)
        if disp:
            print 'pyca Sum: %f, numpy Sum: %f'%(imSum, arrSum)
        self.assertLess(np.abs(imSum-arrSum), eps)

    #
    # Sum2 (Reduction)
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    def test_Sum2(self, disp=False):
        imSum2 = Sum2(self.randIm)
        arrSum2 = np.sum(self.randArr**2)
        if disp:
            print 'pyca Sum2: %f, numpy Sum2: %f'%(imSum2, arrSum2)
        eps = 1e-5*abs(arrSum2)
        self.assertLess(np.abs(imSum2-arrSum2), eps)

    #
    # Exp
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    def test_Exp(self, disp=False):
        Exp_I(self.randIm)
        arr = np.exp(self.randArr)
        self.TestIm(self.randIm, arr, \
                        name='Exp', disp=disp)


    #
    # Pow
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    def test_Pow(self, disp=False):
        # need nonnegative numbers
        Abs_I(self.randIm)
        self.randArr = np.abs(self.randArr)
        # run test
        PowC_I(self.randIm, 0.75)
        arr = np.power(self.randArr, 0.75)
        self.TestIm(self.randIm, arr, \
                        name='Pow', disp=disp)


    #
    # Sin
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    def test_Sin(self, disp=False):
        # run test
        Sin_I(self.randIm)
        arr = np.sin(self.randArr)
        self.TestIm(self.randIm, arr,
                    name='Sin', disp=disp)

    #
    # Asin
    #
    @PyCATest.AddSetUp(randImUnifSetUp, randImTearDown)
    def test_Asin(self, disp=False):
        # run test
        Asin_I(self.randIm)
        arr = np.arcsin(self.randArr)
        self.TestIm(self.randIm, arr,
                    name='Asin', disp=disp)

    #
    # Cos
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    def test_Cos(self, disp=False):
        # run test
        Cos_I(self.randIm)
        arr = np.cos(self.randArr)
        self.TestIm(self.randIm, arr,
                    name='Cos', disp=disp)

    #
    # Acos
    #
    @PyCATest.AddSetUp(randImUnifSetUp, randImTearDown)
    def test_Acos(self, disp=False):
        # run test
        Acos_I(self.randIm)
        arr = np.arccos(self.randArr)
        self.TestIm(self.randIm, arr,
                    name='Acos', disp=disp)

    #
    # Tan
    #
    @PyCATest.AddSetUp(randImUnifSetUp, randImTearDown)
    def test_Tan(self, disp=False):
        # run test
        Tan_I(self.randIm)
        arr = np.tan(self.randArr)
        self.TestIm(self.randIm, arr,
                    name='Tan', disp=disp)

    #
    # Atan
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    def test_Atan(self, disp=False):
        # run test
        Atan_I(self.randIm)
        arr = np.arctan(self.randArr)
        self.TestIm(self.randIm, arr,
                    name='Atan', disp=disp)

    #
    # Csc
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    def test_Csc(self, disp=False):
        # run test
        Csc_I(self.randIm)
        arr = np.sin(self.randArr)
        arr = 1.0/arr
        self.TestIm(self.randIm, arr,
                    name='Csc', disp=disp)

    #
    # Sec
    #
    @PyCATest.AddSetUp(randImUnifSetUp, randImTearDown)
    def test_Sec(self, disp=False):
        # run test
        Sec_I(self.randIm)
        arr = np.cos(self.randArr)
        arr = 1.0/arr
        self.TestIm(self.randIm, arr,
                    name='Sec', disp=disp)

    #
    # Cot
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    def test_Cot(self, disp=False):
        # run test
        Cot_I(self.randIm)
        arr = np.tan(self.randArr)
        arr = 1.0/arr
        self.TestIm(self.randIm, arr,
                    name='Cot', disp=disp)


    #
    # SoftAbs
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    def test_SoftAbs(self, disp=False):
        eps = 0.5
        arr = self.randArr.copy()
        arr2 = self.randArr.copy()
        arr2[arr<-eps] = -arr[arr<-eps]-eps/2.0
        arr2[arr>eps] = arr[arr>eps]-(eps/2.0)
        vals = arr[np.abs(arr)<=eps]
        arr2[np.abs(arr)<=eps] = (vals*vals)/(2*eps)
        SoftAbs_I(self.randIm,eps)

        self.TestIm(self.randIm, arr2, \
                        'SoftAbs', disp=disp)

        if disp:
            x = np.arange(-eps*5, eps*5, eps/10.0)
            xim = common.ImFromNPArr(x, self.mType)
            yim = xim.copy()
            SoftAbs_I(yim,eps)
            plt.figure('SoftAbsFunc')
            xim.toType(MEM_HOST)
            yim.toType(MEM_HOST)
            plt.plot(np.squeeze(xim.asnp()), np.squeeze(yim.asnp()), 'b')
            plt.plot(np.squeeze(xim.asnp()), np.abs(np.squeeze(xim.asnp())),'r')


    #
    # Copy Masked
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(resultImSetUp, resultImTearDown)
    @PyCATest.AddSetUp(randMaskSetUp, randMaskTearDown)
    def test_CopyMasked(self, disp=False):
        SetMem(self.im, 0.0)
        Copy(self.im, self.randIm, self.randMask)
        arr = np.zeros(self.randArr.shape)
        arr[self.randMaskArr != 0] = self.randArr[self.randMaskArr != 0]
        self.TestIm(self.im, arr, 
                    name='Copy_Masked', disp=disp)

    #
    # Abs Masked
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(randMaskSetUp, randMaskTearDown)
    def test_AbsMasked(self, disp=False):
        im = self.randIm.copy()
        Abs_I(im, self.randMask)
        arr = self.randArr.copy()
        arr[self.randMaskArr != 0] = np.abs(arr[self.randMaskArr != 0])
        self.TestIm(im, arr, 
                    name='Abs_I_Masked', disp=disp)
        arr, im = self.genRandImPair()
        Abs(im, self.randIm, self.randMask)
        arr[self.randMaskArr != 0] = np.abs(self.randArr[self.randMaskArr != 0])
        self.TestIm(im, arr, 
                    name='Abs_Masked', disp=disp)

    #
    # Add Masked
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(randMaskSetUp, randMaskTearDown)
    def test_AddMasked(self, disp=False):
        # Add_I
        im = self.randIm.copy()
        arr2, im2 = self.genRandImPair()
        Add_I(im, im2, self.randMask)
        arr = self.randArr.copy()
        arr[self.randMaskArr != 0] += arr2[self.randMaskArr != 0]
        self.TestIm(im, arr, 
                    name='Add_I_Masked', disp=disp)

        # AddC_I
        c = np.random.rand()
        im = self.randIm.copy()
        Add_I(im, c, self.randMask)
        arr = self.randArr.copy()
        arr[self.randMaskArr != 0] += c
        self.TestIm(im, arr, 
                    name='AddC_I_Masked', disp=disp)

        # Add
        arr, im = self.genRandImPair()
        Add(im, self.randIm, im2, self.randMask)
        arr[self.randMaskArr != 0] = \
            self.randArr[self.randMaskArr != 0] + \
            arr2[self.randMaskArr != 0]
        self.TestIm(im, arr, 
                    name='Add_Masked', disp=disp)

        # AddC
        c = np.random.rand()
        arr, im = self.genRandImPair()
        Add(im, self.randIm, c, self.randMask)
        arr[self.randMaskArr != 0] = \
            self.randArr[self.randMaskArr != 0] + c
        self.TestIm(im, arr, 
                    name='AddC_Masked', disp=disp)

    #
    # SoftSgn
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    def test_SoftSgn(self, disp=False):
        eps = 0.5
        arr = self.randArr.copy()
        arr[arr<-eps] = -1.0
        arr[arr>eps] = 1.0
        arr[np.abs(arr)<=eps] = arr[np.abs(arr)<=eps]/eps
        SoftSgn_I(self.randIm,eps)

        self.TestIm(self.randIm, arr, 'SoftSgn', disp=disp)

        if disp:
            x = np.arange(-eps*5, eps*5, eps/10.0)
            xim = common.ImFromNPArr(x, self.mType)
            yim = xim.copy()
            SoftSgn_I(yim,eps)
            plt.figure('SoftSgnFunc')
            xim.toType(MEM_HOST)
            yim.toType(MEM_HOST)
            plt.plot(np.squeeze(xim.asnp()), np.squeeze(yim.asnp()),'b')
            plt.plot(np.squeeze(xim.asnp()), np.sign(np.squeeze(xim.asnp())),'r')

    #
    # Convolve
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    def test_Convolve(self, disp=False):
        im = self.randIm.copy()
        # kernel will just shift everything to the left one pixel
        shift_x_arr = np.array([[[0]],[[0]],[[1]]])
        shift_x = common.ImFromNPArr(shift_x_arr, self.mType)
        self.randArr[:-1,:] = self.randArr[1:,:]
        Convolve(im, self.randIm, shift_x)

        self.TestIm(im, self.randArr, \
                        'convolve', disp=disp)

    #
    # UpwindDiff
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    def test_UpwindDiff(self, disp=False):
        im = self.randIm.copy()
        randSpeedArr = np.random.randn(*self.sz)
        randSpeedIm = common.ImFromNPArr(randSpeedArr, \
                                             mType=self.mType, \
                                             sp=self.imSp)
        forDiff = canp.FiniteDiff(self.randArr, DIM_X, \
                                  DIFF_FORWARD, BC_CLAMP, self.imSp)
        backDiff = canp.FiniteDiff(self.randArr, DIM_X, \
                                   DIFF_BACKWARD, BC_CLAMP, self.imSp)
        self.randArr[randSpeedArr > 0] = backDiff[randSpeedArr > 0]
        self.randArr[randSpeedArr < 0] = forDiff[randSpeedArr < 0]
        UpwindDiff(im, self.randIm, randSpeedIm, DIM_X)

        self.TestIm(im, self.randArr, \
                        'UpwindDiff X', disp=disp)

    #
    # NormalizeSafe
    #
    @PyCATest.AddSetUp(randFieldSetUp, randFieldTearDown)
    def test_NormalizeSafe(self, disp=False):
        eps = 1.0 # very large eps, just for testing
        l = np.sqrt(self.randFArr[:,:,0]**2 + self.randFArr[:,:,1]**2)
        l = np.tile(np.atleast_3d(l),(1,1,2))
        self.randFArr[l>eps] = self.randFArr[l>eps]/l[l>eps]
        self.randFArr[l<=eps] = 0
        NormalizeSafe_I(self.randF, eps)

        self.TestField(self.randF, self.randFArr, \
                           'NormalizeSafe_I', disp=disp)


    #
    # FiniteDiff
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    def test_FiniteDiff(self, disp=False):
        im2 = self.randIm.copy()
        for dim in [DIM_X, DIM_Y]:
            for diffType in [DIFF_FORWARD, DIFF_BACKWARD, DIFF_CENTRAL]:
                for bc in [BC_CLAMP, BC_WRAP, BC_APPROX]:
                    diff = canp.FiniteDiff(self.randArr, dim, diffType, bc, self.imSp)
                    FiniteDiff(im2, self.randIm, dim, diffType, bc)
                    pltname = '%s %s %s'%\
                        (PyCATest.DIMNAMES[dim], \
                             PyCATest.DIFFTNAMES[diffType], \
                             PyCATest.BCNAMES[bc])
                    self.TestIm(im2, diff, name=pltname, disp=disp)

    #
    # Gradient
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(resultFieldSetUp, resultFieldTearDown)
    def test_Gradient(self, disp=False):
        for diffType in [DIFF_FORWARD, DIFF_BACKWARD, DIFF_CENTRAL]:
            for bc in [BC_CLAMP, BC_WRAP, BC_APPROX]:
                gArr = canp.Grad(self.randArr, diffType, bc, self.imSp)
                Gradient(self.field, self.randIm, diffType, bc)
                pltname = 'Gradient %s %s'%\
                    (PyCATest.DIFFTNAMES[diffType], PyCATest.BCNAMES[bc])
                self.TestField(self.field, gArr, name=pltname, disp=disp)

    #
    # Jac. Det.
    #
    @PyCATest.AddSetUp(randFieldSetUp, randFieldTearDown)
    @PyCATest.AddSetUp(resultImSetUp, resultImTearDown)
    def test_JacDet(self, disp=False):
        arr = canp.JacDet(self.randFArr, DIFF_CENTRAL, BC_APPROX, self.imSp)
        JacDetH(self.im, self.randF, DIFF_CENTRAL, BC_APPROX)
        self.TestIm(self.im, arr, name='JacDet', disp=disp)

    #
    # Masked Gradient
    #
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(maskSetUp, maskTearDown)
    @PyCATest.AddSetUp(resultFieldSetUp, resultFieldTearDown)
    def test_GradMask(self, disp=False):
        diffType = DIFF_CENTRAL
        bc = BC_CLAMP
        gArr = canp.GradMask(self.randArr, self.maskArr, diffType, bc, self.imSp)
        GradientMask(self.field, self.randIm, self.mask, diffType, bc)
        pltname = 'GradientMask %s %s'%\
            (PyCATest.DIFFTNAMES[diffType], PyCATest.BCNAMES[bc])
        self.TestField(self.field, gArr, name=pltname, disp=disp)


    #
    # Gaussian Blur
    #
    #unittest.skip("We know this will fail")
    @PyCATest.SkipIfNotDisp("We know this test fails")
    @PyCATest.AddSetUp(randImSetUp, randImTearDown)
    @PyCATest.AddSetUp(resultImSetUp, resultImTearDown)
    def test_GaussianBlur(self, disp=False):
        import scipy.ndimage as ndimage
        tmp = self.randIm.copy()
        sig = 2.0
        # mode controls the border handling
        #arr = ndimage.filters.gaussian_filter(arr,sig,mode='reflect')
        arr = ndimage.filters.gaussian_filter(self.randArr,sig,mode='nearest')
        #arr = ndimage.filters.gaussian_filter(arr,sig,mode='constant', cval=0)
        if self.mType == MEM_HOST:
            gaussFilter = GaussianFilterCPU()
        else:
            gaussFilter = GaussianFilterGPU()
        sigVec = Vec3Df(sig, sig, sig)
        kVec = Vec3Di(int(6*sig),int(6*sig),int(6*sig))
        gaussFilter.updateParams(self.imSz, sigVec, kVec)
        gaussFilter.Filter(self.im, self.randIm, tmp)

        self.TestIm(self.im, arr, 'Gaussian Filter', disp=disp)


    #
    # Simple math operations
    #

    def RunOp(self, disp=False,
              arg1Type='Image3D', arg2Type='Image3D', outType='Image3D',
              pycaOpName='Add', npOpName='__add__',
              inPlace=False):
        """Test simple math operations (+,-,*=,etc.) on Imag3D and Field3D
        """

        if disp:
            print 'arg1Type=', arg1Type, ', arg2Type=', arg2Type,\
                ', outType=', outType,\
                ', pycaOpName=', pycaOpName, ', npOpName=', npOpName,\
                ', inPlace=', inPlace

        #
        # set up data
        #

        # if we're testing division, avoid divide-by-zero (or
        # near-zero)
        NonNeg = False
        if npOpName == '__div__' or npOpName == '__idiv__':
            NonNeg = True

        if disp:
            print 'NonNeg:', NonNeg

        # Create first argument
        if arg1Type == 'Image3D':
            npArg1, pycaArg1 = self.genRandImPair(NonNeg=NonNeg)
        elif arg1Type == 'Field3D':
            npArg1, pycaArg1 = self.genRandFieldPair(NonNeg=NonNeg)
        elif arg1Type == 'Vec3Df':
            npArg1, pycaArg1 = self.genRandVecfPair(NonNeg=NonNeg)
        elif arg1Type == 'Vec3Di':
            npArg1, pycaArg1 = self.genRandVeciPair(NonNeg=NonNeg)
        else:
            raise Exception('Unknown first arg type '+arg1Type)

        # Create second argument
        if arg2Type == 'Image3D':
            npArg2, pycaArg2 = self.genRandImPair(NonNeg=NonNeg)
        elif arg2Type == 'Field3D':
            npArg2, pycaArg2 = self.genRandFieldPair(NonNeg=NonNeg)
        elif arg2Type == 'Vec3Df':
            npArg2, pycaArg2 = self.genRandVecfPair(NonNeg=NonNeg)
            # handle Field3D (op) Vec3Df case
            if arg1Type == 'Field3D' and outType == 'Field3D':
                pycaArg2.z = 0
                npArg2 = npArg2[:2]
        elif arg2Type == 'Vec3Di':
            npArg2, pycaArg2 = self.genRandVeciPair(NonNeg=NonNeg)
        elif arg2Type == 'float':
            npArg2 = float(np.random.rand())
            if NonNeg:
                npArg2 = abs(npArg2) + self.minNonNeg
            pycaArg2 = npArg2
        elif arg2Type == 'int':
            npArg2 = int(np.random.random_integers(self.minRandInt,
                                                   self.maxRandInt,
                                                   1))
            if NonNeg:
                npArg2 = int(np.maximum(npArg2,1))
            pycaArg2 = npArg2
        else:
            raise Exception('Unknown second arg type '+arg2Type)

        #
        # execute test functions
        #
        if pycaOpName.startswith('__'):
            pycaOp = pycaArg1.__getattribute__(pycaOpName)
            if inPlace:
                # e.g. pycaArg1.__iadd__(pycaArg2)
                #      (equivalent: pycaArg1 += pycaArg2)
                pycaOp(pycaArg2)
            else:
                # e.g. pycaRtnVal = pycaArg1.__add__(pycaArg2)
                #      (equivalent: pycaRtnVal = pycaArg1 + pycaArg2)
                pycaRtnVal = pycaOp(pycaArg2)
                if pycaRtnVal is NotImplemented:
                    common.DebugHere()
        else:
            pycaOp = eval(pycaOpName)
            if inPlace:
                # e.g. Add_I(pycaArg1, pycaArg2)
                pycaOp(pycaArg1, pycaArg2)
            else:
                # e.g. Add(pycaRtnVal, pycaArg1, pycaArg2)
                if outType == 'Image3D':
                    pycaRtnVal = Image3D(self.grid, self.mType)
                elif outType == 'Field3D':
                    pycaRtnVal = Field3D(self.grid, self.mType)
                else:
                    raise Exception('invalid first arg type for '
                                    'argument return value'
                                    + arg1Type)

                pycaOp(pycaRtnVal, pycaArg1, pycaArg2)

        if npOpName.startswith('__'):
            npOp = npArg1.__getattribute__(npOpName)
            if inPlace:
                npOp(npArg2)
            else:
                npRtnVal = npOp(npArg2)
        else:
            npOp = eval(npOpName)
            if inPlace:
                npOp(npArg1, npArg2)
            else:
                npRtnVal = npOp(npArg1, npArg2)

        # test results
        if inPlace:
            pycaRtnVal = pycaArg1
            npRtnVal = npArg1

        if outType == 'Image3D':
            self.TestIm(pycaRtnVal, npRtnVal, \
                        name=pycaOpName, disp=disp)
        elif outType == 'Field3D':
            self.TestField(pycaRtnVal, npRtnVal, \
                           name=pycaOpName, disp=disp)
        elif outType == 'Vec3Df':
            diff = np.sum(np.abs(np.array(pycaRtnVal.tolist())-npRtnVal))
            self.assertLess(diff, self.MaxEps)
        elif outType == 'Vec3Di':
            diff = np.sum(np.abs(np.array(pycaRtnVal.tolist())-npRtnVal))
            self.assertTrue(diff == 0)
        else:
            raise Exception('Unknown first arg type '+firstArgType)


    # Image3D math tests
    def test_BinaryOps(self, disp=False):
        for optest in self.BinaryOperTests:
            self.RunOp(disp=disp,
                       arg1Type=optest[0], arg2Type=optest[1],
                       outType=optest[2],
                       pycaOpName=optest[3],
                       npOpName=optest[4],
                       inPlace=optest[5])




    # runTest is only added so that the class can be instantiated
    # directly in order to call individual tests
    def runTest():
        print 'No tests to run directly, all are member functions'

# end NumpyTestCase

#
# \min_{u,b} \|I_s-I(x+u(x))\|^2  + \alpha<Lu,u>
#
#

from PyCA.Core import *

import PyCA.Common as common
import PyCA.Display as display

import numpy as np
import matplotlib.pyplot as plt

import os, errno

def GreedyReg\
        (I0Orig, \
         I1Orig, \
         scales = [1], \
         nIters = [1000], \
         ustep = [0.25], \
         fluidParams = [0.1, 0.1, 0.001], \
         plotEvery = 100):

    mType = I0Orig.memType()
    origGrid = I0Orig.grid()

    # allocate vars
    I0 = Image3D(origGrid, mType)
    I1 = Image3D(origGrid, mType)
    h = Field3D(origGrid, mType)
    Idef = Image3D(origGrid, mType)
    diff = Image3D(origGrid, mType)
    gI = Field3D(origGrid, mType)
    gU = Field3D(origGrid, mType)
    scratchI = Image3D(origGrid, mType)
    scratchV = Field3D(origGrid, mType)

    # allocate diffOp
    if mType == MEM_HOST:
        diffOp = FluidKernelFFTCPU()
    else:
        diffOp = FluidKernelFFTGPU()

    # initialize some vars
    zerovec = Vec3Df(0.0, 0.0, 0.0)

    nScales = len(scales)
    scaleManager = MultiscaleManager(origGrid)
    for s in scales:
        scaleManager.addScaleLevel(s)

    # Initalize the thread memory manager (needed for resampler)
    # num pools is 2 (images) + 2*3 (fields)
    ThreadMemoryManager.init(origGrid, mType, 8)

    if mType == MEM_HOST:
        resampler = MultiscaleResamplerGaussCPU(origGrid)
    else:
        resampler = MultiscaleResamplerGaussGPU(origGrid)

    def setScale(scale):
        global curGrid

        scaleManager.set(scale)
        curGrid = scaleManager.getCurGrid()
        # since this is only 2D:
        curGrid.spacing().z = 1.0;

        resampler.setScaleLevel(scaleManager)

        diffOp.setAlpha(fluidParams[0])
        diffOp.setBeta(fluidParams[1])
        diffOp.setGamma(fluidParams[2])
        diffOp.setGrid(curGrid)

        # downsample images
        I0.setGrid(curGrid)
        I1.setGrid(curGrid)
        if scaleManager.isLastScale():
            Copy(I0, I0Orig)
            Copy(I1, I1Orig)
        else:
            resampler.downsampleImage(I0,I0Orig)
            resampler.downsampleImage(I1,I1Orig)

        # initialize / upsample deformation
        if scaleManager.isFirstScale():
            h.setGrid(curGrid)
            SetToIdentity(h)
        else:
            resampler.updateHField(h)

        # set grids
        gI.setGrid(curGrid)
        Idef.setGrid(curGrid)
        diff.setGrid(curGrid)
        gU.setGrid(curGrid)
        scratchI.setGrid(curGrid)
        scratchV.setGrid(curGrid)
    # end function

    energy = [[] for _ in xrange(3)]

    for scale in range(len(scales)):

        setScale(scale)

        for it in range(nIters[scale]):
            print 'iter %d'%it

            # compute deformed image
            ApplyH(Idef, I0, h)

            # update gradient
            Gradient(gI, Idef)

            # update u
            Sub(diff, I1, Idef)

            gI *= diff

            diffOp.applyInverseOperator(gU, gI)

            gU *= ustep[scale]
            # ApplyV(scratchV, h, gU, BACKGROUND_STRATEGY_PARTIAL_ID)
            ComposeHV(scratchV, h, gU)
            h.swap(scratchV)

            # compute energy
            energy[0].append(Sum2(diff))

            if it % plotEvery == 0 or it == nIters[scale]-1:
                clrlist = ['r','g','b','m','c','y','k']
                plt.figure('energy')
                for i in range(len(energy)):
                    plt.plot(energy[i],clrlist[i])
                    if i == 0:
                        plt.hold(True)
                plt.hold(False)
                plt.draw()

                plt.figure('results')
                plt.clf()
                plt.subplot(3,2,1)
                display.DispImage(I0, 'I0', newFig=False)
                plt.subplot(3,2,2)
                display.DispImage(I1, 'I1', newFig=False)
                plt.subplot(3,2,3)
                display.DispImage(Idef, 'def', newFig=False)
                plt.subplot(3,2,4)
                display.DispImage(diff, 'diff', newFig=False)
                plt.colorbar()
                plt.subplot(3,2,5)
                display.GridPlot(h, every=4, isVF=False)
                plt.draw()
                plt.show()

            # end plot
        # end iteration
    # end scale
    return (Idef, h, energy)
# end function

if __name__ == '__main__':

    plt.close('all')

    if GetNumberOfCUDADevices() > 0:
        mType = MEM_DEVICE
    else:
        print "No CUDA devices found, running on CPU"
        mType = MEM_HOST

    imagedir='./Images/'

    #
    # Run lena images
    #

    I0 = common.LoadPNGImage(imagedir + 'lena_deformed.png', mType)
    I1 = common.LoadPNGImage(imagedir + 'lena_orig.png', mType)

    (Idef, h, energy) = \
           GreedyReg(I0, \
                    I1, \
                    scales = [2,1], \
                    nIters = [1000,1000], \
                    ustep = [0.1, 0.1], \
                    fluidParams = [0.5, 0.5, 0.001], \
                    plotEvery = 500)


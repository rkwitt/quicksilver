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

def ElastReg\
        (I0, \
             I1, \
             nIters = 1000, \
             ustep = 0.25, \
             fluidParams = [0.1, 0.1, 0.001], \
             VFC = 0.2, \
             Mask = None, \
             plotEvery = 100):

    mType = I0.memType()
    grid = I0.grid()
 
    # allocate vars
    u = Field3D(grid, mType)
    SetMem(u,0.0)
    Idef = Image3D(grid, mType)
    diff = Image3D(grid, mType)
    gI = Field3D(grid, mType)
    gU = Field3D(grid, mType)
    scratchI = Image3D(grid, mType)
    scratchV = Field3D(grid, mType)

    # allocate diffOp 
    if mType == MEM_HOST:
        diffOp = FluidKernelFFTCPU()
    else:
        diffOp = FluidKernelFFTGPU()
    diffOp.setAlpha(fluidParams[0])
    diffOp.setBeta(fluidParams[1])
    diffOp.setGamma(fluidParams[2])
    diffOp.setGrid(grid)

    energy = [[] for _ in xrange(3)]

    # update gradient
    Gradient(gI, I0)

    for it in range(nIters):
        print 'iter %d'%it

        # compute deformed image
        ApplyV(Idef, I0, u, 1.0)

        # update u
        Sub(diff, I1, Idef)

        if Mask != None:
            Mul_I(diff, Mask)

        ApplyV(scratchV, gI, u, BACKGROUND_STRATEGY_CLAMP)
        Mul_I(scratchV, diff)

        diffOp.applyInverseOperator(gU, scratchV)

        # for computing energy
        diffOp.applyOperator(scratchV, u)
        
        # u =  u*(1-VFC*ustep) + (-2.0*ustep)*gU
        MulC_Add_MulC_I(u, (1-VFC*ustep),
                               gU, 2.0*ustep)

        # compute energy
        energy[0].append(Sum2(diff))
        energy[1].append(0.5*VFC*Dot(scratchV, u))
        energy[2].append(energy[0][-1]+\
                         energy[1][-1])

        if plotEvery > 0 and \
               ((it+1) % plotEvery == 0 or \
                it == nIters-1):
            print 'plotting'
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
            display.GridPlot(u, every=4)
            plt.draw()
            plt.show()

        # end plot
    # end iteration
    return (Idef, u, energy)
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

    (Idef, u, energy) = \
           ElastReg(I0, \
                    I1, \
                    nIters = 1000, \
                    ustep = 0.2, \
                    fluidParams = [0.5, 0.5, 0.001], \
                    VFC = 0.2, \
                    plotEvery = 500)


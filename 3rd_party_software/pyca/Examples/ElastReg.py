#
# \argmin{u} \|I_s-I(x+u(x))\|^2  + 0.5*\alpha<Lu,u>
#
#

import PyCA.Core as ca

import PyCA.Common as common
import PyCA.Display as display

import numpy as np
import matplotlib.pyplot as plt

def ElastReg(I0Orig, 
             I1Orig, 
             scales = [1], 
             nIters = [1000], 
             maxPert = [0.2], 
             fluidParams = [0.1, 0.1, 0.001], 
             VFC = 0.2, 
             Mask = None, 
             plotEvery = 100):

    mType = I0Orig.memType()
    origGrid = I0Orig.grid()

    # allocate vars
    I0 = ca.Image3D(origGrid, mType)
    I1 = ca.Image3D(origGrid, mType)
    u = ca.Field3D(origGrid, mType)
    Idef = ca.Image3D(origGrid, mType)
    diff = ca.Image3D(origGrid, mType)
    gI = ca.Field3D(origGrid, mType)
    gU = ca.Field3D(origGrid, mType)
    scratchI = ca.Image3D(origGrid, mType)
    scratchV = ca.Field3D(origGrid, mType)

    # mask
    if Mask != None:
        MaskOrig = Mask.copy()

    # allocate diffOp
    if mType == ca.MEM_HOST:
        diffOp = ca.FluidKernelFFTCPU()
    else:
        diffOp = ca.FluidKernelFFTGPU()

    # initialize some vars
    nScales = len(scales)
    scaleManager = ca.MultiscaleManager(origGrid)
    for s in scales:
        scaleManager.addScaleLevel(s)

    # Initalize the thread memory manager (needed for resampler)
    # num pools is 2 (images) + 2*3 (fields)
    ca.ThreadMemoryManager.init(origGrid, mType, 8)

    if mType == ca.MEM_HOST:
        resampler = ca.MultiscaleResamplerGaussCPU(origGrid)
    else:
        resampler = ca.MultiscaleResamplerGaussGPU(origGrid)

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
            ca.Copy(I0, I0Orig)
            ca.Copy(I1, I1Orig)
        else:
            resampler.downsampleImage(I0,I0Orig)
            resampler.downsampleImage(I1,I1Orig)

        if Mask != None:
            if scaleManager.isLastScale():
                Mask.setGrid(curGrid)
                ca.Copy(Mask, MaskOrig)
            else:
                resampler.downsampleImage(Mask,MaskOrig)

        # initialize / upsample deformation
        if scaleManager.isFirstScale():
            u.setGrid(curGrid)
            ca.SetMem(u, 0.0)
        else:
            resampler.updateVField(u)

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
        ustep = None
        # update gradient
        ca.Gradient(gI, I0)

        for it in range(nIters[scale]):
            print 'iter %d'%it

            # compute deformed image
            ca.ApplyV(Idef, I0, u, 1.0)

            # update u
            ca.Sub(diff, I1, Idef)

            if Mask != None:
                ca.Mul_I(diff, Mask)

            ca.ApplyV(scratchV, gI, u, ca.BACKGROUND_STRATEGY_CLAMP)
            ca.Mul_I(scratchV, diff)

            diffOp.applyInverseOperator(gU, scratchV)

            vfcEn = VFC*ca.Dot(scratchV, gU)

            # why is this negative necessary?
            ca.MulC_I(gU, -1.0)

            # u =  u*(1-VFC*ustep) + (-2.0*ustep)*gU
            # MulC_Add_MulC_I(u, (1-VFC*ustep),
            #                        gU, 2.0*ustep)

            # u =  u - ustep*(VFC*u + 2.0*gU)
            ca.MulC_I(gU, 2.0)

            # subtract average if gamma is zero (result of nullspace
            # of L for K(L(u)))
            if fluidParams[2] == 0:
                av = ca.SumComp(u)
                av /= scratchI.nVox()
                ca.SubC(scratchV, u, av)
            # continue computing gradient
            ca.MulC(scratchV, u, VFC)
            ca.Add_I(gU, scratchV)

            ca.Magnitude(scratchI, gU)
            gradmax = ca.Max(scratchI)
            if ustep is None or ustep*gradmax > maxPert:
                ustep = maxPert[scale]/gradmax
                print 'step is %f'%ustep
            
            ca.MulC_I(gU, ustep)
            # apply gradient
            ca.Sub_I(u, gU)

            # compute energy
            energy[0].append(ca.Sum2(diff))
            diffOp.applyOperator(scratchV, u)
            energy[1].append(0.5*VFC*ca.Dot(u, scratchV))
            energy[2].append(energy[0][-1]+\
                             energy[1][-1])

            if plotEvery > 0 and \
                   ((it+1) % plotEvery == 0 or \
                    (scale == nScales-1 and it == nIters[scale]-1)):
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
                plt.subplot(3,2,6)
                display.JacDetPlot(u)
                plt.colorbar()
                plt.draw()
                plt.show()

            # end plot
        # end iteration
    # end scale
    return (Idef, u, energy)
# end function

if __name__ == '__main__':

    plt.close('all')

    if ca.GetNumberOfCUDADevices() > 0:
        mType = ca.MEM_DEVICE
    else:
        print "No CUDA devices found, running on CPU"
        mType = ca.MEM_HOST


    imagedir='./Images/'

    #
    # Run lena images
    #

    I0 = common.LoadPNGImage(imagedir + 'lena_deformed.png', mType)
    I1 = common.LoadPNGImage(imagedir + 'lena_orig.png', mType)

    (Idef, u, energy) = \
           ElastReg(I0, 
                    I1, 
                    scales = [2, 1], 
                    nIters = [1000, 1000], 
                    maxPert = [0.1, 0.1], 
                    fluidParams = [0.5, 0.5, 0.001], 
                    # fluidParams = [0.5, 0.0, 0.001], 
                    # fluidParams = [0.5, 0.5, 0.0], 
                    VFC = 0.2, 
                    plotEvery = 500)


from PyCA.Core import *

import PyCA.Common as common
import PyCA.Display as display

import numpy as np
import matplotlib.pyplot as plt

def PrimalDualTV(I0, \
                 DataFidC, \
                 TVC = 1.0, \
                 nIters = 5000, \
                 stepP = None, \
                 stepI = None, \
                 disp = False, \
                 dispEvery = 0):


    #
    # Initialize data
    #
    mType = I0.memType()
    grid = I0.grid().copy()

    if stepP == None:
        stepP = 1.0/8.0

    if stepI == None:
        stepI = min(stepP,1.0/DataFidC)

    bc = BC_CLAMP
    # bc = BC_WRAP

    # primal var
    I = I0.copy()
    # I = Image3D(grid, mType)
    # SetMem(I, 0.0)


    # dual var
    p = Field3D(grid, mType)
    
    # zerovec = Vec3Df(0.0,0.0,0.0)
    # SetMem(p, zerovec)
    Gradient(p, I0, DIFF_FORWARD, bc)
    ReprojectToUnitVec(p)

    # Initialize other data
    energy = [[] for _ in xrange(2)]

    #
    # Allocate all necessary data
    #
    scratchI = Image3D(grid, mType)
    scratchI2 = Image3D(grid, mType)
    scratchV = Field3D(grid, mType)

    EnergyFig = plt.figure('PrimalDual Energy');
    plt.clf();
    ResultFig = plt.figure('PrimalDual Results');
    plt.clf();

    # overwrites LDefSum
    def plotResults(fig,cmap='gray',rng=[0,1]):
        plt.figure(fig)
        plt.subplot(1,3,1)
        display.DispImage(I0, 'Orig', cmap=cmap, \
                                   newFig=False, rng=rng, t=False)
        plt.subplot(1,3,2)
        display.DispImage(I, 'Denoised', cmap=cmap, \
                                   newFig=False, rng=rng, t=False)
        Sub(scratchI, I, I0)
        plt.subplot(1,3,3)
        display.DispImage(scratchI, 'diff', cmap=cmap, \
                                   newFig=False, rng=None, t=False)
        plt.draw()
        plt.show()

    def plotEnergy(en, fig):
        plt.figure(fig)
        plt.plot(en[0][1:],'r')
        plt.hold(True)
        plt.plot(en[1][1:],'g')
        plt.hold(False)
        plt.draw()
        plt.show()

    for k in range(nIters+1):

        print 'iteration %d...'%k

        #
        # Display images
        #
        if disp and dispEvery > 0 and k%dispEvery == 0:
            plotResults(ResultFig.number)

        #
        # Compute energy
        #

        # primal energy

        Sub(scratchI, I, I0)
        primalEnergy = (DataFidC/2.0)*Sum2(scratchI)
        GradientMag(scratchI, I, DIFF_FORWARD, bc)
        primalEnergy += TVC*Sum(scratchI)

        # dual energy

        Divergence(scratchI, p, DIFF_BACKWARD, bc)
        MulC_I(scratchI, TVC/DataFidC)
        Sqr_I(scratchI)
        Divergence(scratchI2, p, DIFF_BACKWARD, bc)
        MulC_I(scratchI2, 2.0*(TVC/DataFidC))
        Mul_I(scratchI2, I0)
        Add_I(scratchI, scratchI2)
        dualEnergy = (-DataFidC/2.0)*Sum(scratchI)

        energy[0].append(primalEnergy)
        energy[1].append(dualEnergy)

        if disp and dispEvery > 0 and k%dispEvery == 0:
            plotEnergy(energy, EnergyFig.number)

        # just compute energy on final iteration
        if k >= nIters:
            break

        # primal step

        # scratchI = I - I0 - (TVC/DataFidC)*div(p)
        Divergence(scratchI, p, DIFF_BACKWARD, bc)
        MulC_I(scratchI, -TVC/DataFidC)
        Sub(scratchI2, I, I0)
        Add_I(scratchI, scratchI2)
        # I = I - stepI*gI
        Add_MulC_I(I, scratchI, -stepI)

        # dual step
        Gradient(scratchV, I, DIFF_FORWARD, bc)
        # weighting update by 1/TVC to speed convergence
        #Add_MulC_I(p, scratchV, stepP*TVC)
        Add_MulC_I(p, scratchV, stepP)
        # reproject onto constraint
        ReprojectToUnitVec(p)

    if disp:
        plotResults(ResultFig.number)
        plotEnergy(energy, EnergyFig.number)
    
    return (I, energy)

#
# End function
#

if __name__ == '__main__':
    plt.close('all')

    # number of iterations
    nIters = 2000

    disp = True
    dispEvery = 1000

    if GetNumberOfCUDADevices() > 0:
        mType = MEM_DEVICE
    else:
        print "No CUDA devices found, running on CPU"
        mType = MEM_HOST

    # data fidelity modifier
    DataFidC = 1.0
    TVC = 0.05

    imagedir='./Images/'

    #
    # Run lena images
    #
    
    I0 = common.LoadPNGImage(imagedir + 'lena_orig.png', mType)

    imSz = I0.size()
    sz = imSz.tolist()[0:2]

    (I,energy) = \
        PrimalDualTV(I0, \
                     DataFidC, \
                     TVC = TVC, \
                     nIters = nIters, \
                     stepP = 1.0, \
                     stepI = 1.0/16.0, \
                     disp = disp, \
                     dispEvery = dispEvery)

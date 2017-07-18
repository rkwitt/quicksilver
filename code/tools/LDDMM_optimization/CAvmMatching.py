#!/usr/bin/python2

# configuration files related modules
from Configs import Config, Optim, Compute, VMConfig
import os.path

# pyca modules
import PyCA.Core as ca
import PyCA.Common as common
import PyCA.Display as display

# vector momentum modules
from Libraries import CAvmCommon

# others
import numpy as np
import matplotlib.pyplot as plt
import os, errno

import logging
import sys
import copy
import math
import time

StudySpec = {
    'I0':
    Config.Param(default='subject1.mhd', required=True,
                    comment="Initial (moving) image file"),
    'I1':
    Config.Param(default='subject2.mhd', required=True,
                    comment="Target (fixed) image file")}

MatchingConfigSpec = {
    'compute': Compute.ComputeConfigSpec,
    'vectormomentum': VMConfig.VMConfigSpec,
    'study': StudySpec,
    'optim': Optim.OptimConfigSpec,
    'io': {
        'plotEvery':
        Config.Param(default=10,
                     comment="Update plots every N iterations"),
        'plotSlice':
        Config.Param(default=None,
                     comment="Slice to plot.  Defaults to mid axial"),
        'quiverEvery':
        Config.Param(default=1,
                     comment="How much to downsample for quiver plots"),
        'outputPrefix':
        Config.Param(default="./",
                     comment="Where to put output.  Don't forget trailing "
                     + "slash")},
    '_resource': 'VectorMomentum_Matching'}


class MatchingVariables:
    """
    Parameters for doing matching or regression.

    The constructor actually does the allocation of all scratch variables
    required for computation of PDiff matching gradients
    """
    def __init__(self,I0,I1, sigma, t, cpinds, cpstates, alpha, beta, gamma, nIter, stepSize, maxPert, nTimeSteps, integMethod='RK4', optMethod='FIXEDGD', nInv=10, plotEvery=None, plotSlice=None, quiverEvery=None, outputPrefix='./'):
        """
        Initialize everything with the size and type given
        """
        print I0
        self.I0 = I0
        self.I1 = I1
        self.grid = I0.grid()
        print(self.grid)
        self.memtype = I0.memType()

        # matching param
        self.sigma = sigma

        # initial conditions
        self.m0 = ca.Field3D(self.grid, self.memtype)
        ca.SetMem(self.m0,0.0)

        # state variables
        self.g = ca.Field3D(self.grid, self.memtype)
        self.ginv = ca.Field3D(self.grid, self.memtype)
        self.m = ca.Field3D(self.grid, self.memtype)
        self.I = ca.Image3D(self.grid, self.memtype)
        self.residualIm = ca.Image3D(self.grid, self.memtype)

        # adjoint variables
        self.madj = ca.Field3D(self.grid, self.memtype)
        self.Iadj = ca.Image3D(self.grid, self.memtype)
        self.madjtmp = ca.Field3D(self.grid, self.memtype)
        self.Iadjtmp = ca.Image3D(self.grid, self.memtype)

        # time array
        self.t = t

        # checkpointing variables
        self.checkpointinds = cpinds
        self.checkpointstates = cpstates

        # set up diffOp
        if self.memtype == ca.MEM_HOST:
            self.diffOp = ca.FluidKernelFFTCPU()
        else:
            self.diffOp = ca.FluidKernelFFTGPU()
        self.diffOp.setAlpha(alpha)
        self.diffOp.setBeta(beta)
        self.diffOp.setGamma(gamma)
        self.diffOp.setGrid(self.grid)

        # energy
        self.Energy = None

        # optimization stuff
        self.nIter = nIter
        self.stepSize = stepSize
        self.maxPert = maxPert
        self.nTimeSteps = nTimeSteps
        self.optMethod = optMethod
        self.integMethod = integMethod
        self.nInv = nInv # for interative update to inverse deformation

        # plotting variables
        self.plotEvery = plotEvery
        self.plotSlice = plotSlice
        self.quiverEvery = quiverEvery
        self.outputPrefix = outputPrefix
        # scratch variables
        self.scratchV1 = ca.Field3D(self.grid, self.memtype)
        self.scratchV2 = ca.Field3D(self.grid, self.memtype)
        self.scratchV3 = ca.Field3D(self.grid, self.memtype)

        self.scratchV4 = ca.Field3D(self.grid, self.memtype)
        self.scratchV5 = ca.Field3D(self.grid, self.memtype)
        self.scratchV6 = ca.Field3D(self.grid, self.memtype)
        self.scratchV7 = ca.Field3D(self.grid, self.memtype)
        self.scratchV8 = ca.Field3D(self.grid, self.memtype)
        self.scratchV9 = ca.Field3D(self.grid, self.memtype)
        self.scratchV10 = ca.Field3D(self.grid, self.memtype)
        self.scratchV11 = ca.Field3D(self.grid, self.memtype)


    # end __init__
# end MatchingVariables

def MatchingGradient(p):
    # shoot the geodesic forward
    CAvmCommon.IntegrateGeodesic(p.m0,p.t,p.diffOp, \
                      p.m, p.g, p.ginv,\
                      p.scratchV1, p.scratchV2,p. scratchV3,\
                      p.checkpointstates, p.checkpointinds,\
                      Ninv=p.nInv, integMethod = p.integMethod, RK4=p.scratchV4,scratchG=p.scratchV5)


    endidx = p.checkpointinds.index(len(p.t)-1)
    # compute residual image
    ca.ApplyH(p.residualIm,p.I0,p.ginv)
    ca.Sub_I(p.residualIm, p.I1)
    # while we have residual, save the image energy
    IEnergy = ca.Sum2(p.residualIm)/(2*p.sigma*p.sigma*float(p.I0.nVox()))

    ca.DivC_I(p.residualIm, p.sigma*p.sigma) # gradient at measurement

    # integrate backward
    CAvmCommon.IntegrateAdjoints(p.Iadj,p.madj,\
                      p.I,p.m,p.Iadjtmp, p.madjtmp,p.scratchV1,\
                      p.scratchV2,p.scratchV3,\
                      p.I0,p.m0,\
                      p.t, p.checkpointstates, p.checkpointinds,\
                      [p.residualIm], [endidx],\
                      p.diffOp,
                      p.integMethod, p.nInv, \
                      scratchV3=p.scratchV7, scratchV4=p.g,scratchV5=p.ginv,scratchV6=p.scratchV8, scratchV7=p.scratchV9, \
                      scratchV8=p.scratchV10,scratchV9=p.scratchV11,\
                      RK4=p.scratchV4, scratchG=p.scratchV5, scratchGinv=p.scratchV6)


    # compute gradient
    ca.Copy(p.scratchV1, p.m0)
    p.diffOp.applyInverseOperator(p.scratchV1)
    # while we have velocity, save the vector energy
    VEnergy = 0.5*ca.Dot(p.m0,p.scratchV1)/float(p.I0.nVox())

    ca.Sub_I(p.scratchV1, p.madj)
    #p.diffOp.applyOperator(p.scratchV1)
    return (p.scratchV1, VEnergy, IEnergy)

# end MatchingGradient

def MatchingIteration(p):
    # compute gradient for matching
    (grad_m, VEnergy, IEnergy) = MatchingGradient(p)

    if p.optMethod == 'FIXEDGD':
        # take fixed stepsize gradient step
        ca.Add_MulC_I(p.m0,grad_m,-p.stepSize)
    else:
        raise Exception("Unknown optimization scheme: "+p.optMethod)
    # end if
    return (VEnergy, IEnergy)
# end MatchingIteration


def MatchingPlots(p):
    """
    Do some summary plots for image matching
    """
    #reset all variables by shooting once, may have been overwritten
    CAvmCommon.IntegrateGeodesic(p.m0,p.t,p.diffOp,\
                          p.m, p.g, p.ginv,\
                          p.scratchV1, p.scratchV2,p. scratchV3,\
                          p.checkpointstates, p.checkpointinds,\
                          Ninv=p.nInv, integMethod = p.integMethod)

    # plot the images
    fig = plt.figure('images')
    plt.clf()
    fig.patch.set_facecolor('white')
    plt.subplot(2,2,1)
    display.DispImage(p.I0, 'I0', newFig=False, cmap='gray', sliceIdx=p.plotSlice)
    plt.subplot(2,2,2)
    indx_of_last_tp = p.checkpointinds.index(len(p.t)-1)
    (g, ginv) = p.checkpointstates[indx_of_last_tp]
    ca.ApplyH(p.I,p.I0,ginv)
    display.DispImage(p.I, '\phi.I0', newFig=False, cmap='gray', sliceIdx=p.plotSlice)
    plt.subplot(2,2,3)
    display.DispImage(p.I1, 'I1', newFig=False, cmap='gray', sliceIdx=p.plotSlice)
    plt.subplot(2,2,4)
    ca.ApplyH(p.I,p.I1,g)
    display.DispImage(p.I, '\phi^{-1}.I1', newFig=False, cmap='gray', sliceIdx=p.plotSlice)

    plt.draw()
    plt.show()
    if p.outputPrefix != None: plt.savefig(p.outputPrefix+'images.pdf')
    fig = plt.figure('def')
    plt.clf()
    fig.patch.set_facecolor('white')
    plt.subplot(2,2,1)
    display.GridPlot(ginv,every=p.quiverEvery,color='k', sliceIdx=p.plotSlice, isVF=False)
    plt.axis('equal')
    plt.axis('off')
    plt.title('\phi^{-1}')
    plt.subplot(2,2,2)
    display.GridPlot(g,every=p.quiverEvery,color='k', sliceIdx=p.plotSlice, isVF=False)
    plt.axis('equal')
    plt.axis('off')
    plt.title('\phi')
    plt.subplot(2,2,3)
    ca.JacDetH(p.I,ginv) #p.I used as scratch variable to compute jacobian
    display.DispImage(p.I, '|D\phi^{-1}|', newFig=False, sliceIdx=p.plotSlice)
    plt.subplot(2,2,4)
    ca.MulC_I(p.residualIm, p.sigma*p.sigma)
    display.DispImage(p.residualIm, '\phi.I0-I1', newFig=False, sliceIdx=p.plotSlice)
    plt.colorbar()
    plt.draw()
    plt.show()
    if p.outputPrefix != None: plt.savefig(p.outputPrefix+'def.pdf')
    fig = plt.figure('energy')
    fig.patch.set_facecolor('white')

    TE = [sum(x) for x in p.Energy]
    VE = [row[0] for row in p.Energy]
    IE = [row[1] for row in p.Energy]
    plt.subplot(1,3,1)
    plt.plot(TE)
    plt.title('Total Energy')
    plt.hold(False)
    plt.subplot(1,3,2)
    plt.plot(VE)
    plt.title('Vector Energy')
    plt.hold(False)
    plt.subplot(1,3,3)
    plt.plot(IE)
    plt.title('Image Energy')
    plt.hold(False)
    plt.draw()
    plt.show()
    if p.outputPrefix != None: plt.savefig(p.outputPrefix+'energy.pdf')
# end MatchingPlots

def RunMatching(p):
    p.Energy = []
    for it in xrange(p.nIter):
        (VEnergy, IEnergy) = MatchingIteration(p)
        p.Energy.append([VEnergy, IEnergy])
        print VEnergy+IEnergy, '(Total) = ',VEnergy, '(Vector)+',IEnergy,'(Image)'
# plot some stuff
        # if p.plotEvery > 0 and (((it+1) % p.plotEvery) == 0 or it == p.nIter-1):
            # MatchingPlots(p)
        # end if
    # end for
# end RunMatching

def Matching(cf):

    if cf.compute.useCUDA and cf.compute.gpuID is not None:
        ca.SetCUDADevice(cf.compute.gpuID)
    if os.path.isfile(cf.io.outputPrefix+'m0.mhd'):
        print cf.io.outputPrefix;
        return();

    # if os.path.isfile(cf.io.outputPrefix+'m0.mhd'):
    #     return();
    # else:
    #     print cf.io.outputPrefix;
    #     return();        
    # prepare output directory
    common.Mkdir_p(os.path.dirname(cf.io.outputPrefix))

    # Output loaded config
    if cf.io.outputPrefix is not None:
        cfstr = Config.ConfigToYAML(MatchingConfigSpec, cf)
        with open(cf.io.outputPrefix + "parsedconfig.yaml", "w") as f:
            f.write(cfstr)

    mType = ca.MEM_DEVICE if cf.compute.useCUDA else ca.MEM_HOST


    I0 = common.LoadITKImage(cf.study.I0, mType)
    I1 = common.LoadITKImage(cf.study.I1, mType)
    #ca.DivC_I(I0,255.0)
    #ca.DivC_I(I1,255.0)
    grid = I0.grid()

    It = ca.Image3D(grid,mType)

    ca.ThreadMemoryManager.init(grid, mType, 1)

    #common.DebugHere()
    # TODO: need to work on these
    t = [x*1./cf.optim.nTimeSteps for x in range(cf.optim.nTimeSteps+1)]
    checkpointinds = range(1,len(t))
    checkpointstates =  [(ca.Field3D(grid,mType),ca.Field3D(grid,mType)) for idx in checkpointinds]

    p = MatchingVariables(I0,I1, cf.vectormomentum.sigma, t,checkpointinds, checkpointstates, cf.vectormomentum.diffOpParams[0], cf.vectormomentum.diffOpParams[1], cf.vectormomentum.diffOpParams[2], cf.optim.Niter, cf.optim.stepSize, cf.optim.maxPert, cf.optim.nTimeSteps, integMethod = cf.optim.integMethod, optMethod=cf.optim.method, nInv=cf.optim.NIterForInverse,plotEvery=cf.io.plotEvery, plotSlice = cf.io.plotSlice, quiverEvery = cf.io.quiverEvery, outputPrefix = cf.io.outputPrefix)

    print(p.stepSize)

    RunMatching(p)

    # write output
    if cf.io.outputPrefix is not None:
        # reset all variables by shooting once, may have been overwritten
        CAvmCommon.IntegrateGeodesic(p.m0,p.t,p.diffOp,\
                          p.m, p.g, p.ginv,\
                          p.scratchV1, p.scratchV2,p. scratchV3,\
                          p.checkpointstates, p.checkpointinds,\
                          Ninv=p.nInv, integMethod = p.integMethod)
        ca.ApplyH(It,I0,p.ginv)
        common.SaveITKField(p.m0, cf.io.outputPrefix+"m0.mhd")
        common.SaveITKField(p.ginv, cf.io.outputPrefix+"phiinv.mhd")
        #common.SaveITKField(p.g, cf.io.outputPrefix+"phi.mhd")
        common.SaveITKImage(It, cf.io.outputPrefix+"I1.mhd")
    # end if
# end Matching

if __name__ == '__main__':

    try:
        usercfg = Config.Load(spec=MatchingConfigSpec, argv=sys.argv)
    except Config.MissingConfigError:
        # Don't freak out, this should have printed a sample config to STDOUT
        sys.exit(1)

    Matching(usercfg)

#!/usr/bin/python2
import os.path

# pyca modules
import PyCA.Core as ca
import PyCA.Common as common

# vector momentum modules
from Libraries import CAvmCommon

# others
import numpy as np
import os, errno

import logging
import sys
import copy
import math


def geodesic_shooting_diffOp(moving, target, m0, steps, mType, config):
	grid = moving.grid()
	m0.setGrid(grid)
	ca.ThreadMemoryManager.init(grid, mType, 1)
	if mType == ca.MEM_HOST:
		diffOp = ca.FluidKernelFFTCPU()
	else:
		diffOp = ca.FluidKernelFFTGPU()
	diffOp.setAlpha(config['deformation_params']['diffOpParams'][0])
	diffOp.setBeta(config['deformation_params']['diffOpParams'][1])
	diffOp.setGamma(config['deformation_params']['diffOpParams'][2])
	diffOp.setGrid(grid)

	g = ca.Field3D(grid, mType)
	ginv = ca.Field3D(grid, mType)
	mt = ca.Field3D(grid, mType)
	It = ca.Image3D(grid, mType)
	It_inv = ca.Image3D(grid, mType)

	if (steps <= 0):
		time_steps = config['deformation_params']['timeSteps'];
	else:
		time_steps = steps;

	t = [x*1./time_steps for x in range(time_steps+1)]
	checkpointinds = range(1,len(t))
	checkpointstates =  [(ca.Field3D(grid,mType),ca.Field3D(grid,mType)) for idx in checkpointinds]

	scratchV1 = ca.Field3D(grid,mType)
	scratchV2 = ca.Field3D(grid,mType)
	scratchV3 = ca.Field3D(grid,mType)    

	CAvmCommon.IntegrateGeodesic(m0,t,diffOp, mt, g, ginv,\
								 scratchV1,scratchV2,scratchV3,\
                                 keepstates=checkpointstates,keepinds=checkpointinds,
                                 Ninv=config['deformation_params']['NIterForInverse'], integMethod = config['deformation_params']['integMethod'])
	ca.ApplyH(It,moving,ginv)
	ca.ApplyH(It_inv,target,g)

	output={'I1':It, 'I1_inv': It_inv, "phiinv":ginv}
	return output
#enddef



def geodesic_shooting(moving, target, m0, steps, mType, config):
	return geodesic_shooting_diffOp(moving, target, m0, steps, mType, config)
#enddef 
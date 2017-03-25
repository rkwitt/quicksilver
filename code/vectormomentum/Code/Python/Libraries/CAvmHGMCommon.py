import PyCA.Core as ca
import PyCA.Common as common
from CAvmCommon import *

def HGMComputeVhatResidual(out_v, rhohat, p, phat, diffOp):
     '''
     '''
     CoAdInf(out_v, phat ,p)
     MulC_I(out_v, -1.0)
     Add_I(out_v, rhohat)
     diffOp.applyInverseOperator(out_v)
     return out_v
# end HGMComputeVhatResidual

def HGMEvaluateRHSBwdResidual(out_v, scratchV, v, p, phat, rhohat, rho, diffOp):
    '''
    '''

    #first compute vhat_t
    HGMComputeVhatResidual(out_v, rhohat, p, phat, diffOp)

    #compute Dvmhat
    JacobianXY(scratchV, v, phat)

    #Add both
    Add_I(scratchV, out_v)

    #deform
    ApplyH(out_v, scratchV, rho, BACKGROUND_STRATEGY_PARTIAL_ZERO);

    return out_v
#end HGMEvaluateRHSBwdResidual

class HGMTimeDependentStateGroup:
    def __init__(self, t, grid, mType=MEM_HOST, Jt=None, nt=None, checkpoint=False,subjectId=None):
        self.t = t # time
        self.J = Jt # intercept data (if we have it here)
        self.n = nt # slope data (if we have it here)

        if Jt is not None: # if we have individual data here, we'll also need initial momenta for residual geodesic and store its unique id (name) for writing output
            self.p0 =  Field3D(grid, mType) # initial momenta for residual geodesic
            SetMem(self.p0,0.0)
            self.Energy=[] # energy history for this subject
            if subjectId is not None: #this should always evaluate to true if Jt is not None
                 self.subjectId = subjectId
        else:
            self.p0 = None
        # end if

        # set up checkpointing when requested, and force checkpointing at datapoints
        if checkpoint or Jt is not None:
            self.cpg = Field3D(grid, mType)
            self.cpginv = Field3D(grid, mType)
            SetToIdentity(self.cpg)
            SetToIdentity(self.cpginv)
        else:
            self.cpg = self.cgpinv = None
        # end if
# end TimeDependentStateGroup

class HGMTimeDependentStateResidual:
    def __init__(self, s, grid, mType=MEM_HOST, checkpoint=False):
        self.s = s # time

        # set up checkpointing when requested
        if checkpoint:
            self.cpg = Field3D(grid, mType)
            self.cpginv = Field3D(grid, mType)
            SetToIdentity(self.cpg)
            SetToIdentity(self.cpginv)
        else:
            self.cpg = self.cpginv = None
        # end if
# end TimeDependentStateResidual


class HGMGroupState:
    def __init__(self, grid, mType, alpha, beta, gamma,t, Ninv, SigmaIntercept, SigmaSlope, Sigma, StepSize, integMethod='EULER'):
        self.grid = grid
        self.mType = mType

        # initial conditions
        self.I0 = Image3D(grid,mType)
        SetMem(self.I0,0.0)
        self.m0 = Field3D(grid,mType)
        SetMem(self.m0,0.0)
        self.EnergyHistory=[]

        # state at time t
        self.g = Field3D(grid,mType)
        self.ginv = Field3D(grid,mType)
        self.m = Field3D(self.grid, self.mType)
        self.I = Image3D(self.grid, self.mType)

        # adjoint variables at time t
        self.madj = Field3D(self.grid, self.mType)
        self.Iadj = Image3D(self.grid, self.mType)

        # image variables for closed-form template update
        self.sumSplatI =  Image3D(self.grid, self.mType)
        self.sumJac =  Image3D(self.grid, self.mType)

        # differential operator (size same as a Field3D)
        if mType == MEM_HOST: self.diffOp = FluidKernelFFTCPU()
        else: self.diffOp = FluidKernelFFTGPU()
        self.diffOp.setAlpha(alpha)
        self.diffOp.setBeta(beta)
        self.diffOp.setGamma(gamma)
        self.diffOp.setGrid(grid)

        # some extras
        self.t = t
        self.Ninv = Ninv
        self.integMethod = integMethod
        self.SigmaIntercept = SigmaIntercept
        self.SigmaSlope = SigmaSlope
        self.Sigma = Sigma
        self.StepSize = StepSize
# end HGMGroupState

class HGMResidualState:
    def __init__(self, J0, p0, grid, mType, alpha, beta, gamma,s, Ninv, SigmaIntercept, SigmaSlope, Sigma, StepSize, integMethod='EULER'):
        self.grid = grid
        self.mType = mType

        # initial conditions, do not create new memory these are always references
        self.J0 = J0
        self.p0 = p0

        # state at time s
        self.rho = Field3D(self.grid,self.mType)
        self.rhoinv = Field3D(self.grid,self.mType)
        self.p = Field3D(self.grid, self.mType)
        self.J = Image3D(self.grid, self.mType)

        # adjoint variables at time s
        self.padj = Field3D(self.grid, self.mType)
        self.rhoadj = Field3D(self.grid, self.mType)

        # differential operator (size same as a Field3D)
        if mType == MEM_HOST: self.diffOp = FluidKernelFFTCPU()
        else: self.diffOp = FluidKernelFFTGPU()
        self.diffOp.setAlpha(alpha)
        self.diffOp.setBeta(beta)
        self.diffOp.setGamma(gamma)
        self.diffOp.setGrid(grid)

        # some extras
        self.s = s
        self.Ninv = Ninv
        self.integMethod = integMethod
        self.SigmaIntercept = SigmaIntercept
        self.SigmaSlope = SigmaSlope
        self.Sigma = Sigma
        self.StepSize = StepSize
# end HGMResidualState

def HGMSetupTimeDiscretizationGroup(t,J,n, Jind,checkpointinds,cpMemType=None,subjectIds=None):
    """
    Set up a tdisc object for group geodesic (list of HGMTimeDependentStatesGroup)
    """
    grid = J[0].grid()
    if cpMemType is None:
        cpMemType = J[0].memType()
    td = []
    for i in xrange(len(t)):
        if i in Jind:
            td.append(HGMTimeDependentStateGroup(t[i],grid,cpMemType,J[Jind.index(i)],n[Jind.index(i)],i in checkpointinds,subjectIds[Jind.index(i)]))
        else:
            td.append(HGMTimeDependentStateGroup(t[i],grid,cpMemType,None,None,i in checkpointinds))
    return td
# end HGMSetupTimeDiscretizationGroup

def HGMSetupTimeDiscretizationResidual(s, checkpointinds, cpGrid, cpMemType=None):
    """
    Set up a tdisc object for residual geodesic (list of HGMTimeDependentStatesResidual)
    """

    td = []
    for i in xrange(len(s)):
        td.append(HGMTimeDependentStateResidual(s[i],cpGrid,cpMemType,i in checkpointinds))
    return td
# end HGMSetupTimeDiscretizationResidual

def HGMResidualIter(gradIntercept, gradSlope, residualState, tDiscResidual, mt, J1, n1):
    """
    This integrates the residual geodesic forward and then integrates the residual geodesics adjoint backward.
    Takes a gradient step to update residualState.p0 and returns the intercept and slope gradient.
    """

    # shoot residual geodesic forward
    HGMIntegrateGeodesic(residualState.p0, residualState.s,residualState.diffOp, residualState.p, residualState.rho,residualState.rhoinv, tDiscResidual, residualState.Ninv,residualState.integMethod)

    # integrate residual geodesic backward
    HGMIntegrateAdjointsResidual(residualState, tDiscResidual, mt, J1, n1)

    # compute gradients to return to group geodesic
    (pEnergy, iEnergy, slopeEnergy) = HGMComputeJumpGradientsFromResidual(gradIntercept, gradSlope, residualState, mt, J1, n1)

    # gradient descent step for residualState.p0
    HGMTakeGradientStepResidual(residualState)

    return (pEnergy, iEnergy, slopeEnergy)
# end HGMResidualIter

def HGMComputeJumpGradientsFromResidual(gradIntercept, gradSlope, residualState, mt, J1, n1):
     interceptEnergy = 0.0
     slopeEnergy = 0.0

     grid = residualState.J0.grid()
     mType = residualState.J0.memType()

     imdiff = ca.ManagedImage3D(grid, mType)
     vecdiff = ca.ManagedField3D(grid, mType)

     # 1. Gradient for Intercept
     ApplyH(imdiff,residualState.J0,residualState.rhoinv)
     Sub_I(imdiff, J1)
     iEnergy = 0.5*ca.Sum2(imdiff)/(float(residualState.p0.nVox())*residualState.Sigma*residualState.Sigma*residualState.SigmaIntercept*residualState.SigmaIntercept) # save for use in intercept energy term
     MulC_I(imdiff, 1.0/(residualState.SigmaIntercept*residualState.SigmaIntercept*residualState.Sigma*residualState.Sigma))
     #Splat(gradIntercept, residualState.rhoinv, imdiff,False)
     SplatSafe(gradIntercept, residualState.rhoinv, imdiff)


     # 2. Gradient for Slope
     CoAd(residualState.p,residualState.rhoinv, mt)
     Sub_I(residualState.p, n1)
     Copy(vecdiff, residualState.p) # save for use in slope energy term
     residualState.diffOp.applyInverseOperator(residualState.p)
     # apply Ad and get the result in gradSlope
     Ad(gradSlope,residualState.rhoinv,residualState.p)
     MulC_I(gradSlope, 1.0/(residualState.SigmaSlope*residualState.SigmaSlope))
     # energy computation here
     # slope term
     slopeEnergy = ca.Dot(residualState.p,vecdiff)/(float(residualState.p0.nVox())*residualState.SigmaSlope*residualState.SigmaSlope)

     # intercept term. p is used as scratch variable
     Copy(residualState.p, residualState.p0)
     residualState.diffOp.applyInverseOperator(residualState.p)
     pEnergy = ca.Dot(residualState.p0, residualState.p)/(float(residualState.p0.nVox())*residualState.SigmaIntercept*residualState.SigmaIntercept)

     return (pEnergy,iEnergy, slopeEnergy)

# end HGMComputeJumpGradientsFromResidual
def HGMWriteListToFile(listname,filename):
    try:
        os.remove(filename)
    except OSError:
        pass
    theFile = open(filename, 'w')
    for item in listname:
        print>>theFile, item
    theFile.close()
# end HGMWriteListToFile

def HGMTakeGradientStepResidual(residualState):
     # residualState.p is used as scratch variable
     v = residualState.p
     # Gradient for Momenta
     Copy(v,residualState.p0)
     residualState.diffOp.applyInverseOperator(v)
     MulC_I(v, 1.0/(residualState.SigmaIntercept*residualState.SigmaIntercept))
     Sub_I(v, residualState.padj)
     #residualState.diffOp.applyOperator(v)
     Add_MulC_I(residualState.p0, v,-residualState.StepSize)

# end HGMTakeGradientStepResidual

def HGMIntegrateGeodesic(m0,t,diffOp,\
                             m,g,ginv,\
                             tdisc=None,\
                             Ninv=5,integMethod='EULER', startIndex=0, endIndex=0, initializePhi=True,\
                             scratchV1=None, scratchV2=None, scratchV3=None,\
                             scratchG=None, RK4=None):
    '''
    Resulted g and ginv are  diffeomorphism and its inverse at end point of shooting. and keepstates is
    populated with g and ginv at checkpoints in tdisc. m is NOT the momentum at the end point. Must call CoAd(m,ginv,m0) after the call to get it.
    If startTime is anything other than 0, it assumes g and ginv are appropriately initialized.
    '''
    grid = m0.grid()
    mType = m0.memType()
    # scratch variables
    if scratchV1 is None:
        scratchV1 = ca.ManagedField3D(grid, mType)
    if scratchV2 is None:
        scratchV2 = ca.ManagedField3D(grid, mType)
    if scratchV3 is None:
        scratchV3 = ca.ManagedField3D(grid, mType)
    # RK4 needs two extra scratch variables
    if integMethod == "RK4":
        if scratchG is None:
            scratchG = ca.ManagedField3D(grid, mType)
        if RK4 is None:
            RK4 = ca.ManagedField3D(grid, mType)
    # end if

    # initial conditions
    if endIndex == 0:
         endIndex=len(t)-1

    if startIndex == 0:
         if initializePhi==True: # if t=0 is not the identity diffeomorphisms
              SetToIdentity(g)
              SetToIdentity(ginv)
              Copy(m,m0)
         else:
              CoAd(m,ginv,m0)
    else:
         CoAd(m,ginv,m0) # assumes ginv is initialized when the function was called
    # end if

    # Smooth to get velocity at time t
    diffOp.applyInverseOperator(m)
    if tdisc is not None:
        if (tdisc[0].cpg is not None) & (startIndex == 0):
            # remember to copy everything
            Copy(tdisc[0].cpg,g)
            Copy(tdisc[0].cpginv,ginv)
         # end if
    # end if

    # do integration
    for i in range(startIndex+1,endIndex+1,1):

         dt = t[i]-t[i-1] # time step
         if integMethod == "EULER":
              # print 'Running Euler integration for shooting'
              # Compute forward step, w for g
              EvaluateRHSFwd(scratchV1, m, g)
              MulC_I(scratchV1, dt)

              # Take the fwd step
              Add_I(g,scratchV1)

              # Update ginv corresponding to this fwd step
              UpdateInverse(scratchV3, scratchV2, ginv, scratchV1, Ninv)
              Copy(ginv,scratchV3)

         elif integMethod == "RK4":
              # RK4 integration two extra scratch fields
              EvaluateRHSFwd(scratchV1, m, g); MulC_I(scratchV1, dt) # k1 computed
              Copy(RK4,scratchV1)

              # for k2
              Copy(scratchG, g); MulC_I(scratchV1,0.5); Add_I(scratchG,scratchV1); UpdateInverse(scratchV3, scratchV2, ginv, scratchV1, Ninv); CoAd(m,scratchV3,m0); diffOp.applyInverseOperator(m);
              EvaluateRHSFwd(scratchV1, m, scratchG); MulC_I(scratchV1, dt) # k2 computed
              Add_MulC_I(RK4,scratchV1,2.0)

              # for k3
              Copy(scratchG, g); MulC_I(scratchV1,0.5); Add_I(scratchG,scratchV1); UpdateInverse(scratchV3, scratchV2, ginv, scratchV1, Ninv); CoAd(m,scratchV3,m0); diffOp.applyInverseOperator(m);
              EvaluateRHSFwd(scratchV1, m, scratchG); MulC_I(scratchV1, dt) # k3 computed
              Add_MulC_I(RK4,scratchV1,2.0)

              # for k4
              Copy(scratchG, g);  Add_I(scratchG,scratchV1); UpdateInverse(scratchV3, scratchV2, ginv, scratchV1, Ninv); CoAd(m,scratchV3,m0); diffOp.applyInverseOperator(m);
              EvaluateRHSFwd(scratchV1, m, scratchG); MulC_I(scratchV1, dt) # k4 computed
              Add_I(RK4,scratchV1)

              # final update
              MulC_I(RK4,1.0/float(6.0))
              Add_I(g,RK4)
              UpdateInverse(scratchV3, scratchV2, ginv, RK4, Ninv)
              Copy(ginv,scratchV3)
         else:
              raise Exception('Unknown integration method: '+integMethod)
         # end if

         # common lines of code executed regardless of the integration scheme
         # check whether we should store this state
         if tdisc is not None:
             if (tdisc[i].cpg is not None):
                 # remember to copy everything
                 Copy(tdisc[i].cpg,g)
                 Copy(tdisc[i].cpginv,ginv)
             # end if
         # end if

         # skip if it is the last iteration.
         if i<endIndex:
              # Coadjoint action gives us momentum at time t
              CoAd(m,ginv,m0)
              # Smooth to get velocity at time t
              diffOp.applyInverseOperator(m)
         # end if
    # end for
    '''
    del ScratchV1
    del ScratchV2
    del ScratchV3
    if integMethod == "RK4":
        del ScratchG
        del RK4
    # end if
    '''
# end IntegrateGeodesic

def HGMIntegrateGeodesicBwdIteration(t,i,m0,g1,ginv1,m1,bwdG,bwdGinv,
                                  gprev,ginvprev,
                                  m1initialized,prev_was_checkpoint,
                                  diffOp,
                                  m,
                                  scratchV1, scratchV2,scratchV3,
                                  Ninv=5,integMethod='EULER',
                                  RK4=None, scratchG=None):
     if m1 is None:
         print 'm1 variable is initialized in HGMIntegrateGeodesicBwdIteration'
         m1 = ManagedField3D(m0.grid(),m0.memType())
     if bwdG is None:
         print 'bwdG variable is initialized in HGMIntegrateGeodesicBwdIteration'
         bwdG = ManagedField3D(m0.grid(),m0.memType())
     if bwdGinv is None:
         print 'bwdGinv variable is initialized in HGMIntegrateGeodesicBwdIteration'
         bwdGinv = ManagedField3D(m0.grid(),m0.memType())

     if ( m1initialized == False ):
         SetToIdentity(bwdG)
         SetToIdentity(bwdGinv)
         CoAd(m1,ginv1,m0)
         m1initialized = True
     # end if

     # If previous iteration had a checkpoint bwdG, bwdGinv would have been updated. If not need an updated ones
     if (prev_was_checkpoint == True) & (i != (len(t)-2)):
          ComposeHH(bwdG,gprev,ginv1)
          ComposeHH(bwdGinv,g1,ginvprev)
     # end if

     HGMIntegrateGeodesic(m1,[0,t[i]-t[i+1]],diffOp,\
                              m,bwdG,bwdGinv,\
                              tdisc=None,\
                              Ninv=Ninv, integMethod=integMethod, initializePhi=False,\
                              scratchG=scratchG,RK4=RK4)

     ComposeHH(gprev,bwdG,g1)
     ComposeHH(ginvprev, ginv1,bwdGinv)
     prev_was_checkpoint = False

# end HGMIntegrateGeodesicBwdIteration


def HGMIntegrateAdjoints(groupState, tDiscGroup, residualState, tDiscResidual):
    '''
    groupState and tDiscGroup is for the group geodesic while residualState and tDiscResidual are re-used for all residuals as scratch data structures
    '''

    grid = groupState.I0.grid()
    mType = groupState.I0.memType()

    # allocate scratch memory for this function
    scratchV1 = ca.ManagedField3D(grid, mType)
    scratchV2 = ca.ManagedField3D(grid, mType)

    v = ca.ManagedField3D(grid,mType)
    mhat = ca.ManagedField3D(grid,mType)
    Ihat = ca.ManagedImage3D(grid,mType)
    imOnes = ca.ManagedImage3D(grid,mType)
    SetMem(imOnes,1.0)

    # extra reference names used
    m1 = None
    bwdG = None
    bwdGinv = None
    g = None
    ginv = None

    SetMem(groupState.madj,0.0)
    SetMem(groupState.Iadj,0.0)
    indx_of_last_tp = len(groupState.t)-1
    g = tDiscGroup[indx_of_last_tp].cpg
    ginv = tDiscGroup[indx_of_last_tp].cpginv

    # initialize variables for closed form template update
    SetMem(groupState.sumSplatI,0.0)
    SetMem(groupState.sumJac,0.0)

    #I(t) and m(t) at end point
    CoAd(groupState.m,ginv,groupState.m0)
    ApplyH(groupState.I,groupState.I0,ginv)
    Copy(v,groupState.m)
    groupState.diffOp.applyInverseOperator(v) # has v(t)
    SetMem(mhat,0.0) # will be used for hat version of madj
    SetMem(Ihat,0.0) # will be used for hat version of Iadj
    # initial conditions for backward integration
    # check if there is an individual at last timepoint

    if tDiscGroup[indx_of_last_tp].J is not None:
        residualState.J0 = groupState.I
        residualState.p0 = tDiscGroup[indx_of_last_tp].p0
        # get slope gradients and intercept gradients and also update initial momenta p0 for residual geodesic rho as per gradient descent
        tDiscGroup[indx_of_last_tp].Energy.append(HGMResidualIter(Ihat, mhat, residualState,tDiscResidual, groupState.m, tDiscGroup[indx_of_last_tp].J,tDiscGroup[indx_of_last_tp].n))
        # accumulate splatI and jacI (groupState.ginv used as scratch)
        ComposeHH(groupState.ginv,ginv,residualState.rhoinv)
        scratchI = ca.ManagedImage3D(grid, mType)
        #Splat(scratchI, groupState.ginv, tDiscGroup[indx_of_last_tp].J,False)
        SplatSafe(scratchI, groupState.ginv, tDiscGroup[indx_of_last_tp].J)
        Add_I(groupState.sumSplatI, scratchI)
        #Splat(scratchI, groupState.ginv, imOnes, False)
        SplatSafe(scratchI, groupState.ginv, imOnes)
        Add_I(groupState.sumJac, scratchI)
        del scratchI

        MulC_I(Ihat,-1)
        MulC_I(mhat,-1)
        # splat the negative of the intercept gradient backward
        #Splat(groupState.Iadj, ginv, Ihat,False)
        SplatSafe(groupState.Iadj, ginv, Ihat)
        # deform the negative of the slope gradient forward
        ApplyH(groupState.madj, mhat, g, BACKGROUND_STRATEGY_PARTIAL_ZERO);
    # end if

    prev_was_checkpoint = True
    m1initialized = False
    # backward integrate
    for i in range(len(groupState.t)-2,-1,-1):
        dt = groupState.t[i] - groupState.t[i+1]
        if groupState.integMethod == "EULER":
            EvaluateRHSBwd(scratchV1, scratchV2, v, groupState.I, groupState.m, mhat, Ihat, g, groupState.diffOp)

            # Euler step for madj
            Add_MulC_I(groupState.madj, scratchV1, dt)
            indx_of_cur_tp = i
            if (tDiscGroup[i].cpg is not None):
                 g = tDiscGroup[indx_of_cur_tp].cpg
                 ginv = tDiscGroup[indx_of_cur_tp].cpginv
                 prev_was_checkpoint = True
            elif i > 0:# i=0 need not be calculated
                 # compute g and ginv by backward integrating geodesic
                 # oops we are going to need a lot of scratch variables
                 scratchV3 = ca.ManagedField3D(grid, mType) # for m1
                 scratchV4 = ca.ManagedField3D(grid, mType) # for bwdG
                 scratchV5 = ca.ManagedField3D(grid, mType) # for bwdGinv
                 scratchV6 = ca.ManagedField3D(grid, mType) # for g
                 scratchV7 = ca.ManagedField3D(grid, mType) # for ginv
                 Copy(scratchV6,g)
                 Copy(scratchV7,ginv)
                 # update reference
                 g=scratchV6
                 ginv=scratchV7

                 # mhat and v are used as scratch variables in below call
                 m1 = scratchV3; bwdG = scratchV4; bwdGinv = scratchV5; # update references for ease of reading

                 HGMIntegrateGeodesicBwdIteration(groupState.t,i,groupState.m0, tDiscGroup[indx_of_last_tp].cpg, tDiscGroup[indx_of_last_tp].cpginv, m1, bwdG, bwdGinv,\
                                                      g,ginv,\
                                                      m1initialized,prev_was_checkpoint,\
                                                      groupState.diffOp,\
                                                      mhat,\
                                                      scratchV1,scratchV2,v,\
                                                      Ninv=groupState.Ninv,integMethod=groupState.integMethod)


            # end if

            # check if there is an individual at this timepoint
            if i>0:
                 CoAd(groupState.m,ginv,groupState.m0)
            else:
                 Copy(groupState.m,groupState.m0)
            if tDiscGroup[indx_of_cur_tp].J is not None:
                 if i>0:
                      ApplyH(groupState.I,groupState.I0,ginv)
                 else:
                      Copy(groupState.I,groupState.I0)
                 residualState.J0 = groupState.I
                 residualState.p0 = tDiscGroup[indx_of_cur_tp].p0
                 # get slope gradients and intercept gradients and also update initial momenta p0 for residual geodesic rho as per gradient descent
                 tDiscGroup[indx_of_cur_tp].Energy.append(HGMResidualIter(Ihat, mhat, residualState,tDiscResidual,groupState.m, tDiscGroup[indx_of_cur_tp].J,tDiscGroup[indx_of_cur_tp].n))

                 # accumulate splatI and jacI (groupState.ginv used as scratch)

                 ComposeHH(groupState.ginv,ginv,residualState.rhoinv)
                 scratchI = ca.ManagedImage3D(grid, mType)
                 #Splat(scratchI, groupState.ginv, tDiscGroup[indx_of_cur_tp].J,False)
                 SplatSafe(scratchI, groupState.ginv, tDiscGroup[indx_of_cur_tp].J)
                 Add_I(groupState.sumSplatI, scratchI)
                 #Splat(scratchI, groupState.ginv, imOnes, False)
                 SplatSafe(scratchI, groupState.ginv, imOnes)
                 Add_I(groupState.sumJac, scratchI)
                 del scratchI

                 if i > 0:
                     # splat the the intercept gradient backward and subtract
                     #Splat(groupState.I, ginv, Ihat, False)
                     SplatSafe(groupState.I, ginv, Ihat)
                     Sub_I(groupState.Iadj,groupState.I)
                     # deform the slope gradient forward and subtract
                     ApplyH(scratchV1, mhat, g, BACKGROUND_STRATEGY_PARTIAL_ZERO);
                     Sub_I(groupState.madj,scratchV1)
                 else:
                     Sub_I(groupState.Iadj,Ihat)
                     Sub_I(groupState.madj,mhat)
                 # end if
            # end if

            # update variables for next iteration
            if i > 0: # last iteration skipped
                 ApplyH(groupState.I,groupState.I0,ginv)
                 Copy(v,groupState.m)
                 groupState.diffOp.applyInverseOperator(v) # has v(t)
                 ApplyH(mhat,groupState.madj,ginv,BACKGROUND_STRATEGY_PARTIAL_ZERO) # hat version of madj
                 #Splat(Ihat, g, groupState.Iadj,False)  # hat version of Iadj
                 SplatSafe(Ihat, g, groupState.Iadj)  # hat version of Iadj
            # end if
        elif groupState.integMethod == "RK4":
            raise Exception('RK4 not yet implemented for HGM group geodesic backward integration.')
        else:
            raise Exception('Unknown integration method: '+integMethod)
        #end if
    # end for

    # Gather all individual energies
    TotalInterceptEnergy = 0.0
    TotalSlopeEnergy = 0.0
    for tdisc in tDiscGroup:
         if tdisc.J is not None:
              TotalInterceptEnergy += tdisc.Energy[-1][0]+ tdisc.Energy[-1][1]
              TotalSlopeEnergy += tdisc.Energy[-1][2]

    # compute vector energy    (m is used as scratch)
    Copy(groupState.m,groupState.m0)
    groupState.diffOp.applyInverseOperator(groupState.m)
    VectorEnergy = Dot(groupState.m0,groupState.m)/float(groupState.m0.nVox())

    # append all energies in history
    groupState.EnergyHistory.append([VectorEnergy, TotalInterceptEnergy, TotalSlopeEnergy])

#end HGMIntegrateAdjoints

##########
def HGMIntegrateAdjointsResidual(residualState, tDiscResidual, mt ,J1, n1):
    '''
    '''
    grid = residualState.J0.grid()
    mType = residualState.J0.memType()

    # allocate scratch memory for this function
    scratchV1 = ca.ManagedField3D(grid, mType)
    scratchV2 = ca.ManagedField3D(grid, mType)
    v = ca.ManagedField3D(grid,mType)
    phat = ca.ManagedField3D(grid,mType)
    rhohat = ca.ManagedField3D(grid,mType)
    if residualState.integMethod == "RK4":
         RK4 = ManagedField3D(grid,mType)
         scratchG = ManagedField3D(grid,mType)
         scratchGinv = ManagedField3D(grid,mType)

    # extra reference names used
    p1 = None
    bwdRho = None
    bwdRhoinv = None
    rho = None
    rhoinv = None

    SetMem(residualState.padj,0.0)
    SetMem(residualState.rhoadj,0.0)
    indx_of_last_tp = len(residualState.s)-1
    rho = tDiscResidual[indx_of_last_tp].cpg
    rhoinv = tDiscResidual[indx_of_last_tp].cpginv

    #I(t) and m(t) at end point
    CoAd(residualState.p,rhoinv, residualState.p0)
    ApplyH(residualState.J,residualState.J0,rhoinv)
    Copy(v,residualState.p)
    residualState.diffOp.applyInverseOperator(v) # has v(t)
    SetMem(phat,0.0) # will be used for hat version of padj
    SetMem(rhohat,0.0) # will be used for hat version of rhoadj

    # initial conditions for backward integration
    # first term
    Gradient(rhohat, residualState.J)
    Sub_I(residualState.J,J1) # imagediff
    MulMulC_I(rhohat, residualState.J, -1.0/(residualState.SigmaIntercept*residualState.SigmaIntercept*residualState.Sigma*residualState.Sigma))

    # second term
    CoAd(scratchV1,rhoinv, mt)
    Sub(scratchV2, scratchV1, n1)
    residualState.diffOp.applyInverseOperator(scratchV2)
    CoAdInf(residualState.rhoadj, scratchV2, scratchV1)  # residualState.rhoadj used as scratch variable
    MulC_I(residualState.rhoadj, 1.0/(residualState.SigmaSlope*residualState.SigmaSlope))
    # add the two terms
    Add_I(rhohat,residualState.rhoadj)
    CoAd(residualState.rhoadj, rho, rhohat)

    prev_was_checkpoint = True
    m1initialized = False
    for i in range(len(residualState.s)-2,-1,-1):
        ds = residualState.s[i] - residualState.s[i+1]
        if residualState.integMethod == "EULER":
             # print 'Running Euler integration for adjoint integration'
             HGMEvaluateRHSBwdResidual(scratchV1, scratchV2, v, residualState.p, phat, rhohat, rho, residualState.diffOp)
             # Euler step for padj
             Add_MulC_I(residualState.padj, scratchV1, ds)
             indx_of_cur_tp = i
             if (tDiscResidual[i].cpg is not None):
                 rho = tDiscResidual[indx_of_cur_tp].cpg
                 rhoinv = tDiscResidual[indx_of_cur_tp].cpginv
                 prev_was_checkpoint = True
             elif i > 0:# i=0 need not be calculated
                 # compute g and ginv by backward integrating geodesic
                 # oops we are going to need a lot of scratch variables
                 scratchV3 = ca.ManagedField3D(grid, mType) # for m1
                 scratchV4 = ca.ManagedField3D(grid, mType) # for bwdG
                 scratchV5 = ca.ManagedField3D(grid, mType) # for bwdGinv
                 scratchV6 = ca.ManagedField3D(grid, mType) # for g
                 scratchV7 = ca.ManagedField3D(grid, mType) # for ginv
                 Copy(scratchV6,rho)
                 Copy(scratchV7,rhoinv)
                 # update reference
                 rho=scratchV6
                 rhoinv=scratchV7
                 # mhat and v are used as scratch variables in below call
                 p1 = scratchV3; bwdRho = scratchV4; bwdRhoinv = scratchV5; # update references for ease of reading
                 HGMIntegrateGeodesicBwdIteration(residualState.s,i,residualState.p0, tDiscResidual[indx_of_last_tp].cpg, tDiscResidual[indx_of_last_tp].cpginv, p1, bwdRho, bwdRhoinv,\
                                                      rho,rhoinv,\
                                                      m1initialized,prev_was_checkpoint,\
                                                      residualState.diffOp,\
                                                      phat,\
                                                      scratchV1,scratchV2,v,\
                                                      Ninv=residualState.Ninv,integMethod=residualState.integMethod)

             # end if

             # update variables for next iteration
             if i > 0: # last iteration skipped
                 CoAd(residualState.p,rhoinv,residualState.p0)
                 Copy(v,residualState.p)
                 residualState.diffOp.applyInverseOperator(v) # has v(t)
                 ApplyH(phat,residualState.padj,rhoinv,BACKGROUND_STRATEGY_PARTIAL_ZERO) # hat version of madj
                 CoAd(rhohat,rhoinv,residualState.rhoadj)  # hat version of rhoadj
             # end if
        elif residualState.integMethod == "RK4":

             # first get g and ginv for current timepoint
             if (tDiscResidual[i].cpg is not None) or (i==0): # assuming g and g inv points to prev
                 indx_of_cur_tp = i
                 # just assign the references
                 if i>0:
                      rhocur = tDiscResidual[indx_of_cur_tp].cpg
                      rhoinvcur = tDiscResidual[indx_of_cur_tp].cpginv
                 # end if note if i==0, gcur and ginvcur are treated as identity

                 prev_was_checkpoint = True

                 # begin rk4 integration for adjoint
                 HGMEvaluateRHSBwdResidual(scratchV1, scratchV2, v,  residualState.p, phat, rhohat, rho, residualState.diffOp); MulC_I(scratchV1, ds) # k1 computed
                 Copy(RK4,scratchV1)

                 # compute and store phi0t, phit0, v_t, i_t, m_t and i_hat at t+1/2 for computing k2 and k3
                 if i>0:
                      SubMulC(scratchG,rho,rhocur,0.5); #scratchG has w
                      UpdateInverse(scratchGinv, scratchV2, rhoinvcur, scratchG, residualState.Ninv)
                      Add_I(scratchG,rhocur)

                 else:# add 0.5 times identity
                      HtoV(scratchG,rho);MulC_I(scratchG,0.5) #scratchG has w
                      # iteratively update inverse as, g_{1,0} = Id - w\circ g_{1,0}
                      SetToIdentity(scratchGinv)
                      for k in range(residualState.Ninv):
                           ApplyH(scratchV2, scratchG, scratchGinv, BACKGROUND_STRATEGY_PARTIAL_ZERO)
                           HtoV(scratchGinv,scratchV2); MulC_I(scratchGinv,-1.0)
                      VtoH_I(scratchG)
                 # end if

                 CoAd(residualState.p,scratchGinv,residualState.p0); Copy(v,residualState.p); residualState.diffOp.applyInverseOperator(v);
                 CoAd(rhohat,scratchGinv, residualState.rhoadj)  # hat version of rhoadj

                 # for k2
                 # mhat_t at t+1/2 for k2
                 Add_MulC(scratchV2,residualState.padj,scratchV1,0.5) # mtilde at t+1/2
                 ApplyH(phat,scratchV2,scratchGinv, BACKGROUND_STRATEGY_PARTIAL_ZERO); #mhat at t+1/2
                 HGMEvaluateRHSBwdResidual(scratchV1, scratchV2, v, residualState.p, phat, rhohat, scratchG, residualState.diffOp); MulC_I(scratchV1, ds) # k2 computed
                 Add_MulC_I(RK4, scratchV1, 2.0)

                 # for k3
                 # mhat_t at t+1/2 for k3
                 Add_MulC(scratchV2,residualState.padj,scratchV1,0.5) # mtilde at t+1/2
                 ApplyH(phat,scratchV2,scratchGinv, BACKGROUND_STRATEGY_PARTIAL_ZERO); #mhat at t+1/2

                 HGMEvaluateRHSBwdResidual(scratchV1, scratchV2, v,  residualState.p, phat, rhohat, scratchG, residualState.diffOp); MulC_I(scratchV1, ds) # k3 computed
                 Add_MulC_I(RK4, scratchV1, 2.0)

                 # compute v_t, i_t, m_t, i_hat at t for computing k4
                 if i>0:
                      CoAd(residualState.p,rhoinvcur,residualState.p0)
                      CoAd(rhohat,rhoinvcur, residualState.rhoadj)
                 else:
                      Copy(residualState.p,residualState.p0)
                      Copy(rhohat, residualState.rhoadj)

                 Copy(v,residualState.p); residualState.diffOp.applyInverseOperator(v);
                 # for k4
                 # mhat_t at t for k4
                 Add(scratchV2,residualState.padj,scratchV1) # mtilde at t
                 if i>0:
                      ApplyH(phat,scratchV2,rhoinvcur, BACKGROUND_STRATEGY_PARTIAL_ZERO);
                 else:
                      Copy(phat,scratchV2)
                 # end if #mhat at t

                 if i>0:
                      HGMEvaluateRHSBwdResidual(scratchV1, scratchV2, v, residualState.p, phat, rhohat, rhocur, residualState.diffOp); MulC_I(scratchV1, ds) # k4 computed
                 else:
                      SetToIdentity(scratchG)
                      HGMEvaluateRHSBwdResidual(scratchV1, scratchV2, v, residualState.p, phat, rhohat, scratchG, residualState.diffOp); MulC_I(scratchV1, ds) # k4 computed
                 Add_I(RK4, scratchV1)

                 # final update
                 MulC_I(RK4,1.0/float(6.0))
                 Add_I(residualState.padj,RK4)

                 #FOR NEXT ITERATION:
                 # compute mhat_t, ihat_t at t to use in k1 computation. Note v_t, i_t and m_t are still stored from this iteration.

                 if i > 0: # last iteration skipped
                      ApplyH(phat,residualState.padj,rhoinvcur,BACKGROUND_STRATEGY_PARTIAL_ZERO) # hat version of madj
                      CoAd(rhohat,rhoinvcur,residualState.rhoadj)  # hat version of rhoadj
                      # assign g, ginv appropriately for next iteration
                      rho = rhocur
                      rhoinv = rhoinvcur
                 # end if

             else:
                  raise Exception('RK4 without all checkpoints for HGMIntegrateAdjointsResidual not yet implemented!')
             # end if

        else:
            raise Exception('Unknown integration method: '+integMethod)
        #end if

        # common lines of code executed regardless of the integration scheme

    # end for

# end HGMIntegrateAdjointsResidual

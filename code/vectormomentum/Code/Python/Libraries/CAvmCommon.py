from PyCA.Core import *
import PyCA.Common as common

import numpy as np
import matplotlib.pyplot as plt

def SplatSafe(outMass, g, mass):
     mType = mass.memType()
     if mType == MEM_DEVICE:
          minmaxl = MinMax(mass)
          maxval = max([abs(x) for x in minmaxl])
          if maxval > 2000.00:
               # print 'Warning, range too big for splatting. Range: ',minmaxl
               # print 'Temporary downscaling values for splatting. Will scale it back after splatting.'
               scalefactor =  float(100.00)/maxval
               MulC_I(mass, scalefactor)
               Splat(outMass, g, mass, False)
               MulC_I(mass, 1.0/scalefactor)
               MulC_I(outMass, 1.0/scalefactor)
          else:
               Splat(outMass, g, mass, False)
     else:
          Splat(outMass,g, mass, False)
# end SplatSafe

def ComputeVhat(out_v, scratchV, I, Ihat, m, mhat, diffOp):
     '''
     '''
     Gradient(out_v, I)
     MulMulC_I(out_v, Ihat, 1.0)
     CoAdInf(scratchV, mhat ,m)
     Sub_I(out_v, scratchV)
     diffOp.applyInverseOperator(out_v)
     return out_v
# end ComputeVhat

def EvaluateRHSFwd(out_v, v, phi):
     '''
     Evaluate RHS for forward integration of \frac{d\phi}{dt}=v\circ\phi
     '''
     ApplyH(out_v, v, phi, BACKGROUND_STRATEGY_PARTIAL_ZERO)
     return out_v
#end EvaluateRHSFwd

def EvaluateRHSBwd(out_v, scratchV, v, I, m, mhat, Ihat, g, diffOp):
    '''
    '''
    #first compute vhat_t
    ComputeVhat(out_v, scratchV, I, Ihat, m, mhat, diffOp)

    #compute Dvmhat
    JacobianXY(scratchV, v, mhat)

    #Add both
    Add_I(scratchV, out_v)

    #deform
    ApplyH(out_v, scratchV, g, BACKGROUND_STRATEGY_PARTIAL_ZERO);

    return out_v
#end EvaluateRHSBwd

def IntegrateGeodesic(m0,t,diffOp,\
                           m,g,ginv,\
                           scratchV1,scratchV2,scratchV3,\
                           keepstates=None,keepinds=None,\
                           Ninv=5,integMethod='EULER', startIndex=0, endIndex=0, initializePhi=True,\
                           RK4=None, scratchG=None):
    '''
    Resulted g and ginv are  diffeomorphism and its inverse at end point of shooting. and keepstates is
    populated with g and ginv at requested timepoints mentioned in keepinds. m is NOT the momentum at the end point. Must call CoAd(m,ginv,m0) after the call to get it.
    If startTime is anything other than 0, it assumes g and ginv are appropriately initialized.
    '''
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
    if (keepinds!=None) & (keepstates!=None):
         if (0 in keepinds) & (startIndex == 0):
              # remember to copy everything
              indx_of_cur_tp = keepinds.index(0)
              Copy(keepstates[indx_of_cur_tp][0],g)
              Copy(keepstates[indx_of_cur_tp][1],ginv)
         # end if
    # end if

    # do integration
    for i in range(startIndex+1,endIndex+1,1):
        #sys.stdout.write(',')
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
              # RK4 integration two extra fields
              if scratchG is None:
                   print 'scratchG variable is initialized in geodesic shooting'
                   scratchG = Field3D(m0.grid(),m0.memType())
              if RK4 is None:
                   print 'RK4 variable is initialized in geodesic shooting'
                   RK4 = Field3D(m0.grid(),m0.memType())

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
         if (keepinds!=None) & (keepstates!=None):
              if i in keepinds:
                   # remember to copy everything
                   indx_of_cur_tp = keepinds.index(i)
                   Copy(keepstates[indx_of_cur_tp][0],g)
                   Copy(keepstates[indx_of_cur_tp][1],ginv)
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
# end IntegrateGeodesic

def IntegrateGeodesicBwdIteration(t,i,m0,g1,ginv1,m1,bwdG,bwdGinv,
                                  gprev,ginvprev,
                                  m1initialized,prev_was_checkpoint,
                                  diffOp,
                                  m,
                                  scratchV1, scratchV2,scratchV3,
                                  Ninv=5,integMethod='EULER',
                                  RK4=None, scratchG=None):
     if m1 is None:
          print 'm1 variable is initialized in IntegrateGeodesicBwdIteration'
          m1 = Field3D(m0.grid(),m0.memType())
     if bwdG is None:
          print 'bwdG variable is initialized in IntegrateGeodesicBwdIteration'
          bwdG = Field3D(m0.grid(),m0.memType())
     if bwdGinv is None:
          print 'bwdGinv variable is initialized in IntegrateGeodesicBwdIteration'
          bwdGinv = Field3D(m0.grid(),m0.memType())

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

     IntegrateGeodesic(m1,[0,t[i]-t[i+1]],diffOp,\
                            m,bwdG,bwdGinv,\
                            scratchV1,scratchV2,scratchV3,\
                            keepstates=None,keepinds=None,\
                            Ninv=Ninv,integMethod=integMethod,initializePhi=False,RK4=RK4,scratchG=scratchG)

     ComposeHH(gprev,bwdG,g1)
     ComposeHH(ginvprev, ginv1,bwdGinv)
     prev_was_checkpoint = False

# end IntegrateGeodesicBwd
def IntegrateAdjoints(Iadj,madj,
                      I,m,Iadjtmp, madjtmp,v,
                      scratchV1,scratchV2,
                      I0,m0,
                      t, checkpointstates, checkpointinds,
                      IGradAtMsmts, msmtInds,
                      diffOp,
                      integMethod='EULER',Ninv=5,
                      scratchV3=None, scratchV4=None,scratchV5=None,scratchV6=None, scratchV7=None, # used when all timepoints are not checkpointed or with RK4
                      scratchV8=None, scratchV9=None, # used with RK4 only when all are not checkpointed
                      RK4 = None, scratchG = None,scratchGinv = None, # used only with RK4
                      scratchI = None
                      ):
    '''
    '''
    if len(t)-1 not in checkpointinds:
        raise Exception('Endpoint must be one of the checkpoints passed to IntegrateAdjoints')
    else:
        indx_of_last_tp = checkpointinds.index(len(t)-1)

    # extra reference names used
    m1 = None
    bwdG = None
    bwdGinv = None

    SetMem(madj,0.0)
    SetMem(Iadj,0.0)
    (g, ginv) = checkpointstates[indx_of_last_tp]

    #I(t) and m(t) at end point
    CoAd(m,ginv,m0)
    ApplyH(I,I0,ginv)

    Copy(v,m)
    diffOp.applyInverseOperator(v) # has v(t)
    SetMem(madjtmp,0.0) # will be used for hat version of madj
    SetMem(Iadjtmp,0.0) # will be used for hat version of Iadj

    # initial conditions
    for k in range(len(msmtInds)):
        if checkpointinds[msmtInds[k]] == (len(t)-1):
            # there is a measurement at the last time point which will always be the case for matching but not necessarily regression
            MulC(Iadjtmp,IGradAtMsmts[k],-1)
            #Splat(Iadj,ginv, Iadjtmp,False) #Also, ginv = checkpointstates[msmtInds[k]][1]
            SplatSafe(Iadj, ginv, Iadjtmp)
        # end if
    # end for

    prev_was_checkpoint = True
    m1initialized = False
    for i in range(len(t)-2,-1,-1):
        dt = t[i] - t[i+1]
        if integMethod == "EULER":
             # print 'Running Euler integration for adjoint integration'
             EvaluateRHSBwd(scratchV1, scratchV2, v, I, m, madjtmp, Iadjtmp, g, diffOp)

             # Euler step for madj
             Add_MulC_I(madj, scratchV1, dt)

             if i in checkpointinds:
                 indx_of_cur_tp = checkpointinds.index(i)
                 (g, ginv) = checkpointstates[indx_of_cur_tp]
                 prev_was_checkpoint = True
             elif i > 0: # i=0 need not be calculated
                      # compute g and ginv by backward integrating geodesic
                      # oops we are going to need a lot of scratch variables
                      if scratchV6 is None:
                           print 'scratchV6 variable is initialized in IntegrateAdjoints'
                           scratchV6 = Field3D(I0.grid(),I0.memType())
                      if scratchV7 is None:
                           print 'scratchV7 variable is initialized in IntegrateAdjoints'
                           scratchV7 = Field3D(I0.grid(),I0.memType())
                      if  (prev_was_checkpoint == True): # so that we do not update checkpointed states
                           Copy(scratchV6,g)
                           Copy(scratchV7,ginv)
                           # update reference
                           g=scratchV6
                           ginv=scratchV7
                      # madjtmp and v are used as scratch variables in below call
                      m1 = scratchV3; bwdG = scratchV4; bwdGinv = scratchV5; # update references for ease of reading
                      IntegrateGeodesicBwdIteration(t,i,m0, checkpointstates[indx_of_last_tp][0], checkpointstates[indx_of_last_tp][1],m1, bwdG, bwdGinv,
                                                    g,ginv,
                                                    m1initialized,prev_was_checkpoint,
                                                    diffOp,
                                                    madjtmp,
                                                    scratchV1,scratchV2,v,
                                                    Ninv=Ninv,integMethod=integMethod)

             # end if

             # if there is a measurement at this time point (for regression)

             for k in range(len(msmtInds)):
                  if i>0:
                       if checkpointinds[msmtInds[k]] == i:
                            # I is used as scratch variable
                            #Splat(I, ginv, IGradAtMsmts[k],False)
                            SplatSafe(I, ginv, IGradAtMsmts[k])
                            Sub_I(Iadj,I)
                            # end if
                  elif msmtInds[k]==-1:# if there is a measurement at time t=0 it won't be checkpointed but HARDCODED to have msmstInds == -1. Note this will be checked only for t==0
                       Sub_I(Iadj, IGradAtMsmts[k])
                  # end if
             # end for

             # update variables for next iteration
             if i > 0: # last iteration skipped
                 CoAd(m,ginv,m0)
                 ApplyH(I,I0,ginv)
                 Copy(v,m)
                 diffOp.applyInverseOperator(v) # has v(t)
                 ApplyH(madjtmp,madj,ginv,BACKGROUND_STRATEGY_PARTIAL_ZERO) # hat version of madj
                 #Splat(Iadjtmp, g, Iadj,False)  # hat version of Iadj
                 SplatSafe(Iadjtmp, g, Iadj)  # hat version of Iadj
             # end if
        elif integMethod == "RK4":
             if RK4 is None:
                  print 'RK4 variable is initialized'
                  RK4 = Field3D(I0.grid(),I0.memType())
             if scratchG is None:
                  print 'scratchG variable is initialized'
                  scratchG = Field3D(I0.grid(),I0.memType())
             if scratchGinv is None:
                  print 'scratchGinv variable is initialized'
                  scratchGinv = Field3D(I0.grid(),I0.memType())
             if scratchI is None:
                  print 'scratchI variable is initialized'
                  scratchI = Image3D(I0.grid(),I0.memType())

             # first get g and ginv for current timepoint
             if (i in checkpointinds) or (i==0): # assuming g and g inv points to prev
                 # just assign the references
                 if i>0:
                      indx_of_cur_tp = checkpointinds.index(i)
                      (gcur, ginvcur) = checkpointstates[indx_of_cur_tp]
                 # end if note if i==0, gcur and ginvcur are treated as identity

                 prev_was_checkpoint = True

                 # begin rk4 integration for adjoint
                 EvaluateRHSBwd(scratchV1, scratchV2, v, I, m, madjtmp, Iadjtmp, g, diffOp); MulC_I(scratchV1, dt) # k1 computed
                 Copy(RK4,scratchV1)

                 # compute and store phi0t, phit0, v_t, i_t, m_t and i_hat at t+1/2 for computing k2 and k3
                 if i>0:
                      SubMulC(scratchG,g,gcur,0.5); #scratchG has w
                      UpdateInverse(scratchGinv, scratchV2, ginvcur, scratchG, Ninv)
                      Add_I(scratchG,gcur)

                 else:# add 0.5 times identity
                      HtoV(scratchG,g);MulC_I(scratchG,0.5) #scratchG has w
                      # iteratively update inverse as, g_{1,0} = Id - w\circ g_{1,0}
                      SetToIdentity(scratchGinv)
                      for k in range(Ninv):
                           ApplyH(scratchV2, scratchG, scratchGinv, BACKGROUND_STRATEGY_PARTIAL_ZERO)
                           HtoV(scratchGinv,scratchV2); MulC_I(scratchGinv,-1.0)
                      VtoH_I(scratchG)
                 # end if

                 CoAd(m,scratchGinv,m0); Copy(v,m); diffOp.applyInverseOperator(v); ApplyH(I,I0,scratchGinv)
                 #Splat(Iadjtmp, scratchG, Iadj,False)
                 SplatSafe(Iadjtmp, scratchG, Iadj)

                 # for k2
                 # mhat_t at t+1/2 for k2
                 Add_MulC(scratchV2,madj,scratchV1,0.5) # mtilde at t+1/2
                 ApplyH(madjtmp,scratchV2,scratchGinv, BACKGROUND_STRATEGY_PARTIAL_ZERO); #mhat at t+1/2

                 EvaluateRHSBwd(scratchV1, scratchV2, v, I, m, madjtmp, Iadjtmp, scratchG, diffOp); MulC_I(scratchV1, dt) # k2 computed
                 Add_MulC_I(RK4, scratchV1, 2.0)

                 # for k3
                 # mhat_t at t+1/2 for k3
                 Add_MulC(scratchV2,madj,scratchV1,0.5) # mtilde at t+1/2
                 ApplyH(madjtmp,scratchV2,scratchGinv, BACKGROUND_STRATEGY_PARTIAL_ZERO); #mhat at t+1/2

                 EvaluateRHSBwd(scratchV1, scratchV2, v, I, m, madjtmp, Iadjtmp, scratchG, diffOp); MulC_I(scratchV1, dt) # k3 computed
                 Add_MulC_I(RK4, scratchV1, 2.0)

                 # compute v_t, i_t, m_t, i_hat at t for computing k4
                 if i>0:
                      CoAd(m,ginvcur,m0)
                      #Splat(Iadjtmp, gcur, Iadj,False)
                      SplatSafe(Iadjtmp, gcur, Iadj)
                      ApplyH(I,I0,ginvcur)
                 else:
                      Copy(m,m0)
                      Copy(Iadjtmp,Iadj)
                      Copy(I,I0)
                 Copy(v,m); diffOp.applyInverseOperator(v);
                 # for k4
                 # mhat_t at t for k4
                 Add(scratchV2,madj,scratchV1) # mtilde at t
                 if i>0:
                      ApplyH(madjtmp,scratchV2,ginvcur, BACKGROUND_STRATEGY_PARTIAL_ZERO);
                 else:
                      Copy(madjtmp,scratchV2)
                 # end if #mhat at t
                 if i>0:
                      EvaluateRHSBwd(scratchV1, scratchV2, v, I, m, madjtmp, Iadjtmp, gcur, diffOp); MulC_I(scratchV1, dt) # k4 computed
                 else:
                      SetToIdentity(scratchG)
                      EvaluateRHSBwd(scratchV1, scratchV2, v, I, m, madjtmp, Iadjtmp, scratchG, diffOp); MulC_I(scratchV1, dt) # k4 computed

                 Add_I(RK4, scratchV1)

                 # final update
                 MulC_I(RK4,1.0/float(6.0))
                 Add_I(madj,RK4)

                 #FOR NEXT ITERATION:
                 # compute mhat_t, ihat_t at t to use in k1 computation. Note v_t, i_t and m_t are still stored from this iteration.

                 # if there is a measurement at this time point (for regression)
                 for k in range(len(msmtInds)):
                      if i>0:
                           if checkpointinds[msmtInds[k]] == i:
                                #Splat(scratchV1, ginvcur, IGradAtMsmts[k],False)
                                SplatSafe(scratchI, ginvcur, IGradAtMsmts[k])
                                Sub_I(Iadj,scratchI)
                                #Splat(Iadjtmp, gcur, Iadj,False)  # hat version of Iadj
                                SplatSafe(Iadjtmp, gcur, Iadj)  # hat version of Iadj
                           # end if
                      elif msmtInds[k]==-1: # if there is a measurement at time t=0 it won't be checkpointed but HARDCODED to have msmstInds == -1. Note this will be checked only for t==0
                           Sub_I(Iadj, IGradAtMsmts[k])
                      # end if
                 # end for

                 if i > 0: # last iteration skipped
                      ApplyH(madjtmp,madj,ginvcur, BACKGROUND_STRATEGY_PARTIAL_ZERO) # hat version of madj
                      # assign g, ginv appropriately for next iteration
                      g = gcur
                      ginv = ginvcur
                 # end if

             else:
                 # raise Exception('RK4 integration without all checkpoints not yet implemented')

                 # compute gcur and ginvcur by backward integrating geodesic
                 if scratchV6 is None:
                      print 'scratchV6 variable is initialized in IntegrateAdjoints'
                      scratchV6 = Field3D(I0.grid(),I0.memType())
                 if scratchV7 is None:
                      print 'scratchV7 variable is initialized in IntegrateAdjoints'
                      scratchV7 = Field3D(I0.grid(),I0.memType())
                 if scratchV8 is None:
                      print 'scratchV8 variable is initialized in IntegrateAdjoints'
                      scratchV8 = Field3D(I0.grid(),I0.memType())
                 if scratchV9 is None:
                      print 'scratchV9 variable is initialized in IntegrateAdjoints'
                      scratchV9 = Field3D(I0.grid(),I0.memType())

                 # initialize with previous
                 if prev_was_checkpoint == True:
                      gcur=scratchV8
                      ginvcur=scratchV9
                      Copy(gcur,g)
                      Copy(ginvcur,ginv)
                 # endif --previous was not checkpoint scratchV6 and scratch V8 should both have g and scratchV7 and scratch V9 should both have ginv. so no need to copy

                 # scratchG, scratchGinv and v are used as scratch variables in below call
                 m1 = scratchV3; bwdG = scratchV4; bwdGinv = scratchV5; # update references for ease of reading
                 IntegrateGeodesicBwdIteration(t,i,m0, checkpointstates[indx_of_last_tp][0], checkpointstates[indx_of_last_tp][1],m1, bwdG, bwdGinv,
                                                    gcur,ginvcur,
                                                    m1initialized,prev_was_checkpoint,
                                                    diffOp,
                                                    scratchGinv,
                                                    scratchV1,scratchV2,v,
                                                    Ninv=Ninv,integMethod=integMethod,
                                                    RK4=RK4,scratchG=scratchG)

                 # begin rk4 integration for adjoint
                 Copy(v,m); diffOp.applyInverseOperator(v);
                 EvaluateRHSBwd(scratchV1, scratchV2, v, I, m, madjtmp, Iadjtmp, g, diffOp); MulC_I(scratchV1, dt) # k1 computed
                 Copy(RK4,scratchV1)

                 # compute and store phi0t, phit0, v_t, i_t, m_t and i_hat at t+1/2 for computing k2 and k3
                 SubMulC(scratchG,g,gcur,0.5); #scratchG has w
                 UpdateInverse(scratchGinv, scratchV2, ginvcur, scratchG, Ninv)
                 Add_I(scratchG,gcur)

                 CoAd(m,scratchGinv,m0); Copy(v,m); diffOp.applyInverseOperator(v); ApplyH(I,I0,scratchGinv)
                 #Splat(Iadjtmp, scratchG, Iadj,False)
                 SplatSafe(Iadjtmp, scratchG, Iadj)

                 # for k2
                 # mhat_t at t+1/2 for k2
                 Add_MulC(scratchV2,madj,scratchV1,0.5) # mtilde at t+1/2
                 ApplyH(madjtmp,scratchV2,scratchGinv, BACKGROUND_STRATEGY_PARTIAL_ZERO); #mhat at t+1/2

                 EvaluateRHSBwd(scratchV1, scratchV2, v, I, m, madjtmp, Iadjtmp, scratchG, diffOp); MulC_I(scratchV1, dt) # k2 computed
                 Add_MulC_I(RK4, scratchV1, 2.0)

                 # for k3
                 # mhat_t at t+1/2 for k3
                 Add_MulC(scratchV2,madj,scratchV1,0.5) # mtilde at t+1/2
                 ApplyH(madjtmp,scratchV2,scratchGinv, BACKGROUND_STRATEGY_PARTIAL_ZERO); #mhat at t+1/2

                 EvaluateRHSBwd(scratchV1, scratchV2, v, I, m, madjtmp, Iadjtmp, scratchG, diffOp); MulC_I(scratchV1, dt) # k3 computed
                 Add_MulC_I(RK4, scratchV1, 2.0)

                 # compute v_t, i_t, m_t, i_hat at t for computing k4
                 CoAd(m,ginvcur,m0)
                 #Splat(Iadjtmp, gcur, Iadj,False)
                 SplatSafe(Iadjtmp, gcur, Iadj)
                 ApplyH(I,I0,ginvcur)
                 Copy(v,m); diffOp.applyInverseOperator(v);

                 # for k4
                 # mhat_t at t for k4
                 Add(scratchV2,madj,scratchV1) # mtilde at t
                 ApplyH(madjtmp,scratchV2,ginvcur, BACKGROUND_STRATEGY_PARTIAL_ZERO);

                 EvaluateRHSBwd(scratchV1, scratchV2, v, I, m, madjtmp, Iadjtmp, gcur, diffOp); MulC_I(scratchV1, dt) # k4 computed
                 Add_I(RK4, scratchV1)

                 # final update
                 MulC_I(RK4,1.0/float(6.0))
                 Add_I(madj,RK4)
                 #FOR NEXT ITERATION:
                 # compute mhat_t, ihat_t at t to use in k1 computation. Note v_t, i_t and m_t are still stored from this iteration.

                 # if there is a measurement at this time point (for regression)
                 for k in range(len(msmtInds)):
                    if i>0:
                         if checkpointinds[msmtInds[k]] == i:

                              #Splat(scratchV1, ginvcur, IGradAtMsmts[k],False)
                              SplatSafe(scratchI, ginvcur, IGradAtMsmts[k])
                              Sub_I(Iadj,scratchI)
                              #Splat(Iadjtmp, gcur, Iadj,False)  # hat version of Iadj
                              SplatSafe(Iadjtmp, gcur, Iadj)  # hat version of Iadj
                         # end if
                    elif msmtInds[k]==-1: # if there is a measurement at time t=0 it won't be checkpointed but HARDCODED to have msmstInds == -1. Note this will be checked only for t==0
                         Sub_I(Iadj,IGradAtMsmts[k])

                    # end if
                 # end for

                 ApplyH(madjtmp,madj,ginvcur, BACKGROUND_STRATEGY_PARTIAL_ZERO) # hat version of madj
                 # assign g, ginv appropriately for next iteration
                 # first assign references
                 g = scratchV6
                 ginv = scratchV7
                 # then copy memory
                 Copy(g,gcur)
                 Copy(ginv,ginvcur)
             # end if

        else:
            raise Exception('Unknown integration method: '+integMethod)
        #end if

        # common lines of code executed regardless of the integration scheme

    # end for
    return m
#end IntegrateAdjoints

def ParallelTransport(m1, n0, m0, nTimeSteps, diffOp, Ninv=10, integMethod='RK4',saveOutput=False, mtArray=None, ginvArray=None):
     '''
     Parallel translation of vector momentum, m0 along the geodesic
     denoted by initial condition, n0 from t=0 to t=1 using given time
     descritization, nTimeSteps.  Returns the parallely tranlated
     momentum in argument variable m1 at the end point of geodesic at
     t=1.  Also, returns list of norm of m(t), norm(n) and inner
     product between m and n at all timepoints along the geodesic.
     '''

     t = [x*1./nTimeSteps for x in range(nTimeSteps+1)]
     mGrid = n0.grid()
     mType = n0.memType()
     v = Field3D(mGrid, mType)
     #m = m1
     m = Field3D(mGrid, mType)
     w = Field3D(mGrid, mType)
     n = Field3D(mGrid, mType)

     g = Field3D(mGrid, mType)
     ginv = Field3D(mGrid, mType)

     k_n = Field3D(mGrid, mType)
     k_m = Field3D(mGrid, mType)
     scratchV1 = Field3D(mGrid, mType)
     scratchV2 = Field3D(mGrid, mType)


     scratchG = Field3D(mGrid, mType)
     scratchM = Field3D(mGrid, mType)
     RK4_m = Field3D(mGrid, mType)
     RK4_n = Field3D(mGrid, mType)

     SetToIdentity(g)
     SetToIdentity(ginv)
     Copy(n, n0)
     Copy(w,n)
     diffOp.applyInverseOperator(w)

     Copy(m, m0)
     Copy(v,m)
     diffOp.applyInverseOperator(v)
     #common.DebugHere()
     mtArray=[]
     ginvArray=[]
     if saveOutput:
          #common.DebugHere()
          mtArray.append(m.copy())
          ginvArray.append(ginv.copy())

     norm2_m=[Dot(m,v)]
     norm2_n=[Dot(n,w)]
     inner_m_n=[Dot(m,w)]

     for i in range(1,len(t),1):
          dt = t[i] - t[i-1]
          if integMethod == "EULER":
               raise Exception('Euler integration for parallel transport  not implemented: '+integMethod)
          elif integMethod == "RK4":
               EvaluateRHSFwd(k_n, w, g); MulC_I(k_n, dt) # k1_n computed
               Copy(RK4_n,k_n)

               EvaluateRHSParallelTransport(k_m, scratchV1, n, w, m, v, diffOp); MulC_I(k_m, dt) # k1_m computed
               Copy(RK4_m,k_m)

               # for k2_n
               Copy(scratchG, g); MulC_I(k_n,0.5); Add_I(scratchG,k_n); UpdateInverse(scratchV2, scratchV1, ginv, k_n, Ninv); CoAd(n,scratchV2, n0); Copy(w,n); diffOp.applyInverseOperator(w);
               EvaluateRHSFwd(k_n, w, scratchG); MulC_I(k_n, dt) # k2 computed
               Add_MulC_I(RK4_n, k_n, 2.0)

               # for k2_m
               Copy(scratchM, m); MulC_I(k_m,0.5); Add_I(scratchM,k_m); Copy(v,scratchM); diffOp.applyInverseOperator(v);
               EvaluateRHSParallelTransport(k_m, scratchV1, n, w, scratchM, v, diffOp); MulC_I(k_m, dt) # k2_m computed
               Add_MulC_I(RK4_m, k_m, 2.0)

               # for k3_n
               Copy(scratchG, g); MulC_I(k_n,0.5); Add_I(scratchG,k_n); UpdateInverse(scratchV2, scratchV1, ginv, k_n, Ninv); CoAd(n,scratchV2, n0);Copy(w,n); diffOp.applyInverseOperator(w);
               EvaluateRHSFwd(k_n, w, scratchG); MulC_I(k_n, dt) # k3_n computed
               Add_MulC_I(RK4_n, k_n, 2.0)

               # for k3_m
               Copy(scratchM, m); MulC_I(k_m,0.5); Add_I(scratchM,k_m); Copy(v,scratchM); diffOp.applyInverseOperator(v);
               EvaluateRHSParallelTransport(k_m, scratchV1, n, w, scratchM, v, diffOp); MulC_I(k_m, dt) # k3_m computed
               Add_MulC_I(RK4_m, k_m, 2.0)

               # for k4_n
               Copy(scratchG, g);  Add_I(scratchG,k_n); UpdateInverse(scratchV2, scratchV1, ginv, k_n, Ninv); CoAd(n,scratchV2, n0);Copy(w,n); diffOp.applyInverseOperator(w);
               EvaluateRHSFwd(k_n, w, scratchG); MulC_I(k_n, dt) # k4_n computed
               Add_I(RK4_n, k_n)

               # for k4_m
               Copy(scratchM, m);  Add_I(scratchM,k_m); Copy(v,scratchM); diffOp.applyInverseOperator(v);
               EvaluateRHSParallelTransport(k_m, scratchV1, n, w, scratchM, v, diffOp); MulC_I(k_m, dt) # k4_m computed
               Add_I(RK4_m, k_m)

               # final update
               MulC_I(RK4_n,1.0/float(6.0))
               Add_I(g,RK4_n)
               UpdateInverse(scratchV2, scratchV1, ginv, RK4_n, Ninv)
               Copy(ginv,scratchV2)

               Add_MulC_I(m, RK4_m, 1.0/float(6.0))
          else:
               raise Exception('Unknown integration method: '+integMethod)
          # end if

          # Coadjoint action gives us momentum at time t
          CoAd(n, ginv, n0)
          Copy(w,n)
          # Smooth to get velocity at time t
          diffOp.applyInverseOperator(w)
          Copy(v,m)
          diffOp.applyInverseOperator(v)
          if saveOutput:
               mtArray.append(m.copy())
               ginvArray.append(ginv.copy())
          # save norms and inner product
          norm2_m.append(Dot(m,v))
          norm2_n.append(Dot(n,w))
          inner_m_n.append(Dot(m,w))
     # end for

     Copy(m1,m)
     return norm2_m, norm2_n, inner_m_n, mtArray, ginvArray

def EvaluateRHSParallelTransport(out_m, scratchV, n, w, m, v, diffOp):
     '''
     Evaluate RHS for parallel transport
     '''
     AdInf(out_m, w, v)
     diffOp.applyOperator(out_m)

     CoAdInf(scratchV, w, m)
     Sub_I(out_m, scratchV)

     CoAdInf(scratchV, v, n)
     Sub_I(out_m, scratchV)

     MulC_I(out_m, 0.5)
     return out_m
#end EvaluateRHSFwd
def MyGridPlot(vf, sliceIdx=None, dim='z', every=1, isVF=True,
               color='g', plotBase=True, colorbase='#A0A0FF'):
    sliceArr = common.ExtractSliceArrVF(vf, sliceIdx, dim)
    sz = sliceArr.shape
    hID = np.mgrid[0:sz[0], 0:sz[1]]
    d1 = np.squeeze(hID[1, ::every, ::every])
    d2 = np.squeeze(hID[0, ::every, ::every])
    sliceArr = sliceArr[::every, ::every, :]
    if plotBase:
        plt.plot(d1, d2, colorbase)
        plt.hold(True)
        plt.plot(d1.T, d2.T, colorbase)
    if not isVF:
        d1 = np.zeros(d1.shape)
        d2 = np.zeros(d2.shape)

    if dim=='z':
        plt.plot(d1+np.squeeze(sliceArr[:,:,1]),
                 d2+np.squeeze(sliceArr[:,:,0]), color)
        plt.hold(True)
        plt.plot((d1+np.squeeze(sliceArr[:,:,1])).T,
                 (d2+np.squeeze(sliceArr[:,:,0])).T, color)
        plt.hold(False)
    elif dim=='x':
        plt.plot(d1+np.squeeze(sliceArr[:,:,2]),
                 d2+np.squeeze(sliceArr[:,:,1]), color)
        plt.hold(True)
        plt.plot((d1+np.squeeze(sliceArr[:,:,2])).T,
                 (d2+np.squeeze(sliceArr[:,:,1])).T, color)
        plt.hold(False)
    elif dim=='y':
        plt.plot(d1+np.squeeze(sliceArr[:,:,2]),
                 d2+np.squeeze(sliceArr[:,:,0]), color)
        plt.hold(True)
        plt.plot((d1+np.squeeze(sliceArr[:,:,2])).T,
                 (d2+np.squeeze(sliceArr[:,:,0])).T, color)
        plt.hold(False)

    # change axes to match image axes
    if not plt.gca().yaxis_inverted():
        plt.gca().invert_yaxis()
        # force redraw
        plt.draw()


def MyQuiver(vf, sliceIdx=None, dim='z',every=1,thresh=None,scaleArrows=None,arrowCol='r',lineWidth=None, width=None):
    sliceArr = common.ExtractSliceArrVF(vf, sliceIdx, dim)
    if dim=='z':
         vy = np.squeeze(sliceArr[:,:,0])
         vx = np.squeeze(sliceArr[:,:,1])
    elif dim=='x':
         vy = np.squeeze(sliceArr[:,:,1])
         vx = np.squeeze(sliceArr[:,:,2])
    elif dim=='y':
         vy = np.squeeze(sliceArr[:,:,0])
         vx = np.squeeze(sliceArr[:,:,2])
    vxshow = np.zeros(np.shape(vx))
    vyshow = np.zeros(np.shape(vy))

    vxshow[::every,::every] = vx[::every,::every]
    vyshow[::every,::every] = vy[::every,::every]
    valindex = np.zeros(np.shape(vx),dtype=bool)
    valindex[::every,::every] = True

    if thresh is not None:
        valindex[(vx**2+vy**2)<thresh] = False


    #gridX, gridY = np.meshgrid(range(np.shape(vx)[0]),range(np.shape(vx)[1]-1,-1,-1))
    gridX, gridY = np.meshgrid(range(np.shape(vx)[0]),range(np.shape(vx)[1]))
    quiverhandle = plt.quiver(gridX[valindex],gridY[valindex],vx[valindex],vy[valindex],scale=scaleArrows,color=arrowCol,linewidth=lineWidth,width=width,zorder=4)

    # change axes to match image axes
    plt.gca().invert_yaxis()

    # force redraw
    plt.draw()

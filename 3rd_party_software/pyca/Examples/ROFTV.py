from PyCA.Core import *

import PyCA.Common as common
import PyCA.Display as display

import numpy as np
import matplotlib.pyplot as plt

#
# Allocate all necessary data
#
def InitializeData(grid, mType):
    global scratchI
    scratchI = Image3D(grid, mType)
    global dxf
    dxf = Image3D(grid, mType)
    global dyf
    dyf = Image3D(grid, mType)
    global dx
    dx = Image3D(grid, mType)
    global dy
    dy = Image3D(grid, mType)
    global dy_xf
    dy_xf = Image3D(grid, mType)
    global dx_yf
    dx_yf = Image3D(grid, mType)
    global gxf
    gxf = Image3D(grid, mType)
    global gyf
    gyf = Image3D(grid, mType)

    av_x_arr = np.zeros([3,1])
    av_x_arr[1,0] = 0.5
    av_x_arr[2,0] = 0.5
    global av_x
    av_x = Image3D(mType)
    av_x.fromlist(np.atleast_3d(av_x_arr).tolist())

    av_y_arr = np.zeros([1,3])
    av_y_arr[0,1] = 0.5
    av_y_arr[0,2] = 0.5
    global av_y
    av_y = Image3D(mType)
    av_y.fromlist(np.atleast_3d(av_y_arr).tolist())

#
# Set grid (for multiscale)
#
def SetGrid(grid):
    global scratchI
    scratchI.setGrid(grid)
    global dxf
    dxf.setGrid(grid)
    global dyf
    dyf.setGrid(grid)
    global dx
    dx.setGrid(grid)
    global dy
    dy.setGrid(grid)
    global dy_xf
    dy_xf.setGrid(grid)
    global dx_yf
    dx_yf.setGrid(grid)
    global gxf
    gxf.setGrid(grid)
    global gyf
    gyf.setGrid(grid)

def calcEnergy(I, Data, DataFidC, TVC, TVPow=1.0, Mask=None):
    energy = []
    Sub(scratchI, I, Data)
    if Mask != None:
        Mul_I(scratchI, Mask)
    imEn = (DataFidC/2.0)*Sum2(scratchI)
    energy.append(imEn)
    GradForMag(scratchI, I)
    #GradientMag(scratchI, I)
    if TVPow != 1.0:
        PowC_I(scratchI, TVPow)
    tvEn = TVC*Sum(scratchI)
    energy.append(tvEn)
    energy.append(imEn+tvEn)
    return energy

def calcTVGradient(I, TVC, Beta, gTV, TVPow=1.0):
    # Compute TV gradient
    FiniteDiff(dxf, I, DIM_X, DIFF_FORWARD)
    FiniteDiff(dyf, I, DIM_Y, DIFF_FORWARD)
    FiniteDiff(dx, I, DIM_X, DIFF_CENTRAL)
    FiniteDiff(dy, I, DIM_Y, DIFF_CENTRAL)
    
    Convolve(dy_xf, dy, av_x)
    Convolve(dx_yf, dx, av_y)
    
    # dy_xf = dy_xf*dy_xf + dxf*dxf
    Mul(gxf, dy_xf, dy_xf)
    Add_Mul_I(gxf, dxf, dxf)

    # dx_yf = dx_yf*dx_yf + dyf*dyf
    Mul(gyf, dx_yf, dx_yf)
    Add_Mul_I(gyf, dyf, dyf)

    # add regularizer
    AddC_I(gxf, Beta)
    AddC_I(gyf, Beta)

    if TVPow == 1.0:
        Sqrt_I(gxf)
        Sqrt_I(gyf)
    else:
        PowC_I(gxf, 1.0-(TVPow/2.0))
        PowC_I(gyf, 1.0-(TVPow/2.0))

    Div_I(dxf, gxf)
    Div_I(dyf, gyf)

    FiniteDiff(dx, dxf, DIM_X, DIFF_BACKWARD)
    FiniteDiff(dy, dyf, DIM_Y, DIFF_BACKWARD)

    AddMulC(gTV, dx, dy, -(TVPow*TVC))

# overwrites LDefSum
def plotResults(I, Data, fig,cmap='gray',rng=[0,1]):
    plt.figure(fig)
    plt.subplot(1,3,1)
    display.DispImage(Data, 'Orig', cmap=cmap, \
                               newFig=False, rng=rng, t=False)
    plt.subplot(1,3,2)
    display.DispImage(I, 'Denoised', cmap=cmap, \
                               newFig=False, rng=rng, t=False)
    Sub(scratchI, I, Data)
    plt.subplot(1,3,3)
    display.DispImage(scratchI, 'diff', cmap=cmap, \
                               newFig=False, rng=None, t=False)
    plt.draw()
    plt.show()

def plotEnergy(en, fig):
    plt.figure(fig)
    plt.plot(en[0][:],'r')
    plt.hold(True)
    plt.plot(en[1][:],'g')
    plt.plot(en[2][:],'b')
    plt.hold(False)
    plt.draw()
    plt.show()


def RunROFTV(Data, \
                 I0=None, \
                 DataFidC=25.0, \
                 TVC=1.0, \
                 Beta=1e-5, \
                 TVPow=1.0, \
                 stepI=None, \
                 nIters=5000, \
                 dispEvery=1000, \
                 disp=True, \
                 Mask=None):

    if stepI == None:
        stepI = 1.0/(4.0*DataFidC)

    grid = Data.grid().copy()
    mType = Data.memType()

    #
    # Initialize data
    #
    
    # TV-denoised image
    if I0 is None:
        I = Data.copy()
    else:
        I = I0.copy()

    diff = Image3D(grid, mType)
    speed = Image3D(grid, mType)
    gI = Image3D(grid, mType)

    # other vars
    InitializeData(grid, mType)
    
    # Initialize other data
    energy = [[] for _ in range(3)]

    if disp:
        EnergyFig = plt.figure('ROF Energy');
        plt.clf();
        ResultFig = plt.figure('ROF Results');
        plt.clf();

    for k in range(nIters+1):

        print 'iteration %d...'%k

        #
        # Display images
        #
        if disp and k%dispEvery == 0:
            plotResults(I, Data, ResultFig.number)

        #
        # Compute energy
        #

        en = calcEnergy(I, Data, DataFidC, TVC, TVPow, Mask)
        # primal energy

        energy[0].append(en[0])
        energy[1].append(en[1])
        energy[2].append(en[2])

        if disp and k%dispEvery == 0:
            plotEnergy(energy, EnergyFig.number)

        # just compute energy on final iteration
        if k >= nIters:
            break

        # gradient stored in gI
        calcTVGradient(I, TVC, Beta, gI, TVPow)

        # compute the image difference term
        Sub(diff, I, Data)
        MulC_I(diff, DataFidC)
        if Mask != None:
            Mul_I(diff, Mask)

        Add_I(gI, diff)

        if False:
            MulC(speed, gI, -1.0)
            UpwindGradMag(gI, I, speed)
            MulC_I(gI, -1.0)

        # apply update
        Add_MulC_I(I, gI, -stepI)

    if disp:
        plotResults(I, Data, ResultFig.number)
        plotEnergy(energy, EnergyFig.number)

    return (I, energy)

def RunTest():

    # number of iterations
    nIters = 2000
    #nIters = 0

    disp = True
    dispEvery = 1000

    if GetNumberOfCUDADevices() > 0:
        mType = MEM_DEVICE
    else:
        print "No CUDA devices found, running on CPU"
        mType = MEM_HOST

    # data fidelity modifier
    DataFidC = 20.0
    # TV modifier
    TVC=1.0

    TVPow=1.0

    UseMask = True

    # regularization term to avoid zero denominator
    Beta = 1e-5

    stepI = 0.001

    imagedir='./Images/'

    #
    # Run lena images
    #
    
    Data = common.LoadPNGImage(imagedir + 'lena_orig.png', mType)

    imSz = Data.size()
    sz = imSz.tolist()[0:2]

    if True:
        I0 = Data.copy()
    else:
        I0 = common.RandImage(nSig=1.0, gSig=5.0, mType=mType)

    Mask = None
    if UseMask:
        bdr = 10
        MaskArr = np.zeros(sz)
        MaskArr[bdr:-bdr,bdr:-bdr] = 1.0
        Mask = common.ImFromNPArr(MaskArr, mType)

    (I, energy) = \
        RunROFTV(Data=Data, \
                     I0 = I0, \
                     DataFidC=DataFidC, \
                     TVC=TVC, \
                     TVPow=TVPow, \
                     stepI=stepI, \
                     Beta=Beta, \
                     nIters=nIters, \
                     dispEvery=dispEvery, \
                     disp=disp, \
                     Mask=Mask)

    print 'final energy: {ttl:n} = {im:n} + {tv:n}'\
        .format(ttl=energy[2][-1], \
                    im=energy[0][-1], \
                    tv=energy[1][-1])

if __name__ == '__main__':

    plt.close('all')

    RunTest()

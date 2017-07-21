#
# General additions to unittest for testing PyCA
#
import sys
import unittest

from PyCA.Core import DIM_X,DIM_Y,DIM_Z
from PyCA.Core import DIFF_FORWARD, DIFF_BACKWARD, DIFF_CENTRAL
from PyCA.Core import BC_APPROX, BC_WRAP, BC_CLAMP

#
# decorator func to add setup/teardown code to member function
#
def AddSetUp(setUp, tearDown=None):
    def AddSetUpFunc(f):
        def AddSetUpWrappedFunc(*args,**kwargs):
            selfOb = args[0]
            if setUp != None:
                setUp(selfOb)
            f(*args,**kwargs)
            if tearDown != None:
                tearDown(selfOb)
        return AddSetUpWrappedFunc
    return AddSetUpFunc

#
# decorator func to skip test if disp=False (assumes disp is first arg
# for member function, ie second overall argument after 'self')
#
def SkipIfNotDisp(msg=""):
    def SkipIfNotDispFunc(f):
        def SkipIfNotDispWrappedFunc(*args,**kwargs):
            if (len(args) > 1 and args[1] == True) or \
                    'disp' in kwargs and kwargs['disp'] == True:
                f(*args,**kwargs)
            else:
                sys.stdout.write(msg+'...')
                sys.stdout.flush()
        return SkipIfNotDispWrappedFunc
    return SkipIfNotDispFunc

# stringification of constants
DIMNAMES={DIM_X : 'DIM_X', DIM_Y : 'DIM_Y', DIM_Z : 'DIM_Z'}
DIFFTNAMES={DIFF_FORWARD : 'DIFF_FORWARD', DIFF_BACKWARD : 'DIFF_BACKWARD', DIFF_CENTRAL : 'DIFF_CENTRAL'}
BCNAMES={BC_APPROX : 'BC_APPROX', BC_WRAP : 'BC_WRAP', BC_CLAMP : 'BC_CLAMP'}

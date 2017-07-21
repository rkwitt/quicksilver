import PyCA.Core as core

import inspect

# Debugging, just add DebugHere() to start debugger at that line in
# ipython
try:
    from IPython.core.debugger import Tracer
except ImportError:
    ipy_debug_available = False
else:
    ipy_debug_available = True

def DebugHere():
    if ipy_debug_available:
        t = Tracer()
        t()
    else:
        print 'ipython not available for debugging'
        #raise Exception('ipython not available for debugging')
        
def GetOuterFrame():
    curframe = inspect.currentframe()
    upframe = inspect.getouterframes(curframe)[2][0]
    info = inspect.getframeinfo(upframe)
    return info[:3]


def CheckCUDAError(msg=''):
    (file, line, func) = GetOuterFrame()
    funcmsg = func + ": " + msg
    core.CheckCUDAError(file, line, funcmsg)

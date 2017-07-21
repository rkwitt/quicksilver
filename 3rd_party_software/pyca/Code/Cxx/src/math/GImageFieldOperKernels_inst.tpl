def IOperBG = [
BACKGROUND_STRATEGY_CLAMP,
BACKGROUND_STRATEGY_WRAP,
BACKGROUND_STRATEGY_ZERO,
BACKGROUND_STRATEGY_PARTIAL_ZERO]

def Interp = [
INTERP_NN,
INTERP_LINEAR,
INTERP_CUBIC]

def Fwd = [
true, 
false]

def SqrRoot = [
true, 
false]

def DiffT = [
DIFF_FORWARD, 
DIFF_BACKWARD, 
DIFF_CENTRAL]

template void PyCA::
ApplyH<${IOperBG}, ${Interp}>(float*, const float*, const float*, const float*, const float*, int, int, int, StreamT);

template void PyCA::
ApplyV<${Fwd}, ${IOperBG}, ${Interp}>(float*, const float*, const float*, const float*, const float*, const float&, int, int, int, float, float, float, StreamT, bool);

template void PyCA::
ComposeTranslation<${IOperBG}, ${Interp}>(float*, const float*, const Vec3Df&, const Vec3Di&, StreamT, bool);

template void PyCA::
g_gradient<${DiffT}>(float*, float*, float*, const float*, int, int, int, float, float, float, BoundaryCondT, StreamT);

template void PyCA::
g_gradientMask<${DiffT}>(float*, float*, float*, const float*, const float*, int, int, int, float, float, float, BoundaryCondT, StreamT);

template void PyCA::
g_gradientMag<float, ${DiffT}>(float*, const float*, int, int, int, float, float, float, BoundaryCondT, StreamT);

template void PyCA::
g_divergence<${DiffT}>(float*, const float*, const float*, const float*, int, int, int, float, float, float, BoundaryCondT, StreamT);

template void PyCA::
JacDetV<${Fwd}>(float*, const float*, const float*, const float*, const float*, const float*, const float*, const float*, const float*, const float*, size_t, StreamT);


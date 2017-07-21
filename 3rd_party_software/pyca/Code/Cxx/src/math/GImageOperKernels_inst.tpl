def IOperBG = [
BACKGROUND_STRATEGY_CLAMP,
BACKGROUND_STRATEGY_WRAP,
BACKGROUND_STRATEGY_ZERO,
BACKGROUND_STRATEGY_PARTIAL_ZERO]

def Interp = [
INTERP_NN,
INTERP_LINEAR,
INTERP_CUBIC]

def BoolArg=[true, false]

// template instantiations
template void Resample<${IOperBG}, ${Interp}, ${BoolArg}>(float*, const Vec3Di&, const float*, const Vec3Di &, StreamT);

template void ResampleWorld<${IOperBG}, ${Interp}>(float*, const Vec3Di&, const Vec3Df&, const Vec3Df&, const float*, const Vec3Di&, const Vec3Df&, const Vec3Df&, StreamT);

template void SplatWorld<${IOperBG}>(float*, const Vec3Di&, const Vec3Df&, const Vec3Df&, const float*, const Vec3Di&, const Vec3Df&, const Vec3Df&, StreamT);

template void SplatWorld<${IOperBG}>(float*, const Vec3Di&, const Vec3Df&, const Vec3Df&, const float*, const Vec3Di&, const Vec3Df&, const Vec3Df&, float*, StreamT);

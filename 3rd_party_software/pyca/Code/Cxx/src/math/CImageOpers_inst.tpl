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
template void CImageOpers::Resample<${IOperBG}, ${Interp}, ${BoolArg}>(Image3D&, const Image3D&, StreamT);

// template instantiations
template void CImageOpers::ResampleWorld<${IOperBG}, ${Interp}>(Image3D&, const Image3D&, StreamT);

template void CImageOpers::SplatWorld<${IOperBG}>(Image3D&, const Image3D&, StreamT);
template void CImageOpers::SplatWorld<${IOperBG}>(Image3D&, const Image3D&, Image3D&, StreamT);

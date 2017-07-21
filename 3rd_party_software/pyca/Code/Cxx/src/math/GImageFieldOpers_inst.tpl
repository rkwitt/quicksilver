def IOperBG = [
BACKGROUND_STRATEGY_CLAMP,
BACKGROUND_STRATEGY_WRAP,
BACKGROUND_STRATEGY_ZERO,
BACKGROUND_STRATEGY_PARTIAL_ZERO]

def Interp = [
INTERP_NN,
INTERP_LINEAR,
INTERP_CUBIC]

template void GImageFieldOpers::
ApplyH<${IOperBG}, ${Interp}>(Image3D&, const Image3D&, const Field3D&, StreamT);

template void GImageFieldOpers::
ApplyV<${IOperBG}, ${Interp}>(Image3D&, const Image3D&, const Field3D&, const float&, StreamT, bool);

template void GImageFieldOpers::
ApplyVInv<${IOperBG}, ${Interp}>(Image3D&, const Image3D&, const Field3D&, const float&, StreamT, bool);

template void GImageFieldOpers::
ComposeTranslation<${IOperBG}, ${Interp}>(Image3D&, const Image3D&, const Vec3Df&, StreamT, bool);

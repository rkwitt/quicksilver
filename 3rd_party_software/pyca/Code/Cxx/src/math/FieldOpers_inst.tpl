def VOperBG = [
BACKGROUND_STRATEGY_CLAMP,
BACKGROUND_STRATEGY_WRAP,
BACKGROUND_STRATEGY_ZERO,
BACKGROUND_STRATEGY_PARTIAL_ZERO]

def HOperBG = [
BACKGROUND_STRATEGY_CLAMP,
BACKGROUND_STRATEGY_ID,
BACKGROUND_STRATEGY_PARTIAL_ID]

def AllBG = [
BACKGROUND_STRATEGY_CLAMP,
BACKGROUND_STRATEGY_WRAP,
BACKGROUND_STRATEGY_ZERO,
BACKGROUND_STRATEGY_PARTIAL_ZERO,
BACKGROUND_STRATEGY_ID,
BACKGROUND_STRATEGY_PARTIAL_ID]

def IDBG = [
BACKGROUND_STRATEGY_ID,
BACKGROUND_STRATEGY_PARTIAL_ID]

def BoolArg = [true, false]

#ifdef CUDA_ENABLED
def ExecMode = [EXEC_CPU, EXEC_GPU]
#else
def ExecMode = [EXEC_CPU]
#endif

template void FieldOpers<${ExecMode}>::
ComposeHH<${HOperBG}>(Field3D&, const Field3D&, const Field3D&, StreamT);

template void FieldOpers<${ExecMode}>::
ComposeVH<${VOperBG}>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);

template void FieldOpers<${ExecMode}>::
ComposeVInvH<${VOperBG}>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);

template void FieldOpers<${ExecMode}>::
ComposeHV<${HOperBG}>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);

template void FieldOpers<${ExecMode}>::
ComposeHVInv<${HOperBG}>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT,bool);

template void FieldOpers<${ExecMode}>::
ComposeTranslation<${AllBG}>(Field3D&, const Field3D&, const Vec3Df&, StreamT,bool);

template void FieldOpers<${ExecMode}>::
ApplyH<${AllBG}>(Field3D&, const Field3D&, const Field3D&, StreamT);

template void FieldOpers<${ExecMode}>::
ApplyV<${AllBG}>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT, bool);

template void FieldOpers<${ExecMode}>::
ApplyVInv<${AllBG}>(Field3D&, const Field3D&, const Field3D&, const float&, StreamT, bool);

template void FieldOpers<${ExecMode}>::
Resample<${AllBG}, ${BoolArg}>(Field3D&, const Field3D&, StreamT);

template void FieldOpers<${ExecMode}>::
FixedPointInverse<${IDBG}>(Field3D&, const Field3D&, unsigned int, StreamT, bool);

template void FieldOpers<${ExecMode}>::
Ad<${VOperBG}>(Field3D&, const Field3D&, const Field3D&, StreamT, bool);

template void FieldOpers<${ExecMode}>::
CoAd<${VOperBG}>(Field3D&, const Field3D&, const Field3D&, StreamT, bool);

template class FieldOpers<${ExecMode}>;

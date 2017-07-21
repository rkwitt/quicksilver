def IOperBG = [
BACKGROUND_STRATEGY_CLAMP,
BACKGROUND_STRATEGY_WRAP,
BACKGROUND_STRATEGY_ZERO,
BACKGROUND_STRATEGY_PARTIAL_ZERO]

def Interp = [
INTERP_NN,
INTERP_LINEAR,
INTERP_CUBIC]

#ifdef CUDA_ENABLED
def ExecMode = [EXEC_CPU, EXEC_GPU]
#else
def ExecMode = [EXEC_CPU]
#endif

template void ImageFieldOpers<${ExecMode}>::
ApplyH<${IOperBG}, ${Interp}>(Image3D&, const Image3D&, const Field3D&, StreamT);

template void ImageFieldOpers<${ExecMode}>::
ApplyV<${IOperBG}, ${Interp}>(Image3D&, const Image3D&, const Field3D&, const float&, StreamT, bool);

template void ImageFieldOpers<${ExecMode}>::
ApplyVInv<${IOperBG}, ${Interp}>(Image3D&, const Image3D&, const Field3D&, const float&, StreamT, bool);

template void ImageFieldOpers<${ExecMode}>::
ComposeTranslation<${IOperBG}, ${Interp}>(Image3D&, const Image3D&, const Vec3Df&, StreamT, bool);

template class ImageFieldOpers<${ExecMode}>;

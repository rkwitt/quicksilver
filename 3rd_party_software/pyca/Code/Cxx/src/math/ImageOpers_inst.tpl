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

def BoolArg = [true, false]

template
void ImageOpers<${ExecMode}>::Resample<${IOperBG}, ${Interp}, ${BoolArg}>(Image3D&, const Image3D&, StreamT);

// function template instantiation

template
void ImageOpers<${ExecMode}>::ResampleWorld<${IOperBG}, ${Interp}>(Image3D&, const Image3D&, StreamT);

template
void ImageOpers<${ExecMode}>::SplatWorld<${IOperBG}>(Image3D&, const Image3D&, StreamT);
template
void ImageOpers<${ExecMode}>::SplatWorld<${IOperBG}>(Image3D&, const Image3D&, Image3D&, StreamT);

template class ImageOpers<${ExecMode}>;


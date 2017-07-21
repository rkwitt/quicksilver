#ifdef CUDA_ENABLED
def ExecMode = [EXEC_CPU, EXEC_GPU]
#else
def ExecMode = [EXEC_CPU]
#endif

template class VFieldOpers<${ExecMode}>;


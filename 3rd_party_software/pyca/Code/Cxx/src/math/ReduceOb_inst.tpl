
#ifdef CUDA_ENABLED
def ExecModeAll = [EXEC_CPU, EXEC_GPU, EXEC_GPU_PARAM]
#else
def ExecModeAll = [EXEC_CPU]
#endif

template class ReduceOb<${ExecModeAll}>;
   

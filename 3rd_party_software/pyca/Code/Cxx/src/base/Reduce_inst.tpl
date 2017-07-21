/* ================================================================
 *
 * PyCA Project
 *
 * Copyright (c) J. Samuel Preston, Linh K. Ha, Sarang C. Joshi. All
 * rights reserved.  See Copyright.txt or for details.
 *
 * This software is distributed WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the above copyright notice for more information.
 *
 * ================================================================ */

#ifdef CUDA_ENABLED
def ExecModeAll = [EXEC_CPU, EXEC_GPU, EXEC_GPU_PARAM]
#else
def ExecModeAll = [EXEC_CPU]
#endif

def ScalarType = [float, int, uint]

template
void Reduce<${ExecModeAll}>::Max(${ScalarType}& h_o, const ${ScalarType}* h_i, size_t n, bool update, StreamT stream);

template
void Reduce<${ExecModeAll}>::Min(${ScalarType}& h_o, const ${ScalarType}* h_i, size_t n, bool update, StreamT stream);

template
void Reduce<${ExecModeAll}>::Sum(${ScalarType}& h_o, const ${ScalarType}* h_i, size_t n, bool update, StreamT stream);

template
void Reduce<${ExecModeAll}>::LInf(${ScalarType}& h_o, const ${ScalarType}* h_i, size_t n, bool update, StreamT stream);

template
void Reduce<${ExecModeAll}>::L1(${ScalarType}& h_o, const ${ScalarType}* h_i, size_t n, bool update, StreamT stream);

template
void Reduce<${ExecModeAll}>::Sum2(${ScalarType}& h_o, const ${ScalarType}* h_i, size_t n, bool update, StreamT stream);

template
void Reduce<${ExecModeAll}>::MaxMin(Vec2D<${ScalarType}>& h_o, const ${ScalarType}* h_i, size_t n, bool update, StreamT stream);

// have to repeat three times so we don't get all permutations of two ScalarType
template
void Reduce<${ExecModeAll}>::Dot(${ScalarType}& h_o, const ${ScalarType}* h_i, const ${ScalarType}* h_i1, size_t n, bool update, StreamT stream);

template class Reduce<${ExecModeAll}>;


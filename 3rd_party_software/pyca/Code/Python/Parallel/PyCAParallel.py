"""Parallelization helper functions for PyCA"""

import PyCA.Core as ca

import numpy as np

try:
   from mpi4py import MPI
   hasMPI = True
except:
   hasMPI = False

def Reduce(A, hA, op=None):
    """Reduce PyCA Image3D or Field3D over MPI

    A can live anywhere but hA needs to be of mType MEM_HOST

    """
    if not hasMPI:
      raise Exception("mpi4py required for Reduce operations: not found")

    if op is None:
      op = MPI.SUM

    if A.memType() == ca.MEM_HOST:
        # can do this in place without using hA
        if isinstance(A, ca.Image3D):
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, A.asnp(), op=op)
        elif isinstance(A, ca.Field3D):
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, A.x_asnp(), op=op)
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, A.y_asnp(), op=op)
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, A.z_asnp(), op=op)
        else:
            raise Exception('Can only reduce Image3D and Field3D')
    else:
        # will need some uploading and downloading
        assert(hA.memType() == ca.MEM_HOST)  # make sure we can actually use hA
        ca.Copy(hA, A)  # download

        if isinstance(A, ca.Image3D):
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, hA.asnp(), op=op)
        elif isinstance(A, ca.Field3D):
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, hA.x_asnp(), op=op)
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, hA.y_asnp(), op=op)
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, hA.z_asnp(), op=op)
        else:
            raise Exception('Can only reduce Image3D and Field3D')

        ca.Copy(A, hA)  # upload


def ReduceFloat(f, op=None):
    """Reduce a single float value over MPI"""
    if not hasMPI:
      raise Exception("mpi4py required for Reduce operations: not found")

    if op is None:
      op = MPI.SUM

    fa = np.array([f])  # can only reduce over numpy arrays
    MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE,
                             fa,
                             op=MPI.SUM)
    return fa[0]

"""Distributed computation convenience functions for use with PyCA"""

import Config

import PyCA.Core as ca

from mpi4py import MPI
import os
import sys

# This is meant to be a general computation configuration spec
# Nothing in here should be specific to the algorithm being run
# The idea is that subsections conforming to this spec could be held in a
# central spot such as ~/.pyca/ and be included on the fly by command-line
# overrides such as
# ./Matching.py study.yaml computeConfig="~/.pyca/mpicluster.yaml"
# or by simply changing a single line within study.yaml
ComputeConfigSpec = {
    'useCUDA':
    Config.Param(default=False,
                 comment="Use GPU if available"),
    'gpuID':
    Config.Param(default=0,
                 comment="Used when useCUDA is true. This is when a node has multiple GPUs and \
                 user needs to specify which GPU to use. Can be any integer in [0,#GPUs-1]. Typically not useful with MPI processes."),
    'useMPI':
    Config.Param(default=False,
                 comment="Use MPI for multiple nodes. Is overriden and set to True if you spawn the mpi processes from outside."),
    'interactive':
    Config.Param(default=False,
                 comment="Run in interactive mode (i.e. do plotting etc)"),
    'numProcesses':
    Config.Param(default=1,
                 comment="Number of MPI processes to spawn. Not used if started process already is an MPI process."),
    'hostFile':
    Config.Param(default="",
                 comment="Path to the file that list of hosts to use as MPI worker nodes. If left blank, localhost is used. Not used if started process already is an MPI process.")}

sys_excepthook = sys.excepthook
def mpi_excepthook(type, value, traceback): 
    sys_excepthook(type, value, traceback) 
    MPI.COMM_WORLD.Abort(1) 

def Spawn(nproc=1, hostList=None):
    """Spawn mpiexec if not already done"""
    # detect if we've started this process in MPI yet or do that now
    if not 'OMPI_COMM_WORLD_LOCAL_RANK' in os.environ:
        if hostList is None or len(hostList) == 0:
            hostList = ['localhost']

        print "Spawning " + str(nproc) + " MPI processes..."
        # TODO: use mpi4py for spawning instead of mpiexec
        #comm = MPI.COMM_SELF.Spawn(sys.executable,
                                    #args=sys.argv, maxprocs=nproc)
        CMD = "mpiexec -x PYTHONPATH -n " + str(nproc) + " --host \"" +\
              ','.join(hostList) + "\" " + ' '.join(sys.argv)
        print(CMD)
        os.system(CMD)
        sys.exit(0)            

def Reduce(A, hA, op=MPI.SUM):
    """Reduce PyCA Image3D or Field3D over MPI
    A can live anywhere but hA needs to be of mType MEM_HOST
    """
    if A.memType() == ca.MEM_HOST:
        # can do this in place without using hA
        if isinstance(A, ca.Field3D):
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, A.asnp()[0], op=op)
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, A.asnp()[1], op=op)
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, A.asnp()[2], op=op)
        else:
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, A.asnp(), op=op)
    else:
        # will need some uploading and downloading
        assert(hA.memType() == ca.MEM_HOST)  # make sure we can actually use hA
        ca.Copy(hA, A)  # download
        if isinstance(A, ca.Field3D):
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, hA.asnp()[0], op=op)
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, hA.asnp()[1], op=op)
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, hA.asnp()[2], op=op)
        else:
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, hA.asnp(), op=op)
        ca.Copy(A, hA)  # upload


def GetMPIInfo():
    """Just get the size, rank, name, localrank of MPI proces

    Note these are safe even when MPI is not running

    """
    return {
        "size": MPI.COMM_WORLD.Get_size(),
        "rank": MPI.COMM_WORLD.Get_rank(),
        "name": MPI.Get_processor_name(),
        "local_rank": (int(os.environ['OMPI_COMM_WORLD_LOCAL_RANK']) if
                       'OMPI_COMM_WORLD_LOCAL_RANK' in os.environ.keys()
                       else 0)}


    
def ReadHostFile(hostFilePath):
    hostList = []
    if hostFilePath is not None and hostFilePath != "":
        with open(hostFilePath, 'r') as f:
            for line in f:
                li=line.strip()
                if not li.startswith("#"):
                    hostList.append(line.rstrip())
    else:
        hostList.append("localhost")
    return hostList
    

def Compute(workerFn, cf):
    """Run parallel workerFn using compute configuration"""

    if cf.compute.useMPI and 'OMPI_COMM_WORLD_LOCAL_RANK' not in os.environ:        
        # is it is not already started as an MPI process
        # get hostlists from hostfile and pass it to Spawn 
        Spawn(nproc=min(cf.compute.numProcesses,
                        cf.study.numSubjects),
              hostList=ReadHostFile(cf.compute.hostFile))
    elif 'OMPI_COMM_WORLD_LOCAL_RANK' in os.environ:
        # if already running as an MPI process, may be from outside python
        # behave as if started from inside
        cf.compute.useMPI = True
        cf.compute.numProcesses = GetMPIInfo()['size']        

    if cf.compute.useMPI:
        # so that exceptions raised by mpi threads are
        # properly handled
        sys.excepthook = mpi_excepthook 

    if cf.compute.useMPI:        
        # Simple check that we're not using more host processes than GPUs
        mpiLocalRank = GetMPIInfo()['local_rank']
        if cf.compute.useCUDA:
            if mpiLocalRank >= ca.GetNumberOfCUDADevices():
                #MPI.COMM_WORLD.Abort(1) 
                raise Exception("Please don't use more host processes "+
                                "than GPUs")


            ca.SetCUDADevice(mpiLocalRank)

    # run the worker process
    workerFn(cf)

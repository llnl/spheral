#-------------------------------------------------------------------------------
# mpi
#
# This module reproduces the pyMPI interface using mpi4py.
#-------------------------------------------------------------------------------
import SpheralConfigs
import sys
from SpheralTestUtilities import globalFrame

# NOTE: this logic for disabling recv_mprobe seems to be necessary with newer
# mpi4py versions, since the LC MPI implementations apparently report matched_probes
# as supported, but seem to be broken.
import mpi4py
mpi4py.rc.recv_mprobe = False
if ("LEOS" in SpheralConfigs.component_configs()):
    mpi4py.rc.finalize = False

# Now go on as usual...
from mpi4py import MPI

#-------------------------------------------------------------------------------
# Communicator and geometry.
#-------------------------------------------------------------------------------
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
procs = comm.Get_size()

#-------------------------------------------------------------------------------
# Define the operations.
#-------------------------------------------------------------------------------
MIN = MPI.MIN
MAX = MPI.MAX
SUM = MPI.SUM
MINLOC = MPI.MINLOC
MAXLOC = MPI.MAXLOC

#-------------------------------------------------------------------------------
# Prepare files to keep the stdout and stderr streams in.
# The pyMPI defaults are only rank 0 writes stdout, but all
# processes write stderr.
#-------------------------------------------------------------------------------
globalscope = globalFrame().f_globals
if rank > 0:
    exec("""
import sys
__mpi_stdoutfile__ = open("/dev/null", "w")
sys.stdout = __mpi_stdoutfile__
""", globalscope)

#-------------------------------------------------------------------------------
# Allow us to be discriminated against as not the fake MPI package
#-------------------------------------------------------------------------------
def is_fake_mpi():
    return False

#-------------------------------------------------------------------------------
# A common helper to convert vector_of_* types to lists for communication
#-------------------------------------------------------------------------------
def __listify(obj):
    if hasattr(obj, "__qualname__") and "vector_of" in obj.__qualname__:
        return list(obj)
    else:
        return obj

#-------------------------------------------------------------------------------
# send
#-------------------------------------------------------------------------------
def send(obj, dest=0, tag=100):
    comm.send(obj=__listify(obj), dest=dest, tag=tag)

#-------------------------------------------------------------------------------
# recv
#-------------------------------------------------------------------------------
def recv(source=0, tag=100):
    return (comm.recv(source=source, tag=tag), )

#-------------------------------------------------------------------------------
# isend
#-------------------------------------------------------------------------------
def isend(obj, dest=0, tag=100):
    return comm.isend(obj=__listify(obj), dest=dest, tag=tag)

#-------------------------------------------------------------------------------
# reduce
#-------------------------------------------------------------------------------
def reduce(obj, op=SUM, root=0):
    return comm.reduce(sendobj=__listify(obj), op=op, root=root)

#-------------------------------------------------------------------------------
# allreduce
#-------------------------------------------------------------------------------
def allreduce(obj, op=SUM):
    return comm.allreduce(sendobj=__listify(obj), op=op)

#-------------------------------------------------------------------------------
# gather
#-------------------------------------------------------------------------------
def gather(obj, root=0):
    return comm.gather(sendobj=__listify(obj), root=root)

#-------------------------------------------------------------------------------
# allgather
#-------------------------------------------------------------------------------
def allgather(obj):
    return comm.allgather(sendobj=__listify(obj))

#-------------------------------------------------------------------------------
# bcast
#-------------------------------------------------------------------------------
def bcast(obj, root=0):
    return comm.bcast(__listify(obj), root=root)

#-------------------------------------------------------------------------------
# barrier
#-------------------------------------------------------------------------------
def barrier():
    comm.barrier()

#-------------------------------------------------------------------------------
# synchronizeQueuedOutput
#-------------------------------------------------------------------------------
def synchronizeQueuedOutput(stdoutfile = None,
                            stderrfile = None):
    if stdoutfile == None:
        exec("import sys; sys.stdout = sys.__stdout__", globalscope)
    else:
        exec("__mpi_stdoutfile__ = open(%s, 'w'); sys.stdout = __mpi_stdoutfile__" % stdoutfile,
             globalscope)

    if stderrfile == None:
        exec("import sys; sys.stderr = sys.__stderr__", globalscope)
    else:
        exec("__mpi_stderrfile__ = open(%s, 'w'); sys.stderr = __mpi_stderrfile__" % stderrfile,
             globalscope)

    return

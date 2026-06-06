#-------------------------------------------------------------------------------
# A fake version of the pyMPI mpi module for use in serial runs.
#-------------------------------------------------------------------------------
print("Invoking fake mpi module.")
rank = 0
procs = 1
MIN = -1
MAX = -2
SUM = -3
MINLOC = -4
MAXLOC = -5

def is_fake_mpi():
    return True

def reduce(var, op):
    return var

def allreduce(var, op):
    return var

def gather(obj, root=0):
    return [obj,]

def allgather(obj, op):
    return [obj,]

def bcast(obj, root=0):
    return obj

def barrier():
    return


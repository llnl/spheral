from PYB11Generator import *
from LinearSolver import *

@PYB11holder("std::shared_ptr")
class HypreLinearSolver(LinearSolver):
    def pyinit(self,
               options = "std::shared_ptr<HypreOptions>"):
        "Hypre solver"

PYB11inject(LinearSolverAbstractMethods, HypreLinearSolver, virtual = True)

from PYB11Generator import *
from LinearSolver import *

@PYB11holder("std::shared_ptr")
class EigenLinearSolver(LinearSolver):
    def pyinit(self,
               options = "std::shared_ptr<EigenOptions>"):
        "Eigen solver"
        
PYB11inject(LinearSolverAbstractMethods, EigenLinearSolver, virtual = True)

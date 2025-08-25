from PYB11Generator import *

@PYB11ignore
class LinearSolverAbstractMethods:
    @PYB11cppname("initialize")
    def initializePtr(self,
                    numLocVal = "const unsigned",
                    firstGlobInd = "const unsigned",
                    numValsPerRow = "const unsigned*"):
        "Create matrices and vectors"
        return "void"

    def beginFill(self):
        "Begin matrix fill"
        return "void"

    @PYB11cppname("setMatRow")
    def setMatRowPtr(self,
                     numVals = "const unsigned",
                     globRowInd = "const unsigned",
                     globColInd = "const unsigned*",
                     colVal = "const double*"):
        "Set values in matrix"
        return "void"

    @PYB11cppname("setMatRows")
    def setMatRowsPtr(self,
                      numRows = "const unsigned",
                      numColsPerRow = "const unsigned*",
                      globRowInd = "const unsigned*",
                      globColInd = "const unsigned*",
                      colVal = "const double*"):
        "Set multiple rows in matrix"
        return "void"

    def assemble(self):
        "Assemble matrix after fill"
        return "void"
      
    def finalize(self):
        "Create solver/preconditioner"
        return "void"

    @PYB11cppname("set")
    def setPtr(self,
               c = "const LinearSolver::Component",
               numVals = "const unsigned",
               firstGlobInd = "const unsigned",
               val = "const double*"):
        "Set values from external vector"
        return "void"

    @PYB11cppname("get")
    def getPtr(self,
               c = "const LinearSolver::Component",
               numVals = "const unsigned",
               firstGlobInd = "const unsigned",
               val = "double*"):
        "Get values from internal vector"
        return "void"

    def solve(self):
        "Solve the system of equations in place"
        return "void"

    def multiply(self):
        "Multiply the matrix by a vector"
        return "void"

@PYB11holder("std::shared_ptr")
class LinearSolver:
    Component = PYB11enum(("LHS",
                           "RHS"),
                          export_values = True,
                          doc = "Choose which vector to use when applying a function")
    
    def pyinit(self):
        "Pure virtual linear solver class"

    @PYB11implementation("""[](LinearSolver& self,
                               const unsigned numLocVal,
                               const unsigned firstGlobInd,
                               const std::vector<unsigned>& numValsPerRow) {
                               self.initialize(numLocVal, firstGlobInd, &numValsPerRow[0]); }""")
    @PYB11pyname("initialize")
    def initializeVec(self,
                      numLocVal = "const unsigned",
                      firstGlobInd = "const unsigned",
                      numValsPerRow = "const std::vector<unsigned>&"):
        "Create matrices and vectors"
        return "void"

    @PYB11implementation("""[](LinearSolver& self,
                               const unsigned numVals,
                               const unsigned globRowInd,
                               const std::vector<unsigned>& globColInd,
                               const std::vector<double>& colVal) {
                               self.setMatRow(numVals, globRowInd, &globColInd[0], &colVal[0]); }""")
    @PYB11pyname("setMatRow")
    def setMatRowVec(self,
                     numVals = "const unsigned",
                     globRowInd = "const unsigned",
                     globColInd = "const vector<unsigned>&",
                     colVal = "const std::vector<double>&"):
        "Set values in matrix"
        return "void"

    @PYB11implementation("""[](LinearSolver& self,
                               const unsigned numRows,
                               const std::vector<unsigned>& numColsPerRow,
                               const std::vector<unsigned>& globRowInd,
                               const std::vector<unsigned>& globColInd,
                               const std::vector<double>& colVal) {
                               self.setMatRows(numRows, &numColsPerRow[0], &globRowInd[0], &globColInd[0], &colVal[0]); }""")
    @PYB11pyname("setMatRows")
    def setMatRowsVec(self,
                      numRows = "const unsigned",
                      numColsPerRow = "const vector<unsigned>&",
                      globRowInd = "const vector<unsigned>&",
                      globColInd = "const vector<unsigned>&",
                      colVal = "const std::vector<double>&"):
        "Set multiple rows in matrix"
        return "void"

    @PYB11implementation("""[](LinearSolver& self,
                               const LinearSolver::Component c,
                               const unsigned numVals,
                               const unsigned firstGlobInd,
                               const std::vector<double>& val) {
                               self.set(c, numVals, firstGlobInd, &val[0]); }""")
    @PYB11pyname("set")
    def setVec(self,
               c = "const LinearSolver::Component",
               numVals = "const unsigned",
               firstGlobInd = "const unsigned",
               val = "const std::vector<double>&"):
        "Set values in RHS"
        return "void"
    
    @PYB11implementation("""[](LinearSolver& self,
                               const LinearSolver::Component c,
                               const unsigned numVals,
                               const unsigned firstGlobInd) {
                               std::vector<double> result(numVals);
                               self.get(c, numVals, firstGlobInd, &result[0]);
                               return result; }""")
    @PYB11pyname("get")
    def getVec(self,
               c = "const LinearSolver::Component",
               numVals = "const unsigned",
               firstGlobInd = "const unsigned"):
        "Get values from RHS"
        return "std::vector<double>"
    
    description = PYB11property(returnType = "std::string",
                                getter = "getDescription",
                                setter = "setDescription")

    @PYB11virtual
    @PYB11const
    def statistics(self):
        return "std::vector<std::shared_ptr<IncrementalStatistic<double>>>"

PYB11inject(LinearSolverAbstractMethods, LinearSolver, pure_virtual = True)

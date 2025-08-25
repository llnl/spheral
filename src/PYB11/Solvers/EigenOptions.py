from PYB11Generator import *

@PYB11holder("std::shared_ptr")
class EigenOptions:
    def pyinit(self):
        "Holds the options for Eigen solvers"

    qr = PYB11readwrite()

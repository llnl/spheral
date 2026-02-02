#-------------------------------------------------------------------------------
# The functors for integration and such
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11namespace("Spheral::PythonBoundFunctors")
@PYB11holder("std::shared_ptr")
@PYB11template("argT", "retT")
class SpheralFunctor:
    def pyinit(self):
        return

    @PYB11pure_virtual
    @PYB11const
    def __call__(self, x1="const %(argT)s&"):
        "Required operator() to map %(argT)s --> %(retT)s"
        return "%(retT)s"

@PYB11namespace("Spheral::PythonBoundFunctors")
@PYB11holder("std::shared_ptr")
@PYB11template("argT1", "argT2", "retT")
class Spheral2ArgFunctor:
    def pyinit(self):
        return

    @PYB11pure_virtual
    @PYB11const
    def __call__(self,
                 x1 = "const %(argT1)s&",
                 x2 = "const %(argT2)s&"):
        "Required operator() to map %(argT1)s %(argT2)s --> %(retT)s"
        return "%(retT)s"

@PYB11namespace("Spheral::PythonBoundFunctors")
@PYB11holder("std::shared_ptr")
@PYB11template("argT1", "argT2", "argT3", "retT")
class Spheral3ArgFunctor:
    def pyinit(self):
        return

    @PYB11pure_virtual
    @PYB11const
    def __call__(self,
                 x1 = "const %(argT1)s&",
                 x2 = "const %(argT2)s&",
                 x3 = "const %(argT3)s&"):
        "Required operator() to map %(argT1)s %(argT2)s %(argT3)s --> %(retT)s"
        return "%(retT)s"

@PYB11namespace("Spheral::PythonBoundFunctors")
@PYB11holder("std::shared_ptr")
@PYB11template("argT1", "argT2", "argT3", "argT4", "retT")
class Spheral4ArgFunctor:
    def pyinit(self):
        return

    @PYB11pure_virtual
    @PYB11const
    def __call__(self,
                 x1 = "const %(argT1)s&",
                 x2 = "const %(argT2)s&",
                 x3 = "const %(argT3)s&",
                 x4 = "const %(argT4)s&"):
        "Required operator() to map %(argT1)s %(argT2)s %(argT3)s %(argT4)s --> %(retT)s"
        return "%(retT)s"

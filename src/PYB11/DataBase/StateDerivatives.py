#-------------------------------------------------------------------------------
# StateDerivatives
#-------------------------------------------------------------------------------
from PYB11Generator import *
from StateBase import *

@PYB11template("Dimension")
class StateDerivatives(StateBase):

    PYB11typedefs = """
    using Scalar =typename %(Dimension)s::Scalar;
    using Vector =typename %(Dimension)s::Vector;
    using Tensor =typename %(Dimension)s::Tensor;
    using SymTensor =typename %(Dimension)s::SymTensor;
    using KeyType =typename StateBase<%(Dimension)s>::KeyType;
    using FieldName =typename StateBase<%(Dimension)s>::FieldName;
    using MeshPtr =typename StateBase<%(Dimension)s>::MeshPtr;
    using PackageList =typename StateDerivatives<%(Dimension)s>::PackageList;
"""

    #...........................................................................
    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                dataBase = "DataBase<%(Dimension)s>&",
                physicsPackages = "PackageList&"):
        "Construct using the Physics::registerDerivatives methods of the packages"

    #...........................................................................
    # Operators
    def __eq__(self):
        return

    #...........................................................................
    # Methods
    def Zero(self):
        "Set all derivative Fields to zero"
        return "void"

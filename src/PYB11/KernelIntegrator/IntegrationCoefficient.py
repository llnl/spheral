from PYB11Generator import *

@PYB11template("Dimension", "CoefficientType")
@PYB11holder("std::shared_ptr")
class IntegrationCoefficient:
    def pyinit(self):
        "Coefficient for integrals"

    @PYB11const
    @PYB11pure_virtual
    def evaluateCoefficient(self,
                            kid = "const KernelIntegrationData<%(Dimension)s>&"):
        "Return a value for the coefficient"
        return "%(CoefficientType)s"

@PYB11template("Dimension", "CoefficientType")
@PYB11holder("std::shared_ptr")
class ConstantIntegrationCoefficient(IntegrationCoefficient):
    def pyinit(self):
        "Returns a constant coefficient"
    def pyinit1(self,
                coeff = "%(CoefficientType)s"):
        "Returns a constant coefficient"

    data = PYB11property(doc = "The coefficient",
                         getter = "getData",
                         setter = "setData")
    
    @PYB11const
    @PYB11virtual
    def evaluateCoefficient(self,
                            kid = "const KernelIntegrationData<%(Dimension)s>&"):
        "Return a value for the coefficient"
        return "%(CoefficientType)s"
    
@PYB11template("Dimension", "CoefficientType")
@PYB11holder("std::shared_ptr")
class FieldListIntegrationCoefficient(IntegrationCoefficient):
    def pyinit(self):
        "Returns the value of the FieldList interpolated in the integration region"
    def pyinit1(self,
                data = "const FieldList<%(Dimension)s, %(CoefficientType)s>&"):
        "Returns the value of the FieldList interpolated in the integration region"

    data = PYB11property(doc = "The FieldList",
                         getter = "getData",
                         setter = "setData")

    @PYB11virtual
    def setDataPoint(self,
                     nodeListi = "int",
                     nodei = "int",
                     data = "const %(CoefficientType)s&"):
        "Set a point of data"
        return "void"
    
    @PYB11const
    @PYB11virtual
    def evaluateCoefficient(self,
                            kid = "const KernelIntegrationData<%(Dimension)s>&"):
        "Return a value for the coefficient"
        return "%(CoefficientType)s"


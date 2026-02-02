from PYB11Generator import *

@PYB11template("Dimension", "CoefficientType")
@PYB11holder("std::shared_ptr")
class IntegralDependsOnCoefficient:
    def pyinit(self):
        "Choose coefficient for integral (coefficient defaults to one)"

    coefficient = PYB11property(doc="The coefficient",
                                getter = "getCoefficient",
                                setter = "setCoefficient")

@PYB11template("Dimension", "CoefficientType")
@PYB11holder("std::shared_ptr")
class IntegralDependsOnFieldListCoefficient:
    def pyinit(self):
        "Choose coefficient for integral"

    coefficient = PYB11property(doc="The coefficient",
                                getter = "getCoefficient",
                                setter = "setCoefficient")

@PYB11template("Dimension")
@PYB11holder("std::shared_ptr")
class KernelIntegralBase:
    def pyinit(self):
        "Base integral class"

    @PYB11pure_virtual
    @PYB11const
    def bilinear(self):
        "Does integral need bilinear indexing?"
        return "bool"

    @PYB11pure_virtual
    @PYB11const
    def volume(self):
        "Does integral have volume integral terms?"
        return "bool"

    @PYB11pure_virtual
    @PYB11const
    def surface(self):
        "Does integral have surface integral terms?"
        return "bool"
    
    @PYB11pure_virtual
    def initialize(self,
                   flatConnectivity = "const FlatConnectivity<%(Dimension)s>&"):
        "Initialize the integral"
        return "void"
    
    @PYB11virtual
    def finalize(self,
                 flatConnectivity = "const FlatConnectivity<%(Dimension)s>&"):
        "Finalize the integral"
        return "void"
    
    @PYB11virtual
    def addToIntegral(self,
                      data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"

    @PYB11virtual
    def addToSurfaceIntegral(self,
                             data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"
    
@PYB11template("Dimension", "DataType")
@PYB11holder("std::shared_ptr")
class KernelIntegral(KernelIntegralBase):
    PYB11typedefs = '''
    typedef std::vector<%(DataType)s> StorageType;
'''
    def pyinit(self):
        "Integral of kernels with return type"

    @PYB11pure_virtual
    @PYB11const
    def bilinear(self):
        "Does integral need bilinear indexing?"
        return "bool"

    @PYB11pure_virtual
    @PYB11const
    def volume(self):
        "Does integral have volume integral terms?"
        return "bool"

    @PYB11pure_virtual
    @PYB11const
    def surface(self):
        "Does integral have surface integral terms?"
        return "bool"
    
    @PYB11const
    @PYB11pure_virtual
    def values(self):
        "Return integral values"
        return "const StorageType&"
        
    @PYB11pure_virtual
    def initialize(self,
                   flatConnectivity = "const FlatConnectivity<%(Dimension)s>&"):
        "Initialize the integral"
        return "void"

    @PYB11virtual
    def addToIntegral(self,
                      data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"

    @PYB11virtual
    def addToSurfaceIntegral(self,
                             data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"
    
@PYB11template("Dimension", "DataType")
@PYB11holder("std::shared_ptr")
class LinearIntegral(KernelIntegral):
    PYB11typedefs = '''
    typedef std::vector<%(DataType)s> StorageType;
'''
    
    def pyinit(self):
        "Linear integral"

    @PYB11virtual
    @PYB11const
    def bilinear(self):
        "Does integral need bilinear indexing?"
        return "bool"

    @PYB11pure_virtual
    @PYB11const
    def volume(self):
        "Does integral have volume integral terms?"
        return "bool"

    @PYB11pure_virtual
    @PYB11const
    def surface(self):
        "Does integral have surface integral terms?"
        return "bool"
    
    @PYB11const
    @PYB11virtual
    def values(self):
        "Return integral values"
        return "const StorageType&"
        
    @PYB11virtual
    def initialize(self,
                   flatConnectivity = "const FlatConnectivity<%(Dimension)s>&"):
        "Initialize the integral"
        return "void"
    
    @PYB11virtual
    def addToIntegral(self,
                      data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"
    
@PYB11template("Dimension", "BaseDataType")
@PYB11template_dict({"DataType" : "std::vector<%(BaseDataType)s>"})
@PYB11holder("std::shared_ptr")
class BilinearIntegral(KernelIntegral):
    PYB11typedefs = '''
    typedef std::vector<std::vector<%(BaseDataType)s>> StorageType;
'''
    def pyinit(self):
        "Bilinear integral"

    @PYB11virtual
    @PYB11const
    def bilinear(self):
        "Does integral need bilinear indexing?"
        return "bool"

    @PYB11pure_virtual
    @PYB11const
    def volume(self):
        "Does integral have volume integral terms?"
        return "bool"

    @PYB11pure_virtual
    @PYB11const
    def surface(self):
        "Does integral have surface integral terms?"
        return "bool"
    
    @PYB11const
    @PYB11virtual
    def values(self):
        "Return integral values"
        return "const StorageType&"
        
    @PYB11virtual
    def initialize(self,
                   flatConnectivity = "const FlatConnectivity<%(Dimension)s>&"):
        "Initialize the integral"
        return "void"

    @PYB11virtual
    def addToIntegral(self,
                      data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"

@PYB11template("Dimension", "BaseDataType")
@PYB11template_dict({"DataType" : "std::vector<%(BaseDataType)s>"})
@PYB11holder("std::shared_ptr")
class LinearSurfaceDependentIntegral(KernelIntegral):
    PYB11typedefs = '''
    typedef std::vector<std::vector<%(BaseDataType)s>> StorageType;
'''
    def pyinit(self):
        "Linear surface integral"

    @PYB11virtual
    @PYB11const
    def bilinear(self):
        "Does integral need bilinear indexing?"
        return "bool"

    @PYB11virtual
    @PYB11const
    def volume(self):
        "Does integral have volume integral terms?"
        return "bool"

    @PYB11virtual
    @PYB11const
    def surface(self):
        "Does integral have surface integral terms?"
        return "bool"
    
    @PYB11const
    @PYB11virtual
    def values(self):
        "Return integral values"
        return "const StorageType&"
        
    @PYB11virtual
    def initialize(self,
                   flatConnectivity = "const FlatConnectivity<%(Dimension)s>&"):
        "Initialize the integral"
        return "void"

    @PYB11virtual
    def addToSurfaceIntegral(self,
                             data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"

@PYB11template("Dimension", "BaseDataType")
@PYB11template_dict({"DataType" : "std::vector<%(BaseDataType)s>"})
@PYB11holder("std::shared_ptr")
class BilinearSurfaceDependentIntegral(KernelIntegral):
    PYB11typedefs = '''
    typedef std::vector<std::vector<%(BaseDataType)s>> StorageType;
'''
    def pyinit(self):
        "Bilinear surface integral"

    @PYB11virtual
    @PYB11const
    def bilinear(self):
        "Does integral need bilinear indexing?"
        return "bool"

    @PYB11virtual
    @PYB11const
    def volume(self):
        "Does integral have volume integral terms?"
        return "bool"

    @PYB11virtual
    @PYB11const
    def surface(self):
        "Does integral have surface integral terms?"
        return "bool"
    
    @PYB11const
    @PYB11virtual
    def values(self):
        "Return integral values"
        return "const StorageType&"
        
    @PYB11virtual
    def initialize(self,
                   flatConnectivity = "const FlatConnectivity<%(Dimension)s>&"):
        "Initialize the integral"
        return "void"

    @PYB11virtual
    def addToSurfaceIntegral(self,
                             data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"
        
@PYB11template("Dimension")
@PYB11template_dict({"DataType" : "typename %(Dimension)s::Scalar",
                     "CoefficientType" : "typename %(Dimension)s::Scalar"})
@PYB11holder("std::shared_ptr")
class LinearKernel(LinearIntegral,
                   IntegralDependsOnCoefficient):
    def pyinit(self):
        "Integral of kernel with a functional coefficient"

    @PYB11virtual
    @PYB11const
    def volume(self):
        "Does integral have volume integral terms?"
        return "bool"

    @PYB11virtual
    @PYB11const
    def surface(self):
        "Does integral have surface integral terms?"
        return "bool"
    
    @PYB11virtual
    def addToIntegral(self,
                      data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"

@PYB11template("Dimension")
@PYB11template_dict({"DataType" : "typename %(Dimension)s::Vector",
                     "CoefficientType" : "typename %(Dimension)s::Scalar"})
@PYB11holder("std::shared_ptr")
class LinearGrad(LinearIntegral,
                 IntegralDependsOnCoefficient):
    def pyinit(self):
        "Integral of kernel with a functional coefficient"

    @PYB11virtual
    @PYB11const
    def volume(self):
        "Does integral have volume integral terms?"
        return "bool"

    @PYB11virtual
    @PYB11const
    def surface(self):
        "Does integral have surface integral terms?"
        return "bool"
    
    @PYB11virtual
    def addToIntegral(self,
                      data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"

@PYB11template("Dimension")
@PYB11template_dict({"DataType" : "std::vector<typename %(Dimension)s::Scalar>",
                     "CoefficientType" : "std::vector<typename %(Dimension)s::Scalar>"})
@PYB11holder("std::shared_ptr")
class LinearKernelStdVector(LinearIntegral,
                            IntegralDependsOnCoefficient):
    def pyinit(self,
               size = "int"):
        "Integral of kernel with a functional coefficient"

    @PYB11virtual
    @PYB11const
    def volume(self):
        "Does integral have volume integral terms?"
        return "bool"

    @PYB11virtual
    @PYB11const
    def surface(self):
        "Does integral have surface integral terms?"
        return "bool"
    
    @PYB11virtual
    def addToIntegral(self,
                      data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"
    
    @PYB11virtual
    def initialize(self,
                   flatConnectivity = "const FlatConnectivity<%(Dimension)s>&"):
        "Initialize the integral"
        return "void"

@PYB11template("Dimension")
@PYB11template_dict({"DataType" : "std::vector<typename %(Dimension)s::Vector>",
                     "CoefficientType" : "std::vector<typename %(Dimension)s::Scalar>"})
@PYB11holder("std::shared_ptr")
class LinearGradStdVector(LinearIntegral,
                          IntegralDependsOnCoefficient):
    def pyinit(self,
               size = "int"):
        "Integral of kernel with a functional coefficient"

    @PYB11virtual
    @PYB11const
    def volume(self):
        "Does integral have volume integral terms?"
        return "bool"

    @PYB11virtual
    @PYB11const
    def surface(self):
        "Does integral have surface integral terms?"
        return "bool"
    
    @PYB11virtual
    def addToIntegral(self,
                      data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"
    
    @PYB11virtual
    def initialize(self,
                   flatConnectivity = "const FlatConnectivity<%(Dimension)s>&"):
        "Initialize the integral"
        return "void"
    
@PYB11template("Dimension")
@PYB11template_dict({"DataType" : "typename %(Dimension)s::Vector",
                     "CoefficientType" : "typename %(Dimension)s::Vector"})
@PYB11holder("std::shared_ptr")
class LinearKernelVector(LinearIntegral,
                         IntegralDependsOnCoefficient):
    def pyinit(self):
        "Integral of kernel with a functional coefficient"

    @PYB11virtual
    @PYB11const
    def volume(self):
        "Does integral have volume integral terms?"
        return "bool"

    @PYB11virtual
    @PYB11const
    def surface(self):
        "Does integral have surface integral terms?"
        return "bool"
    
    @PYB11virtual
    def addToIntegral(self,
                      data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"
    
@PYB11template("Dimension")
@PYB11template_dict({"BaseDataType" : "typename %(Dimension)s::Scalar",
                     "CoefficientType" : "typename %(Dimension)s::Scalar"})
@PYB11holder("std::shared_ptr")
class BilinearKernelKernel(BilinearIntegral,
                           IntegralDependsOnCoefficient):
    def pyinit(self):
        "Bilinear integral of kernel and kernel"
        
    @PYB11virtual
    @PYB11const
    def volume(self):
        "Does integral have volume integral terms?"
        return "bool"

    @PYB11virtual
    @PYB11const
    def surface(self):
        "Does integral have surface integral terms?"
        return "bool"
    
    @PYB11virtual
    def addToIntegral(self,
                      data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"

@PYB11template("Dimension")
@PYB11template_dict({"BaseDataType" : "typename %(Dimension)s::Vector",
                     "CoefficientType" : "typename %(Dimension)s::Scalar"})
@PYB11holder("std::shared_ptr")
class BilinearGradKernel(BilinearIntegral,
                         IntegralDependsOnCoefficient):
    def pyinit(self):
        "Bilinear integral of grad kernel and kernel"

    @PYB11virtual
    @PYB11const
    def volume(self):
        "Does integral have volume integral terms?"
        return "bool"

    @PYB11virtual
    @PYB11const
    def surface(self):
        "Does integral have surface integral terms?"
        return "bool"
    
    @PYB11virtual
    def addToIntegral(self,
                      data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"

@PYB11template("Dimension")
@PYB11template_dict({"BaseDataType" : "typename %(Dimension)s::Vector",
                     "CoefficientType" : "typename %(Dimension)s::Scalar"})
@PYB11holder("std::shared_ptr")
class BilinearKernelGrad(BilinearIntegral,
                         IntegralDependsOnCoefficient):
    def pyinit(self):
        "Bilinear integral of grad kernel and kernel"

    @PYB11virtual
    @PYB11const
    def volume(self):
        "Does integral have volume integral terms?"
        return "bool"

    @PYB11virtual
    @PYB11const
    def surface(self):
        "Does integral have surface integral terms?"
        return "bool"
    
    @PYB11virtual
    def addToIntegral(self,
                      data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"

@PYB11template("Dimension")
@PYB11template_dict({"BaseDataType" : "typename %(Dimension)s::Scalar",
                     "CoefficientType" : "typename %(Dimension)s::Scalar"})
@PYB11holder("std::shared_ptr")
class BilinearGradDotGrad(BilinearIntegral,
                          IntegralDependsOnCoefficient):
    def pyinit(self):
        "Bilinear integral of grad kernel and grad kernel"

    @PYB11virtual
    @PYB11const
    def volume(self):
        "Does integral have volume integral terms?"
        return "bool"

    @PYB11virtual
    @PYB11const
    def surface(self):
        "Does integral have surface integral terms?"
        return "bool"
    
    @PYB11virtual
    def addToIntegral(self,
                      data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"

@PYB11template("Dimension")
@PYB11template_dict({"BaseDataType" : "typename %(Dimension)s::Tensor",
                     "CoefficientType" : "typename %(Dimension)s::Scalar"})
@PYB11holder("std::shared_ptr")
class BilinearGradProdGrad(BilinearIntegral,
                           IntegralDependsOnCoefficient):
    def pyinit(self):
        "Bilinear integral of grad kernel and grad kernel"

    @PYB11virtual
    @PYB11const
    def volume(self):
        "Does integral have volume integral terms?"
        return "bool"

    @PYB11virtual
    @PYB11const
    def surface(self):
        "Does integral have surface integral terms?"
        return "bool"
    
    @PYB11virtual
    def addToIntegral(self,
                      data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"

@PYB11template("Dimension")
@PYB11template_dict({"BaseDataType" : "typename %(Dimension)s::Vector",
                     "CoefficientType" : "typename %(Dimension)s::Scalar"})
@PYB11holder("std::shared_ptr")
class BilinearSurfaceNormalKernelKernelFromGrad(BilinearIntegral,
                                                IntegralDependsOnCoefficient):
    def pyinit(self):
        "Bilinear surface integral of kernel and kernel times normal, calculated using gradients"

    @PYB11virtual
    @PYB11const
    def volume(self):
        "Does integral have volume integral terms?"
        return "bool"

    @PYB11virtual
    @PYB11const
    def surface(self):
        "Does integral have surface integral terms?"
        return "bool"
    
    @PYB11virtual
    def addToIntegral(self,
                      data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"

@PYB11template("Dimension")
@PYB11template_dict({"DataType" : "typename %(Dimension)s::Scalar",
                     "CoefficientType" : "typename %(Dimension)s::Scalar"})
@PYB11holder("std::shared_ptr")
class LinearSurfaceKernel(LinearIntegral,
                          IntegralDependsOnCoefficient):
    def pyinit(self):
        "Integral of kernel with a functional coefficient over a surface"

    @PYB11virtual
    @PYB11const
    def volume(self):
        "Does integral have volume integral terms?"
        return "bool"

    @PYB11virtual
    @PYB11const
    def surface(self):
        "Does integral have surface integral terms?"
        return "bool"
    
    @PYB11virtual
    def addToSurfaceIntegral(self,
                             data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"
    
@PYB11template("Dimension")
@PYB11template_dict({"BaseDataType" : "typename %(Dimension)s::Vector",
                     "CoefficientType" : "typename %(Dimension)s::Scalar"})
@PYB11holder("std::shared_ptr")
class LinearSurfaceNormalKernel(LinearSurfaceDependentIntegral,
                                IntegralDependsOnCoefficient):
    def pyinit(self):
        "Linear surface integral of kernel and normal"

    @PYB11virtual
    def addToIntegral(self,
                      data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"

@PYB11template("Dimension")
@PYB11template_dict({"BaseDataType" : "std::vector<typename %(Dimension)s::Vector>",
                     "CoefficientType" : "std::vector<typename %(Dimension)s::Scalar>"})
@PYB11holder("std::shared_ptr")
class LinearSurfaceNormalKernelStdVector(LinearSurfaceDependentIntegral,
                                         IntegralDependsOnCoefficient):
    def pyinit(self,
               size = "int"):
        "Linear surface integral of kernel and normal"

    @PYB11virtual
    def addToIntegral(self,
                      data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"
    
    @PYB11virtual
    def initialize(self,
                   flatConnectivity = "const FlatConnectivity<%(Dimension)s>&"):
        "Initialize the integral"
        return "void"
    
@PYB11template("Dimension")
@PYB11template_dict({"BaseDataType" : "typename %(Dimension)s::Scalar",
                     "CoefficientType" : "typename %(Dimension)s::Scalar"})
@PYB11holder("std::shared_ptr")
class BilinearSurfaceKernelKernel(BilinearIntegral,
                                  IntegralDependsOnCoefficient):
    def pyinit(self):
        "Bilinear surface integral of kernel, kernel, and normal"

    @PYB11virtual
    @PYB11const
    def volume(self):
        "Does integral have volume integral terms?"
        return "bool"

    @PYB11virtual
    @PYB11const
    def surface(self):
        "Does integral have surface integral terms?"
        return "bool"
    
    @PYB11virtual
    def addToIntegral(self,
                      data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"
    
@PYB11template("Dimension")
@PYB11template_dict({"BaseDataType" : "typename %(Dimension)s::Vector",
                     "CoefficientType" : "typename %(Dimension)s::Scalar"})
@PYB11holder("std::shared_ptr")
class BilinearSurfaceNormalKernelKernel(BilinearSurfaceDependentIntegral,
                                        IntegralDependsOnCoefficient):
    def pyinit(self):
        "Bilinear surface integral of kernel, kernel, and normal"

    @PYB11virtual
    def addToIntegral(self,
                      data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"
    
@PYB11template("Dimension")
@PYB11template_dict({"BaseDataType" : "typename %(Dimension)s::Scalar",
                     "CoefficientType" : "typename %(Dimension)s::Scalar"})
@PYB11holder("std::shared_ptr")
class BilinearSurfaceNormalKernelDotGrad(BilinearIntegral,
                                         IntegralDependsOnCoefficient):
    def pyinit(self):
        "Bilinear surface integral of kernel and normal dotted into grad"

    @PYB11virtual
    @PYB11const
    def volume(self):
        "Does integral have volume integral terms?"
        return "bool"

    @PYB11virtual
    @PYB11const
    def surface(self):
        "Does integral have surface integral terms?"
        return "bool"
    
    @PYB11virtual
    def addToIntegral(self,
                      data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"

@PYB11template("Dimension")
@PYB11template_dict({"DataType" : "typename %(Dimension)s::Scalar",
                     "CoefficientType" : "typename %(Dimension)s::Scalar"})
@PYB11holder("std::shared_ptr")
class CellCoefficient(LinearIntegral,
                      IntegralDependsOnCoefficient):
    def pyinit(self):
        "Integral of kernel with a functional coefficient"

    @PYB11virtual
    @PYB11const
    def volume(self):
        "Does integral have volume integral terms?"
        return "bool"

    @PYB11virtual
    @PYB11const
    def surface(self):
        "Does integral have surface integral terms?"
        return "bool"
    
    @PYB11virtual
    def addToIntegral(self,
                      data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"

@PYB11template("Dimension")
@PYB11template_dict({"BaseDataType" : "typename %(Dimension)s::Vector",
                     "CoefficientType" : "typename %(Dimension)s::Scalar"})
@PYB11holder("std::shared_ptr")
class SurfaceNormalCoefficient(LinearSurfaceDependentIntegral,
                               IntegralDependsOnCoefficient):
    def pyinit(self):
        "Linear surface integral of kernel and normal"

    @PYB11virtual
    def addToIntegral(self,
                      data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"
    
@PYB11template("Dimension", "BaseDataType")
@PYB11template_dict({"CoefficientType" : "typename %(Dimension)s::Scalar"})
@PYB11holder("std::shared_ptr")
class BilinearMultiplyByFieldList(BilinearIntegral,
                                  IntegralDependsOnFieldListCoefficient):
    PYB11typedefs = '''
    typedef std::vector<std::vector<%(BaseDataType)s>> StorageType;
'''
    def pyinit(self,
               baseIntegral = "std::shared_ptr<BilinearIntegral<%(Dimension)s, %(BaseDataType)s>>"):
        "Bilinear integral multiplied by a FieldList"

    @PYB11pure_virtual
    @PYB11const
    def volume(self):
        "Does integral have volume integral terms?"
        return "bool"

    @PYB11pure_virtual
    @PYB11const
    def surface(self):
        "Does integral have surface integral terms?"
        return "bool"
    
    @PYB11virtual
    def finalize(self,
                 flatConnectivity = "const FlatConnectivity<%(Dimension)s>&"):
        "Finalize the integral"
        return "void"

    @PYB11virtual
    def addToIntegral(self,
                      data = "const KernelIntegrationData<%(Dimension)s>&"):
        "Add value for a quadrature point to the integral"
        return "void"


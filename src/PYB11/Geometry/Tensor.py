from PYB11Generator import *

#-------------------------------------------------------------------------------
# Tensor (rank 2) template
#-------------------------------------------------------------------------------
@PYB11template("ndim")
class Tensor:
    "Spheral geometric tensor (rank 2: %(ndim)sx%(ndim)s) class"

    # constexpr attributes
    nDimensions = PYB11property(constexpr=True, static=True, doc="Number of dimensions")
    numElements = PYB11property(constexpr=True, static=True, doc="Number of elements stored in the type")
    zero = PYB11property(constexpr=True, static=True, doc="The zero value equivalent")
    one = PYB11property(constexpr=True, static=True, doc="The one value equivalent")

    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                rhs = "const Dim<%(ndim)s>::Tensor"):
        "Copy constructor (tensor)"

    def pyinit2(self,
                rhs = "const Dim<%(ndim)s>::SymTensor"):
        "Copy constructor (symmetric tensor)"

    def pyinit3(self,
                xx="double"):
        "Construct with element values."

    def pyinit4(self,
                xx="double", xy="double",
                yx="double", yy="double"):
        "Construct with element values."

    def pyinit5(self,
                xx="double", xy="double", xz="double",
                yx="double", yy="double", yz="double",
                zx="double", zy="double", zz="double"):
        "Construct with element values."

    # Sequence methods
    @PYB11implementation("[](const Dim<%(ndim)s>::Tensor&) { return Dim<%(ndim)s>::Tensor::numElements(); }")
    def __len__(self):
        "The size (number of elements) of the Tensor."

    @PYB11implementation("[](const Dim<%(ndim)s>::Tensor &s, size_t i) -> double { if (i >= Dim<%(ndim)s>::Tensor::numElements()) throw py::index_error(); return s[i]; }") 
    def __getitem__(self):
        "Python indexing to get an element."

    @PYB11implementation("[](Dim<%(ndim)s>::Tensor &s, size_t i, double v) { if (i >= Dim<%(ndim)s>::Tensor::numElements()) throw py::index_error(); s[i] = v; }") 
    def __setitem__(self):
        "Python indexing to set an element."

    @PYB11implementation("[](const Dim<%(ndim)s>::Tensor &s) { return py::make_iterator(s.begin(), s.end()); }, py::keep_alive<0,1>()")
    def __iter__(self):
        "Python iteration through a Tensor."

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def __call__(self,
                 row="Dim<%(ndim)s>::Tensor::size_type", 
                 col="Dim<%(ndim)s>::Tensor::size_type"):
        "Extract the (row, column) element."
        return "double"

    @PYB11pycppname("__call__")
    @PYB11implementation("[](Dim<%(ndim)s>::Tensor& self, Dim<%(ndim)s>::Tensor::size_type row, Dim<%(ndim)s>::Tensor::size_type col, double val) { self(row,col) = val; }")
    def assignCall(self,
                   row="Dim<%(ndim)s>::Tensor::size_type", 
                   col="Dim<%(ndim)s>::Tensor::size_type",
                   val="double"):
        return "void"

    # String representation
    @PYB11implementation("""
[](const Dim<%(ndim)s>::Tensor& self) {
  std::string result = "Tensor" + std::to_string(%(ndim)s) + "d(";
  for (auto val: self) result += (" " + std::to_string(val) + " ");
  result += ")";
  return result;
}""")
    def __repr__(self):
        return

    # Operators
    def __neg__(self):
        return
    def __add__(self, rhs="Dim<%(ndim)s>::Tensor()"):
        return
    def __sub__(self, rhs="Dim<%(ndim)s>::Tensor()"):
        return
    def __mul__(self, rhs="Dim<%(ndim)s>::Tensor()"):
        return
    def __iadd__(self, rhs="Dim<%(ndim)s>::Tensor()"):
        return
    def __isub__(self, rhs="Dim<%(ndim)s>::Tensor()"):
        return

    @PYB11pycppname("__add__")
    def __add__ST(self, rhs="Dim<%(ndim)s>::SymTensor()"):
        return
    @PYB11pycppname("__sub__")
    def __sub__ST(self, rhs="Dim<%(ndim)s>::SymTensor()"):
        return
    @PYB11pycppname("__mul__")
    def __mul__ST(self, rhs="Dim<%(ndim)s>::SymTensor()"):
        return
    @PYB11pycppname("__iadd__")
    def __iadd__ST(self, rhs="Dim<%(ndim)s>::SymTensor()"):
        return
    @PYB11pycppname("__isub__")
    def __isub__ST(self, rhs="Dim<%(ndim)s>::SymTensor()"):
        return

    @PYB11pycppname("__mul__")
    def __mul__f(self, rhs="double()"):
        return
    @PYB11pycppname("__rmul__")
    def __rmul__f(self, rhs="double()"):
        return
    @PYB11pycppname("__truediv__")
    def __truediv__f(self, rhs="double()"):
        return
    @PYB11pycppname("__imul__")
    def __imul__f(self, rhs="double()"):
        return
    @PYB11pycppname("__itruediv__")
    def __itruediv__f(self, rhs="double()"):
        return

    @PYB11pycppname("__mul__")
    def __mul__V(self, rhs="Dim<%(ndim)s>::Vector()"):
        return

    # Comparison
    def __eq__(self):
        return
    def __ne__(self):
        return
    def __lt__(self):
        return
    def __gt__(self):
        return
    def __le__(self):
        return
    def __ge__(self):
        return

    # Methods
    def getRow(self):
        "Extract the indexed row as a Vector%(ndim)sd."
    def getColumn(self):
        "Extract the indexed column as a Vector%(ndim)sd."
    def Zero(self):
        "Zero out the elements."
    def Identity(self):
        "Set equal to the I tensor."
    def Symmetric(self):
        "Return a SymTensor with the symmetric portion of this tensor."
    def SkewSymmetric(self):
        "Return a Tensor with the skew-symmetric portion of this tensor."
    def Transpose(self):
        "Return the transpose of this Tensor."
    def Inverse(self):
        "Return the inverse of the tensor."
    def diagonalElements(self):
        "Return a Vector%(ndim)sd with the diagonal elements of this tensor."
    def Trace(self):
        "Compute the trace (sum of diagonal elements)."
    def Determinant(self):
        "Compute the determinant of the tensor."
    @PYB11const
    def dot(self, rhs="const Dim<%(ndim)s>::Vector&"):
        "Product with a Vector%(ndim)sd."
        return "Dim<%(ndim)s>::Vector"
    @PYB11const
    @PYB11pycppname("dot")
    def dot2(self, rhs="const Dim<%(ndim)s>::Tensor&"):
        "Product with a Tensor%(ndim)sd."
        return "Dim<%(ndim)s>::Tensor"
    @PYB11const
    @PYB11pycppname("dot")
    def dot3(self, rhs="const Dim<%(ndim)s>::SymTensor&"):
        "Product with a SymTensor%(ndim)sd."
        return "Dim<%(ndim)s>::Tensor"
    @PYB11const
    def doubledot(self, rhs="const Dim<%(ndim)s>::Tensor&"):
        "Double dot contraction with another tensor (returns a scalar)."
        return "double"
    @PYB11const
    @PYB11pycppname("doubledot")
    def doubledot2(self, rhs="const Dim<%(ndim)s>::SymTensor&"):
        "Double dot contraction with a SymTensor (returns a scalar)."
        return "double"
    def selfDoubledot(self):
        "Double dot contraction with ourself (returns a scalar)."
    def square(self):
        "Compute the product with ourself as a matrix product."
    def squareElements(self):
        "Returns the element-wise square of the tensor."
    def eigenValues(self):
        "Return a Vector%(ndim)sd with the eigenvalues of this tensor."
    def rotationalTransform(self):
        "Apply the given rotational transform to the tensor."
    def maxAbsElement(self):
        "Return the maximum of the absolute values of the elements."

    # Properties
    xx = PYB11property("double", "xx", "xx", doc="The xx element.")
    xy = PYB11property("double", "xy", "xy", doc="The xy element.")
    xz = PYB11property("double", "xz", "xz", doc="The xz element.")
    yx = PYB11property("double", "yx", "yx", doc="The yx element.")
    yy = PYB11property("double", "yy", "yy", doc="The yy element.")
    yz = PYB11property("double", "yz", "yz", doc="The yz element.")
    zx = PYB11property("double", "zx", "zx", doc="The zx element.")
    zy = PYB11property("double", "zy", "zy", doc="The zy element.")
    zz = PYB11property("double", "zz", "zz", doc="The zz element.")

#-------------------------------------------------------------------------------
# Tensor instantiations.
#-------------------------------------------------------------------------------
Tensor1d = PYB11TemplateClass(Tensor,
                              template_parameters = ("1"),
                              cppname = "Dim<1>::Tensor",
                              pyname = "Tensor1d")
Tensor2d = PYB11TemplateClass(Tensor,
                              template_parameters = ("2"),
                              cppname = "Dim<2>::Tensor",
                              pyname = "Tensor2d")
Tensor3d = PYB11TemplateClass(Tensor,
                              template_parameters = ("3"),
                              cppname = "Dim<3>::Tensor",
                              pyname = "Tensor3d")

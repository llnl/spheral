from PYB11Generator import *

#-------------------------------------------------------------------------------
# SymTensor (rank 2) template
#-------------------------------------------------------------------------------
@PYB11template("ndim")
class SymTensor:
    "Spheral geometric symmetric tensor (rank 2: %(ndim)sx%(ndim)s) class"

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

    def pyinit_double(self,
                      a = "double"):
        "Construct from a single scalar (results in a diagonal tensor with 'a' along the diagonal)"

    @PYB11ignoreTest(lambda m_attrs, k_attrs: k_attrs["template_dict"]["ndim"] != "1")
    def pyinit3(self,
                xx = "double",
                yy = ("double", 0.0),
                zz = ("double", 0.0)):
        "Construct with element values (1D)"

    @PYB11ignoreTest(lambda m_attrs, k_attrs: k_attrs["template_dict"]["ndim"] != "2")
    def pyinit4(self,
                xx="double", xy="double",
                yx="double", yy="double",
                zz=("double", 0.0)):
        "Construct with element values (2D)"

    @PYB11ignoreTest(lambda m_attrs, k_attrs: k_attrs["template_dict"]["ndim"] != "3")
    def pyinit5(self,
                xx="double", xy="double", xz="double",
                yx="double", yy="double", yz="double",
                zx="double", zy="double", zz="double"):
        "Construct with element values (3D)"

    # Sequence methods
    @PYB11implementation("[](const Dim<%(ndim)s>::SymTensor&) { return Dim<%(ndim)s>::SymTensor::numElements(); }")
    def __len__(self):
        "The size (number of elements) of the SymTensor."

    @PYB11implementation("[](const Dim<%(ndim)s>::SymTensor &s, size_t i) { if (i >= Dim<%(ndim)s>::SymTensor::numElements()) throw py::index_error(); return s[i]; }") 
    @PYB11returnpolicy("reference_internal")
    def __getitem__(self):
        "Python indexing to get an element."

    @PYB11implementation("[](Dim<%(ndim)s>::SymTensor &s, size_t i, double v) { if (i >= Dim<%(ndim)s>::SymTensor::numElements()) throw py::index_error(); s[i] = v; }") 
    def __setitem__(self):
        "Python indexing to set an element."

    @PYB11implementation("[](const Dim<%(ndim)s>::SymTensor &s) { return py::make_iterator(s.begin(), s.end()); }, py::keep_alive<0,1>()")
    def __iter__(self):
        "Python iteration through a SymTensor."

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def __call__(self,
                 row="Dim<%(ndim)s>::SymTensor::size_type", 
                 col="Dim<%(ndim)s>::SymTensor::size_type"):
        "Extract the (row, column) element."
        return "double"

    @PYB11pycppname("__call__")
    @PYB11implementation("[](Dim<%(ndim)s>::SymTensor& self, Dim<%(ndim)s>::SymTensor::size_type row, Dim<%(ndim)s>::SymTensor::size_type col, double val) { self(row,col) = val; }")
    def assignCall(self,
                   row="Dim<%(ndim)s>::SymTensor::size_type", 
                   col="Dim<%(ndim)s>::SymTensor::size_type",
                   val="double"):
        return "void"

    # String representation
    @PYB11implementation("""
[](const Dim<%(ndim)s>::SymTensor& self) {
  std::string result = "SymTensor" + std::to_string(%(ndim)s) + "d(";
  for (auto val: self) result += (" " + std::to_string(val) + " ");
  result += ")";
  return result;
}""")
    def __repr__(self):
        return

    # Operators
    def __neg__(self):
        return
    def __add__(self, rhs="Dim<%(ndim)s>::SymTensor()"):
        return
    def __sub__(self, rhs="Dim<%(ndim)s>::SymTensor()"):
        return
    def __mul__(self, rhs="Dim<%(ndim)s>::SymTensor()"):
        return
    def __iadd__(self, rhs="Dim<%(ndim)s>::SymTensor()"):
        return
    def __isub__(self, rhs="Dim<%(ndim)s>::SymTensor()"):
        return

    @PYB11pycppname("__add__")
    def __add__T(self, rhs="Dim<%(ndim)s>::Tensor()"):
        return
    @PYB11pycppname("__sub__")
    def __sub__T(self, rhs="Dim<%(ndim)s>::Tensor()"):
        return
    @PYB11pycppname("__mul__")
    def __mul__T(self, rhs="Dim<%(ndim)s>::Tensor()"):
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

    # Methods special to symmetric tensor (not in tensor).
    def sqrt(self):
        "Sqrt of the symmetric tensor."
    def pow(self):
        "Raise the symmetric tensor to an arbitrary power."
    def eigenVectors(self):
        "Return an EigenStruct with the eigenvalues and eigenvectors."

    # A few methods using the full 3D elements
    @PYB11const
    def diagonalElements3D(self):
        "Return a Vector3d with the 3D diagonal"
    @PYB11const
    def Trace3D(self):
        "The full 3D trace"
    @PYB11const
    def Determinant3D(self):
        "The full 3D determinant"
    @PYB11const
    def doubledot3D(self, rhs="const Dim<%(ndim)s>::Tensor&"):
        "3D double dot contraction (a_ij * b_ji)"
        return "double"
    @PYB11cppname("doubledot3D")
    @PYB11const
    def doubledot3D_sym(self, rhs="const Dim<%(ndim)s>::SymTensor&"):
        "3D double dot contraction (a_ij * b_ji)"
        return "double"
    @PYB11const
    def selfDoubledot3D(self):
        "3D double dot contraction with ourself (a_ij * a_ji)"
        return "double"
    @PYB11const
    def maxAbsElement3D(self):
        "Return the maximum of the absolute valuses of the elements (3D)"
        return "double"
    @PYB11const
    def eigenValues3D(self):
        "Compute the full 3D eigenvalues (always returns a Vector3d)"
        return
    @PYB11const
    def eigenVectors3D(self):
        "Compute the full 3D eigenvalues & eigenvectors (always returns an EigenStruct3d)"
        return

    @PYB11ignoreTest(lambda m_attrs, k_attrs: k_attrs["template_dict"]["ndim"] == "3")
    @PYB11pycppname("dot")
    def dot3D(self, rhs="const Dim<3>::Vector&"):
        "Product with a 3D vector"
        return

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
# SymTensor instantiations.
#-------------------------------------------------------------------------------
SymTensor1d = PYB11TemplateClass(SymTensor,
                                 template_parameters = ("1"),
                                 cppname = "Dim<1>::SymTensor",
                                 pyname = "SymTensor1d")
SymTensor2d = PYB11TemplateClass(SymTensor,
                                 template_parameters = ("2"),
                                 cppname = "Dim<2>::SymTensor",
                                 pyname = "SymTensor2d")
SymTensor3d = PYB11TemplateClass(SymTensor,
                                 template_parameters = ("3"),
                                 cppname = "Dim<3>::SymTensor",
                                 pyname = "SymTensor3d")

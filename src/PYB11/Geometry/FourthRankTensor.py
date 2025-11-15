from PYB11Generator import *

#-------------------------------------------------------------------------------
# FourthRankTensor template
#-------------------------------------------------------------------------------
@PYB11template("ndim")
class FourthRankTensor:
    "Spheral fourth rank tensor (%(ndim)sx%(ndim)sx%(ndim)sx%(ndim)s) class"

    # constexpr attributes
    nrank = PYB11property(doc="Rank of the tensor")
    nDimensions = PYB11property(doc="Number of dimensions")
    numElements = PYB11property(doc="Number of elements stored in the type")
    zero = PYB11property(doc="The zero value equivalent")

    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                rhs = "const Dim<%(ndim)s>::FourthRankTensor"):
        "Copy constructor"

    def pyinit2(self,
                rhs="double"):
        "Construct setting the element values to a constant value."

    # Sequence methods
    @PYB11implementation("[](const Dim<%(ndim)s>::FourthRankTensor&) { return Dim<%(ndim)s>::FourthRankTensor::numElements(); }")
    def __len__(self):
        "The size (number of elements) of the FourthRankTensor."

    @PYB11implementation("[](const Dim<%(ndim)s>::FourthRankTensor &s, size_t i) { if (i >= Dim<%(ndim)s>::FourthRankTensor::numElements()) throw py::index_error(); return s[i]; }") 
    @PYB11returnpolicy("reference_internal")
    def __getitem__(self):
        "Python indexing to get an element."

    @PYB11implementation("[](Dim<%(ndim)s>::FourthRankTensor &s, size_t i, double v) { if (i >= Dim<%(ndim)s>::FourthRankTensor::numElements()) throw py::index_error(); s[i] = v; }") 
    def __setitem__(self):
        "Python indexing to set an element."

    @PYB11implementation("[](const Dim<%(ndim)s>::FourthRankTensor &s) { return py::make_iterator(s.begin(), s.end()); }, py::keep_alive<0,1>()")
    def __iter__(self):
        "Python iteration through a FourthRankTensor."

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def __call__(self,
                 i="Dim<%(ndim)s>::FourthRankTensor::size_type", 
                 j="Dim<%(ndim)s>::FourthRankTensor::size_type",
                 k="Dim<%(ndim)s>::FourthRankTensor::size_type",
                 m="Dim<%(ndim)s>::FourthRankTensor::size_type"):
        "Extract the (i,j,k,m) element."
        return "double"

    @PYB11pycppname("__call__")
    @PYB11implementation("""[](Dim<%(ndim)s>::FourthRankTensor& self, 
                               Dim<%(ndim)s>::FourthRankTensor::size_type i,
                               Dim<%(ndim)s>::FourthRankTensor::size_type j,
                               Dim<%(ndim)s>::FourthRankTensor::size_type k,
                               Dim<%(ndim)s>::FourthRankTensor::size_type m,
                               double val) { self(i,j,k,m) = val; }""")
    def assignCall(self,
                   i="Dim<%(ndim)s>::FourthRankTensor::size_type", 
                   j="Dim<%(ndim)s>::FourthRankTensor::size_type",
                   k="Dim<%(ndim)s>::FourthRankTensor::size_type",
                   m="Dim<%(ndim)s>::FourthRankTensor::size_type",
                   val="double"):
        return "void"

    # Methods
    def Zero(self):
        "Zero out the elements"
        return "void"

    @PYB11const
    def doubledot(self, rhs="const RankNTensor<%(ndim)s, 4, GeomFourthRankTensor<%(ndim)s>>& rhs"):
        return "double"

    @PYB11const
    def squareElements(self):
        return "const Dim<%(ndim)s>::FourthRankTensor"

    @PYB11const
    def maxAbsElement(self):
        return "double"

    # Operators
    def __neg__(self):
        return
    def __iadd__(self):
        return
    def __isub__(self):
        return
    def __add__(self):
        return
    def __sub__(self):
        return
    def __imul__(self, rhs="double()"):
        return
    def __itruediv__(self, rhs="double()"):
        return
    def __mul__(self, rhs="double()"):
        return
    def __truediv__(self, rhs="double()"):
        return
                 
    # Comparisons
    def __eq__(self):
        return
    def __ne__(self):
        return

    # String representation
    @PYB11implementation("""
[](const Dim<%(ndim)s>::FourthRankTensor& self) {
  std::string result = "FourthRankTensor" + std::to_string(%(ndim)s) + "d(";
  for (auto val: self) result += (" " + std::to_string(val) + " ");
  result += ")";
  return result;
}""")
    def __repr__(self):
        return

#-------------------------------------------------------------------------------
# FourthRankTensor instantiations.
#-------------------------------------------------------------------------------
FourthRankTensor1d = PYB11TemplateClass(FourthRankTensor,
                                       template_parameters = ("1"),
                                       cppname = "Dim<1>::FourthRankTensor",
                                       pyname = "FourthRankTensor1d")
FourthRankTensor2d = PYB11TemplateClass(FourthRankTensor,
                                       template_parameters = ("2"),
                                       cppname = "Dim<2>::FourthRankTensor",
                                       pyname = "FourthRankTensor2d")
FourthRankTensor3d = PYB11TemplateClass(FourthRankTensor,
                                       template_parameters = ("3"),
                                       cppname = "Dim<3>::FourthRankTensor",
                                       pyname = "FourthRankTensor3d")

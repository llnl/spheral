from PYB11Generator import *

#-------------------------------------------------------------------------------
# Vector template
#-------------------------------------------------------------------------------
@PYB11template("ndim")
class Vector:
    "Spheral Geometric Vector class"

    # constexpr attributes
    nDimensions = PYB11property(constexpr=True, static=True, doc="Number of dimensions")
    numElements = PYB11property(constexpr=True, static=True, doc="Number of elements stored in the type")
    zero = PYB11property(constexpr=True, static=True, doc="The zero value equivalent")
    one = PYB11property(constexpr=True, static=True, doc="The one value equivalent")

    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                rhs = "const Dim<%(ndim)s>::Vector"):
        "Copy constructor"

    def pyinit2(self,
                x = "double",
                y = ("double", "0.0"),
                z = ("double", "0.0")):
        "Construct with element values."

    # Methods
    def Zero(self):
        "Zero out the vector elements"

    @PYB11const
    def dot(self):
        "Dot (inner) product with a Vector."

    @PYB11const
    def cross(self):
        "Cross product with a Vector."

    @PYB11const
    def dyad(self):
        "Dyadic (outer) product with a Vector."

    @PYB11const
    def selfdyad(self):
        "Dyadic (outer) product with ourself."

    @PYB11const
    def unitVector(self):
        "Unit vector in the direction of this one."

    @PYB11const
    def magnitude(self):
        "The magnitude of the Vector."

    @PYB11const
    def magnitude2(self):
        "The square of the magnitude of the Vector."
        return

    @PYB11const
    def minElement(self):
        "Minimum (x,y,z) in the Vector."

    @PYB11const
    def maxElement(self):
        "Maximum (x,y,z) in the Vector."

    @PYB11const
    def maxAbsElement(self):
        "Maximum absolute element (|x|,|y|,|z|) in the Vector."

    @PYB11const
    def sumElements(self):
        "Sum of the elements (x+y+z) in the Vector."

    # Operators
    def __neg__(self):
        return
    def __add__(self):
        return
    def __sub__(self):
        return
    def __mul__(self):
        return
    def __iadd__(self):
        return
    def __isub__(self):
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

    # Sequence methods
    @PYB11implementation("[](const Dim<%(ndim)s>::Vector&) { return Dim<%(ndim)s>::nDim; }")
    def __len__(self):
        "The size (in number of coordinates) of the Vector."

    @PYB11implementation("[](const Dim<%(ndim)s>::Vector &s, size_t i) -> double { if (i >= Dim<%(ndim)s>::Vector::numElements()) throw py::index_error(); return s[i]; }") 
    @PYB11returnpolicy("reference_internal")
    def __getitem__(self):
        "Python indexing to get a coordinate."
        return "double"

    @PYB11implementation("[](Dim<%(ndim)s>::Vector &s, size_t i, double v) { if (i >= Dim<%(ndim)s>::Vector::numElements()) throw py::index_error(); s[i] = v; }") 
    def __setitem__(self):
        "Python indexing to set a coordinate."
        return "void"

    @PYB11implementation("[](const Dim<%(ndim)s>::Vector &s) { return py::make_iterator(s.begin(), s.end()); }, py::keep_alive<0,1>()")
    def __iter__(self):
        "Python iteration through a Vector."

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def __call__(self, i="Dim<%(ndim)s>::Vector::size_type"):
        "Index for a coordinate using parens."
        return "double"

    # Comparison
    @PYB11const
    def compare(self,
                rhs = "const Dim<%(ndim)s>::Vector&"):
        "Compare (-1,0,1) with a Vector."
        return "int"

    @PYB11const
    @PYB11pyname("compare")
    @PYB11cppname("compare")
    def compare1(self,
                 rhs = "const double"):
        "Compare (-1,0,1) with a double."
        return "int"

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

    # String representation
    @PYB11implementation("""
[](const Dim<%(ndim)s>::Vector& self) {
  std::string result = "Vector" + std::to_string(%(ndim)s) + "d(";
  for (auto val: self) result += (" " + std::to_string(val) + " ");
  result += ")";
  return result;
}""")
    def __repr__(self):
        return

    # # Pickle support
    # @PYB11implementation("[](const Dim<%(ndim)s>::Vector& self) { return py::make_tuple(self.x(), self.y(), self.z()); }")
    # def __getnewargs__(self):
    #     return "py::tuple"

    # @PYB11implementation("[](const Dim<%(ndim)s>::Vector& self) { return py::make_tuple(self.x(), self.y(), self.z()); }")
    # def __getstate__(self):
    #     return "py::tuple"

    # @PYB11implementation("[](Dim<%(ndim)s>::Vector& self, const py::tuple& state) { self.x(state[0]); self.y(state[1]); self.z(state[2]); }")
    # def __setstate__(self):
    #     return

    # Properties
    x = PYB11property("double", getter="x", setter="x", doc="The x coordinate.")
    y = PYB11property("double", getter="y", setter="y", doc="The y coordinate.")
    z = PYB11property("double", getter="z", setter="z", doc="The z coordinate.")

#-------------------------------------------------------------------------------
# Vector instantiations.
#-------------------------------------------------------------------------------
Vector1d = PYB11TemplateClass(Vector,
                              template_parameters = ("1"),
                              cppname = "Dim<1>::Vector",
                              pyname = "Vector1d",
                              docext = " (1D).")
Vector2d = PYB11TemplateClass(Vector,
                              template_parameters = ("2"),
                              cppname = "Dim<2>::Vector",
                              pyname = "Vector2d",
                              docext = " (2D).")
Vector3d = PYB11TemplateClass(Vector,
                              template_parameters = ("3"),
                              cppname = "Dim<3>::Vector",
                              pyname = "Vector3d",
                              docext = " (3D).")


#-------------------------------------------------------------------------------
# PythonScalarArtificialViscosity
#-------------------------------------------------------------------------------
from PYB11Generator import *
from ArtificialViscosity import *

@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralArtificialViscosity")
class PythonScalarArtificialViscosity(ArtificialViscosity):
    """Python-overridable artificial viscosity (Scalar version, CPU-only).

    Allows rapid prototyping of new viscosity models in Python.

    WARNING: Significantly slower than C++ implementations (100-1000x slower).
             Use for prototyping only, not production runs.
             Forces serial execution (no OpenMP) due to Python GIL.

    Example usage:
        from Spheral1d import *

        class MyViscosity(PythonScalarArtificialViscosity1d):
            def __init__(self, Cl, Cq, kernel):
                PythonScalarArtificialViscosity1d.__init__(self, Cl, Cq, kernel)

            def computeQPiij(self, QPiij, QPiji, Qij, Qji,
                           xi, vi, rhoi, csi,
                           xj, vj, rhoj, csj,
                           etai, etaj,
                           fCli, fCqi, fClj, fCqj):
                # Your custom viscosity implementation here
                vij = vi - vj
                xij = xi - xj
                if vij.dot(xij) < 0.0:  # Compression only
                    mui = vij.dot(etai) / (etai.magnitude2() + self.epsilon2())
                    # ... compute QPiij, QPiji, Qij, Qji
                else:
                    QPiij = 0.0
                    QPiji = 0.0
                    Qij = 0.0
                    Qji = 0.0

            def label(self):
                return "MyViscosity"

        # Use it
        av = MyViscosity(1.0, 2.0, kernel)
        hydro = SPH1d(dataBase, av, kernel)
    """

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using ThirdRankTensor = typename %(Dimension)s::ThirdRankTensor;
    using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
    using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
    using VolumeRequirements = typename Physics<%(Dimension)s>::VolumeRequirements;
    using RKRequirements = typename Physics<%(Dimension)s>::RKRequirements;
    using ConnectivityRequirements = typename Physics<%(Dimension)s>::ConnectivityRequirements;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               Clinear = "const Scalar",
               Cquadratic = "const Scalar",
               kernel = "const TableKernel<%(Dimension)s>&"):
        "PythonScalarArtificialViscosity constructor"

    #...........................................................................
    # Methods
    @PYB11pure_virtual
    @PYB11const
    def computeQPiij(self,
                     # Particle i state
                     xi    = "const Vector&",
                     vi    = "const Vector&",
                     rhoi  = "const Scalar",
                     csi   = "const Scalar",
                     # Particle j state
                     xj    = "const Vector&",
                     vj    = "const Vector&",
                     rhoj  = "const Scalar",
                     csj   = "const Scalar",
                     # Pre-computed H*dx (eta)
                     etai  = "const Vector&",
                     etaj  = "const Vector&",
                     # Viscosity coefficient multipliers
                     fCli   = "const Scalar",
                     fCqi   = "const Scalar",
                     fClj   = "const Scalar",
                     fCqj   = "const Scalar"):
        """Compute artificial viscosity for particle pair (i,j).

        Override this method in Python to implement custom viscosity models.

        Should return results a tuple: (QPiij, QPiji, Qij, Qji)

        This is a SIMPLIFIED interface with only 14 arguments (vs 20+ in full QPiij).
        All FieldList lookups and complex data extraction is handled by C++.

        Returns a tuple:
            QPiij, QPiji: Q/rho^2 viscous pressure outputs (modify in place)
            Qij, Qji: Viscous pressure Q outputs (modify in place)

        Args:
            xi, vi, rhoi, csi: Position, velocity, density, sound speed for particle i
            xj, vj, rhoj, csj: Position, velocity, density, sound speed for particle j
            etai, etaj: Pre-computed Hi*(xi-xj) and Hj*(xj-xi)
                        Use these instead of recomputing from positions
            fCli, fCqi: Linear and quadratic viscosity multipliers for particle i
            fClj, fCqj: Linear and quadratic viscosity multipliers for particle j

        Typical implementation pattern:
            vij = vi - vj
            xij = xi - xj

            if vij.dot(xij) < 0.0:  # Compression
                mui = vij.dot(etai) / (etai.magnitude2() + self.epsilon2)
                muj = vij.dot(etaj) / (etaj.magnitude2() + self.epsilon2)

                Clij = 0.5 * (fCli + fClj) * self.Cl
                Cqij = 0.5 * (fCqi + fCqj) * self.Cq

                ei = -Clij * csi * min(0.0, mui) + Cqij * min(0.0, mui)**2
                ej = -Clij * csj * min(0.0, muj) + Cqij * min(0.0, muj)**2

                return (ei / rhoi, 
                        ej / rhoj,
                        rhoi * ei,
                        rhoj * ej)
            else:
                return (0.0, 0.0, 0.0, 0.0)
        """
        return "std::tuple<Scalar, Scalar, Scalar, Scalar>"

    @PYB11virtual
    @PYB11protected
    def updateManagedPtr(self):
        "Update member data for managed pointer."
        return "void"

    @PYB11virtual
    @PYB11const
    def QPiTypeIndex(self):
        "Require ArtificialViscosities to specify the type_index of the descendant QPiType"
        return "std::type_index"

    @PYB11virtual
    @PYB11const
    def label(self):
        return "std::string"

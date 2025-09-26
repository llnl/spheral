"""
Spheral Integrator module.

Provides the base classes and implementations for time integrators in Spheral.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

HAVE_SOLVERS = "@SPHERAL_ENABLE_SUNDIALS@"

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"DataBase/DataBase.hh"',
                  '"DataBase/State.hh"',
                  '"DataBase/StateDerivatives.hh"',
                  '"Physics/Physics.hh"',
                  '"Boundary/Boundary.hh"',
                  '"FileIO/FileIO.hh"',
                  '"Integrator/Integrator.hh"',
                  '"Integrator/ImplicitIntegrator.hh"',
                  '"Integrator/PredictorCorrector.hh"',
                  '"Integrator/ForwardEuler.hh"',
                  '"Integrator/SynchronousRK2.hh"',
                  '"Integrator/SynchronousRK4.hh"',
                  '"Integrator/CheapSynchronousRK2.hh"',
                  '"Integrator/Verlet.hh"',
                  '"Integrator/BackwardEuler.hh"',
                  '"Integrator/CrankNicolson.hh"',
                  '"Integrator/ImplicitIntegrationVectorOperator.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
from Integrator import *
from ImplicitIntegrator import *
from PredictorCorrectorIntegrator import *
from ForwardEulerIntegrator import *
from SynchronousRK2Integrator import *
from SynchronousRK4Integrator import *
from CheapSynchronousRK2Integrator import *
from VerletIntegrator import *
from CrankNicolsonIntegrator import *

if HAVE_SOLVERS:
    from ImplicitIntegrationVectorOperator import *
    from BackwardEulerIntegrator import *

for ndim in dims:
    Dimension = f"Dim<{ndim}>"
    exec(f'''
Integrator{ndim}d = PYB11TemplateClass(Integrator, template_parameters="{Dimension}")
ImplicitIntegrator{ndim}d = PYB11TemplateClass(ImplicitIntegrator, template_parameters="{Dimension}")
PredictorCorrectorIntegrator{ndim}d = PYB11TemplateClass(PredictorCorrectorIntegrator, template_parameters="{Dimension}")
ForwardEulerIntegrator{ndim}d = PYB11TemplateClass(ForwardEulerIntegrator, template_parameters="{Dimension}")
SynchronousRK2Integrator{ndim}d = PYB11TemplateClass(SynchronousRK2Integrator, template_parameters="{Dimension}")
SynchronousRK4Integrator{ndim}d = PYB11TemplateClass(SynchronousRK4Integrator, template_parameters="{Dimension}")
CheapSynchronousRK2Integrator{ndim}d = PYB11TemplateClass(CheapSynchronousRK2Integrator, template_parameters="{Dimension}")
VerletIntegrator{ndim}d = PYB11TemplateClass(VerletIntegrator, template_parameters="{Dimension}")
CrankNicolsonIntegrator{ndim}d = PYB11TemplateClass(CrankNicolsonIntegrator, template_parameters="{Dimension}")

if HAVE_SOLVERS:
    ImplicitIntegrationVectorOperator{ndim}d = PYB11TemplateClass(ImplicitIntegrationVectorOperator, template_parameters="{Dimension}")
    BackwardEulerIntegrator{ndim}d = PYB11TemplateClass(BackwardEulerIntegrator, template_parameters="{Dimension}")
''')

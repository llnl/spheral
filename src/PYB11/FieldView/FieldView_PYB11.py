"""
Spheral FieldView module.

Provides the FieldView classes.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from FieldView import *
from ArithmeticFieldView import *
from MinMaxFieldView import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Geometry/Dimension.hh"',
                  '"Field/FieldView.hh"',
                  '"Utilities/FieldDataTypeTraits.hh"',
                  '"Utilities/DomainNode.hh"',
                  '"Geometry/CellFaceFlag.hh"',
                  '<vector>']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Do our dimension dependent instantiations.
#-------------------------------------------------------------------------------
for ndim in dims:
    Dimension = f"Dim<{ndim}>"
    Scalar = f"{Dimension}::Scalar"
    Vector = f"{Dimension}::Vector"
    Tensor = f"{Dimension}::Tensor"
    SymTensor = f"{Dimension}::SymTensor"
    ThirdRankTensor = f"{Dimension}::ThirdRankTensor"
    FourthRankTensor = f"{Dimension}::FourthRankTensor"
    FifthRankTensor = f"{Dimension}::FifthRankTensor"

    #...........................................................................
    # non-numeric types
    for (value, label) in ((f"{Dimension}::FacetedVolume",       "FacetedVolume"),
                           ( "std::vector<int>",                 "VectorInt"),
                           ( "std::vector<unsigned>",            "VectorUnsigned"),
                           ( "std::vector<uint64_t>",            "VectorULL"),
                           ( "std::vector<double>",              "VectorDouble"),
                           (f"std::vector<{Vector}>",            "VectorVector"),
                           (f"std::vector<{Tensor}>",            "VectorTensor"),
                           (f"std::vector<{SymTensor}>",         "VectorSymTensor"),
                           ( "std::vector<CellFaceFlag>",        "vector_of_CellFaceFlag"),
                           (f"DomainNode<{Dimension}>",          "DomainNode"),
                           (f"RKCoefficients<{Dimension}>",      "RKCoefficients"),
                           (ThirdRankTensor,  "ThirdRankTensor"),
                           (FourthRankTensor, "FourthRankTensor"),
                           (FifthRankTensor,  "FifthRankTensor")):
        exec(f'''
{label}FieldView{ndim}d = PYB11TemplateClass(FieldView, template_parameters=("{Dimension}", "{value}"))
''')

    #...........................................................................
    # arithmetic fields
    for (value, label) in (("int",            "Int"),
                           ("unsigned",       "Unsigned"),
                           ("uint64_t",       "ULL"),
                           (Vector,           "Vector"),
                           (Tensor,           "Tensor")):
        exec(f'''
{label}FieldView{ndim}d = PYB11TemplateClass(ArithmeticFieldView, template_parameters=("{Dimension}", "{value}"))
''')

    #...........................................................................
    # A few fields can apply the min/max with a scalar addtionally
    for (value, label) in ((Scalar,           "Scalar"),
                           (SymTensor,        "SymTensor")):
        exec(f'''
{label}FieldView{ndim}d = PYB11TemplateClass(MinMaxFieldView, template_parameters=("{Dimension}", "{value}"))
''')

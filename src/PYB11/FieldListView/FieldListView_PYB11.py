"""
Spheral FieldListView module.

Provides the FieldListView classes.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from FieldListView import *
from ArithmeticFieldListView import *
from MinMaxFieldListView import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Geometry/Dimension.hh"',
                  '"Field/FieldView.hh"',
                  '"Field/FieldListView.hh"',
                  '"Utilities/FieldDataTypeTraits.hh"',
                  '"Utilities/DomainNode.hh"',
                  '"Geometry/CellFaceFlag.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Do our dimension dependent instantiations.
#-------------------------------------------------------------------------------
for ndim in dims:
    Dimension = f"Dim<{ndim}>"
    Vector = f"{Dimension}::Vector"
    Tensor = f"{Dimension}::Tensor"
    SymTensor = f"{Dimension}::SymTensor"
    ThirdRankTensor = f"{Dimension}::ThirdRankTensor"
    FourthRankTensor = f"{Dimension}::FourthRankTensor"
    FifthRankTensor = f"{Dimension}::FifthRankTensor"
    FacetedVolume = f"{Dimension}::FacetedVolume"

    #...........................................................................
    # FieldListView -- non-numeric types 
    for (value, label) in (( FacetedVolume,                 "FacetedVolume"), 
                           ( "std::vector<int>",            "VectorInt"),
                           ( "std::vector<unsigned>",       "VectorUnsigned"),
                           ( "std::vector<uint64_t>",       "VectorULL"),
                           ( "std::vector<double>",         "VectorDouble"),
                           (f"std::vector<{Vector}>",       "VectorVector"),
                           (f"std::vector<{Tensor}>",       "VectorTensor"),
                           (f"std::vector<{SymTensor}>",    "VectorSymTensor"),
                           ( "std::vector<CellFaceFlag>",   "vector_of_CellFaceFlag"),
                           (f"DomainNode<{Dimension}>",     "DomainNode"),
                           (f"RKCoefficients<{Dimension}>", "RKCoefficients")):
        exec(f'''
{label}FieldListView{ndim}d = PYB11TemplateClass(FieldListView, template_parameters=("{Dimension}", "{value}"))
''')

    #...........................................................................
    # arithmetic FieldListViews
    for (value, label) in (("int",            "Int"),
                           ("unsigned",       "Unsigned"),
                           ("uint64_t",       "ULL"),
                           (Vector,           "Vector"),
                           (Tensor,           "Tensor"),
                           (ThirdRankTensor,  "ThirdRankTensor"),
                           (FourthRankTensor, "FourthRankTensor"),
                           (FifthRankTensor,  "FifthRankTensor")):
        exec(f'''
{label}FieldListView{ndim}d = PYB11TemplateClass(ArithmeticFieldListView, template_parameters=("{Dimension}", "{value}"))
''')

    #...........................................................................
    # A few FieldListView types can apply the min/max with a scalar additionally
    for (value, label) in (("double",          "Scalar"),
                           (SymTensor,        "SymTensor")):
        exec(f'''
{label}FieldListView{ndim}d = PYB11TemplateClass(MinMaxFieldListView, template_parameters=("{Dimension}", "{value}"))
''')

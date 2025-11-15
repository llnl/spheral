"""
Spheral FieldList module.

Provides the FieldList classes.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from FieldListBase import *
from FieldList import *
from ArithmeticFieldList import *
from MinMaxFieldList import *
from FieldListSet import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Geometry/Dimension.hh"',
                  '"Field/FieldBase.hh"',
                  '"Field/Field.hh"',
                  '"Field/FieldList.hh"',
                  '"Field/FieldListSet.hh"',
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
    Vector = f"{Dimension}::Vector"
    Tensor = f"{Dimension}::Tensor"
    SymTensor = f"{Dimension}::SymTensor"
    ThirdRankTensor = f"{Dimension}::ThirdRankTensor"
    FourthRankTensor = f"{Dimension}::FourthRankTensor"
    FifthRankTensor = f"{Dimension}::FifthRankTensor"
    FacetedVolume = f"{Dimension}::FacetedVolume"

    #...........................................................................
    # FieldListBase, FieldListSet
    exec(f'''
FieldListBase{ndim}d = PYB11TemplateClass(FieldListBase, template_parameters="{Dimension}")
FieldListSet{ndim}d = PYB11TemplateClass(FieldListSet, template_parameters="{Dimension}")
''')

    #...........................................................................
    # FieldList -- non-numeric types
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
{label}FieldList{ndim}d = PYB11TemplateClass(FieldList, template_parameters=("{Dimension}", "{value}"))
''')

    #...........................................................................
    # arithmetic FieldLists
    for (value, label) in ((ThirdRankTensor,  "ThirdRankTensor"),
                           (FourthRankTensor, "FourthRankTensor"),
                           (FifthRankTensor,  "FifthRankTensor")):
        exec(f'''
{label}FieldList{ndim}d = PYB11TemplateClass(ArithmeticFieldList, template_parameters=("{Dimension}", "{value}"))
''')

    #...........................................................................
    # A few FieldLists types can apply the min/max with a scalar additionally
    for (value, label) in (("int",            "Int"),
                           ("unsigned",       "Unsigned"),
                           ("uint64_t",       "ULL"),
                           ("double",         "Scalar"),
                           (Vector,           "Vector"),
                           (Tensor,           "Tensor"),
                           (SymTensor,        "SymTensor")):
        exec(f'''
{label}FieldList{ndim}d = PYB11TemplateClass(MinMaxFieldList, template_parameters=("{Dimension}", "{value}"))
''')

    #...........................................................................
    # STL collections of FieldList types
    for value, label in (("int",     "Int"),
                         ("double",  "Scalar"),
                         (Vector,    "Vector"),
                         (Tensor,    "Tensor"),
                         (SymTensor, "SymTensor")):
        exec(f'''
vector_of_{label}FieldList{ndim}d = PYB11_bind_vector("FieldList<{Dimension}, {value}>", opaque=True, local=False)
vector_of_{label}FieldListPtr{ndim}d = PYB11_bind_vector("FieldList<{Dimension}, {value}>*", opaque=True, local=False)
vector_of_vector_of_{label}FieldList{ndim}d = PYB11_bind_vector("std::vector<FieldList<{Dimension}, {value}>>", opaque=True, local=False)
''')

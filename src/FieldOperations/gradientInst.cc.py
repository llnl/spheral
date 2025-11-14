text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "FieldOperations/gradient.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
"""

for Value in ("Scalar", "Vector"):
    text += """
template 
FieldList<%%(Dim)s, MathTraits<%%(Dim)s, %(Value)s>::GradientType> 
gradient<%%(Dim)s, %(Value)s>(const FieldList<%%(Dim)s, %(Value)s>& fieldList,
                              const FieldList<%%(Dim)s, %%(Vector)s>& position,
                              const FieldList<%%(Dim)s, %%(Scalar)s>& weight,
                              const FieldList<%%(Dim)s, %%(Scalar)s>& mass,
                              const FieldList<%%(Dim)s, %%(Scalar)s>& density,
                              const FieldList<%%(Dim)s, %%(SymTensor)s>& Hfield,
                              const TableKernel< %%(Dim)s >& kernel);

template 
FieldList<%%(Dim)s, std::vector<MathTraits<%%(Dim)s, %(Value)s>::GradientType>>
gradient<%%(Dim)s, %(Value)s>(const FieldList<%%(Dim)s, std::vector<%(Value)s>>& fieldList,
                              const FieldList<%%(Dim)s, %%(Vector)s>& position,
                              const FieldList<%%(Dim)s, %%(Scalar)s>& weight,
                              const FieldList<%%(Dim)s, %%(Scalar)s>& mass,
                              const FieldList<%%(Dim)s, %%(Scalar)s>& density,
                              const FieldList<%%(Dim)s, %%(SymTensor)s>& Hfield,
                              const TableKernel< %%(Dim)s >& kernel);

template 
void
gradientPairs<%%(Dim)s, %(Value)s>(FieldList<%%(Dim)s, MathTraits<%%(Dim)s, %(Value)s>::GradientType>& result,
                                   const FieldList<%%(Dim)s, %(Value)s>& field,
                                   const FieldList<%%(Dim)s, %%(Vector)s>& position,
                                   const FieldList<%%(Dim)s, %%(Scalar)s>& weight,
                                   const FieldList<%%(Dim)s, %%(SymTensor)s>& Hfield,
                                   const ConnectivityMap<%%(Dim)s>& conn,
                                   const TableKernel<%%(Dim)s>& kernel);
""" % {"Value" : "%(" + Value + ")s"}

text += """
//============================== limiter() ==============================
template 
FieldList<%(Dim)s, %(SymTensor)s> 
limiter<%(Dim)s, %(Scalar)s>(const FieldList<%(Dim)s, %(Scalar)s>& fieldList,
                             const FieldList<%(Dim)s, MathTraits<%(Dim)s, %(Scalar)s>::GradientType>& gradient,
                             const FieldList<%(Dim)s, %(Vector)s>& position,
                             const FieldList<%(Dim)s, %(SymTensor)s>& Hfield,
                             const TableKernel< %(Dim)s >& kernel);
template 
FieldList<%(Dim)s, %(SymTensor)s> 
limiter<%(Dim)s, %(Vector)s>(const FieldList<%(Dim)s, %(Vector)s>& fieldList,
                             const FieldList<%(Dim)s, MathTraits<%(Dim)s, %(Vector)s>::GradientType>& gradient,
                             const FieldList<%(Dim)s, %(Vector)s>& position,
                             const FieldList<%(Dim)s, %(SymTensor)s>& Hfield,
                             const TableKernel< %(Dim)s >& kernel);

}
"""

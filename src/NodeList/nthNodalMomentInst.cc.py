text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeList/nthNodalMoment.cc"

namespace Spheral {
  template
  FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Scalar> 
  nthNodalMoment<Dim<%(ndim)s>, 0U>(const DataBase<Dim<%(ndim)s>>&,
                                    const TableKernel<Dim<%(ndim)s>>&,
                                    const bool);
  template
  FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Vector> 
  nthNodalMoment<Dim<%(ndim)s>, 1U>(const DataBase<Dim<%(ndim)s>>&,
                                    const TableKernel<Dim<%(ndim)s>>&,
                                    const bool);
  template
  FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::SymTensor> 
  nthNodalMoment<Dim<%(ndim)s>, 2U>(const DataBase<Dim<%(ndim)s>>&,
                                    const TableKernel<Dim<%(ndim)s>>&,
                                    const bool);

  template
  void
  zerothAndFirstNodalMoments<Dim<%(ndim)s>>(const DataBase<Dim<%(ndim)s>>&,
                                            const TableKernel<Dim<%(ndim)s>>&,
                                            const bool,
                                            FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Scalar>&,
                                            FieldList<Dim<%(ndim)s>, Dim<%(ndim)s>::Vector>&);

}
"""

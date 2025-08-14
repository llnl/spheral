text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "RK/RKUtilities.cc"
#include "Geometry/Dimension.hh"
"""

for order in ["ZerothOrder", "LinearOrder", "QuadraticOrder"]:
    text += """
namespace Spheral {
template class RKUtilities<Dim<%(ndim)s>, """
    text += """RKOrder::%(order)s>;
}
""" % {"order" : order}

text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "RK/RKUtilities.cc"
#include "Geometry/Dimension.hh"
"""

for order in ["QuarticOrder"]:
    text += """
namespace Spheral {
template class RKUtilities<Dim<%(ndim)s>, """
    text += """RKOrder::%(order)s>;
}
""" % {"order" : order}

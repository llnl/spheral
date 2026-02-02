//---------------------------------Spheral++----------------------------------//
// Flag the nodes commonly connected within the given H radius as distinct
// fragments.  Result returned as an integer Field of fragement IDs.
//------------------------------------------------------------------------------
#ifndef __Spheral_identifyFragment__
#define __Spheral_identifyFragment__

namespace Spheral {

template<typename Dimension> class NodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
Field<Dimension, int>
computeFragmentField(const NodeList<Dimension>& nodeList,
                     const double linkRadius,
                     const Field<Dimension, typename Dimension::Scalar>& density,
                     const Field<Dimension, typename Dimension::SymTensor>& damage,
                     const Field<Dimension, int>& mask,
                     const double densityThreshold,
                     const double damageThreshold,
                     const bool assignDustToFragments);

}

#endif

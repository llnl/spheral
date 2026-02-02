//---------------------------------Spheral++----------------------------------//
// ASPHRadialFunctor
//
// Provides user-overridable hooks to modify how the ASPH object computes
// radial normals and magnitude
//
// Created by JMO, Thu Oct 10 14:12:37 PDT 2024
//----------------------------------------------------------------------------//
#ifndef __Spheral_ASPHRadialFunctor__
#define __Spheral_ASPHRadialFunctor__

namespace Spheral {

template<typename Dimension>
class ASPHRadialFunctor {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;

  // Constructors, destructor.
  ASPHRadialFunctor()           {}
  virtual ~ASPHRadialFunctor()  {}

  // Compute the outward pointing radial unit vector
  virtual Vector radialUnitVector(const size_t nodeListi,
                                  const size_t i,
                                  const Vector& posi) const { return posi.unitVector(); }

  // Compute the radial coordinate
  virtual Scalar radialCoordinate(const size_t nodeListi,
                                  const size_t i,
                                  const Vector& posi) const { return posi.magnitude(); }
};

}

#endif

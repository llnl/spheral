//---------------------------------Spheral++----------------------------------//
// Modified form of the standard SPH pair-wise viscosity due to Monaghan &
// Gingold.  This form is modified to use the velocity gradient to limit the
// velocity jump at the mid-point between points.
//
// Created by JMO, Thu Nov 20 14:13:18 PST 2014
//
// CURRENTLY UNUSABLE
//----------------------------------------------------------------------------//
#ifndef LimitedMonaghanGingoldViscosityRZ_HH
#define LimitedMonaghanGingoldViscosityRZ_HH

#include "LimitedMonaghanGingoldViscosity.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

class LimitedMonaghanGingoldViscosityRZ final : public LimitedMonaghanGingoldViscosity<Dim<2> > {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<2> Dimension;
  typedef Dimension::Scalar Scalar;
  typedef Dimension::Vector Vector;
  typedef Dimension::Tensor Tensor;
  typedef Dimension::SymTensor SymTensor;
  typedef Dimension::ThirdRankTensor ThirdRankTensor;
  typedef Dimension::FourthRankTensor FourthRankTensor;
  typedef Dimension::FifthRankTensor FifthRankTensor;
  typedef ArtificialViscosity<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  LimitedMonaghanGingoldViscosityRZ(const Scalar Clinear,
                                    const Scalar Cquadratic,
                                    const bool linearInExpansion,
                                    const bool quadraticInExpansion,
                                    const Scalar etaCritFrac,
                                    const Scalar etaFoldFrac);

  // Destructor.
  virtual ~LimitedMonaghanGingoldViscosityRZ();

  // The required method to compute the artificial viscous P/rho^2.
  virtual std::pair<Tensor, Tensor> Piij(const unsigned nodeListi, const unsigned i, 
                                         const unsigned nodeListj, const unsigned j,
                                         const Vector& xi,
                                         const Vector& etai,
                                         const Vector& vi,
                                         const Scalar rhoi,
                                         const Scalar csi,
                                         const SymTensor& Hi,
                                         const Vector& xj,
                                         const Vector& etaj,
                                         const Vector& vj,
                                         const Scalar rhoj,
                                         const Scalar csj,
                                         const SymTensor& Hj) const;

  // Restart methods.
  virtual std::string label() const { return "LimitedMonaghanGingoldViscosityRZ"; }

private:
  //--------------------------- Private Interface ---------------------------//
  LimitedMonaghanGingoldViscosityRZ();
  LimitedMonaghanGingoldViscosityRZ(const LimitedMonaghanGingoldViscosityRZ&);
  LimitedMonaghanGingoldViscosityRZ& operator=(const LimitedMonaghanGingoldViscosityRZ&) const;
};

}

#endif

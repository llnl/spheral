//---------------------------------Spheral++----------------------------------//
// A simple form for the artificial viscosity due to Monaghan & Gingold.
// References: 
//   Monaghan, J. J, & Gingold, R. A. 1983, J. Comput. Phys., 52, 374
//   Monaghan, J. J. 1992, ARA&A, 30, 543
//
// This version is adapted for cylindrical (RZ) coordinates based on the work
//   García-Senz, D., Relaño, A., Cabezón, R. M., & Bravo, E. (2009).
//   Axisymmetric smoothed particle hydrodynamics with self-gravity.
//   Monthly Notices of the Royal Astronomical Society, 392(1), 346–360.
//   http://doi.org/10.1111/j.1365-2966.2008.14058.x
// 
// Created by JMO, Mon Nov 20 15:50:29 PST 2006
//
// CURRENTLY UNUSABLE
//----------------------------------------------------------------------------//
#ifndef __Spheral_MonaghanGingoldViscosityGSRZ_hh__
#define __Spheral_MonaghanGingoldViscosityGSRZ_hh__

#include "MonaghanGingoldViscosity.hh"

namespace Spheral {

class MonaghanGingoldViscosityGSRZ: public MonaghanGingoldViscosity<Dim<2> > {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<2> Dimension;
  typedef Dimension::Scalar Scalar;
  typedef Dimension::Vector Vector;
  typedef Dimension::Tensor Tensor;
  typedef Dimension::SymTensor SymTensor;

  // Constructors.
  MonaghanGingoldViscosityGSRZ(const Scalar Clinear,
                               const Scalar Cquadratic,
                               const bool linearInExpansion,
                               const bool quadraticInExpansion);

  // Destructor.
  virtual ~MonaghanGingoldViscosityGSRZ();

  // The required method to compute the artificial viscous P/rho^2.
  virtual std::pair<Tensor, Tensor> Piij(const size_t nodeListi, const size_t i, 
                                         const size_t nodeListj, const size_t j,
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
  virtual std::string label() const { return "MonaghanGingoldViscosityGSRZ"; }

private:
  //--------------------------- Private Interface ---------------------------//
  MonaghanGingoldViscosityGSRZ();
};

}

#endif

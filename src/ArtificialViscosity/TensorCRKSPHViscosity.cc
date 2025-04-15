//---------------------------------Spheral++----------------------------------//
// A specialized form of the TensorMonaghanGingoldViscosity for use with CRKSPH.
//
// Created by J. Michael Owen, Wed Nov  5 23:51:31 PST 2014
//----------------------------------------------------------------------------//
#include "TensorCRKSPHViscosity.hh"
#include "FileIO/FileIO.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Boundary/Boundary.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/FluidNodeList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/rotationMatrix.hh"
#include "Utilities/GeometricUtilities.hh"
#include "RK/RKFieldNames.hh"
#include "RK/gradientRK.hh"
#include "Utilities/Timer.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
TensorCRKSPHViscosity<Dimension>::
TensorCRKSPHViscosity(const Scalar Clinear,
                      const Scalar Cquadratic,
                      const TableKernel<Dimension>& kernel,
                      const RKOrder order):
  TensorMonaghanGingoldViscosity<Dimension>(Clinear, Cquadratic, kernel),
  mOrder(order) {
}

//------------------------------------------------------------------------------
// Override the base method of computing the velocity gradient
//------------------------------------------------------------------------------
template<typename Dimension>
void
TensorCRKSPHViscosity<Dimension>::
updateVelocityGradient(const DataBase<Dimension>& dataBase,
                       const State<Dimension>& state,
                       const StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("TensorCRKSPHViscosity_updateVelocityGradient");

  // Get the necessary state.
  const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto  velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto  rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto  WR = state.template get<ReproducingKernel<Dimension>>(RKFieldNames::reproducingKernel(mOrder));
  const auto  corrections = state.fields(RKFieldNames::rkCorrections(mOrder), RKCoefficients<Dimension>());
  const auto& connectivityMap = dataBase.connectivityMap();

  auto DvDx_AV = state.fields(HydroFieldNames::ArtificialViscosityVelocityGradient, Tensor::zero);

  // Compute the basic velocity gradient.
  const auto vol = mass/rho;
  DvDx_AV.assignFields(gradientRK(velocity, position, vol, H, connectivityMap, WR, corrections,  NodeCoupling()));
  for (auto* fptr: DvDx_AV) fptr->name(HydroFieldNames::ArtificialViscosityVelocityGradient);

  TIME_END("TensorCRKSPHViscosity_updateVelocityGradient");
}

}

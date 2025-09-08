//---------------------------------Spheral++----------------------------------//
// PalphaPorosity
// 
// See header for references and such.
//----------------------------------------------------------------------------//
#include "Porosity/PalphaPorosity.hh"
#include "Material/EquationOfState.hh"
#include "FileIO/FileIO.hh"
#include "Field/Field.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "Utilities/globalNodeIDs.hh"

#include <string>

using std::vector;
using std::string;
using std::to_string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PalphaPorosity<Dimension>::
PalphaPorosity(const SolidNodeList<Dimension>& nodeList,
               const double phi0,
               const double Pe,
               const double Pt,
               const double Ps,
               const double alphae,
               const double alphat,
               const double n1,
               const double n2,
               const double cS0,
               const double c0,
               const double rhoS0,
               const bool jutziStateUpdate):
  PorosityModel<Dimension>(nodeList, phi0, cS0, c0, rhoS0, jutziStateUpdate),
  mPe(Pe),
  mPt(Pt),
  mPs(Ps),
  mAlphae(alphae),
  mAlphat(alphat),
  mn1(n1),
  mn2(n2),
  mdPdU(HydroFieldNames::partialPpartialEps, nodeList),
  mdPdR(HydroFieldNames::partialPpartialRho, nodeList) {
  VERIFY2((mPe <= mPt) and (mPt <= mPs),
          "PalphaPorosity input ERROR : require Pe <= Pt <= Ps: (Pe, Pt, Ps) = " << mPe << ", Pt = " << mPt << ", " << mPs);
  const auto nglobal = numGlobalNodes(nodeList);
  if (nglobal > 0) {
    const auto alpha0_max = mAlpha0.max();
    VERIFY2((1.0 <= mAlphae) and (mAlphat <= mAlphae) and (mAlphae <= alpha0_max),
            "PalphaPorosity input ERROR : require 1.0 <= alphat <= alphae <= alpha0, (alphat, alphae, alpha0) = " << mAlphat << ", " << mAlphae << ", " << alpha0_max);
  }
}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PalphaPorosity<Dimension>::
PalphaPorosity(const SolidNodeList<Dimension>& nodeList,
               const Field<Dimension, Scalar>& phi0,
               const double Pe,
               const double Pt,
               const double Ps,
               const double alphae,
               const double alphat,
               const double n1,
               const double n2,
               const double cS0,
               const Field<Dimension, Scalar>& c0,
               const double rhoS0,
               const bool jutziStateUpdate):
  PorosityModel<Dimension>(nodeList, phi0, cS0, c0, rhoS0, jutziStateUpdate),
  mPe(Pe),
  mPt(Pt),
  mPs(Ps),
  mAlphae(alphae),
  mAlphat(alphat),
  mn1(n1),
  mn2(n2),
  mdPdU(HydroFieldNames::partialPpartialEps, nodeList),
  mdPdR(HydroFieldNames::partialPpartialRho, nodeList) {
  VERIFY2((mPe <= mPt) and (mPt <= mPs),
          "PalphaPorosity input ERROR : require Pe <= Pt <= Ps: (Pe, Pt, Ps) = " << mPe << ", Pt = " << mPt << ", " << mPs);
  const auto nglobal = numGlobalNodes(nodeList);
  if (nglobal > 0) {
    const auto alpha0_max = mAlpha0.max();
    VERIFY2((1.0 <= mAlphae) and (mAlphat <= mAlphae) and (mAlphae <= alpha0_max),
            "PalphaPorosity input ERROR : require 1.0 <= alphat <= alphae <= alpha0, (alphat, alphae, alpha0) = " << mAlphat << ", " << mAlphae << ", " << alpha0_max);
    const auto n = nodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      mc0[i] = c0[i];
    }
  }
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PalphaPorosity<Dimension>::
~PalphaPorosity() {
}

//------------------------------------------------------------------------------
// Evaluate derivatives.
// For this model we evaluate the derivative of the alpha field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PalphaPorosity<Dimension>::
evaluateDerivatives(const Scalar time,
                    const Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  // Get the state fields.
  const auto  buildKey = [&](const std::string& fkey) { return StateBase<Dimension>::buildFieldKey(fkey, mNodeList.name()); };
  const auto& rho = state.field(buildKey(HydroFieldNames::massDensity), 0.0);
  const auto& P = state.field(buildKey(HydroFieldNames::pressure), 0.0);
  const auto& dPdu = state.field(buildKey(HydroFieldNames::partialPpartialEps), 0.0);
  const auto& dPdr = state.field(buildKey(HydroFieldNames::partialPpartialRho), 0.0);
  const auto& alpha = state.field(buildKey(SolidFieldNames::porosityAlpha), 0.0);
  const auto& DrhoDt = derivs.field(buildKey(IncrementBoundedState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity), 0.0);
  const auto& DuDt = derivs.field(buildKey(IncrementBoundedState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy), 0.0);
  auto&       DalphaDt = derivs.field(buildKey(IncrementBoundedState<Dimension, Scalar>::prefix() + SolidFieldNames::porosityAlpha), 0.0);
  auto&       fDSnew = derivs.field(buildKey(ReplaceBoundedState<Dimension, Scalar>::prefix() + SolidFieldNames::fDSjutzi), 0.0);

  // Walk the nodes.
  const auto n = mNodeList.numInternalNodes();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    const auto rhoi = rho(i);
    const auto Pi = P(i);
    const auto dPsdui = dPdu(i);
    const auto dPsdri = dPdr(i);
    const auto alphai = alpha(i);
    const auto c0i = mc0(i);         // Initial sound speed with porosity
    const auto DrhoDti = DrhoDt(i);
    const auto DuDti = DuDt(i);

    // If fully compressed we're done
    if (alphai == 1.0) {

      DalphaDt(i) = 0.0;
      fDSnew(i) = 1.0;

    } else {

      // Time evolution of the solid pressure
      const auto dPsdti = dPsdri*DrhoDti + dPsdui*DuDti;

      // If we're above the solid transition pressure we should compress out all porosity
      if (Pi >= mPs) {

        DalphaDt(i) = (1.0 - alphai)*safeInvVar(dt);
        fDSnew(i) = 1.0;

      } else {

        // First compute the derivative with respect to pressure
        auto DalphaDpi = 0.0;
        if (alphai > 1.0) {
          if (Pi < mPe or dPsdti < 0.0) {

            // Elastic
            if (c0i != mcS0) {  // If initial porous sound speed is the same as solid phase, no elastic evolution
              const auto halpha = 1.0 + (alphai - 1.0)*(c0i - mcS0)*safeInvVar(mcS0*(mAlphae - 1.0));
              DalphaDpi = alphai*alphai/mKS0*(1.0 - safeInvVar(halpha*halpha));
            }

          } else {

            // Plastic
            DalphaDpi = (Pi < mPt ?
                         mn1*(mAlphat - mAlphae)*pow((mPt - Pi)/(mPt - mPe), mn1)*safeInv(mPt - Pi) + mn2*(1.0 - mAlphat)*pow((mPs - Pi)/(mPs - mPe), mn2)*safeInv(mPs - Pi) :
                         mn2*(1.0 - mAlphat)*pow((mPs - Pi)/(mPs - mPe), mn2)*safeInv(mPs - Pi));

          }
        }

        // Now we can compute the final time derivative
        DalphaDpi = min(0.0, DalphaDpi);  // Keep things physical
        const auto Ainv = safeInvVar(alphai + DalphaDpi*(Pi - rhoi*dPsdri));
        const auto dPdti = (alphai*dPsdri*DrhoDti + dPsdui*DuDti)*Ainv;
        DalphaDt(i) = DalphaDpi*dPdti;

        // Optionally update the deviatoric stress scaling term
        if (mJutziStateUpdate) {
          auto DalphaDrhoi = (Pi/(rhoi*rhoi)*dPsdui + alphai*dPsdri)*Ainv * DalphaDpi;
          fDSnew(i) = std::max(0.0, std::min(1.0, 1.0 + DalphaDrhoi*rhoi/alphai));
        } else {
          fDSnew(i) = 1.0;
        }
      }
    }

#pragma omp critical
    {
      mMaxAbsDalphaDt = max(mMaxAbsDalphaDt, abs(DalphaDt(i)));    // For use in the timestep calculation
    }
  }
}

//------------------------------------------------------------------------------
// Register our state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PalphaPorosity<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  // PorosityModel does a lot of the work
  PorosityModel<Dimension>::registerState(dataBase, state);

  // We need the pressure derivatives
  state.enroll(mdPdU);
  state.enroll(mdPdR);
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PalphaPorosity<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  PorosityModel<Dimension>::dumpState(file, pathName);
  file.write(mdPdU, pathName + "/dPdU");
  file.write(mdPdR, pathName + "/dPdR");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PalphaPorosity<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  PorosityModel<Dimension>::restoreState(file, pathName);
  file.read(mdPdU, pathName + "/dPdU");
  file.read(mdPdR, pathName + "/dPdR");
}

}

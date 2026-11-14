//---------------------------------Spheral++----------------------------------//
// AxisBoundary -- A specialized Reflecting boundary condition for use along
// the axis in RZ calculations.
// Note: this boundary is automatically constructred by the RZ hydro objects,
//       so the user should not explicitly add this boundary.
//
// Created by JMO, Sun Aug 14 10:20:53 PDT 2016
//----------------------------------------------------------------------------//
#ifndef __Spheral_AxisBoundary__
#define __Spheral_AxisBoundary__

#include "Geometry/Dimension.hh"
#include "ReflectingBoundary.hh"

namespace Spheral {

class AxisBoundaryRZ: public ReflectingBoundary<Dim<2> > {

public:
  //--------------------------- Public Interface ---------------------------//
  using Dimension = Dim<2>;
  using Scalar = Dimension::Scalar;
  using Vector = Dimension::Vector;
  using Tensor = Dimension::Tensor;
  using SymTensor = Dimension::SymTensor;
  using ThirdRankTensor = Dimension::ThirdRankTensor;

  // Constructors and destructors.
  AxisBoundaryRZ(const DataBase<Dimension>& dataBase,
                 const Scalar etamin);
  virtual ~AxisBoundaryRZ() = default;

  // Find the set of nodes in violation of this boundary in the given NodeList.
  virtual void setViolationNodes(NodeList<Dimension>& nodeList) override;
  virtual void updateViolationNodes(NodeList<Dimension>& nodeList) override;

  // Access the fuzz from the axis we're using to enforce the BC.
  Scalar etamin() const;
  void etamin(const Scalar x);

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "AxisBoundaryRZ"; }
  //****************************************************************************

  // Prevent the Boundary virtual methods from being hidden
  using Boundary<Dimension>::applyGhostBoundary;
  using Boundary<Dimension>::enforceBoundary;

private:
  //--------------------------- Private Interface ---------------------------//
  double mEtaMin;
};

}

#endif

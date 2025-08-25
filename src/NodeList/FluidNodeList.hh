//---------------------------------Spheral++----------------------------------//
// FluidNodeList -- An abstract base class for the NodeLists to represent
//                  fluids.
//
// Created by JMO, Sun Sep 12 17:12:24 PDT 1999
//----------------------------------------------------------------------------//
#ifndef __Spheral_FluidNodeList__
#define __Spheral_FluidNodeList__

#include "NodeList.hh"

#include <float.h>
#include <string>

namespace Spheral {

template<typename Dimension> class NodeIteratorBase;
template<typename Dimension> class AllNodeIterator;
template<typename Dimension> class InternalNodeIterator;
template<typename Dimension> class GhostNodeIterator;
template<typename Dimension> class MasterNodeIterator;
template<typename Dimension> class CoarseNodeIterator;
template<typename Dimension> class RefineNodeIterator;
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class Neighbor;
template<typename Dimension> class ConnectivityMap;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class EquationOfState;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;
class FileIO;

template<typename Dimension>
class FluidNodeList: public NodeList<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  
  // Constructors.
  FluidNodeList(std::string name,
                EquationOfState<Dimension>& eos,
                const size_t numInternal,
                const size_t numGhost,
                const Scalar hmin,
                const Scalar hmax,
                const Scalar hminratio,
                const Scalar nPerh,
                const size_t maxNumNeighbors,
                const Scalar rhoMin,
                const Scalar rhoMax);

  // Destructor
  virtual ~FluidNodeList() = default;

  // Access the fluid state variables.
  Field<Dimension, Scalar>& massDensity();
  Field<Dimension, Scalar>& specificThermalEnergy();

  const Field<Dimension, Scalar>& massDensity() const;
  const Field<Dimension, Scalar>& specificThermalEnergy() const;

  void massDensity(const Field<Dimension, Scalar>& rho);
  void specificThermalEnergy(const Field<Dimension, Scalar>& eps);

  // These are quantities which are not stored, but can be computed.
  virtual void pressure(Field<Dimension, Scalar>& field) const;
  virtual void temperature(Field<Dimension, Scalar>& field) const;
  virtual void soundSpeed(Field<Dimension, Scalar>& field) const;
  virtual void volume(Field<Dimension, Scalar>& field) const;
  virtual void linearMomentum(Field<Dimension, Vector>& field) const;
  virtual void totalEnergy(Field<Dimension, Scalar>& field) const;

  // Access the equation of state.
  const EquationOfState<Dimension>& equationOfState() const;
  void equationOfState(const EquationOfState<Dimension>& eos);

  // Optional bounding mass densities for use when time integrating
  // the density.
  Scalar rhoMin() const;
  Scalar rhoMax() const;

  void rhoMin(Scalar x);
  void rhoMax(Scalar x);

  //****************************************************************************
  // Methods required for restarting.
  // Dump and restore the NodeList state.
  virtual std::string label() const { return "FluidNodeList"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

  // No default constructor or copying.
  FluidNodeList() = delete;
  FluidNodeList(const FluidNodeList& nodes) = delete;
  FluidNodeList& operator=(const FluidNodeList& rhs) = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  // Min/max mass densities.
  Scalar mRhoMin, mRhoMax;

  // Fields that define the fluid's current state.
  Field<Dimension, Scalar> mMassDensity;
  Field<Dimension, Scalar> mSpecificThermalEnergy;

  // Equation of state.
  const EquationOfState<Dimension>* mEosPtr;
};

}

#include "FluidNodeListInline.hh"

#endif


//---------------------------------Spheral++----------------------------------//
// SubPointPressureHourglassControl
//
// Impose additional forces on each point using subdivisions of the Voronoi
// control volume to constrain the unphysical degrees of freedom in our hydro
// discretization and avoid spurious so-called houglass modes.
//----------------------------------------------------------------------------//
#ifndef __Spheral_SubPointPressureHourglassControl__
#define __Spheral_SubPointPressureHourglassControl__

#include "Field/FieldList.hh"
#include "Physics/Physics.hh"
#include "DataOutput/registerWithRestart.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;
template<typename Dimension> class Boundary;

template<typename Dimension>
class SubPointPressureHourglassControl : public Physics<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using FacetedVolume = typename Dimension::FacetedVolume;
  
  using BoundaryIterator = typename std::vector<Boundary<Dimension>*>::iterator;
  using ConstBoundaryIterator = typename std::vector<Boundary<Dimension>*>::const_iterator;
  using TimeStepType = typename std::pair<double, std::string>;

  // Constructor
  SubPointPressureHourglassControl(const Scalar fHG,
                                   const Scalar xfilter);

  // Destructor.
  virtual ~SubPointPressureHourglassControl();

  //******************************************************************************//
  // Stuff to do once on problem startup
  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  // Evaluate derivatives
  virtual void evaluateDerivatives(const Scalar time,
                                   const Scalar dt,
                                   const DataBase<Dimension>& dataBase,
                                   const State<Dimension>& state,
                                   StateDerivatives<Dimension>& derivatives) const override;
  
  // Vote on a time step.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar time) const override;

  // Register the state
  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;

  // Register the state derivatives
  virtual void registerDerivatives(DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs) override;

  // Label
  virtual std::string label() const override { return "SubPointPressureHourglassControl"; }

  // Does this package require Voronoi-like cells per point?
  virtual bool requireVoronoiCells() const override { return true; }

  // Access parameters
  Scalar fHG() const                                   { return mfHG; }
  Scalar xfilter() const                               { return mxfilter; }
  void fHG(const Scalar x)                             { mfHG = x; }
  void xfilter(const Scalar x)                         { mxfilter = x; }
  const FieldList<Dimension, Vector>&    DvDt() const  { return mDvDt; }

  //****************************************************************************
  // Methods required for restarting.
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

  // No default constructor, copying, or assignment.
  SubPointPressureHourglassControl() = delete;
  SubPointPressureHourglassControl(const SubPointPressureHourglassControl&) = delete;
  SubPointPressureHourglassControl& operator=(const SubPointPressureHourglassControl&) = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  Scalar mfHG, mxfilter;
  mutable FieldList<Dimension, Vector> mDvDt;

  // The restart registration.
  RestartRegistrationType mRestart;
};

}

#endif

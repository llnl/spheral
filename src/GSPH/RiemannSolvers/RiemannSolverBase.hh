//---------------------------------Spheral++----------------------------------//
// RiemannSolverBase -- base class for our riemann solvers 
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//
#ifndef __Spheral_RiemannSolverBase_hh__
#define __Spheral_RiemannSolverBase_hh__

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;
template<typename Dimension> class WaveSpeedBase;
template<typename Dimension> class LimiterBase;
template<typename Dimension> class Boundary;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
class RiemannSolverBase {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename std::vector<Boundary<Dimension>*>::const_iterator ConstBoundaryIterator;

public:

  RiemannSolverBase(LimiterBase<Dimension>& slopeLimiter,
                    WaveSpeedBase<Dimension>& waveSpeedBase,
                    const bool linearReconstruction);

  ~RiemannSolverBase();

  virtual 
  void initialize(const DataBase<Dimension>& dataBase,
                  const State<Dimension>& state,
                  const StateDerivatives<Dimension>& derivs,
                  ConstBoundaryIterator boundaryBegin,
                  ConstBoundaryIterator boundaryEnd,
                  const Scalar time,
                  const Scalar dt,
                  const TableKernel<Dimension>& W);


    virtual
  void interfaceState(const Vector& ri,
                      const Vector& rj,
                      const SymTensor& Hi,
                      const SymTensor& Hj,
                      const Scalar& rhoi,   
                      const Scalar& rhoj,
                      const Scalar& epsi,
                      const Scalar& epsj, 
                      const Scalar& ci,   
                      const Scalar& cj, 
                      const Scalar& Pi,    
                      const Scalar& Pj,
                      const Vector& vi,    
                      const Vector& vj,
                      const Vector& DrhoDxi,
                      const Vector& DrhoDxj,
                      const Vector& DepsDxi,
                      const Vector& DepsDxj,
                      const Vector& DpDxi,    
                      const Vector& DpDxj,
                      const Tensor& DvDxi,    
                      const Tensor& DvDxj,
                            Scalar& Pstar,
                            Vector& vstar,
                            Scalar& rhostari,
                            Scalar& rhostarj,
                            Scalar& epsstari,
                            Scalar& epsstarj) const;


  virtual
  void interfaceState(const Vector& ri,
                      const Vector& rj,
                      const SymTensor& Hi,
                      const SymTensor& Hj,
                      const Scalar& rhoi,   
                      const Scalar& rhoj, 
                      const Scalar& ci,   
                      const Scalar& cj,
                      const Scalar& Pi,    
                      const Scalar& Pj,
                      const Vector& vi,    
                      const Vector& vj,
                      const SymTensor& Si,    
                      const SymTensor& Sj,
                      const Tensor& Di,    
                      const Tensor& Dj,
                            Vector& Tstar,
                            Vector& vstar) const;

  LimiterBase<Dimension>& limiter() const;
  WaveSpeedBase<Dimension>& waveSpeed() const;

  bool linearReconstruction() const;
  void linearReconstruction(bool x);

  virtual 
  void linearReconstruction(const Vector& ri,
                            const Vector& rj,
                            const Scalar& yi,
                            const Scalar& yj,
                            const Vector& DyDxi,
                            const Vector& DyDxj,
                                  Scalar& ytildei,
                                  Scalar& ytildej) const;

  virtual 
  void linearReconstruction(const Vector& ri,
                            const Vector& rj,
                            const Vector& yi,
                            const Vector& yj,
                            const Tensor& DyDxi,
                            const Tensor& DyDxj,
                                  Vector& ytildei,
                                  Vector& ytildej) const;



private:
  
  LimiterBase<Dimension>& mSlopeLimiter;   
  WaveSpeedBase<Dimension>& mWaveSpeed;
  bool mLinearReconstruction;

};

}

#include "RiemannSolverBaseInline.hh"

#endif


//
//  HelmholtzEquationOfState.h
//  
//
//  Created by Raskin, Cody Dantes on 8/28/14.
//
//

#ifndef ____HelmholtzEquationOfState_hh__
#define ____HelmholtzEquationOfState_hh__

#include "EquationOfState.hh"
#include "Field/FieldList.hh"

namespace Spheral {
        
template<typename Dimension>
class HelmholtzEquationOfState: public EquationOfState<Dimension> {
      
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
          
  // Constructors, destructors.
  HelmholtzEquationOfState(const PhysicalConstants& constants,
                           const double minimumPressure,
                           const double maximumPressure,
                           const double minimumTemperature,
                           const MaterialPressureMinType minPressureType,
                           const Scalar abar0,
                           const Scalar zbar0,
                           const double externalPressure);
  ~HelmholtzEquationOfState();

          
  // We require any equation of state to define the following properties.
  virtual void setPressure(Field<Dimension, Scalar>& Pressure,
                           const Field<Dimension, Scalar>& massDensity,
                           const Field<Dimension, Scalar>& specificThermalEnergy) const override;
          
  virtual void setPressureAndDerivs(Field<Dimension, Scalar>& Pressure,           // set pressure
                                    Field<Dimension, Scalar>& dPdu,               // set (\partial P)/(\partial u) (specific thermal energy)
                                    Field<Dimension, Scalar>& dPdrho,             // set (\partial P)/(\partial rho) (density)
                                    const Field<Dimension, Scalar>& massDensity,
                                    const Field<Dimension, Scalar>& specificThermalEnergy) const override;

  virtual void setTemperature(Field<Dimension, Scalar>& temperature,
                              const Field<Dimension, Scalar>& massDensity,
                              const Field<Dimension, Scalar>& specificThermalEnergy) const override;
          
  virtual void setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                                        const Field<Dimension, Scalar>& massDensity,
                                        const Field<Dimension, Scalar>& temperature) const override;
          
  virtual void setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                               const Field<Dimension, Scalar>& massDensity,
                               const Field<Dimension, Scalar>& temperature) const override;
          
  virtual void setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
                             const Field<Dimension, Scalar>& massDensity,
                             const Field<Dimension, Scalar>& specificThermalEnergy) const override;
          
  virtual void setGammaField(Field<Dimension, Scalar>& gamma,
                             const Field<Dimension, Scalar>& massDensity,
                             const Field<Dimension, Scalar>& specificThermalEnergy) const override;
                      
  virtual void setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
                              const Field<Dimension, Scalar>& massDensity,
                              const Field<Dimension, Scalar>& specificThermalEnergy) const override;
          
  virtual void setEntropy(Field<Dimension, Scalar>& entropy,
                          const Field<Dimension, Scalar>& massDensity,
                          const Field<Dimension, Scalar>& specificThermalEnergy) const override;
          
  // Some of the following methods are disabled
  virtual Scalar pressure(const Scalar /*massDensity*/,
                          const Scalar /*specificThermalEnergy*/) const { VERIFY2(false, "HelmholtzEquationOfState does not support individual state calls."); return 0; }
          
  virtual Scalar temperature(const Scalar /*massDensity*/,
                             const Scalar /*specificThermalEnergy*/) const { VERIFY2(false, "HelmholtzEquationOfState does not support individual state calls."); return 0; }
          
  virtual Scalar specificThermalEnergy(const Scalar /*massDensity*/,
                                       const Scalar /*temperature*/) const { VERIFY2(false, "HelmholtzEquationOfState does not support individual state calls."); return 0; }
          
  virtual Scalar specificHeat(const Scalar /*massDensity*/,
                              const Scalar /*temperature*/) const { VERIFY2(false, "HelmholtzEquationOfState does not support individual state calls."); return 0; }
          
  virtual Scalar soundSpeed(const Scalar /*massDensity*/,
                            const Scalar /*specificThermalEnergy*/) const { VERIFY2(false, "HelmholtzEquationOfState does not support individual state calls."); return 0; }
          
  // Get the effective gamma (ratio of specific heats) for this eos.
  virtual Scalar gamma(const Scalar /*massDensity*/,
                       const Scalar /*specificThermalEnergy*/) const { VERIFY2(false, "HelmholtzEquationOfState does not support individual state calls."); return 0; }
          
  // Get the bulk modulus.
  virtual Scalar bulkModulus(const Scalar /*massDensity*/,
                             const Scalar /*specificThermalEnergy*/) const { VERIFY2(false, "HelmholtzEquationOfState does not support individual state calls."); return 0; }
          
  virtual Scalar entropy(const Scalar /*massDensity*/,
                         const Scalar /*specificThermalEnergy*/) const { VERIFY2(false, "HelmholtzEquationOfState does not support individual state calls."); return 0; }

          
  const Field<Dimension, Scalar>& abar() const;
  const Field<Dimension, Scalar>& zbar() const;
                      
  bool getUpdateStatus() const;
  void setUpdateStatus(bool bSet);
          
  virtual bool valid() const override;
          
private:
  //--------------------------- Private Interface ---------------------------//
          
  mutable std::shared_ptr<Field<Dimension, Scalar> > myAbar;
  mutable std::shared_ptr<Field<Dimension, Scalar> > myZbar;
  mutable std::shared_ptr<Field<Dimension, Scalar> > mySpecificThermalEnergy;
  mutable std::shared_ptr<Field<Dimension, Scalar> > myMassDensity;
  mutable std::shared_ptr<Field<Dimension, Scalar> > myTemperature;
  mutable std::shared_ptr<Field<Dimension, Scalar> > myPressure;
  mutable std::shared_ptr<Field<Dimension, Scalar> > mySoundSpeed;
  mutable std::shared_ptr<Field<Dimension, Scalar> > myGamma;
  mutable std::shared_ptr<Field<Dimension, Scalar> > myEntropy;
          
  Scalar mabar0, mzbar0, mPmin, mPmax;
  mutable Scalar mTmin;
  bool needUpdate;
                      
  const PhysicalConstants& mConstants;
  Scalar mDistincm, mMassing, mEnergyinergpg, mTimeins, mPressureinbarye, mDensingpccm, mVelincmps;
          
  void storeFields(const Field<Dimension, Scalar>& thisMassDensity, const Field<Dimension, Scalar>& thisSpecificThermalEnergy) const;

};

}

#endif

//---------------------------------Spheral++----------------------------------//
// SolidNodeList -- A form of the FluidNodeList appropriate for use with 
// solid materials.
//
// Created by JMO, Tue Sep 7 22:44:37 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_SolidNodeList_hh__
#define __Spheral_SolidNodeList_hh__

#include "NodeList/FluidNodeList.hh"

namespace Spheral {

template<typename Dimension> class StrengthModel;

template<typename Dimension>
class SolidNodeList: public FluidNodeList<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors.
  SolidNodeList(std::string name,
                EquationOfState<Dimension>& eos,
                StrengthModel<Dimension>& strength,
                const size_t numInternal,
                const size_t numGhost,
                const Scalar hmin,
                const Scalar hmax,
                const Scalar hminratio,
                const Scalar nPerh,
                const size_t maxNumNeighbors,
                const Scalar rhoMin,
                const Scalar rhoMax);

  // Destructor.
  virtual ~SolidNodeList() = default;

  // Override the base method for calculating the sound speed.
  virtual void soundSpeed(Field<Dimension, Scalar>& field) const;

  // Convenience methods to compute derived moduli and such.
  virtual void bulkModulus(Field<Dimension, Scalar>& field) const;
  virtual void shearModulus(Field<Dimension, Scalar>& field) const;
  virtual void yieldStrength(Field<Dimension, Scalar>& field) const;
  virtual void YoungsModulus(Field<Dimension, Scalar>& field,
                             const Field<Dimension, Scalar>& K,          // bulk modulus
                             const Field<Dimension, Scalar>& mu) const;  // shear modulus
  virtual void longitudinalSoundSpeed(Field<Dimension, Scalar>& field,
                                      const Field<Dimension, Scalar>& rho,        // mass density
                                      const Field<Dimension, Scalar>& K,          // bulk modulus
                                      const Field<Dimension, Scalar>& mu) const;  // shear modulus

  // Access the stored deviatoric stress and plastic strain.
  Field<Dimension, SymTensor>& deviatoricStress();
  Field<Dimension, Scalar>& plasticStrain();
  Field<Dimension, Scalar>& plasticStrainRate();

  const Field<Dimension, SymTensor>& deviatoricStress() const;
  const Field<Dimension, Scalar>& plasticStrain() const;
  const Field<Dimension, Scalar>& plasticStrainRate() const;

  // The tensor damage
  Field<Dimension, SymTensor>& damage();
  const Field<Dimension, SymTensor>& damage() const;

  // The fragment ID field.
  Field<Dimension, int>& fragmentIDs();
  const Field<Dimension, int>& fragmentIDs() const;

  // The particle type field.
  Field<Dimension, int>& particleTypes();
  const Field<Dimension, int>& particleTypes() const;

  // The strength model this solid is using.
  const StrengthModel<Dimension>& strengthModel() const;

  //****************************************************************************
  // Methods required for restarting.
  // Dump and restore the NodeList state.
  virtual std::string label() const { return "SolidNodeList"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

  // No default constructor or copying.
  SolidNodeList() = delete;
  SolidNodeList(const SolidNodeList& solidNodeList) = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  Field<Dimension, SymTensor> mDeviatoricStress;
  Field<Dimension, Scalar> mPlasticStrain;
  Field<Dimension, Scalar> mPlasticStrainRate;
  Field<Dimension, SymTensor> mDamage;
  Field<Dimension, int> mFragmentIDs;
  Field<Dimension, int> mParticleTypes;

  // Pointer to the associated strength object.
  StrengthModel<Dimension>& mStrength;
};

}

#include "SolidNodeListInline.hh"

#endif

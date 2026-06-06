//---------------------------------Spheral++----------------------------------//
// IntegrationCoefficient
//
// Represents an integral to be integrated by KernelIntegrator
//----------------------------------------------------------------------------//
#ifndef __Spheral_IntegrationCoefficient_hh__
#define __Spheral_IntegrationCoefficient_hh__

#include <memory>
#include <vector>
#include <unordered_map>
#include "BilinearIndex.hh"
#include "KernelIntegrationData.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// IntegrationCoefficient
//
// Lets us do integrals of arbitrary coefficients
//------------------------------------------------------------------------------
template<typename Dimension, typename CoefficientType>
class IntegrationCoefficient {
public:
  IntegrationCoefficient() { }
  virtual ~IntegrationCoefficient() { }
  virtual CoefficientType evaluateCoefficient(const KernelIntegrationData<Dimension>& kid) const = 0;
};
    
//------------------------------------------------------------------------------
// ConstantIntegrationCoefficient
//
// Returns one
//------------------------------------------------------------------------------
template<typename Dimension, typename CoefficientType>
class ConstantIntegrationCoefficient : public IntegrationCoefficient<Dimension, CoefficientType> {
public:
  ConstantIntegrationCoefficient();
  virtual ~ConstantIntegrationCoefficient() {}
  ConstantIntegrationCoefficient(CoefficientType coeff);
  virtual CoefficientType evaluateCoefficient(const KernelIntegrationData<Dimension>& kid) const override;
  virtual const CoefficientType& getData() const;
  virtual void setData(CoefficientType coeff);
private:
  bool mDataSet;
  CoefficientType mCoeff;
};

//------------------------------------------------------------------------------
// DefaultIntegrationCoefficient
//
// Returns one
//------------------------------------------------------------------------------
template<typename Dimension, typename CoefficientType>
class DefaultIntegrationCoefficient : public IntegrationCoefficient<Dimension, CoefficientType> {
public:
  DefaultIntegrationCoefficient() { }
  virtual ~DefaultIntegrationCoefficient() { }
  virtual CoefficientType evaluateCoefficient(const KernelIntegrationData<Dimension>& kid) const override;
};

//------------------------------------------------------------------------------
// FieldListIntegrationCoefficient
//
// Interpolates the FieldList at the integration point
// f\left(x\right)=\sum_{i}V_{i}f_{i}W_{i}\left(x\right)
//------------------------------------------------------------------------------
template<typename Dimension, typename CoefficientType>
class FieldListIntegrationCoefficient : public IntegrationCoefficient<Dimension, CoefficientType> {
public:
  FieldListIntegrationCoefficient();
  virtual ~FieldListIntegrationCoefficient() {};
  FieldListIntegrationCoefficient(const FieldList<Dimension, CoefficientType>& data);
  virtual const FieldList<Dimension, CoefficientType>& getData() const;
  virtual void setData(const FieldList<Dimension, CoefficientType>& data);
  virtual void setDataPoint(int nodeListi, int nodei, const CoefficientType& d);
  virtual CoefficientType evaluateCoefficient(const KernelIntegrationData<Dimension>& kid) const override;
  
private:
  bool mDataSet;
  FieldList<Dimension, CoefficientType> mData;
};

//------------------------------------------------------------------------------
// IntegralDependsOnCoefficient
//
// Coefficient defaults to DefaultIntegrationCoefficient
//------------------------------------------------------------------------------
template<typename Dimension, typename CoefficientType>
class IntegralDependsOnCoefficient {
public:
  IntegralDependsOnCoefficient()  {
    mCoefficient = std::make_shared<DefaultIntegrationCoefficient<Dimension, CoefficientType>>();
  }

  virtual ~IntegralDependsOnCoefficient() { }

  // Give the coefficient to the integal
  virtual void setCoefficient(std::shared_ptr<IntegrationCoefficient<Dimension, CoefficientType>> coeff) {
    mCoefficient = coeff;
  }

  // Return a value for the coefficient
  virtual std::shared_ptr<IntegrationCoefficient<Dimension, CoefficientType>>
  getCoefficient() const {
    return mCoefficient;
  }
  
protected:
  std::shared_ptr<IntegrationCoefficient<Dimension, CoefficientType>> mCoefficient;
};

//------------------------------------------------------------------------------
// IntegralDependsOnFieldListCoefficient
//
// Coefficient defaults to DefaultIntegrationCoefficient
//------------------------------------------------------------------------------
template<typename Dimension, typename CoefficientType>
class IntegralDependsOnFieldListCoefficient {
public:
  IntegralDependsOnFieldListCoefficient()  {
    mCoefficient = std::make_shared<FieldListIntegrationCoefficient<Dimension, CoefficientType>>();
  }

  virtual ~IntegralDependsOnFieldListCoefficient() { }
  
  // Give the coefficient to the integal
  virtual void setCoefficient(std::shared_ptr<FieldListIntegrationCoefficient<Dimension, CoefficientType>> coeff) {
    mCoefficient = coeff;
  }

  // Return a value for the coefficient
  virtual std::shared_ptr<FieldListIntegrationCoefficient<Dimension, CoefficientType>>
  getCoefficient() const {
    return mCoefficient;
  }
  
protected:
  std::shared_ptr<FieldListIntegrationCoefficient<Dimension, CoefficientType>> mCoefficient;
};

} // end namespace Spheral

#include "IntegrationCoefficientInline.hh"

#endif

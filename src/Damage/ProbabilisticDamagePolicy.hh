//---------------------------------Spheral++----------------------------------//
// ProbabilisticDamagePolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent scalar damage state.
//
// References:
//   Benz, W. & Asphaug, E., 1995 "Computer Physics Comm.", 87, 253-265.
//   Benz, W. & Asphaug, E., 1994 "Icarus", 107, 98-116.
//   Randles, P.W. & Libersky, L.D., 1996, "Comput. Methods Appl. Engrg, 
//     139, 375-408
//
// Created by JMO, Sun Sep 26 16:15:17 PDT 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_ProbabilisticDamagePolicy_hh__
#define __Spheral_ProbabilisticDamagePolicy_hh__

#include "DataBase/FieldUpdatePolicy.hh"

#include <string>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class FluidNodeList;
template<typename Dimension, typename DataType> class Field;

template<typename Dimension>
class ProbabilisticDamagePolicy: public FieldUpdatePolicy<Dimension, typename Dimension::SymTensor> {
public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using KeyType = typename FieldUpdatePolicy<Dimension, SymTensor>::KeyType;

  // Constructors, destructor.
  ProbabilisticDamagePolicy(const bool damageInCompression,  // allow damage in compression
                            const double kWeibull,           // coefficient in Weibull power-law
                            const double mWeibull);          // exponenent in Weibull power-law
  virtual ~ProbabilisticDamagePolicy() = default;
  
  // Overload the methods describing how to update Fields.
  virtual void update(const KeyType& key,
                      State<Dimension>& state,
                      StateDerivatives<Dimension>& derivs,
                      const double multiplier,
                      const double t,
                      const double dt);

  // Equivalence.
  virtual bool operator==(const UpdatePolicyBase<Dimension>& rhs) const;

  static const std::string prefix() { return "delta "; }

  // Forbidden methods
  ProbabilisticDamagePolicy() = delete;
  ProbabilisticDamagePolicy(const ProbabilisticDamagePolicy& rhs) = delete;
  ProbabilisticDamagePolicy& operator=(const ProbabilisticDamagePolicy& rhs) = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  bool mDamageInCompression;
  double mkWeibull, mmWeibull;
};

}

#endif

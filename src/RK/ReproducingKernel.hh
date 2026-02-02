//---------------------------------Spheral++----------------------------------//
// ReproducingKernel
//
// Front-end class to wrap the Reproducing Kernel corrected methods
//----------------------------------------------------------------------------//
#ifndef __Spheral_ReproducingKernel__
#define __Spheral_ReproducingKernel__

#include "RK/ReproducingKernelMethods.hh"

namespace Spheral {

template<typename Dimension>
class ReproducingKernel: public ReproducingKernelMethods<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructor
  ReproducingKernel(const TableKernel<Dimension>& W,
                    const RKOrder order);
  ReproducingKernel();
  ReproducingKernel(const ReproducingKernel& rhs);
  ReproducingKernel& operator=(const ReproducingKernel& rhs);
  virtual ~ReproducingKernel() {}
  bool operator==(const ReproducingKernel& rhs) const;

  // Base kernel calls
  Scalar evaluateBaseKernel(const Vector& x,
                            const SymTensor& H) const;
  Vector evaluateBaseGradient(const Vector& x,
                              const SymTensor& H) const;
  SymTensor evaluateBaseHessian(const Vector& x,
                                const SymTensor& H) const;
  std::pair<Scalar, Vector> evaluateBaseKernelAndGradient(const Vector& x,
                                                          const SymTensor& H) const;

  // Corrected kernel calls
  Scalar evaluateKernel(const Vector& x,
                        const SymTensor& H,
                        const RKCoefficients<Dimension>& corrections) const;
  Vector evaluateGradient(const Vector& x,
                          const SymTensor& H,
                          const RKCoefficients<Dimension>& corrections) const;
  SymTensor evaluateHessian(const Vector& x,
                            const SymTensor& H,
                            const RKCoefficients<Dimension>& corrections) const;
  std::pair<Scalar, Vector> evaluateKernelAndGradient(const Vector& x,
                                                      const SymTensor& H,
                                                      const RKCoefficients<Dimension>& corrections) const;
  std::tuple<Scalar, Vector, Scalar> evaluateKernelAndGradients(const Vector& x,
                                                                const SymTensor& H,
                                                                const RKCoefficients<Dimension>& corrections) const;

  // Compute corrections and normals
  void computeCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                          const FieldList<Dimension, Scalar>& volume,
                          const FieldList<Dimension, Vector>& position,
                          const FieldList<Dimension, SymTensor>& H,
                          const bool needHessian,
                          FieldList<Dimension, RKCoefficients<Dimension>>& zerothCorrections,
                          FieldList<Dimension, RKCoefficients<Dimension>>& corrections) const;
  void computeNormal(const ConnectivityMap<Dimension>& connectivityMap,
                     const FieldList<Dimension, Scalar>& volume,
                     const FieldList<Dimension, Vector>& position,
                     const FieldList<Dimension, SymTensor>& H,
                     const FieldList<Dimension, RKCoefficients<Dimension>>& corrections,
                     FieldList<Dimension, Scalar>& surfaceArea,
                     FieldList<Dimension, Vector>& normal) const;

  // Access the internal state
  const TableKernel<Dimension>& kernel() const;

private:
  //--------------------------- Private Interface ---------------------------//
  const TableKernel<Dimension>* mWptr;   // The base interpolation kernel
};

}

#include "RK/ReproducingKernelInline.hh"

#endif

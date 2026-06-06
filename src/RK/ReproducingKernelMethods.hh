//---------------------------------Spheral++----------------------------------//
// ReproducingKernelMethods
//
// Front-end class to wrap the Reproducing Kernel corrected methods
//----------------------------------------------------------------------------//
#ifndef __Spheral_ReproducingKernelMethods__
#define __Spheral_ReproducingKernelMethods__

#include "RK/RKCoefficients.hh"
#include "RK/RKUtilities.hh"

namespace Spheral {

template<typename Dimension>
class ReproducingKernelMethods {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Eigen::SparseMatrix<double> TransformationMatrix;

  // Constructors, destructor
  ReproducingKernelMethods(const RKOrder order);
  ReproducingKernelMethods();
  ReproducingKernelMethods(const ReproducingKernelMethods& rhs);
  ReproducingKernelMethods& operator=(const ReproducingKernelMethods& rhs);
  virtual ~ReproducingKernelMethods() {}
  bool operator==(const ReproducingKernelMethods& rhs) const;

  // Build a transformation operator
  TransformationMatrix transformationMatrix(const Tensor& T,
                                            const bool needHessian) const;

  // Apply a transformation operator to a corrections vector
  void applyTransformation(const TransformationMatrix& T,
                           RKCoefficients<Dimension>& corrections) const;

  // Access the internal state
  RKOrder order() const;
  int gradCorrectionsSize() const;
  unsigned int hessCorrectionsSize() const;

protected:
  //--------------------------- Protected Interface ---------------------------//
  RKOrder mOrder;
  int mGradCorrectionsSize, mHessCorrectionsSize;

  // Pointers to the correct member methods of RKUtilities
  Scalar                    (*mEvaluateBaseKernel)(const TableKernel<Dimension>&, const Vector&, const SymTensor&);
  Vector                    (*mEvaluateBaseGradient)(const TableKernel<Dimension>&, const Vector&, const SymTensor&);
  SymTensor                 (*mEvaluateBaseHessian)(const TableKernel<Dimension>&, const Vector&, const SymTensor&);
  std::pair<Scalar, Vector> (*mEvaluateBaseKernelAndGradient)(const TableKernel<Dimension>&, const Vector&, const SymTensor&);

  Scalar                             (*mEvaluateKernel)(const TableKernel<Dimension>&, const Vector&, const SymTensor&, const RKCoefficients<Dimension>&);
  Vector                             (*mEvaluateGradient)(const TableKernel<Dimension>&, const Vector&, const SymTensor&, const RKCoefficients<Dimension>&);
  SymTensor                          (*mEvaluateHessian)(const TableKernel<Dimension>&, const Vector&, const SymTensor&, const RKCoefficients<Dimension>&);
  std::pair<Scalar, Vector>          (*mEvaluateKernelAndGradient)(const TableKernel<Dimension>&, const Vector&, const SymTensor&, const RKCoefficients<Dimension>&);
  std::tuple<Scalar, Vector, Scalar> (*mEvaluateKernelAndGradients)(const TableKernel<Dimension>&, const Vector&, const SymTensor&, const RKCoefficients<Dimension>&);

  void (*mComputeCorrections)(const ConnectivityMap<Dimension>&,
                              const TableKernel<Dimension>&,
                              const FieldList<Dimension, Scalar>&,
                              const FieldList<Dimension, Vector>&,
                              const FieldList<Dimension, SymTensor>&,
                              const bool,
                              FieldList<Dimension, RKCoefficients<Dimension>>&,
                              FieldList<Dimension, RKCoefficients<Dimension>>&);
  void (*mComputeNormal)(const ConnectivityMap<Dimension>&,
                         const TableKernel<Dimension>&,
                         const FieldList<Dimension, Scalar>&,
                         const FieldList<Dimension, Vector>&,
                         const FieldList<Dimension, SymTensor>&,
                         const FieldList<Dimension, RKCoefficients<Dimension>>&,
                         FieldList<Dimension, Scalar>&,
                         FieldList<Dimension, Vector>&);
  void (*mGetTransformationMatrix)(const Tensor&,
                                   const bool,
                                   TransformationMatrix&);
  void (*mApplyTransformation)(const TransformationMatrix& T,
                               RKCoefficients<Dimension>& corrections);
};

}

#include "RK/ReproducingKernelMethodsInline.hh"

#endif

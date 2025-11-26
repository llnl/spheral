namespace Spheral {

//------------------------------------------------------------------------------
// The reflection operator appropriate between the given master and slave 
// positions.
//------------------------------------------------------------------------------
inline
Dim<3>::Tensor
SphericalBoundary::
reflectOperator(const Dim<3>::Vector& r0,
                const Dim<3>::Vector& r1) const {
  const Vector normal = (r0 - r1).unitVector();
  const Tensor result = Tensor::one() - 2.0*normal.selfdyad();
  ENSURE(fuzzyEqual(std::abs(result.Determinant()), 1.0, 1.0e-10));
  return result;
}

}

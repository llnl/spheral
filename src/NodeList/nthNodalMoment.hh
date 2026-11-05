//---------------------------------Spheral++----------------------------------//
// nthNodalMoment
// nodalMoments
//
// Compute the nth moment of the local nodal distribution in \eta space:
//
//    \sum_j  (\eta_i)^n W_ij
//    -----------------------
//         \sum_j W_ij
//
// Nodal moments is specialized to simultaneously compute the non-normalized 
// zeroth and normalized first moments.
//
// Created by JMO, Mon May  9 16:20:18 PDT 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral_nthNodalMoment__
#define __Spheral_nthNodalMoment__

namespace Spheral {

// Forward declarations.
template<typename Dimension, typename Value> class FieldList;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;

// A trait class to figure out the result type.
template<typename Dimension, unsigned n> struct MomentTraits;
template<typename Dimension> struct MomentTraits<Dimension, 0U> { typedef typename Dimension::Scalar Moment; };
template<typename Dimension> struct MomentTraits<Dimension, 1U> { typedef typename Dimension::Vector Moment; };
template<typename Dimension> struct MomentTraits<Dimension, 2U> { typedef typename Dimension::SymTensor Moment; };

// Compute a particular moment.
template<typename Dimension, unsigned moment>
FieldList<Dimension, typename MomentTraits<Dimension, moment>::Moment>
nthNodalMoment(const DataBase<Dimension>& dataBase,
               const TableKernel<Dimension>& W,
               const bool renormalize);

// Compute the non-normalized zeroth and normalized first moment.
template<typename Dimension>
void
zerothAndFirstNodalMoments(const DataBase<Dimension>& dataBase,
                           const TableKernel<Dimension>& W,
                           const bool useKernelAsGradient,
                           FieldList<Dimension, typename Dimension::Scalar>& zerothMoment,
                           FieldList<Dimension, typename Dimension::Vector>& firstMoment);

}

#endif

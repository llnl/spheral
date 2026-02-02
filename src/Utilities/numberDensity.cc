//---------------------------------Spheral++----------------------------------//
// numberDensity
//
// Compute the number density using the ASPH sum of the nodes.
//
// Created by JMO, Tue Feb  9 13:51:33 PST 2010
//----------------------------------------------------------------------------//
#include "numberDensity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Kernel/TableKernel.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

//------------------------------------------------------------------------------
// Return the number density for all nodes in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Scalar>
numberDensity(const DataBase<Dimension>& dataBase,
              const TableKernel<Dimension>& W) {

  typedef typename Dimension::Scalar Scalar;

  // Get the state.
  const auto  position = dataBase.globalPosition();
  const auto  H = dataBase.globalHfield();
  const auto& cm = dataBase.connectivityMap();
  const auto  numNodeLists = dataBase.numNodeLists();

  // Prepare the result.
  FieldList<Dimension, Scalar> result = dataBase.newGlobalFieldList(0.0, "number density");

  // Some useful variables.
  const auto W0 = W.kernelValue(0.0, 1.0);

  // The set of interacting node pairs.
  const auto& pairs = cm.nodePairList();
  const auto  npairs = pairs.size();

  // First the self contribution.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = result[nodeListi]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      result(nodeListi, i) += Hdeti*W0;
    }
  }

  // Now the pair contributions.
#pragma omp parallel
  {
    int i, j, nodeListi, nodeListj;
    auto result_thread = result.threadCopy();

#pragma omp for
    for (auto k = 0u; k < npairs; ++k) {
      i = pairs[k].i_node;
      j = pairs[k].j_node;
      nodeListi = pairs[k].i_list;
      nodeListj = pairs[k].j_list;

      const auto& ri = position(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();

      const auto& rj = position(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();

      const auto rij = ri - rj;
      const auto etai = (Hi*rij).magnitude();
      const auto etaj = (Hj*rij).magnitude();
      const auto Wi = W.kernelValue(etai, Hdeti);
      const auto Wj = W.kernelValue(etaj, Hdetj);

      result(nodeListi, i) += Wi;
      result(nodeListj, j) += Wj;
    }

#pragma omp critical
    {
      result_thread.threadReduce();
    }
  }

  // That's it.
  return result;
}

}


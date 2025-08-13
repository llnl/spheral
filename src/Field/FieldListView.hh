#ifndef __Spheral_FieldListView_hh__
#define __Spheral_FieldListView_hh__

#include "chai/ManagedArray.hpp"
#include "config.hh"

namespace Spheral {

template <typename Dimension, typename DataType>
class FieldListView {

  using ElementType = FieldView<Dimension, DataType>;
  using ContainerType = chai::ManagedArray<ElementType>;
  ContainerType mData;

public:

  FieldListView(const ContainerType& d) : mData(d) {}

  SPHERAL_HOST_DEVICE DataType& operator()(const size_t fieldIndex,
                                           const size_t nodeIndx) {
    return mData[fieldIndex][nodeIndx];
  }

  SPHERAL_HOST_DEVICE DataType& operator()(const size_t fieldIndex,
                                           const size_t nodeIndx) const {
    return mData[fieldIndex][nodeIndx];
  }

  SPHERAL_HOST_DEVICE ElementType operator[](const size_t index) {
    return mData[index];
  }

  SPHERAL_HOST_DEVICE ElementType operator[](const size_t index) const {
    return mData[index];
  }

  void move(chai::ExecutionSpace space, bool recursive = true) {
    mData.move(space);
    if (recursive) {
      for (auto d: mData) {
        d.move(space);
      }
    }
  }

  SPHERAL_HOST_DEVICE size_t size() const { return mData.size(); }

  SPHERAL_HOST_DEVICE ElementType* data() const { return mData.data(); }

  SPHERAL_HOST
  ElementType* data(chai::ExecutionSpace space, bool do_move = true) const {
    return mData.data(space, do_move);
  }

  SPHERAL_HOST
  void touch(chai::ExecutionSpace space, bool recursive = true) {
    mData.registerTouch(space);
    if (recursive) {
      for (auto d : mData) {
        d.touch(space);
      }
    }
  }
};
} // namespace Spheral
#endif // __Spheral_FieldListView_hh__

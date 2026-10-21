//---------------------------------Spheral++----------------------------------//
// DEMNodeList -- 
// 3/10/2021
//----------------------------------------------------------------------------//
#ifndef __Spheral_DEMNodeList__
#define __Spheral_DEMNodeList__

#include "NodeList.hh"
#include "Geometry/Dimension.hh"

#include <float.h>
#include <string>

namespace Spheral {

template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;

class FileIO;

template<typename Dimension>
class DEMNodeList: public NodeList<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  // Constructors.
  DEMNodeList(std::string name,
              const size_t numInternal,
              const size_t numGhost,
              const Scalar hmin,
              const Scalar hmax,
              const Scalar hminratio,
              const Scalar nPerh,
              const Scalar neighborSearchBuffer,
              const size_t maxNumNeighbors);

  // Destructor
  virtual ~DEMNodeList() = default;

  Field<Dimension, Scalar>& particleRadius();
  const Field<Dimension, Scalar>& particleRadius() const;
  void particleRadius(const Field<Dimension, Scalar>& particleRadius);

  Field<Dimension, int>& compositeParticleIndex();
  const Field<Dimension, int>& compositeParticleIndex() const;
  void compositeParticleIndex(const Field<Dimension, int>& compositeParticleIndex);

  Field<Dimension, size_t>& uniqueIndex();
  const Field<Dimension, size_t>& uniqueIndex() const;
  void uniqueIndex(const Field<Dimension, size_t>& uniqueIndex);

  Scalar neighborSearchBuffer() const;
  void   neighborSearchBuffer(Scalar x);

  void setHfieldFromParticleRadius(const size_t startUniqueIndex);

  //****************************************************************************
  // Methods required for restarting.
  // Dump and restore the NodeList state.
  virtual std::string label() const { return "DEMNodeList"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

  // No default constructor or copying.
  DEMNodeList() = delete;
  DEMNodeList(const DEMNodeList& nodes) = delete;
  DEMNodeList& operator=(const DEMNodeList& rhs) = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  Scalar mNeighborSearchBuffer;
  Field<Dimension, Scalar> mParticleRadius;
  Field<Dimension, int> mCompositeParticleIndex;
  Field<Dimension, size_t> mUniqueIndex;
};

}

#include "DEMNodeListInline.hh"

#endif


//------------------------------------------------------------------------------
// GeometryRegistrar
//
// A singleton class for maintaining global information about the type of
// coordinate Spheral is running in.
//
// Created by JMO, Thu Mar 11 16:41:33 PST 2021
//------------------------------------------------------------------------------

#include "Geometry/GeometryRegistrar.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Instance
//------------------------------------------------------------------------------
GeometryRegistrar&
GeometryRegistrar::instance() {
  static GeometryRegistrar theInstance;
  return theInstance;
}

//------------------------------------------------------------------------------
// Constructor (private)
//------------------------------------------------------------------------------
GeometryRegistrar::GeometryRegistrar(){
}

//------------------------------------------------------------------------------
// Destructor (private)
//------------------------------------------------------------------------------
GeometryRegistrar::~GeometryRegistrar() {
}

//------------------------------------------------------------------------------
// Initializations
//------------------------------------------------------------------------------
CoordinateType GeometryRegistrar::mCoords = CoordinateType::Cartesian;

}

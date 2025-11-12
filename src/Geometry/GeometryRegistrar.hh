//------------------------------------------------------------------------------
// GeometryRegistrar
//
// A singleton class for maintaining global information about the type of
// coordinate Spheral is running in.
//
// Created by JMO, Thu Mar 11 16:41:33 PST 2021
//------------------------------------------------------------------------------
#ifndef __Spheral_GeometryRegistrar__
#define __Spheral_GeometryRegistrar__

namespace Spheral {

// Enum to capture cartesian vs. curvilinear coordinates
enum class CoordinateType {
  Cartesian = 0,
  Spherical = 1,
  RZ = 2,
};

class GeometryRegistrar {
public:
  //--------------------------- Public Interface ---------------------------//
  // Get the instance.
  static GeometryRegistrar& instance();

  // The attribute we hang onto defining the local coordinate system.
  static CoordinateType coords()             { return mCoords; }
  static void coords(const CoordinateType x) { mCoords = x; }

private:
  //--------------------------- Private Interface --------------------------//
  static CoordinateType mCoords;

  // No public constructors, destructor, or assignment.
  GeometryRegistrar();
  GeometryRegistrar(const GeometryRegistrar&);
  GeometryRegistrar& operator=(const GeometryRegistrar&);
  ~GeometryRegistrar();
};

}

#endif

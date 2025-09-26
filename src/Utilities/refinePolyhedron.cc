//---------------------------------Spheral++----------------------------------//
// refinePolyhedron
//
// Refine/smooth a polyhedron.
//
// We use the Pixar OpenSubdiv package for our work here.  This implementation
// is based on an example I found in their source:  tutorials/far/tutorial_0/far_tutorial_0.cpp
//
// Created by JMO, Tue Mar 18 15:23:04 PDT 2014
//   revised, JMO, Thu Jan 26 20:56:45 PST 2017
//----------------------------------------------------------------------------//
#include "refinePolyhedron.hh"
#include "Geometry/Dimension.hh"
#include "config.hh"

// We use Pixar's opensubdiv package to do the refinement.
#ifdef SPHERAL_ENABLE_OPENSUBDIV
#include "opensubdiv/far/topologyDescriptor.h"
#include "opensubdiv/far/primvarRefiner.h"
#endif

#include <vector>
#include <algorithm>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

namespace {   // anonymous 

//------------------------------------------------------------------------------
// Vertex container implementation.
// Directly cribbed from the Opensubdiv example.
struct Vertex {

    // Minimal required interface ----------------------
    Vertex() { }

    Vertex(Vertex const & src) {
        _position[0] = src._position[0];
        _position[1] = src._position[1];
        _position[2] = src._position[2];
    }

    // This method added by JMO
    Vertex(Dim<3>::Vector const & src) {
        _position[0] = src[0];
        _position[1] = src[1];
        _position[2] = src[2];
    }

    void Clear( void * =0 ) {
        _position[0]=_position[1]=_position[2]=0.0f;
    }

    void AddWithWeight(Vertex const & src, float weight) {
        _position[0]+=weight*src._position[0];
        _position[1]+=weight*src._position[1];
        _position[2]+=weight*src._position[2];
    }

    // Public interface ------------------------------------
    void SetPosition(float x, float y, float z) {
        _position[0]=x;
        _position[1]=y;
        _position[2]=z;
    }

    const float * GetPosition() const {
        return _position;
    }

private:
    float _position[3];
};

//------------------------------------------------------------------------------

}             // anonymous

GeomPolyhedron refinePolyhedron(const GeomPolyhedron& poly0,
                                const unsigned numLevels) {
  CONTRACT_VAR(poly0);
  CONTRACT_VAR(numLevels);

#ifdef SPHERAL_ENABLE_OPENSUBDIV

  using namespace OpenSubdiv;
  typedef Dim<3>::Vector Vector;
  typedef GeomPolyhedron::Facet Facet;
  typedef Far::TopologyDescriptor Descriptor;

  // Start by reading out the Spheral polyhedral data to an Opensubdiv topology
  // to shove into opensubdiv.
  const vector<Vector>& verts0 = poly0.vertices();
  const vector<vector<unsigned> >& facetVerts0 = poly0.facetVertices();
  const unsigned numVertices0 = verts0.size();
  const unsigned numFaces0 = facetVerts0.size();
  // float g_verts[numVertices0][3];
  vector<int> g_vertsperface(numFaces0);
  unsigned vertsPerFaceSum = 0;
  {
    // for (unsigned i = 0; i != numVertices0; ++i) {
    //   g_verts[i][0] = verts0[i][0];
    //   g_verts[i][1] = verts0[i][1];
    //   g_verts[i][2] = verts0[i][2];
    // }
    for (unsigned i = 0; i != numFaces0; ++i) {
      g_vertsperface[i] = facetVerts0[i].size();
      vertsPerFaceSum += facetVerts0[i].size();
    }
  }
  vector<int> g_vertIndices(vertsPerFaceSum);
  {
    unsigned j = 0;
    for (const auto& inds: facetVerts0) {
      for (const auto i: inds) {
        CHECK(j < vertsPerFaceSum);
        g_vertIndices[j] = i;
        ++j;
      }
    }
    CHECK(j == vertsPerFaceSum);
  }

  vector<float> positions0(3*numVertices0);
  vector<int> faceVertices0, numVertsPerFace0(numFaces0);
  {
    const vector<Vector>& verts0 = poly0.vertices();
    const vector<Facet>& facets0 = poly0.facets();
    for (unsigned i = 0; i != numVertices0; ++i) {
      positions0[3*i  ] = verts0[i].x();
      positions0[3*i+1] = verts0[i].y();
      positions0[3*i+2] = verts0[i].z();
    }
    for (unsigned i = 0; i != numFaces0; ++i) {
      const vector<unsigned>& ipoints = facets0[i].ipoints();
      numVertsPerFace0[i] = ipoints.size();
      for (unsigned j = 0; j != ipoints.size(); ++j) {
        faceVertices0.push_back(ipoints[j]);
      }
    }
  }
    
  // Initialize the OSD topology.
  Sdc::SchemeType type = OpenSubdiv::Sdc::SCHEME_CATMARK;
  Sdc::Options options;
  options.SetVtxBoundaryInterpolation(Sdc::Options::VTX_BOUNDARY_EDGE_ONLY);

  Descriptor desc;
  desc.numVertices  = numVertices0;
  desc.numFaces     = numFaces0;
  desc.numVertsPerFace = &g_vertsperface[0];
  desc.vertIndicesPerFace  = &g_vertIndices[0];

  // Instantiate a FarTopologyRefiner from the descriptor
  Far::TopologyRefiner * refiner = Far::TopologyRefinerFactory<Descriptor>::Create(desc,
                                                                                   Far::TopologyRefinerFactory<Descriptor>::Options(type, options));

  // Uniformly refine the topology up to numLevels.
  refiner->RefineUniform(Far::TopologyRefiner::UniformOptions(numLevels));

  // Allocate a buffer for vertex primvar data. The buffer length is set to
  // be the sum of all children vertices up to the highest level of refinement.
  std::vector<Vertex> vbuffer(refiner->GetNumVerticesTotal());
  Vertex * verts = &vbuffer[0];

  // Initialize coarse mesh positions
  for (unsigned i = 0; i != numVertices0; ++i) {
    verts[i].SetPosition(verts0[i][0], verts0[i][1], verts0[i][2]);
  }

  // Interpolate vertex primvar data
  Far::PrimvarRefiner primvarRefiner(*refiner);

  Vertex * src = verts;
  for (unsigned level = 1; level <= numLevels; ++level) {
    Vertex * dst = src + refiner->GetLevel(level-1).GetNumVertices();
    primvarRefiner.Interpolate(level, src, dst);
    src = dst;
  }

  // Construct the final refined polyhedron.
  Far::TopologyLevel const & refLastLevel = refiner->GetLevel(numLevels);
  const unsigned numVertices1 = refLastLevel.GetNumVertices();
  const unsigned numFacets1 = refLastLevel.GetNumFaces();
  vector<Vector> vertices1;
  vertices1.reserve(numVertices1);
  const unsigned firstOfLastVerts = refiner->GetNumVerticesTotal() - numVertices1;
  for (unsigned i = 0; i != numVertices1; ++i) {
    float const * pos = verts[firstOfLastVerts + i].GetPosition();
    vertices1.push_back(Vector(pos[0],
                               pos[1],
                               pos[2]));
  }
  CHECK(vertices1.size() == numVertices1);
  vector<vector<unsigned> > facets1(numFacets1);
  for (unsigned i = 0; i != numFacets1; ++i) {
    Far::ConstIndexArray fverts = refLastLevel.GetFaceVertices(i);
    const unsigned nvf = fverts.size();
    for (unsigned j = 0; j != nvf; ++j) {
      facets1[i].push_back(fverts[j]);
    }
  }
  GeomPolyhedron poly1(vertices1, facets1);

  // That's it.
  return poly1;
#else
  VERIFY2(false, "ERROR: attempt to call refinePolyhedron, but OpenSubdiv has not been compiled into Spheral.");
  return GeomPolyhedron();
#endif
}

}

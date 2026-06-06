//---------------------------------Spheral++----------------------------------//
// PolyClipperUtilities
//
// Methods for talking with PolyClipper:
//       https://github.com/LLNL/PolyClipper
//
// Created by JMO, Thu Feb  4 16:24:44 PST 2021
//----------------------------------------------------------------------------//

#include "Geometry/PolyClipperUtilities.hh"
#include "Utilities/DBC.hh"
#include "Utilities/Timer.hh"


namespace Spheral {

//------------------------------------------------------------------------------
// Spheral::GeomPolygon -> PolyClipper::Polygon.
//------------------------------------------------------------------------------
void convertToPolyClipper(PolyClipperPolygon& polygon,
                          const Dim<2>::FacetedVolume& Spheral_polygon) {
  ADV_TIME_BEGIN("PC2d_convertto");

  // Construct the vertices without connectivity first.
  const auto& coords = Spheral_polygon.vertices();
  const auto  nverts0 = coords.size();
  polygon.resize(nverts0);
  for (auto i = 0u; i < nverts0; ++i) polygon[i].position = coords[i];

  // Build the connectivity.
  const auto& facets = Spheral_polygon.facets();
  int v1, v2;
  for (const auto& facet: facets) {
    v1 = facet.ipoint1();
    v2 = facet.ipoint2();
    CHECK(v1 < (int)nverts0);
    CHECK(v2 < (int)nverts0);
    polygon[v1].neighbors.second = v2;
    polygon[v2].neighbors.first  = v1;
  }

  CHECK(polygon.size() == nverts0);
  ADV_TIME_END("PC2d_convertto");
}

//------------------------------------------------------------------------------
// PolyClipper::Polygon -> Spheral::GeomPolygon.
// The return value is the set of plane IDs responsible for each vertex.
//------------------------------------------------------------------------------
vector<set<int>> convertFromPolyClipper(Dim<2>::FacetedVolume& Spheral_polygon,
                                        const PolyClipperPolygon& polygon) {
  ADV_TIME_BEGIN("PC2d_convertfrom");

  // Useful types.
  using FacetedVolume = Spheral::Dim<2>::FacetedVolume;
  using Vector = Spheral::Dim<2>::Vector;

  vector<set<int>> vertexPlanes;

  if (polygon.empty()) {

    Spheral_polygon = FacetedVolume();

  } else {

    // Numbers of vertices.
    const auto nverts = polygon.size();
    const auto nactive = count_if(polygon.begin(), polygon.end(),
                                  [](const PolyClipperVertex2d& x) { return x.comp >= 0; });
    set<int> usedVertices;

    // Go until we hit all the active vertices.
    vector<Vector> coords(nactive);
    vector<vector<unsigned>> facets(nactive, vector<unsigned>(2));
    vertexPlanes.resize(nactive);
    auto k = 0, loopStart = 0;
    while ((int)usedVertices.size() < nactive) {

      // Look for the first active unused vertex.
      unsigned vstart = 0;
      while (vstart < nverts and
             (polygon[vstart].comp < 0 or usedVertices.find(vstart) != usedVertices.end())) vstart++;
      CHECK(vstart < nverts);
      auto vnext = vstart;

      // Read out this loop.
      auto force = true;
      while (force or vnext != vstart) {
        CHECK2(k < nactive, polygon2string(polygon));
        coords[k] = polygon[vnext].position;
        facets[k][0] = k;
        facets[k][1] = k + 1;
        vertexPlanes[k].insert(polygon[vnext].clips.begin(), polygon[vnext].clips.end());
        ++k;
        force = false;
        usedVertices.insert(vnext);
        vnext = polygon[vnext].neighbors.second;
      }
      facets[k-1][1] = loopStart;
      loopStart = k;
    }
    CHECK(k == nactive);

    Spheral_polygon = FacetedVolume(coords, facets);

  }

  // Return the set of planes responsible for each vertex.
  ENSURE(vertexPlanes.size() == Spheral_polygon.vertices().size());
  ADV_TIME_END("PC2d_convertfrom");
  return vertexPlanes;
}

//------------------------------------------------------------------------------
// Spheral::GeomPolyhedron -> PolyClipper::Polyhedron.
//------------------------------------------------------------------------------
void convertToPolyClipper(PolyClipperPolyhedron& polyhedron,
                          const Dim<3>::FacetedVolume& Spheral_polyhedron) {
  ADV_TIME_BEGIN("PC3d_convertto");

  const auto& vertPositions = Spheral_polyhedron.vertices();
  const auto& facets = Spheral_polyhedron.facets();
  const auto  nverts = vertPositions.size();

  // Build the PolyClipper Vertex3d's, but without connectivity yet.
  polyhedron.resize(nverts);
  for (auto k = 0u; k < nverts; ++k) {
    polyhedron[k] = PolyClipperVertex3d(vertPositions[k], 1);
  }

  // Note all the edges associated with each vertex.
  vector<vector<pair<int, int>>> vertexPairs(nverts);
  int v0, vprev, vnext;
  for (const auto& facet: facets) {
    const auto& ipoints = facet.ipoints();
    const int   n = ipoints.size();
    for (int k = 0; k < n; ++k) {
      v0 = ipoints[k];
      vprev = ipoints[(k - 1 + n) % n];
      vnext = ipoints[(k + 1 + n) % n];
      vertexPairs[v0].push_back(make_pair(vnext, vprev));   // Gotta reverse!
    }
  }

  // Now we can build vertex->vertex connectivity.
  for (auto k = 0u; k < nverts; ++k) {

    // Sort the edges associated with this vertex.
    const auto n = vertexPairs[k].size();
    CHECK(n >= 3);
    for (auto i = 0u; i < n - 1; ++i) {
      auto j = i + 1;
      while (j < n and vertexPairs[k][i].second != vertexPairs[k][j].first) ++j;
      CHECK(j < n);
      swap(vertexPairs[k][i + 1], vertexPairs[k][j]);
    }

    // Now that they're in order, create the neighbors.
    for (auto i = 0u; i < n; ++i) {
      polyhedron[k].neighbors.push_back(vertexPairs[k][i].first);
    }
    CHECK(polyhedron[k].neighbors.size() == vertexPairs[k].size());
  }

  CHECK(polyhedron.size() == nverts);
  ADV_TIME_END("PC3d_convertto");
}

//------------------------------------------------------------------------------
// PolyClipper::Polyhedron -> Spheral::GeomPolyhedron.
//------------------------------------------------------------------------------
vector<set<int>> convertFromPolyClipper(Dim<3>::FacetedVolume& Spheral_polyhedron,
                                        const PolyClipperPolyhedron& polyhedron) {
  ADV_TIME_BEGIN("PC3d_convertfrom");

  // Useful types.
  using FacetedVolume = Spheral::Dim<3>::FacetedVolume;
  using Vector = Spheral::Dim<3>::Vector;

  vector<set<int>> vertexPlanes;

  if (polyhedron.empty()) {

    Spheral_polyhedron = FacetedVolume();

  } else {

    // extractFaces actually does most of the work.
    const auto faces = extractFaces(polyhedron);

    // Number and extract the active vertices.
    vector<Vector> coords;
    {
      int i = 0;
      for (auto itr = polyhedron.begin(); itr != polyhedron.end(); ++itr) {
        if (itr->comp >= 0) {
          coords.push_back(itr->position);
          vertexPlanes.push_back(itr->clips);
          itr->ID = i++;
        }
      }
    }
    CHECK((int)coords.size() == count_if(polyhedron.begin(), polyhedron.end(),
                                         [](const PolyClipperVertex3d& x) { return x.comp >= 0; }));

    // Extract the faces as integer vertex index loops.
    vector<vector<unsigned>> facets(faces.size());
    for (auto k = 0u; k < faces.size(); ++k) {
      facets[k].resize(faces[k].size());
      transform(faces[k].begin(), faces[k].end(), facets[k].begin(),
                [&](const int k) { return polyhedron[k].ID; });
    }

    // Now we can build the Spheral::Polyhedron.
    Spheral_polyhedron = FacetedVolume(coords, facets);

  }

  // Return the set of planes responsible for each vertex.
  ENSURE(vertexPlanes.size() == Spheral_polyhedron.vertices().size());
  ADV_TIME_END("PC3d_convertfrom");

  return vertexPlanes;
}

}

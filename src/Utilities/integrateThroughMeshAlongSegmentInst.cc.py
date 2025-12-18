text = """
//------------------------------------------------------------------------------
// Explicit instantiations.
//------------------------------------------------------------------------------
#include "Utilities/integrateThroughMeshAlongSegment.cc"

//------------------------------------------------------------------------------
// Find the points of intersection with the mesh planes for the given segment.
//------------------------------------------------------------------------------

namespace Spheral {

    template Dim< %(ndim)s >::Scalar integrateThroughMeshAlongSegment<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>(const vector<vector<Dim< %(ndim)s >::Scalar> >& values, const Dim< %(ndim)s >::Vector& xmin, const Dim< %(ndim)s >::Vector& xmax, const vector<unsigned>& ncells, const Dim< %(ndim)s >::Vector& s0, const Dim< %(ndim)s >::Vector& s1);

}
"""

specializations = {
    1:
"""
    vector<Dim<1>::Vector>
    findIntersections(const Dim<1>::Vector& xmin,
            const Dim<1>::Vector& xmax,
            const vector<unsigned>& ncells,
            const Dim<1>::Vector& s0,
            const Dim<1>::Vector& s1) {
      REQUIRE(xmin.x() < xmax.x());
      REQUIRE(ncells.size() == 1);
      REQUIRE(ncells[0] > 0);

      // Find the min and max bounding mesh planes indicies.
      typedef Dim<1>::Vector Vector;
      const Vector smin = elementWiseMin(s0, s1);
      const Vector smax = elementWiseMax(s0, s1);
      const unsigned ixmin = rightCellBoundary(smin(0), xmin(0), xmax(0), ncells[0]);
      const unsigned ixmax =  leftCellBoundary(smax(0), xmin(0), xmax(0), ncells[0]);
      CHECK(ixmin <= ncells[0]);
      CHECK(ixmax <= ncells[0]);

      // The intersections are just the intermediate planes.
      vector<Vector> result;
      const double xstep = (xmax.x() - xmin.x())/ncells[0];
      for (unsigned iplane = ixmin; iplane < ixmax; ++iplane) {
        result.push_back(xmin + Dim<1>::Vector(iplane*xstep));
      }

      // Post-conditions.
      BEGIN_CONTRACT_SCOPE
      {
        if (result.size() > 0) {
          ENSURE(result.front().x() >= smin.x());
          ENSURE(result.back().x() <= smax.x());
        }
      }
      END_CONTRACT_SCOPE

      return result;
    }
""",
    2:
"""
    vector<Dim<2>::Vector>
    findIntersections(const Dim<2>::Vector& xmin,
                      const Dim<2>::Vector& xmax,
                      const vector<unsigned>& ncells,
                      const Dim<2>::Vector& s0,
                      const Dim<2>::Vector& s1) {
      REQUIRE(xmin.x() < xmax.x() and xmin.y() < xmax.y());
      REQUIRE(ncells.size() == 2);
      REQUIRE(ncells[0] > 0 and ncells[1] > 0);

      // Find the min and max bounding mesh planes indicies.
      typedef Dim<2>::Vector Vector;
      const Vector smin = elementWiseMin(s0, s1);
      const Vector smax = elementWiseMax(s0, s1);
      vector<unsigned> ixmin, ixmax;
      for (size_t idim = 0; idim != 2; ++idim) {
        ixmin.push_back(rightCellBoundary(smin(idim), xmin(idim), xmax(idim), ncells[idim]));
        ixmax.push_back( leftCellBoundary(smax(idim), xmin(idim), xmax(idim), ncells[idim]));
        CHECK(ixmin.back() <= ncells[idim]);
        CHECK(ixmax.back() <= ncells[idim]);
      }
      CHECK(ixmin.size() == 2);
      CHECK(ixmax.size() == 2);

      // The intersections are just the intermediate planes.
      vector<Vector> result;
      const double xstep = (xmax.x() - xmin.x())/ncells[0];
      const double ystep = (xmax.y() - xmin.y())/ncells[1];
      for (unsigned iplane = ixmin[0]; iplane < ixmax[0]; ++iplane) {
        const double xseg = xmin.x() + iplane*xstep;
        const Vector meshSeg0(xseg, -1.0e10);
        const Vector meshSeg1(xseg,  1.0e10);
        Vector intersect1, intersect2;
        const char test = segmentSegmentIntersection(s0, s1, meshSeg0, meshSeg1, intersect1, intersect2);
        CONTRACT_VAR(test);
        CHECK(test != '0');
        result.push_back(intersect1);
      }
      for (unsigned iplane = ixmin[1]; iplane < ixmax[1]; ++iplane) {
        const double yseg = xmin.y() + iplane*ystep;
        const Vector meshSeg0(-1.0e10, yseg);
        const Vector meshSeg1( 1.0e10, yseg);
        Vector intersect1, intersect2;
        const char test = segmentSegmentIntersection(s0, s1, meshSeg0, meshSeg1, intersect1, intersect2);
        CONTRACT_VAR(test);
        CHECK(test != '0');
        result.push_back(intersect1);
      }

      // Post-conditions.
      BEGIN_CONTRACT_SCOPE
      {
        const Vector shat = (s1 - s0).unitVector();
        CONTRACT_VAR(shat);
        const double segLen = (s1 - s0).magnitude();
        for (vector<Vector>::const_iterator itr = result.begin();
             itr != result.end();
             ++itr) {
          CONTRACT_VAR(segLen);
          ENSURE((*itr - s0).dot(shat) >= 0.0);
          ENSURE(fuzzyEqual(abs((*itr - s0).dot(shat)), (*itr - s0).magnitude(), 1.0e-10));
          ENSURE((*itr - s0).dot(shat)*safeInv(segLen) <= (1.0 + 1.0e-8));
        }
      }
      END_CONTRACT_SCOPE

      return result;
    }
""",
    3:
"""
    vector<Dim<3>::Vector>
    findIntersections(const Dim<3>::Vector& xmin,
                      const Dim<3>::Vector& xmax,
                      const vector<unsigned>& ncells,
                      const Dim<3>::Vector& s0,
                      const Dim<3>::Vector& s1) {
      REQUIRE(xmin.x() < xmax.x() and xmin.y() < xmax.y() and xmin.z() < xmax.z());
      REQUIRE(ncells.size() == 3);
      REQUIRE(ncells[0] > 0 and ncells[1] > 0 and ncells[2] > 0);

      // Find the min and max bounding mesh planes indicies.
      typedef Dim<3>::Vector Vector;
      const Vector smin = elementWiseMin(s0, s1);
      const Vector smax = elementWiseMax(s0, s1);
      vector<unsigned> ixmin, ixmax;
      for (size_t idim = 0; idim != 3; ++idim) {
        ixmin.push_back(rightCellBoundary(smin(idim), xmin(idim), xmax(idim), ncells[idim]));
        ixmax.push_back( leftCellBoundary(smax(idim), xmin(idim), xmax(idim), ncells[idim]));
        CHECK(ixmin.back() <= ncells[idim]);
        CHECK(ixmax.back() <= ncells[idim]);
      }
      CHECK(ixmin.size() == 3);
      CHECK(ixmax.size() == 3);

      // The intersections are just the intermediate planes.
      vector<Vector> result;
      const double xstep = (xmax.x() - xmin.x())/ncells[0];
      const double ystep = (xmax.y() - xmin.y())/ncells[1];
      const double zstep = (xmax.z() - xmin.z())/ncells[2];
      for (unsigned iplane = ixmin[0]; iplane < ixmax[0]; ++iplane) {
        const double xplane = xmin.x() + iplane*xstep;
        const Vector meshPlane0(xplane, 0.0, 0.0);
        const Vector meshPlane1(xplane, 0.0, 1.0);
        const Vector meshPlane2(xplane, 1.0, 1.0);
        Vector intersect;
        const char test = segmentPlaneIntersection(s0, s1, meshPlane0, meshPlane1, meshPlane2, intersect);
        CONTRACT_VAR(test);
        CHECK(test != '0');
        result.push_back(intersect);
      }
      for (unsigned iplane = ixmin[1]; iplane < ixmax[1]; ++iplane) {
        const double yplane = xmin.y() + iplane*ystep;
        const Vector meshPlane0(0.0, yplane, 0.0);
        const Vector meshPlane1(0.0, yplane, 1.0);
        const Vector meshPlane2(1.0, yplane, 1.0);
        Vector intersect;
        const char test = segmentPlaneIntersection(s0, s1, meshPlane0, meshPlane1, meshPlane2, intersect);
        CONTRACT_VAR(test);
        CHECK(test != '0');
        result.push_back(intersect);
      }
      for (unsigned iplane = ixmin[2]; iplane < ixmax[2]; ++iplane) {
        const double zplane = xmin.z() + iplane*zstep;
        const Vector meshPlane0(0.0, 0.0, zplane);
        const Vector meshPlane1(0.0, 1.0, zplane);
        const Vector meshPlane2(1.0, 1.0, zplane);
        Vector intersect;
        const char test = segmentPlaneIntersection(s0, s1, meshPlane0, meshPlane1, meshPlane2, intersect);
        CONTRACT_VAR(test);
        CHECK(test != '0');
        result.push_back(intersect);
      }

      // Post-conditions.
      BEGIN_CONTRACT_SCOPE
      {
        const Vector shat = (s1 - s0).unitVector();
        CONTRACT_VAR(shat);
        const double segLen = (s1 - s0).magnitude();
        for (vector<Vector>::const_iterator itr = result.begin();
             itr != result.end();
             ++itr) {
          CONTRACT_VAR(segLen);
          ENSURE((*itr - s0).dot(shat) >= 0.0);
          ENSURE(fuzzyEqual(abs((*itr - s0).dot(shat)), (*itr - s0).magnitude(), 1.0e-10));
          ENSURE((*itr - s0).dot(shat)*safeInv(segLen) <= 1.0);
        }
      }
      END_CONTRACT_SCOPE

      return result;
    }
"""
}

// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/extrude.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/MultiSolid.h"
#include "SFCGAL/NURBSCurve.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Solid.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"

#include "SFCGAL/Exception.h"

#include "SFCGAL/algorithm/area.h"
#include "SFCGAL/algorithm/force3D.h"
#include "SFCGAL/algorithm/intersection.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/algorithm/normal.h"
#include "SFCGAL/algorithm/plane.h"
#include "SFCGAL/algorithm/translate.h"
#include "SFCGAL/numeric.h"
#include "SFCGAL/triangulate/triangulatePolygon.h"

#include "SFCGAL/detail/tools/Log.h"

#include <CGAL/Line_3.h>
#include <CGAL/Ray_3.h>

#include <utility>

namespace SFCGAL::algorithm {

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

auto
extrude(const Point &g, const Kernel::Vector_3 &v) -> LineString *;
auto
extrude(const LineString &g, const Kernel::Vector_3 &v) -> PolyhedralSurface *;
auto
extrude(const Polygon &g, const Kernel::Vector_3 &v, bool addTop = true)
    -> Solid *;
auto
extrude(const Triangle &g, const Kernel::Vector_3 &v) -> Solid *;

auto
extrude(const MultiPoint &g, const Kernel::Vector_3 &v) -> MultiLineString *;
auto
extrude(const MultiLineString &g, const Kernel::Vector_3 &v)
    -> PolyhedralSurface *;
auto
extrude(const MultiPolygon &g, const Kernel::Vector_3 &v) -> MultiSolid *;

auto
extrude(const TriangulatedSurface &g, const Kernel::Vector_3 &v) -> Solid *;
auto
extrude(const PolyhedralSurface &g, const Kernel::Vector_3 &v) -> Solid *;

auto
extrude(const GeometryCollection &g, const Kernel::Vector_3 &v)
    -> GeometryCollection *;

auto
extrude(const Point &g, const Kernel::Vector_3 &v) -> LineString *
{
  if (g.isEmpty()) {
    return new LineString();
  }

  Kernel::Point_3 const a = g.toPoint_3();
  Kernel::Point_3 const b = a + v;

  return new LineString(Point(a), Point(b));
}

auto
extrude(const LineString &g, const Kernel::Vector_3 &v) -> PolyhedralSurface *
{

  std::unique_ptr<PolyhedralSurface> polyhedralSurface(new PolyhedralSurface());

  if (g.isEmpty()) {
    return polyhedralSurface.release();
  }

  for (size_t i = 0; i < g.numPoints() - 1; i++) {
    std::unique_ptr<LineString> ring(new LineString);

    Kernel::Point_3 const a = g.pointN(i).toPoint_3();
    Kernel::Point_3 const b = g.pointN(i + 1).toPoint_3();
    ring->addPoint(new Point(a));
    ring->addPoint(new Point(b));
    ring->addPoint(new Point(b + v));
    ring->addPoint(new Point(a + v));
    ring->addPoint(new Point(a));

    polyhedralSurface->addPatch(new Polygon(ring.release()));
  }

  return polyhedralSurface.release();
}

auto
extrude(const Polygon &g, const Kernel::Vector_3 &v, bool addTop) -> Solid *
{
  if (g.isEmpty()) {
    return new Solid();
  }

  bool const reverseOrientation = (v * normal3D<Kernel>(g)) > 0;

  // resulting shell
  PolyhedralSurface polyhedralSurface;

  // "bottom"
  Polygon bottom(g);
  force3D(bottom);

  if (reverseOrientation) {
    bottom.reverse();
  }

  polyhedralSurface.addPatch(bottom);

  // "top"
  if (addTop) {
    Polygon top(bottom);
    top.reverse();
    translate(top, v);
    polyhedralSurface.addPatch(top);
  }
  // exterior ring and interior rings extruded
  for (size_t i = 0; i < bottom.numRings(); i++) {
    std::unique_ptr<PolyhedralSurface> boundaryExtruded(
        extrude(bottom.ringN(i), v));

    for (size_t j = 0; j < boundaryExtruded->numPatches(); j++) {
      boundaryExtruded->patchN(j).reverse();
      polyhedralSurface.addPatch(boundaryExtruded->patchN(j));
    }
  }

  return new Solid(polyhedralSurface);
}

auto
extrude(const Triangle &g, const Kernel::Vector_3 &v) -> Solid *
{
  return extrude(g.toPolygon(), v);
}

auto
extrude(const MultiPoint &g, const Kernel::Vector_3 &v) -> MultiLineString *
{
  std::unique_ptr<MultiLineString> result(new MultiLineString());

  if (g.isEmpty()) {
    return result.release();
  }

  for (size_t i = 0; i < g.numGeometries(); i++) {
    result->addGeometry(extrude(g.pointN(i), v));
  }

  return result.release();
}

auto
extrude(const MultiLineString &g, const Kernel::Vector_3 &v)
    -> PolyhedralSurface *
{
  std::unique_ptr<PolyhedralSurface> result(new PolyhedralSurface());

  if (g.isEmpty()) {
    return result.release();
  }

  for (size_t i = 0; i < g.numGeometries(); i++) {
    std::unique_ptr<PolyhedralSurface> extruded(extrude(g.lineStringN(i), v));

    for (size_t j = 0; j < extruded->numPatches(); j++) {
      result->addPatch(extruded->patchN(j));
    }
  }

  return result.release();
}

auto
extrude(const MultiPolygon &g, const Kernel::Vector_3 &v) -> MultiSolid *
{
  std::unique_ptr<MultiSolid> result(new MultiSolid());

  if (g.isEmpty()) {
    return result.release();
  }

  for (size_t i = 0; i < g.numGeometries(); i++) {
    result->addGeometry(extrude(g.polygonN(i), v));
  }

  return result.release();
}

auto
extrude(const TriangulatedSurface &g, const Kernel::Vector_3 &v) -> Solid *
{
  std::unique_ptr<Solid> result(new Solid());

  if (g.isEmpty()) {
    return result.release();
  }

  // bottom and top
  for (size_t i = 0; i < g.numPatches(); i++) {
    Triangle bottomPart(g.patchN(i));
    force3D(bottomPart);
    bottomPart.reverse();
    result->exteriorShell().addPatch(bottomPart);

    Triangle topPart(g.patchN(i));
    force3D(topPart);
    translate(topPart, v);
    result->exteriorShell().addPatch(topPart);
  }

  // boundary
  std::unique_ptr<Geometry> boundary(g.boundary());
  BOOST_ASSERT(boundary.get() != NULL);

  // closed surface extruded
  if (!boundary->isEmpty()) {
    std::unique_ptr<Geometry> extrudedBoundary(extrude(*boundary, v));
    BOOST_ASSERT(extrudedBoundary->is<PolyhedralSurface>());
    result->exteriorShell().addPolygons(
        extrudedBoundary->as<PolyhedralSurface>());
  }

  return result.release();
}

auto
extrude(const PolyhedralSurface &g, const Kernel::Vector_3 &v) -> Solid *
{
  if (g.isEmpty()) {
    return new Solid();
  }

  TriangulatedSurface triangulatedSurface;
  triangulate::triangulatePolygon3D(g, triangulatedSurface);
  return extrude(triangulatedSurface, v);
}

auto
extrude(const GeometryCollection &g, const Kernel::Vector_3 &v)
    -> GeometryCollection *
{
  std::unique_ptr<GeometryCollection> result(new GeometryCollection());

  if (g.isEmpty()) {
    return result.release();
  }

  for (size_t i = 0; i < g.numGeometries(); i++) {
    result->addGeometry(extrude(g.geometryN(i), v).release());
  }

  return result.release();
}

/// @private
auto
extrude(const Geometry &inputGeometry, const Kernel::Vector_3 &vector)
    -> std::unique_ptr<Geometry>
{
  switch (inputGeometry.geometryTypeId()) {
  case TYPE_POINT:
    return std::unique_ptr<Geometry>(
        extrude(inputGeometry.as<Point>(), vector));

  case TYPE_LINESTRING:
    return std::unique_ptr<Geometry>(
        extrude(inputGeometry.as<LineString>(), vector));

  case TYPE_NURBSCURVE: {
    auto lineString =
        inputGeometry.as<NURBSCurve>().toLineString(); // default parameters
    if (!lineString || lineString->isEmpty()) {
      return std::make_unique<PolyhedralSurface>(); // empty result
    }
    return std::unique_ptr<Geometry>(extrude(*lineString, vector));
  }

  case TYPE_POLYGON:
    return std::unique_ptr<Geometry>(
        extrude(inputGeometry.as<Polygon>(), vector));

  case TYPE_TRIANGLE:
    return std::unique_ptr<Geometry>(
        extrude(inputGeometry.as<Triangle>(), vector));

  case TYPE_GEOMETRYCOLLECTION:
    return std::unique_ptr<Geometry>(
        extrude(inputGeometry.as<GeometryCollection>(), vector));

  case TYPE_MULTIPOINT:
    return std::unique_ptr<Geometry>(
        extrude(inputGeometry.as<MultiPoint>(), vector));

  case TYPE_MULTILINESTRING:
    return std::unique_ptr<Geometry>(
        extrude(inputGeometry.as<MultiLineString>(), vector));

  case TYPE_MULTIPOLYGON:
    return std::unique_ptr<Geometry>(
        extrude(inputGeometry.as<MultiPolygon>(), vector));

  case TYPE_TRIANGULATEDSURFACE:
    return std::unique_ptr<Geometry>(
        extrude(inputGeometry.as<TriangulatedSurface>(), vector));

  case TYPE_POLYHEDRALSURFACE:
    return std::unique_ptr<Geometry>(
        extrude(inputGeometry.as<PolyhedralSurface>(), vector));

  case TYPE_SOLID:
  case TYPE_MULTISOLID:
    // extrusion not available
    break;
  }

  BOOST_THROW_EXCEPTION(InappropriateGeometryException(
      (boost::format("Extrusion of %s is not supported") %
       inputGeometry.geometryType())
          .str()));
}

/// @private
auto
extrude(const Geometry &inputGeom, const Kernel::FT &displacementX,
        const Kernel::FT &displacementY, const Kernel::FT &displacementZ,
        NoValidityCheck /*unused*/) -> std::unique_ptr<Geometry>
{
  return extrude(inputGeom,
                 Kernel::Vector_3(displacementX, displacementY, displacementZ));
}

/// @private
auto
extrude(const Geometry &geometry, const Kernel::FT &deltaX,
        const Kernel::FT &deltaY, const Kernel::FT &deltaZ)
    -> std::unique_ptr<Geometry>
{
  SFCGAL_ASSERT_GEOMETRY_VALIDITY(geometry);
  std::unique_ptr<Geometry> result(
      extrude(geometry, deltaX, deltaY, deltaZ, NoValidityCheck()));
  propagateValidityFlag(*result, true);
  return result;
}

/// @private
SFCGAL_API auto
extrude(const Geometry &geom, const double &displacementX,
        const double &displacementY, const double &displacementZ)
    -> std::unique_ptr<Geometry>
{
  if (!std::isfinite(displacementX) || !std::isfinite(displacementY) ||
      !std::isfinite(displacementZ)) {
    BOOST_THROW_EXCEPTION(NonFiniteValueException(
        "trying to extrude with non finite value in direction"));
  }

  return extrude(geom, Kernel::FT(displacementX), Kernel::FT(displacementY),
                 Kernel::FT(displacementZ));
}

/// @private
SFCGAL_API auto
extrude(const Polygon &polygon, const double &height)
    -> std::unique_ptr<Geometry>
{

  if (!std::isfinite(height)) {
    BOOST_THROW_EXCEPTION(NonFiniteValueException(
        "trying to extrude with non finite value in direction"));
  }

  return std::unique_ptr<Geometry>(
      extrude(polygon, Kernel::Vector_3(0.0, 0.0, height), false));
}

/// @private
// NOLINTBEGIN(readability-function-cognitive-complexity)
SFCGAL_API auto
extrudeUntil(const Polygon &footprint, const PolyhedralSurface &roof)
    -> std::unique_ptr<Solid>
{
  if (footprint.isEmpty() || roof.isEmpty()) {
    return std::make_unique<Solid>();
  }

  // Normalize footprint to CCW orientation (positive signed area)
  Polygon normalizedFootprint = footprint;
  if (signedArea(normalizedFootprint.exteriorRing()) < 0) {
    normalizedFootprint.reverse();
  }

  std::unique_ptr<PolyhedralSurface> shell =
      std::make_unique<PolyhedralSurface>();

  // Note: Bottom face will be added later after collecting all boundary
  // vertices

  // Helper to lift a 2D polygon to 3D using a plane
  auto liftPolygon = [](const Polygon               &poly2D,
                        const CGAL::Plane_3<Kernel> &plane) -> Polygon {
    Polygon poly3D = poly2D;
    for (size_t i = 0; i < poly3D.numRings(); ++i) {
      LineString &ring = poly3D.ringN(i);
      for (size_t j = 0; j < ring.numPoints(); ++j) {
        Point     &point = ring.pointN(j);
        Kernel::FT x     = point.x();
        Kernel::FT y     = point.y();
        // z = -(ax + by + d) / c
        if (plane.c() != 0) {
          Kernel::FT z =
              -(plane.a() * x + plane.b() * y + plane.d()) / plane.c();
          point = Point(x, y, z);
        }
      }
    }
    return poly3D;
  };

  std::vector<Polygon> topPolygons;

  // 2. Process roof polygons
  for (size_t i = 0; i < roof.numPolygons(); ++i) {
    const Polygon &roofPoly = roof.polygonN(i);

    // Get plane
    CGAL::Plane_3<Kernel> plane;
    try {
      plane = algorithm::plane3D<Kernel>(roofPoly);
    } catch (...) {
      continue;
    }

    // Ignore vertical faces (normal perpendicular to Z)
    if (plane.c() == 0) {
      continue;
    }

    // Project roof polygon to 2D
    Polygon roofPoly2D = roofPoly;
    for (size_t ringIdx = 0; ringIdx < roofPoly2D.numRings(); ++ringIdx) {
      LineString &ring = roofPoly2D.ringN(ringIdx);
      for (size_t pointIdx = 0; pointIdx < ring.numPoints(); ++pointIdx) {
        ring.pointN(pointIdx) =
            Point(ring.pointN(pointIdx).x(), ring.pointN(pointIdx).y(), 0.0);
      }
    }

    // Intersect with footprint
    std::unique_ptr<Geometry> intersectionResult;
    try {
      intersectionResult = intersection(normalizedFootprint, roofPoly2D);
    } catch (...) {
      continue;
    }

    if (!intersectionResult || intersectionResult->isEmpty()) {
      continue;
    }

    // Collect result polygons
    std::vector<Polygon> currentPolys;
    if (intersectionResult->is<Polygon>()) {
      currentPolys.push_back(intersectionResult->as<Polygon>());
    } else if (intersectionResult->is<Triangle>()) {
      currentPolys.emplace_back(intersectionResult->as<Triangle>());
    } else if (intersectionResult->is<MultiPolygon>()) {
      const MultiPolygon &mp = intersectionResult->as<MultiPolygon>();
      for (size_t k = 0; k < mp.numGeometries(); ++k) {
        currentPolys.push_back(mp.geometryN(k).as<Polygon>());
      }
    } else if (intersectionResult->is<GeometryCollection>()) {
      const GeometryCollection &gc =
          intersectionResult->as<GeometryCollection>();
      for (size_t k = 0; k < gc.numGeometries(); ++k) {
        if (gc.geometryN(k).is<Polygon>()) {
          currentPolys.push_back(gc.geometryN(k).as<Polygon>());
        }
      }
    }

    // Lift and add
    for (const auto &poly : currentPolys) {
      Polygon poly3D = liftPolygon(poly, plane);
      topPolygons.push_back(poly3D);
      shell->addPolygon(poly3D);
    }
  }

  if (topPolygons.empty()) {
    return std::make_unique<Solid>();
  }

  // 3. Create lateral faces from boundary of top surface
  struct Edge {
    Point point1, point2;
    auto
    operator<(const Edge &other) const -> bool
    {
      if (point1.x() != other.point1.x()) {
        return point1.x() < other.point1.x();
      }
      if (point1.y() != other.point1.y()) {
        return point1.y() < other.point1.y();
      }
      if (point1.z() != other.point1.z()) {
        return point1.z() < other.point1.z();
      }
      if (point2.x() != other.point2.x()) {
        return point2.x() < other.point2.x();
      }
      if (point2.y() != other.point2.y()) {
        return point2.y() < other.point2.y();
      }
      return point2.z() < other.point2.z();
    }
  };

  std::map<Edge, int> edgeCount;

  for (const auto &poly : topPolygons) {
    for (size_t ringIdx = 0; ringIdx < poly.numRings(); ++ringIdx) {
      const LineString &ring = poly.ringN(ringIdx);
      for (size_t i = 0; i < ring.numPoints() - 1; ++i) {
        Point edgePoint1 = ring.pointN(i);
        Point edgePoint2 = ring.pointN(i + 1);

        Edge edge;
        if (edgePoint1.x() < edgePoint2.x() ||
            (edgePoint1.x() == edgePoint2.x() &&
             edgePoint1.y() < edgePoint2.y()) ||
            (edgePoint1.x() == edgePoint2.x() &&
             edgePoint1.y() == edgePoint2.y() &&
             edgePoint1.z() < edgePoint2.z())) {
          edge.point1 = edgePoint1;
          edge.point2 = edgePoint2;
        } else {
          edge.point1 = edgePoint2;
          edge.point2 = edgePoint1;
        }
        edgeCount[edge]++;
      }
    }
  }

  // Collect oriented boundary edges
  std::vector<std::pair<Point, Point>> boundaryEdges;

  for (const auto &[edge, count] : edgeCount) {
    if (count == 1) {
      // Boundary edge. Find orientation.
      Point edgeStart    = edge.point1;
      Point edgeEnd      = edge.point2;
      bool  foundForward = false;

      for (const auto &poly : topPolygons) {
        for (size_t ringIdx = 0; ringIdx < poly.numRings(); ++ringIdx) {
          const LineString &ring = poly.ringN(ringIdx);
          for (size_t i = 0; i < ring.numPoints() - 1; ++i) {
            if (ring.pointN(i) == edgeStart && ring.pointN(i + 1) == edgeEnd) {
              foundForward = true;
              goto found;
            }
          }
        }
      }
    found:

      Point start = foundForward ? edgeStart : edgeEnd;
      Point end   = foundForward ? edgeEnd : edgeStart;

      boundaryEdges.emplace_back(start, end);
    }
  }

  if (boundaryEdges.empty()) {
    return std::make_unique<Solid>();
  }

  // Build connected boundary loops from edges (handles multiple loops for holes)
  std::vector<std::vector<std::pair<Point, Point>>> allLoops;

  while (!boundaryEdges.empty()) {
    std::vector<std::pair<Point, Point>> currentLoop;
    currentLoop.push_back(boundaryEdges[0]);
    boundaryEdges.erase(boundaryEdges.begin());

    // Build connected loop
    bool loopClosed = false;
    while (!loopClosed && !boundaryEdges.empty()) {
      Point currentEnd = currentLoop.back().second;
      Point loopStart  = currentLoop.front().first;

      // Check if loop is closed (using 2D comparison)
      if (currentEnd.x() == loopStart.x() && currentEnd.y() == loopStart.y()) {
        loopClosed = true;
        break;
      }

      // Find edge starting at currentEnd
      bool foundNext = false;
      for (auto it = boundaryEdges.begin(); it != boundaryEdges.end(); ++it) {
        if (it->first.x() == currentEnd.x() &&
            it->first.y() == currentEnd.y()) {
          currentLoop.push_back(*it);
          boundaryEdges.erase(it);
          foundNext = true;
          break;
        }
      }
      if (!foundNext) {
        break;
      }
    }

    if (!currentLoop.empty()) {
      allLoops.push_back(std::move(currentLoop));
    }
  }

  if (allLoops.empty()) {
    return std::make_unique<Solid>();
  }

  // Find the exterior loop (largest area) and interior loops (holes)
  size_t exteriorIdx = 0;
  if (allLoops.size() > 1) {
    // Calculate signed area to find exterior (positive/CCW) vs holes
    // (negative/CW) For the bottom face, we need to reverse orientation
    auto calcSignedArea2D =
        [](const std::vector<std::pair<Point, Point>> &loop) -> Kernel::FT {
      Kernel::FT area(0);
      for (const auto &[start, end] : loop) {
        area += (end.x() - start.x()) * (end.y() + start.y());
      }
      return area / Kernel::FT(2);
    };

    Kernel::FT maxArea(0);
    for (size_t i = 0; i < allLoops.size(); ++i) {
      Kernel::FT area = CGAL::abs(calcSignedArea2D(allLoops[i]));
      if (area > maxArea) {
        maxArea     = area;
        exteriorIdx = i;
      }
    }
  }

  // Create bottom face with all rings (reversed for inward normal)
  // Exterior ring first, then interior rings (holes)
  LineString bottomExtRing;
  for (auto it = allLoops[exteriorIdx].rbegin();
       it != allLoops[exteriorIdx].rend(); ++it) {
    bottomExtRing.addPoint(
        Point(it->second.x(), it->second.y(), Kernel::FT(0)));
  }
  if (!bottomExtRing.isEmpty()) {
    bottomExtRing.addPoint(bottomExtRing.pointN(0)); // Close the ring
  }

  Polygon bottomFace(bottomExtRing);

  // Add interior rings (holes) - reversed again (so CW for holes in bottom)
  for (size_t i = 0; i < allLoops.size(); ++i) {
    if (i == exteriorIdx) {
      continue;
    }
    LineString holeRing;
    for (auto it = allLoops[i].rbegin(); it != allLoops[i].rend(); ++it) {
      holeRing.addPoint(Point(it->second.x(), it->second.y(), Kernel::FT(0)));
    }
    if (!holeRing.isEmpty()) {
      holeRing.addPoint(holeRing.pointN(0)); // Close the ring
      bottomFace.addRing(holeRing);
    }
  }

  shell->addPolygon(bottomFace);

  // Create walls from all boundary edges
  for (const auto &loop : allLoops) {
    for (const auto &[start, end] : loop) {
      Point startBase(start.x(), start.y(), Kernel::FT(0));
      Point endBase(end.x(), end.y(), Kernel::FT(0));

      LineString wallRing;
      wallRing.addPoint(start);
      wallRing.addPoint(startBase);
      wallRing.addPoint(endBase);
      wallRing.addPoint(end);
      wallRing.addPoint(start);
      Polygon wall(wallRing);

      shell->addPolygon(wall);
    }
  }

  auto solid = std::make_unique<Solid>();
  solid->setExteriorShell(std::move(shell));
  return solid;
}
// NOLINTEND(readability-function-cognitive-complexity)

SFCGAL_API auto
extrudeUntil(const Polygon &footprint, const Geometry &roof)
    -> std::unique_ptr<Solid>
{
  // Convert roof geometry to PolyhedralSurface based on type
  PolyhedralSurface roofSurface;

  switch (roof.geometryTypeId()) {
  case TYPE_POLYHEDRALSURFACE:
    roofSurface = roof.as<PolyhedralSurface>();
    break;

  case TYPE_POLYGON:
    roofSurface.addPolygon(roof.as<Polygon>());
    break;

  case TYPE_TRIANGLE:
    roofSurface.addPatch(roof.as<Triangle>());
    break;

  case TYPE_TRIANGULATEDSURFACE: {
    const auto &tin = roof.as<TriangulatedSurface>();
    for (size_t i = 0; i < tin.numPatches(); ++i) {
      roofSurface.addPatch(tin.patchN(i));
    }
    break;
  }

  default:
    BOOST_THROW_EXCEPTION(
        Exception("extrudeUntil: roof must be Polygon, Triangle, "
                  "PolyhedralSurface, or TriangulatedSurface"));
  }

  return extrudeUntil(footprint, roofSurface);
}

} // namespace SFCGAL::algorithm

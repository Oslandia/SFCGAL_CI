// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/split.h"

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"

#include "SFCGAL/Exception.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/algorithm/area.h"
#include "SFCGAL/algorithm/covers.h"
#include "SFCGAL/algorithm/intersection.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/detail/SegmentStore.h"

#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

#include <vector>

namespace SFCGAL::algorithm {

// ----------------------------------------------------------------------------------
// -- private interface
// ----------------------------------------------------------------------------------
/// @{
/// @privatesection

using Arr_traits  = CGAL::Arr_segment_traits_2<Kernel>;
using Arrangement = CGAL::Arrangement_2<Arr_traits>;

/**
 * @brief Insert a SFCGAL LineString into a CGAL arrangement.
 *
 * Converts each segment of the linestring and inserts it into the arrangement.
 * Handles degenerate segments (zero-length) by skipping them.
 *
 * @param arr The arrangement to insert into (modified in-place)
 * @param linestring The linestring whose segments to insert
 */
static void
insertLineStringIntoArrangement(Arrangement &arr, const LineString &linestring)
{
  if (linestring.isEmpty() || linestring.numPoints() < 2) {
    return;
  }

  for (size_t i = 0; i < linestring.numPoints() - 1; ++i) {
    const Point      &point1  = linestring.pointN(i);
    const Point      &point2  = linestring.pointN(i + 1);
    const Kernel::Point_2 arrPoint1 = point1.toPoint_2();
    const Kernel::Point_2 arrPoint2 = point2.toPoint_2();

    if (arrPoint1 == arrPoint2) {
      continue;
    }

    CGAL::insert(arr, Kernel::Segment_2(arrPoint1, arrPoint2));
  }
}

/**
 * @brief Insert a SFCGAL Polygon boundary into a CGAL arrangement.
 *
 * Inserts all edges of the exterior ring and interior rings (holes).
 * Duplicate consecutive points are skipped to avoid degenerate segments.
 *
 * @param arr The arrangement to insert into (modified in-place)
 * @param polygon The polygon whose boundary to insert
 */
static void
insertPolygonBoundaryIntoArrangement(Arrangement &arr, const Polygon &polygon)
{
  if (polygon.isEmpty()) {
    return;
  }

  const LineString &exterior    = polygon.exteriorRing();
  const size_t      numSegments = exterior.isClosed() ? exterior.numPoints() - 1
                                                      : exterior.numPoints();

  if (exterior.numPoints() >= 2) {
    for (size_t i = 0; i < numSegments; ++i) {
      const Point          &point1    = exterior.pointN(i);
      const Point          &point2    = exterior.pointN((i + 1) % exterior.numPoints());
      const Kernel::Point_2 arrPoint1 = point1.toPoint_2();
      const Kernel::Point_2 arrPoint2 = point2.toPoint_2();

      if (arrPoint1 == arrPoint2) {
        continue;
      }

      CGAL::insert(arr, Kernel::Segment_2(arrPoint1, arrPoint2));
    }
  }

  for (size_t ringIndex = 0; ringIndex < polygon.numInteriorRings(); ++ringIndex) {
    const LineString &hole       = polygon.interiorRingN(ringIndex);
    const size_t      numHoleSegs = hole.isClosed() ? hole.numPoints() - 1
                                                     : hole.numPoints();

    if (hole.numPoints() >= 2) {
      for (size_t i = 0; i < numHoleSegs; ++i) {
        const Point          &point1    = hole.pointN(i);
        const Point          &point2    = hole.pointN((i + 1) % hole.numPoints());
        const Kernel::Point_2 arrPoint1 = point1.toPoint_2();
        const Kernel::Point_2 arrPoint2 = point2.toPoint_2();

        if (arrPoint1 == arrPoint2) {
          continue;
        }

        CGAL::insert(arr, Kernel::Segment_2(arrPoint1, arrPoint2));
      }
    }
  }
}

/**
 * @brief Extract a SFCGAL Polygon from an arrangement face.
 *
 * Constructs a polygon by traversing the outer CCB (counter-clockwise boundary)
 * of the face and any inner CCBs (holes). Interpolates Z and M values from
 * the original geometry segments.
 *
 * @param face An arrangement face (must be bounded)
 * @param segmentStore Segment store for Z and M interpolation
 * @param dimension Coordinate dimension type for output points
 * @return SFCGAL polygon representing the face
 */
static auto
faceToPolygon(Arrangement::Face_const_handle      face,
              const detail::SegmentStore         &segmentStore,
              CoordinateType                      dimension) -> std::unique_ptr<Polygon>
{
  if (face->is_unbounded()) {
    return std::make_unique<Polygon>();
  }

  auto               outerCcb = face->outer_ccb();
  auto               current  = outerCcb;
  std::vector<Point> outerPoints;

  do {
    const Kernel::Point_2 &arrPoint = current->source()->point();
    const double           x        = CGAL::to_double(arrPoint.x());
    const double           y        = CGAL::to_double(arrPoint.y());
    outerPoints.push_back(segmentStore.createPoint(x, y, dimension));
    ++current;
  } while (current != outerCcb);

  if (!outerPoints.empty()) {
    outerPoints.push_back(outerPoints.front());
  }

  auto polygon = std::make_unique<Polygon>(LineString(outerPoints));

  for (auto holeIt = face->holes_begin(); holeIt != face->holes_end(); ++holeIt) {
    auto               holeCcb = *holeIt;
    auto               holeCurr = holeCcb;
    std::vector<Point> holePoints;

    do {
      const Kernel::Point_2 &arrPoint = holeCurr->source()->point();
      const double           x        = CGAL::to_double(arrPoint.x());
      const double           y        = CGAL::to_double(arrPoint.y());
      holePoints.push_back(segmentStore.createPoint(x, y, dimension));
      ++holeCurr;
    } while (holeCurr != holeCcb);

    if (!holePoints.empty()) {
      holePoints.push_back(holePoints.front());
    }

    polygon->addInteriorRing(LineString(holePoints));
  }

  return polygon;
}

/**
 * @brief Check if a polygon face is contained within the original polygon.
 *
 * Uses the centroid of the face and checks if it's covered by the original
 * polygon. This works because the arrangement ensures faces don't overlap.
 *
 * @param face The arrangement face to check
 * @param originalPolygon The original polygon before splitting
 * @return true if face is inside originalPolygon, false otherwise
 */
static auto
isFaceInsideOriginalPolygon(Arrangement::Face_const_handle face,
                            const Polygon                 &originalPolygon,
                            const detail::SegmentStore    &segmentStore,
                            CoordinateType                 dimension) -> bool
{
  if (face->is_unbounded()) {
    return false;
  }

  auto facePoly = faceToPolygon(face, segmentStore, dimension);
  if (facePoly->isEmpty()) {
    return false;
  }

  try {
    auto intersectionResult =
        intersection(*facePoly, originalPolygon, NoValidityCheck());

    if (!intersectionResult || intersectionResult->isEmpty()) {
      return false;
    }

    const double faceArea = area(*facePoly);
    const double intArea  = area(*intersectionResult);

    return std::abs(faceArea - intArea) < 1e-10;
  } catch (...) {
    return false;
  }
}

/**
 * @brief Count intersection points between linestring and polygon boundary.
 *
 * This is used to determine if the linestring actually crosses the polygon
 * (requires at least 2 intersection points for a proper split).
 *
 * @param polygon The polygon
 * @param linestring The linestring
 * @return Number of distinct intersection points
 */
static auto
countBoundaryIntersections(const Polygon &polygon, const LineString &linestring)
    -> size_t
{
  Arrangement temp_arr;
  insertPolygonBoundaryIntoArrangement(temp_arr, polygon);
  insertLineStringIntoArrangement(temp_arr, linestring);

  size_t intersection_count = 0;

  for (auto vit = temp_arr.vertices_begin(); vit != temp_arr.vertices_end();
       ++vit) {
    if (vit->degree() >= 3) {
      intersection_count++;
    }
  }

  return intersection_count;
}

/// @} end of private section

// ----------------------------------------------------------------------------------
// -- public interface
// ----------------------------------------------------------------------------------
/// @publicsection

auto
split(const Polygon &polygon, const LineString &linestring,
      NoValidityCheck /*unused*/) -> std::unique_ptr<Geometry>
{
  if (polygon.isEmpty()) {
    return std::make_unique<GeometryCollection>();
  }

  if (linestring.isEmpty() || linestring.numPoints() < 2) {
    auto collection = std::make_unique<GeometryCollection>();
    collection->addGeometry(polygon.clone().release());
    return collection;
  }

  const size_t numIntersections = countBoundaryIntersections(polygon, linestring);

  if (numIntersections < 2) {
    auto collection = std::make_unique<GeometryCollection>();
    collection->addGeometry(polygon.clone().release());
    return collection;
  }

  detail::SegmentStore segmentStore;
  segmentStore.extractSegments(polygon);
  segmentStore.extractSegments(linestring);

  CoordinateType dimension = COORDINATE_XY;
  if (segmentStore.hasZ() && segmentStore.hasM()) {
    dimension = COORDINATE_XYZM;
  } else if (segmentStore.hasZ()) {
    dimension = COORDINATE_XYZ;
  } else if (segmentStore.hasM()) {
    dimension = COORDINATE_XYM;
  }

  Arrangement arr;

  try {
    insertPolygonBoundaryIntoArrangement(arr, polygon);
    insertLineStringIntoArrangement(arr, linestring);

    std::vector<std::unique_ptr<Polygon>> resultPolygons;

    for (auto faceIt = arr.faces_begin(); faceIt != arr.faces_end(); ++faceIt) {
      if (!faceIt->is_unbounded()) {
        if (isFaceInsideOriginalPolygon(faceIt, polygon, segmentStore, dimension)) {
          auto facePoly = faceToPolygon(faceIt, segmentStore, dimension);

          if (!facePoly->isEmpty() && area(*facePoly) > 1e-10) {
            resultPolygons.push_back(std::move(facePoly));
          }
        }
      }
    }

    auto result = std::make_unique<GeometryCollection>();

    if (resultPolygons.empty() || resultPolygons.size() == 1) {
      result->addGeometry(polygon.clone().release());
    } else {
      for (auto &poly : resultPolygons) {
        result->addGeometry(poly.release());
      }
    }

    return result;

  } catch (const std::exception &) {
    auto collection = std::make_unique<GeometryCollection>();
    collection->addGeometry(polygon.clone().release());
    return collection;
  }
}

auto
split(const Polygon &polygon, const LineString &linestring)
    -> std::unique_ptr<Geometry>
{
  // Validate inputs
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(polygon);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(linestring);

  return split(polygon, linestring, NoValidityCheck());
}

auto
split(const Geometry &geometry, const LineString &linestring,
      NoValidityCheck /*unused*/) -> std::unique_ptr<Geometry>
{
  auto result = std::make_unique<GeometryCollection>();

  switch (geometry.geometryTypeId()) {
  case TYPE_POLYGON:
    return split(geometry.as<Polygon>(), linestring, NoValidityCheck());

  case TYPE_TRIANGLE: {
    const auto         &triangle = geometry.as<Triangle>();
    std::vector<Point> points    = {triangle.vertex(0), triangle.vertex(1),
                                 triangle.vertex(2), triangle.vertex(0)};
    LineString         ring(points);
    Polygon            polygon(ring);
    return split(polygon, linestring, NoValidityCheck());
  }

  case TYPE_MULTIPOLYGON: {
    const auto &multiPolygon = geometry.as<MultiPolygon>();
    for (size_t i = 0; i < multiPolygon.numGeometries(); ++i) {
      auto splitResult =
          split(multiPolygon.polygonN(i), linestring, NoValidityCheck());

      if (splitResult && !splitResult->isEmpty()) {
        const auto &gc = splitResult->as<GeometryCollection>();
        for (size_t j = 0; j < gc.numGeometries(); ++j) {
          result->addGeometry(gc.geometryN(j).clone().release());
        }
      }
    }
    break;
  }

  case TYPE_POLYHEDRALSURFACE: {
    const auto &surface = geometry.as<PolyhedralSurface>();
    for (size_t i = 0; i < surface.numPolygons(); ++i) {
      auto splitResult =
          split(surface.polygonN(i), linestring, NoValidityCheck());

      if (splitResult && !splitResult->isEmpty()) {
        const auto &gc = splitResult->as<GeometryCollection>();
        for (size_t j = 0; j < gc.numGeometries(); ++j) {
          result->addGeometry(gc.geometryN(j).clone().release());
        }
      }
    }
    break;
  }

  case TYPE_TRIANGULATEDSURFACE: {
    const auto &surface = geometry.as<TriangulatedSurface>();
    for (size_t i = 0; i < surface.numTriangles(); ++i) {
      const auto         &triangle = surface.triangleN(i);
      std::vector<Point> points    = {triangle.vertex(0), triangle.vertex(1),
                                   triangle.vertex(2), triangle.vertex(0)};
      LineString         ring(points);
      Polygon            polygon(ring);
      auto               splitResult = split(polygon, linestring, NoValidityCheck());

      if (splitResult && !splitResult->isEmpty()) {
        const auto &gc = splitResult->as<GeometryCollection>();
        for (size_t j = 0; j < gc.numGeometries(); ++j) {
          result->addGeometry(gc.geometryN(j).clone().release());
        }
      }
    }
    break;
  }

  case TYPE_GEOMETRYCOLLECTION: {
    const auto &collection = geometry.as<GeometryCollection>();
    for (size_t i = 0; i < collection.numGeometries(); ++i) {
      auto splitResult =
          split(collection.geometryN(i), linestring, NoValidityCheck());

      if (splitResult && !splitResult->isEmpty()) {
        const auto &gc = splitResult->as<GeometryCollection>();
        for (size_t j = 0; j < gc.numGeometries(); ++j) {
          result->addGeometry(gc.geometryN(j).clone().release());
        }
      }
    }
    break;
  }

  default:
    BOOST_THROW_EXCEPTION(Exception(
        "split: unsupported geometry type " + geometry.geometryType()));
  }

  if (result->isEmpty()) {
    result->addGeometry(geometry.clone().release());
  }

  return result;
}

auto
split(const Geometry &geometry, const LineString &linestring)
    -> std::unique_ptr<Geometry>
{
  // Validate inputs
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(geometry);
  SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(linestring);

  return split(geometry, linestring, NoValidityCheck());
}

} // namespace SFCGAL::algorithm

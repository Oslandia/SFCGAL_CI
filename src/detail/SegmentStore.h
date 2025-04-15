// Copyright (c) 2025-2025, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_DETAIL_SEGMENTSTORE_H_
#define SFCGAL_DETAIL_SEGMENTSTORE_H_

#include <vector>
#include <limits>
#include <cmath>

#include "SFCGAL/Segment.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/numeric.h"

namespace SFCGAL {
namespace detail {

/**
 * @brief Collection of segments from a geometry for interpolation
 */
class SegmentStore {
private:
  std::vector<Segment> segments;
  bool                 hasZCoord = false;
  bool                 hasMCoord = false;

public:
  SegmentStore() {}

  void
  addSegment(const Segment &segment)
  {
    segments.push_back(segment);

    // Update dimension flags
    hasZCoord = segment.source().is3D() && segment.target().is3D();
    hasMCoord = segment.source().isMeasured() && segment.target().isMeasured();
  }

  // Check if store contains segments with Z coordinates
  bool
  hasZ() const
  {
    return hasZCoord;
  }

  // Check if store contains segments with M values
  bool
  hasM() const
  {
    return hasMCoord;
  }

  // Find the nearest segment to a point
  Segment
  findNearestSegment(double x, double y) const
  {
    if (segments.empty()) {
      // Return a dummy segment if store is empty
      return Segment(Point(), Point());
    }

    double minDist = std::numeric_limits<double>::max();
    size_t bestIdx = 0;

    for (size_t i = 0; i < segments.size(); ++i) {
      double dist = segments[i].distanceToPoint(x, y);
      if (dist < minDist) {
        minDist = dist;
        bestIdx = i;
      }
    }

    return segments[bestIdx];
  }

  // Interpolate Z value for a point
  double
  interpolateZ(double x, double y) const
  {
    if (segments.empty() || !hasZCoord) {
      return NaN();
    }

    Segment nearest = findNearestSegment(x, y);
    double  t       = nearest.interpolationParameter(x, y);

    if (!nearest.source().is3D() || !nearest.target().is3D()) {
      return NaN();
    }

    double z1 = CGAL::to_double(nearest.source().z());
    double z2 = CGAL::to_double(nearest.target().z());

    return z1 + t * (z2 - z1);
  }

  // Interpolate M value for a point
  double
  interpolateM(double x, double y) const
  {
    if (segments.empty() || !hasMCoord) {
      return NaN();
    }

    Segment nearest = findNearestSegment(x, y);
    double  t       = nearest.interpolationParameter(x, y);

    if (!nearest.source().isMeasured() || !nearest.target().isMeasured()) {
      return NaN();
    }

    double m1 = nearest.source().m();
    double m2 = nearest.target().m();

    return m1 + t * (m2 - m1);
  }
};

/**
 * @brief Extract segments from a LineString for interpolation
 */
void
extractSegments(const LineString &lineString, SegmentStore &store)
{
  if (lineString.numPoints() < 2) {
    return;
  }

  for (size_t i = 1; i < lineString.numPoints(); ++i) {
    const Point &p1 = lineString.pointN(i - 1);
    const Point &p2 = lineString.pointN(i);
    store.addSegment(Segment(p1, p2));
  }
}

/**
 * @brief Extract segments from a Polygon for interpolation
 */
void
extractSegments(const Polygon &polygon, SegmentStore &store)
{
  // Extract from exterior ring
  extractSegments(polygon.exteriorRing(), store);

  // Extract from interior rings
  for (size_t i = 0; i < polygon.numInteriorRings(); ++i) {
    extractSegments(polygon.interiorRingN(i), store);
  }
}

/**
 * @brief Extract segments from all geometry types
 */
void
extractSegments(const Geometry &geometry, SegmentStore &store)
{
  switch (geometry.geometryTypeId()) {
  case TYPE_POINT:
    // Points don't have segments
    break;

  case TYPE_LINESTRING:
    extractSegments(static_cast<const LineString &>(geometry), store);
    break;

  case TYPE_POLYGON:
    extractSegments(static_cast<const Polygon &>(geometry), store);
    break;

  case TYPE_MULTIPOINT:
    // MultiPoints don't have segments
    break;

  case TYPE_MULTILINESTRING: {
    const auto &multiLine = static_cast<const MultiLineString &>(geometry);
    for (size_t i = 0; i < multiLine.numGeometries(); ++i) {
      extractSegments(multiLine.lineStringN(i), store);
    }
    break;
  }

  case TYPE_MULTIPOLYGON: {
    const auto &multiPolygon = static_cast<const MultiPolygon &>(geometry);
    for (size_t i = 0; i < multiPolygon.numGeometries(); ++i) {
      extractSegments(multiPolygon.polygonN(i), store);
    }
    break;
  }

  case TYPE_GEOMETRYCOLLECTION: {
    const auto &collection = static_cast<const GeometryCollection &>(geometry);
    for (size_t i = 0; i < collection.numGeometries(); ++i) {
      extractSegments(collection.geometryN(i), store);
    }
    break;
  }

  case TYPE_POLYHEDRALSURFACE: {
    const auto &surface = static_cast<const PolyhedralSurface &>(geometry);
    for (size_t i = 0; i < surface.numPolygons(); ++i) {
      extractSegments(surface.polygonN(i), store);
    }
    break;
  }

  default:
    // Other types not supported
    break;
  }
}

/**
 * @brief Create a point with interpolated Z and M values
 */
Point
createPoint(double x, double y, double z,
                              const SegmentStore &store,
                              CoordinateType dimension)
{
  // Determine if we need Z and/or M values
  bool needsZ = (dimension == CoordinateType::COORDINATE_XYZ ||
                 dimension == CoordinateType::COORDINATE_XYZM);

  bool needsM = (dimension == CoordinateType::COORDINATE_XYM ||
                 dimension == CoordinateType::COORDINATE_XYZM);

  // If no interpolation needed, return point with original values
  if (!needsZ && !needsM) {
    return Point(x, y);
  }

  double finalZ = z;
  double finalM = NaN();

  // Interpolate Z if needed
  if (needsZ && store.hasZ()) {
    double interpolatedZ = store.interpolateZ(x, y);
    if (!std::isnan(interpolatedZ)) {
      finalZ = interpolatedZ;
    }
  }

  // Interpolate M if needed
  if (needsM && store.hasM()) {
    finalM = store.interpolateM(x, y);
  }

  // Create point with appropriate dimensions
  if (needsZ && needsM) {
    return Point(x, y, finalZ, finalM);
  } else if (needsZ) {
    return Point(x, y, finalZ);
  } else if (needsM) {
    Point p(x, y);
    p.setM(finalM);
    return p;
  } else {
    return Point(x, y);
  }
}

} // namespace detail
} // namespace SFCGAL

#endif // SFCGAL_DETAIL_SEGMENTSTORE_H_

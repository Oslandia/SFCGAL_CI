// Copyright (c) 2025-2025, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_DETAIL_SEGMENTSTORE_H_
#define SFCGAL_DETAIL_SEGMENTSTORE_H_

#include <cmath>
#include <limits>
#include <tuple>
#include <vector>

#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPoint.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Segment.h"
#include "SFCGAL/numeric.h"

namespace SFCGAL::detail {

/**
 * @brief Collection of segments from a geometry for interpolation
 */
class SegmentStore {
private:
  std::vector<Segment>
       segments;          ///< Collection of segments stored for interpolation
  bool hasZCoord = false; ///< Flag indicating if any segment has Z coordinates
  bool hasMCoord = false; ///< Flag indicating if any segment has M values

public:
  /**
   * @brief Default constructor
   */
  SegmentStore() {}

  /**
   * @brief Add a segment to the store
   * @param segment The segment to add
   */
  void
  addSegment(const Segment &segment)
  {
    segments.push_back(segment);

    // Update dimension flags
    hasZCoord =
        hasZCoord || (segment.source().is3D() && segment.target().is3D());
    hasMCoord = hasMCoord || (segment.source().isMeasured() &&
                              segment.target().isMeasured());
  }

  /**
   * @brief Check if store contains segments with Z coordinates
   * @return True if any segment has Z coordinates
   */
  bool
  hasZ() const
  {
    return hasZCoord;
  }

  /**
   * @brief Check if store contains segments with M values
   * @return True if any segment has M values
   */
  bool
  hasM() const
  {
    return hasMCoord;
  }

  /**
   * @brief Find the nearest segment to a point
   * @param x X-coordinate of the point
   * @param y Y-coordinate of the point
   * @return The nearest segment to the given point
   */
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

  /**
   * @brief Interpolate Z and M values for a point
   * @param x X-coordinate of the point
   * @param y Y-coordinate of the point
   * @return Tuple with interpolated (z, m) values, NaN if not applicable
   */
  std::tuple<double, double>
  interpolateZM(double x, double y) const
  {
    double z = NaN();
    double m = NaN();

    if (segments.empty()) {
      return std::make_tuple(z, m);
    }

    Segment nearest = findNearestSegment(x, y);
    double  t       = nearest.interpolationParameter(x, y);

    // Interpolate Z if available
    if (hasZCoord && nearest.source().is3D() && nearest.target().is3D()) {
      double z1 = CGAL::to_double(nearest.source().z());
      double z2 = CGAL::to_double(nearest.target().z());
      z         = z1 + t * (z2 - z1);
    }

    // Interpolate M if available
    if (hasMCoord && nearest.source().isMeasured() &&
        nearest.target().isMeasured()) {
      double m1 = nearest.source().m();
      double m2 = nearest.target().m();
      m         = m1 + t * (m2 - m1);
    }

    return std::make_tuple(z, m);
  }

  /**
   * @brief Extract segments from a LineString for interpolation
   * @param lineString The linestring to extract segments from
   */
  void
  extractSegments(const LineString &lineString)
  {
    if (lineString.numPoints() < 2) {
      return;
    }

    for (size_t i = 1; i < lineString.numPoints(); ++i) {
      const Point &p1 = lineString.pointN(i - 1);
      const Point &p2 = lineString.pointN(i);
      addSegment(Segment(p1, p2));
    }
  }

  /**
   * @brief Extract segments from a Polygon for interpolation
   * @param polygon The polygon to extract segments from
   */
  void
  extractSegments(const Polygon &polygon)
  {
    if (polygon.isEmpty()) {
      return;
    }

    // Extract from exterior ring
    extractSegments(polygon.exteriorRing());

    // Extract from interior rings
    for (size_t i = 0; i < polygon.numInteriorRings(); ++i) {
      extractSegments(polygon.interiorRingN(i));
    }
  }

  /**
   * @brief Extract segments from all geometry types
   * @param geometry The geometry to extract segments from
   */
  void
  extractSegments(const Geometry &geometry)
  {
    switch (geometry.geometryTypeId()) {
    case TYPE_POINT:
      // Points don't have segments
      break;

    case TYPE_LINESTRING:
      extractSegments(static_cast<const LineString &>(geometry));
      break;

    case TYPE_POLYGON:
      extractSegments(static_cast<const Polygon &>(geometry));
      break;

    case TYPE_MULTIPOINT:
      // MultiPoints don't have segments
      break;

    case TYPE_MULTILINESTRING: {
      const auto &multiLine = static_cast<const MultiLineString &>(geometry);
      for (size_t i = 0; i < multiLine.numGeometries(); ++i) {
        extractSegments(multiLine.lineStringN(i));
      }
      break;
    }

    case TYPE_MULTIPOLYGON: {
      const auto &multiPolygon = static_cast<const MultiPolygon &>(geometry);
      for (size_t i = 0; i < multiPolygon.numGeometries(); ++i) {
        extractSegments(multiPolygon.polygonN(i));
      }
      break;
    }

    case TYPE_GEOMETRYCOLLECTION: {
      const auto &collection =
          static_cast<const GeometryCollection &>(geometry);
      for (size_t i = 0; i < collection.numGeometries(); ++i) {
        extractSegments(collection.geometryN(i));
      }
      break;
    }

    case TYPE_POLYHEDRALSURFACE: {
      const auto &surface = static_cast<const PolyhedralSurface &>(geometry);
      for (size_t i = 0; i < surface.numPatches(); ++i) {
        extractSegments(surface.patchN(i));
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
   * @param x X-coordinate of the point
   * @param y Y-coordinate of the point
   * @param dimension The coordinate type for the point
   * @return Point with interpolated Z and M values
   */
  Point
  createPoint(double x, double y, CoordinateType dimension) const
  {
    auto [interpZ, interpM] = interpolateZM(x, y);
    return Point(x, y, interpZ, interpM, dimension);
  }
};

} // namespace SFCGAL::detail

#endif // SFCGAL_DETAIL_SEGMENTSTORE_H_

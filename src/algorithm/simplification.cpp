// Copyright (c) 2025-2025, SFCGAL Team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/simplification.h"
#include "SFCGAL/GeometryCollection.h"
#include "SFCGAL/LineString.h"
#include "SFCGAL/MultiLineString.h"
#include "SFCGAL/MultiPolygon.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/TriangulatedSurface.h"
#include "SFCGAL/algorithm/isValid.h"
#include "SFCGAL/detail/SegmentStore.h"
#include "SFCGAL/detail/algorithm/simplification.h"

namespace SFCGAL {
namespace algorithm {

/**
 * @brief Main entry point for geometry simplification
 */
auto
simplify(const Geometry &geometry, double threshold, bool preserveTopology)
    -> std::unique_ptr<Geometry>
{
  switch (geometry.geometryTypeId()) {
  case TYPE_SOLID: {
    // Is it make sense?
    break;
  }
  default:
    SFCGAL_ASSERT_GEOMETRY_VALIDITY_2D(geometry);
    break;
  }

  std::unique_ptr<Geometry> result(
      simplify(geometry, threshold, preserveTopology, NoValidityCheck()));
  propagateValidityFlag(*result, true);
  return result;
}

auto
simplify(const Geometry &geometry, double threshold, bool preserveTopology,
         NoValidityCheck /*unused*/) -> std::unique_ptr<Geometry>
{
  if (geometry.isEmpty()) {
    return std::unique_ptr<Geometry>(geometry.clone());
  }

  const double squaredThreshold = threshold * threshold;

  switch (geometry.geometryTypeId()) {
  case TYPE_LINESTRING: {
    // Extract segments for interpolation
    detail::SegmentStore store;
    store.extractSegments(geometry);

    return detail::simplifyLineString(static_cast<const LineString &>(geometry),
                                      squaredThreshold, preserveTopology,
                                      store);
  }

  case TYPE_MULTILINESTRING:
    return detail::simplifyMultiLineString(
        static_cast<const MultiLineString &>(geometry), squaredThreshold,
        preserveTopology);

  case TYPE_POLYGON:
    return detail::simplifyPolygon(static_cast<const Polygon &>(geometry),
                                   squaredThreshold, preserveTopology);

  case TYPE_MULTIPOLYGON:
    return detail::simplifyMultiPolygon(
        static_cast<const MultiPolygon &>(geometry), squaredThreshold,
        preserveTopology);

  case TYPE_POLYHEDRALSURFACE:
    return detail::simplifyPolyhedralSurface(
        static_cast<const PolyhedralSurface &>(geometry), squaredThreshold,
        preserveTopology);

  case TYPE_GEOMETRYCOLLECTION:
    return detail::simplifyGeometryCollection(
        static_cast<const GeometryCollection &>(geometry), squaredThreshold,
        preserveTopology);

  default:
    // For unsupported types, return a copy
    return std::unique_ptr<Geometry>(geometry.clone());
  }
}

} // namespace algorithm
} // namespace SFCGAL

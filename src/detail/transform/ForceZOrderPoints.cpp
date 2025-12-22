// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/transform/ForceZOrderPoints.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/algorithm/orientation.h"

#include <utility>

namespace SFCGAL::transform {

ForceZOrderPoints::ForceZOrderPoints(Kernel::FT defaultZ)
    : _defaultZ(std::move(defaultZ))
{
}

void
ForceZOrderPoints::transform(Point &point)
{
  if (point.isEmpty()) {
    return;
  }

  if (!point.is3D()) {
    point = Point(point.x(), point.y(), _defaultZ);
  }
}

void
ForceZOrderPoints::visit(Triangle &triangle)
{
  if (triangle.isEmpty()) {
    return;
  }

  if (!triangle.is3D()) {
    if (!SFCGAL::algorithm::isCounterClockWiseOriented(triangle)) {
      // not pointing up, reverse
      triangle.reverse();
    }

    Transform::visit(triangle);
  }
}

void
ForceZOrderPoints::visit(Polygon &polygon)
{
  if (polygon.isEmpty()) {
    return;
  }

  if (!polygon.is3D()) {
    LineString &ext = polygon.exteriorRing();

    if (!SFCGAL::algorithm::isCounterClockWiseOriented(
            polygon.exteriorRing())) {
      // exterior ring not pointing up, reverse
      ext.reverse();
    }

    for (size_t i = 0; i < polygon.numInteriorRings(); ++i) {
      LineString &inter = polygon.interiorRingN(i);

      if (SFCGAL::algorithm::isCounterClockWiseOriented(inter)) {
        // interior ring is pointing up, reverse
        inter.reverse();
      }
    }

    Transform::visit(polygon);
  }
}

} // namespace SFCGAL::transform

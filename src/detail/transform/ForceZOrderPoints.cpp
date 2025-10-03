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
ForceZOrderPoints::transform(Point &p)
{
  if (!p.is3D()) {
    p = Point(p.x(), p.y(), _defaultZ);
  }
}

void
ForceZOrderPoints::visit(Triangle &t)
{
  if (!t.is3D()) {
    if (!SFCGAL::algorithm::isCounterClockWiseOriented(t)) {
      // not pointing up, reverse
      t.reverse();
    }

    Transform::visit(t);
  }
}

void
ForceZOrderPoints::visit(Polygon &p)
{
  if (!p.is3D()) {
    LineString &ext = p.exteriorRing();

    if (!SFCGAL::algorithm::isCounterClockWiseOriented(p.exteriorRing())) {
      // exterior ring not pointing up, reverse
      ext.reverse();
    }

    for (size_t i = 0; i < p.numInteriorRings(); ++i) {
      LineString &inter = p.interiorRingN(i);

      if (SFCGAL::algorithm::isCounterClockWiseOriented(inter)) {
        // interior ring is pointing up, reverse
        inter.reverse();
      }
    }

    Transform::visit(p);
  }
}

} // namespace SFCGAL::transform

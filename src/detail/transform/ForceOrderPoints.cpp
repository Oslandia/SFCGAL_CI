// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/transform/ForceOrderPoints.h"
#include "SFCGAL/Point.h"
#include "SFCGAL/Polygon.h"
#include "SFCGAL/Triangle.h"
#include "SFCGAL/algorithm/orientation.h"

namespace SFCGAL::transform {

ForceOrderPoints::ForceOrderPoints(bool orientCCW) : _orientCCW(orientCCW) {}

void
ForceOrderPoints::transform(Point & /*point*/)
{
}

void
ForceOrderPoints::visit(Triangle &t)
{
  if (!SFCGAL::algorithm::isCounterClockWiseOriented(t)) {
    // not pointing up, reverse
    if (_orientCCW) {
      t.reverse();
    }
  } else {
    if (!_orientCCW) {
      t.reverse();
    }
  }

  Transform::visit(t);
}

void
ForceOrderPoints::visit(Polygon &p)
{
  LineString &ext = p.exteriorRing();

  if (!SFCGAL::algorithm::isCounterClockWiseOriented(p.exteriorRing())) {
    // exterior ring not pointing up, reverse
    if (_orientCCW) {
      ext.reverse();
    }
  } else {
    if (!_orientCCW) {
      ext.reverse();
    }
  }
  const bool isCCWO{SFCGAL::algorithm::isCounterClockWiseOriented(ext)};
  for (size_t i = 0; i < p.numInteriorRings(); ++i) {
    LineString &inter = p.interiorRingN(i);

    if (SFCGAL::algorithm::isCounterClockWiseOriented(inter) == isCCWO) {
      inter.reverse();
    }
  }

  Transform::visit(p);
}

} // namespace SFCGAL::transform

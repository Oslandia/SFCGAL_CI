// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/algorithm/orientation.h>
#include <SFCGAL/detail/transform/ForceOrderPoints.h>

namespace SFCGAL {
namespace transform {

///
///
///
ForceOrderPoints::ForceOrderPoints(bool orientCCW) : _orientCCW(orientCCW) {}

///
///
///
void
ForceOrderPoints::transform(Point &)
{
}

///
///
///
void
ForceOrderPoints::visit(Triangle &t)
{
  if (!t.is3D()) {
    if (!algorithm::isCounterClockWiseOriented(t)) {
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
}

void
ForceOrderPoints::visit(Polygon &p)
{
  if (!p.is3D()) {
    LineString &ext = p.exteriorRing();

    if (!algorithm::isCounterClockWiseOriented(p.exteriorRing())) {
      // exterior ring not pointing up, reverse
      if (_orientCCW) {
        ext.reverse();
      }
    } else {
      if (!_orientCCW) {
        ext.reverse();
      }
    }

    for (size_t i = 0; i < p.numInteriorRings(); ++i) {
      LineString inter = p.interiorRingN(i);

      if (algorithm::isCounterClockWiseOriented(inter)) {
        // interior ring is pointing up, reverse
        if (_orientCCW) {
          inter.reverse();
        }
      } else {
        if (!_orientCCW) {
          inter.reverse();
        }
      }
    }

    Transform::visit(p);
  }
}

} // namespace transform
} // namespace SFCGAL

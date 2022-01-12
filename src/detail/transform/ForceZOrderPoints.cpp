// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/Point.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/algorithm/orientation.h>
#include <SFCGAL/detail/transform/ForceZOrderPoints.h>

#include <utility>

namespace SFCGAL {
namespace transform {

///
///
///
ForceZOrderPoints::ForceZOrderPoints(Kernel::FT defaultZ)
    : _defaultZ(std::move(defaultZ))
{
}

///
///
///
void
ForceZOrderPoints::transform(Point &p)
{
  if (!p.is3D()) {
    p = Point(p.x(), p.y(), _defaultZ);
  }
}

///
///
///
void
ForceZOrderPoints::visit(Triangle &t)
{
  if (!t.is3D()) {
    if (!algorithm::isCounterClockWiseOriented(t)) {
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

    if (!algorithm::isCounterClockWiseOriented(p.exteriorRing())) {
      // exterior ring not pointing up, reverse
      ext.reverse();
    }

    for (size_t i = 0; i < p.numInteriorRings(); ++i) {
      LineString &inter = p.interiorRingN(i);

      if (algorithm::isCounterClockWiseOriented(inter)) {
        // interior ring is pointing up, reverse
        inter.reverse();
      }
    }

    Transform::visit(p);
  }
}

} // namespace transform
} // namespace SFCGAL

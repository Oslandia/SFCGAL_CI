// Copyright (c) 2025-2025, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/transform/ForceM.h"
#include "SFCGAL/Point.h"

#include <utility>

namespace SFCGAL::transform {

ForceM::ForceM(const double &defaultM) : _defaultM(defaultM) {}

void
ForceM::transform(Point &point)
{
  if (!point.isEmpty() && !point.isMeasured()) {
    if (point.is3D()) {
      point = Point(point.x(), point.y(), point.z(), _defaultM);
    } else {
      point = Point(point.x(), point.y());
      point.setM(_defaultM);
    }
  }
}

} // namespace SFCGAL::transform

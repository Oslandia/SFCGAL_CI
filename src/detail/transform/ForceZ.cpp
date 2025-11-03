// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/transform/ForceZ.h"
#include "SFCGAL/Point.h"

#include <utility>

namespace SFCGAL::transform {

ForceZ::ForceZ(Kernel::FT defaultZ) : _defaultZ(std::move(defaultZ)) {}

void
ForceZ::transform(Point &point)
{
  if (!point.isEmpty() && !point.is3D()) {
    Point pt(point.x(), point.y(), _defaultZ);
    if (point.isMeasured()) {
      pt.setM(point.m());
    }
    point = pt;
  }
}

} // namespace SFCGAL::transform

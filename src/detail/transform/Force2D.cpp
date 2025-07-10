// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/detail/transform/Force2D.h"
#include "SFCGAL/Point.h"

namespace SFCGAL::transform {

void
Force2D::transform(Point &p)
{
  if (!p.isEmpty() && p.is3D()) {
    Point pt(p.x(), p.y());
    if (p.isMeasured()) {
      pt.setM(p.m());
    }
    p = pt;
  }
}

} // namespace SFCGAL::transform

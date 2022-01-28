// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include <SFCGAL/Point.h>
#include <SFCGAL/detail/transform/Force2D.h>

namespace SFCGAL {
namespace transform {

///
///
///
void
Force2D::transform(Point &p)
{
  if (!p.isEmpty() && p.is3D()) {
    Point pt(p.x(), p.y());
    if (p.isMeasured())
      pt.setM(p.m());
    p = pt;
  }
}

} // namespace transform
} // namespace SFCGAL

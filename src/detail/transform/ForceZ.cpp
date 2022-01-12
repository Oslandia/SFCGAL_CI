// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/Point.h>
#include <SFCGAL/detail/transform/ForceZ.h>

#include <utility>

namespace SFCGAL {
namespace transform {

///
///
///
ForceZ::ForceZ(const Kernel::FT defaultZ) : _defaultZ(std::move(defaultZ)) {}

///
///
///
void
ForceZ::transform(Point &p)
{
  if (!p.isEmpty() && !p.is3D()) {
    Point pt(p.x(), p.y(), _defaultZ);
    if (p.isMeasured())
      pt.setM(p.m());
    p = pt;
  }
}

} // namespace transform
} // namespace SFCGAL

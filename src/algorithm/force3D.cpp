// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#include <SFCGAL/algorithm/force3D.h>
#include <SFCGAL/detail/transform/ForceZ.h>

namespace SFCGAL {
namespace algorithm {

///
///
///
void
force3D(Geometry &g, const Kernel::FT &defaultZ)
{
  transform::ForceZ t(defaultZ);
  g.accept(t);
}

} // namespace algorithm
} // namespace SFCGAL

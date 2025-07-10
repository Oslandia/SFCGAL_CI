// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/force3D.h"
#include "SFCGAL/detail/transform/ForceZ.h"

namespace SFCGAL::algorithm {

void
force3D(Geometry &g, const Kernel::FT &defaultZ)
{
  transform::ForceZ t(defaultZ);
  g.accept(t);
}

} // namespace SFCGAL::algorithm

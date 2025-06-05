// Copyright (c) 2025-2025, Oslandia.
// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#include "SFCGAL/algorithm/forceMeasured.h"
#include "SFCGAL/detail/transform/ForceM.h"

namespace SFCGAL::algorithm {

void
forceMeasured(Geometry &geometry, const double &defaultM)
{
  transform::ForceM transform(defaultM);
  geometry.accept(transform);
}

} // namespace SFCGAL::algorithm

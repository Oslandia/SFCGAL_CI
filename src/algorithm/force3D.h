// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_FORCE3D_H_
#define SFCGAL_ALGORITHM_FORCE3D_H_

#include "SFCGAL/config.h"

#include "SFCGAL/Kernel.h"

namespace SFCGAL {
class Geometry;
}

namespace SFCGAL::algorithm {

/**
 * @brief force a 2D geometry to be 3D (replace undefined Z by defaultZ,
 * existing Z values remains unchanged)
 * @param g the geometry to force to 3D
 * @param defaultZ the default Z value to use (default: 0)
 * @warning ignore empty geometries
 */
SFCGAL_API void
force3D(Geometry &g, const Kernel::FT &defaultZ = 0);

} // namespace SFCGAL::algorithm

#endif

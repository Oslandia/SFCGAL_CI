// Copyright (c) 2025-2025, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_FORCE_MEASURED_H_
#define SFCGAL_ALGORITHM_FORCE_MEASURED_H_

#include "SFCGAL/config.h"

namespace SFCGAL {
class Geometry;
}

namespace SFCGAL::algorithm {

/**
 * @brief force a 2D or M geometry to be M (replace undefined M by defaultM,
 * existing M values remains unchanged)
 * @warning ignore empty geometries
 */
SFCGAL_API void
forceMeasured(Geometry &geometry, const double &defaultM = 0);

} // namespace SFCGAL::algorithm

#endif

// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_COVERS_POINTS_ALGORITHM
#define SFCGAL_COVERS_POINTS_ALGORITHM

#include "SFCGAL/config.h"

namespace SFCGAL {
class Geometry;
namespace detail {
namespace algorithm {

/**
 * Pseudo cover test on 2D geometries. Collect points of gb and tests if no
 * points of gb is outside ga
 * @param ga The covering geometry
 * @param gb The geometry to test coverage for
 * @return True if ga covers all points of gb
 */
SFCGAL_API auto
coversPoints(const Geometry &ga, const Geometry &gb) -> bool;

/**
 * Pseudo cover test on 3D geometries. Collect points of gb and tests if no
 * points of gb is outside ga
 * @param ga The covering geometry
 * @param gb The geometry to test coverage for
 * @return True if ga covers all points of gb in 3D
 */
SFCGAL_API auto
coversPoints3D(const Geometry &ga, const Geometry &gb) -> bool;
} // namespace algorithm
} // namespace detail
} // namespace SFCGAL

#endif

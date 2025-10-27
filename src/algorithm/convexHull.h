// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_CONVEXHULL_H_
#define SFCGAL_ALGORITHM_CONVEXHULL_H_

#include "SFCGAL/config.h"

#include "SFCGAL/Geometry.h"

namespace SFCGAL::algorithm {

/**
 * Compute the 2D convex hull for a geometry
 * @param geometry the input geometry
 * @return the 2D convex hull as a unique_ptr<Geometry>
 */
SFCGAL_API std::unique_ptr<Geometry>
           convexHull(const Geometry &geometry);

/**
 * Compute the 3D convex hull for a geometry
 * @param geometry the input geometry
 * @return the 3D convex hull as a unique_ptr<Geometry>
 * @todo improve to handle collinear points and coplanar points
 */
SFCGAL_API std::unique_ptr<Geometry>
           convexHull3D(const Geometry &geometry);

} // namespace SFCGAL::algorithm

#endif

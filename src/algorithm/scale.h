// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_SCALE_H
#define SFCGAL_ALGORITHM_SCALE_H

#include "SFCGAL/Geometry.h"
#include "SFCGAL/config.h"

namespace SFCGAL::algorithm {

/**
 * Scale a geometry by a given factor
 * @param geometry input geometry
 * @param scaleFactor scale factor
 */
SFCGAL_API void
scale(Geometry &geometry, double scaleFactor);

/**
 * Scale a geometry by different factors for each dimension
 * @param geometry input geometry
 * @param sx scale factor for x dimension
 * @param sy scale factor for y dimension
 * @param sz scale factor for z dimension
 */
SFCGAL_API void
scale(Geometry &geometry, double sx, double sy, double sz = 0.0);

/**
 * Scale a geometry by different factors for each dimension around a center
 * point
 * @param geometry input geometry
 * @param sx scale factor for x dimension
 * @param sy scale factor for y dimension
 * @param sz scale factor for z dimension
 * @param cx x-coordinate of the center point
 * @param cy y-coordinate of the center point
 * @param cz z-coordinate of the center point
 */
SFCGAL_API void
scale(Geometry &geometry, double sx, double sy, double sz, double cx, double cy,
      double cz);

} // namespace SFCGAL::algorithm

#endif

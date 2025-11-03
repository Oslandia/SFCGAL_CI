// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_GENERATOR_DISC_H_
#define SFCGAL_GENERATOR_DISC_H_

#include "SFCGAL/config.h"

#include <memory>

namespace SFCGAL {
class Point;
class Polygon;
} // namespace SFCGAL

namespace SFCGAL {
namespace generator {

/**
 * @brief Generate a discrete circle
 * @param center The center point of the circle
 * @param radius The radius of the circle
 * @param nQuadrantSegments Number of segments per quadrant
 * @return A unique pointer to the generated polygon
 * @todo unittest
 */
SFCGAL_API auto
disc(const Point &center, const double &radius,
     const unsigned int &nQuadrantSegments = 8U) -> std::unique_ptr<Polygon>;

} // namespace generator
} // namespace SFCGAL

#endif

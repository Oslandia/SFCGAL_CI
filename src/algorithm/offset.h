// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_OFFSET_H_
#define SFCGAL_ALGORITHM_OFFSET_H_

#include "SFCGAL/config.h"

#include <memory>

namespace SFCGAL {
class Geometry;
class MultiPolygon;
} // namespace SFCGAL

namespace SFCGAL::algorithm {
struct NoValidityCheck;

/**
 * @brief [experimental]compute polygon offset
 *
 * @param geometry The input geometry
 * @param radius The offset radius
 * @return The offset polygon as a MultiPolygon
 * @warning test in order to compare with minkowski sum
 * @pre geometry is a valid Geometry
 * @note When applied to NURBSCurve geometries, the offset
 *   is internally computed on a LineString obtained via
 *   toLineString() with its default parameters.
 */
SFCGAL_API std::unique_ptr<MultiPolygon>
           offset(const Geometry &geometry, const double &radius);

/**
 * @brief [experimental]compute polygon offset
 *
 * @param geometry The input geometry
 * @param radius The offset radius
 * @return The offset polygon as a MultiPolygon
 * @warning test in order to compare with minkowski sum
 * @pre geometry is a valid Geometry
 * @warning No actual validity check is done.
 */
SFCGAL_API std::unique_ptr<MultiPolygon>
offset(const Geometry &geometry, const double &radius, NoValidityCheck);

} // namespace SFCGAL::algorithm

#endif

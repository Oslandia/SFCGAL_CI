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

namespace SFCGAL {
namespace algorithm {
struct NoValidityCheck;

/**
 * @brief [experimental]compute polygon offset
 *
 * @warning test in order to compare with minkowski sum
 * @pre g is a valid Geometry
 * @note When applied to NURBSCurve geometries, the offset
 *   is internally computed on a LineString obtained via
 *   toLineString() with its default parameters.
 */
SFCGAL_API std::unique_ptr<MultiPolygon>
           offset(const Geometry &g, const double &r);

/**
 * @brief [experimental]compute polygon offset
 *
 * @warning test in order to compare with minkowski sum
 * @pre g is a valid Geometry
 * @warning No actual validity check is done.
 */
SFCGAL_API std::unique_ptr<MultiPolygon>
           offset(const Geometry &g, const double &r, NoValidityCheck);

} // namespace algorithm
} // namespace SFCGAL

#endif

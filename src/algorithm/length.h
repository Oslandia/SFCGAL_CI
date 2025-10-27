// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_LENGTH_H_
#define SFCGAL_ALGORITHM_LENGTH_H_

#include "SFCGAL/config.h"

namespace SFCGAL {
class Geometry;
class LineString;
class GeometryCollection;
} // namespace SFCGAL

namespace SFCGAL::algorithm {

/**
 * @brief Compute the 2D length for a Geometry (0 for incompatible types)
 * @param g The geometry to compute length for
 * @return The 2D length of the geometry
 * @note When applied to NURBSCurve geometries, the 2D length
 *   is internally computed on a LineString obtained via
 *   toLineString() with its default parameters.
 */
SFCGAL_API auto
length(const Geometry &g) -> double;
/**
 * @brief Compute the 2D length for a LineString
 * @param g The linestring to compute length for
 * @return The 2D length of the linestring
 */
SFCGAL_API auto
length(const LineString &g) -> double;
/**
 * @brief Compute the 2D length for a GeometryCollection
 * @param g The geometry collection to compute length for
 * @return The total 2D length of all geometries in the collection
 */
SFCGAL_API auto
length(const GeometryCollection &g) -> double;

/**
 * @brief Compute the 3D length for a geometry
 * @param g The geometry to compute 3D length for
 * @return the 3D length of the Geometry, 0 for incompatible types
 */
SFCGAL_API auto
length3D(const Geometry &g) -> double;
/**
 * @brief Compute the 3D length for a LineString
 * @param g The linestring to compute 3D length for
 * @return The 3D length of the linestring
 */
SFCGAL_API auto
length3D(const LineString &g) -> double;
/**
 * @brief Compute the 3D length for a GeometryCollection
 * @param g The geometry collection to compute 3D length for
 * @return The total 3D length of all geometries in the collection
 */
SFCGAL_API auto
length3D(const GeometryCollection &g) -> double;

} // namespace SFCGAL::algorithm

#endif

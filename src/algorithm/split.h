// Copyright (c) 2025-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_SPLIT_H_
#define SFCGAL_ALGORITHM_SPLIT_H_

#include "SFCGAL/config.h"

#include <memory>

namespace SFCGAL {
class Geometry;
class Polygon;
class LineString;

namespace algorithm {
struct NoValidityCheck;

/**
 * @brief Split a polygon using a linestring.
 *
 * This function uses CGAL's Arrangement_2 with exact predicates and exact
 * constructions (EPECK kernel) to robustly split a polygon along a linestring.
 *
 * Algorithm:
 * 1. Insert polygon boundary edges into a planar arrangement
 * 2. Insert linestring segments into the arrangement
 * 3. Extract bounded faces from the arrangement
 * 4. Filter faces that lie within the original polygon
 * 5. Convert back to SFCGAL geometries with interpolated Z/M values
 *
 * Dimension Handling (Union Strategy):
 * The output dimension is the union of input dimensions with linear
 * interpolation of Z and M values from the nearest segment:
 *
 * | Polygon | LineString | Result | Interpolation Strategy |
 * |---------|------------|--------|------------------------|
 * | XY      | XY         | XY     | No interpolation needed |
 * | XYZ     | XY         | XYZ    | Z from polygon segments |
 * | XY      | XYZ        | XYZ    | Z from linestring for new points |
 * | XYZ     | XYZ        | XYZ    | Z from nearest segment |
 * | XYM     | XYM        | XYM    | M from nearest segment |
 * | XYZ     | XYM        | XYZM   | Z from polygon, M from linestring |
 * | XYZM    | XY         | XYZM   | Z and M from polygon |
 * | XYZM    | XYZM       | XYZM   | Z and M from nearest segment |
 *
 * For each new point on the split line, Z and M values are interpolated
 * linearly from the nearest segment in the original geometries. If a dimension
 * is not present in any input, it will be NaN in the output for those points.
 *
 * Edge Cases Handled:
 * - Linestring doesn't intersect polygon: returns original polygon
 * - Linestring only touches boundary (single point): returns original polygon
 * - Linestring collinear with polygon edge: returns original polygon
 * - Linestring passes through vertices: handled by arrangement
 * - Multiple crossings creating 3+ pieces: all pieces returned
 * - Polygon with holes: holes are preserved in appropriate pieces
 *
 * @param polygon The polygon to split (may have holes)
 * @param linestring The cutting linestring
 * @return GeometryCollection containing the split polygons, or original polygon
 * if split is not valid
 *
 * @pre polygon and linestring are valid geometries
 * @pre Arrangement operates on 2D projection (XY plane)
 *
 * @note Uses exact arithmetic (EPECK kernel) for robustness
 * @note Performance: O((n+k) log(n+k)) where n=polygon vertices, k=line
 * segments
 * @note Z and M interpolation uses nearest segment distance in XY plane
 *
 * Examples:
 * \code{.cpp}
 * // Simple vertical split (XY)
 * Polygon rect({{0,0}, {10,0}, {10,6}, {0,6}, {0,0}});
 * LineString line({{5,-1}, {5,7}});
 * auto result = split(rect, line);
 * // Returns GeometryCollection with 2 XY polygons
 *
 * // Split with Z interpolation (XYZ polygon + XY line)
 * Polygon rect3d({{0,0,10}, {10,0,20}, {10,6,30}, {0,6,15}, {0,0,10}});
 * LineString line2d({{5,-1}, {5,7}});
 * auto result2 = split(rect3d, line2d);
 * // Returns GeometryCollection with 2 XYZ polygons, Z interpolated on split
 *
 * // No split (line doesn't cross)
 * LineString external({{20,0}, {20,10}});
 * auto result3 = split(rect, external);
 * // Returns GeometryCollection with original polygon
 * \endcode
 */
SFCGAL_API auto
split(const Polygon &polygon, const LineString &linestring)
    -> std::unique_ptr<Geometry>;

/**
 * @brief Split a polygon using a linestring (no validity check variant).
 *
 * Same as split() but skips input geometry validity checks.
 * Use only when geometries are known to be valid.
 *
 * @param polygon The polygon to split
 * @param linestring The cutting linestring
 * @param nvc NoValidityCheck marker
 * @return GeometryCollection containing the split polygons
 *
 * @warning No input validity check is performed
 */
SFCGAL_API auto
split(const Polygon &polygon, const LineString &linestring, NoValidityCheck nvc)
    -> std::unique_ptr<Geometry>;

} // namespace algorithm
} // namespace SFCGAL

#endif // SFCGAL_ALGORITHM_SPLIT_H_

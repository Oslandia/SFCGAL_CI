// Copyright (c) 2012-2024, SFCGAL Contributors and Oslandia
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_MINKOWSKISUM3D_H_
#define SFCGAL_ALGORITHM_MINKOWSKISUM3D_H_

/**
 * @file minkowskiSum3D.h
 * @brief 3D Minkowski sum algorithm using CGAL Nef polyhedra.
 *
 * Provides functions to compute the Minkowski sum of two 3D geometries.
 * The implementation uses CGAL's Nef_polyhedron_3 which supports all
 * combinations of 0D (points), 1D (linestrings), 2D (surfaces), and
 * 3D (solids) geometry types.
 *
 * @see https://doc.cgal.org/latest/Minkowski_sum_3/index.html
 */

#include "SFCGAL/Geometry.h"
#include "SFCGAL/config.h"

#include <memory>

namespace SFCGAL::algorithm {

struct NoValidityCheck;

/**
 * @brief Computes the 3D Minkowski sum of two geometries (A + B).
 *
 * The Minkowski sum of two point sets A and B is defined as:
 *   A + B = { a + b : a in A, b in B }
 *
 * This function uses CGAL's Nef polyhedron implementation which supports
 * all combinations of the following geometry types:
 *
 * - **Point**: Singular vertex in Nef polyhedron
 * - **LineString**: Polyline using CGAL's native Polylines_tag constructor
 * - **NURBSCurve**: Converted to LineString first (using default parameters)
 * - **Triangle**: Singular convex facet
 * - **Polygon**: Triangulated to handle non-convex shapes and interior rings
 * - **TriangulatedSurface**: Converted directly to Nef polyhedron
 * - **PolyhedralSurface**: Converted directly to Nef polyhedron
 * - **Solid**: Exterior shell converted, interior shells (voids) subtracted
 * - **Multi* types**: Components combined via union
 * - **GeometryCollection**: Components combined via union
 *
 * @note The result is typically a PolyhedralSurface (closed shell) when both
 *       inputs have 3D extent, or may be empty when the result is degenerate
 *       (e.g., 2D+2D coplanar geometries).
 *
 * @note Complexity: O(n^3 * m^3) in worst case, where n and m are the
 *       complexities of the two input Nef polyhedra.
 *
 * @param gA The first geometry (validity checked)
 * @param gB The second geometry (validity checked)
 *
 * @return The 3D Minkowski sum as a unique_ptr<Geometry>:
 *         - PolyhedralSurface for volumetric results
 *         - GeometryCollection (empty) if the result is empty or degenerate
 *
 * @throws GeometryInvalidityException if input geometries are invalid
 * @throws GeometryInvalidityException for unsupported geometry types
 *
 * @pre gA and gB are valid 3D geometries
 * @post Result has validity flag set to true
 *
 * @see https://doc.cgal.org/latest/Minkowski_sum_3/index.html
 * @see minkowskiSum3D(const Geometry&, const Geometry&, NoValidityCheck) for
 *      unchecked version
 */
SFCGAL_API auto
minkowskiSum3D(const Geometry &gA, const Geometry &gB)
    -> std::unique_ptr<Geometry>;

/**
 * @brief Computes the 3D Minkowski sum without validity checking.
 *
 * Same as minkowskiSum3D(gA, gB) but skips input validation for performance.
 * Use when you know the inputs are already valid geometries.
 *
 * @param gA The first geometry (no validity check)
 * @param gB The second geometry (no validity check)
 *
 * @return The 3D Minkowski sum as a unique_ptr<Geometry>:
 *         - PolyhedralSurface for volumetric results
 *         - GeometryCollection (empty) if the result is empty or degenerate
 *
 * @warning Calling with invalid geometries may cause undefined behavior
 *          or incorrect results.
 *
 * @see minkowskiSum3D(const Geometry&, const Geometry&) for validated version
 */
SFCGAL_API auto
minkowskiSum3D(const Geometry &gA, const Geometry &gB, NoValidityCheck)
    -> std::unique_ptr<Geometry>;

} // namespace SFCGAL::algorithm

#endif // SFCGAL_ALGORITHM_MINKOWSKISUM3D_H_

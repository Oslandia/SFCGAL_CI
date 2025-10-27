// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_TRIANGULATE_TRIANGULATE2DZ_H_
#define SFCGAL_TRIANGULATE_TRIANGULATE2DZ_H_

#include "SFCGAL/config.h"

#include "SFCGAL/Geometry.h"
#include "SFCGAL/detail/triangulate/ConstraintDelaunayTriangulation.h"

namespace SFCGAL::triangulate {
/**
 * @brief Constraint 2DZ Delaunay Triangulation for geometry collections (keep Z
 * if defined, a projectionPlane may be provided)
 * @param g The input geometry collection to triangulate
 * @param triangulation The triangulation object to fill
 */
SFCGAL_API void
triangulateCollection2DZ(const Geometry                  &g,
                         ConstraintDelaunayTriangulation &triangulation);

/**
 * @brief Constraint 2DZ Delaunay Triangulation (keep Z if defined, a
 * projectionPlane may be provided)
 * @param g The input geometry to triangulate
 * @param triangulation The triangulation object to fill
 */
SFCGAL_API void
triangulate2DZ(const Geometry                  &g,
               ConstraintDelaunayTriangulation &triangulation);

/**
 * @brief Constraint 2DZ Delaunay Triangulation (keep Z if defined, project
 * points in OXY plane)
 * @param g The input geometry to triangulate
 * @return The constraint Delaunay triangulation
 */
SFCGAL_API ConstraintDelaunayTriangulation
triangulate2DZ(const Geometry &g);

} // namespace SFCGAL::triangulate

#endif

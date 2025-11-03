// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_ORIENTATION_H_
#define SFCGAL_ALGORITHM_ORIENTATION_H_

#include "SFCGAL/config.h"

#include "SFCGAL/Geometry.h"
#include "SFCGAL/Kernel.h"

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>

namespace SFCGAL::algorithm {

/**
 * Make valid 2D orientation
 * @param polygon the CGAL polygon to fix orientation for
 */
SFCGAL_API auto
makeValidOrientation(CGAL::Polygon_2<Kernel> &polygon) -> void;
/**
 * Make valid 2D orientation
 * @param polygon the CGAL polygon with holes to fix orientation for
 */
SFCGAL_API auto
makeValidOrientation(CGAL::Polygon_with_holes_2<Kernel> &polygon) -> void;
/**
 * Make valid 2D orientation
 * @param polygon the polygon to fix orientation for
 */
SFCGAL_API auto
makeValidOrientation(Polygon &polygon) -> void;

/**
 * Test if a Geometry has a consistent orientation
 * @param g the triangulated surface to test
 * @return true if the surface has consistent orientation, false otherwise
 */
SFCGAL_API auto
hasConsistentOrientation3D(const TriangulatedSurface &g) -> bool;
/**
 * Test if a PolyhedralSurface has a consistent orientation
 * @param g the polyhedral surface to test
 * @return true if the surface has consistent orientation, false otherwise
 */
SFCGAL_API auto
hasConsistentOrientation3D(const PolyhedralSurface &g) -> bool;

/**
 * Try to make consistent orientation in a TriangulatedSurface
 * @param g the triangulated surface to fix orientation for
 */
SFCGAL_API auto
makeConsistentOrientation3D(TriangulatedSurface &g) -> void;

/**
 * Test if a 2D surface is oriented counter clockwise
 * @param polygon the polygon to test
 * @return true if counter clockwise oriented, false otherwise
 */
SFCGAL_API auto
isCounterClockWiseOriented(const Polygon &polygon) -> bool;

/**
 * Test if a 2D surface is oriented counter clockwise
 * @param triangle the triangle to test
 * @return true if counter clockwise oriented, false otherwise
 */
SFCGAL_API auto
isCounterClockWiseOriented(const Triangle &triangle) -> bool;

/**
 * Test if a 2D surface is oriented counter clockwise
 * @param lineString the line string to test
 * @return true if counter clockwise oriented, false otherwise
 */
SFCGAL_API auto
isCounterClockWiseOriented(const LineString &lineString) -> bool;

} // namespace SFCGAL::algorithm

#endif

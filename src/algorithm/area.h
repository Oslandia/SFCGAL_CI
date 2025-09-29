// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_AREA_H_
#define SFCGAL_ALGORITHM_AREA_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/Kernel.h"

namespace SFCGAL {
namespace algorithm {
struct NoValidityCheck;

/**
 * @brief Compute the 2D area of a Geometry.
 *
 * @warning Z component is ignored, there is no 2D projection for 3D geometries
 * @pre The input geometry @p geom must be a valid geometry instance.
 *
 * @param geom The input geometry.
 *
 * @return The 2D area as a double precision value.
 */
SFCGAL_API auto
area(const Geometry &geom) -> double;

/**
 * @brief Compute the planar (2D) area of a Geometry without performing validity
 * checks.
 *
 * Dispatches on the runtime geometry type:
 * - Returns 0 for POINT, LINESTRING, NURBSCURVE, SOLID and MULTISOLID.
 * - Delegates to specialized overloads for POLYGON, TRIANGLE,
 * TRIANGULATEDSURFACE, POLYHEDRALSURFACE and collections (MULTI and
 * GEOMETRYCOLLECTION types).
 *
 * @warning Z component is ignored, there is no 2D projection for 3D geometries
 * @warning No actual validity check is done
 *
 * @pre The input geometry @p geom must be valid.
 *
 * @param geom The geometry whose planar area is computed.
 * @param noCheck Unused tag parameter indicating that geometry validity
 * is not checked.
 * @return double The computed planar area (>= 0).
 * @throws Exception If the geometry type is not handled by this dispatch.
 *
 * @ingroup detail
 */
SFCGAL_API auto
area(const Geometry &geom, NoValidityCheck noCheck) -> double;

/**
 * @brief Computes the 2D signed area of a Triangle.
 *
 * This function calculates the oriented area of a triangle in the plane.
 * The sign of the area depends on the orientation of the triangle's vertices:
 * - Positive if the points are ordered counterclockwise (CCW).
 * - Negative if the points are ordered clockwise (CW).
 * - Zero if the points are collinear.
 *
 * @param triangle The input Triangle whose signed area is to be computed.
 * @return Kernel::FT The signed area of the Triangle.
 *
 * @ingroup detail
 */
SFCGAL_API auto
signedArea(const Triangle &triangle) -> Kernel::FT;

/**
 * @brief Computes the 2D signed area of a LineString.
 *
 * This function calculates the oriented area of a lineString in the plane.
 * The sign of the area depends on the orientation of the lineString's vertices:
 * - Positive if the points are ordered counterclockwise (CCW).
 * - Negative if the points are ordered clockwise (CW).
 * - Zero if the points are collinear.
 *
 * @param lineString The input LineString whose signed area is to be computed.
 * @return Kernel::FT The signed area of the LineString.
 *
 * @ingroup detail
 */
SFCGAL_API auto
signedArea(const LineString &lineString) -> Kernel::FT;

/**
 * @brief Computes the 2D area of a Triangle.
 *
 * @param triangle The input Triangle.
 * @return The total area of the triangle.
 *
 * @ingroup detail
 */
SFCGAL_API auto
area(const Triangle &triangle) -> double;

/**
 * @brief Computes the 2D area of a Polygon.
 *
 * @param polygon The input Polygon.
 * @return The total area of the polygon.
 *
 * @ingroup detail
 */
SFCGAL_API auto
area(const Polygon &polygon) -> double;

/**
 * @brief Computes the 2D area of a GeometryCollection.
 *
 * @param collection The input GeometryCollection.
 * @return The total area of the collection.
 *
 * @ingroup detail
 */
SFCGAL_API auto
area(const GeometryCollection &collection) -> double;

/**
 * @brief Computes the 2D area of a TriangulatedSurface.
 *
 * @param tin The input TriangulatedSurface.
 * @return The total area of the TriangulatedSurface.
 *
 * @ingroup detail
 */
SFCGAL_API auto
area(const TriangulatedSurface &tin) -> double;

/**
 * @brief Computes the 2D area of a PolyhedralSurface.
 *
 * @param surface The input PolyhedralSurface.
 * @return The total area of the PolyhedralSurface.
 *
 * @ingroup detail
 */
SFCGAL_API auto
area(const PolyhedralSurface &surface) -> double;

/**
 * @brief Compute the 3D area of a Geometry.
 *
 * @warning Solid area is set to 0 (might be defined as the area of the surface)
 * @pre The input geometry @p geom must be a valid geometry instance.
 *
 * @param geom The input geometry.
 *
 * @return The 3D area as a double precision value.
 * @ingroup detail
 */
SFCGAL_API auto
area3D(const Geometry &geom) -> double;

/**
 * @brief Compute the 3D surface area of a geometry, dispatching by geometry
 * type.
 *
 * This overload assumes validity checks have already been performed (the
 * NoValidityCheck parameter is present to select this unchecked path) and
 * returns the 3D area in the geometry's own coordinate space.
 *
 * @warning Solid area is set to 0 (might be defined as the area of the surface)
 * @warning No actual validity check is done
 * @pre geom is a valid geometry
 *
 * Behavior by geometry type:
 * - POINT, LINESTRING, NURBSCURVE, SOLID, MULTISOLID: return 0.
 * - POLYGON, TRIANGLE, TRIANGULATEDSURFACE, POLYHEDRALSURFACE: delegate to
 *   the corresponding area3D(...) function for that type.
 * - MULTI types and GEOMETRYCOLLECTION: sum the 3D area of contained
 * geometries.
 *
 * @param geom Geometry to measure.
 * @param noCheck Unused tag parameter selecting the unchecked path.
 * @return double 3D surface area of the geometry.
 *
 * @ingroup detail
 */
SFCGAL_API auto
area3D(const Geometry &geom, NoValidityCheck noCheck) -> double;

/**
 * @brief Compute the 3D area of a Polygon.
 *
 * @pre The input geometry @p polygon must be a valid geometry instance.
 * @param polygon The input polygon.
 * @return The 3D area as a double precision value.
 * @ingroup detail
 */
SFCGAL_API auto
area3D(const Polygon &polygon) -> double;

/**
 * @brief Compute the 3D area of a GeometryCollection.
 *
 * @pre The input geometry @p collection must be a valid geometry instance.
 * @param collection The input geometryCollection.
 * @return The 3D area as a double precision value.
 * @ingroup detail
 */
SFCGAL_API auto
area3D(const GeometryCollection &collection) -> double;

/**
 * @brief Compute the 3D area of a PolyhedralSurface.
 *
 * @pre The input geometry @p surface must be a valid geometry instance.
 * @param surface The input PolyhedralSurface
 * @return The 3D area as a double precision value.
 * @ingroup detail
 */
SFCGAL_API auto
area3D(const PolyhedralSurface &surface) -> double;

/**
 * @brief Compute the 3D area of a Triangle.
 *
 * @pre The input geometry @p triangle must be a valid geometry instance.
 * @param triangle The input Triangle
 * @return The 3D area as a double precision value.
 * @ingroup detail
 */
SFCGAL_API auto
area3D(const Triangle &triangle) -> double;

/**
 * @brief Compute the 3D area of a TriangulatedSurface.
 *
 * @pre The input geometry @p tin must be a valid geometry instance.
 * @param tin The input TriangulatedSurface
 * @return The 3D area as a double precision value.
 * @ingroup detail
 */
SFCGAL_API auto
area3D(const TriangulatedSurface &tin) -> double;

} // namespace algorithm
} // namespace SFCGAL

#endif

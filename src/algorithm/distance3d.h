// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_DISTANCE3D_H_
#define SFCGAL_ALGORITHM_DISTANCE3D_H_

#include "SFCGAL/Geometry.h"
#include "SFCGAL/Kernel.h"
#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/config.h"

namespace SFCGAL {
namespace algorithm {
struct NoValidityCheck;

/**
 * @brief Dispatch 3D distance between two Geometries.
 *
 * @todo complete with solid
 *
 * @pre gA is a valid geometry
 * @pre gB is a valid geometry
 *
 * @param gA The first Geometry.
 * @param gB The second Geometry.
 * @return The 3D distance between the two Geometries.
 */
SFCGAL_API auto
distance3D(const Geometry &gA, const Geometry &gB) -> double;

/**
 * @brief Dispatch 3D distance between two Geometries.
 *
 *
 * @pre gA is a valid geometry
 * @pre gB is a valid geometry
 *
 * @warning No actual validity check is done
 *
 * @param gA The first Geometry.
 * @param gB The second Geometry.
 * @param noCheck validity check, unused parameter
 * @return The 3D distance between the two Geometries.
 */
SFCGAL_API auto
distance3D(const Geometry &gA, const Geometry &gB, NoValidityCheck noCheck)
    -> double;

/**
 * @brief Dispatch 3D distance from Point to Geometry
 *
 * @param gA The input Point.
 * @param gB The second Geometry.
 * @return The 3D distance between the two Geometries.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distancePointGeometry3D(const Point &gA, const Geometry &gB) -> double;

/**
 * @brief Computes the 3D distance between two Points.
 *
 * @param gA The first Point.
 * @param gB The second Point.
 * @return The 3D distance between the two Points.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distancePointPoint3D(const Point &gA, const Point &gB) -> double;

/**
 * @brief Computes the 3D distance between a Point and a LineString.
 *
 * @param gA The input Point.
 * @param gB The input LineString.
 * @return The 3D distance between the point and the LineString.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distancePointLineString3D(const Point &gA, const LineString &gB) -> double;

/**
 * @brief Computes the 3D distance between a Point and a Triangle.
 *
 * @param gA The input Point.
 * @param gB The input Triangle.
 * @return The 3D distance between the Point and the Triangle.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distancePointTriangle3D(const Point &gA, const Triangle &gB) -> double;

/**
 * @brief Computes the 3D distance between a Point and a Polygon.
 *
 * @param gA The input Point.
 * @param gB The input Polygon.
 * @return The 3D distance between the Point and the Polygon.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distancePointPolygon3D(const Point &gA, const Polygon &gB) -> double;

/**
 * @brief Computes the 3D distance between a Point and a PolyhedralSurface.
 *
 * @param pointA The input Point.
 * @param polySurfaceB The input PolyhedralSurface.
 * @return The 3D distance between the Point and the PolyhedralSurface.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distancePointPolyhedralSurface3D(const Point             &pointA,
                                 const PolyhedralSurface &polySurfaceB)
    -> double;

/**
 * @brief Computes the 3D distance between a Point and a TriangulatedSurface.
 *
 * @param pointA The input Point.
 * @param triangulatedSurfaceB The input TriangulatedSurface.
 * @return The 3D distance between the Point and the TriangulatedSurface.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distancePointTriangulatedSurface3D(
    const Point &pointA, const TriangulatedSurface &triangulatedSurfaceB)
    -> double;

/**
 * @brief Computes the 3D distance between a Point and a Solid.
 *
 * @param gA The input Point.
 * @param gB The input Solid.
 * @return The 3D distance between the Point and the Solid.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distancePointSolid3D(const Point &gA, const Solid &gB) -> double;

/**
 * @brief Dispatch 3D distance from LineString to Geometry
 *
 * @param gA The input LineString.
 * @param gB The input Geometry.
 * @return The 3D distance between the LineString and the Geometry
 *
 * @ingroup detail
 */
SFCGAL_API auto
distanceLineStringGeometry3D(const LineString &gA, const Geometry &gB)
    -> double;

/**
 * @brief Computes the 3D distance between two LineStrings.
 *
 * @param gA The first LineString.
 * @param gB The second LineString.
 * @return The 3D distance between the two LineStrings.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distanceLineStringLineString3D(const LineString &gA, const LineString &gB)
    -> double;

/**
 * @brief Computes the 3D distance between a LineString and a Triangle.
 *
 * @param gA The input LineString.
 * @param gB The input Triangle.
 * @return The 3D distance between the LineString and the Triangle.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distanceLineStringTriangle3D(const LineString &gA, const Triangle &gB)
    -> double;

/**
 * @brief Computes the 3D distance between a LineString and a Polygon.
 *
 * @param gA The input LineString.
 * @param gB The input Polygon.
 * @return The 3D distance between the LineString and the Polygon.
 *
 * @todo same method than distancePointPolygon3D (unify if triangulate is
 * available)
 *
 * @ingroup detail
 */
SFCGAL_API auto
distanceLineStringPolygon3D(const LineString &gA, const Polygon &gB) -> double;

/**
 * @brief Computes the 3D distance between a LineString and a PolyhedralSurface.
 *
 * @param lineA The input LineString.
 * @param polySurfaceB The input PolyhedralSurface.
 * @return The 3D distance between the LineString and the PolyhedralSurface.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distanceLineStringPolyhedralSurface3D(const LineString        &lineA,
                                      const PolyhedralSurface &polySurfaceB)
    -> double;

/**
 * @brief Computes the 3D distance between a LineString and a
 * TriangulatedSurface.
 *
 * @param lineA The input LineString.
 * @param triangulatedSurfaceB The input TriangulatedSurface.
 * @return The 3D distance between the LineString and the TriangulatedSurface.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distanceLineStringTriangulatedSurface3D(
    const LineString &lineA, const TriangulatedSurface &triangulatedSurfaceB)
    -> double;

/**
 * @brief Computes the 3D distance between a LineString and a Solid.
 *
 * @param gA The input LineString.
 * @param gB The input Solid.
 * @return The 3D distance between the LineString and the Solid.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distanceLineStringSolid3D(const LineString &gA, const Solid &gB) -> double;

/**
 * @brief Dispatch 3D distance from a Triangle to a Geometry
 *
 * @param gA The Triangle.
 * @param gB The input Geometry.
 * @return The 3D distance between the Triangle and the Geometry
 *
 * @ingroup detail
 */
SFCGAL_API auto
distanceTriangleGeometry3D(const Triangle &gA, const Geometry &gB) -> double;

/**
 * @brief Computes the 3D distance between a Triangle and a PolyhedralSurface.
 *
 * @param triangleA The input Triangle.
 * @param polySurfaceB The input PolyhedralSurface.
 * @return The 3D distance between the Triangle and the PolyhedralSurface.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distanceTrianglePolyhedralSurface3D(const Triangle          &triangleA,
                                    const PolyhedralSurface &polySurfaceB)
    -> double;

/**
 * @brief Computes the 3D distance between a Triangle and a Solid.
 *
 * @param gA The input Triangle.
 * @param gB The input Solid.
 * @return The 3D distance between the Triangle and the Solid.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distanceTriangleSolid3D(const Triangle &gA, const Solid &gB) -> double;

/**
 * @brief Dispatch 3D distance from a Polygon to a Geometry
 *
 * @param gA The Polygon.
 * @param gB The input Geometry.
 * @return The 3D distance between the Polygon and the Geometry
 *
 * @ingroup detail
 */
SFCGAL_API auto
distancePolygonGeometry3D(const Polygon &gA, const Geometry &gB) -> double;

/**
 * @brief Dispatch 3D distance from a PolyhedralSurface to a Geometry
 *
 * @param polySurfaceA The PolyhedralSurface.
 * @param geomB The input Geometry.
 * @return The 3D distance between the PolyhedralSurface and the Geometry
 *
 * @ingroup detail
 */
SFCGAL_API auto
distancePolyhedralSurfaceGeometry3D(const PolyhedralSurface &polySurfaceA,
                                    const Geometry          &geomB) -> double;

/**
 * @brief Dispatch 3D distance from a TriangulatedSurface to a Geometry
 *
 * @param triangulatedSurfaceA The TriangulatedSurface.
 * @param geomB The input Geometry.
 * @return The 3D distance between the TriangulatedSurface and the Geometry
 *
 * @ingroup detail
 */
SFCGAL_API auto
distanceTriangulatedSurfaceGeometry3D(
    const TriangulatedSurface &triangulatedSurfaceA, const Geometry &geomB)
    -> double;

/**
 * @brief Dispatch 3D distance from a Solid to a Geometry
 *
 * @param gA The Solid.
 * @param gB The input Geometry.
 * @return The 3D distance between the Solid and the Geometry
 *
 * @ingroup detail
 */
SFCGAL_API auto
distanceSolidGeometry3D(const Solid &gA, const Geometry &gB) -> double;

/**
 * @brief Computes the 3D distance between two Solids.
 *
 * @param gA The first Solid.
 * @param gB The second Solid.
 * @return The 3D distance between the two Solids.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distanceSolidSolid3D(const Solid &gA, const Solid &gB) -> double;

/**
 * @brief Dispatch 3D distance from a GeometryCollection to a Geometry
 *
 * @param gA The GeometryCollection.
 * @param gB The input Geometry.
 * @return The 3D distance between the GeometryCollection and the Geometry
 *
 * @ingroup detail
 */
SFCGAL_API auto
distanceGeometryCollectionToGeometry3D(const Geometry &gA, const Geometry &gB)
    -> double;

/**
 * @brief Computes the 3D distance between a Point and a Segment.
 *
 * @param p The input Point.
 * @param a The start Point of the segment.
 * @param b The end Point of the segment.
 * @return The 3D distance between the Point and the Segment.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distancePointSegment3D(const Point &p, const Point &a, const Point &b)
    -> double;

/**
 * @brief Computes the 3D distance between a Point and a Triangle.
 *
 * @param p The input Point.
 * @param a The first Vertex of the Triangle
 * @param b The second Vertex of the Triangle
 * @param c The third Vertex of the Triangle
 * @return The 3D distance between the Point and the Segment.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distancePointTriangle3D(const Point &p, const Point &a, const Point &b,
                        const Point &c) -> double;

/**
 * @brief Computes the 3D distance between two Segments.
 *
 * @param a The start Point of the first segment.
 * @param b The end Point of the first segment.
 * @param c The start Point of the second segment.
 * @param d The end Point of the second segment.
 * @return The 3D distance between the two Segments.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distanceSegmentSegment3D(const Point &a, const Point &b, const Point &c,
                         const Point &d) -> double;

/**
 * @brief Computes the 3D distance between a Segment and a Triangle.
 *
 * @param sA The start Point of the Segment
 * @param sB The end Point of the Segment
 * @param tA The first Vertex of the Triangle
 * @param tB The second Vertex of the Triangle
 * @param tC The third Vertex of the Triangle
 * @return The 3D distance between the Segment and the Triangle.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distanceSegmentTriangle3D(const Point &sA, const Point &sB, const Point &tA,
                          const Point &tB, const Point &tC) -> double;

/**
 * @brief Computes the 3D distance between two Triangles.
 *
 * @param gA The first Triangle
 * @param gB The second Triangle
 * @return The 3D distance between the Segment and the Triangle.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distanceTriangleTriangle3D(const Triangle &gA, const Triangle &gB) -> double;

} // namespace algorithm
} // namespace SFCGAL

#endif

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

namespace SFCGAL::algorithm {
struct NoValidityCheck;

/**
 * @brief Dispatch 3D distance between two Geometries.
 *
 * @param geometry1 First geometry
 * @param geometry2 Second geometry
 * @return 3D distance between the geometries
 * @note When applied to NURBSCurve geometries, the distance
 *   is internally computed on a LineString obtained via
 *   toLineString() with its default parameters.
 * @todo complete with solid
 */
SFCGAL_API auto
distance3D(const Geometry &geometry1, const Geometry &geometry2) -> double;

/**
 * @brief Dispatch 3D distance between two Geometries.
 *
 * @pre geometry1 is a valid geometry
 * @pre geometry2 is a valid geometry
 *
 * @param geometry1 First geometry
 * @param geometry2 Second geometry
 * @param noCheck Validity check parameter
 * @return 3D distance between the geometries
 * @warning No actual validity check is done
 *
 */
SFCGAL_API auto
distance3D(const Geometry &geometry1, const Geometry &geometry2,
           NoValidityCheck noCheck) -> double;

/**
 * @brief Dispatch 3D distance from Point to Geometry
 *
 * @param point The input Point.
 * @param geometry The second Geometry.
 * @return The 3D distance between the two Geometries.
 *
 */
SFCGAL_API auto
distancePointGeometry3D(const Point &point, const Geometry &geometry) -> double;

/**
 * @brief Computes the 3D distance between two Points.
 *
 * @param point1 The first Point.
 * @param point2 The second Point.
 * @return The 3D distance between the two Points.
 *
 */
SFCGAL_API auto
distancePointPoint3D(const Point &point1, const Point &point2) -> double;

/**
 * @brief Computes the 3D distance between a Point and a LineString.
 *
 * @param point The input Point.
 * @param lineString The input LineString.
 * @return The 3D distance between the point and the LineString.
 *
 */
SFCGAL_API auto
distancePointLineString3D(const Point &point, const LineString &lineString)
    -> double;

/**
 * @brief Computes the 3D distance between a Point and a Triangle.
 *
 * @param point The input Point.
 * @param triangle The input Triangle.
 * @return The 3D distance between the Point and the Triangle.
 *
 */
SFCGAL_API auto
distancePointTriangle3D(const Point &point, const Triangle &triangle) -> double;

/**
 * @brief Computes the 3D distance between a Point and a Polygon.
 *
 * @param point The input Point.
 * @param polygon The input Polygon.
 * @return The 3D distance between the Point and the Polygon.
 *
 */
SFCGAL_API auto
distancePointPolygon3D(const Point &point, const Polygon &polygon) -> double;

/**
 * @brief Computes the 3D distance between a Point and a PolyhedralSurface.
 *
 * @param pointA The input Point.
 * @param polySurfaceB The input PolyhedralSurface.
 * @return The 3D distance between the Point and the PolyhedralSurface.
 *
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
 */
SFCGAL_API auto
distancePointTriangulatedSurface3D(
    const Point &pointA, const TriangulatedSurface &triangulatedSurfaceB)
    -> double;

/**
 * @brief Computes the 3D distance between a Point and a Solid.
 *
 * @param point The input Point.
 * @param solid The input Solid.
 * @return The 3D distance between the Point and the Solid.
 *
 */
SFCGAL_API auto
distancePointSolid3D(const Point &point, const Solid &solid) -> double;

/**
 * @brief Dispatch 3D distance from LineString to Geometry
 *
 * @param lineString The input LineString.
 * @param geometry The input Geometry.
 * @return The 3D distance between the LineString and the Geometry
 *
 */
SFCGAL_API auto
distanceLineStringGeometry3D(const LineString &lineString,
                             const Geometry   &geometry) -> double;

/**
 * @brief Computes the 3D distance between two LineStrings.
 *
 * @param lineString1 The first LineString.
 * @param lineString2 The second LineString.
 * @return The 3D distance between the two LineStrings.
 *
 */
SFCGAL_API auto
distanceLineStringLineString3D(const LineString &lineString1,
                               const LineString &lineString2) -> double;

/**
 * @brief Computes the 3D distance between a LineString and a Triangle.
 *
 * @param lineString The input LineString.
 * @param triangle The input Triangle.
 * @return The 3D distance between the LineString and the Triangle.
 *
 */
SFCGAL_API auto
distanceLineStringTriangle3D(const LineString &lineString,
                             const Triangle   &triangle) -> double;

/**
 * @brief Computes the 3D distance between a LineString and a Polygon.
 *
 * @param lineString The input LineString.
 * @param polygon The input Polygon.
 * @return The 3D distance between the LineString and the Polygon.
 *
 * @todo same method than distancePointPolygon3D (unify if triangulate is
 * available)
 *
 */
SFCGAL_API auto
distanceLineStringPolygon3D(const LineString &lineString,
                            const Polygon    &polygon) -> double;

/**
 * @brief Computes the 3D distance between a LineString and a PolyhedralSurface.
 *
 * @param lineA The input LineString.
 * @param polySurfaceB The input PolyhedralSurface.
 * @return The 3D distance between the LineString and the PolyhedralSurface.
 *
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
 */
SFCGAL_API auto
distanceLineStringTriangulatedSurface3D(
    const LineString &lineA, const TriangulatedSurface &triangulatedSurfaceB)
    -> double;

/**
 * @brief Computes the 3D distance between a LineString and a Solid.
 *
 * @param lineString The input LineString.
 * @param solid The input Solid.
 * @return The 3D distance between the LineString and the Solid.
 *
 */
SFCGAL_API auto
distanceLineStringSolid3D(const LineString &lineString, const Solid &solid)
    -> double;

/**
 * @brief Dispatch 3D distance from a Triangle to a Geometry
 *
 * @param triangle The Triangle.
 * @param geometry The input Geometry.
 * @return The 3D distance between the Triangle and the Geometry
 *
 */
SFCGAL_API auto
distanceTriangleGeometry3D(const Triangle &triangle, const Geometry &geometry)
    -> double;

/**
 * @brief Computes the 3D distance between a Triangle and a PolyhedralSurface.
 *
 * @param triangleA The input Triangle.
 * @param polySurfaceB The input PolyhedralSurface.
 * @return The 3D distance between the Triangle and the PolyhedralSurface.
 *
 */
SFCGAL_API auto
distanceTrianglePolyhedralSurface3D(const Triangle          &triangleA,
                                    const PolyhedralSurface &polySurfaceB)
    -> double;

/**
 * @brief Computes the 3D distance between a Triangle and a Solid.
 *
 * @param triangle The input Triangle.
 * @param solid The input Solid.
 * @return The 3D distance between the Triangle and the Solid.
 *
 */
SFCGAL_API auto
distanceTriangleSolid3D(const Triangle &triangle, const Solid &solid) -> double;

/**
 * @brief Dispatch 3D distance from a Polygon to a Geometry
 *
 * @param polygon The Polygon.
 * @param geometry The input Geometry.
 * @return The 3D distance between the Polygon and the Geometry
 *
 */
SFCGAL_API auto
distancePolygonGeometry3D(const Polygon &polygon, const Geometry &geometry)
    -> double;

/**
 * @brief Dispatch 3D distance from a PolyhedralSurface to a Geometry
 *
 * @param polySurfaceA The PolyhedralSurface.
 * @param geomB The input Geometry.
 * @return The 3D distance between the PolyhedralSurface and the Geometry
 *
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
 */
SFCGAL_API auto
distanceTriangulatedSurfaceGeometry3D(
    const TriangulatedSurface &triangulatedSurfaceA, const Geometry &geomB)
    -> double;

/**
 * @brief Dispatch 3D distance from a Solid to a Geometry
 *
 * @param solid The Solid.
 * @param geometry The input Geometry.
 * @return The 3D distance between the Solid and the Geometry
 *
 */
SFCGAL_API auto
distanceSolidGeometry3D(const Solid &solid, const Geometry &geometry) -> double;

/**
 * @brief Computes the 3D distance between two Solids.
 *
 * @param solid1 The first Solid.
 * @param solid2 The second Solid.
 * @return The 3D distance between the two Solids.
 *
 */
SFCGAL_API auto
distanceSolidSolid3D(const Solid &solid1, const Solid &solid2) -> double;

/**
 * @brief Dispatch 3D distance from a GeometryCollection to a Geometry
 *
 * @param geometryCollection The GeometryCollection.
 * @param geometry The input Geometry.
 * @return The 3D distance between the GeometryCollection and the Geometry
 *
 */
SFCGAL_API auto
distanceGeometryCollectionToGeometry3D(const Geometry &geometryCollection,
                                       const Geometry &geometry) -> double;

/**
 * @brief Computes the 3D distance between a Point and a Segment.
 *
 * @param point The input Point.
 * @param startPoint The start Point of the segment.
 * @param endPoint The end Point of the segment.
 * @return The 3D distance between the Point and the Segment.
 *
 */
SFCGAL_API auto
distancePointSegment3D(const Point &point, const Point &startPoint,
                       const Point &endPoint) -> double;

/**
 * @brief Computes the 3D distance between a Point and a Triangle.
 *
 * @param point The input Point.
 * @param vertex1 The first Vertex of the Triangle
 * @param vertex2 The second Vertex of the Triangle
 * @param vertex3 The third Vertex of the Triangle
 * @return The 3D distance between the Point and the Segment.
 *
 */
SFCGAL_API auto
distancePointTriangle3D(const Point &point, const Point &vertex1,
                        const Point &vertex2, const Point &vertex3) -> double;

/**
 * @brief Computes the 3D distance between two Segments.
 *
 * @param startPoint1 The start Point of the first segment.
 * @param endPoint1 The end Point of the first segment.
 * @param startPoint2 The start Point of the second segment.
 * @param endPoint2 The end Point of the second segment.
 * @return The 3D distance between the two Segments.
 *
 */
SFCGAL_API auto
distanceSegmentSegment3D(const Point &startPoint1, const Point &endPoint1,
                         const Point &startPoint2, const Point &endPoint2)
    -> double;

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
 */
SFCGAL_API auto
distanceSegmentTriangle3D(const Point &sA, const Point &sB, const Point &tA,
                          const Point &tB, const Point &tC) -> double;

/**
 * @brief Computes the 3D distance between two Triangles.
 *
 * @param triangle1 The first Triangle
 * @param triangle2 The second Triangle
 * @return The 3D distance between the Segment and the Triangle.
 *
 */
SFCGAL_API auto
distanceTriangleTriangle3D(const Triangle &triangle1, const Triangle &triangle2)
    -> double;

} // namespace SFCGAL::algorithm

#endif

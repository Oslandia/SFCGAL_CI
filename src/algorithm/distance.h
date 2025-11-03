// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_DISTANCE_H_
#define SFCGAL_ALGORITHM_DISTANCE_H_

#include "SFCGAL/config.h"

#include "SFCGAL/Geometry.h"

namespace SFCGAL::algorithm {
struct NoValidityCheck;

/**
 * @brief Compute the distance between two Geometries.
 *
 * @param geometry1 First geometry
 * @param geometry2 Second geometry
 * @return Distance between the geometries
 * @note When applied to NURBSCurve geometries, the distance
 *   is internally computed on a LineString obtained via
 *   toLineString() with its default parameters.
 * @warning No actual validity check is done
 *
 */
SFCGAL_API auto
distance(const Geometry &geometry1, const Geometry &geometry2) -> double;

/**
 * @brief Dispatch distance between two Geometries.
 *
 * @pre geometry1 is a valid geometry
 * @pre geometry2 is a valid geometry
 *
 * @param geometry1 First geometry
 * @param geometry2 Second geometry
 * @param noCheck Validity check parameter
 * @return Distance between the geometries
 * @warning No actual validity check is done
 *
 */
SFCGAL_API auto
distance(const Geometry &geometry1, const Geometry &geometry2,
         NoValidityCheck noCheck) -> double;

/**
 * @brief Dispatch distance from Point to Geometry
 *
 * @param point The input Point.
 * @param geometry The second Geometry.
 * @return The distance between the two Geometries.
 *
 */
SFCGAL_API auto
distancePointGeometry(const Point &point, const Geometry &geometry) -> double;

/**
 * @brief Computes the distance between two Points.
 *
 * @param point1 The first Point.
 * @param point2 The second Point.
 * @return The distance between the two Points.
 *
 */
SFCGAL_API auto
distancePointPoint(const Point &point1, const Point &point2) -> double;

/**
 * @brief Computes the distance between a Point and a LineString.
 *
 * @param point The input Point.
 * @param lineString The input LineString.
 * @return The distance between the point and the LineString.
 *
 */
SFCGAL_API auto
distancePointLineString(const Point &point, const LineString &lineString)
    -> double;

/**
 * @brief Computes the distance between a Point and a Polygon.
 *
 * @param point The input Point.
 * @param polygon The input Polygon.
 * @return The distance between the Point and the Polygon.
 *
 */
SFCGAL_API auto
distancePointPolygon(const Point &point, const Polygon &polygon) -> double;

/**
 * @brief Computes the distance between a Point and a Triangle.
 *
 * @param point The input Point.
 * @param triangle The input Triangle.
 * @return The distance between the Point and the Triangle.
 *
 */
SFCGAL_API auto
distancePointTriangle(const Point &point, const Triangle &triangle) -> double;

/**
 * @brief Dispatch distance from LineString to Geometry
 *
 * @param lineString The input LineString.
 * @param geometry The input Geometry.
 * @return The distance between the LineString and the Geometry
 *
 */
SFCGAL_API auto
distanceLineStringGeometry(const LineString &lineString,
                           const Geometry   &geometry) -> double;

/**
 * @brief Computes the distance between two LineStrings.
 *
 * @param lineString1 The first LineString.
 * @param lineString2 The second LineString.
 * @return The distance between the two LineStrings.
 *
 */
SFCGAL_API auto
distanceLineStringLineString(const LineString &lineString1,
                             const LineString &lineString2) -> double;

/**
 * @brief Computes the distance between a LineString and a Polygon.
 *
 * @param lineString The input LineString.
 * @param polygon The input Polygon.
 * @return The distance between the LineString and the Polygon.
 *
 */
SFCGAL_API auto
distanceLineStringPolygon(const LineString &lineString, const Polygon &polygon)
    -> double;

/**
 * @brief Computes the distance between a LineString and a Triangle.
 *
 * @param lineString The input LineString.
 * @param triangle The input Triangle.
 * @return The distance between the LineString and the Triangle.
 *
 */
SFCGAL_API auto
distanceLineStringTriangle(const LineString &lineString,
                           const Triangle   &triangle) -> double;

/**
 * @brief Dispatch distance from Polygon to Geometry
 *
 * @param polygon The Polygon.
 * @param geometry The input Geometry.
 * @return The distance between the Polygon and the Geometry
 *
 */
SFCGAL_API auto
distancePolygonGeometry(const Polygon &polygon, const Geometry &geometry)
    -> double;

/**
 * @brief Computes the distance between two Polygons.
 *
 * @param polygon1 The first Polygon.
 * @param polygon2 The second Polygon.
 * @return The distance between the two Polygons.
 *
 */
SFCGAL_API auto
distancePolygonPolygon(const Polygon &polygon1, const Polygon &polygon2)
    -> double;

/**
 * @brief Computes the distance between a Polygon and a Triangle.
 *
 * @param polygon The input Polygon.
 * @param triangle The input Triangle.
 * @return The distance between the Polygon and the Triangle.
 *
 */
SFCGAL_API auto
distancePolygonTriangle(const Polygon &polygon, const Triangle &triangle)
    -> double;

/**
 * @brief Dispatch distance from a Triangle to a Geometry
 *
 * @param triangle The Triangle.
 * @param geometry The input Geometry.
 * @return The distance between the Triangle and the Geometry
 *
 */
SFCGAL_API auto
distanceTriangleGeometry(const Triangle &triangle, const Geometry &geometry)
    -> double;

/**
 * @brief Computes the distance from a GeometryCollection to a Geometry
 *
 * @param geometryCollection The GeometryCollection.
 * @param geometry The input Geometry.
 * @return The distance between the GeometryCollection and the Geometry
 *
 */
SFCGAL_API auto
distanceGeometryCollectionToGeometry(const Geometry &geometryCollection,
                                     const Geometry &geometry) -> double;

/**
 * @brief Computes the distance between a Point and a Segment.
 *
 * @param point The input Point.
 * @param segmentStart The start Point of the segment.
 * @param segmentEnd The end Point of the segment.
 * @return The distance between the Point and the Segment.
 *
 */
SFCGAL_API auto
distancePointSegment(const Point &point, const Point &segmentStart,
                     const Point &segmentEnd) -> double;

/**
 * @brief Computes the distance between two Segments.
 *
 * @param point1 The start Point of the first segment.
 * @param point2 The end Point of the first segment.
 * @param point3 The start Point of the second segment.
 * @param point4 The end Point of the second segment.
 * @return The distance between the two Segments.
 *
 */
SFCGAL_API auto
distanceSegmentSegment(const Point &point1, const Point &point2,
                       const Point &point3, const Point &point4) -> double;

} // namespace SFCGAL::algorithm

#endif

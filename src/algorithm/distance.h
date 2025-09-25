// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_DISTANCE_H_
#define SFCGAL_ALGORITHM_DISTANCE_H_

#include "SFCGAL/config.h"

#include "SFCGAL/Geometry.h"

namespace SFCGAL {
namespace algorithm {
struct NoValidityCheck;

/**
 * @brief Compute the distance between two Geometries.
 *
 * @param gA First geometry
 * @param gB Second geometry
 * @return Distance between the geometries
 * @warning No actual validity check is done
 *
 * @ingroup detail
 */
SFCGAL_API auto
distance(const Geometry &gA, const Geometry &gB) -> double;

/**
 * @brief Dispatch distance between two Geometries.
 *
 * @pre gA is a valid geometry
 * @pre gB is a valid geometry
 *
 * @param gA First geometry
 * @param gB Second geometry
 * @param noCheck Validity check parameter
 * @return Distance between the geometries
 * @warning No actual validity check is done
 *
 * @ingroup detail
 */
SFCGAL_API auto
distance(const Geometry &gA, const Geometry &gB, NoValidityCheck noCheck)
    -> double;

/**
 * @brief Dispatch distance from Point to Geometry
 *
 * @param gA The input Point.
 * @param gB The second Geometry.
 * @return The distance between the two Geometries.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distancePointGeometry(const Point &gA, const Geometry &gB) -> double;

/**
 * @brief Computes the distance between two Points.
 *
 * @param gA The first Point.
 * @param gB The second Point.
 * @return The distance between the two Points.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distancePointPoint(const Point &gA, const Point &gB) -> double;

/**
 * @brief Computes the distance between a Point and a LineString.
 *
 * @param gA The input Point.
 * @param gB The input LineString.
 * @return The distance between the point and the LineString.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distancePointLineString(const Point &gA, const LineString &gB) -> double;

/**
 * @brief Computes the distance between a Point and a Polygon.
 *
 * @param gA The input Point.
 * @param gB The input Polygon.
 * @return The distance between the Point and the Polygon.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distancePointPolygon(const Point &gA, const Polygon &gB) -> double;

/**
 * @brief Computes the distance between a Point and a Triangle.
 *
 * @param gA The input Point.
 * @param gB The input Triangle.
 * @return The distance between the Point and the Triangle.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distancePointTriangle(const Point &gA, const Triangle &gB) -> double;

/**
 * @brief Dispatch distance from LineString to Geometry
 *
 * @param gA The input LineString.
 * @param gB The input Geometry.
 * @return The distance between the LineString and the Geometry
 *
 * @ingroup detail
 */
SFCGAL_API auto
distanceLineStringGeometry(const LineString &gA, const Geometry &gB) -> double;

/**
 * @brief Computes the distance between two LineStrings.
 *
 * @param gA The first LineString.
 * @param gB The second LineString.
 * @return The distance between the two LineStrings.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distanceLineStringLineString(const LineString &gA, const LineString &gB)
    -> double;

/**
 * @brief Computes the distance between a LineString and a Polygon.
 *
 * @param gA The input LineString.
 * @param gB The input Polygon.
 * @return The distance between the LineString and the Polygon.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distanceLineStringPolygon(const LineString &gA, const Polygon &gB) -> double;

/**
 * @brief Computes the distance between a LineString and a Triangle.
 *
 * @param gA The input LineString.
 * @param gB The input Triangle.
 * @return The distance between the LineString and the Triangle.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distanceLineStringTriangle(const LineString &gA, const Triangle &gB) -> double;

/**
 * @brief Dispatch distance from Polygon to Geometry
 *
 * @param gA The Polygon.
 * @param gB The input Geometry.
 * @return The distance between the Polygon and the Geometry
 *
 * @ingroup detail
 */
SFCGAL_API auto
distancePolygonGeometry(const Polygon &gA, const Geometry &gB) -> double;

/**
 * @brief Computes the distance between two Polygons.
 *
 * @param gA The first Polygon.
 * @param gB The second Polygon.
 * @return The distance between the two Polygons.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distancePolygonPolygon(const Polygon &gA, const Polygon &gB) -> double;

/**
 * @brief Computes the distance between a Polygon and a Triangle.
 *
 * @param gA The input Polygon.
 * @param gB The input Triangle.
 * @return The distance between the Polygon and the Triangle.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distancePolygonTriangle(const Polygon &gA, const Triangle &gB) -> double;

/**
 * @brief Dispatch distance from a Triangle to a Geometry
 *
 * @param gA The Triangle.
 * @param gB The input Geometry.
 * @return The distance between the Triangle and the Geometry
 *
 * @ingroup detail
 */
SFCGAL_API auto
distanceTriangleGeometry(const Triangle &gA, const Geometry &gB) -> double;

/**
 * @brief Computes the distance from a GeometryCollection to a Geometry
 *
 * @param gA The GeometryCollection.
 * @param gB The input Geometry.
 * @return The distance between the GeometryCollection and the Geometry
 *
 * @ingroup detail
 */
SFCGAL_API auto
distanceGeometryCollectionToGeometry(const Geometry &gA, const Geometry &gB)
    -> double;

/**
 * @brief Computes the distance between a Point and a Segment.
 *
 * @param p The input Point.
 * @param a The start Point of the segment.
 * @param b The end Point of the segment.
 * @return The distance between the Point and the Segment.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distancePointSegment(const Point &p, const Point &a, const Point &b) -> double;

/**
 * @brief Computes the distance between two Segments.
 *
 * @param a The start Point of the first segment.
 * @param b The end Point of the first segment.
 * @param c The start Point of the second segment.
 * @param d The end Point of the second segment.
 * @return The distance between the two Segments.
 *
 * @ingroup detail
 */
SFCGAL_API auto
distanceSegmentSegment(const Point &a, const Point &b, const Point &c,
                       const Point &d) -> double;

} // namespace algorithm
} // namespace SFCGAL

#endif

// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _SFCGAL_ALGORITHM_DISTANCE_H_
#define _SFCGAL_ALGORITHM_DISTANCE_H_

#include <SFCGAL/config.h>

#include <SFCGAL/Geometry.h>

namespace SFCGAL {
namespace algorithm {
struct NoValidityCheck;

/**
 * Compute the distance between two Geometries
 * @ingroup public_api
 * @pre gA is a valid geometry
 * @pre gB is a valid geometry
 */
SFCGAL_API double
distance(const Geometry &gA, const Geometry &gB);

/**
 * Compute the distance between two Geometries
 * @ingroup detail
 * @pre gA is a valid geometry
 * @pre gB is a valid geometry
 * @warning No actual validity check is done
 */
SFCGAL_API double
distance(const Geometry &gA, const Geometry &gB, NoValidityCheck);

/**
 * dispatch distance from Point to Geometry
 * @ingroup detail
 */
SFCGAL_API double
distancePointGeometry(const Point &gA, const Geometry &gB);
/**
 * distance between two Points
 * @ingroup detail
 */
SFCGAL_API double
distancePointPoint(const Point &gA, const Point &gB);
/**
 * distance between a Point and a LineString
 * @ingroup detail
 */
SFCGAL_API double
distancePointLineString(const Point &gA, const LineString &gB);
/**
 * distance between a Point and a Polygon
 * @ingroup detail
 */
SFCGAL_API double
distancePointPolygon(const Point &gA, const Polygon &gB);
/**
 * distance between a Point and a Triangle
 * @ingroup detail
 */
SFCGAL_API double
distancePointTriangle(const Point &gA, const Triangle &gB);

/**
 * dispatch distance from LineString to Geometry
 * @ingroup detail
 */
SFCGAL_API double
distanceLineStringGeometry(const LineString &gA, const Geometry &gB);
/**
 * distance between two LineStrings
 * @ingroup detail
 */
SFCGAL_API double
distanceLineStringLineString(const LineString &gA, const LineString &gB);
/**
 * distance between a LineString and a Polygon
 * @ingroup detail
 */
SFCGAL_API double
distanceLineStringPolygon(const LineString &gA, const Polygon &gB);
/**
 * distance between a LineString and a Triangle
 * @ingroup detail
 */
SFCGAL_API double
distanceLineStringTriangle(const LineString &gA, const Triangle &gB);

/**
 * dispatch distance from Polygon to Geometry
 * @ingroup detail
 */
SFCGAL_API double
distancePolygonGeometry(const Polygon &gA, const Geometry &gB);
/**
 * distance between two Polygons
 * @ingroup detail
 */
SFCGAL_API double
distancePolygonPolygon(const Polygon &gA, const Polygon &gB);
/**
 * distance between a Polygon and a Triangle
 * @ingroup detail
 */
SFCGAL_API double
distancePolygonTriangle(const Polygon &gA, const Triangle &gB);

/**
 * dispatch distance from a Triangle to a Geometry
 * @ingroup detail
 */
SFCGAL_API double
distanceTriangleGeometry(const Triangle &gA, const Geometry &gB);

/**
 * dispatch distance from a collection of geometry (gA) to a Geometry (gB)
 * @ingroup detail
 */
SFCGAL_API double
distanceGeometryCollectionToGeometry(const Geometry &gA, const Geometry &gB);

/**
 * @ingroup detail
 */
SFCGAL_API double
distancePointSegment(const Point &p, const Point &a, const Point &b);

/**
 * @ingroup detail
 */
SFCGAL_API double
distanceSegmentSegment(const Point &a, const Point &b, const Point &c,
                       const Point &d);

} // namespace algorithm
} // namespace SFCGAL

#endif

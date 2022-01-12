// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: GPL-2.0-or-later

#ifndef _SFCGAL_ALGORITHM_DISTANCE3D_H_
#define _SFCGAL_ALGORITHM_DISTANCE3D_H_

#include <SFCGAL/config.h>

#include <SFCGAL/Geometry.h>

namespace SFCGAL {
namespace algorithm {
struct NoValidityCheck;

/**
 * dispatch distance between two Geometries
 * @todo complete with solid
 * Compute distance between two 3D Geometries
 * @ingroup public_api
 * @pre gA is a valid geometry
 * @pre gB is a valid geometry
 */
SFCGAL_API double
distance3D(const Geometry &gA, const Geometry &gB);

/**
 * Compute distance between two 3D Geometries
 * @ingroup detail
 * @pre gA is a valid geometry
 * @pre gB is a valid geometry
 * @warning No actual validity check is done
 */
SFCGAL_API double
distance3D(const Geometry &gA, const Geometry &gB, NoValidityCheck);

/**
 * dispatch distance from Point to Geometry
 * @ingroup detail
 */
SFCGAL_API double
distancePointGeometry3D(const Point &gA, const Geometry &gB);
/**
 * distance between two Points
 * @ingroup detail
 */
SFCGAL_API double
distancePointPoint3D(const Point &gA, const Point &gB);
/**
 * distance between a Point and a LineString
 * @ingroup detail
 */
SFCGAL_API double
distancePointLineString3D(const Point &gA, const LineString &gB);
/**
 * distance between a Point and a Triangle
 * @ingroup detail
 */
SFCGAL_API double
distancePointTriangle3D(const Point &gA, const Triangle &gB);
/**
 * distance between a Point and a Triangle
 * @ingroup detail
 */
SFCGAL_API double
distancePointPolygon3D(const Point &gA, const Polygon &gB);
/**
 * distance between a Point and a Solid
 * @ingroup detail
 */
SFCGAL_API double
distancePointSolid3D(const Point &gA, const Solid &gB);

/**
 * dispatch distance between a LineString and a Geometry
 * @ingroup detail
 */
SFCGAL_API double
distanceLineStringGeometry3D(const LineString &gA, const Geometry &gB);
/**
 * distance between two LineStrings
 * @ingroup detail
 */
SFCGAL_API double
distanceLineStringLineString3D(const LineString &gA, const LineString &gB);
/**
 * distance between a LineString and a Triangle
 * @ingroup detail
 */
SFCGAL_API double
distanceLineStringTriangle3D(const LineString &gA, const Triangle &gB);
/**
 * distance between a LineString and a Polygon
 * @todo same method than distancePointPolygon3D (unify if triangulate is
 * available)
 * @ingroup detail
 */
SFCGAL_API double
distanceLineStringPolygon3D(const LineString &gA, const Polygon &gB);
/**
 * distance between a LineString and a Solid
 * @ingroup detail
 */
SFCGAL_API double
distanceLineStringSolid3D(const LineString &gA, const Solid &gB);

/**
 * dispatch distance between a Triangle and a Geometry
 * @ingroup detail
 */
SFCGAL_API double
distanceTriangleGeometry3D(const Triangle &gA, const Geometry &gB);
/**
 * distance between a Triangle and a Solid
 * @ingroup detail
 */
SFCGAL_API double
distanceTriangleSolid3D(const Triangle &gA, const Solid &gB);

/**
 * dispatch distance between a Polygon and a Geometry
 * @ingroup detail
 */
SFCGAL_API double
distancePolygonGeometry3D(const Polygon &gA, const Geometry &gB);

/**
 * dispatch distance between a Solid and a Geometry
 * @ingroup detail
 */
SFCGAL_API double
distanceSolidGeometry3D(const Solid &gA, const Geometry &gB);
/**
 * distance between two Solids
 * @ingroup detail
 */
SFCGAL_API double
distanceSolidSolid3D(const Solid &gA, const Solid &gB);

/**
 * dispatch distance from a collection of geometry (gA) to a Geometry (gB)
 * @ingroup detail
 */
SFCGAL_API double
distanceGeometryCollectionToGeometry3D(const Geometry &gA, const Geometry &gB);

/**
 * @ingroup detail
 */
SFCGAL_API double
distancePointSegment3D(const Point &p, const Point &a, const Point &b);
/**
 * @ingroup detail
 */
SFCGAL_API double
distancePointTriangle3D(const Point &p, const Point &a, const Point &b,
                        const Point &c);
/**
 * @ingroup detail
 */
SFCGAL_API double
distanceSegmentSegment3D(const Point &a, const Point &b, const Point &c,
                         const Point &d);
/**
 * @ingroup detail
 */
SFCGAL_API double
distanceSegmentTriangle3D(const Point &sA, const Point &sB, const Point &tA,
                          const Point &tB, const Point &tC);
/**
 * distance between two Triangles
 * @ingroup detail
 */
SFCGAL_API double
distanceTriangleTriangle3D(const Triangle &gA, const Triangle &gB);

} // namespace algorithm
} // namespace SFCGAL

#endif

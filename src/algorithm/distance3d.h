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
 * dispatch distance between two Geometries
 * @todo complete with solid
 * Compute distance between two 3D Geometries
 * @pre gA is a valid geometry
 * @pre gB is a valid geometry
 */
SFCGAL_API double
distance3D(const Geometry &gA, const Geometry &gB);

/**
 * Compute distance between two 3D Geometries
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
 * distance between a Point and a PolyhedralSurface
 * @ingroup detail
 */
SFCGAL_API auto
distancePointPolyhedralSurface3D(const Point             &pointA,
                                 const PolyhedralSurface &polySurfaceB)
    -> double;
/**
 * distance between a Point and a TriangulatedSurface
 * @ingroup detail
 */
SFCGAL_API auto
distancePointTriangulatedSurface3D(
    const Point &pointA, const TriangulatedSurface &triangulatedSurfaceB)
    -> double;
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
 * distance between a LineString and a PolyhedralSurface
 * @ingroup detail
 */
SFCGAL_API auto
distanceLineStringPolyhedralSurface3D(const LineString        &lineA,
                                      const PolyhedralSurface &polySurfaceB)
    -> double;
/**
 * distance between a LineString and a TriangulatedSurface
 * @ingroup detail
 */
SFCGAL_API auto
distanceLineStringTriangulatedSurface3D(
    const LineString &lineA, const TriangulatedSurface &triangulatedSurfaceB)
    -> double;
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
 * distance between a Triangle and a PolyhedralSurface
 * @ingroup detail
 */
SFCGAL_API double
distanceTrianglePolyhedralSurface3D(const Triangle          &triangleA,
                                    const PolyhedralSurface &polySurfaceB);
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
 * dispatch distance between a PolyhedralSurface and a Geometry
 * @ingroup detail
 */
SFCGAL_API auto
distancePolyhedralSurfaceGeometry3D(const PolyhedralSurface &polySurfaceA,
                                    const Geometry          &geomB) -> double;

/**
 * dispatch distance between a TriangulatedSurface and a Geometry
 * @ingroup detail
 */
SFCGAL_API auto
distanceTriangulatedSurfaceGeometry3D(
    const TriangulatedSurface &triangulatedSurfaceA, const Geometry &geomB)
    -> double;

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

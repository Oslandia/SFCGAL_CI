// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_ALGORITHM_STRAIGHTSKELETON_H_
#define SFCGAL_ALGORITHM_STRAIGHTSKELETON_H_

#include "SFCGAL/PolyhedralSurface.h"
#include "SFCGAL/algorithm/extrude.h"
#include "SFCGAL/algorithm/union.h"
#include "SFCGAL/config.h"
#include <memory>

namespace SFCGAL {
class Geometry;
class Polygon;
class MultiPolygon;
class MultiLineString;
} // namespace SFCGAL

namespace SFCGAL {
namespace algorithm {
struct NoValidityCheck;

/**
 * @brief build an approximate medial axis for a Polygon
 * @param geom input geometry
 * @pre geom is a valid geometry
 * @throws NotImplementedException If geom is a Polygon with point touching
 * rings.
 */
SFCGAL_API std::unique_ptr<MultiLineString>
           approximateMedialAxis(const Geometry &geom);

/**
 * @brief build a 2D straight skeleton for a Polygon
 * @todo add supports for TriangulatedSurface and PolyhedralSurface
 * @param geom input geometry
 * @param autoOrientation check and fix polygon orientation
 * @param innerOnly Skip non-inner edges if requested
 * @param outputDistanceInM whether to output the distance to border as M
 * @param toleranceAbs Distance tolerance between returned points. A line must
 * have a maximum distance of toleranceAbs.
 * @pre geom is a valid geometry
 * @throws NotImplementedException If geom is a Polygon with point touching
 * rings.
 */
SFCGAL_API std::unique_ptr<MultiLineString>
           straightSkeleton(const Geometry &geom, bool autoOrientation = true,
                            bool innerOnly = false, bool outputDistanceInM = false,
                            const double &toleranceAbs = EPSILON);

/**
 * @brief build a 2D straight skeleton for a Polygon
 * @param geom input geometry
 * @param autoOrientation check and fix polygon orientation
 * @param innerOnly Skip non-inner edges if requested
 * @param outputDistanceInM whether to output the distance to border as M
 * @param toleranceAbs Distance tolerance between returned points. A line must
 * have a maximum distance of toleranceAbs.
 * @pre geom is a valid geometry
 * @warning No actual validity check is done
 * @throws NotImplementedException If geom is a Polygon with point touching
 * rings.
 */
SFCGAL_API std::unique_ptr<MultiLineString>
straightSkeleton(const Geometry &geom, bool autoOrientation, NoValidityCheck,
                 bool innerOnly = false, bool outputDistanceInM = false,
                 const double &toleranceAbs = EPSILON);

/**
 * @brief build a 2D straight skeleton for a Polygon
 * @ingroup detail
 * @throws NotImplementedException If geom is a Polygon with point touching
 * rings.
 */
SFCGAL_API std::unique_ptr<MultiLineString>
           straightSkeleton(const Polygon &geom, bool autoOrientation = true,
                            bool innerOnly = false, bool outputDistanceInM = false,
                            const double &toleranceAbs = EPSILON);

/**
 * @brief build a 2D straight skeleton for a Polygon
 * @ingroup detail
 * @throws NotImplementedException If geom is a Polygon with point touching
 * rings.
 */
SFCGAL_API std::unique_ptr<MultiLineString>
straightSkeleton(const MultiPolygon &geom, bool autoOrientation = true,
                 bool innerOnly = false, bool outputDistanceInM = false,
                 const double &toleranceAbs = EPSILON);

/**
 * @brief build a 3D straight skeleton extruded for a Polygon
 * @ingroup detail
 * @throws NotImplementedException If geom is a Polygon with point touching
 * rings.
 */
SFCGAL_API auto
extrudedStraightSkeleton(const Polygon &geom, double height)
    -> std::unique_ptr<PolyhedralSurface>;

SFCGAL_API auto
extrudeStraightSkeleton(const Geometry &geom, double height)
    -> std::unique_ptr<PolyhedralSurface>;

SFCGAL_API auto
extrudeStraightSkeleton(const Geometry &geom, double building_height,
                        double roof_height)
    -> std::unique_ptr<PolyhedralSurface>;

/**
 * @brief Build a 2D straight skeleton partition for a Geometry
 * @param[in] geom The input geometry
 * @param[in] autoOrientation Check and fix polygon orientation
 * @return A unique pointer to a MultiPolygon representing the partitioned
 * geometry
 * @throws Exception If CGAL fails to create the straight skeleton
 * @note Only Triangle, Polygon, and MultiPolygon geometries are supported
 *
 * This function creates a partition of the input geometry based on its straight
 * skeleton. For unsupported geometry types, an empty MultiPolygon is returned.
 */
SFCGAL_API auto
straightSkeletonPartition(const Geometry &geom, bool autoOrientation = true)
    -> std::unique_ptr<PolyhedralSurface>;

/**
 * @brief Build a 2D straight skeleton partition for a Polygon
 * @ingroup detail
 * @param[in] geom The input polygon
 * @param[in] autoOrientation Check and fix polygon orientation (not used in
 * this implementation)
 * @return A unique pointer to a MultiPolygon representing the partitioned
 * polygon
 * @throws Exception If CGAL fails to create the straight skeleton
 *
 * This function creates a partition of the input polygon based on its straight
 * skeleton. It uses CGAL's Arrangement_2 to handle the intersection of skeleton
 * and polygon edges.
 */
SFCGAL_API auto
straightSkeletonPartition(const Polygon &geom, bool autoOrientation = true)
    -> std::unique_ptr<PolyhedralSurface>;

/**
 * @brief Build a 2D straight skeleton partition for a MultiPolygon
 * @ingroup detail
 * @param[in] geom The input multi-polygon
 * @param[in] autoOrientation Check and fix polygon orientation
 * @return A unique pointer to a MultiPolygon representing the partitioned
 * multi-polygon
 *
 * This function applies the straight skeleton partition to each polygon in the
 * input multi-polygon and combines the results into a single MultiPolygon.
 */
SFCGAL_API auto
straightSkeletonPartition(const MultiPolygon &geom, bool autoOrientation = true)
    -> std::unique_ptr<PolyhedralSurface>;

} // namespace algorithm
} // namespace SFCGAL

#endif

// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2025, SFCGAL team.
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

namespace SFCGAL::algorithm {
struct NoValidityCheck;

/**
 * @brief build an approximate medial axis for a Polygon
 * @param geom input geometry
 * @param projectToEdges if true, project free endpoints to polygon boundary
 *                       using edge midpoint method
 * @return approximate medial axis as a MultiLineString
 * @pre geom is a valid geometry
 * @throws NotImplementedException If geom is a Polygon with point touching
 * rings.
 *
 * When projectToEdges is true, free endpoints of the medial axis are extended
 * to the polygon boundary. The projection uses the midpoint of the common
 * defining edge from the straight skeleton, providing a geometrically exact
 * and symmetric result.
 */
SFCGAL_API auto
approximateMedialAxis(const Geometry &geom, bool projectToEdges = false)
    -> std::unique_ptr<MultiLineString>;

/**
 * @brief build a 2D straight skeleton for a Polygon
 * @todo add supports for TriangulatedSurface and PolyhedralSurface
 * @param geom input geometry
 * @param autoOrientation check and fix polygon orientation
 * @param innerOnly Skip non-inner edges if requested
 * @param outputDistanceInM whether to output the distance to border as M
 * @param toleranceAbs Distance tolerance between returned points. A line must
 * have a maximum distance of toleranceAbs.
 * @return 2D straight skeleton as a MultiLineString
 * @pre geom is a valid geometry
 * @throws NotImplementedException If geom is a Polygon with point touching
 * rings.
 */
SFCGAL_API auto
straightSkeleton(const Geometry &geom, bool autoOrientation = true,
                 bool innerOnly = false, bool outputDistanceInM = false,
                 const double &toleranceAbs = EPSILON)
    -> std::unique_ptr<MultiLineString>;

/**
 * @brief build a 2D straight skeleton for a Polygon
 * @param geom input geometry
 * @param autoOrientation check and fix polygon orientation
 * @param innerOnly Skip non-inner edges if requested
 * @param outputDistanceInM whether to output the distance to border as M
 * @param toleranceAbs Distance tolerance between returned points. A line must
 * have a maximum distance of toleranceAbs.
 * @return 2D straight skeleton as a MultiLineString
 * @pre geom is a valid geometry
 * @warning No actual validity check is done
 * @throws NotImplementedException If geom is a Polygon with point touching
 * rings.
 */
SFCGAL_API auto
straightSkeleton(const Geometry &geom, bool autoOrientation, NoValidityCheck,
                 bool innerOnly = false, bool outputDistanceInM = false,
                 const double &toleranceAbs = EPSILON)
    -> std::unique_ptr<MultiLineString>;

/**
 * @brief build a 2D straight skeleton for a Polygon
 * @param geom input polygon
 * @param autoOrientation check and fix polygon orientation
 * @param innerOnly Skip non-inner edges if requested
 * @param outputDistanceInM whether to output the distance to border as M
 * @param toleranceAbs Distance tolerance between returned points
 * @return 2D straight skeleton as a MultiLineString
 * @throws NotImplementedException If geom is a Polygon with point touching
 * rings.
 */
SFCGAL_API auto
straightSkeleton(const Polygon &geom, bool autoOrientation = true,
                 bool innerOnly = false, bool outputDistanceInM = false,
                 const double &toleranceAbs = EPSILON)
    -> std::unique_ptr<MultiLineString>;

/**
 * @brief build a 2D straight skeleton for a Polygon
 * @param geom input multi-polygon
 * @param autoOrientation check and fix polygon orientation
 * @param innerOnly Skip non-inner edges if requested
 * @param outputDistanceInM whether to output the distance to border as M
 * @param toleranceAbs Distance tolerance between returned points
 * @return 2D straight skeleton as a MultiLineString
 * @throws NotImplementedException If geom is a Polygon with point touching
 * rings.
 */
SFCGAL_API auto
straightSkeleton(const MultiPolygon &geom, bool autoOrientation = true,
                 bool innerOnly = false, bool outputDistanceInM = false,
                 const double &toleranceAbs = EPSILON)
    -> std::unique_ptr<MultiLineString>;

/**
 * @brief build a 3D straight skeleton extruded for a Polygon
 * @param geom the input polygon
 * @param height extrusion height
 * @param angles vector of vector of angles for each polygon ring edge
 * @return extruded straight skeleton as a PolyhedralSurface
 * @throws NotImplementedException If geom is a Polygon with point touching
 * rings.
 */
SFCGAL_API auto
extrudeStraightSkeleton(const Polygon &geom, double height,
                        std::vector<std::vector<Kernel::FT>> angles = {{}})
    -> std::unique_ptr<PolyhedralSurface>;

/**
 * @brief build a 3D straight skeleton extruded for a Geometry
 * @param geom input geometry
 * @param height extrusion height
 * @param angles vector of vector of angles for each polygon ring edge
 * @return extruded straight skeleton as a PolyhedralSurface
 * @throws NotImplementedException If geom is a Polygon with point touching
 * rings.
 */
SFCGAL_API auto
extrudeStraightSkeleton(const Geometry &geom, double height,
                        std::vector<std::vector<Kernel::FT>> angles = {{}})
    -> std::unique_ptr<PolyhedralSurface>;

/**
 * @brief build a 3D straight skeleton extruded for a Geometry with building and
 * roof heights
 * @param geom input geometry
 * @param building_height building height
 * @param roof_height roof height
 * @param angles vector of vector of angles for each polygon ring edge
 * @return extruded straight skeleton as a PolyhedralSurface
 * @throws NotImplementedException If geom is a Polygon with point touching
 * rings.
 */
SFCGAL_API auto
extrudeStraightSkeleton(const Geometry &geom, double building_height,
                        double                               roof_height,
                        std::vector<std::vector<Kernel::FT>> angles = {{}})
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

} // namespace SFCGAL::algorithm

#endif

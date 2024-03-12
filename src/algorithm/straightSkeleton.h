// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2022, Oslandia.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef _SFCGAL_ALGORITHM_STRAIGHTSKELETON_H_
#define _SFCGAL_ALGORITHM_STRAIGHTSKELETON_H_

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
 * @param g input geometry
 * @ingroup public_api
 * @pre g is a valid geometry
 * @throws NotImplementedException If g is a Polygon with point touching rings.
 */
SFCGAL_API std::unique_ptr<MultiLineString>
           approximateMedialAxis(const Geometry &g);

/**
 * @brief build a 2D straight skeleton for a Polygon
 * @todo add supports for TriangulatedSurface and PolyhedralSurface
 * @param g input geometry
 * @param autoOrientation check and fix polygon orientation
 * @param outputM whether to output the distance to border as M
 * @param toleranceAbs Distance tolerance between returned points. A line must
 * have a maximum distance of toleranceAbs.
 * @ingroup public_api
 * @pre g is a valid geometry
 * @throws NotImplementedException If g is a Polygon with point touching rings.
 */
SFCGAL_API std::unique_ptr<MultiLineString>
           straightSkeleton(const Geometry &g, bool autoOrientation = true,
                            bool innerOnly = false, bool outputDistanceInM = false,
                            const double &toleranceAbs = 1e-8);

/**
 * @brief build a 2D straight skeleton for a Polygon
 * @param g input geometry
 * @param autoOrientation check and fix polygon orientation
 * @param outputM whether to output the distance to border as M
 * @param toleranceAbs Distance tolerance between returned points. A line must
 * have a maximum distance of toleranceAbs.
 * @ingroup public_api
 * @pre g is a valid geometry
 * @warning No actual validity check is done
 * @throws NotImplementedException If g is a Polygon with point touching rings.
 */
SFCGAL_API std::unique_ptr<MultiLineString>
straightSkeleton(const Geometry &g, bool autoOrientation, NoValidityCheck,
                 bool innerOnly = false, bool outputDistanceInM = false,
                 const double &toleranceAbs = 1e-8);

/**
 * @brief build a 2D straight skeleton for a Polygon
 * @ingroup detail
 * @throws NotImplementedException If g is a Polygon with point touching rings.
 */
SFCGAL_API std::unique_ptr<MultiLineString>
           straightSkeleton(const Polygon &g, bool autoOrientation = true,
                            bool innerOnly = false, bool outputDistanceInM = false,
                            const double &toleranceAbs = 1e-8);

/**
 * @brief build a 2D straight skeleton for a Polygon
 * @ingroup detail
 * @throws NotImplementedException If g is a Polygon with point touching rings.
 */
SFCGAL_API std::unique_ptr<MultiLineString>
           straightSkeleton(const MultiPolygon &g, bool autoOrientation = true,
                            bool innerOnly = false, bool outputDistanceInM = false,
                            const double &toleranceAbs = 1e-8);

/**
 * @brief build a 3D straight skeleton extruded for a Polygon
 * @ingroup detail
 * @throws NotImplementedException If g is a Polygon with point touching rings.
 */
SFCGAL_API auto
extrudedStraightSkeleton(const Polygon &g, double height)
    -> std::unique_ptr<PolyhedralSurface>;

SFCGAL_API auto
extrudeStraightSkeleton(const Geometry &g, double height)
    -> std::unique_ptr<PolyhedralSurface>;

SFCGAL_API auto
extrudeStraightSkeleton(const Geometry &g, double building_height,
                        double roof_height)
    -> std::unique_ptr<PolyhedralSurface>;

} // namespace algorithm
} // namespace SFCGAL

#endif
